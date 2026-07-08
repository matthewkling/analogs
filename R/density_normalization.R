# Density normalization helpers ------------------------------------
#
# Implements the global D_max scalar that turns raw `sum_weights` (kernel-
# weighted analog count) and `tabulate` (kernel-weighted per-class count)
# into normalized fractions on roughly [0, 1]. The framing -- worked out at
# length in earlier design discussions -- is "Flavor C, idealized grid":
#
#   D_max = (1 / Ā) · ∫_0^max_geog  K(0, r) · 2 π r  dr
#
# where Ā is the mean cell area (in units that match max_geog for the
# index's coord_type; see `.compute_cell_area_weights()` for unit details)
# and K(0, r) is the kernel evaluated at climate-distance = 0 and the given
# geographic distance.
#
# Conceptually this is the kernel-weight-per-cell that a hypothetical focal
# would see at the interior of an infinite uniform raster with the same
# cell size and projection, where every cell within max_geog is a perfect
# climate match. It's a property of the kernel + max_geog + mean cell area
# only, so it's a single global scalar per query (not per-focal). Discrete
# real focals fluctuate around D / D_max ≤ 1 due to finite-cell
# discretization noise; this is the same noise inherent in the raw
# sum_weights values themselves and is not introduced by normalization.
#
# All integrals below are 1D and have closed forms. They evaluate at zero
# environmental distance, so:
#   - environment-only kernels (gaussian_env, inverse_env) reduce to a
#     constant-K integral over the disk
#   - geo-only kernels (gaussian_geog, inverse_geog) keep the geographic
#     factor unchanged at env_dist = 0
#   - joint kernels (gaussian_joint, inverse_joint) collapse to their
#     geo-only form at env_dist = 0
#   - uniform: D_max = pi * max_geog^2 / mean_cell_area
#
# All kernel types are supported. For environment-aware kernels (gaussian_env,
# inverse_env, gaussian_joint, inverse_joint), D / D_max measures
# eivironment-quality-weighted analog availability. For uniform and geo-only
# kernels, D / D_max measures the fraction of the disk satisfying max_env
# (when max_env is active); without max_env, those reduce to D == D_max
# trivially.


# Compute the global D_max scalar for a given kernel/theta/max_geog/mean
# cell area configuration. Returns a positive finite scalar.
#
# Validation of compatibility (kernel non-NULL, finite max_geog,
# raster-derived mean_cell_area available, etc.) is the caller's
# responsibility -- typically done via .validate_normalize_compat() before
# this is reached.
.compute_global_dmax <- function(kernel_env, kernel_geog,
                                 theta_env, theta_geog,
                                 max_geog, mean_cell_area) {

      if (!is.finite(max_geog) || max_geog <= 0) {
            stop("Internal error: .compute_global_dmax called with invalid max_geog.",
                 call. = FALSE)
      }
      if (!is.finite(mean_cell_area) || mean_cell_area <= 0) {
            stop("Internal error: .compute_global_dmax called with invalid mean_cell_area.",
                 call. = FALSE)
      }

      # D_max is the theoretical maximum aggregate kernel weight for a focal
      # whose entire geographic disk (radius max_geog) consists of perfect
      # climate matches. Per cell the weight is
      #     w_env(env_dist = 0) * w_geog(geo_dist),
      # and w_env(0) = 1 for every kernel shape (uniform, Gaussian exp(0)=1,
      # and the reparameterized inverse 1/(1 + 0/theta) = 1). So the environmental
      # kernel contributes a constant factor of 1 and D_max depends only on the
      # GEOGRAPHIC kernel, integrated over the disk:
      #     D_max = [ \int_0^G w_geog(r) 2*pi*r dr ] / mean_cell_area.
      # This holds for any composition of per-family kernels (the environmental shape
      # is irrelevant at env_dist = 0), so no fused-kernel reconstruction is
      # needed. The disk integral (with closed forms per shape) lives in the
      # shared helper .kernel_disk_integral().
      geo_theta <- if (is.null(theta_geog) ||
                       (length(theta_geog) == 1L && is.na(theta_geog))) {
            NULL
      } else {
            theta_geog
      }

      kernel_integral <- .kernel_disk_integral(kernel_geog, geo_theta, max_geog)

      # Convert "kernel weight integrated over physical area" into "kernel
      # weight per cell" by dividing by mean cell area. After this division
      # D_max is on the same dimensionless "per cell" scale as `sum_weights`
      # for cells with cell_area_weight (mean 1) applied.
      D_max <- kernel_integral / mean_cell_area

      if (!is.finite(D_max) || D_max <= 0) {
            stop("Internal error: computed D_max is non-positive or non-finite ",
                 "(got ", D_max, "). Check kernel_geog/theta_geog/max_geog/",
                 "mean_cell_area.", call. = FALSE)
      }

      D_max
}


# Helper: would `normalize = TRUE` succeed AND have any effect given these
# query parameters and this index? Used by `auto` resolution to decide
# silently between TRUE and FALSE without erroring. Returns logical(1).
#
# Mirrors the index-side and query-side preconditions enforced by
# .validate_normalize_compat(), but as a predicate rather than a validator.
# Returns FALSE when there are no normalizable stats so that `auto`
# resolves to FALSE for those queries -- this keeps the result attribute
# (`normalize = FALSE`) honest and skips downstream D_max computation.
.normalize_is_feasible <- function(stat, max_geog, index) {
      # Stat-relevance: if no column would be normalized, normalize would
      # be a no-op. Treat as not-feasible so `auto` resolves to FALSE.
      if (!any(stat %in% c("sum_weights", "tabulate"))) return(FALSE)

      # Index-side prereqs: built from a raster pool with cell-area
      # weighting active.
      if (is.null(index$mean_cell_area)) return(FALSE)
      if (is.null(index$cell_area_weight)) return(FALSE)

      # Query-side prereq: finite max_geog (bounds the disk over which D_max is
      # integrated). Any kernel configuration is supported -- including a
      # uniform kernel, where D_max is the disk area and D / D_max is the
      # fraction of the disk satisfying the constraints. D_max depends only on
      # the geographic kernel (the climate factor is 1 at env_dist = 0), so a
      # non-uniform kernel is NOT required.
      if (is.null(max_geog) || !is.finite(max_geog) || max_geog <= 0) return(FALSE)

      TRUE
}


# Validate that the requested `normalize = TRUE` configuration is
# well-defined for this index/query. Throws an informative error if not.
#
# Preconditions checked:
#   - kernel is non-NULL and not "uniform" (D_max would equal D)
#   - kernel is climate-aware: gaussian_env, inverse_env, gaussian_joint,
#     inverse_joint. Geo-only kernels (gaussian_geog, inverse_geog) are
#     rejected because for them K(0, r) = K(d_env, r) for all d_env,
#     so D_max == D up to discretization, which is meaningless.
#   - max_geog is finite and positive (a finite disk is required for a
#     finite D_max)
#   - the index has a non-NULL `mean_cell_area` (built from a raster pool
#     with cell-area weighting active)
#   - the index has cell_area_weight active (otherwise the normalization
#     denominator and the sum_weights numerator are on different scales)
#   - at least one of the "sum_weights" / "tabulate" stats is requested
#     (other stats are unaffected by normalize; warn but don't error)
#
# Returns invisible(TRUE) when the configuration is valid; otherwise
# stop()s with a user-facing message.
.validate_normalize_compat <- function(normalize, stat, max_geog,
                                       index) {
      if (!isTRUE(normalize)) return(invisible(TRUE))

      # 0. Stat-relevance short-circuit. If `stat` doesn't include any
      # column that normalization would apply to, normalize is a no-op
      # and we don't validate any further preconditions. This matches the
      # documented behavior ("silently ignored when the requested `stat`
      # includes none of `sum_weights` or `tabulate`") and avoids spurious
      # errors for callers that pass normalize = TRUE through generic
      # plumbing while requesting unrelated stats (e.g. analog_velocity's
      # stat = "none", analog_availability's stat = "count").
      if (!any(stat %in% c("sum_weights", "tabulate"))) {
            return(invisible(TRUE))
      }

      # 1. max_geog must be finite and positive. (Any kernel configuration is
      # supported, including uniform: D_max depends only on the geographic
      # kernel and is well-defined for all shapes, so no non-uniform kernel is
      # required.)
      if (is.null(max_geog) || !is.finite(max_geog) || max_geog <= 0) {
            stop("`normalize = TRUE` requires a finite positive `max_geog`. ",
                 "Without a geographic constraint, the D_max integral is ",
                 "undefined for some kernels and uninterpretable for others.",
                 call. = FALSE)
      }

      # 2. Index must have mean_cell_area available (from a raster pool).
      if (is.null(index$mean_cell_area)) {
            stop("`normalize = TRUE` requires that the index was built from a ",
                 "SpatRaster pool, so that a mean cell area is available for ",
                 "the D_max integral. Rebuild your index from a raster (or ",
                 "supply `mean_cell_area` explicitly to `build_analog_index()`).",
                 call. = FALSE)
      }

      # 3. Index must have cell-area weighting active. Otherwise the numerator
      # and denominator are on inconsistent scales (sum_weights computed
      # without per-cell area normalization, but D_max derived assuming
      # mean-1 area weights). We require the user to opt in explicitly.
      if (is.null(index$cell_area_weight)) {
            stop("`normalize = TRUE` requires `cell_area_weight` to be active ",
                 "on the index (it is by default for raster pools). Rebuild ",
                 "the index with `cell_area_weight = 'auto'` or `TRUE`, or ",
                 "use `normalize = FALSE`.", call. = FALSE)
      }

      # 4. Stat compatibility. Only sum_weights and tabulate get normalized;
      # for other stats `normalize = TRUE` is a silent no-op (per
      # documentation). We don't warn here, since `normalize = TRUE` is the
      # default for some helpers and users may legitimately request stats
      # that don't happen to include a normalizable column.

      invisible(TRUE)
}


# Apply D_max normalization in-place to the relevant columns of an
# aggregation-mode result data.frame. Currently this means dividing
# `sum_weights` and any tabulate-class columns (named `n_<lev>` or
# `<varname>_n_<lev>`) by D_max.
#
# `out` is the post-rebuild data.frame produced by query_analog_index()
# right before .format_output(). It has columns `index`, `x`, `y` plus one
# column per stat (with multi-y suffixes where applicable).
#
# `tabulate_col_names` is the character vector of column names produced for
# the `tabulate` stat (assembled by query_analog_index() during column
# layout); we receive it explicitly because tabulate column names follow a
# convention rather than a fixed pattern.
#
# Returns the modified data.frame.
.apply_dmax_normalization <- function(out, D_max, stat, tabulate_col_names) {
      if (!is.finite(D_max) || D_max <= 0) {
            stop("Internal error: cannot normalize by non-positive D_max (",
                 D_max, ").", call. = FALSE)
      }

      if ("sum_weights" %in% stat && "sum_weights" %in% names(out)) {
            out[["sum_weights"]] <- out[["sum_weights"]] / D_max
      }

      if ("tabulate" %in% stat && length(tabulate_col_names) > 0L) {
            present <- intersect(tabulate_col_names, names(out))
            for (cn in present) {
                  out[[cn]] <- out[[cn]] / D_max
            }
      }

      out
}
