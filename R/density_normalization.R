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
# climate distance, so:
#   - climate-only kernels (gaussian_clim, inverse_clim) reduce to a
#     constant-K integral over the disk
#   - geo-only kernels (gaussian_geog, inverse_geog) keep the geographic
#     factor unchanged at clim_dist = 0
#   - joint kernels (gaussian_joint, inverse_joint) collapse to their
#     geo-only form at clim_dist = 0
#   - uniform: D_max = pi * max_geog^2 / mean_cell_area
#
# All kernel types are supported. For climate-aware kernels (gaussian_clim,
# inverse_clim, gaussian_joint, inverse_joint), D / D_max measures
# climate-quality-weighted analog availability. For uniform and geo-only
# kernels, D / D_max measures the fraction of the disk satisfying max_clim
# (when max_clim is active); without max_clim, those reduce to D == D_max
# trivially.


# Compute the global D_max scalar for a given kernel/theta/max_geog/mean
# cell area configuration. Returns a positive finite scalar.
#
# Validation of compatibility (kernel non-NULL, finite max_geog,
# raster-derived mean_cell_area available, etc.) is the caller's
# responsibility -- typically done via .validate_normalize_compat() before
# this is reached.
.compute_global_dmax <- function(kernel, theta, max_geog, mean_cell_area) {

      if (!is.finite(max_geog) || max_geog <= 0) {
            stop("Internal error: .compute_global_dmax called with invalid max_geog.",
                 call. = FALSE)
      }
      if (!is.finite(mean_cell_area) || mean_cell_area <= 0) {
            stop("Internal error: .compute_global_dmax called with invalid mean_cell_area.",
                 call. = FALSE)
      }

      # Defaults for theta (mirroring C++-side defaults documented in
      # ?analog_search). For kernels that don't use theta in the integrand
      # below (the climate-only and uniform cases), this is a no-op.
      eps_clim_default <- 1e-12
      eps_geog_default <- 1e-6
      sigma_default    <- 1.0   # for single-dim Gaussians; not used at clim=0
      theta_clim_default <- 1.0
      theta_geog_default <- 1.0

      # Resolve theta into named scalars per kernel for clarity below.
      th <- if (is.null(theta) || (length(theta) == 1L && is.na(theta))) {
            NULL
      } else {
            as.numeric(theta)
      }

      # Integrand-specific calculations.
      # All integrals are over r in [0, max_geog] in the planar disk, so the
      # 2-D-area element is 2 π r dr.
      G <- max_geog                # alias for readability
      twopi <- 2 * pi

      kernel_integral <- switch(
            kernel,

            "uniform" = {
                  # K(0, r) = 1 for all r (no kernel weighting at all).
                  # ∫₀^G 1 · 2π r dr = π G². With max_clim filtering active,
                  # D / D_max gives the fraction of the disk satisfying
                  # max_clim; without max_clim, D == D_max trivially.
                  pi * G * G
            },

            "gaussian_clim" = {
                  # K(0, r) = exp(-0 / 2σ²) = 1 for all r.
                  # ∫₀^G 1 · 2π r dr = π G².
                  pi * G * G
            },

            "inverse_clim" = {
                  # K(0, r) = 1 / (0 + eps) = 1 / eps (theta = epsilon).
                  eps <- if (!is.null(th)) th[1] else eps_clim_default
                  if (!is.finite(eps) || eps <= 0) {
                        stop("Internal error: invalid epsilon for inverse_clim.",
                             call. = FALSE)
                  }
                  (1 / eps) * pi * G * G
            },

            "gaussian_geog" = {
                  # K(0, r) = exp(-r² / 2σ²). With max_clim filtering active,
                  # D / D_max gives the geo-weighted fraction of the disk
                  # satisfying max_clim; without max_clim, D == D_max
                  # trivially since this kernel doesn't depend on climate
                  # distance.
                  sigma <- if (!is.null(th)) th[1] else sigma_default
                  if (!is.finite(sigma) || sigma <= 0) {
                        stop("Internal error: invalid sigma for gaussian_geog.",
                             call. = FALSE)
                  }
                  # ∫₀^G exp(-r²/(2σ²)) · 2π r dr
                  #   = 2πσ² · (1 - exp(-G²/(2σ²)))
                  twopi * sigma^2 * (1 - exp(-G^2 / (2 * sigma^2)))
            },

            "inverse_geog" = {
                  # K(0, r) = 1 / (r + eps). Closed form for the disk
                  # integral:
                  #   ∫₀^G (r / (r + eps)) · 2π dr
                  #     = 2π · (G - eps · ln((G + eps) / eps))
                  eps <- if (!is.null(th)) th[1] else eps_geog_default
                  if (!is.finite(eps) || eps <= 0) {
                        stop("Internal error: invalid epsilon for inverse_geog.",
                             call. = FALSE)
                  }
                  twopi * (G - eps * log((G + eps) / eps))
            },

            "gaussian_joint" = {
                  # K(0, r) = exp(-(0 + r²/(2 σ_g²))) -> same as gaussian_geog.
                  # theta = c(theta_clim, theta_geog).
                  sigma_g <- if (!is.null(th)) th[2] else theta_geog_default
                  if (!is.finite(sigma_g) || sigma_g <= 0) {
                        stop("Internal error: invalid theta_geog for gaussian_joint.",
                             call. = FALSE)
                  }
                  twopi * sigma_g^2 * (1 - exp(-G^2 / (2 * sigma_g^2)))
            },

            "inverse_joint" = {
                  # K(clim_dist, r) = 1 / (sqrt(clim_dist² + r²) + eps).
                  # At clim_dist = 0: 1 / (r + eps). theta = c(theta_clim,
                  # theta_geog); the geographic eps lives in theta[2] per
                  # the C++ convention (matches inverse_geog scale).
                  eps <- if (!is.null(th)) th[2] else theta_geog_default
                  if (!is.finite(eps) || eps <= 0) {
                        stop("Internal error: invalid theta_geog for inverse_joint.",
                             call. = FALSE)
                  }
                  twopi * (G - eps * log((G + eps) / eps))
            },

            stop("Internal error: unrecognized kernel '", kernel,
                 "' in .compute_global_dmax.", call. = FALSE)
      )

      # Convert "kernel weight integrated over physical area" into "kernel
      # weight per cell" by dividing by mean cell area. After this division
      # D_max is on the same dimensionless "per cell" scale as `sum_weights`
      # for cells with cell_area_weight (mean 1) applied.
      D_max <- kernel_integral / mean_cell_area

      if (!is.finite(D_max) || D_max <= 0) {
            stop("Internal error: computed D_max is non-positive or non-finite ",
                 "(got ", D_max, "). Check kernel/theta/max_geog/mean_cell_area.",
                 call. = FALSE)
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
.normalize_is_feasible <- function(stat, kernel, max_geog, index) {
      # Stat-relevance: if no column would be normalized, normalize would
      # be a no-op. Treat as not-feasible so `auto` resolves to FALSE.
      if (!any(stat %in% c("sum_weights", "tabulate"))) return(FALSE)

      # Index-side prereqs: built from a raster pool with cell-area
      # weighting active.
      if (is.null(index$mean_cell_area)) return(FALSE)
      if (is.null(index$cell_area_weight)) return(FALSE)

      # Query-side prereqs: a kernel is set + finite max_geog. All kernel
      # types are supported -- including uniform and geo-only kernels --
      # though for those, D / D_max measures fraction of the disk
      # satisfying max_clim rather than climate-quality-weighted analog
      # availability.
      if (is.null(kernel)) return(FALSE)
      if (is.null(max_geog) || !is.finite(max_geog) || max_geog <= 0) return(FALSE)

      TRUE
}


# Validate that the requested `normalize = TRUE` configuration is
# well-defined for this index/query. Throws an informative error if not.
#
# Preconditions checked:
#   - kernel is non-NULL and not "uniform" (D_max would equal D)
#   - kernel is climate-aware: gaussian_clim, inverse_clim, gaussian_joint,
#     inverse_joint. Geo-only kernels (gaussian_geog, inverse_geog) are
#     rejected because for them K(0, r) = K(d_clim, r) for all d_clim,
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
.validate_normalize_compat <- function(normalize, stat, kernel, max_geog,
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

      # 1. A kernel must be set. All kernel types are supported -- even
      # uniform and geo-only kernels yield a well-defined D_max via the
      # closed-form integral, though for those D/D_max reduces to a
      # fraction-of-disk-satisfying-max_clim measurement when max_clim is
      # active.
      if (is.null(kernel)) {
            stop("`normalize = TRUE` requires a kernel. ",
                 "Set `kernel` to one of: 'uniform', 'gaussian_clim', ",
                 "'inverse_clim', 'gaussian_geog', 'inverse_geog', ",
                 "'gaussian_joint', or 'inverse_joint'.", call. = FALSE)
      }

      # 2. max_geog must be finite and positive.
      if (is.null(max_geog) || !is.finite(max_geog) || max_geog <= 0) {
            stop("`normalize = TRUE` requires a finite positive `max_geog`. ",
                 "Without a geographic constraint, the D_max integral is ",
                 "undefined for some kernels and uninterpretable for others.",
                 call. = FALSE)
      }

      # 3. Index must have mean_cell_area available (from a raster pool).
      if (is.null(index$mean_cell_area)) {
            stop("`normalize = TRUE` requires that the index was built from a ",
                 "SpatRaster pool, so that a mean cell area is available for ",
                 "the D_max integral. Rebuild your index from a raster (or ",
                 "supply `mean_cell_area` explicitly to `build_analog_index()`).",
                 call. = FALSE)
      }

      # 4. Index must have cell-area weighting active. Otherwise the numerator
      # and denominator are on inconsistent scales (sum_weights computed
      # without per-cell area normalization, but D_max derived assuming
      # mean-1 area weights). We require the user to opt in explicitly.
      if (is.null(index$cell_area_weight)) {
            stop("`normalize = TRUE` requires `cell_area_weight` to be active ",
                 "on the index (it is by default for raster pools). Rebuild ",
                 "the index with `cell_area_weight = 'auto'` or `TRUE`, or ",
                 "use `normalize = FALSE`.", call. = FALSE)
      }

      # 5. Stat compatibility. Only sum_weights and tabulate get normalized;
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
