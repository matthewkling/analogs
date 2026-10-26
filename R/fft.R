# Binned-convolution (FFT) execution method ------------------------------------
#
# An alternative backend for `select = "all"` kernel-weighted aggregation.
# Instead of visiting every focal-analog pair, it bins the pool's weighted
# moments onto a regular grid over the active axes (environment and/or
# geography), convolves each moment grid with the kernel via FFTs, and reads
# the convolved grids off at the focals. The C++ side lives in
# src/fft_grid.cpp; this file decides the grid geometry, kernel arrays, and
# which moment "fields" each requested statistic needs.


#' Fast Fourier Transform (FFT) execution method
#'
#' Configures the binned-convolution backend for [analog_search()] and its
#' wrappers. For kernel-weighted summary stats (regression, density, impact,
#' etc.), the FFT method computes an approximate result without visiting
#' every focal-analog pair, which can make it orders of magnitude
#' faster for large problems. Use this function by passing the result as
#' `method = fft_method(...)`, or use `method = "fft"` for the defaults.
#'
#' ## How it works
#'
#' FFT is only compatible with kernel-weighted summary statistics. Rather than
#' separately visiting every focal-analog pair, as in the `"exact"` method,
#' FFT (a multivariate version of binned kernel estimation; Wand
#' 1994) evaluates them for all focals at once using an coarse internal grid.
#' Nodes are spaced at a fraction (`bin`) of each kernel's bandwidth.
#' Each sum (or related value liek regression coefficient) is then computed
#' in three steps:
#'
#' \enumerate{
#'   \item **Bin.** Each pool point's contribution is split among the grid
#'     nodes surrounding it, in proportion to its closeness to each (linear
#'     binning). The grid then holds how much pool "mass" lies near each node.
#'   \item **Convolve.** Each node is assigned the kernel-weighted sum of the
#'     mass at every node, done via the fast Fourier transform by transforming
#'     the grid and the kernel, multiplying them, and transforming back.
#'     The grid is padded by one kernel width of empty nodes along each axis,
#'     so that mass near one edge does not wrap around to the opposite edge.
#'   \item **Read off.** Each focal's value is interpolated from the grid nodes
#'     surrounding its exact location and environment (multilinear
#'     interpolation).
#' }
#'
#' Because the convolved surface is a sum of kernel-shaped bumps, it changes
#' little between nodes spaced well within the kernel's width, so the coarse
#' grid loses little accuracy. The output keeps the resolution of the focal
#' data: two neighboring focal raster cells with different climates read off
#' different parts of the environmental axes and can get quite different
#' results. Statistics that need several sums (for example, a regression
#' needs weighted cross-products of the covariates and `y`) repeat the bin and
#' convolve steps once per sum. For leave-one-out (`exclude_self = TRUE`),
#' each focal's own binned contribution is subtracted after read-off.
#'
#' The grid covers only the region where focals and pool points are within
#' each other's reach, so focals confined to part of the pool's extent cost
#' less. A family's nodes are spaced at `bin * theta` (`bin * max` for a
#' uniform kernel, or `bin * min` if a geographic `min` is smaller than
#' `theta`), shrunk slightly so that nodes fall exactly on the ends of the
#' range.
#'
#' ## Accuracy and cost
#'
#' `bin` controls the trade-off between accuracy and cost. For the Gaussian
#' kernel, binning error shrinks with the square of `bin`; for the inverse
#' kernel (which has a sharp peak at zero distance) and for hard edges
#' (uniform kernels, `max`/`min` cutoffs where the weight is still
#' appreciable, `count`) it shrinks roughly in proportion to `bin`. Hard limits
#' are anti-aliased on the grid (each node counts in proportion to the
#' fraction of its cell inside the limit), which reduces but does not remove
#' this edge error. Halving `bin` multiplies the number of grid nodes by `2^d`
#' for `d` active axes. The default `bin = 0.25` typically gives errors well
#' under 1%, but error at a given `bin` is larger in sparse neighborhoods and
#' for ratio statistics such as regression coefficients; use `check` to
#' measure it.
#'
#' ## When to use it
#'
#' The FFT method's cost is set by the size of its grid rather than by the
#' number of focal-analog pairs. Each axis has roughly `range / (bin * theta)`
#' nodes, where `range` is the data's extent along that axis, and grid size is
#' the product over axes. So what matters is how many axes are active and how
#' narrow the kernels are relative to the data's range:
#'
#' \itemize{
#'   \item **Two active axes** (e.g. two environmental variables, or a
#'     geographic kernel alone): the grid stays small even for narrow kernels,
#'     and the FFT method is typically one to three orders of magnitude faster
#'     than the exact method for large numbers of focals.
#'   \item **Three or more active axes**: grid size grows with the product of
#'     the range-to-bandwidth ratios across axes. The FFT method pays off when
#'     the kernels are wide relative to the data's range. When a kernel is
#'     narrow, analog neighborhoods are small, the exact method's search index
#'     is efficient, and the FFT grid may exceed `max_cells`.
#' }
#'
#' Speedups also shrink with more threads, since the exact method
#' parallelizes over focals more completely than the FFT method. When in
#' doubt, time both on a subset of focals. If `max_cells` is exceeded, use a
#' coarser `bin` or the exact method; for raster inputs,
#' [tiled_analog_search()] accepts `method = "fft"` and builds a separate,
#' smaller grid for each tile and its halo (spanning only the environments
#' found there), which is the most effective way to cut peak memory when
#' geographic axes are active.
#'
#' ## Requirements
#'
#' Geographic kernels require projected coordinates (ideally an equal-area
#' projection), since on a longitude/latitude grid a fixed-distance kernel
#' changes shape with latitude. The FFT method requires `select = "all"` and
#' a raw (non-index) `pool`, and does not support `stat = "none"`, `x_cov`,
#' `downsample < 1`, or design-based standard errors for `"regression"`.
#' Gaussian kernels without a `max` are truncated at 7 `theta`.
#'
#' @param bin Grid spacing, as a fraction of each kernel's scale: `theta` for
#'   Gaussian and inverse kernels, `max` for uniform kernels (see Details).
#'   Default `0.25`.
#' @param check Number of randomly chosen focals to recompute with the exact
#'   engine for an accuracy report, attached to the output as the
#'   `fft_check` attribute (see [metadata()]). Default `0` (no check).
#' @param max_cells Upper limit on the number of grid nodes, as a guard
#'   against memory exhaustion (each node of each grid takes 8 bytes).
#'   Default `1e8`.
#'
#' @return An object of class `"analog_fft_method"`.
#'
#' @references Wand, M. P. (1994). Fast computation of multivariate kernel
#'   estimators. *Journal of Computational and Graphical Statistics*, 3(4),
#'   433-445.
#'
#' @seealso [analog_search()]
#'
#' @examples
#' fft_method()
#' fft_method(bin = 0.1, check = 500)
#'
#' @export
fft_method <- function(bin = 0.25, check = 0, max_cells = 1e8) {
      if (!is.numeric(bin) || length(bin) != 1L || !is.finite(bin) ||
          bin <= 0 || bin > 1) {
            stop("`bin` must be a single number in (0, 1].", call. = FALSE)
      }
      if (!is.numeric(check) || length(check) != 1L || !is.finite(check) ||
          check < 0) {
            stop("`check` must be a single non-negative number.", call. = FALSE)
      }
      if (!is.numeric(max_cells) || length(max_cells) != 1L ||
          !is.finite(max_cells) || max_cells < 1) {
            stop("`max_cells` must be a single positive number.", call. = FALSE)
      }
      structure(list(bin = bin, check = as.integer(check), max_cells = max_cells),
                class = "analog_fft_method")
}

#' @export
print.analog_fft_method <- function(x, ...) {
      cat("<analog_fft_method> bin = ", x$bin, ", check = ", x$check,
          ", max_cells = ", format(x$max_cells, big.mark = ","), "\n", sep = "")
      invisible(x)
}

# Resolve the `method` argument: NULL for the exact engine, otherwise an
# analog_fft_method object.
.resolve_method <- function(method) {
      if (identical(method, "exact")) return(NULL)
      if (identical(method, "fft")) return(fft_method())
      if (inherits(method, "analog_fft_method")) return(method)
      stop('`method` must be "exact", "fft", or an fft_method() object.',
           call. = FALSE)
}


# Query ------------------------------------------------------------------------

# FFT replacement for .query_cpp_chunked(): returns a matrix in the exact
# engine's column layout (one row per non-NA focal). Called from
# query_analog_index() after all inputs are validated and formatted.
# `exact_fn(rows)` runs the exact engine on a subset of focal rows (for
# `check`).
.query_fft <- function(fft, focal, ref, stat, aggregate_codes, se, se_code,
                       values_mat, covariates_mat, n_classes_per_var,
                       area_weight, user_weight, lambda,
                       kernel_env, theta_env, max_env,
                       kernel_geog, theta_geog, max_geog, min_geog,
                       x_cov, exclude_self, coord_type, exact_fn) {

      n_env <- ncol(ref) - 2L
      env <- .fft_family(kernel_env, theta_env, max_env)
      geog <- .fft_family(kernel_geog, theta_geog, max_geog, min_geog)

      # Preconditions (select = "all", pool not an index, and no downsampling
      # are enforced in analog_search()).
      if (identical(stat, "none")) {
            stop('`method = "fft"` does not support `stat = "none"`; ',
                 "use the exact method for analog pairs.", call. = FALSE)
      }
      if (!is.null(x_cov)) {
            stop('`method = "fft"` does not support `x_cov` (focal-specific ',
                 "covariances make the kernel vary between focals).",
                 call. = FALSE)
      }
      if ("regression" %in% stat && se == "design") {
            stop('`method = "fft"` does not support `se = "design"` for ',
                 '"regression"; use `se = "ess"` or the exact method.',
                 call. = FALSE)
      }
      if (geog$active && !identical(coord_type, "projected")) {
            stop('`method = "fft"` with a geographic kernel or limit needs ',
                 "projected coordinates: on a lon/lat grid the kernel's shape ",
                 "varies with latitude. Project the data to an equal-area CRS ",
                 "(e.g. with terra::project()), or use the exact method.",
                 call. = FALSE)
      }
      if (env$active && n_env > 3L) {
            stop('`method = "fft"` supports at most 3 environmental variables ',
                 "(found ", n_env, ").", call. = FALSE)
      }

      # Grid axes: geography first, then environment (each only if active).
      fams <- list()
      if (geog$active) fams$geog <- c(geog, list(cols = 1:2))
      if (env$active)  fams$env  <- c(env,  list(cols = 2L + seq_len(n_env)))

      # Fields, then kernels for the kernel types the fields use.
      fields <- .fft_build_fields(stat, se, values_mat, covariates_mat,
                                  n_classes_per_var, area_weight, user_weight,
                                  n_pool = nrow(ref))
      axes <- lapply(fams, .fft_family_axes,
                     pool = ref, focal = focal, bin = fft$bin,
                     kernel_types = fields$kernel_types)

      n_nodes <- unlist(lapply(axes, `[[`, "n"), use.names = FALSE)
      pad <- unlist(lapply(axes, `[[`, "pad"), use.names = FALSE)
      if (prod(as.numeric(n_nodes)) > fft$max_cells) {
            stop("The FFT grid would have ", format(prod(as.numeric(n_nodes)), big.mark = ","),
                 " nodes, more than `max_cells` (",
                 format(fft$max_cells, big.mark = ","), "). Increase `bin` in ",
                 "fft_method(), raise `max_cells`, run raster inputs through ",
                 "tiled_analog_search(), or use the exact method.",
                 call. = FALSE)
      }
      cols <- unlist(lapply(fams, `[[`, "cols"), use.names = FALSE)
      kernels <- lapply(fields$kernel_types, function(kt) {
            lapply(axes, function(ax) ax$kernels[[kt]])
      })

      vals <- .fft_fields_cpp(
            pool_pts     = ref[, cols, drop = FALSE],
            focal_pts    = focal[, cols, drop = FALSE],
            origin       = unlist(lapply(axes, `[[`, "origin"), use.names = FALSE),
            delta        = unlist(lapply(axes, `[[`, "delta"), use.names = FALSE),
            n_nodes      = as.integer(n_nodes),
            group_sizes  = as.integer(vapply(fams, function(f) length(f$cols), 1L)),
            pad_sizes    = as.integer(pad),
            kernels      = kernels,
            masses       = fields$masses,
            field_kernel = fields$field_kernel,
            exclude_self = isTRUE(exclude_self)
      )

      uw_mean <- if (is.null(user_weight) || !any(user_weight > 0)) 1 else
            mean(user_weight[user_weight > 0])

      res <- .fft_finalize_cpp(
            field_vals        = vals,
            roles             = fields$roles,
            aggregate_codes   = as.integer(aggregate_codes),
            se_code           = as.integer(se_code),
            n_vars            = if (is.null(values_mat)) 0L else ncol(values_mat),
            n_covs            = if (is.null(covariates_mat)) 0L else ncol(covariates_mat),
            n_classes_per_var = as.integer(n_classes_per_var),
            lambda            = lambda,
            has_values        = !is.null(values_mat),
            values_have_na    = !is.null(values_mat) && anyNA(values_mat),
            user_weight_mean  = uw_mean
      )

      attr(res, "n_x") <- nrow(focal)
      attr(res, "n_pool") <- nrow(ref)
      attr(res, "n_env") <- n_env
      attr(res, "coord_type") <- coord_type
      attr(res, "binning_method") <- "fft"
      attr(res, "fft_grid") <- data.frame(
            axis    = unlist(lapply(names(fams), function(nm) {
                  if (nm == "geog") c("x", "y") else colnames(ref)[fams$env$cols] %||%
                        paste0("env", seq_len(n_env))
            }), use.names = FALSE),
            n_nodes = n_nodes,
            spacing = unlist(lapply(axes, `[[`, "delta"), use.names = FALSE),
            fft_len = pad
      )

      if (fft$check > 0L) {
            rows <- sort(sample.int(nrow(focal), min(fft$check, nrow(focal))))
            exact <- exact_fn(rows)
            err <- res[rows, , drop = FALSE] - exact
            attr(res, "fft_check") <- data.frame(
                  column      = seq_len(ncol(res)),
                  n           = colSums(is.finite(err)),
                  max_abs_err = suppressWarnings(apply(abs(err), 2, max, na.rm = TRUE)),
                  rms_err     = sqrt(colMeans(err^2, na.rm = TRUE)),
                  rms_exact   = sqrt(colMeans(exact^2, na.rm = TRUE))
            )
      }

      res
}


# Families and axes ------------------------------------------------------------

# Per-family settings. A family is active (gets grid axes) when it weights by
# distance or applies a distance limit.
.fft_family <- function(kernel, theta, max, min = 0) {
      max <- if (is.null(max)) Inf else as.numeric(max)
      list(kernel = kernel,
           theta  = if (is.null(theta) || is.na(theta)) 1 else as.numeric(theta),
           max    = max,
           min    = if (is.null(min)) 0 else as.numeric(min),
           active = kernel != "uniform" || any(is.finite(max)) ||
                 (!is.null(min) && min > 0))
}

# Grid geometry and kernel arrays for one family's axes.
#
# Node spacing is `bin` times the kernel's scale (theta, or max for a uniform
# kernel, or min if smaller), shrunk slightly so nodes span the grid's range
# exactly. Along each
# axis the grid covers only the overlap of the pool and focal ranges, each
# widened by the kernel's support: pool points beyond it reach no focal, and
# focals beyond it have no analogs (the C++ side skips / zeroes them).
#
# The kernel arrays hold each kernel type's values at integer node offsets,
# out to the largest support any used type needs:
#   "K"  distance weight inside the family's limits (and within 7 theta for a
#        Gaussian);
#   "K2" squared distance weight inside the limits, for sums of squared weights;
#   "I"  indicator of the limits, for counts and unweighted sums.
# Hard limits (max, min, per-variable boxes) are anti-aliased: a node's value
# is scaled by the fraction of its grid cell that lies inside the limits,
# which removes most of the area bias a hard edge otherwise picks up on a
# grid.
.fft_family_axes <- function(fam, pool, focal, bin, kernel_types) {
      d <- length(fam$cols)
      pc <- pool[, fam$cols, drop = FALSE]
      fc <- focal[, fam$cols, drop = FALSE]
      p_lo <- apply(pc, 2, min); p_hi <- apply(pc, 2, max)
      f_lo <- apply(fc, 2, min); f_hi <- apply(fc, 2, max)

      box <- length(fam$max) > 1L                 # per-variable limits
      lim <- rep_len(fam$max, d)                  # Inf when unlimited
      limited <- any(is.finite(lim)) || fam$min > 0

      # Support radius (distance units) per kernel type
      r_K <- if (fam$kernel == "gaussian") pmin(lim, 7 * fam$theta) else lim
      r_need <- rep(0, d)
      if (any(kernel_types %in% c("K", "K2"))) r_need <- pmax(r_need, r_K)
      if ("I" %in% kernel_types) r_need <- pmax(r_need, lim)

      # Grid range: overlap of pool and focal ranges, widened by the support
      lo <- pmax(pmin(p_lo, f_lo), pmax(p_lo, f_lo) - r_need)
      hi <- pmin(pmax(p_hi, f_hi), pmin(p_hi, f_hi) + r_need)
      hi <- pmax(hi, lo)

      # Length scale the grid must resolve: theta (max for a uniform kernel),
      # or the inner radius of an annulus if that is smaller.
      scale <- if (fam$kernel == "uniform") lim else rep(fam$theta, d)
      if (fam$min > 0) scale <- pmin(scale, fam$min)
      step <- bin * scale
      n <- ifelse(hi > lo, ceiling((hi - lo) / step) + 1, 1)
      delta <- ifelse(n > 1, (hi - lo) / (n - 1), step)
      # Half a cell extra so anti-aliased edge cells are kept
      h <- as.integer(pmin(ceiling(r_need / delta + 0.5), n - 1))
      dims <- 2L * h + 1L
      offs <- lapply(seq_len(d), function(k) seq(-h[k], h[k]) * delta[k])

      # Squared distance at each node offset, built with outer() so arrays come
      # out column-major with dims 2h + 1.
      sq_dist <- function(shift = rep(0, d)) {
            d2 <- 0
            for (k in seq_len(d)) d2 <- outer(d2, (offs[[k]] + shift[k])^2, `+`)
            array(d2, dim = dims)
      }
      dist <- sqrt(sq_dist())

      # Fraction of each node's cell inside the limits
      cover <- array(1, dim = dims)
      if (limited) {
            if (box) {
                  # Per-axis overlap of [o - delta/2, o + delta/2] with [-lim, lim]
                  cover <- 1
                  for (k in seq_len(d)) {
                        a <- offs[[k]] - delta[k] / 2
                        b <- offs[[k]] + delta[k] / 2
                        frac <- pmax(0, pmin(b, lim[k]) - pmax(a, -lim[k])) / delta[k]
                        cover <- outer(cover, frac, `*`)
                  }
                  cover <- array(cover, dim = dims)
                  if (fam$min > 0) cover <- cover * (dist >= fam$min)
            } else {
                  # Supersample each cell on an s^d lattice
                  s <- if (d <= 2L) 8L else 4L
                  sub <- (seq_len(s) - 0.5) / s - 0.5
                  shifts <- as.matrix(expand.grid(rep(list(sub), d)))
                  cover <- array(0, dim = dims)
                  for (q in seq_len(nrow(shifts))) {
                        r <- sqrt(sq_dist(shifts[q, ] * delta))
                        cover <- cover + (r <= fam$max & r >= fam$min)
                  }
                  cover <- cover / nrow(shifts)
            }
      }

      kern <- list()
      if (any(kernel_types %in% c("K", "K2"))) {
            w <- .kernel_weight_fn(fam$kernel, fam$theta)(dist)
            if (fam$kernel == "gaussian") w[dist > 7 * fam$theta] <- 0
            kern$K <- array(w * cover, dim = dims)
            kern$K2 <- array(w^2 * cover, dim = dims)
      }
      if ("I" %in% kernel_types) kern$I <- cover

      list(origin = lo, delta = delta, n = as.integer(n), h = h,
           pad = as.integer(vapply(n + h, stats::nextn, numeric(1), factors = c(2, 3, 5))),
           kernels = kern)
}


# Fields -----------------------------------------------------------------------

# Decide which kernel-weighted sums ("fields") the requested stats need, build
# their per-pool-point masses, and map each accumulator slot of the C++
# finalizer to its field (0-based; -1 = absent). Each field is a kernel type
# ("K", "K2", "I") and a product of named per-point factors:
#   pw    point weight (cell area x user weight)
#   y<v>  response v (missing as 0)      o<v>  1 if response v is observed
#   x<p>  covariate p                    c<v>_<k>  1 if response v is class k
# Identical (kernel, factors) pairs share one field, so e.g. sum_w serves
# sum_weights, weighted_mean's denominator, and the regression's X'WX[0, 0].
.fft_build_fields <- function(stat, se, values_mat, covariates_mat,
                              n_classes_per_var, area_weight, user_weight,
                              n_pool) {
      pw <- rep(1, n_pool)
      if (!is.null(area_weight)) pw <- pw * area_weight
      if (!is.null(user_weight)) pw <- pw * user_weight
      n_vars <- if (is.null(values_mat)) 0L else ncol(values_mat)
      n_covs <- if (is.null(covariates_mat)) 0L else ncol(covariates_mat)

      factor_vec <- function(tok) {
            type <- substr(tok, 1, 1)
            if (tok == "pw") return(pw)
            if (type == "y") {
                  y <- values_mat[, as.integer(substring(tok, 2))]
                  y[is.na(y)] <- 0
                  return(y)
            }
            if (type == "o") return(as.numeric(!is.na(values_mat[, as.integer(substring(tok, 2))])))
            if (type == "x") return(covariates_mat[, as.integer(substring(tok, 2))])
            if (type == "c") {
                  vk <- as.integer(strsplit(substring(tok, 2), "_")[[1]])
                  y <- values_mat[, vk[1]]
                  return(as.numeric(!is.na(y) & y == vk[2]))
            }
            stop("Internal error: unknown field factor '", tok, "'.", call. = FALSE)
      }

      keys <- character(0)
      masses <- list()
      kern_of <- character(0)
      add <- function(kern, toks) {
            toks <- sort(toks)
            key <- paste0(kern, ":", paste(toks, collapse = "*"))
            hit <- match(key, keys)
            if (!is.na(hit)) return(hit - 1L)
            keys <<- c(keys, key)
            kern_of <<- c(kern_of, kern)
            masses[[length(masses) + 1L]] <<- if (length(toks) == 0L) rep(1, n_pool) else
                  Reduce(`*`, lapply(toks, factor_vec))
            length(keys) - 1L
      }
      yv <- function(v) paste0("y", v)
      xt <- function(a) if (a == 0L) character(0) else paste0("x", a)  # a = 0: intercept

      weighted_stats <- c("sum_weights", "mean_weights", "ess", "weighted_sum",
                          "weighted_mean", "regression", "tabulate")
      value_stats <- stat[stat %in% c("sum", "mean", "weighted_sum", "weighted_mean")]
      want_se <- se != "none"
      wm_se <- want_se && "weighted_mean" %in% stat && n_vars > 0L
      reg <- "regression" %in% stat && n_vars > 0L

      roles <- list(
            count = if (any(stat %in% c("count", "mean_weights", "mean"))) add("I", character(0)) else -1L,
            sw    = if (any(stat %in% weighted_stats)) add("K", "pw") else -1L,
            sw2   = if ("ess" %in% stat || (want_se && any(stat %in% c("weighted_mean", "regression"))))
                  add("K2", c("pw", "pw")) else -1L
      )

      roles$value <- unlist(lapply(seq_len(n_vars), function(v) {
            vapply(value_stats, function(st) {
                  add(if (st %in% c("sum", "mean")) "I" else "K", c("pw", yv(v)))
            }, integer(1))
      }), use.names = FALSE)
      roles$wm_wy2 <- if (wm_se) vapply(seq_len(n_vars), function(v) add("K", c("pw", yv(v), yv(v))), integer(1))
      roles$wm_w2y <- roles$wm_w2y2 <- NULL
      if (wm_se && se == "design") {
            roles$wm_w2y  <- vapply(seq_len(n_vars), function(v) add("K2", c("pw", "pw", yv(v))), integer(1))
            roles$wm_w2y2 <- vapply(seq_len(n_vars), function(v) add("K2", c("pw", "pw", yv(v), yv(v))), integer(1))
      }
      # Per-variable denominators when any response is missing (see
      # focal_stats.h): fields over the analogs observing each variable, or the
      # shared fields for variables with no missing values.
      if (n_vars > 0L && anyNA(values_mat) && length(value_stats) > 0L) {
            has_na <- colSums(is.na(values_mat)) > 0
            per_var <- function(shared, kern, toks) {
                  if (shared < 0L) return(integer(0))
                  vapply(seq_len(n_vars), function(v) {
                        if (has_na[v]) add(kern, c(toks, paste0("o", v))) else shared
                  }, integer(1))
            }
            roles$obs_count <- per_var(roles$count, "I", character(0))
            roles$obs_sw    <- per_var(roles$sw, "K", "pw")
            roles$obs_sw2   <- per_var(roles$sw2, "K2", c("pw", "pw"))
      }
      roles$tab <- if ("tabulate" %in% stat) unlist(lapply(seq_len(n_vars), function(v) {
            vapply(seq_len(n_classes_per_var[v]), function(k) add("K", c("pw", paste0("c", v, "_", k))), integer(1))
      }), use.names = FALSE)

      if (reg) {
            dim <- n_covs + 1L
            grid_ac <- expand.grid(c = seq_len(dim) - 1L, a = seq_len(dim) - 1L)  # row-major (a, c)
            roles$xtwx <- mapply(function(a, c) add("K", c("pw", xt(a), xt(c))), grid_ac$a, grid_ac$c)
            roles$xtwy <- unlist(lapply(seq_len(n_vars), function(v) {
                  vapply(seq_len(dim) - 1L, function(a) add("K", c("pw", yv(v), xt(a))), integer(1))
            }), use.names = FALSE)
            if (se == "ess") {
                  roles$reg_wy2 <- vapply(seq_len(n_vars), function(v) add("K", c("pw", yv(v), yv(v))), integer(1))
            }
            # Own systems for responses with missing values (see ridge_solve.h)
            own_vars <- which(colSums(is.na(values_mat)) > 0)
            roles$own <- rep(-1L, n_vars)
            roles$own[own_vars] <- seq_along(own_vars) - 1L
            roles$own_xtwx <- unlist(lapply(own_vars, function(v) {
                  mapply(function(a, c) add("K", c("pw", paste0("o", v), xt(a), xt(c))), grid_ac$a, grid_ac$c)
            }), use.names = FALSE)
            roles$own_sw <- vapply(own_vars, function(v) add("K", c("pw", paste0("o", v))), integer(1))
            if (se == "ess") {
                  roles$own_sw2 <- vapply(own_vars, function(v) add("K2", c("pw", "pw", paste0("o", v))), integer(1))
            }
      }

      slots <- c("count", "sw", "sw2", "obs_count", "obs_sw", "obs_sw2",
                 "value", "wm_wy2", "wm_w2y", "wm_w2y2", "tab",
                 "xtwx", "xtwy", "reg_wy2", "own", "own_xtwx", "own_sw", "own_sw2")
      roles <- lapply(stats::setNames(slots, slots), function(s) as.integer(roles[[s]] %||% integer(0)))

      kernel_types <- unique(kern_of)
      list(masses       = do.call(cbind, masses) %||% matrix(0, n_pool, 0),
           field_kernel = as.integer(match(kern_of, kernel_types) - 1L),
           kernel_types = kernel_types,
           roles        = roles)
}
