#' Find Climate Analogs
#'
#' Identifies locations in a reference dataset that are climatically similar to
#' focal locations, with optional constraints on climate distance and geographic
#' distance. This function supports multiple use cases including climate velocity
#' analysis, analog availability mapping, and climate impact assessment.
#'
#' The function uses a spatial indexing structure (lattice-based) to quickly
#' search through large reference datasets. Climate similarity is measured
#' using Euclidean distance in climate space (ideally pre-whitened; see Details).
#' Geographic distance can be computed for lon/lat coordinates (great-circle
#' distance) or projected coordinates (planar distance).
#'
#' @param focal Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers. Each row represents a focal
#'   location for which analogs will be found.
#'
#' @param ref Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers. The reference dataset to search
#'   for analogs.
#'
#' @param mode Character string specifying the analog search mode. One of:
#'   \itemize{
#'     \item \code{"knn_clim"}: For each focal, return up to \code{k} analogs
#'       with smallest climate distance, subject to \code{max_clim} and
#'       \code{max_geog} filters.
#'     \item \code{"knn_geog"}: For each focal, return up to \code{k} analogs
#'       with smallest geographic distance, subject to \code{max_clim} and
#'       \code{max_geog} filters.
#'     \item \code{"all"}: Return all analogs that satisfy the filters.
#'     \item \code{"count"}: For each focal, count how many analogs satisfy
#'       the filters.
#'     \item \code{"sum"}: For each focal, sum weights of all analogs that
#'       satisfy the filters (see \code{weight} and \code{theta}).
#'     \item \code{"mean"}: For each focal, mean of weights of all analogs that
#'       satisfy the filters.
#'   }
#'
#' @param max_clim Maximum climate distance constraint (default: NULL = no
#'   climate constraint). Can be either:
#'   \itemize{
#'     \item A scalar: Euclidean radius in climate space (e.g., 0.5)
#'     \item A vector: Per-variable absolute differences (length must equal
#'       number of climate variables)
#'   }
#'   Only reference locations within this climate distance are considered.
#'
#' @param max_geog Maximum geographic distance constraint in km (default:
#'   NULL = no geographic constraint). When specified, only reference locations
#'   within this distance are considered.
#'
#' @param k Number of nearest analogs to return per focal location for kNN
#'   modes. Required when \code{mode} is \code{"knn_geog"} or \code{"knn_clim"};
#'   must be \code{NULL} for other modes.
#'
#' @param weight Weighting function for matches, used only when
#'   \code{mode} is \code{"sum"} or \code{"mean"}. One of:
#'   \itemize{
#'     \item \code{"uniform"}: All matches weighted equally (weight = 1.0).
#'     \item \code{"inverse_clim"}: Weight = 1 / (climate_distance + epsilon),
#'       with epsilon given by \code{theta} (or a small default if \code{theta}
#'       is \code{NULL}).
#'     \item \code{"inverse_geog"}: Weight = 1 / (geographic_distance + epsilon),
#'       with epsilon given by \code{theta} (or a small default if \code{theta}
#'       is \code{NULL}).
#'   }
#'   For \code{mode} in \code{"knn_geog"}, \code{"knn_clim"}, \code{"count"},
#'   or \code{"all"}, \code{weight} must be \code{NULL}.
#'
#' @param theta Optional numeric parameter used by some weighting kernels
#'   when \code{mode} is \code{"sum"} or \code{"mean"} and \code{weight} is
#'   not \code{"uniform"}. Currently interpreted as:
#'   \itemize{
#'     \item For \code{"gaussian_clim"}: sigma bandwidth for climate distance.
#'     \item For \code{"gaussian_geog"}: sigma bandwidth for geographic distance.
#'     \item For \code{"gaussian_joint"}: length-two vector of sigma bandwidths
#'       for climate and geographic distances, respectively.
#'     \item For \code{"inverse_clim"}: epsilon added to climate distance.
#'     \item For \code{"inverse_geog"}: epsilon added to geographic distance.
#'   }
#'   If \code{theta} is \code{NULL}, a small default epsilon is used. For
#'   \code{weight = "uniform"} or for non-aggregating modes, \code{theta}
#'   must be \code{NULL}.
#'
#' @param report_dist Logical; if TRUE (default), include distance columns in
#'   output when \code{mode} is \code{"knn_geog"}, \code{"knn_clim"} or
#'   \code{"all"}. Set to FALSE for more compact output.
#'
#' @param coord_type Coordinate system type (default: "auto"):
#'   \itemize{
#'     \item \code{"auto"}: Automatically detect from coordinate ranges.
#'     \item \code{"lonlat"}: Unprojected lon/lat coordinates (uses great-circle distance).
#'     \item \code{"projected"}: Projected XY coordinates (uses planar distance).
#'   }
#'
#' @param lattice_res Tuning parameter giving the number of bins per dimension
#'   of the internally-used lattice search index. Either:
#'   \itemize{
#'     \item A positive integer.
#'     \item \code{"auto"} (the default): Automatically tune the lattice resolution
#'       by optimizing compute time on a subsample of focal points. If focal has
#'       relatively few rows, auto-tuning is skipped and a default resolution of
#'       16 is used.
#'   }
#'   If "auto" and if there are more than 2000
#'
#' @param n_threads Optional integer number of threads to use for the
#'   computation. If \code{NULL} (default), the global RcppParallel setting
#'   is used (see \code{RcppParallel::setThreadOptions}).
#'
#' @return
#' The return value depends on the \code{mode} parameter:
#'
#' **For mode = "knn_geog", "knn_clim" or "all"**:
#' A data.frame with one row per focal-analog pair:
#' \itemize{
#'   \item \code{focal_index}: Index of focal location (1-based).
#'   \item \code{focal_x, focal_y}: Coordinates of focal location.
#'   \item \code{analog_index}: Index of analog location in reference dataset (1-based).
#'   \item \code{analog_x, analog_y}: Coordinates of analog location.
#'   \item \code{clim_dist}: Climate distance (if \code{report_dist = TRUE}).
#'   \item \code{geog_dist}: Geographic distance in km (if \code{report_dist = TRUE}).
#' }
#'
#' **For mode = "sum", "mean", or "count"**:
#' A data.frame with one row per focal location:
#' \itemize{
#'   \item \code{focal_index}: Index of focal location (1-based).
#'   \item \code{focal_x, focal_y}: Coordinates of focal location.
#'   \item \code{value}: Aggregated value (count, sum of weights, or mean of weights).
#' }
#'
#' All outputs include diagnostic attributes propagated from the C++ core,
#' including:
#' \itemize{
#'   \item \code{total_bins}: Number of spatial bins in the lattice index.
#'   \item \code{avg_bin_occupancy}: Average points per bin.
#'   \item \code{min_bin_occupancy, max_bin_occupancy}: Range of bin occupancy.
#'   \item \code{binning_method}: Method used ("multi_dim_lattice" or "none").
#'   \item \code{n_ref, n_clim}: Size of reference dataset and number of climate variables.
#' }
#'
#' @details
#' **Common Use Cases:**
#'
#' \strong{Climate Velocity} (nearest geographic neighbor with similar climate):
#' \preformatted{
#' find_analogs(
#'   focal   = clim$clim1,
#'   ref     = clim$clim2,
#'   mode    = "knn_geog",
#'   max_clim = 0.5,
#'   max_geog = NULL,
#'   k        = 1
#' )
#' }
#'
#' \strong{Climate Impact} (climatically similar locations within dispersal range):
#' \preformatted{
#' find_analogs(
#'   focal   = clim$clim1,
#'   ref     = clim$clim2,
#'   mode    = "knn_clim",
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   k        = 20
#' )
#' }
#'
#' \strong{Analog Availability} (count of suitable locations):
#' \preformatted{
#' find_analogs(
#'   focal   = clim$clim1,
#'   ref     = clim$clim1,
#'   mode    = "count",
#'   max_clim = 0.5,
#'   max_geog = 100
#' )
#' }
#'
#' \strong{Weighted Analog Intensity} (e.g., distance-weighted availability):
#' \preformatted{
#' find_analogs(
#'   focal   = clim$clim1,
#'   ref     = clim$clim1,
#'   mode    = "sum",
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   weight   = "inverse_geog",
#'   theta    = 1e-6
#' )
#' }
#'
#' @export
find_analogs <- function(
            focal,
            ref,
            mode = c("knn_clim", "knn_geog", "count", "sum", "mean", "all"),
            max_clim = NULL,
            max_geog = NULL,
            k = NULL,
            weight = c("uniform",
                       "gaussian_clim", "gaussian_geog", "gaussian_joint",
                       "inverse_clim", "inverse_geog", "inverse_joint"),
            theta = NULL,
            report_dist = TRUE,
            coord_type = c("auto", "lonlat", "projected"),
            lattice_res = "auto",
            n_threads = NULL
) {

      # ---- Input validation --------------------------------------------------
      coord_type <- match.arg(coord_type)
      mode <- match.arg(mode)

      weight <- match.arg(weight)
      if(!weight %in% c("inverse_clim", "inverse_geog")) weight <- NULL

      # Validate combination of mode, k, weight, theta
      if (mode %in% c("knn_clim", "knn_geog")) {
            # kNN modes: require k, disallow weight/theta
            if (is.null(k)) {
                  k <- 1L
            }
            k <- as.integer(k)
            if (length(k) != 1L || k <= 0L) {
                  stop("For mode '", mode, "', k must be a positive integer.")
            }
            if (!is.null(weight)) {
                  stop("For mode '", mode, "', weight must be NULL.")
            }
            if (!is.null(theta)) {
                  stop("For mode '", mode, "', theta must be NULL.")
            }
      } else {
            # Non-kNN modes: k must be NULL
            if (!is.null(k)) {
                  stop("For mode '", mode, "', k must be NULL.")
            }
            k <- 0L
      }

      if (mode %in% c("all", "count")) {
            # No weighting allowed
            if (!is.null(weight)) {
                  stop("For mode '", mode, "', weight must be NULL.")
            }
            if (!is.null(theta)) {
                  stop("For mode '", mode, "', theta must be NULL.")
            }
      }

      if (mode %in% c("sum", "mean")) {
            # Aggregation modes: weight is required, theta optional
            valid_weights <- c("uniform", "inverse_clim", "inverse_geog")
            if (is.null(weight)) {
                  weight <- "uniform"
            }
            if (!weight %in% valid_weights) {
                  stop(
                        "For mode '",
                        mode,
                        "', weight must be one of: ",
                        paste(valid_weights, collapse = ", ")
                  )
            }
            if (identical(weight, "uniform")) {
                  if (!is.null(theta)) {
                        stop("For weight = 'uniform', theta must be NULL.")
                  }
            } else {
                  # inverse_*: theta is epsilon; if NULL, we'll use a default in aggregators
                  if (!is.null(theta)) {
                        if (
                              !is.numeric(theta) ||
                              length(theta) != 1L ||
                              theta <= 0
                        ) {
                              stop(
                                    "theta must be a single positive numeric value, or NULL."
                              )
                        }
                  }
            }
      } else {
            # For non-aggregation modes, weight/theta must be NULL
            if (!is.null(weight)) {
                  stop("weight must be NULL when mode is not 'sum' or 'mean'.")
            }
            if (!is.null(theta)) {
                  stop("theta must be NULL when mode is not 'sum' or 'mean'.")
            }
      }


      # ---- Data normalization ------------------------------------------------

      focal_mm <- .format_data(focal)
      ref_mm <- .format_data(ref)

      # Detect geographic coordinate system
      geo_mode <- switch(
            coord_type,
            auto = .detect_geo(focal_mm[, 1:2, drop = FALSE],
                               ref_mm[, 1:2, drop = FALSE]),
            lonlat = "lonlat",
            projected = "projected"
      )

      # Parse constraints
      max_geog_num <- if (is.null(max_geog)) Inf else as.numeric(max_geog)[1L]

      max_clim_val <- if (is.null(max_clim)) {
            Inf
      } else {
            max_clim
      }

      # ---- Map mode/weight/theta for C++ -------------------------------------
      mode_code <- switch(
            mode,
            "knn_clim" = 0L,
            "knn_geog" = 1L,
            "count"    = 2L,
            "sum"      = 3L,
            "mean"     = 4L,
            "all"      = 5L
      )

      weight_code <- if (mode %in% c("sum","mean")) {
            switch(
                  weight,
                  "uniform"      = 1L,
                  "inverse_clim" = 2L,
                  "inverse_geog" = 3L
            )
      } else {
            0L
      }

      theta_num <- if (is.null(theta)) NA_real_ else as.numeric(theta)[1L]

      # For kNN modes, k_core = k; for others, k_core = 0 to request "all matches"
      k_core <- if (mode %in% c("knn_clim","knn_geog")) as.integer(k) else 0L


      # ---- Call C++ core ------------------------------------------------------

      # Optional per-call thread control
      if (!is.null(n_threads)) {
            if (!requireNamespace("RcppParallel", quietly = TRUE)) {
                  stop("Package 'RcppParallel' is required to control 'n_threads'. ",
                       "Install it with install.packages('RcppParallel').",
                       call. = FALSE)
            }
            RcppParallel::setThreadOptions(numThreads = as.integer(n_threads)[1L])
      }

      # Internal wrapper so we can reuse it for tuning and the full call
      call_analogs_core <- function(focal_mm_sub, lattice_res_int) {
            .Call(
                  `_analogs_find_analogs_core`,
                  focal_mm_sub,               # matrix of focal sites (xy + climate)
                  ref_mm,                     # matrix of ref sites  (xy + climate)
                  as.integer(k_core),         # k for kNN, 0 for all/aggregates
                  max_clim_val,               # climate filter bandwidth (scalar, vector, or Inf)
                  as.numeric(max_geog_num),   # geographic distance filter (km; Inf if NULL)
                  geo_mode,                   # "lonlat" or "projected"
                  as.integer(mode_code),      # mode
                  as.integer(weight_code),    # weight
                  as.numeric(theta_num),      # theta
                  as.integer(lattice_res_int) # lattice resolution (integer)
            )
      }

      lattice_res_int <- .tune_lattice_res(
            focal_mm       = focal_mm,
            lattice_res    = lattice_res,
            call_analogs_core = call_analogs_core,
            k_core         = k_core,
            default_res    = 16L
      )

      res <- call_analogs_core(focal_mm, lattice_res_int)

      # Capture diagnostic attributes from C++ before post-processing
      cpp_attrs <- attributes(res)
      cpp_attrs$names <- NULL
      cpp_attrs$class <- NULL


      # ---- Post-process results ----------------------------------------------

      if (mode %in% c("knn_clim", "knn_geog", "all")) {
            out <- .emit_pairs_cpp(
                  res,
                  focal_mm,
                  ref_mm,
                  report_dist = report_dist,
                  geo_mode    = geo_mode
            )
            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            attr(out, "mode")   <- mode
            attr(out, "weight") <- weight
            attr(out, "theta")  <- theta
            return(out)
      }

      if (mode %in% c("sum", "mean", "count")) {
            # C++ now returns a NumericVector of length n_focal with the
            # aggregated value (count, sum of weights, or mean of weights)
            values <- as.numeric(res)
            if (length(values) != nrow(focal_mm)) {
                  stop("Internal error: aggregate result length does not match number of focals.")
            }

            out <- data.frame(
                  focal_index = seq_len(nrow(focal_mm)),
                  focal_x     = focal_mm[, 1],
                  focal_y     = focal_mm[, 2],
                  value       = values,
                  stringsAsFactors = FALSE
            )

            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            attr(out, "mode")   <- mode
            attr(out, "weight") <- weight
            attr(out, "theta")  <- theta
            return(out)
      }

      stop("Unreachable code - please report this bug")
}



# ---- Internal Helper Functions ---------------------------------------------

#' Detect coordinate system from data ranges
#' @keywords internal
.detect_geo <- function(focal_xy, ref_xy) {
      lon_rng <- range(c(focal_xy[, 1], ref_xy[, 1]), na.rm = TRUE)
      lat_rng <- range(c(focal_xy[, 2], ref_xy[, 2]), na.rm = TRUE)

      if (
            all(is.finite(c(lon_rng, lat_rng))) &&
            lon_rng[1] >= -180 &&
            lon_rng[2] <= 180 &&
            lat_rng[1] >= -90 &&
            lat_rng[2] <= 90
      ) {
            "lonlat"
      } else {
            "projected"
      }
}

#' Compute great-circle distance using Haversine formula
#'
#' NOT CURRENTLY USED
#'
#' @keywords internal
.haversine_km <- function(xy1, xy2) {
      R <- 6371.0088 # Earth's mean radius in km
      to_rad <- pi / 180

      lon1 <- xy1[, 1] * to_rad
      lat1 <- xy1[, 2] * to_rad
      lon2 <- xy2[, 1] * to_rad
      lat2 <- xy2[, 2] * to_rad

      dlon <- lon2 - lon1
      dlat <- lat2 - lat1

      sdlat <- sin(0.5 * dlat)
      sdlon <- sin(0.5 * dlon)

      a <- sdlat * sdlat + cos(lat1) * cos(lat2) * sdlon * sdlon

      2 * R * asin(pmin(1, sqrt(a)))
}


#' Extract coordinates and climate data from input
#' @keywords internal
.select_xy_climate <- function(obj) {
      nm <- colnames(obj)

      # Try to find x,y columns by name
      if (!is.null(nm) && all(c("x", "y") %in% nm)) {
            xy_idx <- match(c("x", "y"), nm)
      } else {
            # Default to first two columns
            xy_idx <- 1:2
      }

      coords <- as.matrix(obj[, xy_idx, drop = FALSE])
      climate <- as.matrix(obj[,
                               setdiff(seq_len(ncol(obj)), xy_idx),
                               drop = FALSE
      ])

      storage.mode(coords) <- "double"
      storage.mode(climate) <- "double"

      if (ncol(coords) != 2L) {
            stop("Coordinate data must have exactly 2 columns (x, y)")
      }
      if (ncol(climate) < 1L) {
            stop(
                  "No climate variable columns found after extracting coordinates"
            )
      }

      cbind(coords, climate)
}

#' Normalize input to standard format
#' @keywords internal
.format_data <- function(r) {
      if (inherits(r, "SpatRaster")) {
            # Convert SpatRaster to data.frame
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs")
            }
            df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
            .select_xy_climate(df)
      } else if (is.matrix(r) || is.data.frame(r)) {
            .select_xy_climate(r)
      } else {
            stop("Input must be a data.frame, matrix, or SpatRaster")
      }
}

# Internal: auto-tune lattice_res, or resolve user/default value
.tune_lattice_res <- function(focal_mm,
                             lattice_res,
                             call_analogs_core,
                             k_core,
                             default_res = 16L,
                             verbose = FALSE) {

      # Helper: detect monotonic timings with a tolerance
      is_strict_monotonic <- function(x, tol = 0.15) {
            if (length(x) < 3L) return(FALSE)
            dx <- diff(x)
            # strictly decreasing or increasing, with relative tolerance
            all(dx <= -tol * abs(x[-length(x)])) ||
                  all(dx >=  tol * abs(x[-length(x)]))
      }

      # If user supplied a numeric value, just honor it
      if (is.numeric(lattice_res)) {
            return(as.integer(lattice_res))
      }

      # NULL: simple fallback to default
      if (is.null(lattice_res)) {
            return(as.integer(default_res))
      }

      # Only special-case we support is "auto"
      if (!(is.character(lattice_res) && identical(lattice_res, "auto"))) {
            stop("Unsupported value for lattice_res: ", lattice_res)
      }

      ## --- "auto" path below ---

      n_focal <- nrow(focal_mm)

      # Only tune for non-trivial problem sizes
      if (n_focal <= 2000L) {
            return(as.integer(default_res))
      }

      # Subsample focal sites
      n_samp <- min(1000L, max(100L, as.integer(n_focal * 0.01)))
      idx    <- sample.int(n_focal, n_samp)
      focal_mm_samp <- focal_mm[idx, , drop = FALSE]

      # Helper to evaluate timing for a given lattice_res
      eval_time <- function(r) {
            r <- as.integer(max(1L, r))  # ensure valid
            st <- system.time(invisible(call_analogs_core(focal_mm_samp, r)))
            st[["elapsed"]]
      }

      # Start from heuristic center
      r0 <- as.integer(default_res)
      r_vals <- c(r0 %/% 2L, r0, r0 * 2L)
      r_vals <- unique(r_vals[r_vals > 0L])

      # Evaluate initial bracket
      times <- vapply(r_vals, eval_time, numeric(1))

      # One-step adaptive refinement
      best_idx <- which.min(times)
      best_r   <- r_vals[best_idx]

      if (best_idx == 1L && length(r_vals) > 1L) {
            # best is smallest → try halving again
            r_try <- max(1L, best_r %/% 2L)
            t_try <- eval_time(r_try)
            r_vals <- c(r_try, r_vals)
            times  <- c(t_try, times)

      } else if (best_idx == length(r_vals) && length(r_vals) > 1L) {
            # best is largest → try doubling again
            r_try <- best_r * 2L
            t_try <- eval_time(r_try)
            r_vals <- c(r_vals, r_try)
            times  <- c(times, t_try)
      }

      # Recompute best after any refinement
      best_idx <- which.min(times)
      best_r   <- r_vals[best_idx]

      # Sort evaluated points by r for monotonicity check
      o      <- order(r_vals)
      r_eval <- r_vals[o]
      t_eval <- times[o]

      # Optional: warn if timings are monotonic in a meaningful regime
      if (length(t_eval) >= 4L &&
          n_samp >= 300L &&
          !(is.numeric(k_core) && k_core == 1L) &&
          is_strict_monotonic(t_eval)) {

            warning(
                  "Auto-tuning of lattice_res did not detect an interior minimum; ",
                  "elapsed times were monotonic across tested values {",
                  paste(r_eval, collapse = ", "),
                  "}. The optimal lattice_res may lie outside this range. ",
                  "Consider manually specifying lattice_res."
            )
      }

      r <- as.integer(best_r)
      if(verbose) message("Auto-tuning selected lattice resolution of ", r, ".")
      r
}

