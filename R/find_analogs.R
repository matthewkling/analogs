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
      weight = NULL,
      theta = NULL,
      report_dist = TRUE,
      coord_type = c("auto", "lonlat", "projected")
) {
      # ---- Input validation --------------------------------------------------
      coord_type <- match.arg(coord_type)
      mode <- match.arg(mode)

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
            auto = .detect_geo(focal_mm[, 1:2], ref_mm[, 1:2]),
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
      res <- .Call(
            `_analogs_find_analogs_core`,
            focal_mm,                 # matrix of focal sites (xy + climate)
            ref_mm,                   # matrix of ref sites  (xy + climate)
            as.integer(k_core),       # k for kNN, 0 for all/aggregates
            max_clim_val,             # climate filter bandwidth (scalar or vector or Inf)
            as.numeric(max_geog_num), # geographic distance filter (km; Inf if NULL)
            geo_mode,                 # "lonlat" or "projected"
            as.integer(mode_code),    # new
            as.integer(weight_code),  # new
            as.numeric(theta_num)     # new
      )

      # Capture diagnostic attributes from C++ before post-processing
      cpp_attrs <- attributes(res)
      cpp_attrs$names <- NULL
      cpp_attrs$class <- NULL

      # ---- Post-process results ----------------------------------------------
      if (mode %in% c("knn_clim", "knn_geog", "all")) {
            out <- .emit_pairs(
                  res,
                  focal_mm,
                  ref_mm,
                  report_dist = report_dist,
                  geo_mode = geo_mode
            )
            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            attr(out, "mode") <- mode
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

#' Build long table of focal-analog pairs with distances
#' @keywords internal
.emit_pairs <- function(res, focal_mm, ref_mm, report_dist, geo_mode) {
      n_f <- nrow(focal_mm[, 1:2])
      rows <- vector("list", n_f)

      for (i in seq_len(n_f)) {
            idx <- res[[i]]
            if (length(idx) == 0L) {
                  rows[[i]] <- NULL
                  next
            }

            # Focal coordinates (repeated)
            fx <- rep(focal_mm[, 1:2][i, 1], length(idx))
            fy <- rep(focal_mm[, 1:2][i, 2], length(idx))

            # Analog coordinates
            ax <- ref_mm[, 1:2][idx, 1]
            ay <- ref_mm[, 1:2][idx, 2]

            df <- data.frame(
                  focal_index = rep.int(i, length(idx)),
                  focal_x = fx,
                  focal_y = fy,
                  analog_index = idx,
                  analog_x = ax,
                  analog_y = ay,
                  stringsAsFactors = FALSE
            )

            if (isTRUE(report_dist)) {
                  # Climate distances (Euclidean in climate space)
                  v_i <- matrix(focal_mm[, 3:ncol(focal_mm)][i, ], nrow = 1)
                  v_j <- ref_mm[, 3:ncol(ref_mm)][idx, , drop = FALSE]
                  clim_d <- sqrt(rowSums((t(t(v_j) - as.numeric(v_i)))^2))

                  # Geographic distances
                  if (geo_mode == "lonlat") {
                        geog_d <- .haversine_km(cbind(fx, fy), cbind(ax, ay))
                  } else {
                        geog_d <- sqrt((ax - fx)^2 + (ay - fy)^2)
                  }

                  df$clim_dist <- clim_d
                  df$geog_dist <- geog_d
            }

            rows[[i]] <- df
      }

      # Handle case where all focal points have no matches
      if (all(vapply(rows, is.null, logical(1)))) {
            return(data.frame(
                  focal_index = integer(0),
                  focal_x = numeric(0),
                  focal_y = numeric(0),
                  analog_index = integer(0),
                  analog_x = numeric(0),
                  analog_y = numeric(0)
            ))
      }

      do.call(rbind, rows)
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
