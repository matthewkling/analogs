#' Find Climate Analogs
#'
#' Identifies locations in a reference dataset that are climaticallyc similar to
#' focal locations, with optional constraints on climate distance and geographic
#' distance. This function supports multiple use cases including climate velocity
#' analysis, analog availability mapping, and climate impact assessment.
#'
#' The function uses an efficient spatial indexing structure (quantile-binned lattice)
#' to quickly search through large reference datasets. Climate similarity is measured
#' using Euclidean distance in climate space (ideally pre-whitened; see Details).
#' Geographic distance can be computed for lat/lon coordinates (great-circle distance)
#' or projected coordinates (planar distance).
#'
#' @param focal Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers. Each row represents a focal
#'   location for which analogs will be found.
#'
#' @param ref Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers. The reference dataset to search
#'   for analogs.
#'
#' @param max_dist Maximum geographic distance constraint in km (default: NULL =
#'   no constraint). When specified, only reference locations within this distance
#'   are considered. Useful for defining "dispersal-constrained" analogs.
#'
#' @param max_clim Maximum climate distance constraint (default: NULL = no constraint).
#'   Can be either:
#'   \itemize{
#'     \item A scalar: Euclidean radius in climate space (e.g., 0.5)
#'     \item A vector: Per-variable absolute differences (length must equal number of climate variables)
#'   }
#'   Only reference locations within this climate distance are considered.
#'
#' @param weight Weighting function for matches (default: "inverse_clim"):
#'   \itemize{
#'     \item \code{"uniform"}: All matches weighted equally (weight = 1.0)
#'     \item \code{"inverse_clim"}: Weight = 1 / (climate_distance + small_constant)
#'     \item \code{"inverse_dist"}: Weight = 1 / (geographic_distance + small_constant)
#'   }
#'   Used for aggregation functions and optionally reported with matches.
#'
#' @param fun Aggregation function determining what to return (default: "nmax"):
#'   \itemize{
#'     \item \code{"nmax"}: Return top n closest matches (requires \code{n} parameter)
#'     \item \code{"all"}: Return all matches passing constraints
#'     \item \code{"count"}: Count number of matches per focal location
#'     \item \code{"sum"}: Sum of weights across all matches
#'     \item \code{"mean"}: Mean of weights across all matches
#'   }
#'
#' @param n Number of analogs to return per focal location. Required when
#'   \code{fun = "nmax"}, ignored otherwise.
#'
#' @param metric Distance metric for climate space (default: "euclidean"). Currently
#'   only Euclidean distance is supported. Future versions may support additional metrics.
#'
#' @param report_dist Logical; if TRUE (default), include distance and weight columns
#'   in output when \code{fun} is "nmax" or "all". Set to FALSE for more compact output.
#'
#' @param coord_type Coordinate system type (default: "auto"):
#'   \itemize{
#'     \item \code{"auto"}: Automatically detect from coordinate ranges
#'     \item \code{"lonlat"}: Unprojected lon/lat coordinates (uses great-circle distance)
#'     \item \code{"projected"}: Projected XY coordinates (uses planar distance)
#'   }
#'
#' @return
#' The return value depends on the \code{fun} parameter:
#'
#' **For fun = "nmax" or "all"**: A data.frame with one row per focal-analog pair:
#' \itemize{
#'   \item \code{focal_index}: Index of focal location (1-based)
#'   \item \code{focal_x, focal_y}: Coordinates of focal location
#'   \item \code{analog_index}: Index of analog location in reference dataset (1-based)
#'   \item \code{analog_x, analog_y}: Coordinates of analog location
#'   \item \code{clim_dist}: Climate distance (if \code{report_dist = TRUE})
#'   \item \code{geog_dist}: Geographic distance in km (if \code{report_dist = TRUE})
#'   \item \code{weight}: Computed weight value (if \code{report_dist = TRUE})
#' }
#'
#' **For fun = "sum", "mean", or "count"**: A data.frame with one row per focal location:
#' \itemize{
#'   \item \code{focal_index}: Index of focal location (1-based)
#'   \item \code{focal_x, focal_y}: Coordinates of focal location
#'   \item \code{value}: Aggregated value (count, sum, or mean of weights)
#' }
#'
#' All outputs include diagnostic attributes:
#' \itemize{
#'   \item \code{total_bins}: Number of spatial bins created in index
#'   \item \code{avg_bin_occupancy}: Average points per bin
#'   \item \code{min_bin_occupancy, max_bin_occupancy}: Range of bin occupancy
#'   \item \code{binning_method}: Method used ("quantile")
#'   \item \code{n_ref, n_vars}: Size of reference dataset
#' }
#'
#' @details
#' **Common Use Cases:**
#'
#' \strong{Climate Velocity} (nearest geographic neighbor with similar climate):
#' \code{find_analogs(focal, ref, max_clim = 0.5, fun = "nmax", n = 1)}
#'
#' \strong{Climate Impact} (climatically similar locations within dispersal range):
#' \code{find_analogs(focal, ref, max_dist = 500, fun = "nmax", n = 10)}
#'
#' \strong{Analog Availability} (count of suitable locations):
#' \code{find_analogs(focal, ref, max_clim = 1.0, max_dist = 500, fun = "count")}
#'
#' **Climate Whitening:**
#' For best results, climate variables should be whitened (standardized and decorrelated)
#' before analysis. This ensures that Euclidean distance properly represents
#' Mahalanobis distance. Whitening functionality will be added in a future version;
#' for now, users should pre-process data externally.
#'
#' **Performance:**
#' The function uses quantile-based spatial binning to efficiently search large
#' reference datasets. Query time scales approximately O(log n) for top-k searches
#' and O(n^(d-1)/d) for constrained searches, where n is the reference dataset
#' size and d is the number of climate variables.
#'
#' @examples
#' \dontrun{
#' # Example 1: Climate velocity analysis
#' library(terra)
#' current <- rast("climate_1981_2010.tif")
#' future <- rast("climate_2071_2100.tif")
#'
#' focal <- as.data.frame(current, xy = TRUE)[1:100, ]
#' ref <- future
#'
#' velocity <- find_analogs(
#'   focal = focal,
#'   ref = ref,
#'   max_clim = 0.5,      # climate similarity threshold
#'   fun = "nmax",
#'   n = 1                # nearest neighbor
#' )
#'
#' # Example 2: Analog availability (count)
#' availability <- find_analogs(
#'   focal = focal,
#'   ref = ref,
#'   max_clim = 1.0,
#'   max_dist = 500,      # within 500 km
#'   fun = "count"
#' )
#'
#' # Example 3: Multiple analogs for impact assessment
#' analogs <- find_analogs(
#'   focal = focal,
#'   ref = ref,
#'   max_dist = 1000,
#'   fun = "nmax",
#'   n = 20
#' )
#'
#' # Check diagnostic information
#' cat("Index used", attr(analogs, "total_bins"), "bins\n")
#' cat("Average bin occupancy:", attr(analogs, "avg_bin_occupancy"), "\n")
#' }
#'
#' @export
find_analogs <- function(
      focal,
      ref,
      max_dist = NULL,
      max_clim = NULL,
      weight = "inverse_clim",
      fun = "nmax",
      n = NULL,
      metric = "euclidean",
      report_dist = TRUE,
      coord_type = c("auto", "lonlat", "projected")
) {
      # ---- Input validation --------------------------------------------------
      coord_type <- match.arg(coord_type)

      valid_weights <- c("uniform", "inverse_clim", "inverse_dist")
      if (!weight %in% valid_weights) {
            stop(
                  "weight must be one of: ",
                  paste(valid_weights, collapse = ", ")
            )
      }

      # TODO: Rename "nmax" to "topk" and add "argmin" per roadmap
      valid_funs <- c("nmax", "sum", "mean", "count", "all")
      if (!fun %in% valid_funs) {
            stop("fun must be one of: ", paste(valid_funs, collapse = ", "))
      }

      if (fun == "nmax" && is.null(n)) {
            stop("n must be specified when fun = 'nmax'")
      }

      if (!is.null(n)) {
            n <- as.integer(n)
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
      max_dist <- if (is.null(max_dist)) Inf else as.numeric(max_dist)[1L]

      max_clim <- if (is.null(max_clim)) {
            Inf
      } else {
            max_clim
      }

      # Map user's 'fun' parameter to internal k parameter
      k <- if (identical(fun, "nmax")) {
            n
      } else if (identical(fun, "all")) {
            nrow(ref_mm[, 3:ncol(ref_mm)]) # Return all reference points
      } else {
            0L
      }

      # ---- Call C++ core -----------------------------------------------------
      res <- .Call(
            `_analogs_find_analogs_core`,
            focal_mm, # matrix of focal sites, with xy and climate cols
            ref_mm, # matrix of ref sites, with xy and climate cols
            as.integer(k), # number of optima (nearest neighbors) to find
            max_clim, # climate filter bandwidth
            as.numeric(max_dist), # geographic distance filter bandwidth
            geo_mode # either "latlon" or "projected"
      )

      # Capture diagnostic attributes from C++ before post-processing
      cpp_attrs <- attributes(res)
      cpp_attrs$names <- NULL # Remove list element names
      cpp_attrs$class <- NULL # Remove class attribute

      # ---- Post-process results ----------------------------------------------
      if (fun %in% c("nmax", "all")) {
            out <- .emit_pairs(
                  res,
                  focal_mm,
                  ref_mm,
                  report_dist = report_dist,
                  weight = weight,
                  geo_mode = geo_mode
            )
            # Restore diagnostic attributes from C++
            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            return(out)
      }

      if (fun %in% c("sum", "mean", "count")) {
            out <- .emit_aggregates(
                  res,
                  focal_mm,
                  ref_mm,
                  fun = fun,
                  weight = weight,
                  geo_mode = geo_mode
            )
            # Restore diagnostic attributes from C++
            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
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

      # If coordinates fall within plausible lon/lat bounds, assume lonlat
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

#' Build long table of focal-analog pairs with distances and weights
#' @keywords internal
.emit_pairs <- function(res, focal_mm, ref_mm, report_dist, weight, geo_mode) {
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

            # Compute distances if needed for reporting or weighting
            if (isTRUE(report_dist) || !identical(weight, "uniform")) {
                  # Climate distances (Euclidean in whitened space)
                  v_i <- matrix(focal_mm[, 3:ncol(focal_mm)][i, ], nrow = 1)
                  v_j <- ref_mm[, 3:ncol(ref_mm)][idx, , drop = FALSE]
                  clim_d <- sqrt(rowSums((t(t(v_j) - as.numeric(v_i)))^2))

                  # Geographic distances
                  if (geo_mode == "lonlat") {
                        geog_d <- .haversine_km(cbind(fx, fy), cbind(ax, ay))
                  } else {
                        geog_d <- sqrt((ax - fx)^2 + (ay - fy)^2)
                  }

                  # Add distance columns if requested
                  if (isTRUE(report_dist)) {
                        df$clim_dist <- clim_d
                        df$geog_dist <- geog_d
                  }

                  # Compute weights
                  df$weight <- switch(
                        weight,
                        uniform = rep(1.0, length(idx)),
                        inverse_clim = 1.0 / (clim_d + 1e-12),
                        inverse_dist = 1.0 / (geog_d + 1e-12)
                  )
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

#' Aggregate matches per focal location
#' @keywords internal
.emit_aggregates <- function(res, focal_mm, ref_mm, fun, weight, geo_mode) {
      n_f <- nrow(focal_mm[, 1:2])

      out <- data.frame(
            focal_index = seq_len(n_f),
            focal_x = focal_mm[, 1:2][, 1],
            focal_y = focal_mm[, 1:2][, 2],
            value = NA_real_,
            stringsAsFactors = FALSE
      )

      for (i in seq_len(n_f)) {
            idx <- res[[i]]

            if (length(idx) == 0L) {
                  out$value[i] <- if (fun == "count") 0 else 0
                  next
            }

            # Compute distances for weighting
            v_i <- matrix(focal_mm[, 3:ncol(focal_mm)][i, ], nrow = 1)
            v_j <- ref_mm[, 3:ncol(ref_mm)][idx, , drop = FALSE]
            clim_d <- sqrt(rowSums((t(t(v_j) - as.numeric(v_i)))^2))

            if (geo_mode == "lonlat") {
                  fx <- rep(focal_mm[, 1:2][i, 1], length(idx))
                  fy <- rep(focal_mm[, 1:2][i, 2], length(idx))
                  geog_d <- .haversine_km(
                        cbind(fx, fy),
                        ref_mm[, 1:2][idx, 1:2, drop = FALSE]
                  )
            } else {
                  fx <- focal_mm[, 1:2][i, 1]
                  fy <- focal_mm[, 1:2][i, 2]
                  ax <- ref_mm[, 1:2][idx, 1]
                  ay <- ref_mm[, 1:2][idx, 2]
                  geog_d <- sqrt((ax - fx)^2 + (ay - fy)^2)
            }

            # Compute weights
            w <- switch(
                  weight,
                  uniform = rep(1, length(idx)),
                  inverse_clim = 1.0 / (clim_d + 1e-12),
                  inverse_dist = 1.0 / (geog_d + 1e-12)
            )

            # Aggregate
            out$value[i] <- switch(
                  fun,
                  count = length(idx), # Simple count
                  sum = sum(w),
                  mean = mean(w),
                  stop("Unrecognized aggregation function: ", fun)
            )
      }

      out
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
