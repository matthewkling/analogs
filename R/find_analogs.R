#' Find climate analogs
#'
#' @param focal Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers
#' @param ref Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers
#' @param max_dist Maximum geographic distance constraint in km (NULL = no constraint)
#' @param max_clim Maximum climate distance constraint (NULL = no constraint)
#' @param weight Weighting function: "uniform", "inverse_clim", or "inverse_dist"
#' @param fun Aggregation function: "nmax", "sum", "mean", "count", or "all"
#' @param n Number of analogs to return (required if fun = "nmax")
#'
#' @return Data frame with columns depending on fun parameter:
#'   - fun = "nmax" or "all": focal_index, focal_x, focal_y, analog_index,
#'     analog_x, analog_y, clim_dist, geog_dist, weight
#'   - fun = "sum", "mean", "count": focal_index, focal_x, focal_y, value
#'
#' @export
find_analogs <- function(focal,
                         ref,
                         max_dist = NULL,
                         max_clim = NULL,
                         weight = "inverse_clim",
                         fun = "nmax",
                         n = NULL) {

      # Validate parameters
      valid_weights <- c("uniform", "inverse_clim", "inverse_dist")
      if (!weight %in% valid_weights) {
            stop("weight must be one of: ", paste(valid_weights, collapse = ", "))
      }

      valid_funs <- c("nmax", "sum", "mean", "count", "all")
      if (!fun %in% valid_funs) {
            stop("fun must be one of: ", paste(valid_funs, collapse = ", "))
      }

      if (fun == "nmax" && is.null(n)) {
            stop("n must be specified when fun = 'nmax'")
      }

      # Detect input types
      focal_type <- if (inherits(focal, "SpatRaster")) "raster" else "matrix"
      ref_type <- if (inherits(ref, "SpatRaster")) "raster" else "matrix"

      # Dispatch to appropriate C++ function
      if (focal_type == "matrix" && ref_type == "raster") {

            # Matrix focal, Raster ref
            ref_array <- terra::as.array(ref)

            # Get focal coordinates for output
            focal_coords_for_output <- as.matrix(focal[, c("x", "y"), drop = FALSE])

            # Get focal cell indices (snapped to raster grid)
            focal_cells <- terra::cellFromXY(ref, focal[, c("x", "y"), drop = FALSE])
            focal_climate <- as.matrix(focal[, -c(1:2), drop = FALSE])

            # Get candidate cell coordinates for distance calculations
            ref_coords <- terra::xyFromCell(ref, 1:terra::ncell(ref))

            res <- terra::res(ref)
            dims <- dim(ref_array)[1:2]

            result <- find_analogs_matrix_raster(
                  focal_cells = focal_cells,
                  focal_coords = focal_coords_for_output,
                  focal_climate = focal_climate,
                  ref_array = ref_array,
                  ref_coords = ref_coords,
                  xres = res[1],
                  yres = res[2],
                  nrow = dims[1],
                  ncol = dims[2],
                  max_dist = max_dist,
                  max_clim = max_clim,
                  weight = weight,
                  fun = fun,
                  n = n
            )

      } else if (focal_type == "matrix" && ref_type == "matrix") {

            # Matrix focal, Matrix ref
            focal_coords <- as.matrix(focal[, c("x", "y"), drop = FALSE])
            focal_climate <- as.matrix(focal[, -c(1:2), drop = FALSE])
            ref_coords <- as.matrix(ref[, c("x", "y"), drop = FALSE])
            ref_climate <- as.matrix(ref[, -c(1:2), drop = FALSE])

            result <- find_analogs_matrix_matrix(
                  focal_coords = focal_coords,
                  focal_climate = focal_climate,
                  ref_coords = ref_coords,
                  ref_climate = ref_climate,
                  max_dist = max_dist,
                  max_clim = max_clim,
                  weight = weight,
                  fun = fun,
                  n = n
            )

      } else if (focal_type == "raster" && ref_type == "matrix") {

            # Raster focal, Matrix ref
            focal_array <- terra::as.array(focal)
            focal_coords <- terra::crds(focal, na.rm = FALSE)
            ref_coords <- as.matrix(ref[, c("x", "y"), drop = FALSE])
            ref_climate <- as.matrix(ref[, -c(1:2), drop = FALSE])

            result <- find_analogs_raster_matrix(
                  focal_array = focal_array,
                  focal_coords = focal_coords,
                  ref_coords = ref_coords,
                  ref_climate = ref_climate,
                  max_dist = max_dist,
                  max_clim = max_clim,
                  weight = weight,
                  fun = fun,
                  n = n
            )

      } else if (focal_type == "raster" && ref_type == "raster") {
            stop("Raster-raster case not yet implemented. Use matrix focal or matrix ref.")
      }

      return(result)
}
