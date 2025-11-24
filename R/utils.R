
# ---- Internal Helper Functions ---------------------------------------------

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
