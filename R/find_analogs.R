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

  # Normalize both inputs to matrix form
  focal_mm <- format_data(focal)
  ref_mm   <- format_data(ref)

  # Call the single fast C++ path
  find_analogs_matrix_matrix(
    focal_coords  = focal_mm$coords,
    focal_climate = focal_mm$climate,
    ref_coords    = ref_mm$coords,
    ref_climate   = ref_mm$climate,
    max_dist = max_dist,
    max_clim = max_clim,
    weight  = weight,
    fun     = fun,
    n       = n
  )
}


# Helper: select x,y by name if present; otherwise use first 2 cols
.select_xy_climate <- function(obj) {
  nm <- colnames(obj)
  if (!is.null(nm) && all(c("x","y") %in% nm[1:max(2, length(nm))])) {
    x_idx <- match("x", nm)
    y_idx <- match("y", nm)
    xy_idx <- c(x_idx, y_idx)
  } else if (!is.null(nm) && all(c("x","y") %in% nm)) {
    xy_idx <- match(c("x","y"), nm)
  } else {
    xy_idx <- 1:2
  }

  coords  <- as.matrix(obj[, xy_idx, drop = FALSE])
  climate <- as.matrix(obj[, setdiff(seq_len(ncol(obj)), xy_idx), drop = FALSE])

  storage.mode(coords)  <- "double"
  storage.mode(climate) <- "double"

  if (ncol(coords) != 2L) stop("coords must have exactly 2 columns (x,y).")
  if (ncol(climate) < 1L) stop("no climate variable columns found after x,y.")

  list(coords = coords, climate = climate)
}


#' Normalize input (SpatRaster OR data.frame/matrix) to coords/climate matrices
#' @keywords internal
format_data <- function(r) {
  if (inherits(r, "SpatRaster")) {
    df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
    .select_xy_climate(df)

  } else if (is.matrix(r) || is.data.frame(r)) {
    .select_xy_climate(r)

  } else {
    stop("Unsupported input type. Provide a data.frame/matrix with columns x,y,... or a SpatRaster.")
  }
}
