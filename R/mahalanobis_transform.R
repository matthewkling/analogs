#' Transform Climate Data for Global Mahalanobis Distance
#'
#' Transforms climate data to decorrelated, unit-variance space, enabling the
#' use of Euclidean distance as a global Mahalanobis distance. Useful for
#' climate analog analysis when you want to account for spatial correlation
#' structure and/or standardize variables with different units.
#'
#' @param x Climate data for focal/query points. Can be:
#'   \itemize{
#'     \item Matrix/data.frame with columns x, y, and climate variables
#'     \item SpatRaster with climate variable layers
#'   }
#'
#' @param pool Optional reference climate data to transform jointly with \code{x}.
#'   When provided, both datasets are transformed using the same transformation
#'   (computed from their combined covariance structure). Must have the same
#'   number of climate variables as \code{x}. Same format options as \code{x}.
#'
#' @param center Logical; if TRUE (default), center each variable to mean zero
#'   before transformation.
#'
#' @param scale Logical; if TRUE (default), standardize each variable to unit
#'   variance before transformation (correlation-based). If FALSE, use
#'   covariance-based transformation which preserves relative variance magnitudes.
#'   Generally TRUE is recommended when climate variables have different units
#'   (e.g., temperature in °C vs precipitation in mm).
#'
#' @return
#' If \code{pool = NULL}: Returns transformed version of \code{x} in the same
#' format (matrix/data.frame/SpatRaster).
#'
#' If \code{pool} is provided: Returns a list with two components:
#' \itemize{
#'   \item \code{$x}: Transformed version of \code{x}
#'   \item \code{$pool}: Transformed version of \code{pool}
#' }
#' Both use the same format as their respective inputs.
#'
#' @details
#' The transformation uses eigendecomposition of the covariance (or correlation)
#' matrix to decorrelate variables and scale to unit variance in all directions.
#' After transformation, Euclidean distance in the transformed space is
#' equivalent to Mahalanobis distance in the original space.
#'
#' The transformation is: X_transformed = (X - μ) %*% Σ^(-1/2)
#'
#' When \code{pool} is provided, the covariance structure is computed from the
#' combined dataset to ensure both are transformed to the same space.
#'
#' If the covariance matrix is near-singular, a small regularization term
#' (1e-6 * max eigenvalue) is added to stabilize the inversion.
#'
#' @examples
#' \dontrun{
#' # Single dataset transformation
#' clim_transformed <- mahalanobis_transform(climate_data)
#'
#' # Joint transformation for analog analysis
#' transformed <- mahalanobis_transform(
#'   x = baseline_climate,
#'   pool = future_climate
#' )
#'
#' # Use transformed data with Euclidean distance (= global Mahalanobis)
#' analogs <- find_analogs(
#'   x = transformed$x,
#'   pool = transformed$pool,
#'   mode = "knn_geog",
#'   max_clim = 2,  # Now in standardized units
#'   k = 1
#' )
#'
#' # Covariance-based transformation (preserve relative variances)
#' transformed_cov <- mahalanobis_transform(
#'   x = climate_data,
#'   scale = FALSE
#' )
#' }
#'
#' @export
mahalanobis_transform <- function(x, pool = NULL, center = TRUE, scale = TRUE) {

      # Check if terra is needed
      x_is_raster <- inherits(x, "SpatRaster")
      pool_is_raster <- !is.null(pool) && inherits(pool, "SpatRaster")

      if (x_is_raster || pool_is_raster) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs", call. = FALSE)
            }
      }

      # Extract climate data using existing helpers
      x_formatted <- .format_data(x)
      x_coords <- x_formatted[, 1:2, drop = FALSE]
      x_clim <- x_formatted[, -(1:2), drop = FALSE]
      n_clim <- ncol(x_clim)

      pool_formatted <- NULL
      pool_coords <- NULL
      pool_clim <- NULL

      if (!is.null(pool)) {
            pool_formatted <- .format_data(pool)
            pool_coords <- pool_formatted[, 1:2, drop = FALSE]
            pool_clim <- pool_formatted[, -(1:2), drop = FALSE]

            if (ncol(pool_clim) != n_clim) {
                  stop("x and pool must have the same number of climate variables")
            }
      }

      # Combine for shared covariance structure if paired
      if (!is.null(pool_clim)) {
            combined <- rbind(x_clim, pool_clim)
      } else {
            combined <- x_clim
      }

      # Compute centering and scaling parameters from combined data
      clim_means <- if (center) colMeans(combined, na.rm = TRUE) else rep(0, n_clim)
      clim_sds <- if (scale) apply(combined, 2, sd, na.rm = TRUE) else rep(1, n_clim)

      # Check for zero variance
      if (any(clim_sds < 1e-10)) {
            stop("One or more climate variables has zero or near-zero variance")
      }

      # Center and scale
      combined_std <- combined
      for (i in seq_len(n_clim)) {
            combined_std[, i] <- (combined[, i] - clim_means[i]) / clim_sds[i]
      }

      # Compute covariance matrix (correlation matrix if scale=TRUE)
      cov_mat <- stats::cov(combined_std, use = "complete.obs")

      # Eigendecomposition
      eig <- eigen(cov_mat, symmetric = TRUE)

      # Check for near-singular matrix
      max_eval <- max(eig$values)
      min_eval <- min(eig$values)

      if (min_eval < 1e-8 * max_eval) {
            warning("Covariance matrix is near-singular. Adding small regularization term.")
            eig$values <- eig$values + 1e-6 * max_eval
      }

      # Compute whitening matrix: V * D^(-1/2) * V'
      # where cov_mat = V * D * V'
      inv_sqrt_evals <- 1.0 / sqrt(eig$values)
      whitening_matrix <- eig$vectors %*% diag(inv_sqrt_evals, nrow = n_clim) %*% t(eig$vectors)

      # Apply transformation to x
      x_clim_std <- x_clim
      for (i in seq_len(n_clim)) {
            x_clim_std[, i] <- (x_clim[, i] - clim_means[i]) / clim_sds[i]
      }
      x_clim_transformed <- x_clim_std %*% whitening_matrix

      # Apply transformation to pool if provided
      pool_clim_transformed <- NULL
      if (!is.null(pool_clim)) {
            pool_clim_std <- pool_clim
            for (i in seq_len(n_clim)) {
                  pool_clim_std[, i] <- (pool_clim[, i] - clim_means[i]) / clim_sds[i]
            }
            pool_clim_transformed <- pool_clim_std %*% whitening_matrix
      }

      # Reconstruct outputs in original format
      x_result <- .reconstruct_output(x, x_coords, x_clim_transformed, x_is_raster)

      if (!is.null(pool)) {
            pool_result <- .reconstruct_output(pool, pool_coords, pool_clim_transformed, pool_is_raster)
            return(list(x = x_result, pool = pool_result))
      } else {
            return(x_result)
      }
}


#' Reconstruct output in original format
#' @keywords internal
.reconstruct_output <- function(original, coords, clim_transformed, is_raster) {

      if (is_raster) {
            # Reconstruct as SpatRaster
            combined <- cbind(coords, clim_transformed)

            # Get original raster properties
            ncol_orig <- terra::ncol(original)
            nrow_orig <- terra::nrow(original)
            ext_orig <- terra::ext(original)
            crs_orig <- terra::crs(original)

            # Create output raster
            out_rast <- terra::rast(
                  nrows = nrow_orig,
                  ncols = ncol_orig,
                  extent = ext_orig,
                  crs = crs_orig,
                  nlyr = ncol(clim_transformed)
            )

            # Set values
            terra::values(out_rast) <- clim_transformed

            # Set layer names
            orig_names <- names(original)
            if (length(orig_names) == ncol(clim_transformed)) {
                  names(out_rast) <- paste0(orig_names, "_transformed")
            } else {
                  names(out_rast) <- paste0("clim", seq_len(ncol(clim_transformed)), "_transformed")
            }

            return(out_rast)

      } else if (is.matrix(original) || is.data.frame(original)) {
            # Reconstruct as matrix or data.frame
            combined <- cbind(coords, clim_transformed)

            # Try to preserve column names
            coord_names <- colnames(coords)
            if (is.null(coord_names)) {
                  coord_names <- c("x", "y")
            }

            clim_names <- colnames(original)
            if (!is.null(clim_names) && length(clim_names) == ncol(original)) {
                  # Original had names, use them for climate columns (skip x,y)
                  clim_names_subset <- clim_names[-(1:2)]
                  colnames(combined) <- c(coord_names, paste0(clim_names_subset, "_transformed"))
            } else {
                  colnames(combined) <- c(coord_names,
                                          paste0("clim", seq_len(ncol(clim_transformed)), "_transformed"))
            }

            if (is.data.frame(original)) {
                  return(as.data.frame(combined))
            } else {
                  return(combined)
            }

      } else {
            stop("Unsupported input type")
      }
}
