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
#' Both use the same format as their respective inputs. Output variables
#' inherit input names with a \code{_transformed} suffix (e.g. \code{tmean}
#' becomes \code{tmean_transformed}); each transformed axis is the whitened
#' counterpart of its corresponding input variable.
#'
#' Rows / cells with NA values in any coordinate or climate column are
#' excluded from the covariance estimation but preserved as NA in the
#' returned output, so the result has the same shape (row count / cell
#' count) as the input.
#'
#' @details
#' The transformation uses eigendecomposition of the covariance (or correlation)
#' matrix to decorrelate variables and scale to unit variance in all directions.
#' After transformation, Euclidean distance in the transformed space is
#' equivalent to Mahalanobis distance in the original space.
#'
#' The transformation is: X_transformed = (X - μ) %*% Σ^(-1/2), computed via
#' symmetric (ZCA / Mahalanobis) whitening: \eqn{\Sigma^{-1/2} = V D^{-1/2} V'}.
#' This is the unique whitening that keeps the transformed axes as close as
#' possible (in least-squares sense) to the original variables, so each
#' output column corresponds one-to-one with its input counterpart. PCA
#' whitening would give identical pairwise distances but reorient the axes
#' to principal components.
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
#' analogs <- analog_search(
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

      x_is_raster    <- inherits(x, "SpatRaster")
      pool_is_raster <- !is.null(pool) && inherits(pool, "SpatRaster")

      if (x_is_raster || pool_is_raster) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs",
                       call. = FALSE)
            }
      }

      # Prepare each input into a uniform internal representation:
      #   $coords_full : full-length xy coords (NAs preserved)   -- matrix/df only
      #   $clim_full   : full-length climate matrix (NAs preserved) -- matrix/df only
      #   $clim_strip  : NA-stripped climate matrix used for the math
      #   $keep        : logical vector marking kept rows
      #   $n_total     : original row / cell count
      #   $template    : single-layer NA SpatRaster (raster inputs only)
      #   $clim_names  : climate variable names (post-`.select_xy_climate`)
      x_prep    <- .prep_for_mahal(x,    "x")
      pool_prep <- if (!is.null(pool))  .prep_for_mahal(pool, "pool")  else NULL

      if (!is.null(pool_prep) &&
          ncol(x_prep$clim_strip) != ncol(pool_prep$clim_strip)) {
            stop("x and pool must have the same number of climate variables",
                 call. = FALSE)
      }

      n_clim <- ncol(x_prep$clim_strip)

      # Combine the stripped climate data for shared covariance structure.
      combined <- if (!is.null(pool_prep)) {
            rbind(x_prep$clim_strip, pool_prep$clim_strip)
      } else {
            x_prep$clim_strip
      }

      # Centering / scaling parameters from the (stripped) combined data.
      clim_means <- if (center) colMeans(combined) else rep(0, n_clim)
      clim_sds   <- if (scale)  apply(combined, 2, stats::sd) else rep(1, n_clim)

      if (any(clim_sds < 1e-10)) {
            stop("One or more climate variables has zero or near-zero variance",
                 call. = FALSE)
      }

      # Standardize.
      combined_std <- sweep(combined,  2, clim_means, FUN = "-")
      combined_std <- sweep(combined_std, 2, clim_sds, FUN = "/")

      # Covariance (correlation when scale = TRUE) of standardized data.
      # All NAs have been stripped, so plain cov() is fine.
      cov_mat <- stats::cov(combined_std)

      # Eigendecomposition; regularize if near-singular.
      eig <- eigen(cov_mat, symmetric = TRUE)
      max_eval <- max(eig$values)
      min_eval <- min(eig$values)
      if (min_eval < 1e-8 * max_eval) {
            warning("Covariance matrix is near-singular. ",
                    "Adding small regularization term.", call. = FALSE)
            eig$values <- eig$values + 1e-6 * max_eval
      }

      # Whitening matrix: V * D^(-1/2) * V'
      inv_sqrt_evals    <- 1.0 / sqrt(eig$values)
      whitening_matrix  <- eig$vectors %*%
            diag(inv_sqrt_evals, nrow = n_clim) %*%
            t(eig$vectors)

      # Apply transform to each stripped input, then scatter back to full
      # length with NAs at dropped rows.
      x_full <- .apply_and_scatter(x_prep,    clim_means, clim_sds,
                                   whitening_matrix)
      pool_full <- if (!is.null(pool_prep)) {
            .apply_and_scatter(pool_prep, clim_means, clim_sds,
                               whitening_matrix)
      } else NULL

      # Reconstruct outputs in original format.
      x_result <- .reconstruct_mahal_output(x, x_prep, x_full)

      if (!is.null(pool)) {
            pool_result <- .reconstruct_mahal_output(pool, pool_prep, pool_full)
            return(list(x = x_result, pool = pool_result))
      }
      x_result
}


# ---- internal helpers --------------------------------------------------------

# Normalize input and produce both full-length and NA-stripped views of the
# climate data. Centralizes everything `mahalanobis_transform()` needs to know
# about a single input.
#
# `label` is just used in error messages.
#
# Reuses `.select_xy_climate()` (existing helper in utils.R) so column-name
# / xy-detection rules match the rest of the package.
.prep_for_mahal <- function(obj, label) {

      is_raster <- inherits(obj, "SpatRaster")

      if (is_raster) {
            # Pull xy + climate as a full-length data.frame (NAs preserved).
            df <- terra::as.data.frame(obj, xy = TRUE, na.rm = FALSE)
            normalized <- .select_xy_climate(df)
            template   <- stats::setNames(
                  terra::setValues(obj[[1]], NA), "raster"
            )
      } else if (is.matrix(obj) || is.data.frame(obj)) {
            normalized <- .select_xy_climate(obj)
            template   <- NULL
      } else {
            stop("`", label, "` must be a data.frame, matrix, or SpatRaster.",
                 call. = FALSE)
      }

      coords_full <- normalized[, 1:2, drop = FALSE]
      clim_full   <- normalized[, -(1:2), drop = FALSE]

      keep       <- stats::complete.cases(normalized)
      clim_strip <- clim_full[keep, , drop = FALSE]

      if (nrow(clim_strip) == 0L) {
            stop("`", label, "` has no rows with complete (non-NA) values.",
                 call. = FALSE)
      }

      # `.select_xy_climate()` doesn't synthesize column names when the input
      # lacks them (e.g., bare unnamed matrix), so fall back to clim1, clim2…
      clim_names <- colnames(clim_full)
      if (is.null(clim_names) || length(clim_names) != ncol(clim_full)) {
            clim_names <- paste0("clim", seq_len(ncol(clim_full)))
      }

      list(
            coords_full = coords_full,
            clim_full   = clim_full,
            clim_strip  = clim_strip,
            keep        = keep,
            n_total     = nrow(normalized),
            template    = template,
            clim_names  = clim_names,
            is_raster   = is_raster
      )
}


# Apply (center / scale / whiten) to the stripped climate matrix and scatter
# the result back into a full-length matrix with NA at dropped rows.
.apply_and_scatter <- function(prep, clim_means, clim_sds, whitening_matrix) {

      strip <- sweep(prep$clim_strip,  2, clim_means, FUN = "-")
      strip <- sweep(strip,            2, clim_sds,   FUN = "/")
      strip <- strip %*% whitening_matrix

      out <- matrix(NA_real_,
                    nrow = prep$n_total,
                    ncol = ncol(strip))
      out[prep$keep, ] <- strip

      # Preserve column names from the stripped (i.e. post-`.select_xy_climate`)
      # ordering so they line up with the original climate columns.
      colnames(out) <- paste0(prep$clim_names, "_transformed")
      out
}


# Build the final output in the same format as the original input. `clim_full`
# is the full-length, NA-padded transformed climate matrix produced by
# `.apply_and_scatter()`.
.reconstruct_mahal_output <- function(original, prep, clim_full) {

      out_names <- colnames(clim_full)
      n_clim    <- ncol(clim_full)

      if (prep$is_raster) {
            # Build one layer per climate variable using the package's standard
            # template-based rasterization pattern (cf. `.cv_to_raster()`).
            template <- prep$template
            layers <- lapply(seq_len(n_clim), function(j) {
                  stats::setNames(
                        terra::setValues(template, clim_full[, j]),
                        out_names[j]
                  )
            })
            out_rast <- do.call(c, layers)
            terra::varnames(out_rast) <- out_names
            return(out_rast)
      }

      # matrix / data.frame branch: rebuild xy + transformed climate at full
      # length, preserving the input's class.
      coord_names <- colnames(prep$coords_full)
      if (is.null(coord_names)) coord_names <- c("x", "y")

      combined <- cbind(prep$coords_full, clim_full)
      colnames(combined) <- c(coord_names, out_names)

      if (is.data.frame(original)) {
            return(as.data.frame(combined))
      }
      combined
}
