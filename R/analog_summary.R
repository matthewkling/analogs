#' Summarize analog search results and parameters
#' @param x Results from `analog_*()` search functions
#' @export
analog_summary <- function(x) {
      cat("Analog Search Results\n")
      cat("=====================\n\n")

      # Data structure
      is_raster <- inherits(x, "SpatRaster")
      if (is_raster) {
            cat("Format: SpatRaster\n")
            if(!requireNamespace("terra", quietly = TRUE)) return(invisible())
            cat("  Dimensions:", paste(dim(x), collapse = " x "), "\n")
            cat("  Variables:", paste(terra::names(x), collapse = ", "), "\n")
      } else {
            cat("Format: data.frame\n")
            cat("  Dimensions:", paste(dim(x), collapse = " x "), "\n")
            cat("  Variables:", paste(names(x), collapse = ", "), "\n")
      }

      # Input data information
      n_focal <- attr(x, "n_focal", exact = TRUE)
      n_ref <- attr(x, "n_ref", exact = TRUE)
      n_clim <- attr(x, "n_clim", exact = TRUE)
      if (!is.null(n_focal) || !is.null(n_ref) || !is.null(n_clim)) {
            cat("\nInput dataset sizes:\n")
            if (!is.null(n_focal)) {
                  cat("  Focal locations:", format(n_focal, big.mark = ","), "\n")
            }
            if (!is.null(n_ref)) {
                  cat("  Reference pool locations:", format(n_ref, big.mark = ","), "\n")
            }
            if (!is.null(n_clim)) {
                  cat("  Climate variables:", n_clim, "\n")
            }
      }

      # Selection parameters
      select <- attr(x, "select", exact = TRUE)
      max_clim <- attr(x, "max_clim", exact = TRUE)
      max_geog <- attr(x, "max_dist", exact = TRUE)
      k <- attr(x, "k", exact = TRUE)
      cat("\nSelection parameters:\n")
      if (!is.null(select)) {
            cat("  Method:", select, "\n")
      }
      if (!is.null(max_clim)) {
            if (length(max_clim) == 1 && is.finite(max_clim)) {
                  cat("  Max climate distance:", max_clim, "\n")
            } else if (length(max_clim) > 1) {
                  cat("  Max climate distance: per-variable\n")
            } else {
                  cat("  Max climate distance: none\n")
            }
      }
      if (!is.null(max_geog)) {
            if (is.finite(max_geog)) {
                  geo_mode <- attr(x, "geo_mode", exact = TRUE)
                  units <- if (!is.null(geo_mode) && geo_mode == "lonlat") "km" else "units"
                  cat("  Max geographic distance:", max_geog, units, "\n")
            } else {
                  cat("  Max geographic distance: none\n")
            }
      }
      if (!is.null(k)) {
            cat("  k:", k, "nearest neighbors\n")
      }



      # Aggregation
      stat <- attr(x, "stat", exact = TRUE)
      cat("\nAggregation parameters:\n")
      if (!is.null(stat)) {
            cat("  Stats:", paste(stat, collapse = ", "), "\n")
      }

      weight <- attr(x, "weight", exact = TRUE)
      if (!is.null(weight)) {
            cat("  Weighting:", weight, "\n")
            theta <- attr(x, "theta", exact = TRUE)
            if (!is.null(theta)) {
                  cat("    theta:", paste(theta, collapse = ", "), "\n")
            }
      }

      # Coordinate system / units
      geo_mode <- attr(x, "geo_mode", exact = TRUE)
      clim_metric <- attr(x, "x_cov")
      if (!is.null(geo_mode) && !is.null(clim_metric)) {
            units <- if (geo_mode == "lonlat") "km" else "projection units"
            cat("\nDistance metrics:\n")
            cat("  Spatial coordinate system:", geo_mode, "\n")
            cat("  Geographic distance units:", units, "\n")
            cat("  Climate distance measure:",
                if(clim_metric) "site-specific Mahalanobis distance" else "Euclidean distance", "\n")
      }

      # Lattice index diagnostics
      n_cells_nonempty <- attr(x, "n_bins_nonempty", exact = TRUE)
      mean_occ <- signif(attr(x, "avg_nonempty_bin_occupancy", exact = TRUE), 3)
      binning <- attr(x, "binning_method", exact = TRUE)
      total_bins <- attr(x, "total_bins", exact = TRUE)
      if (!is.null(binning) && binning == "lattice_ecef") {
            n_dims <- 3 + n_clim
      } else {
            n_dims <- 2 + n_clim
      }
      index_res <- total_bins ^ (1 / n_dims)

      if (!is.null(index_res)) {
            cat("\nLattice index:\n")
            if (!is.null(index_res)) {
                  cat("  Resolution:", index_res, "bins/dim\n")
            }
            if (!is.null(n_cells_nonempty)) {
                  cat("  Total occupied bins:", format(n_cells_nonempty, big.mark = ","), "\n")
            }
            if (!is.null(mean_occ)) {
                  cat("  Mean occupancy:", mean_occ, "points/bin\n")
            }
      }

      # Data summary
      cat("\nData summary:\n")
      v <- setdiff(names(x), c("x", "y", "index", "analog_x", "analog_y", "analog_index"))
      if(is_raster) print(terra::summary(x[[v]])) else print(summary(x[, v]))

      cat("\nSee `attributes(...)` for complete metadata.\n\n")
      invisible(x)
}
