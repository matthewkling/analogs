#' Tiled analog search for memory-constrained queries
#'
#' Performs analog searches on large raster datasets by dividing the focal
#' region into tiles and processing each tile separately. This reduces memory
#' usage at the cost of increasing compute time. Works with any `analog_*()`
#' function.
#'
#' @param x SpatRaster with focal locations (points to find analogs for).
#' @param pool SpatRaster with reference locations (potential analog pool).
#' @param n_tiles Approximate number of tiles. The function will find a grid
#'   close to this number that creates square-ish tiles. Choosing larger
#'   values for n_tiles will reduce memory usage, but will also reduce
#'   computational efficiency. Choose the smallest n_tiles that fits your
#'   memory constraints.
#' @param fun An analog_* function to apply to each tile (e.g., analog_velocity,
#'   analog_impact).
#' @param geog A [kernel()] object giving the geographic search treatment. Its
#'   `max` (the geographic search radius) is required and finite: it both sizes
#'   the tile halos (so edge focals still see all in-range analogs) and is
#'   forwarded to `fun`. Tiling is only beneficial when the geographic search
#'   is bounded, hence `geog` with a finite `max` is required.
#' @param clim Optional [kernel()] object giving the climate search treatment,
#'   forwarded to `fun`. `NULL` (default) means no climate kernel is passed.
#' @param y Optional SpatRaster with values to aggregate across analogs.
#'   Must have spatial properties matching pool.
#' @param x_cov Optional SpatRaster with covariates for focal points. Must have
#'   spatial properties matching x.
#' @param weight Optional single-layer SpatRaster of per-pool-cell weights.
#'   Must have the same CRS and extent as `pool`. See [analog_search()] for
#'   details.
#' @param cell_area_weight Controls cell-area weighting. One of `"auto"`
#'   (default; on for raster pools, which is always the case here), `TRUE`,
#'   or `FALSE`. Unlike a non-tiled query, where weights are normalized to
#'   mean 1 over whatever pool is passed, here weights are normalized to
#'   mean 1 over the *full pool raster* once at the outer level, and the
#'   resulting raster is sliced per tile. This keeps absolute magnitudes of
#'   weighted stats (e.g. `sum_weights`) consistent across tiles.
#' @param ... Additional arguments passed to fun (e.g. `select`, `stat`, `k`).
#'   May include `covariates` as a SpatRaster matching `pool`'s CRS and
#'   extent (cropped per tile alongside `y`). May also include `normalize`
#'   for helpers that support it (e.g. `analog_density()`); when normalization
#'   is requested, the global mean cell area is computed once over the full
#'   pool and propagated to each tile so per-tile `D_max` values are
#'   consistent.
#' @param output_file Optional filename for disk-based output. If specified and
#'   fun returns a SpatRaster, tiles are written to temporary files during
#'   processing and merged to output_file at the end. This is useful when
#'   results are too large to fit in memory. Ignored for data.frame results.
#' @param progress Logical indicating whether to show progress bar.
#'
#' @return Same type as fun returns (SpatRaster or data.frame). If output_file
#'   is specified, returns a disk-backed SpatRaster.
#'
#' @details
#' Tiled analog searches work by splitting x into a number of smaller tiles and
#' calling the requested analog function on each tile, using an analog pool
#' that is the size of the tile buffered by the geographic radius (`geog`'s
#' `max`). This buffer is necessary
#' for correctness but increases compute time, particularly if the radius is
#' large. The results for each tile are temporarily written to disk, and are
#' merged into a single results raster once all tiles have processed.
#'
#' The function requires `geog` to have a finite `max`, as tiling is only beneficial when
#' geographic distance constraints limit the reference pool size for each
#' focal point. The function will warn if the radius is so large that tiling
#' provides minimal memory benefit.
#'
#' If index_res is specified in ..., all tiles will use the same lattice
#' resolution. If index_res is not specified, each tile will independently
#' auto-tune its lattice resolution based on local data characteristics. This
#' adaptive behavior is generally fine and can even be beneficial when climate
#' distributions vary substantially across the landscape (e.g., mountains vs
#' plains).
#'
#' @export
tiled_analog_search <- function(
            x,
            pool,
            n_tiles,
            fun,
            geog,
            clim = NULL,
            y = NULL,
            x_cov = NULL,
            weight = NULL,
            cell_area_weight = "auto",
            ...,
            output_file = NULL,
            progress = TRUE
) {
      # Validation
      if (!inherits(x, "SpatRaster")) {
            stop("x must be a SpatRaster for tiled_analog_search")
      }
      if (!inherits(pool, "SpatRaster")) {
            stop("pool must be a SpatRaster for tiled_analog_search")
      }
      if (!requireNamespace("terra", quietly = TRUE)) {
            stop("package `terra` is required for tiled_analog_search")
      }

      # Check CRS compatibility
      if (!terra::same.crs(x, pool)) {
            stop("x and pool must have the same CRS")
      }

      # Tiling requires a bounded geographic search radius so that halos can be
      # sized to guarantee edge focals still see all in-range analogs. Extract
      # the numeric radius from the `geog` kernel; the internal tiling geometry
      # below works with this scalar.
      if (!inherits(geog, "analog_kernel")) {
            stop("tiled_analog_search requires `geog` to be a kernel() object ",
                 "with a finite `max` (the geographic search radius).",
                 call. = FALSE)
      }
      max_geog <- geog$max
      if (is.null(max_geog) || !is.finite(max_geog) || max_geog <= 0) {
            stop("tiled_analog_search requires `geog` to have a finite, positive ",
                 "`max`: e.g. geog = kernel(max = 100). Tiling is only beneficial ",
                 "when the geographic search is bounded.", call. = FALSE)
      }
      if (!is.null(clim) && !inherits(clim, "analog_kernel")) {
            stop("`clim` must be a kernel() object or NULL.", call. = FALSE)
      }

      # Validate and check y values if provided
      if (!is.null(y)) {
            if (!inherits(y, "SpatRaster")) {
                  stop("y must be a SpatRaster for tiled_analog_search")
            }
            if (!terra::same.crs(y, pool)) {
                  stop("y must have the same CRS as pool")
            }
            if (!all(terra::ext(y) == terra::ext(pool))) {
                  stop("y must have the same extent as pool")
            }
      }

      # Validate and check x_cov if provided
      if (!is.null(x_cov)) {
            if (!inherits(x_cov, "SpatRaster")) {
                  stop("x_cov must be a SpatRaster for tiled_analog_search")
            }
            if (!terra::same.crs(x_cov, x)) {
                  stop("x_cov must have the same CRS as x")
            }
            if (!all(terra::ext(x_cov) == terra::ext(x))) {
                  stop("x_cov must have the same extent as x")
            }
      }

      # Validate and check weight if provided. Same shape rules as y: a
      # single-layer SpatRaster with the same CRS and extent as pool.
      if (!is.null(weight)) {
            if (!inherits(weight, "SpatRaster")) {
                  stop("weight must be a SpatRaster for tiled_analog_search")
            }
            if (terra::nlyr(weight) != 1L) {
                  stop("weight must have exactly one layer")
            }
            if (!terra::same.crs(weight, pool)) {
                  stop("weight must have the same CRS as pool")
            }
            if (!all(terra::ext(weight) == terra::ext(pool))) {
                  stop("weight must have the same extent as pool")
            }
      }

      # Pull `covariates` out of `...` so we can crop it per tile (it's a
      # pool-side raster, like y). Without this, passing a SpatRaster
      # `covariates` to tiled_analog_search would error out at the helper
      # level: each tile's pool is a halo crop, but `covariates` would be
      # the full-extent raster.
      dot_args   <- list(...)
      covariates <- dot_args$covariates
      dot_args$covariates <- NULL

      if (!is.null(covariates) && inherits(covariates, "SpatRaster")) {
            if (!terra::same.crs(covariates, pool)) {
                  stop("covariates must have the same CRS as pool")
            }
            if (!all(terra::ext(covariates) == terra::ext(pool))) {
                  stop("covariates must have the same extent as pool")
            }
      } else if (!is.null(covariates)) {
            stop("In tiled_analog_search, `covariates` must be a SpatRaster ",
                 "(matching pool's extent and CRS).", call. = FALSE)
      }

      # Resolve cell_area_weight. We need globally consistent normalization
      # across all tiles, so we precompute a normalized cell-area raster on
      # the full pool here, then crop it per tile and feed the cropped
      # values to each tile's index build via the numeric-vector form of
      # `cell_area_weight`. This keeps absolute magnitudes of weighted
      # stats (`sum_weights` etc.) comparable across tiles.
      if (!(identical(cell_area_weight, "auto") ||
            isTRUE(cell_area_weight) ||
            isFALSE(cell_area_weight))) {
            stop('`cell_area_weight` must be "auto", TRUE, or FALSE in ',
                 "tiled_analog_search.", call. = FALSE)
      }
      apply_caw <- if (identical(cell_area_weight, "auto")) {
            TRUE                                     # pool is always a raster here
      } else {
            isTRUE(cell_area_weight)
      }

      # Determine coord_type once, here in the wrapper, using the same
      # detection logic the rest of the package uses (.detect_geo respects
      # CRS metadata first and falls back to coordinate-magnitude only when
      # CRS is unknown). The resolved value is then forwarded explicitly
      # to every per-tile analog_search() call, so all tiles agree even if
      # individual tiles' coordinate ranges would have triggered a
      # different fallback magnitude check than the full raster would.
      #
      # We pass the raster's extent corners (a 2x2 matrix) as the xy
      # fallback. This is consulted by .detect_geo() only when CRS is
      # unknown (terra::is.lonlat returns NA) and gives the magnitude
      # check the full data range without materializing every cell.
      xy_for_detect <- {
            ex <- terra::ext(x)
            matrix(
                  c(ex[1], ex[2], ex[3], ex[4]),
                  nrow = 2, ncol = 2,
                  dimnames = list(NULL, c("x", "y"))
            )
      }
      coord_type <- .detect_geo(xy_for_detect, input = x)
      if (!coord_type %in% c("lonlat", "projected")) {
            stop("Could not determine coord_type from `x`. ",
                 "Provide a SpatRaster with a recognized CRS, or set the ",
                 "CRS explicitly before calling tiled_analog_search().",
                 call. = FALSE)
      }

      # Compute the globally-consistent normalized cell-area raster AND the
      # global mean cell area scalar in one pass. Both must use the same
      # subset of cells (finite cellSize values over the full pool) for
      # consistency between per-tile aggregation and per-tile D_max
      # normalization.
      caw_raster_global   <- NULL
      mean_cell_area_global <- NULL
      if (apply_caw) {
            r_area <- terra::cellSize(pool[[1]], mask = FALSE, unit = "km")
            v <- terra::values(r_area)
            finite <- is.finite(v)
            if (!any(finite)) {
                  stop("Internal error: cellSize returned no finite values.",
                       call. = FALSE)
            }
            mean_area_physical <- mean(v[finite])
            if (mean_area_physical <= 0 || !is.finite(mean_area_physical)) {
                  stop("Internal error: non-positive mean cell area.",
                       call. = FALSE)
            }
            caw_raster_global <- r_area / mean_area_physical  # mean-1 normalized

            # Match the per-coord-type units used by .compute_cell_area_weights()
            # (and consumed by .compute_global_dmax()):
            #   - lonlat:    physical km^2
            #   - projected: planar projection-units^2 (prod(res(pool)))
            mean_cell_area_global <- if (identical(coord_type, "lonlat")) {
                  mean_area_physical
            } else {
                  r <- terra::res(pool)
                  as.numeric(r[1]) * as.numeric(r[2])
            }
      }

      # Get extents
      focal_ext <- terra::ext(x)
      ref_ext <- terra::ext(pool)

      kernel_width <- focal_ext[2] - focal_ext[1]
      kernel_height <- focal_ext[4] - focal_ext[3]
      kernel_aspect <- kernel_width / kernel_height

      # Find optimal tile grid
      grid <- optimize_tile_grid(n_tiles, kernel_aspect)
      n_x <- grid$n_x
      n_y <- grid$n_y

      tile_width <- kernel_width / n_x
      tile_height <- kernel_height / n_y

      # Warning for ineffective tiling
      check_tiling_effectiveness(
            tile_width, tile_height, max_geog,
            ref_ext, coord_type == "lonlat"
      )

      # Calculate pixel-based halo inflation (2 pixels)
      pixel_size <- max(terra::res(pool))
      pixel_buffer <- 2 * pixel_size

      # Main tiling loop
      results <- list()
      tile_files <- character(0)
      if (progress) {
            pb <- utils::txtProgressBar(min = 0, max = n_x * n_y, style = 3)
      }

      tile_idx <- 1
      for (i in 1:n_x) {
            for (j in 1:n_y) {
                  # Calculate tile bounds
                  tile_xmin <- focal_ext[1] + (i - 1) * tile_width
                  tile_xmax <- focal_ext[1] + i * tile_width
                  tile_ymin <- focal_ext[3] + (j - 1) * tile_height
                  tile_ymax <- focal_ext[3] + j * tile_height

                  tile_bbox <- terra::ext(tile_xmin, tile_xmax, tile_ymin, tile_ymax)

                  # Calculate halo bounds using buffer with pixel inflation
                  tile_vect <- terra::vect(tile_bbox, crs = terra::crs(x))

                  if (coord_type == "lonlat") {
                        # Convert both max_geog and pixel buffer to meters
                        pixel_buffer_m <- pixel_buffer * 111000  # degrees to meters (approximate)
                        total_buffer <- max_geog * 1000 + pixel_buffer_m
                  } else {
                        # Both in coordinate units
                        total_buffer <- max_geog + pixel_buffer
                  }

                  halo_bbox <- terra::ext(terra::buffer(tile_vect, total_buffer))

                  # Build argument list for function call
                  fun_args <- list(
                        x = terra::crop(x, tile_bbox),
                        pool = terra::crop(pool, halo_bbox),
                        geog = geog,
                        coord_type = coord_type,
                        progress = FALSE
                  )

                  # Forward the climate kernel if one was supplied.
                  if (!is.null(clim)) {
                        fun_args$clim <- clim
                  }

                  # Add y values if provided
                  if (!is.null(y)) {
                        fun_args$y <- terra::crop(y, halo_bbox)
                  }

                  # Add x_cov if provided
                  if (!is.null(x_cov)) {
                        fun_args$x_cov <- terra::crop(x_cov, tile_bbox)
                  }

                  # Add covariates if provided (cropped to halo, since they're
                  # pool-side). `covariates` was pulled out of `...` upstream.
                  if (!is.null(covariates)) {
                        fun_args$covariates <- terra::crop(covariates, halo_bbox)
                  }

                  # Add weight if provided (cropped to halo, pool-side).
                  if (!is.null(weight)) {
                        fun_args$weight <- terra::crop(weight, halo_bbox)
                  }

                  # Cell-area weights: use the globally-normalized raster's
                  # values within this tile's halo. We pass a numeric vector
                  # (extracted from the cropped raster, ordered by cell) so
                  # that build_analog_index uses these values as-is rather
                  # than recomputing per-tile (which would give per-tile
                  # local normalization).
                  if (apply_caw) {
                        caw_tile <- terra::crop(caw_raster_global, halo_bbox)
                        fun_args$cell_area_weight <- as.numeric(terra::values(caw_tile))

                        # Pass the global mean_cell_area through to the
                        # helper so any per-tile call with `normalize = TRUE`
                        # uses the same D_max denominator across tiles.
                        # Only inject when the helper accepts it -- either
                        # explicitly (named formal) or via `...`. Helpers
                        # that don't take either would error on the unused
                        # argument; in those cases (which can't currently
                        # use normalization anyway) we silently skip.
                        fun_formals <- names(formals(fun))
                        if ("mean_cell_area" %in% fun_formals ||
                            "..." %in% fun_formals) {
                              fun_args$mean_cell_area <- mean_cell_area_global
                        }
                  } else {
                        # Force off explicitly to override any default in helpers
                        fun_args$cell_area_weight <- FALSE
                  }

                  # Add any additional arguments from ... (covariates already
                  # extracted above)
                  fun_args <- c(fun_args, dot_args)

                  # Call function
                  result <- do.call(fun, fun_args)

                  # Remove analog_index column if present (tile-specific, not meaningful in merged result)
                  if (is.data.frame(result) && "analog_index" %in% names(result)) {
                        result$analog_index <- NULL
                  } else if (inherits(result, "SpatRaster") && "analog_index" %in% names(result)) {
                        result <- result[[names(result)[names(result) != "analog_index"]]]
                  }

                  # Handle disk-based output for SpatRasters
                  if (!is.null(output_file) && inherits(result, "SpatRaster")) {
                        tile_file <- tempfile(fileext = ".tif")
                        terra::writeRaster(result, tile_file)
                        tile_files <- c(tile_files, tile_file)
                  } else {
                        results[[tile_idx]] <- result
                  }

                  tile_idx <- tile_idx + 1

                  if (progress) utils::setTxtProgressBar(pb, tile_idx - 1)
            }
      }

      if (progress) close(pb)

      # Combine results
      if (!is.null(output_file) && length(tile_files) > 0) {
            # Disk-based merge for SpatRasters
            tile_rasters <- lapply(tile_files, terra::rast)
            result <- do.call(terra::merge, c(tile_rasters, filename = output_file))
            unlink(tile_files)  # Delete temporary files
      } else if (inherits(results[[1]], "SpatRaster")) {
            result <- do.call(terra::merge, results)
      } else {
            result <- do.call(rbind, results)
      }

      result
}


# Helper function to find optimal tile grid
optimize_tile_grid <- function(n_tiles, kernel_aspect) {
      # Try different factorizations
      max_dim <- ceiling(sqrt(n_tiles) * 2)
      candidates <- expand.grid(
            n_x = 1:max_dim,
            n_y = 1:max_dim
      )

      # Filter to reasonable range around target (allow 50% deviation)
      candidates <- candidates[
            candidates$n_x * candidates$n_y >= n_tiles * 2/3 &
                  candidates$n_x * candidates$n_y <= n_tiles * 3/2,
      ]

      if (nrow(candidates) == 0) {
            # Fallback: closest factorization
            candidates <- expand.grid(n_x = 1:max_dim, n_y = 1:max_dim)
      }

      # Score each candidate
      candidates$total_tiles <- candidates$n_x * candidates$n_y

      # Tile aspect ratio (we want this close to 1)
      candidates$tile_aspect <- (candidates$n_x / kernel_aspect) / candidates$n_y
      candidates$aspect_score <- abs(log(candidates$tile_aspect))

      # Deviation from target count (normalized by target)
      candidates$count_score <- abs(candidates$total_tiles - n_tiles) / n_tiles

      # Combined score: weight aspect more heavily than count
      # Prioritize square tiles, but don't go too far from target count
      candidates$combined_score <- candidates$aspect_score + 0.5 * candidates$count_score

      # Pick best
      best <- candidates[order(candidates$combined_score), ][1, ]

      list(n_x = best$n_x, n_y = best$n_y)
}

# Helper function to check tiling effectiveness
check_tiling_effectiveness <- function(
            tile_width, tile_height, max_geog, ref_ext, is_lonlat
) {
      ref_area <- (ref_ext[2] - ref_ext[1]) * (ref_ext[4] - ref_ext[3])

      # Estimate halo dimensions
      if (is_lonlat) {
            # Rough approximation: 1 degree ≈ 111 km at equator
            max_geog_deg <- max_geog / 111
            halo_width <- tile_width + 2 * max_geog_deg
            halo_height <- tile_height + 2 * max_geog_deg
      } else {
            halo_width <- tile_width + 2 * max_geog
            halo_height <- tile_height + 2 * max_geog
      }

      halo_area <- halo_width * halo_height

      # Warn if halo area is large relative to ref area
      # This means each tile loads most of the reference data, defeating the memory benefit of tiling
      if (halo_area > 0.5 * ref_area) {
            warning(
                  "The combination of tile size and max_geog is large relative to reference kernel. ",
                  "Each tile will load >50% of the reference data, limiting memory benefits. ",
                  "Consider using a smaller max_geog and/or more tiles ",
                  "(or running without tiling if memory allows).",
                  immediate. = TRUE
            )
      }
}
