#' Tiled analog search for memory-constrained queries
#'
#' Performs analog searches on large datasets by dividing the focal region
#' into tiles and processing each tile separately. This reduces memory usage,
#' and with `method = "fft"` it also keeps each tile's FFT grid small. Works
#' with any `analog_*()` function.
#'
#' @param x Focal locations: a SpatRaster with environmental variable layers,
#'   or a matrix/data.frame with columns x, y, and environmental variables.
#' @param pool Reference locations (potential analog pool): a SpatRaster, or a
#'   matrix/data.frame with columns x, y, and environmental variables. An
#'   `analog_index` is not accepted, since each tile needs its own pool subset.
#'   Point tables carry no CRS, so `x` and `pool` must share a coordinate
#'   system whenever either is a point table (not checked), and `coord_type`
#'   should usually be given explicitly.
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
#' @param env Optional [kernel()] object giving the environmental search treatment,
#'   forwarded to `fun`. `NULL` (default) means no environmental kernel is passed.
#' @param y Optional values to aggregate across analogs. A SpatRaster matching
#'   `pool`'s CRS and extent when `pool` is a raster; a vector, factor, or
#'   matrix/data.frame with one row per row of `pool` when `pool` is a point
#'   table.
#' @param x_cov Optional focal covariances. A SpatRaster matching `x`'s CRS
#'   and extent when `x` is a raster; a matrix/data.frame with one row per row
#'   of `x` when `x` is a point table.
#' @param weight Optional per-pool-site weights. A single-layer SpatRaster
#'   matching `pool`'s CRS and extent when `pool` is a raster; a numeric
#'   vector with one value per row of `pool` when `pool` is a point table.
#'   See [analog_search()] for details.
#' @param cell_area_weight Controls cell-area weighting. One of `"auto"`
#'   (default; on for raster pools, off for point pools), `TRUE`, or `FALSE`.
#'   `TRUE` requires a raster pool. Unlike a non-tiled query, where weights
#'   are normalized to mean 1 over whatever pool is passed, here weights are
#'   normalized to mean 1 over the *full pool raster* once at the outer level,
#'   and the resulting raster is sliced per tile. This keeps absolute
#'   magnitudes of weighted stats (e.g. `sum_weights`) consistent across
#'   tiles.
#' @param coord_type One of `"auto"` (default), `"lonlat"`, or `"projected"`;
#'   see [analog_search()]. `"auto"` reads the CRS of whichever of `x` and
#'   `pool` is a raster, and falls back to a coordinate-range check only when
#'   both are point tables. That check reads small projected coordinates (any
#'   within +/-180, +/-90) as lon/lat, so set this explicitly for point
#'   tables. Required to be `"projected"` for geographic kernels with
#'   `method = "fft"`.
#' @param ... Additional arguments passed to fun (e.g. `select`, `stat`, `k`).
#'   May include `covariates`, aligned with `pool` in the same way as `y`
#'   (cropped or subset per tile alongside `y`). May also include `normalize`
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
#'   is specified, returns a disk-backed SpatRaster. When `x` is a point
#'   table, rows come back in the order of `x`, and the `index` column refers
#'   to rows of `x`, exactly as in a non-tiled call.
#'
#' @details
#' Tiled analog searches work by splitting x into a number of smaller tiles and
#' calling the requested analog function on each tile, using an analog pool
#' that is the size of the tile buffered by the geographic radius (`geog`'s
#' `max`). This buffer is necessary
#' for correctness but increases compute time, particularly if the radius is
#' large. Raster results for each tile are merged into a single results raster
#' once all tiles have processed (via temporary files when `output_file` is
#' set).
#'
#' Raster inputs are cropped to each tile and halo; point-table inputs are
#' subset by row. Each focal point belongs to exactly one tile, and focal rows
#' with missing coordinates are passed with the first tile, so they come back
#' as `NA` rows just as they would without tiling.
#'
#' The function requires `geog` to have a finite `max`, as tiling is only beneficial when
#' geographic distance constraints limit the reference pool size for each
#' focal point. The function will warn if the radius is so large that tiling
#' provides minimal memory benefit.
#'
#' If index_res is specified in ..., all tiles will use the same lattice
#' resolution. If index_res is not specified, each tile will independently
#' auto-tune its lattice resolution based on local data characteristics. This
#' adaptive behavior is generally fine and can even be beneficial when environmental
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
            env = NULL,
            y = NULL,
            x_cov = NULL,
            weight = NULL,
            cell_area_weight = "auto",
            coord_type = c("auto", "lonlat", "projected"),
            ...,
            output_file = NULL,
            progress = TRUE
) {
      # Validation
      if (!requireNamespace("terra", quietly = TRUE)) {
            stop("package `terra` is required for tiled_analog_search")
      }
      coord_type <- match.arg(coord_type)
      x_pts    <- .is_point_table(x)
      pool_pts <- .is_point_table(pool)
      if (!x_pts && !inherits(x, "SpatRaster")) {
            stop("x must be a SpatRaster, matrix, or data.frame for tiled_analog_search")
      }
      if (!pool_pts && !inherits(pool, "SpatRaster")) {
            stop("pool must be a SpatRaster, matrix, or data.frame for tiled_analog_search")
      }

      # Check CRS compatibility (point tables carry no CRS to compare)
      if (!x_pts && !pool_pts && !terra::same.crs(x, pool)) {
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
      if (!is.null(env) && !inherits(env, "analog_kernel")) {
            stop("`env` must be a kernel() object or NULL.", call. = FALSE)
      }

      # Pull `covariates` out of `...` so we can subset it per tile (it's
      # pool-side, like y). Without this, each tile's pool would be a halo
      # subset while `covariates` stayed full-extent.
      dot_args   <- list(...)
      covariates <- dot_args$covariates
      dot_args$covariates <- NULL

      # Companion inputs must be aligned with the input they belong to: rasters
      # with a raster (same CRS and extent), row-aligned tables with a point
      # table. y, weight, and covariates are pool-side; x_cov is focal-side.
      .check_aligned(y, pool, "y", "pool")
      .check_aligned(x_cov, x, "x_cov", "x")
      .check_aligned(weight, pool, "weight", "pool")
      .check_aligned(covariates, pool, "covariates", "pool")
      if (!is.null(weight) && !pool_pts && terra::nlyr(weight) != 1L) {
            stop("weight must have exactly one layer")
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
      if (isTRUE(cell_area_weight) && pool_pts) {
            stop("`cell_area_weight = TRUE` requires a raster pool.", call. = FALSE)
      }
      apply_caw <- if (identical(cell_area_weight, "auto")) {
            !pool_pts
      } else {
            isTRUE(cell_area_weight)
      }

      # Determine coord_type once, here in the wrapper, and forward it
      # explicitly to every per-tile call, so all tiles agree even if
      # individual tiles' coordinate ranges would have triggered a different
      # fallback magnitude check than the full data would. .detect_geo()
      # respects CRS metadata first (so we hand it a raster input when there
      # is one) and falls back to coordinate magnitude only when CRS is
      # unknown. The extent corners give that check the full data range
      # without materializing every cell.
      if (coord_type == "auto") {
            detect_from <- if (!x_pts) x else pool
            ex <- .spatial_ext(detect_from)
            xy_for_detect <- matrix(c(ex[1], ex[2], ex[3], ex[4]), nrow = 2, ncol = 2,
                                    dimnames = list(NULL, c("x", "y")))
            coord_type <- .detect_geo(xy_for_detect, input = detect_from)
            if (!coord_type %in% c("lonlat", "projected")) {
                  stop("Could not determine coord_type from `x` or `pool`. ",
                       "Set `coord_type` explicitly, or give a raster input a ",
                       "recognized CRS.", call. = FALSE)
            }
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

      # Get extents. A point table can span zero width or height (e.g. a
      # single transect); floor the span so tile geometry stays finite.
      focal_ext <- .spatial_ext(x)
      ref_ext <- .spatial_ext(pool)

      kernel_width <- max(focal_ext[2] - focal_ext[1], sqrt(.Machine$double.eps))
      kernel_height <- max(focal_ext[4] - focal_ext[3], sqrt(.Machine$double.eps))
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

      # Halo inflation of 2 pixels, guarding against crop() snapping a raster
      # pool's halo inward. Point pools are subset exactly, so need none.
      pixel_buffer <- if (pool_pts) 0 else 2 * max(terra::res(pool))

      # Assign each focal point to exactly one tile (half-open intervals,
      # with the far edges closed). Rows with missing coordinates ride along
      # with tile 1 and come back as NA rows, as in a non-tiled call.
      if (x_pts) {
            ti <- pmin(floor((x[, 1] - focal_ext[1]) / tile_width) + 1, n_x)
            tj <- pmin(floor((x[, 2] - focal_ext[3]) / tile_height) + 1, n_y)
            x_tile <- (ti - 1) * n_y + tj
            x_tile[is.na(x_tile)] <- 1
      }

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

                  # Point focals: skip tiles with none (there is nothing to
                  # query, and tile 1 always exists to carry NA-coordinate rows)
                  if (x_pts) {
                        x_rows <- which(x_tile == tile_idx)
                        if (length(x_rows) == 0) {
                              tile_idx <- tile_idx + 1
                              if (progress) utils::setTxtProgressBar(pb, tile_idx - 1)
                              next
                        }
                  }

                  # Calculate halo bounds using buffer with pixel inflation.
                  # Point tables carry no CRS, so take it from coord_type; the
                  # buffer only needs to know whether distances are geodesic.
                  tile_crs <- if (!x_pts) terra::crs(x) else if (!pool_pts) terra::crs(pool) else
                        if (coord_type == "lonlat") "EPSG:4326" else ""
                  tile_vect <- terra::vect(tile_bbox, crs = tile_crs)

                  if (coord_type == "lonlat") {
                        # Convert both max_geog and pixel buffer to meters
                        pixel_buffer_m <- pixel_buffer * 111000  # degrees to meters (approximate)
                        total_buffer <- max_geog * 1000 + pixel_buffer_m
                  } else {
                        # Both in coordinate units
                        total_buffer <- max_geog + pixel_buffer
                  }

                  halo_bbox <- terra::ext(terra::buffer(tile_vect, total_buffer))

                  # Focal-side inputs are cut to the tile, pool-side inputs to
                  # the halo: rasters by crop(), point tables by row.
                  focal_cut <- if (x_pts) function(v) .subset_rows(v, x_rows) else
                        function(v) terra::crop(v, tile_bbox)
                  pool_cut <- if (pool_pts) {
                        pool_rows <- which(.in_bbox(pool, halo_bbox))
                        function(v) .subset_rows(v, pool_rows)
                  } else {
                        function(v) terra::crop(v, halo_bbox)
                  }

                  # Build argument list for function call
                  fun_args <- list(
                        x = focal_cut(x),
                        pool = pool_cut(pool),
                        geog = geog,
                        coord_type = coord_type,
                        progress = FALSE
                  )

                  # Forward the environmental kernel if one was supplied.
                  if (!is.null(env)) {
                        fun_args$env <- env
                  }

                  # Companion inputs, cut alongside the input they belong to.
                  # `covariates` was pulled out of `...` upstream.
                  if (!is.null(y))          fun_args$y          <- pool_cut(y)
                  if (!is.null(x_cov))      fun_args$x_cov      <- focal_cut(x_cov)
                  if (!is.null(covariates)) fun_args$covariates <- pool_cut(covariates)
                  if (!is.null(weight))     fun_args$weight     <- pool_cut(weight)

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

                  # analog_index is tile-local. A point pool's can be mapped
                  # back to rows of the caller's `pool`; a raster pool's is a
                  # cell of the cropped halo and is dropped.
                  if (is.data.frame(result) && "analog_index" %in% names(result)) {
                        if (pool_pts) {
                              result$analog_index <- pool_rows[result$analog_index]
                        } else {
                              result$analog_index <- NULL
                        }
                  } else if (inherits(result, "SpatRaster") && "analog_index" %in% names(result)) {
                        result <- result[[names(result)[names(result) != "analog_index"]]]
                  }

                  # Point focals: map the tile-local focal index back to rows
                  # of the caller's `x`.
                  if (x_pts && is.data.frame(result) && "index" %in% names(result)) {
                        result$index <- x_rows[result$index]
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
            # Point focals: restore the caller's row order. order() is stable,
            # so multi-row (pair) results keep their within-focal order.
            if (x_pts && "index" %in% names(result)) {
                  result <- result[order(result$index), , drop = FALSE]
                  rownames(result) <- NULL
            }
      }

      result
}


# Spatial input helpers. Inputs are either SpatRasters or point tables
# (matrix/data.frame with x, y in the first two columns).

.is_point_table <- function(d) is.data.frame(d) || is.matrix(d)

.spatial_ext <- function(d) {
      if (inherits(d, "SpatRaster")) return(terra::ext(d))
      terra::ext(c(range(d[, 1], na.rm = TRUE), range(d[, 2], na.rm = TRUE)))
}

.in_bbox <- function(d, bbox) {
      b <- as.vector(bbox)
      !is.na(d[, 1]) & !is.na(d[, 2]) &
            d[, 1] >= b[1] & d[, 1] <= b[2] & d[, 2] >= b[3] & d[, 2] <= b[4]
}

.subset_rows <- function(v, i) if (is.null(dim(v))) v[i] else v[i, , drop = FALSE]

# Check that a companion input (y, weight, x_cov, covariates) is aligned with
# the input it belongs to. NULL passes.
.check_aligned <- function(obj, ref, name, ref_name) {
      if (is.null(obj)) return(invisible())
      if (inherits(ref, "SpatRaster")) {
            if (!inherits(obj, "SpatRaster")) {
                  stop(name, " must be a SpatRaster when ", ref_name,
                       " is a SpatRaster", call. = FALSE)
            }
            if (!terra::same.crs(obj, ref)) {
                  stop(name, " must have the same CRS as ", ref_name)
            }
            if (!all(terra::ext(obj) == terra::ext(ref))) {
                  stop(name, " must have the same extent as ", ref_name)
            }
      } else if (inherits(obj, "SpatRaster") || NROW(obj) != nrow(ref)) {
            stop(name, " must be a vector, matrix, or data.frame with one row ",
                 "per row of ", ref_name, " when ", ref_name, " is a point table",
                 call. = FALSE)
      }
      invisible()
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
