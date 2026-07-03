#' Build Analog Index
#'
#' Pre-builds a reusable lattice index from reference climate data. The index
#' can be queried multiple times with different focal points and parameters,
#' avoiding the need to rebuild the lattice for each query.
#'
#' @param pool The reference dataset to search for analogs. Should be a
#'  matrix/data.frame with columns x, y, and climate variables, or a SpatRaster
#'  with climate variable layers.
#'
#' @param coord_type Coordinate system type:
#'   \itemize{
#'     \item \code{"auto"} (default): Automatically detect from coordinate ranges.
#'     \item \code{"lonlat"}: Unprojected lon/lat coordinates (uses great-circle distance;
#'       assumes \code{max_geog} is in km).
#'     \item \code{"projected"}: Projected XY coordinates (uses planar distance;
#'       assumes \code{max_geog} is in projection units).
#'   }
#'
#' @param clim_res_adj,geog_res_adj Non-negative scalars adjusting the lattice
#'   resolution of the climate and geographic variable families, each as a
#'   multiplier on a data-dependent default. The default allocation targets an
#'   average of roughly 50 pool points per occupied bin and splits that budget
#'   between the two families by their effective dimensionality; `clim_res_adj =
#'   geog_res_adj = 1` (the default) uses that allocation. Larger values give a
#'   finer lattice for that family (more, smaller bins), smaller values a
#'   coarser one, and `0` deactivates the family entirely (one bin per axis),
#'   which is appropriate when a query does not constrain that family. Because
#'   the default scales with pool size, these adjustments travel across datasets
#'   of different sizes. Within each family, bins are allocated in proportion to
#'   each axis' data range so a search radius spans a comparable number of bins
#'   on every axis. [analog_search()] sets these automatically from the query's
#'   constraints (deactivating an unconstrained family). Ignored if `pool` is an
#'   `analog_index`.
#'
#' @param downsample Optional downsampling rate (0-1) indicating the proportion of
#'   points in `pool` to retain. Downsampling reduces memory use and improves
#'   query speed at the cost of some precision; adaptive stratified sampling is
#'   used to minimize loss of precision. The default is 1.0 (no downsampling).
#'   See Details for more info.
#'
#' @param seed Optional random seed for reproducible downsampling. If `NULL`
#'   (default), uses current R random state.
#'
#' @param cell_area_weight Controls cell-area weighting for raster pools.
#'   One of:
#'   \itemize{
#'     \item `"auto"` (default): Compute cell-area weights when `pool` is a
#'       SpatRaster, and skip them otherwise. This corrects aggregation
#'       statistics for non-uniform cell areas (e.g. lonlat grids where cell
#'       area shrinks toward the poles, or projected grids on
#'       non-equal-area projections).
#'     \item `TRUE`: Force cell-area weighting on. Errors if `pool` is not a
#'       SpatRaster.
#'     \item `FALSE`: Force cell-area weighting off; treat all pool points
#'       as having equal weight.
#'     \item A numeric vector of length `nrow(pool)`: Use these
#'       caller-supplied weights as-is, without any further normalization.
#'       This is intended for advanced workflows like `tiled_analog_search()`
#'       that need to maintain a globally consistent normalization across
#'       multiple per-tile index builds; most users should use one of the
#'       three options above.
#'   }
#'   When `"auto"` or `TRUE` triggers computation, weights are computed via
#'   `terra::cellSize()` and normalized to mean 1 over finite values, so
#'   absolute magnitudes of stats like `sum_weights` remain comparable to
#'   the unweighted case. The weights are stored on the returned index and
#'   used during all subsequent queries.
#'
#' @param mean_cell_area Optional scalar mean cell area (in km^2) to attach
#'   to the index, overriding any value auto-computed from the raster pool.
#'   Intended for internal use by `tiled_analog_search()` to propagate a
#'   globally-consistent mean area across per-tile index builds (so that
#'   `analog_density(normalize = TRUE)` produces consistent values across
#'   tiles). Most users should leave this `NULL`.
#'
#' @return An S3 object of class \code{"analog_index"} containing:
#'   \itemize{
#'     \item The compiled lattice index (internal C++ structure)
#'     \item Reference data
#'     \item Metadata: coordinate type, dimensions, ranges, resolution
#'     \item Diagnostics: bin counts, occupancy statistics, and downsampling info
#'   }
#'
#' @details
#' The lattice index is a multidimensional grid of bins, built over both geographic
#' and climate dimensions. This structure enables efficient analog searches by
#' first filtering and sorting bins of similar points before computing exact
#' results. For lon/lat coordinates, the index uses ECEF (Earth-Centered Earth-Fixed)
#' space internally for optimal performance.
#'
#' Index resolution is controlled per family by `geog_res_adj` and
#' `clim_res_adj`, each a multiplier on a pool-size-dependent default (targeting
#' ~50 points per bin, split between the families by effective dimensionality).
#' `1` uses the default for that family, larger values are finer, smaller are
#' coarser, and `0` deactivates a family (one bin per axis). The optimal values
#' depend on your data and query patterns; [analog_search()] sets them
#' automatically from the query's constraints, or accept the defaults of `1`,
#' which work well for many applications.
#'
#' ## Downsampling
#'
#' For very large datasets, downsampling can significantly improve memory usage
#' and query speed, at the cost of some precision. The `downsample` parameter controls
#' the target fraction of the data points in `pool` that are retained in the index.
#' Downsampling uses an adaptive stratified approach: densely-packed bins are thinned more
#' aggressively while sparse bins are preserved, which helps reduce imprecision
#' in sparse regions compared to fully random sampling. Note: The actual rate may be
#' higher than requested if maintaining at least one point per occupied bin requires
#' it (common with sparse data or fine-grained binning); check `index$downsample_actual`.
#'
#' Each remaining analog in the downsampled pool gets a `sample_weight` indicating
#' the number of points it represents in the original pool; this weight is the inverse
#' of the sampling rate in the analog's index bin. For pair queries (`stat = "none"`),
#' results include each analog's `sample_weight`. For aggregation stats (count, sum,
#' mean, etc.), sampling weights are used internally to automatically correct for
#' the downsampling bias.
#'
#' @examples
#' \dontrun{
#' # Build index with default settings
#' index <- build_analog_index(climate_data)
#'
#' # Build with finer geographic resolution than default
#' index <- build_analog_index(climate_data, geog_res_adj = 2)
#'
#' # Build with downsampling for large datasets
#' index <- build_analog_index(
#'   large_climate_data,
#'   geog_res_adj = 1, clim_res_adj = 1,
#'   downsample = 0.1,  # Reduce max bin size to 10%
#'   seed = 123         # Reproducible sampling
#' )
#'
#' # Query the index multiple times
#' v1 <- analog_velocity(sites1, pool = index, max_clim = 0.5)
#' v2 <- analog_velocity(sites2, pool = index, max_clim = 0.3)
#' a1 <- analog_availability(sites3, pool = index, max_clim = 0.5, max_geog = 100)
#' }
#'
#' @export
build_analog_index <- function(pool,
                               coord_type = c("auto", "lonlat", "projected"),
                               geog_res_adj = 1,
                               clim_res_adj = 1,
                               downsample = 1.0,
                               seed = NULL,
                               cell_area_weight = "auto",
                               mean_cell_area = NULL) {

      # Validate inputs
      coord_type <- match.arg(coord_type)

      if (!is.numeric(geog_res_adj) || length(geog_res_adj) != 1L ||
          !is.finite(geog_res_adj) || geog_res_adj < 0) {
            stop("geog_res_adj must be a single non-negative number (0 = deactivate)")
      }

      if (!is.numeric(clim_res_adj) || length(clim_res_adj) != 1L ||
          !is.finite(clim_res_adj) || clim_res_adj < 0) {
            stop("clim_res_adj must be a single non-negative number (0 = deactivate)")
      }

      if (!is.numeric(downsample) || length(downsample) != 1L ||
          downsample <= 0 || downsample > 1) {
            stop("downsample must be a number between 0 and 1 (exclusive of 0)")
      }

      # Validate cell_area_weight: "auto", TRUE, FALSE, or a numeric vector
      # of length n_pool. Numeric vectors are accepted as-is (caller-supplied
      # already-normalized weights) and validated for length once we know
      # n_pool below.
      cell_area_weight_is_vec <- is.numeric(cell_area_weight) &&
            is.null(dim(cell_area_weight))
      if (!(identical(cell_area_weight, "auto") ||
            isTRUE(cell_area_weight) ||
            isFALSE(cell_area_weight) ||
            cell_area_weight_is_vec)) {
            stop('`cell_area_weight` must be "auto", TRUE, FALSE, or a numeric ',
                 "vector of length nrow(pool).",
                 call. = FALSE)
      }

      # Validate optional mean_cell_area override.
      if (!is.null(mean_cell_area)) {
            if (!is.numeric(mean_cell_area) || length(mean_cell_area) != 1L ||
                !is.finite(mean_cell_area) || mean_cell_area <= 0) {
                  stop("`mean_cell_area` must be a single positive finite number, or NULL.",
                       call. = FALSE)
            }
            mean_cell_area <- as.numeric(mean_cell_area)
      }

      # Resolve cell_area_weight to one of three states:
      #   - cell_area_weight_vec is NULL (off)
      #   - cell_area_weight_vec is a numeric vector of length n_pool (on)
      pool_is_raster <- inherits(pool, "SpatRaster")
      if (identical(cell_area_weight, "auto")) {
            compute_from_raster <- pool_is_raster
            user_supplied_vec   <- FALSE
      } else if (isTRUE(cell_area_weight)) {
            if (!pool_is_raster) {
                  stop("`cell_area_weight = TRUE` requires `pool` to be a SpatRaster.",
                       call. = FALSE)
            }
            compute_from_raster <- TRUE
            user_supplied_vec   <- FALSE
      } else if (isFALSE(cell_area_weight)) {
            compute_from_raster <- FALSE
            user_supplied_vec   <- FALSE
      } else {
            # Numeric vector path. We trust whatever the caller hands us
            # (no global re-normalization) — the helper that supplies it is
            # expected to have done its own normalization. Vectors are
            # accepted regardless of pool type since the caller may have a
            # legitimate use case (e.g. tiled_analog_search slicing a global
            # vector for a haloed matrix-form pool).
            compute_from_raster <- FALSE
            user_supplied_vec   <- TRUE

            # Reject NA / negative / non-finite values up front.
            if (any(!is.finite(cell_area_weight) | cell_area_weight < 0)) {
                  stop("`cell_area_weight` numeric vector must contain only ",
                       "non-negative finite values.", call. = FALSE)
            }
      }

      # Handle seed
      if (!is.null(seed)) {
            if (!is.numeric(seed) || length(seed) != 1L) {
                  stop("seed must be a single integer or NULL")
            }
            seed <- as.integer(seed)
      } else {
            # Use R's current RNG state to generate a seed
            seed <- as.integer(sample.int(.Machine$integer.max, 1))
      }

      # Resolve to a single numeric vector (or NULL) using the flags from above.
      # When auto-computing from a raster, we also capture the mean cell area
      # used in the normalization. Caller-supplied numeric vectors are taken
      # as-is; in that case the caller is responsible for also supplying
      # `mean_cell_area` if downstream `analog_density(normalize = TRUE)` is
      # to be supported on this index.
      cell_area_weight_vec <- NULL
      computed_mean_area   <- NULL
      if (compute_from_raster) {
            caw <- .compute_cell_area_weights(pool)
            cell_area_weight_vec <- caw$weights
            computed_mean_area   <- caw$mean_area
      } else if (user_supplied_vec) {
            cell_area_weight_vec <- as.numeric(cell_area_weight)
      }

      # Resolve final mean_cell_area to store on the index.
      # Priority: explicit `mean_cell_area` argument > value computed from
      # the raster (when auto-computing). For pure user-supplied numeric
      # vectors with no `mean_cell_area` argument, this stays NULL — the
      # caller (e.g. tiled_analog_search) must pass it explicitly to
      # support density normalization.
      mean_cell_area_final <- if (!is.null(mean_cell_area)) {
            mean_cell_area
      } else {
            computed_mean_area  # may be NULL
      }

      # Normalize pool to standard matrix format. NA pool rows (e.g. ocean
      # cells in a global raster) are stripped here; `row_map` records the
      # original-pool row index of each kept row so we can translate C++
      # analog indices back to the user's pool indexing in pairs mode.
      ref_mm <- .format_data(pool)
      pool_row_map <- attr(ref_mm, "row_map")  # NULL if no rows stripped

      # Original (pre-strip) pool size, for length-checking user-supplied
      # weight vectors against the input the user actually passed.
      n_pool_original <- if (is.null(pool_row_map)) {
            nrow(ref_mm)
      } else {
            # row_map indices are into the original pool; max index = original size
            # only if the last row was retained, so use the dataset's true row count.
            # We recover that from the input.
            if (inherits(pool, "SpatRaster")) {
                  terra::ncell(pool)
            } else {
                  nrow(pool)
            }
      }

      # Detect coordinate system if auto. We pass the original `pool` so
      # the detector can consult CRS metadata (when pool is a SpatRaster);
      # the xy matrix is a fallback for inputs that lack metadata.
      if (coord_type == "auto") {
            coord_type <- .detect_geo(ref_mm[, 1:2, drop = FALSE], pool)
      }

      # Validate cell-area weight vector length against the *original* pool
      # size (the user supplied it for the pool they passed in), then strip
      # to align with the NA-filtered ref_mm if any rows were dropped.
      if (!is.null(cell_area_weight_vec)) {
            if (length(cell_area_weight_vec) != n_pool_original) {
                  stop(sprintf(
                        "`cell_area_weight` length (%d) does not match pool size (%d).",
                        length(cell_area_weight_vec), n_pool_original
                  ), call. = FALSE)
            }
            if (!is.null(pool_row_map)) {
                  cell_area_weight_vec <- cell_area_weight_vec[pool_row_map]
            }
      }

      # Translate the resolution knob into an absolute total bin-count target.
      #
      # The single source of truth for the occupancy-based default lives here.
      # We target an average of TARGET_OCC usable pool points per occupied bin,
      # giving a default TOTAL bin count B = n_pool_used / TARGET_OCC that scales
      # with pool size. B is split between the geographic and climate families
      # in proportion to their EFFECTIVE dimensionality (so at the defaults the
      # behavior matches a balanced allocation), then each family's share is
      # scaled by its per-family resolution adjustment:
      #
      #   geo_default  = B ^ (n_geo_eff / (n_geo_eff + n_clim))
      #   clim_default = B ^ (n_clim    / (n_geo_eff + n_clim))
      #   geo_target   = geog_res_adj  * geo_default
      #   clim_target  = clim_res_adj * clim_default
      #
      # n_geo_eff is the EFFECTIVE geographic dimensionality: for lon/lat the geo
      # axes are stored as 3 ECEF coords but lie on a 2D manifold, so we use 2;
      # projected geo is genuinely 2. Both cases are 2 here. A *_res_adj of 0
      # (or a target <= 1) deactivates that family (1 bin per axis).
      TARGET_OCC <- 50
      n_pool_used <- nrow(ref_mm)
      n_clim <- ncol(ref_mm) - 2L
      n_geo_eff <- 2  # ECEF (lonlat) manifold dim OR projected geo dim

      B <- n_pool_used / TARGET_OCC
      if (!is.finite(B) || B < 1) B <- 1

      geo_frac  <- n_geo_eff / (n_geo_eff + n_clim)
      clim_frac <- n_clim    / (n_geo_eff + n_clim)
      geo_default  <- B ^ geo_frac
      clim_default <- B ^ clim_frac

      # res_adj = 0 deactivates (target 0 -> C++ leaves family at 1 bin/axis).
      geo_target  <- if (geog_res_adj  == 0) 0 else geog_res_adj  * geo_default
      clim_target <- if (clim_res_adj == 0) 0 else clim_res_adj * clim_default
      geo_target  <- as.numeric(geo_target)
      clim_target <- as.numeric(clim_target)

      # Call C++ to build index
      index_cpp <- build_analog_index_cpp(
            ref_mm = ref_mm,
            coord_type = coord_type,
            geo_target = geo_target,
            clim_target = clim_target,
            downsample = downsample,
            seed = seed
      )

      # Construct S3 object
      obj <- structure(
            list(
                  # C++ components
                  lattice_xptr = index_cpp$lattice_xptr,
                  ecef_xptr = index_cpp$ecef_xptr,

                  # Reference data (needed for distance calculations).
                  # NA-filtered: rows where any coord/clim was NA have been
                  # stripped. See `pool_row_map` for original-pool indexing.
                  ref_data = ref_mm,

                  # Map from stripped (ref_data) row index -> original-pool
                  # row index. NULL when no rows were stripped (identity
                  # mapping). Used by query_analog_index() to translate C++
                  # analog indices back to user-pool indexing in pair mode.
                  pool_row_map = pool_row_map,

                  # Metadata. n_pool reports the original pool size as the
                  # user passed it; n_pool_used reports the count of usable
                  # (non-NA) rows actually held in ref_data and indexed by
                  # the lattice. These differ when the input contained NAs.
                  coord_type = coord_type,
                  n_pool = n_pool_original,
                  n_pool_used = index_cpp$n_pool,
                  n_clim = index_cpp$n_clim,

                  # Resolution parameters. geog_res_adj/clim_res_adj are the
                  # per-family knobs the user set (1 = default, 0 = deactivated);
                  # geo_target/clim_target are the realized per-family bin-count
                  # targets passed to C++; bins_per_axis is the realized per-axis
                  # split (geo axes first, then climate).
                  geog_res_adj = geog_res_adj,
                  clim_res_adj = clim_res_adj,
                  geo_target = index_cpp$geo_target,
                  clim_target = index_cpp$clim_target,
                  bins_per_axis = index_cpp$bins_per_axis,

                  # Coordinate ranges
                  coord_mins = index_cpp$coord_mins,
                  coord_maxs = index_cpp$coord_maxs,

                  # Climate ranges
                  clim_mins = index_cpp$clim_mins,
                  clim_maxs = index_cpp$clim_maxs,

                  # Internal flags
                  use_ecef = index_cpp$use_ecef,

                  # Cell-area weighting (NULL if not applied; otherwise
                  # mean-1-normalized numeric vector of length n_pool_used,
                  # aligned to ref_data rows after any NA stripping).
                  cell_area_weight = cell_area_weight_vec,

                  # Mean cell area (km^2) used to normalize the cell-area
                  # weights. Stored regardless of which pool form was used
                  # (raster auto-compute, user-supplied vector with explicit
                  # `mean_cell_area`, etc.). NULL when neither path supplied
                  # a value — in that case `analog_density(normalize = TRUE)`
                  # cannot be used against this index.
                  mean_cell_area = mean_cell_area_final,

                  # Diagnostics
                  total_bins = index_cpp$total_bins,
                  n_bins_nonempty = index_cpp$n_bins_nonempty,
                  min_bin_occupancy = index_cpp$min_bin_occupancy,
                  max_bin_occupancy = index_cpp$max_bin_occupancy,

                  # Downsampling info
                  downsample_target = downsample,
                  downsample_actual = index_cpp$downsample_actual,
                  seed = seed
            ),
            class = "analog_index"
      )

      return(obj)
}


#' Print Method for Analog Index
#'
#' @param x An \code{analog_index} object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.analog_index <- function(x, ...) {
      cat("Analog Index\n")
      cat("============\n\n")

      cat("Reference data:\n")
      # n_pool_used may be absent on legacy index objects; fall back to n_pool.
      n_used <- x$n_pool_used %||% x$n_pool
      if (!is.null(x$n_pool_used) && x$n_pool_used != x$n_pool) {
            cat(sprintf("  %d locations (%d usable after NA filtering)\n",
                        x$n_pool, n_used))
      } else {
            cat(sprintf("  %d locations\n", x$n_pool))
      }
      cat(sprintf("  %d climate variables\n", x$n_clim))
      cat(sprintf("  Coordinate type: %s\n", x$coord_type))

      cat("\nCoordinate ranges:\n")
      cat(sprintf("  X: [%.3f, %.3f]\n", x$coord_mins[1], x$coord_maxs[1]))
      cat(sprintf("  Y: [%.3f, %.3f]\n", x$coord_mins[2], x$coord_maxs[2]))

      cat("\nClimate ranges:\n")
      for (i in seq_along(x$clim_mins)) {
            cat(sprintf("  Variable %d: [%.3f, %.3f]\n",
                        i, x$clim_mins[i], x$clim_maxs[i]))
      }

      cat("\nIndex structure:\n")
      cat(sprintf("  geog_res_adj: %g   clim_res_adj: %g\n",
                  x$geog_res_adj %||% NA_real_, x$clim_res_adj %||% NA_real_))
      if (!is.null(x$bins_per_axis)) {
            n_ax  <- length(x$bins_per_axis)
            n_clim_ax <- x$n_clim %||% 0L
            n_geo_ax  <- max(0L, n_ax - n_clim_ax)
            geo_b  <- if (n_geo_ax > 0) x$bins_per_axis[seq_len(n_geo_ax)] else integer(0)
            clim_b <- if (n_clim_ax > 0) x$bins_per_axis[(n_geo_ax + 1L):n_ax] else integer(0)
            cat(sprintf("  Bins per axis: geo {%s}  climate {%s}\n",
                        paste(geo_b, collapse = ", "),
                        paste(clim_b, collapse = ", ")))
      }
      cat(sprintf("  Total bins: %s\n", format(x$total_bins, big.mark = ",")))
      cat(sprintf("  Non-empty bins: %s (%.1f%%)\n",
                  format(x$n_bins_nonempty, big.mark = ","),
                  100 * x$n_bins_nonempty / x$total_bins))

      if (x$n_bins_nonempty > 0) {
            avg_occ <- n_used / x$n_bins_nonempty
            cat(sprintf("  Bin occupancy: %d-%d (avg: %.1f) points per bin\n",
                        x$min_bin_occupancy, x$max_bin_occupancy, avg_occ))
      }

      # Show downsampling info if applied
      if (!is.null(x$downsample_actual) && x$downsample_actual < 1.0) {
            cat(sprintf("\nDownsampling:\n"))
            cat(sprintf("  Target rate: %.1f%%\n", x$downsample_target * 100))
            cat(sprintf("  Actual rate: %.1f%% of points retained\n",
                        x$downsample_actual * 100))
            cat(sprintf("  Seed: %d\n", x$seed))
      }

      if (x$use_ecef) {
            cat("\n  Note: Using ECEF coordinates internally for lon/lat data\n")
      }

      if (!is.null(x$cell_area_weight)) {
            cat("\n  Note: Cell-area weights applied to pool points\n")
            if (!is.null(x$mean_cell_area)) {
                  cat(sprintf("        (mean cell area: %.3f km^2)\n",
                              x$mean_cell_area))
            }
      }

      invisible(x)
}


#' Check if Object is an Analog Index
#'
#' @param x Object to test
#' @return Logical indicating whether \code{x} is an \code{analog_index}
#'
#' @keywords internal
is_analog_index <- function(x) {
      inherits(x, "analog_index")
}


#' Validate Analog Index
#'
#' Internal function to validate that an analog index is well-formed and
#' compatible with query data.
#'
#' @param index An \code{analog_index} object
#' @param query_data Optional matrix of query points to validate against
#' @param validate_ranges Logical; if TRUE and query_data provided, check
#'   that query points fall within index bounds
#'
#' @return Invisible TRUE if valid, otherwise throws an error
#' @keywords internal
.validate_analog_index <- function(index,
                                   query_data = NULL,
                                   validate_ranges = FALSE) {

      # Check class
      if (!is_analog_index(index)) {
            stop("Object is not an analog_index")
      }

      # Check required components
      required <- c("lattice_xptr", "ref_data", "coord_type",
                    "n_pool", "n_clim", "geo_target", "clim_target")
      missing <- setdiff(required, names(index))
      if (length(missing) > 0) {
            stop("Invalid analog_index: missing components: ",
                 paste(missing, collapse = ", "))
      }

      # Check Xptr validity
      if (is.null(index$lattice_xptr)) {
            stop("Invalid analog_index: lattice pointer is NULL")
      }

      # Validate against query data if provided
      if (!is.null(query_data)) {
            query_ncol <- ncol(query_data)
            expected_ncol <- 2 + index$n_clim

            if (query_ncol != expected_ncol) {
                  stop(sprintf(
                        "Query data has %d columns but index expects %d (2 coords + %d climate)",
                        query_ncol, expected_ncol, index$n_clim
                  ))
            }

            # Optionally check if query points are within index bounds
            if (validate_ranges) {
                  query_x <- query_data[, 1]
                  query_y <- query_data[, 2]

                  x_out <- query_x < index$coord_mins[1] | query_x > index$coord_maxs[1]
                  y_out <- query_y < index$coord_mins[2] | query_y > index$coord_maxs[2]

                  n_out <- sum(x_out | y_out)
                  if (n_out > 0) {
                        warning(sprintf(
                              "%d query point(s) fall outside index coordinate bounds",
                              n_out
                        ))
                  }

                  # Check climate bounds (warnings only)
                  for (i in seq_len(index$n_clim)) {
                        clim_col <- query_data[, 2 + i]
                        clim_out <- clim_col < index$clim_mins[i] |
                              clim_col > index$clim_maxs[i]
                        n_clim_out <- sum(clim_out)
                        if (n_clim_out > 0) {
                              warning(sprintf(
                                    "%d query point(s) have climate variable %d outside index bounds",
                                    n_clim_out, i
                              ))
                        }
                  }
            }
      }

      invisible(TRUE)
}


#' Detect coordinate system from data ranges
#'
#' Resolves `coord_type = "auto"` to either `"lonlat"` or `"projected"`. When
#' the input has CRS metadata (SpatRaster, SpatVector), that's authoritative.
#' For raw matrix / data.frame inputs we fall back to a coordinate-magnitude
#' heuristic: if all xy values fit within `[-180, 180] x [-90, 90]`, treat as
#' lonlat; otherwise projected.
#'
#' The magnitude heuristic is unavoidably ambiguous (a small projected
#' region in meters can fit in lonlat bounds). Users with matrix-form
#' projected data in that range should pass `coord_type` explicitly.
#'
#' @keywords internal
.detect_geo <- function(xy, input = NULL) {
      # If we have CRS metadata, trust it. terra::is.lonlat() returns
      # TRUE/FALSE/NA; only fall back to the magnitude check on NA (no CRS,
      # or CRS recognized neither as lonlat nor as a known projection).
      if (!is.null(input) &&
          (inherits(input, "SpatRaster") || inherits(input, "SpatVector"))) {
            if (requireNamespace("terra", quietly = TRUE)) {
                  # terra::is.lonlat() emits a "[is.lonlat] unknown crs"
                  # warning when CRS is missing, then returns NA. We treat
                  # NA as a fallback signal, so suppress the warning to
                  # avoid noise on legitimate no-CRS inputs.
                  is_ll <- suppressWarnings(terra::is.lonlat(input))
                  if (isTRUE(is_ll))  return("lonlat")
                  if (isFALSE(is_ll)) return("projected")
                  # is.na(is_ll): unknown CRS; fall through to magnitude check.
            }
      }

      lon_rng <- range(xy[, 1], na.rm = TRUE)
      lat_rng <- range(xy[, 2], na.rm = TRUE)

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
