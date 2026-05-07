#' Generalized analog search
#'
#' Identifies locations in a reference dataset that are climatically similar
#' and/or geographically proximal to focal locations. Analog searches use a
#' two-stage approach: first selecting analogs based on specified criteria,
#' then optionally aggregating the results.
#'
#' @details
#' ## Parameter categories
#'
#' - *Data parameters*
#'   (`x`, `pool`, `x_cov`, `y`, `covariates`, `weight`, `coord_type`)
#'   give attributes of the data on which to operate.
#' - *Selection parameters*
#'   (`select`, `max_clim`, `max_geog`, `k`)
#'   define which analogs to `select` from the `pool` for each `x`.
#' - *Aggregation parameters*
#'   (`stat`, `kernel`, `theta`, `lambda`, `se`)
#'   control how selected analogs are summarized.
#' - *Computation parameters*
#'   (`n_threads`, `index_res`, `downsample`, `seed`, `progress`)
#'   specify behavior for controlling compute performance.
#'
#' ## Distance metrics
#'
#' Geographic distance can be computed for lon/lat coordinates (great-circle
#' distance) or projected coordinates (planar distance).
#'
#' Climate similarity is measured using Euclidean or Mahalanobis distance in
#' climate space. In general, when multiple climate variables are used, it is
#' recommended to use pre-whitened (scaled) climate data, to avoid major artifacts
#' from climate variables with different units. Pre-whitening can be done using
#' `scale()` for dataset-wide Euclidean distances, or `mahalanobis_transform()`
#' for dataset-wide Mahalanobis distances.
#'
#' The function also supports climate distance calculations based on
#' *local temporal* covariance structure at focal locations, via the `x_cov`
#' parameter. These local covariance values need to be pre-calculated.
#'
#' ## Computational optimization
#'
#' The analog search architecture is designed with compute performance in mind:
#'
#' - All internal computations are done in C++.
#' - Searches use a lattice-based indexing structure to efficiently
#'   search through large reference datasets. By default, the lattice
#'   resolution is tuned for optimal performance.
#' - Parallel processing is available via the `n_threads` parameter.
#' - You can `downsample` prohibitively large reference pool datasets
#'   to improve speed and memory, using a stratified sampling
#'   scheme that reduces loss of precision relative to random sampling.
#' - For large datasets, enable `progress = TRUE` to display
#'   a progress bar during computation.
#' - For raster datasets that are too large to fit in memory,
#'   `tiled_analog_search()` offers a memory-safe option.
#'
#' ## Cross-validation
#'
#' For honest prediction error when `x` and `pool` are the same dataset, use
#' [analog_cv()] or set `exclude_self = TRUE` to exclude each focal's own row
#' from its analog neighborhood.
#'
#' @param x Focal locations for which analogs will be found. Should be a
#'   matrix/data.frame with columns x, y, and climate variables, or a
#'   SpatRaster with climate variable layers.
#' @param pool The reference dataset to search for analogs. Either:
#'
#'   - Matrix/data.frame with columns x, y, and climate variables,
#'     or SpatRaster with climate variable layers, OR
#'   - An `analog_index` object created by
#'     [build_analog_index()] (for repeated queries).
#' @param x_cov Optional focal-specific covariance matrices for Mahalanobis
#'   distance calculations. Should be a matrix or data.frame with one row per
#'   focal location and one column per unique covariance component, or a
#'   SpatRaster with a layer for each component. For n climate variables,
#'   there are n*(n+1)/2 unique components, ordered as: variances first
#'   (diagonals), then covariances (upper triangle by row).
#' @param y Optional vector, factor, matrix/data.frame, or SpatRaster giving
#'   values for each reference location (must have same number of rows/cells
#'   as `pool`). Required for stats `"sum"`, `"mean"`, `"weighted_sum"`,
#'   `"weighted_mean"`, `"regression"`, and `"tabulate"`. Numeric for
#'   continuous stats; factor or coercible-to-factor (character, integer,
#'   logical) for `stat = "tabulate"`.
#' @param covariates Optional matrix/data.frame or SpatRaster giving
#'   covariate values for each reference location (must have same number of
#'   rows/cells as `pool`). Required when `stat` includes `"regression"`.
#' @param weight Optional pool site weights for use in aggregation.
#'   Numeric vector, single-column matrix/data.frame, or single-layer
#'   SpatRaster, with one value per row/cell of `pool`. For aggregation
#'   stats like `"weighted_mean"`, `"regression"`, etc., weights multiply
#'   through the weighted aggregation alongside any kernel weighting and
#'   cell-area weighting; they do not influence which analogs are selected
#'   by `knn_*` modes (selection remains distance-only). They are reported
#'   in pair mode as a `user_weight` column. Values must be non-negative;
#'   `NA` is allowed and treated as 0 (the point is excluded from
#'   aggregation). Default `NULL` means no user-supplied weights.
#'
#'   If you want to exclude a static subset of pool sites entirely, masking
#'   `pool` (and any associated `y` / `covariates`) upfront is more
#'   efficient than passing `weight = 0` for those sites, since the
#'   lattice index will not have to scan or distance-compute against them.
#'   Use `weight = 0` for cases where the mask varies per query against a
#'   shared index, or where some sites have a continuous weight and others
#'   should be excluded.
#' @param cell_area_weight Controls cell-area weighting when `pool` is a raster.
#'   One of `"auto"` (default; on for raster pools, off otherwise), `TRUE`
#'   (force on; errors if `pool` is not a SpatRaster), or `FALSE` (force
#'   off). Cell-area weights correct aggregation statistics for non-uniform
#'   cell areas (e.g. lonlat grids near the poles, or projected grids on
#'   non-equal-area projections); they are computed via
#'   `terra::cellSize()` and normalized to mean 1. When `pool` is a
#'   pre-built `analog_index`, this argument must agree with the index's
#'   stored configuration: `cell_area_weight = FALSE` errors if the index
#'   was built with cell-area weighting on (rebuild the index instead).
#' @param max_clim Maximum climate distance constraint (default: NULL = no
#'   climate constraint). Can be either:
#'
#'   - A scalar: Euclidean radius in climate space (e.g., 0.5)
#'   - A vector: Per-variable absolute differences (length must equal
#'     number of climate variables)
#'
#'   Only reference locations within this climate distance are considered.
#'   When `x_cov` is provided, scalar thresholds are interpreted in
#'   Mahalanobis distance units.
#'
#' @param max_geog Maximum geographic distance constraint (default:
#'   NULL = no geographic constraint). When specified, only reference locations
#'   within this distance are considered. Radius units should be specified in
#'   kilometers if `coord_type = "lonlat"`, or in projected coordinate units
#'   if `coord_type = "projected"`.
#'
#' @param select Character string specifying the analog selection strategy.
#'   One of:
#'
#'   - `"all"` (default): Select all analogs that satisfy the
#'     `max_clim` and `max_geog` constraints.
#'   - `"knn_clim"`: For each focal, select up to `k` analogs
#'     with smallest climate distance, subject to filters.
#'   - `"knn_geog"`: For each focal, select up to `k` analogs
#'     with smallest geographic distance, subject to filters.
#'
#' @param k Number of nearest analogs to return per focal location for kNN
#'   selection modes. Required when `select` is `"knn_geog"` or
#'   `"knn_clim"`; must be `NULL` for `select = "all"`.
#'
#' @param stat Statistic(s) used to aggregate selected analogs. Either:
#'
#'   - `NULL` or `"none"`: Return all selected analog pairs as a data.frame.
#'   - `"count"`: For each focal, count the number of selected analogs.
#'   - `"sum_weights"`: For each focal, sum the weights of selected
#'     analogs (see `kernel` and `theta`).
#'   - `"mean_weights"`: For each focal, mean of weights of selected
#'     analogs.
#'   - `"sum"`: Sum of `y` values across analogs (requires `y`).
#'   - `"mean"`: Mean of `y` values across analogs (requires `y`).
#'   - `"weighted_sum"`: Sum of (`y` × kernel weight) across analogs
#'     (requires `y` and `kernel`).
#'   - `"weighted_mean"`: Weighted mean of `y` values across analogs
#'     (requires `y` and `kernel`).
#'   - `"ess"`: Kish's effective sample size (ESS), computed as the
#'     squared sum of weights divided by the sum of squared weights
#'     (requires `kernel`).
#'   - `"regression"`: Weighted least squares (or ridge) regression
#'     of `y` on `covariates` within each analog neighborhood.
#'     Returns intercept and slope coefficients. Requires `y`,
#'     `covariates`, and `kernel`. See `lambda` for regularization.
#'   - `"tabulate"`: For each focal and each level of categorical `y`,
#'     sum the kernel weights of analogs whose `y` falls in that level.
#'     With `kernel = "uniform"` this reduces to a per-class vote count;
#'     with a distance-decay kernel it gives similarity-weighted support
#'     per class. Requires `y` (factor or coercible-to-factor) and
#'     `kernel`. Output has one column per class. `"tabulate"` is
#'     mutually exclusive with `"sum"`, `"mean"`, `"weighted_sum"`,
#'     `"weighted_mean"`, and `"regression"` (different `y` semantics);
#'     it can be combined with `"count"`, `"sum_weights"`,
#'     `"mean_weights"`, and `"ess"`.
#'   - A character vector combining multiple stats (e.g.,
#'     `c("count", "weighted_mean", "regression")`).
#'     Note: `"none"` cannot be combined with other stats.
#'
#' @param kernel Kernel decay function for weighting matches, used only when
#'   `stat` includes a weighted aggregation (`"sum_weights"`, `"mean_weights"`,
#'   `"weighted_sum"`, `"weighted_mean"`, `"ess"`, `"regression"`, or
#'   `"tabulate"`). One of:
#'
#'   - `"uniform"`: All matches weighted equally (kernel weight = 1.0).
#'   - `"inverse_clim"`: Inverse climate distance,
#'     kernel weight = 1 / (climate_distance + eps), with epsilon given by `theta`.
#'   - `"inverse_geog"`: Inverse geographic distance,
#'     kernel weight = 1 / (geographic_distance + eps), with epsilon given by `theta`.
#'   - `"gaussian_clim"`: Gaussian kernel on climate distance,
#'     kernel weight = exp(-climate_distance^2 / (2 sigma^2)), with sigma given by `theta`.
#'   - `"gaussian_geog"`: Gaussian kernel on geographic distance,
#'     kernel weight = exp(-geographic_distance^2 / (2 sigma^2)), with sigma given by `theta`.
#'   - `"gaussian_joint"`: Gaussian kernel on combined distance,
#'     kernel weight = exp(-(clim_dist^2 / (2 sigma_clim^2) + geog_dist^2 / (2 sigma_geog^2))),
#'     with sigmas given by `theta`.
#'   - `"inverse_joint"`: Inverse joint distance,
#'     kernel weight = 1 / (sqrt(clim_dist^2 + geog_dist^2) + eps), with epsilon given by `theta`.
#'
#' @param theta Optional numeric parameter controlling the shape of the
#'   weighting `kernel`, used whenever `kernel` is active (i.e. whenever
#'   `stat` includes a weighted aggregation) and `kernel` is not `"uniform"`.
#'   Interpretation depends on `kernel`:
#'
#'   - For `"inverse_clim"` or `"inverse_geog"`: epsilon value
#'     added to distances (scalar; default: 1e-12 for climate, 1e-6 for geography).
#'   - For `"gaussian_clim"` or `"gaussian_geog"`: sigma bandwidth
#'     parameter (scalar; larger values = slower decay with distance).
#'   - For `"gaussian_joint"` or `"inverse_joint"`: 2-element vector
#'     `c(theta_clim, theta_geog)` (defaults: 1 for climate, 1 for geography).
#'
#' @param lambda Ridge penalty parameter for `stat = "regression"`
#'   (default: 0, giving ordinary weighted least squares). Higher values
#'   shrink covariate coefficients toward zero, with the intercept
#'   approaching the weighted mean as `lambda -> Inf`. Ignored when
#'   `"regression"` is not in `stat`.
#'
#' @param se Standard-error framing to apply to SE-supporting stats
#'   (`"weighted_mean"` and `"regression"`). One of:
#'
#'   - `"none"` (default): no SE columns are returned.
#'   - `"ess"`: effective-sample-size framing. For `weighted_mean`,
#'     `SE = sqrt(var_w(y) / n_eff)`, where `n_eff = (Σw)² / Σw²` is Kish's
#'     effective sample size and `var_w(y) = Σwy²/Σw - ȳ_w²`. For regression,
#'     `Var(β̂) = σ²_ess · (X'WX + λI)⁻¹`, with residual variance corrected
#'     using `n_eff - p` degrees of freedom.
#'   - `"design"`: design-based framing (no assumption that weights are
#'     precisions). For `weighted_mean`,
#'     `SE = sqrt(Σ w²(y - ȳ_w)²) / Σw`.
#'
#' @param exclude_self Logical, default `FALSE`. `TRUE` is typically used
#'   for cross-validation, such as via [analog_cv()], in which case each focal
#'   excludes the pool row at the same index from its analog neighborhood.
#'   This requires `x` and `pool` to be the same R object (checked via
#'   `identical()`), and is incompatible with pre-built indices,
#'   `downsample != 1`, and `progress = TRUE`.
#'
#' @param coord_type Coordinate system type:
#'
#'   - `"auto"` (default): Automatically detect from coordinate ranges.
#'   - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle distance;
#'     assumes `max_geog` is in km).
#'   - `"projected"`: Projected XY coordinates (uses planar distance;
#'     assumes `max_geog` is in projection units).
#'
#' @param downsample Optional downsampling rate (0-1) for the reference pool,
#'   indicating the proportion of points to retain. Values < 1 reduce memory
#'   and improve speed at some cost to precision. Default is 1.0 (no downsampling).
#'   Ignored if `pool` is a pre-built index. When `downsample < 1`, `index_res`
#'   must be set explicitly (auto-tuning is not supported in this case; see
#'   the `index_res` parameter for details).
#'
#' @param seed Optional random seed for reproducible downsampling. If `NULL`
#'   (default), uses current R random state. Ignored if `pool` is a pre-built
#'   index or `downsample = 1`.
#'
#' @param index_res Tuning parameter giving the number of bins per dimension
#'   of the internally-used lattice search index. Either:
#'
#'   - A positive integer.
#'   - `"auto"` (the default): Automatically tune the index resolution
#'     by optimizing compute time on a subsample of focal points. If focal has
#'     relatively few rows, auto-tuning is skipped and a default resolution of
#'     16 is used. Auto-tuning is **not supported** when `downsample < 1`,
#'     because the speed-optimal resolution can sometimes result in higher
#'     uncertainty of stat results under downsampling. In that case set
#'     `index_res` explicitly; finer values (e.g. 32) generally give better
#'     accuracy at the possible cost of query speed.
#'
#'   Ignored if `pool` is an `analog_index` (uses index's resolution).
#'
#' @param n_threads Optional integer number of threads to use for the
#'   computation. If `NULL` (default), the global RcppParallel setting
#'   is used (see `RcppParallel::setThreadOptions`).
#'
#' @param progress Logical; if `TRUE`, display a progress bar during
#'   computation. Progress tracking works by splitting the focal dataset into
#'   chunks and processing them sequentially. Useful for large datasets. Default
#'   is `FALSE`.
#'
#' @return Return type depends on input format and query mode.
#'
#' Returns a data.frame, unless `x` is a SpatRaster and results have exactly one record per
#' input cell (aggregation mode, or pairwise with `k = 1`), in which case returns a
#' SpatRaster with one layer per output variable.
#'
#' Pairwise mode (`stat = NULL` or `"none"`) returns one row per focal-analog pair,
#' with the following variables:
#'
#' - `index`, `x`, `y`: Focal location (1-based index and coordinates) corresponding to input `x`
#' - `analog_index`, `analog_x`, `analog_y`: Analog location corresponding to input `pool`
#' - `clim_dist`: Climate distance (Euclidean or Mahalanobis)
#' - `geog_dist`: Geographic distance (km for lonlat, projection units otherwise)
#' - Value columns (if `y` provided): one per variable
#'
#' Aggregation mode (one or more `stat` values) returns one row per focal location,
#' with the following variables:
#'
#' - `index`, `x`, `y`: Focal location
#' - One column per requested statistic. For `stat` with single `y` variable:
#'   column named by stat (e.g., `sum`, `mean`). For `stat` with multiple `y`
#'   variables: columns named `{stat}_{varname}` (e.g., `sum_biomass`, `mean_richness`)
#' - For `stat = "regression"`: columns for `coef_intercept` and `coef_{covariate}`,
#'   or `coef_intercept_{varname}` and  `coef_{covariate}_{varname}` with multiple `y`
#'   variables.
#' - For `stat = "tabulate"`: one column per level of `y`, named `n_{level}` for a
#'   single unnamed `y`, or `{varname}_n_{level}` when `y` is named or has multiple
#'   columns.
#' - When `se != "none"`: matching SE columns (`se_weighted_mean`,
#'   `se_intercept`, etc.) for each SE-supporting stat.
#'
#' All results include metadata attributes (`select`, `stat`, `kernel`, etc.).
#' Use [metadata()] to retrieve them as a named list, or see
#' `?metadata` for a full reference.
#'
#' @references
#' Hamann A, Roberts DR, Barber QE, Carroll C, Nielsen SE (2015). "Velocity of
#' climate change algorithms for guiding conservation and management."
#' *Global Change Biology*, **21**(2), 997-1004.
#' \doi{10.1111/gcb.12736}
#'
#' Grenier P, Parent A-C, Huard D, Anctil F, Chaumont D (2013). "An assessment
#' of six dissimilarity metrics for climate analogs." *Journal of Applied
#' Meteorology and Climatology*, **52**(4), 733-752.
#' \doi{10.1175/JAMC-D-12-0170.1}
#'
#' @seealso [tiled_analog_search()] offers memory-safe searches on large raster
#'   datasets. Helper functions such as [analog_impact()], [analog_velocity()],
#'   and [analog_density()] offer simpler interfaces for common search types.
#'   [analog_cv()] provides cross-validation workflows.
#'
#' @export
analog_search <- function(

      # data
      x,
      pool,
      x_cov = NULL,
      y = NULL,
      covariates = NULL,
      weight = NULL,

      # candidate filtering
      max_clim = NULL,
      max_geog = NULL,
      select = "all",
      k = NULL,

      # analog aggregation
      stat = NULL,
      kernel = NULL,
      theta = NULL,
      lambda = 0,
      se = c("none", "ess", "design"),

      # cross-validation
      exclude_self = FALSE,

      # args passed to build_lattice_index
      coord_type = c("auto", "lonlat", "projected"),
      downsample = 1.0,
      seed = NULL,
      index_res = "auto",
      cell_area_weight = "auto",

      n_threads = NULL,
      progress = FALSE
) {

      se <- match.arg(se)

      # Validate exclude_self against x, pool, downsample, progress.
      # Must happen BEFORE pool is swapped for a built index below.
      .validate_exclude_self(exclude_self, x, pool, downsample, progress)

      # Validate cell_area_weight format up front so error messages mention
      # the user-facing argument and not internal helpers.
      caw_is_vec <- is.numeric(cell_area_weight) && is.null(dim(cell_area_weight))
      if (!(identical(cell_area_weight, "auto") ||
            isTRUE(cell_area_weight) ||
            isFALSE(cell_area_weight) ||
            caw_is_vec)) {
            stop('`cell_area_weight` must be "auto", TRUE, FALSE, or a numeric ',
                 "vector of length nrow(pool).",
                 call. = FALSE)
      }

      # Check if pool is already an index
      if (is_analog_index(pool)) {
            # Pool is pre-built index - use it directly.
            index <- pool

            # cell_area_weight is determined at index build time and lives on
            # the index. Reconcile the query-time argument with the stored
            # state: silently accept "auto" or matching values, error on
            # explicit disagreement.
            index_has_area <- !is.null(index$cell_area_weight)
            if (isTRUE(cell_area_weight) && !index_has_area) {
                  stop("`cell_area_weight = TRUE` was requested but the supplied ",
                       "`analog_index` was built without cell-area weighting. ",
                       "Rebuild the index with `cell_area_weight = TRUE` (or ",
                       '"auto" for a raster pool) before querying.',
                       call. = FALSE)
            }
            if (isFALSE(cell_area_weight) && index_has_area) {
                  stop("`cell_area_weight = FALSE` was requested but the supplied ",
                       "`analog_index` was built with cell-area weighting on. ",
                       "Rebuild the index with `cell_area_weight = FALSE` to ",
                       "disable it.", call. = FALSE)
            }
            if (caw_is_vec) {
                  stop("A numeric `cell_area_weight` vector cannot be applied ",
                       "to a pre-built `analog_index`; cell-area weights are ",
                       "baked in at build time. Rebuild the index with the ",
                       "desired weights instead.",
                       call. = FALSE)
            }

      } else {
            # Pool is raw data - need to build index

            # Tune resolution if needed
            if (identical(index_res, "auto")) {
                  # Auto-tuning is not safe when downsampling is in effect:
                  # tune_index_res() optimizes for query *speed*, but at
                  # downsample < 1 the speed-optimal resolution can produce
                  # systematically biased and/or high-variance aggregations.
                  # Force the user to choose explicitly.
                  if (downsample < 1) {
                        stop(
                              'Auto-tuning of `index_res` is not supported when ',
                              '`downsample < 1`. The speed-optimal resolution ',
                              'can produce biased and high-variance results ',
                              'under downsampling. Set `index_res` explicitly: ',
                              'finer values (e.g. 32) generally give better ',
                              'accuracy at the cost of query speed.',
                              call. = FALSE
                        )
                  }
                  # Use tune_index_res for auto-tuning
                  index_res_int <- tune_index_res(
                        x = x,
                        pool = pool,
                        x_cov = x_cov,
                        y = y,
                        covariates = covariates,
                        select = select,
                        stat = stat,
                        max_clim = max_clim,
                        max_geog = max_geog,
                        k = k,
                        kernel = kernel,
                        theta = theta,
                        lambda = lambda,
                        se = se,
                        coord_type = coord_type,
                        n_threads = n_threads,
                        verbose = FALSE,
                        downsample = downsample,
                        seed = seed
                  )
            } else if (is.numeric(index_res)) {
                  index_res_int <- as.integer(index_res)
            } else {
                  index_res_int <- 16L  # Default
            }

            # Build index from raw pool data
            index <- build_analog_index(
                  pool = pool,
                  coord_type = coord_type,
                  index_res = index_res_int,
                  downsample = downsample,
                  seed = seed,
                  cell_area_weight = cell_area_weight
            )
      }

      # Query the index
      return(query_analog_index(
            x = x,
            index = index,
            select = select,
            stat = stat,
            max_clim = max_clim,
            max_geog = max_geog,
            x_cov = x_cov,
            y = y,
            covariates = covariates,
            weight = weight,
            k = k,
            kernel = kernel,
            theta = theta,
            lambda = lambda,
            se = se,
            exclude_self = exclude_self,
            n_threads = n_threads,
            show_progress = progress
      ))
}
