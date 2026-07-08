#' Generalized analog search
#'
#' Identifies locations in a reference dataset that are climatically similar
#' and/or geographically proximal to focal locations. Analog searches use a
#' two-stage approach: first selecting analogs based on specified criteria,
#' then optionally aggregating the results.
#'
#'
#' @param x Focal locations for which analogs will be found. Should be a
#'   matrix/data.frame with columns x, y, and climate variables, or a
#'   SpatRaster with climate variable layers.
#'
#' @param pool The reference dataset to search for analogs. Either:
#'
#'   - Matrix/data.frame with columns x, y, and climate variables,
#'     or SpatRaster with climate variable layers, OR
#'   - An `analog_index` object created by
#'     [build_analog_index()] (for repeated queries).
#'
#' @param x_cov Optional focal-specific covariance matrices for Mahalanobis
#'   distance calculations. Should be a matrix or data.frame with one row per
#'   focal location and one column per unique covariance component, or a
#'   SpatRaster with a layer for each component. For n climate variables,
#'   there are n*(n+1)/2 unique components, ordered as: variances first
#'   (diagonals), then covariances (upper triangle by row).
#'
#' @param y Optional vector, factor, matrix/data.frame, or SpatRaster giving
#'   values for each reference location (must have same number of rows/cells
#'   as `pool`). Required for stats `"sum"`, `"mean"`, `"weighted_sum"`,
#'   `"weighted_mean"`, `"regression"`, and `"tabulate"`. Numeric for
#'   continuous stats; factor or coercible-to-factor (character, integer,
#'   logical) for `stat = "tabulate"`.
#'
#' @param covariates Optional matrix/data.frame or SpatRaster giving
#'   covariate values for each reference location (must have same number of
#'   rows/cells as `pool`). Required when `stat` includes `"regression"`.
#'
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
#'
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
#'
#' @param clim,geog Per-family distance treatment, each a [kernel()] object (or
#'   `NULL`). A kernel bundles the hard distance threshold, the weighting kernel
#'   shape, and the kernel's scale for one family: climate (`clim`) or geography
#'   (`geog`). `kernel(weight, theta, max)` where:
#'
#'   - `max`: hard distance threshold — candidates beyond it (in that family's
#'     distance) are excluded. For `clim`, `max` may be a single Euclidean
#'     radius or a per-variable vector of absolute-difference thresholds (length
#'     equal to the number of climate variables); scalar climate thresholds are
#'     in Mahalanobis units when `x_cov` is supplied. For `geog`, `max` is a
#'     single radius (kilometers when `coord_type = "lonlat"`, projected units
#'     otherwise).
#'   - `weight`: kernel shape for weighted aggregations — `"uniform"` (no
#'     distance weighting), `"gaussian"` (`exp(-d^2 / (2 theta^2))`), or
#'     `"inverse"` (`1 / (1 + d / theta)`). The overall kernel weight is the
#'     product of the two families' weights, so shapes may be mixed (e.g. an
#'     inverse climate kernel with a Gaussian geographic kernel).
#'   - `theta`: the kernel's scale (Gaussian bandwidth, or inverse half-weight
#'     distance). See [kernel_params()] for calibrated values.
#'
#'   A `NULL` kernel (the default for both) applies no threshold and no
#'   weighting for that family. See [kernel()] for details.
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
#'     analogs (see `kernel` and `theta`). When `normalize = TRUE`, the
#'     reported value is the normalized density `D / D_max`, on roughly
#'     `[0, 1]`; otherwise it is the raw kernel-weight sum.
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
#'   - `"tabulate"`: if `y` is categorical, separately sum the kernel weights
#'     of analogs matching each level of `y`.
#'     With a uniform kernel (no `clim`/`geog` distance weighting) this
#'     reduces to a per-class vote count; with a distance-decay kernel it
#'     gives similarity-weighted support per class. Requires `y` (factor or
#'     coercible-to-factor).
#'     Output has one column per class. `"tabulate"` is mutually exclusive with
#'     `"sum"`, `"mean"`, `"weighted_sum"`, `"weighted_mean"`, and `"regression"`
#'     (different `y` semantics); it can be combined with `"count"`,
#'     `"sum_weights"`, `"mean_weights"`, and `"ess"`.
#'   - A character vector combining multiple stats (e.g.,
#'     `c("count", "weighted_mean", "regression")`).
#'     Note: `"none"` cannot be combined with other stats.
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
#' @param normalize One of `TRUE`, `FALSE`, or `"auto"` (default). Only used
#'   if `stat` includes `"sum_weights"` or `"tabulate"` and `pool` is a raster.
#'   When active, results for these stats are divided by a global scalar so that
#'   they represent a fraction of a theoretically "perfect" scenario where
#'   the full search area within `max_geog` is occupied wall-to-wall by cells
#'   whose climate exactly matches `x`. See details under [analog_search()] for
#'   more info.
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
#'   Ignored if `pool` is a pre-built index. When `downsample < 1`, resolution
#'   must be set explicitly via `geog_res_adj` / `clim_res_adj` (auto-tuning is
#'   not supported in this case; see those parameters for details).
#'
#' @param seed Optional random seed for reproducible downsampling. If `NULL`
#'   (default), uses current R random state. Ignored if `pool` is a pre-built
#'   index or `downsample = 1`.
#'
#' @param clim_res_adj,geog_res_adj Control the lattice search-index resolution
#'   of the climate and geographic families, each a multiplier on a
#'   data-dependent default (targeting ~50 pool points per occupied bin, split
#'   between families by effective dimensionality, so it scales with pool size).
#'   Each is either:
#'
#'   - A non-negative number: `1` uses the default for that family, larger
#'     values are finer, smaller are coarser, and `0` deactivates it.
#'   - `"auto"` (the default for both): tune a single overall resolution scale
#'     by optimizing compute time on a subsample of focal points. If focal has
#'     relatively few rows, tuning is skipped. Not supported when
#'     `downsample < 1` (set explicit numeric values instead).
#'
#'   A family that the query does not constrain (no corresponding `max_*` and
#'   not the knn sort key) is **automatically deactivated**, overriding any
#'   explicit value (with a message), since binning an unconstrained family
#'   only costs time. Ignored if `pool` is an `analog_index` (uses the index's
#'   resolution).
#'
#' @param mean_cell_area Optional scalar mean cell area to attach to the
#'   index when one is built from raw `pool` data. Mainly intended for
#'   internal use by `tiled_analog_search()` to propagate a globally
#'   consistent value across per-tile index builds; most users should
#'   leave this `NULL` (auto-computed from the raster pool). See
#'   [build_analog_index()] for details.
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
#' Use [metadata()] to retrieve them as a named list, or see  `?metadata` for a full
#' reference.
#'
#' @details
#' ## Parameter categories
#'
#' - *Data parameters*
#'   (`x`, `pool`, `x_cov`, `y`, `covariates`, `weight`, `coord_type`)
#'   give attributes of the data on which to operate.
#' - *Selection parameters*
#'   (`select`, `clim`, `geog`, `k`)
#'   define which analogs to `select` from the `pool` for each `x`.
#' - *Aggregation parameters*
#'   (`stat`, `clim`, `geog`, `lambda`, `se`, `normalize`)
#'   control how selected analogs are summarized. The `clim` / `geog` kernels
#'   carry both the selection thresholds and the weighting kernels.
#' - *Computation parameters*
#'   (`n_threads`, `geog_res_adj`, `clim_res_adj`, `downsample`, `seed`, `progress`)
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
#' ## Normalization
#'
#' Normalization divides `D` (the density result from `sum_weights` or `tabulate`)
#' by the  global scalar `D_max`. `D_max` (which is also attached to the result as
#' an attribute) is the highest `D` you could theoretically expect given `max_geog`,
#' `kernel`, and `theta`, i.e. the density value you'd get if the entire the geographic
#' search radius were filled with grid cells whose climate exactly matches `x`.
#' It is calculated as the analytic integral
#' `(1 / mean_cell_area) * integral_0^max_geog K(0, r) * 2*pi*r dr`,
#' which is the kernel-weighted count an idealized focal would accumulate
#' from a continuous uniform pool of perfect climate matches out to
#' `max_geog`. The resulting columns are unitless availability fractions
#' on roughly `[0, 1]`.
#'
#' Because `D_max` is a continuous-kernel idealization while `D` is a discrete sum over a finite grid, normalized
#' values can occasionally exceed 1 by small amounts (typically a few percent). This
#' is a grid discretization artifact, not an error, and at certain `(cell_size, max_geog)`
#' ratios this is more pronounced. Using a higher-resolution pool grid or choosing
#' a `max_geog` that isn't an integer multiple of the cell size both reduce the effect.
#'
#' `normalize = "auto"` activates normalization if every precondition
#' is met: raster-derived index with cell-area weighting, a kernel set
#' (any of the supported types), and a finite positive `max_geog`.
#' Explicit `TRUE` errors on any missing precondition. `normalize` is silently
#' ignored when no normalizable stat is requested. For non-raster pools, `"auto"` falls
#' back to raw kernel-weighted sums. Pass `normalize = TRUE` to require normalization or
#' `normalize = FALSE` to always return raw sums.
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
#' @seealso [kernel_params()] recommends parameter values calibrated to target
#'   kernel coverage. [tiled_analog_search()] offers memory-safe searches on large raster
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

      # candidate filtering + distance weighting (per-family kernels)
      clim = NULL,
      geog = NULL,
      select = "all",
      k = NULL,

      # analog aggregation
      stat = NULL,
      lambda = 0,
      se = c("none", "ess", "design"),
      normalize = "auto",

      # cross-validation
      exclude_self = FALSE,

      # args passed to build_lattice_index
      coord_type = c("auto", "lonlat", "projected"),
      downsample = 1.0,
      seed = NULL,
      geog_res_adj = "auto",
      clim_res_adj = "auto",
      cell_area_weight = "auto",
      mean_cell_area = NULL,

      n_threads = NULL,
      progress = FALSE
) {

      se <- match.arg(se)

      # Unpack the per-family kernels (clim, geog) into the individual
      # components used throughout: hard thresholds (max_clim/max_geog), and
      # per-family kernel shapes + scales. A NULL kernel means no threshold and
      # no weighting for that family. This is the single point where the
      # kernel() sugar is dissolved; everything downstream consumes the plain
      # per-family values. Validation of the components happened in kernel().
      .unpack_kernel <- function(w, arg) {
            if (is.null(w)) return(list(max = NULL, weight = "uniform", theta = NULL))
            if (!inherits(w, "analog_kernel")) {
                  stop("`", arg, "` must be a kernel() object or NULL.", call. = FALSE)
            }
            list(max    = w$max,
                 weight = w$weight %||% "uniform",
                 theta  = w$theta)
      }
      .clim_w <- .unpack_kernel(clim, "clim")
      .geog_w <- .unpack_kernel(geog, "geog")

      max_clim    <- .clim_w$max
      max_geog    <- .geog_w$max
      kernel_clim <- .clim_w$weight   # "uniform" | "gaussian" | "inverse"
      kernel_geog <- .geog_w$weight
      theta_clim  <- .clim_w$theta
      theta_geog  <- .geog_w$theta

      # Validate normalize argument format up front. Accepts TRUE, FALSE,
      # or the string "auto". Compatibility with kernel/max_geog/index/etc.
      # is checked inside query_analog_index() once those values are
      # resolved (and "auto" is resolved to a logical based on whether the
      # index supplies the structural prerequisites for normalization).
      if (!(identical(normalize, "auto") ||
            (is.logical(normalize) && length(normalize) == 1L &&
             !is.na(normalize)))) {
            stop('`normalize` must be TRUE, FALSE, or "auto".', call. = FALSE)
      }

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

            # mean_cell_area also lives on the index. Silently accept NULL
            # (the common case); error on explicit disagreement so users
            # don't think a passed value is being honored.
            if (!is.null(mean_cell_area)) {
                  stop("`mean_cell_area` cannot be applied to a pre-built ",
                       "`analog_index`; it was set at build time. Rebuild ",
                       "the index with the desired `mean_cell_area`, or ",
                       "leave the argument NULL.", call. = FALSE)
            }

      } else {
            # Pool is raw data - need to build index.

            # Auto-deactivation: a family that this query neither filters (finite
            # max_*) nor sorts on (knn mode) prunes nothing, so we deactivate its
            # lattice (res_adj = 0), directing all resolution to the constraining
            # family. This is derived from the query and OVERRIDES any user value
            # (binning an unconstrained family only wastes time); we message when
            # an explicit non-zero user value is overridden.
            geo_constrained  <- !is.null(max_geog) || identical(select, "knn_geog")
            clim_constrained <- !is.null(max_clim) || identical(select, "knn_clim")

            # Resolve requested adjustments. "auto" (the default for both)
            # requests speed tuning of a single overall resolution scale applied
            # to both families; a numeric value is used directly.
            tune_requested <- identical(geog_res_adj, "auto") ||
                  identical(clim_res_adj, "auto")
            geo_adj_req  <- if (identical(geog_res_adj,  "auto")) 1 else geog_res_adj
            clim_adj_req <- if (identical(clim_res_adj, "auto")) 1 else clim_res_adj

            # Apply deactivation override (with a message if it changes an
            # explicit user value).
            if (!geo_constrained && !identical(geog_res_adj, "auto") &&
                is.numeric(geog_res_adj) && geog_res_adj != 0) {
                  message("Query does not constrain geography; deactivating the ",
                          "geographic lattice (overriding geog_res_adj = ",
                          geog_res_adj, ").")
            }
            if (!clim_constrained && !identical(clim_res_adj, "auto") &&
                is.numeric(clim_res_adj) && clim_res_adj != 0) {
                  message("Query does not constrain climate; deactivating the ",
                          "climate lattice (overriding clim_res_adj = ",
                          clim_res_adj, ").")
            }
            geo_adj_use  <- if (geo_constrained)  geo_adj_req  else 0
            clim_adj_use <- if (clim_constrained) clim_adj_req else 0

            # Resolve tuning (minimal: a single overall scale on both families).
            if (tune_requested) {
                  # Auto-tuning is not safe when downsampling is in effect:
                  # tune_index_res() optimizes for query *speed*, but at
                  # downsample < 1 the speed-optimal resolution can produce
                  # systematically biased and/or high-variance aggregations.
                  if (downsample < 1) {
                        stop(
                              'Auto-tuning of resolution is not supported when ',
                              '`downsample < 1`. The speed-optimal resolution ',
                              'can produce biased and high-variance results ',
                              'under downsampling. Set `geog_res_adj` / ',
                              '`clim_res_adj` explicitly: finer values generally ',
                              'give better accuracy at the cost of query speed.',
                              call. = FALSE
                        )
                  }
                  # Tune per-family resolution (coordinate descent over the
                  # active families; deactivated families are passed through).
                  tuned <- tune_index_res(
                        x = x,
                        pool = pool,
                        x_cov = x_cov,
                        y = y,
                        covariates = covariates,
                        select = select,
                        stat = stat,
                        clim = clim,
                        geog = geog,
                        k = k,
                        lambda = lambda,
                        se = se,
                        coord_type = coord_type,
                        geog_res_adj = geo_adj_use,
                        clim_res_adj = clim_adj_use,
                        n_threads = n_threads,
                        verbose = FALSE,
                        downsample = downsample,
                        seed = seed
                  )
                  geo_adj_use  <- tuned$geog_res_adj
                  clim_adj_use <- tuned$clim_res_adj
            }

            # Build index from raw pool data
            index <- build_analog_index(
                  pool = pool,
                  coord_type = coord_type,
                  geog_res_adj = geo_adj_use,
                  clim_res_adj = clim_adj_use,
                  downsample = downsample,
                  seed = seed,
                  cell_area_weight = cell_area_weight,
                  mean_cell_area = mean_cell_area
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
            kernel_clim = kernel_clim,
            kernel_geog = kernel_geog,
            theta_clim = theta_clim,
            theta_geog = theta_geog,
            lambda = lambda,
            se = se,
            normalize = normalize,
            exclude_self = exclude_self,
            n_threads = n_threads,
            show_progress = progress
      ))
}
