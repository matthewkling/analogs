#' Climate impact assessment via analog impact model
#'
#' Predicts ecological outcome variables based on climatic similarity and geographic
#' proximity. When focal and reference data are from the same time period, performs a
#' climate-informed spatial interpolation. When they are from different eras, implements
#' an analog impact model (AIM) projecting the potential ecological state under climate
#' change. For each focal location's future climate, identifies locations with similar
#' baseline climates within a specified geographic range, then aggregates their ecological
#' characteristics weighted by climate similarity. Aggregation can be a weighted mean
#' for continuous outcomes or a weighted class count for categorical outcomes.
#'
#' @param stat Statistic(s) to compute across analogs (default: c("count",
#'   "sum_weights", "weighted_mean")). See [analog_search()] for the full list
#'   of options. Common choices:
#'   \itemize{
#'     \item `"weighted_mean"`: Expected ecological state under continuous y
#'     \item `"tabulate"`: Per-class similarity-weighted support under
#'       categorical y. Output has one column per class, named
#'       `n_<level>` (single-y) or `<varname>_n_<level>` (multi-y).
#'       `"tabulate"` is mutually exclusive with `weighted_mean`/`sum`/
#'       `mean`/`weighted_sum`/`regression`, but can be combined with
#'       `count`, `sum_weights`, `mean_weights`, and `ess`.
#'     \item `"count"`: Analog availability (number of analogs)
#'     \item `"sum_weights"`: Analog density (sum of climate-similarity
#'       weights across analogs)
#'   }
#'   The default `c("count", "sum_weights", "weighted_mean")` is appropriate
#'   for continuous `y`; for categorical `y`, swap `weighted_mean` for
#'   `tabulate`.
#' @param env,geog Per-family [kernel()] objects giving the distance treatment
#'   for environment and geography. For climate impact models, `env` carries the
#'   climate analog threshold (`max`) and the weighting kernel applied to
#'   climate distance (`weight`, `theta`); it defaults to
#'   `kernel("gaussian", theta = 0.25, max = 1.0)`. `geog` carries the hard
#'   distance constraint (`geog = kernel(max = ...)`) and defaults to
#'   `NULL` (no geographic constraint). A geographic kernel may also be supplied
#'   via `geog` if desired. See [kernel()].
#' @inheritParams analog_search
#' @param ... Additional arguments passed to [analog_search()].
#'
#' @details
#'
#' ## The Analog Impact Model (AIM) Framework
#'
#' When `x` represents future climate and `pool` represent baseline climate,
#' this function implements the "reverse analog" approach from the climate
#' change ecology literature. It addresses the question, "For a location's
#' future climate, what ecological conditions exist in current locations with
#' similar climates that are within dispersal range?"
#'
#' The methodology:
#' \enumerate{
#'   \item For each focal location's future climate conditions
#'   \item Find all current locations with similar climates (within `max_env`)
#'   \item Constrain to dispersal-reachable distance (within `max_geog`)
#'   \item Weight each analog by climate similarity (via `kernel` function)
#'   \item Aggregate ecosystem characteristics across these weighted analogs
#' }
#'
#' Unlike traditional AIM implementations that select k nearest climate neighbors,
#' this function uses all analogs within thresholds combined with climate-distance-based
#' kernel weighting. This approach eliminates arbitrary choice of k, provides smoother,
#' more continuous results, and lets the kernel function (via `theta`) naturally
#' control influence. (Note that the traditional version can be implemented via
#' `analog_search(select = "knn_env", stat = "mean", ...))`.)
#'
#' ## Choosing Parameters
#'
#' \itemize{
#'   \item `geog`'s `max`: Set based on species dispersal ability (e.g., 5-500 km)
#'   \item `env`'s `max`: Defines what counts as an "analog"
#'   \item `env`'s `theta`: Controls kernel decay. The weight should decay to a
#'     small value at the climate/geographic thresholds. If `theta` is too large
#'     relative to the thresholds, the hard cutoffs do most of the filtering and
#'     weighting becomes nearly uniform. For Gaussian kernels with three or fewer
#'     climate variables, a reasonable rule of thumb is to set `theta` to
#'     `max / 3`.
#' }
#'
#' ## Interpreting Results
#'
#' \itemize{
#'   \item `weighted_mean`: Expected ecosystem state if colonized by species
#'     from analog locations.
#'   \item `tabulate`: For each class in `y`, the sum of analog weights that
#'     fall in that class. Each class's column gives the total
#'     climatic-similarity-weighted support among analogs (or, with
#'     `analog_search(stat = "tabulate")` with a uniform kernel, a plain
#'     vote count). The "primary" projection at a focal location is
#'     `which.max()` of these columns; "agreement" is the largest column
#'     value divided by the row sum (or, equivalently, divided by
#'     `sum_weights` when both are requested). Bray-Curtis dissimilarity
#'     between two such weighted-vote matrices (e.g., contemporary vs.
#'     future analogs) gives a measure of compositional vulnerability to
#'     ecological transformation (cf. Hoecker et al. 2026).
#'   \item `count`: How many analogs exist within max_env and max_geog? Low counts
#'     indicate limited analog availability, while zero counts indicate climates
#'     that are novel within the geographic search radius.
#'   \item `sum_weights`: Total analog density. Low values indicate sparse
#'     or distant climate matches. This metric captures both the number and quality
#'     of analogs. Interpretation details vary based on the `kernel` parameter.
#'     Pass `normalize = TRUE` to express `sum_weights` (and any `tabulate`
#'     columns) as a fraction of theoretical maximum analog availability,
#'     on roughly `[0, 1]`. See [analog_search()] for preconditions.
#' }
#'
#' @return A data.frame, or a SpatRaster when `x` is one. Contains
#'   `index`, `x`, `y` plus one or more columns determined by `stat`.
#'   See [analog_search()] for column-naming conventions across stats
#'   and [metadata()] for attached metadata attributes.
#'
#' @seealso [analog_search()] for the underlying flexible analog search function;
#'   [analog_cv()] for cross-validation of AIM fits.
#'
#' @examples
#' \dontrun{
#' # Basic climate impact assessment with continuous response
#' impact <- analog_impact(
#'   x = future_climate,
#'   pool = current_climate,
#'   y = current$habitat,
#'   geog = kernel(max = 100),  # 100 km dispersal range
#'   env = kernel("gaussian", theta = 0.25, max = 0.5)
#' )
#'
#' # With uncertainty quantification on weighted_mean
#' impact_se <- analog_impact(
#'   x = future_climate,
#'   pool = current_climate,
#'   y = current$habitat,
#'   geog = kernel(max = 100),
#'   env = kernel("gaussian", theta = 0.25, max = 0.5),
#'   se = "ess"
#' )
#'
#' # Categorical response (e.g. vegetation type):
#' # compute weighted-vote counts per class within climatic neighborhoods
#' veg_votes <- analog_impact(
#'   x    = future_climate,
#'   pool = current_climate,
#'   y    = factor(current$vegetation_type),
#'   stat = c("count", "sum_weights", "tabulate"),
#'   env = kernel("gaussian", theta = 0.25, max = 1),
#'   geog = kernel(max = 500),
#'   # (env set above)
#' )
#' # Output has one `n_<level>` column per vegetation type.
#' # Primary class per focal:  apply(veg_votes[, grep("^n_", names(veg_votes))],
#' #                                   1, function(r) names(r)[which.max(r)])
#' # Agreement (max share):    apply(votes_mat, 1, max) / rowSums(votes_mat)
#' }
#'
#' @export
analog_impact <- function(
            x,
            pool,
            y,
            weight = NULL,
            covariates = NULL,
            geog = NULL,
            env = kernel("gaussian", theta = 0.25, max = 1.0),
            stat = c("weighted_mean", "count", "sum_weights"),
            lambda = 0,
            se = c("none", "ess", "design"),
            normalize = "auto",
            x_cov = NULL,
            coord_type = "auto",
            env_res_adj = "auto",
            geog_res_adj = "auto",
            cell_area_weight = "auto",
            n_threads = NULL,
            progress = FALSE,
            ...
) {
      se <- match.arg(se)

      # Validate stat values. We don't use match.arg here because stat can be
      # a vector of multiple stats; instead we restrict to the set this
      # wrapper exposes. (analog_search() owns the master valid_stats list.)
      allowed_impact_stats <- c("count", "sum_weights", "mean_weights",
                                "weighted_mean", "ess", "tabulate", "regression")
      bad <- setdiff(stat, allowed_impact_stats)
      if (length(bad) > 0L) {
            stop("`stat` value(s) not supported by analog_impact(): ",
                 paste(sprintf("'%s'", bad), collapse = ", "),
                 ". Allowed: ",
                 paste(sprintf("'%s'", allowed_impact_stats), collapse = ", "),
                 ". For pure pair output, use analog_search(stat = NULL).",
                 call. = FALSE)
      }

      analog_search(
            x           = x,
            pool        = pool,
            select      = "all",
            stat        = stat,
            y           = y,
            weight      = weight,
            covariates  = covariates,
            env         = env,
            geog        = geog,
            k           = NULL,
            lambda      = lambda,
            se          = se,
            normalize   = normalize,
            x_cov       = x_cov,
            coord_type  = coord_type,
            env_res_adj= env_res_adj,
            geog_res_adj = geog_res_adj,
            cell_area_weight = cell_area_weight,
            n_threads   = n_threads,
            progress    = progress,
            ...
      )
}
