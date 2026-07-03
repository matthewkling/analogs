#' Tune Index Resolution
#'
#' Automatically finds fast per-family lattice resolution adjustments
#' (`geog_res_adj`, `clim_res_adj`) for your data and query pattern. Uses
#' alternating coordinate descent (Gibbs-style): each active family's
#' resolution is optimized in turn (holding the other fixed) via an
#' expanding-bracket 1-D search, sweeping back and forth until the selected
#' adjustments stop changing.
#'
#' @inheritParams analog_search
#' @inheritParams build_analog_index
#'
#' @param verbose Logical; if TRUE, print search progress. Default FALSE.
#'
#' @return A list with elements `geog_res_adj` and `clim_res_adj` giving the
#'   recommended per-family resolution adjustments. A family that is inactive
#'   on entry (adjustment of 0, i.e. deactivated) is returned unchanged.
#'
#' @details
#' Each family's 1-D search starts from its current adjustment and expands a
#' multiplicative bracket (halving or doubling) in the direction of decreasing
#' compute time until an interior minimum is bracketed (time decreases then
#' increases) or a bound is reached. Adjustments are constrained to the range
#' \[1/32, 32\]. The outer loop alternates between families until a full sweep
#' leaves both selections unchanged (convergence) or a sweep cap is reached.
#'
#' Only families that are active (non-zero adjustment) on entry are tuned; a
#' deactivated family is skipped and passed through. If neither family is
#' active, or the problem is small (<= 2000 focal points), the inputs are
#' returned unchanged.
#'
#' A subsample of focal points is used for benchmarking to keep tuning fast
#' while still being representative of actual query performance.
#'
#' @examples
#' \dontrun{
#' # Tune per-family resolution for an availability query (both families active)
#' adj <- tune_index_res(
#'   x = sample_sites,
#'   pool = climate_data,
#'   select = "all",
#'   stat = "count",
#'   max_clim = 0.5,
#'   max_geog = 100
#' )
#'
#' index <- build_analog_index(climate_data,
#'                             geog_res_adj = adj$geog_res_adj,
#'                             clim_res_adj = adj$clim_res_adj)
#' }
#'
#' @export
tune_index_res <- function(x,
                           pool,
                           downsample = 1.0,
                           seed = NULL,
                           select = "all",
                           stat = NULL,
                           max_clim = NULL,
                           max_geog = NULL,
                           k = NULL,
                           kernel = NULL,
                           theta = NULL,
                           x_cov = NULL,
                           y = NULL,
                           covariates = NULL,
                           lambda = 0,
                           se = c("none", "ess", "design"),
                           coord_type = c("auto", "lonlat", "projected"),
                           geog_res_adj = 1,
                           clim_res_adj = 1,
                           n_threads = NULL,
                           verbose = FALSE) {

      se <- match.arg(se)

      if (is.null(seed)) seed <- .Random.seed[1]
      set.seed(seed)

      # Bounds on the per-family resolution adjustment during search.
      ADJ_MIN <- 1 / 32
      ADJ_MAX <- 32
      MAX_SWEEPS <- 5L

      # Which families are active (non-zero adj = tunable); deactivated families
      # (adj == 0) are passed through untouched.
      geo_active  <- is.numeric(geog_res_adj) && geog_res_adj > 0
      clim_active <- is.numeric(clim_res_adj) && clim_res_adj > 0

      result <- list(geog_res_adj = geog_res_adj, clim_res_adj = clim_res_adj)

      # Nothing to tune if neither family is active.
      if (!geo_active && !clim_active) {
            if (verbose) message("No active families to tune; returning inputs.")
            return(result)
      }

      focal <- .format_data(x)
      n_focal <- nrow(focal)

      # Only tune for non-trivial problem sizes.
      if (n_focal <= 2000L) {
            if (verbose) {
                  message("Dataset too small for tuning (n=", n_focal,
                          "); returning input adjustments.")
            }
            return(result)
      }

      # Subsample focal sites for faster benchmarking.
      n_samp <- min(1000L, max(100L, as.integer(n_focal * 0.01)))
      idx    <- sample.int(n_focal, n_samp)
      focal_samp <- focal[idx, , drop = FALSE]
      attr(focal_samp, "template") <- attr(focal, "template")

      x_cov_samp <- NULL
      if (!is.null(x_cov)) {
            x_cov_samp <- x_cov[idx, , drop = FALSE]
      }

      # Evaluate elapsed query time for a given (geo_adj, clim_adj) pair.
      eval_time <- function(geo_adj, clim_adj) {
            index <- build_analog_index(
                  pool = pool,
                  coord_type = coord_type,
                  cell_area_weight = FALSE,
                  geog_res_adj = geo_adj,
                  clim_res_adj = clim_adj,
                  downsample = downsample,
                  seed = seed
            )
            st <- system.time({
                  query_analog_index(
                        x = focal_samp,
                        index = index,
                        select = select,
                        stat = stat,
                        max_clim = max_clim,
                        max_geog = max_geog,
                        x_cov = x_cov_samp,
                        y = y,
                        covariates = covariates,
                        lambda = lambda,
                        se = se,
                        k = k,
                        kernel = kernel,
                        theta = theta,
                        n_threads = n_threads
                  )
            })
            st[["elapsed"]]
      }

      # Memoized timing cache keyed on the rounded (geo, clim) adjustment pair,
      # so re-evaluating a point across sweeps costs nothing.
      cache <- new.env(parent = emptyenv())
      key_of <- function(g, c) sprintf("%.6g|%.6g", g, c)
      timed <- function(g, c) {
            key <- key_of(g, c)
            if (!is.null(cache[[key]])) return(cache[[key]])
            val <- eval_time(g, c)
            cache[[key]] <- val
            val
      }

      # 1-D search over ONE family's adjustment, holding the other fixed.
      # `family` is "geo" or "clim"; `cur` is that family's current adjustment;
      # `other` is the fixed adjustment of the other family. Expands a
      # multiplicative bracket from `cur` toward decreasing time until an
      # interior minimum is bracketed or a bound is hit; returns the best adj.
      optimize_family <- function(family, cur, other) {
            eval_at <- function(a) {
                  if (family == "geo") timed(a, other) else timed(other, a)
            }

            cur <- min(ADJ_MAX, max(ADJ_MIN, cur))
            t_cur  <- eval_at(cur)

            # Probe both neighbors (one doubling each way) to pick a direction.
            lo <- max(ADJ_MIN, cur / 2)
            hi <- min(ADJ_MAX, cur * 2)
            t_lo <- if (lo < cur) eval_at(lo) else Inf
            t_hi <- if (hi > cur) eval_at(hi) else Inf

            if (verbose) {
                  message(sprintf("  [%s] center adj=%.3g (%.3fs); down=%.3g (%.3fs) up=%.3g (%.3fs)",
                                  family, cur, t_cur, lo, t_lo, hi, t_hi))
            }

            # If center is already a local min, done.
            if (t_cur <= t_lo && t_cur <= t_hi) {
                  return(cur)
            }

            # Expand in the better direction until time turns back up or bound hit.
            if (t_lo < t_hi) {
                  # Descend by halving.
                  best_a <- lo; best_t <- t_lo; a <- lo
                  while (a > ADJ_MIN) {
                        a_next <- max(ADJ_MIN, a / 2)
                        if (a_next == a) break
                        t_next <- eval_at(a_next)
                        if (verbose) message(sprintf("    down adj=%.3g (%.3fs)", a_next, t_next))
                        if (t_next >= best_t) break     # passed the minimum
                        best_a <- a_next; best_t <- t_next; a <- a_next
                  }
                  best_a
            } else {
                  # Descend by doubling.
                  best_a <- hi; best_t <- t_hi; a <- hi
                  while (a < ADJ_MAX) {
                        a_next <- min(ADJ_MAX, a * 2)
                        if (a_next == a) break
                        t_next <- eval_at(a_next)
                        if (verbose) message(sprintf("    up adj=%.3g (%.3fs)", a_next, t_next))
                        if (t_next >= best_t) break     # passed the minimum
                        best_a <- a_next; best_t <- t_next; a <- a_next
                  }
                  best_a
            }
      }

      # Alternating coordinate descent (Gibbs-style): optimize each active
      # family in turn (geo first, as the higher-leverage lever), holding the
      # other fixed, sweeping until neither selection changes or the cap is hit.
      geo_adj  <- if (geo_active)  min(ADJ_MAX, max(ADJ_MIN, geog_res_adj)) else geog_res_adj
      clim_adj <- if (clim_active) min(ADJ_MAX, max(ADJ_MIN, clim_res_adj)) else clim_res_adj

      converged <- FALSE
      for (sweep in seq_len(MAX_SWEEPS)) {
            if (verbose) message(sprintf("Sweep %d:", sweep))
            changed <- FALSE

            if (geo_active) {
                  new_geo <- optimize_family("geo", geo_adj, clim_adj)
                  if (!isTRUE(all.equal(new_geo, geo_adj))) changed <- TRUE
                  geo_adj <- new_geo
            }
            if (clim_active) {
                  new_clim <- optimize_family("clim", clim_adj, geo_adj)
                  if (!isTRUE(all.equal(new_clim, clim_adj))) changed <- TRUE
                  clim_adj <- new_clim
            }

            if (!changed) { converged <- TRUE; break }
      }

      if (verbose) {
            if (converged) {
                  message(sprintf("Converged: geog_res_adj=%.3g clim_res_adj=%.3g",
                                  geo_adj, clim_adj))
            } else {
                  message(sprintf("Reached sweep cap (%d) without full convergence; ",
                                  MAX_SWEEPS),
                          sprintf("using geog_res_adj=%.3g clim_res_adj=%.3g",
                                  geo_adj, clim_adj))
            }
      }

      list(geog_res_adj = as.numeric(geo_adj),
           clim_res_adj = as.numeric(clim_adj))
}
