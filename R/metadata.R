#' Get metadata from analog search results
#'
#' Returns the metadata attributes attached to results from `analog_search()`,
#' or any `analog_*()` wrapper function that calls it, as a named list. These
#' same values can be accessed individually via `attr(x, "<name>")`;
#' `metadata()` is a convenience that bundles them into a consistent
#' list structure independent of the underlying attribute layout.
#' Use it to inspect or programmatically reconstruct the parameterization
#' that produced a result. Each attribute records one user-facing parameter
#' or input-data property.
#'
#' @param x A result from any `analog_*()` function (data.frame or
#'   SpatRaster).
#'
#' @return A named list. Elements depend on which attributes are present on
#'   `x`; missing attributes are simply absent from the returned list (no
#'   `NA` placeholders). Possible elements include:
#'
#'   ## Selection and aggregation parameters
#'
#'   \describe{
#'     \item{`select`}{Selection strategy: `"all"`, `"knn_clim"`, or
#'       `"knn_geog"`.}
#'     \item{`k`}{Number of neighbors for kNN selection modes; `NULL` for
#'       `select = "all"`.}
#'     \item{`stat`}{Aggregation statistic(s) requested. Character vector;
#'       `"none"` for pair-mode results.}
#'     \item{`kernel`}{Kernel function used for analog weighting (e.g.
#'       `"gaussian_clim"`), or `NULL` if no weighted aggregation was
#'       requested.}
#'     \item{`theta`}{Kernel bandwidth or epsilon parameter; scalar for
#'       most kernels, length-2 numeric for `*_joint` kernels.}
#'     \item{`lambda`}{Ridge penalty for `stat = "regression"`; `0` for
#'       ordinary weighted least squares.}
#'     \item{`se`}{Standard-error framing applied to SE-supporting stats:
#'       `"none"`, `"ess"`, or `"design"`.}
#'     \item{`max_clim`}{Climate-distance threshold used for analog
#'       selection. Scalar for Euclidean / Mahalanobis radius, or
#'       length-`n_clim` for per-variable thresholds.}
#'     \item{`max_geog`}{Geographic-distance threshold used for analog
#'       selection. Units are km when `coord_type = "lonlat"`, projection
#'       units otherwise.}
#'     \item{`exclude_self`}{Logical: was each focal's own pool row
#'       excluded from its analog neighborhood? `TRUE` for results from
#'       [analog_cv()] (LOO) or from `analog_search(exclude_self = TRUE)`.}
#'   }
#'
#'   ## Input-data properties
#'
#'   \describe{
#'     \item{`n_x`}{Number of focal locations in the input `x`.}
#'     \item{`n_pool`}{Number of reference locations in the input `pool`.}
#'     \item{`n_clim`}{Number of climate variables.}
#'     \item{`coord_type`}{Coordinate system: `"lonlat"` (great-circle
#'       distances in km) or `"projected"` (planar distances in projection
#'       units).}
#'     \item{`x_cov_provided`}{Logical: was a per-focal covariance matrix
#'       (`x_cov`) supplied? If `TRUE`, climate distances were computed in
#'       Mahalanobis units using location-specific covariance; otherwise
#'       Euclidean.}
#'     \item{`downsample_actual`}{Actual downsampling rate applied to the
#'       reference pool (1.0 if no downsampling). May exceed the requested
#'       rate when maintaining at least one point per occupied bin
#'       requires it.}
#'   }
#'
#'   ## Index diagnostics (lattice-related)
#'
#'   These attributes describe the internal lattice index used to
#'   accelerate the search; they're useful for performance debugging but
#'   not needed to reproduce the result. See [build_analog_index()] for
#'   detail.
#'
#'   \describe{
#'     \item{`binning_method`}{Internal indexing method; one of
#'       `"lattice"` or `"lattice_ecef"`.}
#'     \item{`total_bins`, `n_bins_nonempty`, `min_bin_occupancy`,
#'       `max_bin_occupancy`, `avg_bin_occupancy`,
#'       `avg_nonempty_bin_occupancy`}{Lattice occupancy summaries.}
#'   }
#'
#'   ## Cross-validation metadata
#'
#'   Present only when `x` is a result from [analog_cv()].
#'
#'   \describe{
#'     \item{`cv_method`}{`"loo"` or `"kfold"`.}
#'     \item{`cv_fun`}{Name of the analog function being cross-validated
#'       (e.g. `"analog_impact"`).}
#'     \item{`cv_n_folds`}{Number of folds (k-fold only).}
#'     \item{`cv_pred_target`}{Mechanism used to extract predictions for
#'       residuals: `"weighted_mean"`, `"regression"`, or `"none"`.}
#'   }
#'
#' @details
#' `metadata()` is the documented entry point for reading these
#' attributes; the attribute names themselves are also part of the public
#' contract, so `attr(x, "select")` etc. are equally valid for individual
#' values. The function's main advantages over raw `attributes()` are: it
#' filters out R-internal attributes (`names`, `class`, `dim`, ...), and
#' the `?metadata` help page gives a single canonical reference for
#' every attribute the package attaches.
#'
#' @examples
#' \dontrun{
#' result <- analog_velocity(
#'   x = future_climate,
#'   pool = current_climate,
#'   max_clim = 0.5
#' )
#'
#' # Get the full parameterization as a list
#' meta <- metadata(result)
#' meta$select
#' meta$max_clim
#' meta$n_pool
#'
#' # Or read individual attributes directly
#' attr(result, "select")
#' }
#'
#' @export
metadata <- function(x) {

      # Set of attribute names this package documents as part of the public
      # contract. Attributes not in this set are filtered out so the result
      # is stable even if R-internal or implementation-detail attrs are
      # attached.
      documented <- c(
            # Selection and aggregation
            "select", "k", "stat", "kernel", "theta",
            "lambda", "se",
            "max_clim", "max_geog",
            "exclude_self",
            # Input-data properties
            "n_x", "n_pool", "n_clim", "coord_type", "x_cov_provided",
            "downsample_actual",
            # Index diagnostics
            "binning_method", "total_bins", "n_bins_nonempty",
            "min_bin_occupancy", "max_bin_occupancy",
            "avg_bin_occupancy", "avg_nonempty_bin_occupancy",
            # Cross-validation
            "cv_method", "cv_fun", "cv_n_folds", "cv_pred_target"
      )

      attrs <- attributes(x)
      keep <- intersect(documented, names(attrs))
      attrs[keep]
}
