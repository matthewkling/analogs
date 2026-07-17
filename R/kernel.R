#' Define the distance treatment for one analog dimension
#'
#' Constructs a specification of how a single distance family — environment or
#' geography — is filtered and weighted in an analog search. A `kernel`
#' bundles the per-family choices: the hard distance thresholds (`max` and
#' `min`), the weighting kernel shape (`weight`), and that kernel's scale
#' parameter (`theta`). Pass one to the `env` argument of [analog_search()]
#' and/or one to the `geog` argument.
#'
#' The overall kernel weight for a candidate is the product of the two
#' families' weights, so the families are specified independently and may use
#' different shapes (e.g. an inverse-distance environmental kernel together with a
#' Gaussian geographic kernel). A family with `weight = "uniform"` (or `NULL`)
#' contributes a constant weight of 1, i.e. it filters (if `max`/`min` are set)
#' but does not down-weight by distance.
#'
#' All four components are optional. Which combinations are valid depends on
#' the operation and is checked by [analog_search()] downstream (for example,
#' climate velocity requires an environmental `max`; a weighted statistic requires a
#' non-uniform `weight` on at least one family). A bare `NULL` passed as the
#' `env` or `geog` argument is equivalent to `kernel()` with all components
#' unset: no threshold and no weighting for that family.
#'
#' @param weight Kernel shape for this family. One of `"uniform"` (constant
#'   weight 1; the default when `NULL`), `"gaussian"`
#'   (`exp(-d^2 / (2 theta^2))`), or `"inverse"` (`1 / (1 + d / theta)`, a
#'   heavy-tailed inverse-distance kernel). The family (environmental vs geographic)
#'   is determined by whether the kernel is passed as `env` or `geog`, so the
#'   shape name is unqualified.
#' @param theta Scale parameter for the `weight` kernel. For `"gaussian"` it is
#'   the bandwidth (sigma); for `"inverse"` it is the half-weight distance (the
#'   weight is 1/2 at `d = theta`). Ignored for `"uniform"`. Defaults to `NULL`,
#'   which lets downstream code apply a default of 1. See [kernel_params()] for
#'   help choosing a value calibrated to a target coverage fraction.
#' @param max Hard upper distance threshold for this family: candidates beyond
#'   `max` (in this family's distance) are excluded. `NULL` (default) means no
#'   upper threshold. Usually a single radius. For the environmental family,
#'   `max` may also be a vector of per-variable absolute-difference thresholds
#'   (length equal to the number of environmental variables); the geographic
#'   family uses a single radius. Supplied to the search as `max_env` /
#'   `max_geog`.
#' @param min Hard lower distance threshold for this family: candidates closer
#'   than `min` (in this family's distance) are excluded, so the retained
#'   candidates form an annulus `min <= d <= max`. A single positive scalar, or
#'   `NULL` (default) for no lower threshold. Currently supported only for the
#'   **geographic** family; setting `min` on an environmental kernel is an
#'   error. The primary use case is buffered spatial cross-validation (e.g. via
#'   [analog_cv()]), where excluding geographically near-duplicate candidates
#'   around each focal gives a less optimistic estimate of predictive skill
#'   than plain leave-one-out. Supplied to the search as `min_geog`.
#'
#' @return An object of class `"analog_kernel"`: a list with elements
#'   `weight`, `theta`, `max`, and `min` (each possibly `NULL`).
#'
#' @seealso [analog_search()], [kernel_params()]
#'
#' @examples
#' # Environment: keep analogs within 2 environmental-distance units, Gaussian-weighted
#' kernel(weight = "gaussian", theta = 0.5, max = 2)
#'
#' # Geography: hard 100 km cutoff, no distance weighting (uniform)
#' kernel(max = 100)
#'
#' # Geography: annulus keeping analogs between 5 km and 100 km of each focal
#' # (e.g. a 5 km buffer for spatial cross-validation)
#' kernel(max = 100, min = 5)
#'
#' # Inverse-distance environmental weighting, no hard cutoff
#' kernel("inverse", theta = 1)
#'
#' @export
kernel <- function(weight = NULL, theta = NULL, max = NULL, min = NULL) {

      # Validate weight (shape). NULL is allowed and means "uniform".
      valid_weights <- c("uniform", "gaussian", "inverse")
      if (!is.null(weight)) {
            if (!is.character(weight) || length(weight) != 1L ||
                !weight %in% valid_weights) {
                  stop("`weight` must be one of ",
                       paste(sprintf('"%s"', valid_weights), collapse = ", "),
                       ", or NULL.", call. = FALSE)
            }
      }

      # Validate theta (kernel scale). Optional; single positive finite number.
      if (!is.null(theta)) {
            if (!is.numeric(theta) || length(theta) != 1L ||
                !is.finite(theta) || theta <= 0) {
                  stop("`theta` must be a single positive number, or NULL.",
                       call. = FALSE)
            }
      }

      # Validate max (upper threshold). Optional. A positive numeric vector; Inf
      # is allowed and means "no bound" for that family (equivalent to omitting
      # the threshold, but accepted so callers can pass max = Inf explicitly).
      # Length 1 is the usual case (a single distance radius). For the climate
      # family, a longer vector specifies per-variable absolute-difference
      # thresholds; the required length (1 or the number of climate variables)
      # is checked downstream against the data. The geographic family uses a
      # single radius. We validate only shape-independent properties here.
      if (!is.null(max)) {
            if (!is.numeric(max) || length(max) < 1L ||
                any(is.na(max)) || any(max <= 0)) {
                  stop("`max` must be a positive number (Inf allowed; or, for the ",
                       "climate family, a vector of positive per-variable ",
                       "thresholds), or NULL.", call. = FALSE)
            }
      }

      # Validate min (lower threshold). Optional. Unlike max, min is always a
      # single positive FINITE scalar: an infinite floor would exclude every
      # candidate, and a per-variable min is not defined (the feature is
      # geographic-only in the engine for now; the environmental-family
      # rejection happens downstream in analog_search()). A per-variable box
      # min would also raise an ambiguous L-infinity-vs-L2 semantics question we
      # deliberately avoid here.
      if (!is.null(min)) {
            if (!is.numeric(min) || length(min) != 1L ||
                !is.finite(min) || min <= 0) {
                  stop("`min` must be a single positive finite number, or NULL.",
                       call. = FALSE)
            }
      }

      structure(
            list(weight = weight, theta = theta, max = max, min = min),
            class = "analog_kernel"
      )
}


#' @export
print.analog_kernel <- function(x, ...) {
      w <- x$weight %||% "uniform"
      cat("<analog_kernel>\n")
      cat(sprintf("  weight: %s\n", w))
      if (!identical(w, "uniform")) {
            cat(sprintf("  theta:  %s\n",
                        if (is.null(x$theta)) "1 (default)" else format(x$theta)))
      }
      cat(sprintf("  max:    %s\n",
                  if (is.null(x$max)) "none"
                  else paste(format(x$max), collapse = ", ")))
      cat(sprintf("  min:    %s\n",
                  if (is.null(x$min)) "none"
                  else format(x$min)))
      invisible(x)
}
