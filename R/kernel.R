#' Define the distance treatment for one analog dimension
#'
#' Constructs a specification of how a single distance family — climate or
#' geography — is filtered and weighted in an analog search. A `kernel`
#' bundles the three per-family choices: the hard distance threshold (`max`),
#' the weighting kernel shape (`weight`), and that kernel's scale parameter
#' (`theta`). Pass one to the `clim` argument of [analog_search()] and/or one
#' to the `geog` argument.
#'
#' The overall kernel weight for a candidate is the product of the two
#' families' weights, so the families are specified independently and may use
#' different shapes (e.g. an inverse-distance climate kernel together with a
#' Gaussian geographic kernel). A family with `weight = "uniform"` (or `NULL`)
#' contributes a constant weight of 1, i.e. it filters (if `max` is set) but
#' does not down-weight by distance.
#'
#' All three components are optional. Which combinations are valid depends on
#' the operation and is checked by [analog_search()] downstream (for example,
#' climate velocity requires a climate `max`; a weighted statistic requires a
#' non-uniform `weight` on at least one family). A bare `NULL` passed as the
#' `clim` or `geog` argument is equivalent to `kernel()` with all components
#' unset: no threshold and no weighting for that family.
#'
#' @param weight Kernel shape for this family. One of `"uniform"` (constant
#'   weight 1; the default when `NULL`), `"gaussian"`
#'   (`exp(-d^2 / (2 theta^2))`), or `"inverse"` (`1 / (1 + d / theta)`, a
#'   heavy-tailed inverse-distance kernel). The family (climate vs geographic)
#'   is determined by whether the kernel is passed as `clim` or `geog`, so the
#'   shape name is unqualified.
#' @param theta Scale parameter for the `weight` kernel. For `"gaussian"` it is
#'   the bandwidth (sigma); for `"inverse"` it is the half-weight distance (the
#'   weight is 1/2 at `d = theta`). Ignored for `"uniform"`. Defaults to `NULL`,
#'   which lets downstream code apply a default of 1. See [kernel_params()] for
#'   help choosing a value calibrated to a target coverage fraction.
#' @param max Hard distance threshold for this family: candidates beyond `max`
#'   (in this family's distance) are excluded. `NULL` (default) means no
#'   threshold. Usually a single radius. For the climate family, `max` may also
#'   be a vector of per-variable absolute-difference thresholds (length equal to
#'   the number of climate variables); the geographic family uses a single
#'   radius. Supplied to the search as `max_clim` / `max_geog`.
#'
#' @return An object of class `"analog_kernel"`: a list with elements
#'   `weight`, `theta`, and `max` (each possibly `NULL`).
#'
#' @seealso [analog_search()], [kernel_params()]
#'
#' @examples
#' # Climate: keep analogs within 2 climate-distance units, Gaussian-weighted
#' kernel(weight = "gaussian", theta = 0.5, max = 2)
#'
#' # Geography: hard 100 km cutoff, no distance weighting (uniform)
#' kernel(max = 100)
#'
#' # Inverse-distance climate weighting, no hard cutoff
#' kernel("inverse", theta = 1)
#'
#' @export
kernel <- function(weight = NULL, theta = NULL, max = NULL) {

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

      # Validate max (threshold). Optional. A positive numeric vector; Inf is
      # allowed and means "no bound" for that family (equivalent to omitting the
      # threshold, but accepted so callers can pass max = Inf explicitly).
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

      structure(
            list(weight = weight, theta = theta, max = max),
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
      invisible(x)
}
