#' Kernel parameter recommendations
#'
#' Returns recommended values for the bandwidth (`theta`) and hard
#' truncation distance (`max`) of a kernel operating on Euclidean
#' distances in `d`-dimensional space. Gives theoretical answers to the
#' questions, "How big should `theta` be in order for my kernel to capture
#' a given fraction of pool sites for the typical focal site?" and "How
#' big should `max_clim` or `max_geog` be to truncate only a given
#' percentage of kernel weight for the typical focal site?" Accounts for
#' the effects of analog space multidimensionality on pairwise distance
#' distributions, which can result in one-dimensional intuitions being
#' incorrect. Supports the kernel types and data distributions used by
#' the `analogs` package.
#'
#' Either `fraction` or `theta` should be provided. When `fraction` is
#' given, the function returns the `theta` that calibrates the kernel
#' to capture, on average, that fraction of pool sites (where partial
#' capture is in proportion to kernel weight). Switching kernel shapes
#' at fixed `fraction` holds the expected total kernel weight constant,
#' so weighted aggregate statistics (e.g. `sum_weights`) remain
#' comparable across kernels.
#'
#' For uniform data, `fraction` is not meaningful because "fraction of
#' space" depends on landscape extent; `theta` must be supplied directly
#' (e.g. a dispersal-derived bandwidth for a geographic kernel).
#'
#' When `loss` is specified, the function additionally returns `max`:
#' the truncation distance beyond which less than `loss` of aggregate
#' kernel weight is discarded. Useful for computational efficiency.
#'
#' Recommendations are averages over the distribution of focal cells;
#' specific focal cells experience effective neighborhoods that vary
#' around these averages, with cells in dense climate regions seeing
#' more neighbors than cells in sparse regions.
#'
#' @param fraction Target fraction of pool sites captured (in
#'   weight-proportional terms) by the kernel for the typical focal.
#'   Use this OR `theta`. Requires `data_dist = "mvn"`.
#' @param theta Bandwidth value. Use this OR `fraction`. For Gaussian
#'   kernels, this is the standard bandwidth parameter; for uniform
#'   kernels, this is the cutoff radius (also returned as `max`); for
#'   inverse-distance kernels, this is the half-weight scale of the
#'   reparameterized kernel `1 / (1 + d / theta)` (weight is 1/2 at
#'   `d = theta`).
#' @param d Dimensionality of the space (e.g., number of climate
#'   variables after Mahalanobis transformation, or 2 for geographic).
#' @param loss Fraction of aggregate kernel weight to discard at the
#'   truncation distance `max`. If `NULL` (default), `max` is not
#'   computed.
#' @param kernel One of `"gaussian"` (default), `"uniform"`, or
#'   `"inverse_distance"`.
#' @param data_dist Distribution of cells in space. Either `"mvn"`
#'   (multivariate standard normal; default; appropriate for
#'   Mahalanobis-transformed climate data) or `"uniform"` (appropriate
#'   for geographic space).
#' @return A named list. For Gaussian and inverse-distance kernels:
#'   element `theta`, and `max` if `loss` is specified. For uniform
#'   kernels: element `max` (the single cutoff radius, which serves as
#'   both bandwidth and truncation distance; supplied in the `analogs`
#'   package as `max_clim` or `max_geog`).
#' @examples
#' # Climate kernel: niche fraction of 5% in 4 climate variables
#' kernel_params(fraction = 0.05, d = 4, loss = 0.01)
#'
#' # Geographic kernel: 500 km dispersal-based bandwidth
#' kernel_params(theta = 500, d = 2, data_dist = "uniform", loss = 0.01)
#'
#' # Switching kernels at fixed niche fraction (matched expected weight)
#' kernel_params(fraction = 0.05, d = 4, kernel = "gaussian")
#' kernel_params(fraction = 0.05, d = 4, kernel = "uniform")
#' kernel_params(fraction = 0.05, d = 4, kernel = "inverse_distance")
#' @export
kernel_params <- function(fraction = NULL,
                          theta = NULL,
                          d,
                          loss = NULL,
                          kernel = c("gaussian", "uniform", "inverse_distance"),
                          data_dist = c("mvn", "uniform")) {

      kernel <- match.arg(kernel)
      data_dist <- match.arg(data_dist)

      has_fraction <- !is.null(fraction)
      has_theta <- !is.null(theta)

      if (has_fraction == has_theta) {
            stop("Provide exactly one of `fraction` or `theta`.")
      }
      if (data_dist == "uniform" && has_fraction) {
            stop("For `data_dist = \"uniform\"`, `theta` must be provided directly.")
      }
      if (has_fraction && (fraction <= 0 || fraction >= 1)) {
            stop("`fraction` must be in (0, 1).")
      }
      if (!is.null(loss) && (loss <= 0 || loss >= 1)) {
            stop("`loss` must be in (0, 1).")
      }
      if (d < 1) {
            stop("`d` must be at least 1.")
      }

      # Chi distribution density (radial distance from origin under MVN)
      chi_density <- function(r) {
            (r^(d - 1) * exp(-r^2 / 2)) / (2^(d / 2 - 1) * gamma(d / 2))
      }

      # ---- Compute theta from fraction ----
      if (!has_theta) {
            if (kernel == "gaussian") {
                  theta <- sqrt(fraction^(2 / d) / (1 - fraction^(2 / d)))
            } else if (kernel == "uniform") {
                  theta <- sqrt(qchisq(fraction, df = d))
            } else if (kernel == "inverse_distance") {
                  # Reparameterized inverse kernel: w(r) = 1 / (1 + r/theta).
                  # Calibrate theta so the expected weight over the chi
                  # distribution equals `fraction`. The integral converges (the
                  # 1 + bounds the integrand at r = 0) and is monotonically
                  # increasing in theta (larger theta -> flatter kernel ->
                  # captures more), so uniroot finds a unique root.
                  target <- function(th) {
                        w <- .kernel_weight_fn("inverse", th)
                        integrate(function(r) chi_density(r) * w(r),
                                  lower = 0, upper = Inf)$value - fraction
                  }
                  theta <- uniroot(target, interval = c(1e-8, 1e8))$root
            }
      }

      # ---- Compute max from loss ----
      max_val <- NULL
      if (!is.null(loss) && kernel != "uniform") {
            if (kernel == "gaussian") {
                  if (data_dist == "mvn") {
                        tau <- theta / sqrt(1 + theta^2)
                        max_val <- tau * sqrt(qchisq(1 - loss, df = d))
                  } else {
                        max_val <- theta * sqrt(qchisq(1 - loss, df = d))
                  }
            } else if (kernel == "inverse_distance") {
                  if (data_dist == "mvn") {
                        w <- .kernel_weight_fn("inverse", theta)
                        total <- integrate(
                              function(r) chi_density(r) * w(r),
                              lower = 0, upper = Inf
                        )$value
                        max_val <- uniroot(
                              function(R) {
                                    integrate(function(r) chi_density(r) * w(r),
                                              lower = R, upper = Inf)$value -
                                          loss * total
                              },
                              interval = c(1e-10, 100)
                        )$root
                  } else {
                        stop("Inverse-distance with uniform data: ",
                             "max depends on landscape extent; specify manually.")
                  }
            }
      }

      # Uniform kernels have only one parameter (the cutoff radius),
      # exposed in the analogs package as `max_clim` / `max_geog`.
      # Return that single parameter as `max` rather than as `theta`.
      if (kernel == "uniform") {
            return(list(max = theta))
      }

      result <- list(theta = theta)
      if (!is.null(max_val)) result$max <- max_val
      result
}
