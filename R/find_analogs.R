#' Find climate analogs
#'
#' Thin, user-facing wrapper that preserves your original parameter names and
#' data structures (\code{focal}, \code{ref} as matrix/data.frame with x,y and
#' climate columns, or \code{SpatRaster}). Internals normalize via
#' \code{format_data()} and call the C++ core. Adds \code{metric}, \code{weight},
#' and \code{report_dist} per the roadmap, plus \code{geo = "auto"} to detect
#' lon/lat vs projected.
#'
#' @param focal Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers.
#' @param ref   Matrix/data.frame with columns x, y, and climate variables,
#'   OR SpatRaster with climate variable layers.
#' @param max_dist Maximum geographic distance constraint in km (NULL = no constraint).
#' @param max_clim Maximum climate constraint. Either a scalar Euclidean radius in
#'   climate space, or a length-p vector of per-variable absolute bands (NULL = no constraint).
#' @param weight Weighting function: "uniform", "inverse_clim", or "inverse_dist".
#' @param fun Aggregation: "nmax", "sum", "mean", "count", or "all".
#' @param n Number of analogs to return (required if fun = "nmax").
#' @param metric Climate distance metric (currently "euclidean"; others TBD per roadmap).
#' @param report_dist Logical; if TRUE, include distance columns in outputs where applicable.
#' @param geo One of "auto" (default), "lonlat", or "projected". "auto" detects from coordinates.
#'
#' @return
#' For fun = "nmax" or "all": data.frame with columns
#'   focal_index, focal_x, focal_y, analog_index, analog_x, analog_y,
#'   and optionally clim_dist, geog_dist, weight (when report_dist = TRUE).
#'
#' For fun = "sum", "mean", or "count": data.frame with
#'   focal_index, focal_x, focal_y, value (aggregated by the chosen weight).
#'
#' @export
find_analogs <- function(
            focal,
            ref,
            max_dist   = NULL,
            max_clim   = NULL,
            weight     = "inverse_clim",
            fun        = "nmax",
            n          = NULL,
            metric     = "euclidean",
            report_dist = TRUE,
            geo        = c("auto", "lonlat", "projected")
) {
      # ---- validate args -----------------------------------------------------
      geo <- match.arg(geo)
      valid_weights <- c("uniform", "inverse_clim", "inverse_dist")
      if (!weight %in% valid_weights) stop("weight must be one of: ", paste(valid_weights, collapse=", "))
      valid_funs <- c("nmax", "sum", "mean", "count", "all")
      if (!fun %in% valid_funs) stop("fun must be one of: ", paste(valid_funs, collapse=", "))
      metric <- tolower(metric)
      if (!metric %in% c("euclidean")) {
            stop("metric='", metric, "' not yet supported; currently only 'euclidean' is implemented in the core")
      }
      if (fun == "nmax" && is.null(n)) stop("n must be specified when fun = 'nmax'")
      if (!is.null(n)) n <- as.integer(n)

      # ---- normalize inputs via user's helpers -------------------------------
      focal_mm <- .format_data(focal)
      ref_mm   <- .format_data(ref)
      if (ncol(focal_mm$climate) != ncol(ref_mm$climate)) stop("number of climate variables must match between focal and ref")

      # geo mode
      geo_mode <- switch(geo,
                         auto      = .detect_geo(focal_mm$coords, ref_mm$coords),
                         lonlat    = "lonlat",
                         projected = "projected"
      )

      # constraints
      radius_km <- if (is.null(max_dist)) Inf else as.numeric(max_dist)[1L]
      climate_band <- if (is.null(max_clim)) Inf else max_clim
      if (length(climate_band) == 1L) {
            climate_band <- as.numeric(climate_band)[1L]
      } else if (length(climate_band) == ncol(ref_mm$climate)) {
            climate_band <- as.numeric(climate_band)
      } else {
            stop("`max_clim` must be NULL, a scalar, or length ncol(climate) vector")
      }

      # choose k behavior to match previous API
      k <- if (identical(fun, "nmax")) n else 0L

      # ---- call C++ core -----------------------------------------------------
      res <- .Call(
            `_analogs_find_analogs_core`,
            focal_mm$climate, ref_mm$climate,
            focal_mm$coords,  ref_mm$coords,
            NA_character_,             # mode (reserved)
            as.integer(k),
            climate_band,
            as.numeric(radius_km),
            geo_mode,
            FALSE                      # compact_bins (reserved)
      )

      # res is list of integer vectors (1-based analog indices per focal)

      # ---- post-process into requested shape --------------------------------
      if (fun %in% c("nmax", "all")) {
            out <- .emit_pairs(res, focal_mm, ref_mm, report_dist = report_dist, weight = weight, geo_mode = geo_mode)
            return(out)
      }

      if (fun %in% c("sum", "mean", "count")) {
            out <- .emit_aggregates(res, focal_mm, ref_mm, fun = fun, weight = weight, geo_mode = geo_mode)
            return(out)
      }

      stop("unreachable")
}

# ---- helpers --------------------------------------------------------------

.detect_geo <- function(focal_xy, ref_xy) {
      # rudimentary heuristic: lon/lat if both appear within plausible bounds
      rng <- range(c(focal_xy[,1], ref_xy[,1], focal_xy[,2], ref_xy[,2]), na.rm = TRUE)
      # More robust: check both axes separately
      lon_rng <- range(c(focal_xy[,1], ref_xy[,1]), na.rm = TRUE)
      lat_rng <- range(c(focal_xy[,2], ref_xy[,2]), na.rm = TRUE)
      if (all(is.finite(c(lon_rng, lat_rng))) &&
          lon_rng[1] >= -180 && lon_rng[2] <= 180 &&
          lat_rng[1] >=  -90 && lat_rng[2] <=  90) {
            "lonlat"
      } else {
            "projected"
      }
}

.haversine_km <- function(xy1, xy2) {
      # xy1, xy2 are 2-column matrices; computes pairwise between matching rows
      R <- 6371.0088
      to_rad <- pi/180
      lon1 <- xy1[,1] * to_rad; lat1 <- xy1[,2] * to_rad
      lon2 <- xy2[,1] * to_rad; lat2 <- xy2[,2] * to_rad
      dlon <- lon2 - lon1; dlat <- lat2 - lat1
      sdlat <- sin(0.5*dlat); sdlon <- sin(0.5*dlon)
      a <- sdlat*sdlat + cos(lat1)*cos(lat2)*sdlon*sdlon
      2 * R * asin(pmin(1, sqrt(a)))
}

# Build long table of focal-analog pairs, with optional distances and weights
.emit_pairs <- function(res, focal_mm, ref_mm, report_dist, weight, geo_mode) {
      n_f <- nrow(focal_mm$coords)
      rows <- vector("list", n_f)
      for (i in seq_len(n_f)) {
            idx <- res[[i]]
            if (length(idx) == 0L) { rows[[i]] <- NULL; next }
            fx <- rep(focal_mm$coords[i,1], length(idx))
            fy <- rep(focal_mm$coords[i,2], length(idx))
            ax <- ref_mm$coords[idx,1]
            ay <- ref_mm$coords[idx,2]

            df <- data.frame(
                  focal_index  = rep.int(i, length(idx)),
                  focal_x = fx, focal_y = fy,
                  analog_index = idx,
                  analog_x = ax, analog_y = ay,
                  stringsAsFactors = FALSE
            )

            if (isTRUE(report_dist) || !identical(weight, "uniform")) {
                  # climate distances
                  v_i <- matrix(focal_mm$climate[i,], nrow = 1)
                  v_j <- ref_mm$climate[idx,,drop=FALSE]
                  clim_d <- sqrt(rowSums((t(t(v_j) - as.numeric(v_i)))^2))

                  # geographic distances
                  if (geo_mode == "lonlat") {
                        geog_d <- .haversine_km(cbind(fx, fy), cbind(ax, ay))
                  } else {
                        geog_d <- sqrt((ax - fx)^2 + (ay - fy)^2)
                  }

                  if (isTRUE(report_dist)) {
                        df$clim_dist <- clim_d
                        df$geog_dist <- geog_d
                  }

                  # weights if requested
                  if (identical(weight, "inverse_clim")) {
                        df$weight <- 1/(clim_d + 1e-12)
                  } else if (identical(weight, "inverse_dist")) {
                        df$weight <- 1/(geog_d + 1e-12)
                  } else {
                        df$weight <- 1.0
                  }
            }

            rows[[i]] <- df
      }
      if (all(vapply(rows, is.null, logical(1)))) return(utils::head(data.frame(), 0L))
      do.call(rbind, rows)
}

.emit_aggregates <- function(res, focal_mm, ref_mm, fun, weight, geo_mode) {
      n_f <- nrow(focal_mm$coords)
      out <- data.frame(
            focal_index = seq_len(n_f),
            focal_x = focal_mm$coords[,1],
            focal_y = focal_mm$coords[,2],
            value = NA_real_,
            stringsAsFactors = FALSE
      )
      for (i in seq_len(n_f)) {
            idx <- res[[i]]
            if (length(idx) == 0L) {
                  out$value[i] <- if (fun == "count") 0 else 0
                  next
            }
            # distances
            v_i <- matrix(focal_mm$climate[i,], nrow = 1)
            v_j <- ref_mm$climate[idx,,drop=FALSE]
            clim_d <- sqrt(rowSums((t(t(v_j) - as.numeric(v_i)))^2))
            if (geo_mode == "lonlat") {
                  geog_d <- .haversine_km(cbind(rep(focal_mm$coords[i,1], length(idx)), rep(focal_mm$coords[i,2], length(idx))),
                                          ref_mm$coords[idx,1:2,drop=FALSE])
            } else {
                  fx <- focal_mm$coords[i,1]; fy <- focal_mm$coords[i,2]
                  ax <- ref_mm$coords[idx,1]; ay <- ref_mm$coords[idx,2]
                  geog_d <- sqrt((ax - fx)^2 + (ay - fy)^2)
            }

            # weights
            w <- switch(weight,
                        uniform      = rep(1, length(idx)),
                        inverse_clim = 1/(clim_d + 1e-12),
                        inverse_dist = 1/(geog_d + 1e-12)
            )

            if (fun == "count") {
                  out$value[i] <- sum(w*0 + 1)  # count of matches
            } else if (fun == "sum") {
                  out$value[i] <- sum(w)
            } else if (fun == "mean") {
                  out$value[i] <- mean(w)
            } else {
                  stop("unrecognized aggregate fun: ", fun)
            }
      }
      out
}



# Helper: select x,y by name if present; otherwise use first 2 cols
.select_xy_climate <- function(obj) {
  nm <- colnames(obj)
  if (!is.null(nm) && all(c("x","y") %in% nm[1:max(2, length(nm))])) {
    x_idx <- match("x", nm)
    y_idx <- match("y", nm)
    xy_idx <- c(x_idx, y_idx)
  } else if (!is.null(nm) && all(c("x","y") %in% nm)) {
    xy_idx <- match(c("x","y"), nm)
  } else {
    xy_idx <- 1:2
  }

  coords  <- as.matrix(obj[, xy_idx, drop = FALSE])
  climate <- as.matrix(obj[, setdiff(seq_len(ncol(obj)), xy_idx), drop = FALSE])

  storage.mode(coords)  <- "double"
  storage.mode(climate) <- "double"

  if (ncol(coords) != 2L) stop("coords must have exactly 2 columns (x,y).")
  if (ncol(climate) < 1L) stop("no climate variable columns found after x,y.")

  list(coords = coords, climate = climate)
}


#' Normalize input (SpatRaster OR data.frame/matrix) to coords/climate matrices
#' @keywords internal
.format_data <- function(r) {
  if (inherits(r, "SpatRaster")) {
    df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
    .select_xy_climate(df)

  } else if (is.matrix(r) || is.data.frame(r)) {
    .select_xy_climate(r)

  } else {
    stop("Unsupported input type. Provide a data.frame/matrix with columns x,y,... or a SpatRaster.")
  }
}
