## Code to create the package hex logo, man/figures/logo.png.
##
## The logo shows real analog search output from the example rasters: a focal point and its
## forward analogs (where its historic climate occurs in the future) and backward analogs (where
## its future climate occurs today). Each analog field is the Gaussian kernel weight of every
## pool cell within reach, drawn as a single filled contour at a fixed weight. Lines connect the focal point to a weight-proportional
## sample of its analogs. The background is the historic climate itself, as a dark bivariate
## map of the two example variables (CWD and AET).
##
## Run from the package root. The title uses Helvetica Neue Condensed Bold (macOS), falling back
## to Helvetica Neue Bold or Nimbus Sans Narrow Bold if it isn't installed. Everything else is
## deterministic, so reruns give the same logo.
##
## The image bleeds to the hexagon's edges, so instead of hexSticker (which doesn't clip its
## subplot to the hexagon), the plot is rendered with ragg and masked to the hexagon with magick,
## giving a transparent background outside it.

library(ggplot2)
devtools::load_all()

out_file <- tempfile(fileext = ".png")   # full-resolution render; use_logo() resizes it
dpi <- 600

# design parameters ---------------------------------------------------------------------------

site <- cbind(-108, 37)                  # focal point (lon, lat)
theta <- 0.1                             # Gaussian climate kernel bandwidth (scaled units)
n_lines <- 40                            # analog lines per direction
R <- 520                                 # hexagon circumradius, km
inch_per_R <- 1.176                      # print scale: inches of logo per R; sets the
# sticker size (sqrt(3) x 2 R) that sizes below assume
glow_dist <- 2.9                         # analog fields enclose analogs within this many
# bandwidths (theta) of the focal's climate
glow_alpha <- 0.5                        # analog field fill opacity
seed <- 7

# palette. The background is a dark bivariate climate map, with hue from CWD and luminance
# from AET. Lines and their endpoints are tints of the glow colors.
col <- list(border = "#07130F",
            fwd = "#FFC23D",                 # forward analogs
            bwd = "#B59BFF",                 # backward analogs
            bg_hue = c(160, 200),
            bg_chroma = 14,
            bg_lum = c(4, 20))
tint <- function(x, f) colorRampPalette(c(x, "white"))(101)[round(100 * f) + 1]
line_tint <- 0.25                        # mix toward white, for lines
point_tint <- 0.5                        # mix toward white, for line endpoints

# title font: first available of these faces
font_choices <- list(c("Helvetica Neue", "Condensed Bold"),
                     c("Helvetica Neue", "Bold"),
                     c("Nimbus Sans Narrow", "Bold"))

# analog search output ------------------------------------------------------------------------

clim <- example_rasters()
lon0 <- site[1]
lat0 <- site[2]

# local km coordinates around the focal point, with north up and true shapes at its latitude
km_lat <- 110.57
km_lon <- 111.32 * cos(lat0 * pi / 180)
to_km <- function(d, x = "x", y = "y"){
      d$X <- (d[[x]] - lon0) * km_lon
      d$Y <- (d[[y]] - lat0) * km_lat
      d
}

# all analogs within reach of the hexagon, with their Gaussian climate kernel weights
analogs <- function(from, pool){
      focal <- data.frame(x = lon0, y = lat0, terra::extract(from, site))
      d <- analog_search(focal, pool, env = kernel(max = 3 * theta),
                         geog = kernel(max = 1.2 * R))
      d$w <- exp(-d$env_dist^2 / (2 * theta^2))
      to_km(d, "analog_x", "analog_y")
}
fwd <- analogs(clim$historic, clim$future)
bwd <- analogs(clim$future, clim$historic)

# analog weight surface, with bilinear disaggregation so its contour is smooth
glow_level <- exp(-glow_dist^2 / 2)
glow <- function(d){
      r <- terra::rast(clim$historic, nlyrs = 1, vals = 0)
      r[terra::cellFromXY(r, cbind(d$analog_x, d$analog_y))] <- d$w
      r <- terra::crop(r, terra::ext(c(lon0 + c(-1.3, 1.3) * R / km_lon,
                                       lat0 + c(-1.3, 1.3) * R / km_lat)))
      r <- terra::disagg(r, 4, method = "bilinear")
      d <- to_km(as.data.frame(r, xy = TRUE))
      names(d)[3] <- "v"
      d
}
glow_fwd <- glow(fwd)
glow_bwd <- glow(bwd)

# lines to a weight-proportional sample of analogs, among those well inside the hexagon
in_hex <- function(X, Y, r) abs(X) <= sqrt(3) / 2 * r & abs(Y) <= r - abs(X) / sqrt(3)
set.seed(seed)
lines <- function(d){
      d <- d[in_hex(d$X, d$Y, 0.85 * R), ]
      d[sample(nrow(d), n_lines, prob = d$w), ]
}
lines_fwd <- lines(fwd)
lines_bwd <- lines(bwd)

# background: historic climate, rank-scaled within the window
bg <- terra::crop(clim$historic, terra::ext(c(lon0 + c(-1.2, 1.2) * R / km_lon,
                                              lat0 + c(-1.2, 1.2) * R / km_lat)))
bg <- to_km(as.data.frame(bg, xy = TRUE))
rank01 <- function(v) (rank(v) - 1) / (length(v) - 1)
bg$fill <- hcl(h = col$bg_hue[1] + diff(col$bg_hue) * rank01(bg$CWD),
               c = col$bg_chroma,
               l = col$bg_lum[1] + diff(col$bg_lum) * rank01(bg$AET))

# title font ----------------------------------------------------------------------------------

fonts <- systemfonts::system_fonts()
family <- NULL
for(fc in font_choices){
      i <- which(fonts$family == fc[1] & fonts$style == fc[2])
      if(length(i)){
            systemfonts::register_font("analogs-logo",
                                       plain = list(fonts$path[i[1]], fonts$index[i[1]]))
            family <- "analogs-logo"
            message("title font: ", fc[1], " ", fc[2])
            break
      }
}
if(is.null(family)) stop("none of the title fonts in `font_choices` is installed")

# plot ----------------------------------------------------------------------------------------

hex_angles <- (90 + 60 * 0:5) * pi / 180  # pointy-top

glow_layers <- function(d, color){
      list(geom_contour_filled(data = d, aes(X, Y, z = v), breaks = c(glow_level, 2),
                               fill = color, alpha = glow_alpha),
           geom_contour(data = d, aes(X, Y, z = v), breaks = glow_level,
                        color = color, linewidth = 0.15))
}

line_layers <- function(d, color){
      list(geom_curve(data = d, aes(0, 0, xend = X, yend = Y), curvature = 0.2,
                      color = tint(color, line_tint), alpha = 0.3, linewidth = 0.22),
           geom_point(data = d, aes(X, Y), color = tint(color, point_tint),
                      alpha = 0.8, size = 0.25))
}

# vignette mask and distance contours
bg$dist <- with(bg, sqrt((X)^2 + (Y)^2))
bg$k <- with(bg, dnorm(dist, sd = 250))

p <- ggplot() +
      geom_raster(data = bg, aes(X, Y, fill = fill)) +
      scale_fill_identity() +
      glow_layers(glow_bwd, col$bwd) +
      glow_layers(glow_fwd, col$fwd) +
      line_layers(lines_bwd, col$bwd) +
      line_layers(lines_fwd, col$fwd) +
      geom_contour(data = bg, aes(X, Y, z = k),
                   color = "white", linewidth = .1, alpha = .5) +
      geom_raster(data = bg, aes(X, Y, alpha = k), fill = "black") +
      scale_alpha_continuous(range = 1:0) +
      annotate("point", 0, 0, size = 1.6, color = "white") +
      annotate("text", 0, -0.3 * R, #-0.47 * R,
               label = "analogs", color = "white",
               size = 7.4, family = family) +
      coord_equal(xlim = c(-1, 1) * sqrt(3) / 2 * R, ylim = c(-1, 1) * R, expand = FALSE) +
      theme_void() +
      theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0),
            plot.background = element_blank(), panel.background = element_blank())

# render, mask to the hexagon, and add a border -----------------------------------------------

width <- sqrt(3) * inch_per_R             # hexagon width sqrt(3) R, height 2 R
height <- 2 * inch_per_R

# The border is a filled ring of uniform width between two hexagons, drawn entirely inside the
# image, so its width is the same on every side and nothing is clipped at the edges. The outer
# hexagon is inset slightly so its points aren't cut off by the image boundary.
border_in <- 0.05                         # border width, inches (perpendicular to each side)
r_out <- inch_per_R * (1 - 0.004)         # outer hexagon circumradius, inches
r_in <- r_out - border_in / (sqrt(3) / 2) # inner: apothem smaller by exactly border_in
hex_npc <- function(r) list(x = 0.5 + r * cos(hex_angles) / width,
                            y = 0.5 + r * sin(hex_angles) / height)
outer <- hex_npc(r_out)
inner <- hex_npc(r_in)

render <- function(draw){
      f <- tempfile(fileext = ".png")
      ragg::agg_png(f, width = width, height = height, units = "in", res = dpi,
                    background = "transparent")
      grid::grid.newpage()
      draw()
      invisible(dev.off())
      magick::image_read(f)
}
art <- render(function() grid::grid.draw(ggplotGrob(p)))
mask <- render(function() grid::grid.polygon(outer$x, outer$y,
                                             gp = grid::gpar(col = NA, fill = "black")))
border <- render(function() grid::grid.path(c(outer$x, inner$x), c(outer$y, inner$y),
                                            id = rep(1:2, each = 6), rule = "evenodd",
                                            gp = grid::gpar(col = NA, fill = col$border)))
logo <- magick::image_composite(art, mask, operator = "DstIn")
logo <- magick::image_composite(logo, border, operator = "Over")
magick::image_write(logo, out_file)

# use_logo() resizes to fit `geometry` (default 240x278); twice that keeps fine detail sharp on
# high-DPI screens, since pkgdown and the README display the logo much smaller than this
usethis::use_logo(out_file, geometry = "480x556")
