
library(terra)
library(tidyverse)
library(patchwork)

scale_rast <- function(x) (x - mean(values(x), na.rm = T)) / sd(values(x), na.rm = T)



bb <- ext(-117, -104, 42, 49)

t1 <- rast("/Volumes/T7/CHELSA/v2/raw/CHELSA_bio1_1981-2010_V.2.1.tif") %>% crop(bb)
t2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
t2 <- t2[grepl("2071", t2) & grepl("ssp370_tas_", t2)] %>% rast() %>% crop(bb) %>% mean()
t <- c(t1, t2) %>% scale_rast()

p1 <- rast("/Volumes/T7/CHELSA/v2/raw/CHELSA_bio12_1981-2010_V.2.1.tif") %>% crop(bb)
p2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
p2 <- p2[grepl("2071", p2) & grepl("ssp370_pr_", p2)] %>% rast() %>% crop(bb) %>% mean() %>% "*"(12)
p <- c(p1, p2) %>% log() %>% scale_rast()

clim1 <- c(t[[1]], p[[1]])
clim2 <- c(t[[2]], p[[2]])




n_focals <- 16

set.seed(123)
f_crds <- crds(clim1)[sample(ncell(clim1), n_focals),, drop = F] # -112.1793 48.79569
f_clim <- terra::extract(clim1, f_crds)
focal <- as.matrix(cbind(f_crds, f_clim))
colnames(focal) <- c("x", "y", "t", "p")



# a %>%
#       ggplot(aes(clim_dist, geog_dist + .5)) +
#       facet_wrap(~focal_index, nrow = 5) +
#       geom_point(color = "darkblue", size = .1) +
#       scale_y_continuous(expand = c(0,0)) +
#       scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
#       theme_bw()


pareto_optimal <- function(
            x # a two-column matrix of clim_dist and geog_dist values for a focal site
){
      x <- as.matrix(x)
      idx <- 1:nrow(x)
      x <- cbind(x, idx)

      frontier_idx <- vector()
      f <- x

      while(nrow(f) > 0){
            i <- which.min(f[,2])
            frontier_idx <- c(frontier_idx, f[i,3])
            f <- f[f[,1] < f[i,1], , drop = F]
      }

      idx %in% frontier_idx
}

on_pareto_front <- function(x) {
      x <- as.matrix(x)
      stopifnot(ncol(x) == 2)

      # Columns: clim_dist, geog_dist
      # Order by geog_dist (col 2), then clim_dist (col 1)
      o <- order(x[, 2], x[, 1])
      xs <- x[o, , drop = FALSE]

      best_clim <- Inf
      on_front_sorted <- logical(nrow(xs))

      for (k in seq_len(nrow(xs))) {
            clim_k <- xs[k, 1]
            if (clim_k <= best_clim) {
                  best_clim <- clim_k
                  on_front_sorted[k] <- TRUE
            }
      }

      # Map back to original order
      on_front <- logical(nrow(x))
      on_front[o] <- on_front_sorted
      on_front
}


pareto_front <- function(x, clim_tol = .01){

      map_dfr(unique(x$focal_index), function(i){
            xi <- filter(x, focal_index == i)

            x_hull <- xi %>%
                  # filter(pareto_optimal(select(., clim_dist, geog_dist))) %>%
                  filter(on_pareto_front(select(., clim_dist, geog_dist))) %>%
                  filter(geog_dist <= min(geog_dist[clim_dist <= min(clim_dist + clim_tol)])) %>%
                  mutate(solution = geog_dist/max(geog_dist) - clim_dist/max(clim_dist)) %>%
                  arrange(solution)

            return(x_hull)
      })
}

pareto_paths <- function(pf, x, ref, combine = TRUE){

      pf$color <- colors3d::colors2d(select(pf, clim_dist, geog_dist),
                                     colors = c("magenta", "yellow", "green", "cyan"),
                                     ytrans = "rank"
      )

      x <- x %>% filter(clim_dist <= max(pf$clim_dist),
                  geog_dist <= max(pf$geog_dist))

      pf_focal <- pf %>% select(focal_index, focal_x, focal_y) %>%
            distinct()

      p_curves <- pf %>%
            ggplot(aes(geog_dist, clim_dist,
                       group = focal_index, color = I(color))) +
            geom_path(linewidth = .25) +
            geom_point(size = .25) +
            theme_dark() +
            theme(panel.background = element_rect(fill = "gray30")) +
            # scale_x_sqrt() +
            labs(x = "migration distance (km)",
                 y = "climate exposure (sigma)")

      p_dist <- ggplot(pf, aes(geog_dist, clim_dist,
                               group = focal_index, color = I(color))) +
            facet_wrap(~focal_index,
                       ncol = round(sqrt(length(unique(pf$focal_index)))),
                       scales = "free_x") +
            geom_point(data = x,
                       color = "black", size = .25) +
            geom_path(size = .25) +
            geom_point(size = .25) +
            theme_minimal() +
            theme(panel.background = element_rect(fill = "gray60"),
                  panel.grid = element_blank(),
                  axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
            labs(x = "migration distance (km)",
                 y = "climate exposure (sigma)")

      if(combine){
            p_map <- ggplot() +
                  geom_raster(data = ref,
                              aes(x, y, fill = z)) +
                  geom_path(data = pf,
                            aes(analog_x, analog_y, group = focal_index, color = I(color)),
                            linewidth = .1) +
                  geom_point(data = pf,
                             aes(analog_x, analog_y, group = focal_index, color = I(color)),
                             size = .5) +
                  geom_text(data = pf_focal,
                            aes(focal_x, focal_y, label = focal_index),
                            color = "white", nudge_x = .1, hjust = 0) +
                  scale_x_continuous(expand = c(0,0)) +
                  scale_y_continuous(expand = c(0,0)) +
                  scale_fill_gradientn(colors = c("black", "gray30")) +
                  theme_void() +
                  theme(legend.position = "none")

      }else{
            p_map <- ggplot() +
                  facet_wrap(~focal_index, scales = "free") +
                  geom_raster(data = ref,
                              aes(x, y, fill = z)) +
                  geom_path(data = pf,
                            aes(analog_x, analog_y, group = focal_index, color = I(color)),
                            linewidth = .1) +
                  geom_point(data = pf,
                             aes(analog_x, analog_y, group = focal_index, color = I(color)),
                             size = .5) +
                  scale_x_continuous(expand = c(0,0)) +
                  scale_y_continuous(expand = c(0,0)) +
                  scale_fill_gradientn(colors = c("black", "gray30")) +
                  theme_void() +
                  theme(legend.position = "none")

      }

      p_curves + p_dist + p_map +
            plot_layout(guides = "collect", widths = c(1, 3),
                        design = "AC
                        BC")
}

# background raster data
temp <- as.data.frame(clim2[[1]], xy = T) %>% rename(z = mean)

# pareto: temperature only
a1 <- find_analogs(focal = focal[, 1:3, drop = FALSE],
                   ref = clim2[[1]],
                   mode = "all")
p1 <- a1 %>%
      pareto_front(clim_tol = .01) %>%
      pareto_paths(x = a1, ref = temp)
ggsave("../pareto_fronts_tmean.png", p1,
       width = 10, height = 6, units = "in")

# pareto: temp and precip
a2 <- find_analogs(focal = focal,
                   ref = clim2,
                   mode = "all")
p2 <- a2 %>%
      pareto_front(clim_tol = .01) %>%
      pareto_paths(x = a2, ref = temp, combine = F)
ggsave("../pareto_fronts_tmean_prec.png", p2,
       width = 10, height = 6, units = "in")


stop("pareto computations complete")

#########

a %>%
      filter(focal_index == 16) %>%
      ggplot(aes(clim_dist, geog_dist)) +
      geom_point(color = "darkblue") +
      scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
      scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
      theme_bw() +
      labs(x = "climate novelty (sigma)",
           y = "geographic distance (km)")


ls_chars <- function(clim_dist, geog_dist, tol = .01, res = 50){
      exposure <- min(clim_dist[geog_dist < (tol * max(geog_dist, na.rm = T))])
      velocity <- min(geog_dist[clim_dist < (tol * max(clim_dist, na.rm = T))])

      # convexity
      clim_bin <- ceiling(clim_dist / max(clim_dist, na.rm = T) * res)
      geog_bin <- ceiling(geog_dist / max(geog_dist, na.rm = T) * res)
      g <- expand_grid(clim = unique(clim_bin),
                       geog = unique(geog_bin)) %>%
            left_join(bind_cols(clim = clim_bin, geog = geog_bin) %>% mutate(present = T),
                      by = join_by(clim, geog)) %>%
            mutate(present = ifelse(is.na(present), F, T))
      fit <- glm(present ~ clim + geog + clim*geog, data = g,
                 family = binomial("logit"))

      # g$fitted <- fitted(fit)
      # ggplot(g, aes(clim, geog, fill =fitted)) + geom_raster()
      curvature <- coef(fit)["clim:geog"]
      data.frame(curvature = curvature, exposure = exposure, velocity = velocity)
}



ls <- a %>%
      group_by(focal_index) %>%
      summarize(ls_chars(clim_dist, geog_dist))

hist(ls$curvature)
pairs(ls[,c("curvature", "velocity", "exposure")])
ggplot(ls, aes(exposure, curvature)) + geom_smooth() + geom_hline(yintercept = 0)

ls <- a %>%
      group_by(focal_index) %>%
      summarize(exposure = min(clim_dist[geog_dist < 3]),
                velocity = min(geog_dist[clim_dist < .05])) %>%
      ungroup() %>%
      mutate(exp = plyr::round_any(rank(exposure), 50, floor),
             vel = plyr::round_any(rank(velocity), 50, floor),
             vel = factor(vel, levels = rev(sort(unique(vel))))) %>%
      group_by(exp, vel) %>%
      filter(focal_index == focal_index[1])

left_join(ls, a) %>%
      ggplot(aes(clim_dist, geog_dist)) +
      facet_grid(vel ~ exp) +
      geom_point(color = "darkblue", size = .25) +
      scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
      scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
      theme_bw() +
      labs(x = "exposure (future climatic difference from focal site baseline)", y = "velocity (distance from focal site)")




ls %>%
      filter(vel == 50 & exp == 400 |
                   vel ==350 & exp == 400 |
                   vel == 350 & exp == 0 |
                   vel == 100 & exp == 0) %>%
      mutate(velocity = ifelse(vel == 350, "high", "low"),
             exp = ifelse(exp == 0, "low", "high"),
             exposure = factor(exp, levels = c("low", "high"))) %>%
      left_join(a) %>%
      ggplot(aes(clim_dist, geog_dist)) +
      facet_grid(velocity ~ exposure, labeller = label_both) +
      geom_point(color = "darkblue", size = 2) +
      annotate(geom = "point", x = 0, y = 0, color = "darkred", size = 100, alpha = .5) +
      scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
      scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
      theme_bw() +
      theme(strip.background = element_rect(fill = "black", color = "black"),
            strip.text = element_text(color = "white", size = 16)) +
      labs(x = "future climatic difference from focal site baseline", y = "distance from focal site")




