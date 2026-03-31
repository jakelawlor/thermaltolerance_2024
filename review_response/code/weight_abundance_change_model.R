# script to model abundance change as a function of thermal affinity

# HIGHLY SAMPLED ONLY

# libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggview)
library(grid)
theme_set(theme_bw())
#theme_update(panel.grid = element_blank())

# data --------------------------------------------------------------------
# upload full df for count species slopes (highly sampled only)
counts <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_counts_cattidalheight_HS.rds"
  )
) %>%
  select(-data, -model)

# upload full df for percent cover species slopes (highly sampled only)
cover <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_cover_cattidalheight_HS.rds"
  )
) %>%
  select(-data, -model)

# upload raw temperature values in the study area
temp <- read.csv(here::here(
  "data-processed",
  "appledore-island-env-data",
  "appledore_temps_1982_2023.csv"
))

gc()
gc()

# make models of temp to get the model-predicted min and max
tempmax_mod_minmax <- broom::augment(lm(data=temp, max ~ sample_year)) %>% filter(sample_year %in% c(1982, 2023)) %>% pull(.fitted)
tempmean_mod_minmax <- broom::augment(lm(data=temp, mean ~ sample_year)) %>% filter(sample_year %in% c(1982, 2023)) %>% pull(.fitted)
tempmin_mod_minmax <- broom::augment(lm(data=temp, min ~ sample_year)) %>% filter(sample_year %in% c(1982, 2023)) %>% pull(.fitted)



# find shared species
shared_spp <- unique(cover$organism)[unique(cover$organism) %in% unique(counts$organism)]


# make annotation coords --------------------------------------------------
ann_y <- .28

# find coordinate limits
xlim <- c(min(c(counts$mean_monthly_mean,cover$mean_monthly_mean),na.rm=T)-.85,
          max(c(counts$mean_monthly_mean,cover$mean_monthly_mean), na.rm=T) + .4)
ylim <- c(-max(c(counts$slope,cover$slope),na.rm=T) - .02,
          max(c(counts$slope,cover$slope), na.rm=T)+.02)


# plot counts -------------------------------------------------------------
# weight by inverse error
counts$weight <- 1/(counts$slope_se^2)

counts_mod <- lm(data = counts %>% select(slope, mean_monthly_mean, weight),
                 "slope ~ mean_monthly_mean",
                 weights = weight)
summary(counts_mod)
counts_p_val <- summary(counts_mod)$coefficients[2,"Pr(>|t|)"]
counts_p_val <- if(counts_p_val < 0.01){"< 0.01"} else
  if(counts_p_val < 0.05) {"< 0.05"} else
    paste0("= ",round(counts_p_val,2))
counts_mod_perf <- performance::check_model(counts_mod)
#plot(counts_mod)

counts_mod_au <- broom::augment(counts_mod,
                                counts,
                                interval = "confidence")

counts_mod_au %>% glimpse()

counts_p_base <- counts_mod_au %>%
  distinct(organism, slope, slope_se, weight, mean_monthly_mean, .upper, .lower, .fitted) %>%
  ggplot(aes(x = mean_monthly_mean,
             y = slope)) +
  annotate(geom = "rect",
           xmin = tempmean_mod_minmax[1],
           xmax  = tempmean_mod_minmax[2],
           # xmin = min(temp$mean),
           # xmax = max(temp$mean),
           ymin = -Inf,
           ymax = Inf,
           fill = "grey85",
           alpha = .5) +
  geom_vline(xintercept = tempmean_mod_minmax[1],
             linetype = "twodash",
             linewidth = .3,
             alpha = .7)+
  geom_vline(xintercept = tempmean_mod_minmax[2],
             linetype = "twodash",
             linewidth = .3,
             alpha = .7)+
  geom_hline(yintercept = 0,
             linetype = "longdash") +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "cyan3",
              alpha = .4) +
  geom_line(aes(y = .fitted)) +
  geom_segment(aes(y = slope + slope_se,
                   yend = slope - slope_se)) + 
  geom_point(shape = 21, 
             fill = "cyan3",
             aes(size = weight/min(weight))) +
  coord_cartesian(xlim = xlim,
                  ylim = ylim) +
  scale_x_continuous(   
    breaks = scales::pretty_breaks(),
    labels = ~paste0(.,"°"))  +
  labs(x = "Thermal Affinity\n(mean monthly average temperature\nof global occurrences)",
       y = "Coefficient of Abundance Change",
       title = "Density Species",
       size = "Weight\n(1/slope_se^2)\nrescaled to multiples") +
  annotate(geom = "text",
           x = xlim[1],
           y = ann_y,
           label = paste(sprintf('\u2191'),"Increasing\n   Abundance") ,
           hjust = 0,
           vjust = 1,
           lineheight = .8) +
  annotate(geom = "text",
           x = xlim[1],
           y = -ann_y,
           label = paste(sprintf('\u2193'),"Decreasing\n   Abundance") ,
           hjust = 0,
           vjust = 0,
           lineheight = .8) +
  geom_rug(data = temp,
           inherit.aes = F,
           aes(x = mean))

counts_p <- 
  counts_p_base  +
  annotate(geom = "text",
           x = xlim[2]-.1,
           y = ann_y,
           hjust = 1,
           vjust = 1,
           label = paste("n = ",nrow(counts %>%
                                       filter(!is.na(slope),
                                              !is.na(mean_monthly_mean)) %>%
                                       distinct(organism)))
  ) +
  annotate(geom = "text",
           x = xlim[2]-.1,
           y = -ann_y,
           hjust = 1,
           vjust = 0,
           lineheight = .8,
           label = paste0("coefficient: ",
                          round(counts_mod$coefficients["mean_monthly_mean"],3),
                          "\np ",counts_p_val)  ) +
  annotate(geom = "text",
           x = 10.88,
           y = -.205,
           hjust = .5,
           label = "Temperature\nChange",
           size = 2.3,
           lineheight = .8,
           vjust = 0) +
  annotate(geom = "segment",
           x = tempmean_mod_minmax[1] + .05,
           xend = tempmean_mod_minmax[2] - .05,
           y = -.215,
           yend = -.215,
           linewidth = .3,
           arrow = arrow(length = unit(4,"pt"))) +
  theme(panel.grid = element_blank(),
        plot.title.position = "plot") +
  scale_size_area(breaks = scales::pretty_breaks())


counts_p

# plot cover -------------------------------------------------------------
cover$weight <- 1/(cover$slope_se^2)

cover_mod <- lm(data = cover,
                "slope ~ mean_monthly_mean",
                weights = weight)
summary(cover_mod)
cover_p_val <- summary(cover_mod)$coefficients[2,"Pr(>|t|)"]
cover_p_val <- if(cover_p_val < 0.01){"< 0.01"} else
  if(cover_p_val < 0.05) {"< 0.05"} else
    paste0("= ",round(cover_p_val,2))
cover_mod_perf <- performance::check_model(cover_mod)

#plot(cover_mod)

cover_mod_au <- broom::augment(cover_mod,
                               cover %>% filter(!is.na(mean_monthly_mean)),
                               interval = "confidence")

cover_mod_au %>% glimpse()

cover_p_base <- cover_mod_au %>%
  distinct(organism, slope, mean_monthly_mean, slope_se,  .upper, .lower, .fitted, weight) %>%
  ggplot(aes(x = mean_monthly_mean,
             y = slope)) +
  annotate(geom = "rect",
           xmin = tempmean_mod_minmax[1],
           xmax  = tempmean_mod_minmax[2],
           # xmin = min(temp$mean),
           # xmax = max(temp$mean),
           ymin = -Inf,
           ymax = Inf,
           fill = "grey85",
           alpha = .5) +
  geom_vline(xintercept = tempmean_mod_minmax[1],
             linetype = "twodash",
             linewidth = .3,
             alpha = .7)+
  geom_vline(xintercept = tempmean_mod_minmax[2],
             linetype = "twodash",
             linewidth = .3,
             alpha = .7)+
  geom_hline(yintercept = 0,
             linetype = "longdash") +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "cyan3",
              alpha = .4) +
  geom_line(aes(y = .fitted)) +
  geom_segment(aes(y = slope + slope_se,
                   yend = slope - slope_se)) + 
  geom_point(shape = 21, 
             fill = "cyan3",
             aes(size = weight/min(weight))) +
  # geom_point(data = . %>% filter(organism %in% shared_spp),
  #            aes(fill = organism),
  #            shape = 21,
  #            size = 3) +
  coord_cartesian(xlim = xlim,
                  ylim = ylim) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(),
    labels = ~paste0(.,"°"))  +
  labs(x = "Thermal Affinity\n(mean monthly average temperature\nof global occurrences)",
       y = "Coefficient of Abundance Change",
       title = "Percent Cover Species",
       fill = "Shared Species",
       size = "Weight\n(1/slope_se^2)\nrescaled to multiples") +
  annotate(geom = "text",
           x = xlim[1],
           y = ann_y,
           label = paste(sprintf('\u2191'),"Increasing\n   Abundance") ,
           hjust = 0,
           vjust = 1,
           lineheight = .8) +
  annotate(geom = "text",
           x = xlim[1],
           y = -ann_y,
           label = paste(sprintf('\u2193'),"Decreasing\n   Abundance") ,
           hjust = 0,
           vjust = 0,
           lineheight = .8) +
  geom_rug(data = temp,
           inherit.aes = F,
           aes(x = mean))

cover_p <- 
  cover_p_base +
  annotate(geom = "text",
           x = xlim[2]-.1,
           y = ann_y,
           hjust = 1,
           vjust = 1,
           label = paste("n = ",nrow(cover %>%
                                       filter(!is.na(slope),
                                              !is.na(mean_monthly_mean)) %>%
                                       distinct(organism)))
  ) +
  annotate(geom = "text",
           x = xlim[2]-.1,
           y = -ann_y,
           hjust = 1,
           vjust = 0,
           lineheight = .8,
           label = paste0("coefficient: ",
                          round(cover_mod$coefficients["mean_monthly_mean"],3),
                          "\np ",cover_p_val)
  ) +
  annotate(geom = "text",
           x = 10.88,
           y = -.205,
           hjust = .5,
           label = "Temperature\nChange",
           size = 2.3,
           lineheight = .8,
           vjust = 0) +
  annotate(geom = "segment",
           x = tempmean_mod_minmax[1] + .05,
           xend = tempmean_mod_minmax[2] - .05,
           y = -.215,
           yend = -.215,
           linewidth = .3,
           arrow = arrow(length = unit(4,"pt"))) +
  theme(panel.grid = element_blank(), plot.title.position = "plot") +
  scale_size_area(breaks = scales::pretty_breaks())

cover_p
range(cover$weight)



# merge plots -------------------------------------------------------------

p_full <- counts_p + cover_p +
  plot_layout(axis_titles = "collect",
  ) +
  # plot_annotation(tag_levels = "a") &
  theme(plot.margin = margin(0,0,0,0)) +
  plot_annotation(theme = theme(plot.margin = margin(0,r=1,0,0)))

p_full +
  canvas(width  = 10,
         height = 4,
         unit = "in")


# make histograms of slope_se and weights for each plot
hist_se_counts <- counts %>%
  ggplot() +
  geom_density(aes(x = slope_se),
               fill = "grey") +
  coord_cartesian(xlim = c(0, .06),
                  expand = F) +
  theme_classic() +
  labs(x = "Standard Error of Abundance Model Slope",
       y = "Density")

cover %>%
  select(organism, 
         slope, 
         slope_se, 
         weight) %>% arrange(desc(weight))

hist_se_cover <- cover %>%
  ggplot() +
  geom_density(aes(x = slope_se),
               fill = "grey") +
  coord_cartesian(xlim = c(0, .06),
                  expand = F) +
  theme_classic() +
  labs(x = "Standard Error of Abundance Model Slope",
       y = "Density")

hist_se <- hist_se_counts + hist_se_cover +
  plot_layout(axis_titles = "collect") 

hist_weights_counts <- counts %>%
  ggplot() +
  geom_density(aes(x = weight/min(weight)),
               fill = "grey") +
  coord_cartesian(expand = F) +
  theme_classic() +
  labs(x = "Weighting Factor\nweight/min(weight)",
       y = "Density")

hist_weights_cover <- cover %>%
  ggplot() +
  geom_density(aes(x = weight/min(weight)),
               fill = "grey") +
  coord_cartesian(expand = F) +
  theme_classic() +
  labs(x = "Weighting Factor\nweight/min(weight)",
       y = "Density")

hist_weights <- hist_weights_counts + hist_weights_cover +
  plot_layout(axis_titles = "collect") 


left_col <- (counts_p + 
               theme(legend.position = "top") +
               guides(size = guide_legend(nrow = 2, byrow = TRUE))
             ) / hist_se_counts / hist_weights_counts +
  plot_layout(nrow = 3, heights = c(1,.2,.2)) +
  plot_annotation(tag_levels = list(c('a', 'b','c')))


right_col <- (cover_p + 
                theme(legend.position = "top") +
                guides(size = guide_legend(nrow = 2, byrow = TRUE))
              ) / hist_se_cover / hist_weights_cover +
  plot_layout(nrow = 3, heights = c(1,.2,.2)) +
  plot_annotation(tag_levels = list(c('d', 'e','f')))




p_full2 <- (wrap_elements(left_col) | wrap_elements(right_col))

p_full2 + ggview::canvas(8,8)


ggsave(p_full2,
       filename = here::here("review_response",
                             "figures",
                             "Fig3_abundance_change_slopes_by_temp_HS_weighted.png"),
       width = 8,
       height = 8,
       units = "in")

