# script to model abundance change as a function of thermal affinity



# libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggview)
library(grid)
theme_set(theme_bw())
#theme_update(panel.grid = element_blank())

# data --------------------------------------------------------------------
# upload full df for count species slopes
counts <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_counts_cattidalheight.rds"
  )
) %>%
  select(-data, -model)

# upload the counts coefficients from tweedie models
counts_tweedie <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_counts_tweedie_cattidalheight.rds"
  )
) %>%
  left_join(counts %>% select(organism, mean_monthly_mean)) %>%
  select(-data, -model)

# upload full df for percent cover species slopes
cover <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_cover_cattidalheight.rds"
  )
) %>%
  select(-data, -model)

# upload raw temperature values in the study area
temp <- read.csv(here::here(
  "data-processed",
  "appledore-island-env-data",
  "appledore_temps_1982_2023.csv"
))

# make models of temp to get the model-predicted min and max
tempmax_mod_minmax <- broom::augment(lm(data=temp, max ~ sample_year)) %>% filter(sample_year %in% c(1982, 2023)) %>% pull(.fitted)
tempmean_mod_minmax <- broom::augment(lm(data=temp, mean ~ sample_year)) %>% filter(sample_year %in% c(1982, 2023)) %>% pull(.fitted)
tempmin_mod_minmax <- broom::augment(lm(data=temp, min ~ sample_year)) %>% filter(sample_year %in% c(1982, 2023)) %>% pull(.fitted)

# find how many species have thermal affinities within appledore mean temp range
cover %>% count(mean_monthly_mean > tempmean_mod_minmax[1] & 
                  mean_monthly_mean < tempmean_mod_minmax[2]) 
# 22 of 48 species for cover, or 46%
cover %>% count(mean_monthly_mean < tempmean_mod_minmax[1])
cover %>% count(mean_monthly_mean > tempmean_mod_minmax[2])
counts %>% count(mean_monthly_mean > tempmean_mod_minmax[1] & 
                  mean_monthly_mean < tempmean_mod_minmax[2]) 
# 13 of 26 for counts, or 50%


# find shared species
shared_spp <- unique(cover$organism)[unique(cover$organism) %in% unique(counts$organism)]


# make annotation coords --------------------------------------------------
ann_y <- .26

# plot counts -------------------------------------------------------------
counts_mod <- lm(data = counts %>% select(slope, mean_monthly_mean),
                 "slope ~ mean_monthly_mean")
summary(counts_mod)
counts_mod_perf <- performance::check_model(counts_mod)
#plot(counts_mod)

counts_mod_au <- broom::augment(counts_mod,
                                counts,
                                interval = "prediction")

counts_mod_au %>% glimpse()

counts_p_base <- counts_mod_au %>%
  distinct(organism, slope, mean_monthly_mean, .upper, .lower, .fitted) %>%
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
 # annotate(geom = "rect",
 #          xmin = quantile(temp$mean,.25),
 #          xmax = quantile(temp$mean,.75),
 #          ymin = -Inf,
 #          ymax = Inf,
 #          fill = "grey80",
 #          alpha = .5) +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "cyan3",
              alpha = .4) +
  geom_line(aes(y = .fitted)) +
  geom_point(shape = 21, 
             size = 3,
             fill = "cyan3") +
 # geom_point(data = . %>% filter(organism %in% shared_spp),
 #            aes(fill = organism),
 #            shape = 21,
 #            size = 3) +
  coord_cartesian(xlim = c(min(c(counts$mean_monthly_mean,cover$mean_monthly_mean),na.rm=T),
                           max(c(counts$mean_monthly_mean,cover$mean_monthly_mean), na.rm=T)),
                  ylim = c(-max(c(counts$slope,cover$slope),na.rm=T),
                           max(c(counts$slope,cover$slope), na.rm=T))) +
  scale_x_continuous(labels = paste0(c(8,10,12,14),"°"))  +
  labs(x = "Thermal Affinity\n(mean monthly average temperature of global occurrences)",
       y = "Coefficient of Abundance Change",
       title = "a) Density Species") +
  annotate(geom = "text",
           x = 7.7,
           y = ann_y,
           label = paste(sprintf('\u2191'),"Increasing\n   Abundance") ,
           hjust = 0,
           vjust = 1,
           lineheight = .8) +
  annotate(geom = "text",
           x = 7.7,
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
           x = 14,
           y = -ann_y,
           hjust = 1,
           vjust = 0,
           label = paste("n = ",nrow(counts %>%
                                       filter(!is.na(slope),
                                              !is.na(mean_monthly_mean)) %>%
                                       distinct(organism)))
           ) +
  annotate(geom = "text",
           x = 14,
           y = ann_y,
           hjust = 1,
           vjust = 1,
           lineheight = .8,
           label = paste0("coefficient: ",
                          round(counts_mod$coefficients["mean_monthly_mean"],3),
                          "\np < 0.01")  ) +
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
        plot.title.position = "plot")


counts_p

# plot cover -------------------------------------------------------------
cover_mod <- lm(data = cover,
                 "slope ~ mean_monthly_mean")
summary(cover_mod)
cover_mod_perf <- performance::check_model(cover_mod)

#plot(cover_mod)

cover_mod_au <- broom::augment(cover_mod,
                                cover %>% filter(!is.na(mean_monthly_mean)),
                                interval = "prediction")

cover_mod_au %>% glimpse()

cover_p_base <- cover_mod_au %>%
  distinct(organism, slope, mean_monthly_mean, .upper, .lower, .fitted) %>%
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
 # annotate(geom = "rect",
 #          xmin = quantile(temp$mean,.25),
 #          xmax = quantile(temp$mean,.75),
 #          ymin = -Inf,
 #          ymax = Inf,
 #          fill = "grey80",
 #          alpha = .5) +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "cyan3",
              alpha = .4) +
  geom_line(aes(y = .fitted)) +
  geom_point(shape = 21, 
             size = 3,
             fill = "cyan3") +
 # geom_point(data = . %>% filter(organism %in% shared_spp),
 #            aes(fill = organism),
 #            shape = 21,
 #            size = 3) +
  coord_cartesian(xlim =c(min(c(counts$mean_monthly_mean,cover$mean_monthly_mean),na.rm=T),
                           max(c(counts$mean_monthly_mean,cover$mean_monthly_mean), na.rm=T)),
                  ylim = c(-max(c(counts$slope,cover$slope),na.rm=T),
                           max(c(counts$slope,cover$slope), na.rm=T))) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°"))  +
  labs(x = "Thermal Affinity\n(mean monthly average temperature of global occurrences)",
       y = "Coefficient of Abundance Change",
       title = "b) Percent Cover Species",
       fill = "Shared Species") +
  annotate(geom = "text",
           x = 7.7,
           y = ann_y,
           label = paste(sprintf('\u2191'),"Increasing\n   Abundance") ,
           hjust = 0,
           vjust = 1,
           lineheight = .8) +
  annotate(geom = "text",
           x = 7.7,
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
           x = 14,
           y = -ann_y,
           hjust = 1,
           vjust = 0,
           label = paste("n = ",nrow(cover %>%
                                       filter(!is.na(slope),
                                              !is.na(mean_monthly_mean)) %>%
                                       distinct(organism)))
  ) +
  annotate(geom = "text",
           x = 14,
           y = ann_y,
           hjust = 1,
           vjust = 1,
           lineheight = .8,
           label = paste0("coefficient: ",
                          round(cover_mod$coefficients["mean_monthly_mean"],3),
                          "\np < 0.05")
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
  theme(panel.grid = element_blank(), plot.title.position = "plot")

cover_p


# merge plots -------------------------------------------------------------

p_full <- counts_p + cover_p +
  plot_layout(axis_titles = "collect",
              guides = "collect",
              ) +
 # plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        plot.margin = margin(0,0,0,0)) +
  plot_annotation(theme = theme(plot.margin = margin(0,r=1,0,0)))

p_full +
  canvas(width  = 8,
         height = 4,
         unit = "in")


ggsave(p_full,
       filename = here::here("outputs",
                             "abundance-change",
                  "Fig3_abundance_change_slopes_by_temp.png"),
       width = 8,
       height = 4,
       units = "in")


library(performance)
library(patchwork)

# Run model checks
check_cover <- plot(check_model(cover_mod, panel=F, base_size = 8,
                                size_title = 10))
check_count <- plot(check_model(counts_mod, panel = F, base_size =8,
                                size_title = 10))




cover_wrap <- wrap_plots(check_cover, nrow = 2) + plot_annotation(title = "Percent Cover Group")
count_wrap <- wrap_plots(check_count, nrow = 2) + plot_annotation(title = "Density Group")

combined <- cowplot::plot_grid(count_wrap, cover_wrap, nrow = 2) 
combined +  ggview::canvas(10,10)
ggsave(combined,
       filename = "outputs/abundance-change/model-assumptions-plot.png",
       width = 10, 
       height = 10)



# highlight shared species ------------------------------------------------
# find how many shared species have same sign coefficient
cover_mod_shared <- lm(data = cover %>% filter(organism %in% shared_spp),
                       "slope ~ mean_monthly_mean")
summary(cover_mod_shared)

counts_mod_shared <- lm(data = counts %>% filter(organism %in% shared_spp),
                        "slope ~ mean_monthly_mean")
summary(counts_mod_shared)


counts_p_shared <- counts_p_base +
  geom_point(data = . %>% filter(organism %in% shared_spp),
             aes(fill = organism),
             size = 3,
             shape = 21)  +
  scale_fill_viridis_d(option = "plasma") +
  
  geom_smooth(data = . %>% filter(organism %in% shared_spp),
              method = "lm",
              linetype = "dashed",
              color = "black",
              linewidth = .5,
              se = F) +
  theme(plot.title.position = "plot", 
        panel.grid = element_blank())

cover_p_shared <- cover_p_base +
  geom_point(data = . %>% filter(organism %in% shared_spp),
             aes(fill = organism),
             size = 3,
             shape = 21)  +
  scale_fill_viridis_d(option = "plasma") +
  geom_smooth(data = . %>% filter(organism %in% shared_spp),
              method = "lm",
              linetype = "dashed",
              color = "black",
              linewidth = .5,
              se=F) +
  theme(plot.title.position = "plot",
        panel.grid = element_blank())

p_shared_full <- (counts_p_shared | cover_p_shared)  +
  plot_layout(guides = "collect",
              axis_titles = "collect") & 
  theme(legend.position = 'right')

# get legend width so we know how much to add to side of plot
legend_grob <- ggplotGrob(p_shared_full)$grobs[which(sapply(ggplotGrob(p_shared_full)$grobs, function(x) x$name) == "guide-box")][[1]]

# Calculate legend dimensions
legend_width <- convertWidth(grobWidth(legend_grob), "in", valueOnly = TRUE)


p_shared_full +
  canvas(width  = 8 + legend_width,
         height = 4,
         unit = "in")

ggsave(p_shared_full,
       filename = here::here("outputs",
                             "abundance-change",
                             "Fig3_abundance_change_slopes_by_temp_SHARED_SPP.png"),
       width = 8 + legend_width,
       height = 4,
       units = "in")






# repeat with tweedie -----------------------------------------------------
counts_mod_tweedie <- lm(data = counts_tweedie,
                 "slope_tweedie ~ mean_monthly_mean")
summary(counts_mod_tweedie)

counts_mod_tweedie_au <- broom::augment(counts_mod_tweedie,
                                counts_tweedie,
                                interval = "prediction")

counts_mod_tweedie_au %>% glimpse()

counts_p_tweedie <- counts_mod_tweedie_au %>%
  distinct(organism, slope_tweedie, mean_monthly_mean, .upper, .lower, .fitted) %>%
  ggplot(aes(x = mean_monthly_mean,
             y = slope_tweedie)) +
  annotate(geom = "rect",
           xmin = tempmean_mod_minmax[1],
           xmax  = tempmean_mod_minmax[2],
          # xmin = min(temp$mean),
          # xmax = max(temp$mean),
           ymin = -Inf,
           ymax = Inf,
           fill = "grey85",
           alpha = .5) +
  # annotate(geom = "rect",
  #          xmin = quantile(temp$mean,.25),
  #          xmax = quantile(temp$mean,.75),
  #          ymin = -Inf,
  #          ymax = Inf,
  #          fill = "grey80",
  #          alpha = .5) +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "cyan3",
              alpha = .4) +
  geom_line(aes(y = .fitted)) +
  geom_point(shape = 21, 
             size = 3,
             fill = "cyan3") +
  # geom_point(data = . %>% filter(organism %in% shared_spp),
  #            aes(fill = organism),
  #            shape = 21,
  #            size = 3) +
  coord_cartesian(xlim = c(min(c(counts_tweedie$mean_monthly_mean,cover$mean_monthly_mean),na.rm=T),
                           max(c(counts_tweedie$mean_monthly_mean,cover$mean_monthly_mean), na.rm=T)),
                  ylim = c(-max(c(counts$slope,cover$slope),na.rm=T),
                           max(c(counts$slope,cover$slope), na.rm=T))) +
 # scale_x_continuous(labels = paste0(c(8,10,12,14),"°"))  +
  labs(x = "Thermal Affinity\n(mean monthly average temperature of global occurrences)",
       y = "Coefficient of Abundance Change",
       title = "a) Density Species",
       fill = "Shared Species") +
  annotate(geom = "text",
           x = 7.7,
           y = ann_y,
           label = paste("Increasing Abundance",sprintf('\u2191')) ,
           hjust = 0,
           vjust = 1) +
  annotate(geom = "text",
           x = 7.7,
           y = -ann_y,
           label = paste("Decreasing Abundance",sprintf('\u2193')) ,
           hjust = 0,
           vjust = 0) +
  annotate(geom = "text",
           x = 14,
           y = -ann_y,
           hjust = 1,
           vjust = 0,
           label = paste("n = ",nrow(counts_tweedie %>%
                                       filter(!is.na(slope_tweedie),
                                              !is.na(mean_monthly_mean)) %>%
                                       distinct(organism)))
  ) +
  annotate(geom = "text",
           x = 14,
           y = ann_y,
           hjust = 1,
           vjust = 1,
           lineheight = .8,
           label = paste0("coefficient: ",
                          round(counts_mod_tweedie$coefficients["mean_monthly_mean"],3),
                          "\np < 0.01") 
           ) +
  # scale_fill_viridis_d(option = "plasma") +
  # theme(legend.position = "bottom") +
  geom_rug(data = temp,
           inherit.aes = F,
           aes(x = mean)) +
  theme(panel.grid = element_blank())

counts_p_tweedie

counts_p_compare <- (counts_p + labs(title = "a) Density Species - Gamma") | 
    counts_p_tweedie +labs(title = "b) Density Species - Tweedie"))+
  plot_layout(axis_titles = "collect",
              guides = "collect") &
  theme(legend.position = "bottom",
        plot.margin = margin(2,2,2,2))


counts_p_compare +
  canvas(width  = 8,
         height = 4,
         unit = "in")

ggsave(counts_p_compare,
       filename = here::here("outputs",
                             "abundance-change",
                             "counts_gamma_tweedie_compare.png"),
       width = 8,
       height = 4,
       units = "in")



# make compare plot in one ------------------------------------------------

counts_p_compare_2 <-
counts_mod_tweedie_au %>%
  select(organism,
         slope_tweedie,
         mean_monthly_mean,
         .fitted,
         .upper,
         .lower) %>%
  rename(slope = slope_tweedie) %>%
  mutate(group = "Tweedie") %>%
  rbind(counts_mod_au %>%
              select(organism,
                     slope,
                     mean_monthly_mean,
                     .fitted,
                     .upper,
                     .lower)%>%
          mutate(group = "Gamma"))  %>%
  mutate(group = forcats::fct_inorder(group)) %>%
  
  # start plot
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
          alpha = .5,
          #stat = "unique"
          ) +
  geom_ribbon( alpha = .2,
               show.legend = F,
               fill = "grey50",
               aes(group = group,
                   color = group,
                   ymin = .lower,
                   ymax = .upper)) +
  geom_line(aes(y = .fitted,
                color = group),
            show.legend = F) +
  
  geom_point(shape = 21,
             size = 3,
             color = "black",
             aes(fill = group)) +
  
 # gghighlight::gghighlight(use_direct_label = F) +
  
  scale_fill_manual(values = c("darkmagenta",
                               "cyan4"))+
  scale_color_manual(values = c("darkmagenta",
                               "cyan4")) +
  
  
  geom_rug(data = temp,
           inherit.aes = F,
           aes(x = mean)) +
  geom_hline(mapping = NULL,
             data = NULL,
             yintercept = 0) +
  labs(x = "Thermal Affinity\n(mean monthly average temperature of global occurrences)",
       y = "Coefficient of Abundance Change",
       title = "Density Species",
       fill = "Abundance Model") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        panel.grid = element_blank())

counts_p_compare_2

counts_p_compare_2 +
  canvas(width= 4.5,
         height = 4,
         unit = "in")


ggsave(counts_p_compare_2,
       filename = here::here("outputs",
                             "abundance-change",
                             "counts_gamma_tweedie_compare_2.png"),
       width = 4.5,
       height = 4,
       units = "in")



# -------------------------------------------------------------------------
# repeat with other thermal affinity metrics ------------------------------------------------
# -------------------------------------------------------------------------

tempmin_mod_minmax
# first for min monthly temps
min_plot_base <- ggplot() +
  annotate(geom = "rect",
           xmin = tempmin_mod_minmax[1],#quantile(temp$min,.1)[[1]],
           xmax  = tempmin_mod_minmax[2],#quantile(temp$min,.9)[[1]],
           ymin = -Inf,
           ymax = Inf,
           fill = "grey85",
           alpha = .5)  +
  geom_rug(data = temp,
           aes(x = min),
           linewidth = .2, 
           length = unit(0.04,"npc")) +
  scale_x_continuous(labels = ~paste0(.,"°")) +
  labs(x = "Thermal Minima",
       y = "Coefficient of Change") +
  coord_cartesian(xlim = c(
    min(counts$perc_monthly_min_05, cover$perc_monthly_min_05, na.rm=T),
    max(counts$perc_monthly_min_05, cover$perc_monthly_min_05, na.rm=T)
   ),
   ylim = c(
     max(counts$slope, cover$slope) *-1,
     max(counts$slope, cover$slope)
   )) +
  theme(panel.grid = element_blank())

# model for count species
counts_mod_min <- lm(data = counts,
                     "slope ~ perc_monthly_min_05")
summary(counts_mod_min)
counts_min_au <- broom::augment(counts_mod_min,
                                counts, 
                                interval = "confidence")
# make plot
plot_min_counts <- min_plot_base + 
  geom_ribbon(data = counts_min_au,
              aes(ymin = .lower,
                  ymax = .upper,
                  x = perc_monthly_min_05),
              fill = "cyan3",
              alpha = .5) +
  geom_line(data = counts_min_au,
            aes(x = perc_monthly_min_05,
                y = .fitted)) +
  geom_point(data = counts_min_au,
              aes(x = perc_monthly_min_05,
                  y = slope),
             shape = 21, 
             size = 1.75,
             stroke = .25,
             fill = "cyan3") +
    geom_hline(yintercept = 0,
             linetype = "longdash") +
 # labs(title = "Density Group") +
  annotate(geom = "text",
           x = 7.5,
           y = -.25,
           hjust = 1,
           vjust = 0,
           #fontface = "bold",
           label = "model p = 0.27",
           size = 3) 
plot_min_counts

# repeat for cover species
cover_mod_min <- lm(data = cover %>% filter(!is.na(perc_monthly_min_05)),
                     "slope ~ perc_monthly_min_05")
summary(cover_mod_min)
cover_min_au <- broom::augment(cover_mod_min,
                                cover %>% filter(!is.na(perc_monthly_min_05)), 
                                interval = "confidence")

# make plot
plot_min_cover <- min_plot_base + 
  geom_ribbon(data = cover_min_au,
              aes(ymin = .lower,
                  ymax = .upper,
                  x = perc_monthly_min_05),
              fill = "cyan3",
              alpha = .5) +
  geom_line(data = cover_min_au,
            aes(x = perc_monthly_min_05,
                y = .fitted)) +
  geom_point(data = cover_min_au,
             aes(x = perc_monthly_min_05,
                 y = slope),
             shape = 21, 
             size = 1.75,
             stroke = .25,
             fill = "cyan3") +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
 # labs(title = "Percent Cover Group") +
  annotate(geom = "text",
           x = 7.5,
           y = -.25,
           hjust = 1,
           vjust = 0,
           fontface = "bold",
           label = "model p < 0.01**",
        size = 3) 
plot_min_cover

( plot_min_counts | plot_min_cover ) + ggview::canvas(8,4)


# repeat with mean ---------------------------------------------------------
mean_plot_base <- ggplot() +
  annotate(geom = "rect",
           xmin = tempmean_mod_minmax[1],#quantile(temp$mean,.1)[[1]],
           xmax  = tempmean_mod_minmax[2],#quantile(temp$mean,.9)[[1]],
           ymin = -Inf,
           ymax = Inf,
           fill = "grey85",
           alpha = .5)  +
  geom_rug(data = temp,
           aes(x = mean),
           linewidth = .2, 
           length = unit(0.04,"npc")) +
  scale_x_continuous(labels = ~paste0(.,"°")) +
  labs(x = "Mean Thermal Affinity",
       y = "Coefficient of Change") +
  coord_cartesian(xlim = c(
    min(counts$mean_monthly_mean, cover$mean_monthly_mean, na.rm=T),
    max(counts$mean_monthly_mean, cover$mean_monthly_mean, na.rm=T)
  ),
  ylim = c(
    max(counts$slope, cover$slope) *-1,
    max(counts$slope, cover$slope)
  )) +
  theme(panel.grid = element_blank())
mean_plot_base

# model for count species
counts_mod_mean <- lm(data = counts,
                     "slope ~ mean_monthly_mean")
summary(counts_mod_mean)
counts_mean_au <- broom::augment(counts_mod_mean,
                                counts, 
                                interval = "confidence")
# make plot
plot_mean_counts <- mean_plot_base + 
  geom_ribbon(data = counts_mean_au,
              aes(ymin = .lower,
                  ymax = .upper,
                  x = mean_monthly_mean),
              fill = "cyan3",
              alpha = .5) +
  geom_line(data = counts_mean_au,
            aes(x = mean_monthly_mean,
                y = .fitted)) +
  geom_point(data = counts_mean_au,
             aes(x = mean_monthly_mean,
                 y = slope),
             shape = 21, 
             size = 1.75,
             stroke = .25,
             fill = "cyan3") +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
 # labs(title = "Count Species") +
  annotate(geom = "text",
           x = 14.1,
           y = -.23,
           hjust = 1,
           vjust = 0,
           fontface = "bold",
           label = "model p < 0.01**",
           size = 3)
plot_mean_counts

# repeat for cover species
cover_mod_mean <- lm(data = cover %>% filter(!is.na(mean_monthly_mean)),
                    "slope ~ mean_monthly_mean")
summary(cover_mod_mean)
cover_mean_au <- broom::augment(cover_mod_mean,
                               cover %>% filter(!is.na(mean_monthly_mean)), 
                               interval = "confidence")

# make plot
plot_mean_cover <- mean_plot_base + 
  geom_ribbon(data = cover_mean_au,
              aes(ymin = .lower,
                  ymax = .upper,
                  x = mean_monthly_mean),
              fill = "cyan3",
              alpha = .5) +
  geom_line(data = cover_mean_au,
            aes(x = mean_monthly_mean,
                y = .fitted)) +
  geom_point(data = cover_mean_au,
             aes(x = mean_monthly_mean,
                 y = slope),
             shape = 21, 
             size = 1.75,
             stroke = .25,
             fill = "cyan3") +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
 # labs(title = "Cover Species") +
  annotate(geom = "text",
           x = 14.1,
           y = -.23,
           hjust = 1,
           vjust = 0,
           fontface = "bold",
           label = "model p < 0.05*",
           size = 3)
plot_mean_cover

( plot_mean_counts | plot_mean_cover ) + ggview::canvas(8,4)



# repeat with max ---------------------------------------------------------
max_plot_base <- ggplot() +
  annotate(geom = "rect",
           xmin = tempmax_mod_minmax[1],#quantile(temp$max,.1)[[1]],
           xmax = tempmax_mod_minmax[2],#quantile(temp$max,.9)[[1]],
           ymin = -Inf,
           ymax = Inf,
           fill = "grey85",
           alpha = .5)  +
  geom_rug(data = temp,
           aes(x = max),
           linewidth = .2, 
           length = unit(0.04,"npc")) +
  scale_x_continuous(labels = ~paste0(.,"°")) +
  labs(x = "Thermal Maxima",
       y = "Coefficient of Change") +
  coord_cartesian(xlim = c(
    min(counts$perc_monthly_max_95, cover$perc_monthly_max_95, na.rm=T),
    max(counts$perc_monthly_max_95, cover$perc_monthly_max_95, na.rm=T)
  ),
  ylim = c(
    max(counts$slope, cover$slope) *-1,
    max(counts$slope, cover$slope)
  )) +
  theme(panel.grid = element_blank()) 
max_plot_base

# model for count species
counts_mod_max <- lm(data = counts,
                     "slope ~ perc_monthly_max_95")
summary(counts_mod_max)
counts_max_au <- broom::augment(counts_mod_max,
                                counts, 
                                interval = "confidence")
# make plot
plot_max_counts <- max_plot_base + 
  geom_ribbon(data = counts_max_au,
              aes(ymin = .lower,
                  ymax = .upper,
                  x = perc_monthly_max_95),
              fill = "cyan3",
              alpha = .5) +
  geom_line(data = counts_max_au,
            aes(x = perc_monthly_max_95,
                y = .fitted)) +
  geom_point(data = counts_max_au,
             aes(x = perc_monthly_max_95,
                 y = slope),
             shape = 21, 
             size = 1.75,
             stroke = .25,
             fill = "cyan3") +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
 # labs(title = "Count Species") +
  annotate(geom = "text",
           x = 29,
           y = -.23,
           hjust = 1,
           vjust = 0,
           fontface = "bold",
           label = "model p < 0.01**",
           size = 3) +
  labs(title = "Density Group")
plot_max_counts

# repeat for cover species
cover_mod_max <- lm(data = cover %>% filter(!is.na(mean_monthly_max)),
                    "slope ~ perc_monthly_max_95")
summary(cover_mod_max)
cover_max_au <- broom::augment(cover_mod_max,
                               cover %>% filter(!is.na(mean_monthly_max)), 
                               interval = "confidence")

# make plot
plot_max_cover <- max_plot_base + 
  geom_ribbon(data = cover_max_au,
              aes(ymin = .lower,
                  ymax = .upper,
                  x = perc_monthly_max_95),
              fill = "cyan3",
              alpha = .5) +
  geom_line(data = cover_max_au,
            aes(x = perc_monthly_max_95,
                y = .fitted)) +
  geom_point(data = cover_max_au,
             aes(x = perc_monthly_max_95,
                 y = slope),
             shape = 21, 
             size = 1.75,
             stroke = .25,
             fill = "cyan3") +
  geom_hline(yintercept = 0,
             linetype = "longdash") +
#  labs(title = "Cover Species") +
  annotate(geom = "text",
           x = 29,
           y = -.23,
           hjust = 1,
           vjust = 0,
           fontface = "bold",
           label = "model p < 0.01**",
           size = 3) +
  labs(title = 'Percent Cover Group')
plot_max_cover

( plot_max_counts | plot_max_cover ) + ggview::canvas(8,4)



# arrange all -------------------------------------------------------------
theme_set(theme_bw(base_size = 10) +
            theme(axis.title.x = element_text(margin = margin(b=-4, t =1,0,0,"pt")),
                  #plot.background = element_blank()
                  ))
full_p <- (((plot_max_counts + theme(panel.border = element_rect(linewidth = 1.7)) |
               plot_max_cover + theme(panel.border = element_rect(linewidth = 1.7))) + 
              plot_layout(axis_title = "collect",
                          axes = "collect")) / 
             ((plot_mean_counts + theme(panel.border = element_rect(linewidth = 1.7)) |
                 plot_mean_cover + theme(panel.border = element_rect(linewidth = 1.7)))+ 
                plot_layout(axis_title = "collect",
                            axes = "collect")) / 
             ((plot_min_counts  |
                 plot_min_cover+ theme(panel.border = element_rect(linewidth = 1.7)))+
                plot_layout(axis_title = "collect",
                            axes = "collect"))) #+
#plot_annotation(tag_levels = "a") #&
# theme(plot.margin = margin(0,0,b=3,0,unit = "pt"))

full_p + ggview::canvas(4,6)
ggsave("outputs/abundance-change/abundance_by_t_all_metrics.png",
       width = 4,
       height = 6,
       unit = "in")

rm(list = ls())

