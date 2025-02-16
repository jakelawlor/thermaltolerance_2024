# abundance change by Topt, by trait group

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggview)
library(grid)
theme_set(theme_bw())
#theme_update(panel.grid = element_blank())
traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed.csv")

# data --------------------------------------------------------------------
# upload full df for count species slopes
counts <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_counts_cattidalheight.rds"
  )
) %>%
  select(-data, -model) %>%
  left_join(traits %>% rename(organism = gen_spp) %>%
              mutate(group = stringr::str_to_title(motility_adult)))

# upload the counts coefficients from tweedie models
counts_tweedie <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_counts_tweedie_cattidalheight.rds"
  )
) %>%
  left_join(counts %>% select(organism, mean_monthly_mean)) %>%
  select(-data, -model) %>%
  left_join(traits %>% rename(organism = gen_spp) %>%
              mutate(group = stringr::str_to_title(motility_adult)))


# upload full df for percent cover species slopes
cover <- readRDS(
  here::here(
    "data-processed",
    "abundance-data",
    "abundance_change_slopes_cover_cattidalheight.rds"
  )
) %>%
  select(-data, -model) %>%
  left_join(traits %>% rename(organism = gen_spp) %>%
              mutate(group = stringr::str_to_title(motility_adult)))

# upload raw temperature values in the study area
temp <- read.csv(here::here(
  "data-processed",
  "appledore_temps_1982_2023.csv"
))

temp_quant <- quantile(temp$mean, c(.1,.9))
# find how many species have thermal affinities within appledore mean temp range
cover %>% count(mean_monthly_mean > temp_quant[1] & 
                  mean_monthly_mean < temp_quant[2]) 
# 28 of 48 species for cover, or 58%
cover %>% count(mean_monthly_mean < temp_quant[1])
cover %>% count(mean_monthly_mean > temp_quant[2])
counts %>% count(mean_monthly_mean > temp_quant[1] & 
                   mean_monthly_mean < temp_quant[2]) 
# 14 of 26 for counts, or 54%


# find shared species
shared_spp <- unique(cover$organism)[unique(cover$organism) %in% unique(counts$organism)]


# make annotation coords --------------------------------------------------
ann_y <- .26

# plot counts -------------------------------------------------------------
counts_mod <- lm(data = counts,
                 "slope ~ mean_monthly_mean * group")
summary(counts_mod)
#plot(counts_mod)

counts_mod_au <- broom::augment(counts_mod,
                                counts,
                                interval = "prediction")

counts_mod_au %>% glimpse()

counts_p_base <- counts_mod_au %>%
  distinct(organism, slope, mean_monthly_mean, .upper, .lower, .fitted, group) %>%
  ggplot(aes(x = mean_monthly_mean,
             y = slope,
             color = group,
             fill = group)) +
  annotate(geom = "rect",
           xmin = temp_quant[1],
           xmax  = temp_quant[2],
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
             # fill = "cyan3",
              alpha = .4,
             color = "transparent") +
  geom_line(aes(y = .fitted),
            color = "black") +
  geom_point(shape = 21, 
             size = 3,
            # fill = "cyan3"
            ) +
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
       title = "Density Species") +
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
  geom_rug(data = temp,
           inherit.aes = F,
           aes(x = mean)) +
  scale_color_manual(values = c("cyan2","darkcyan")) +
  scale_fill_manual(values = c("cyan","darkcyan"))

counts_p_base

# plot cover -------------------------------------------------------------
cover_mod <- lm(data = cover,
                "slope ~ mean_monthly_mean")
summary(cover_mod)
plot(cover_mod)

cover_mod_au <- broom::augment(cover_mod,
                               cover %>% filter(!is.na(mean_monthly_mean)),
                               interval = "prediction")

cover_mod_au %>% glimpse()

cover_p_base <- cover_mod_au %>%
  distinct(organism, slope, mean_monthly_mean, .upper, .lower, .fitted) %>%
  ggplot(aes(x = mean_monthly_mean,
             y = slope)) +
  annotate(geom = "rect",
           xmin = temp_quant[1],
           xmax  = temp_quant[2],
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
  coord_cartesian(xlim =c(min(c(counts$mean_monthly_mean,cover$mean_monthly_mean),na.rm=T),
                          max(c(counts$mean_monthly_mean,cover$mean_monthly_mean), na.rm=T)),
                  ylim = c(-max(c(counts$slope,cover$slope),na.rm=T),
                           max(c(counts$slope,cover$slope), na.rm=T))) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°"))  +
  labs(x = "Thermal Affinity\n(mean monthly average temperature of global occurrences)",
       y = "Coefficient of Abundance Change",
       title = "Percent Cover Species",
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
  theme(panel.grid = element_blank())

cover_p


# merge plots -------------------------------------------------------------

p_full <- counts_p + cover_p +
  plot_layout(axis_titles = "collect",
              guides = "collect") &
  theme(legend.position = "bottom",
        plot.margin = margin(2,2,2,2))

p_full +
  canvas(width  = 8,
         height = 4,
         unit = "in")

