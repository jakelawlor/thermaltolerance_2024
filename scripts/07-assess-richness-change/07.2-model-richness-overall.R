# find island-wide richness change  over time
# to get the overall trend of richness increase on Appledore Island



# libraries ---------------------------------------------------------------
library(ggplot2)
library(dplyr)



# load presence only data -------------------------------------------------
# first, upload presence only data
p_only <- readr::read_csv(
  here::here("data-processed",
             "appledore-survey-data",
             "pa-with-therm",
             "pa-with-therm-highly-sampled.csv")
) %>%
  mutate(tidalheight = (13.5-level)*.3048)  %>%
  arrange(tidalheight) %>%
  mutate(transect_label = paste0("Transect ",transect)) %>%
  arrange(transect) %>%
  mutate(transect_label = forcats::fct_inorder(transect_label))


# find richness in replicates
p_only_rich <- p_only %>%
  group_by(year, transect, transect_label, replicate) %>%
  # find richness in each entire transect/replicate/year
  summarize(richness = n_distinct(organism),
            n_tidal_level = n_distinct(tidalheight, level),
            .groups = "drop")  %>%
  # filter out replicates with 6 or fewer tidal heights sampled
  # note that all transects should have 
  filter(n_tidal_level >= 6) %>%
  select(-n_tidal_level)
p_only_rich %>% glimpse()

# plot
p_only_rich %>%
  ggplot(aes(x = year, 
             y = richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~transect_label)

# rarify replicates -------------------------------------------------------
p_only_rich_rare <- p_only_rich %>%
  group_by(year, transect, transect_label) %>%
  summarise(rich_rarefied = mean(richness),
            rich_min = min(richness),
            rich_max = max(richness),
            n_reps = n_distinct(replicate))

p_only_rich_rare %>%
  ggplot(aes(x = year, 
             y = rich_rarefied)) +
  geom_segment(aes(y = rich_min, 
                   yend = rich_max,
                   xend = year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~transect_label)


# model -------------------------------------------------------------------
richmod <- lm(data = p_only_rich_rare,
              "rich_rarefied ~ year + transect_label")
summary(richmod)

library(performance)
p <- plot(performance::check_model(richmod, panel = T))
p <- p & labs(subtitle= NULL)
p + ggview::canvas(12,8)
ggsave(p, 
       filename = "outputs/richness_change/total_rich_assumptions.png",
       width = 12,
       height = 8)

au <- broom::augment(richmod,
                     p_only_rich_rare,
                     interval = "prediction")

# use emmeans or marginal effects to get only effect of year
p_faceted <- au %>% 
  ggplot(aes(x = year,
             y = rich_rarefied,
             group = transect_label)) +
  geom_point() +
  # gghighlight::gghighlight(use_direct_label = F,
  #                          unhighlighted_params = list(alpha = .2)) +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              alpha = .4,
              fill = "cyan3") +
  geom_line(aes(y = .fitted)) +
  facet_wrap(~transect_label, nrow = 2) +
  labs(y = "Full-Transect Species Richness",
       x = NULL) +
  scale_x_continuous(breaks = c(1980, 2000, 2020)) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1,
                                   size = 7),
        panel.spacing.x = unit(.6,"lines"),
        strip.text = element_text(size = 8))

p_faceted + ggview::canvas(4,4)




# plot marginal effect ----------------------------------------------------
library(effects)

# Get the effect of year
year_effect <- effect("year", richmod, 
                      xlevels = max(p_only_rich_rare$year)-min(p_only_rich_rare$year)+1)

# Convert to data frame for ggplot
year_effect_df <- as.data.frame(year_effect)

# Plot
total_rich_eff <- 
  ggplot(year_effect_df, aes(x = year, y = fit)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.4, 
              fill = "cyan3") +
  geom_line(color = "black") +
  ggthemes::theme_few() +
  labs(
    title = "a) Overall Richness Change",
    x = "Year",
    y = "Species Richness per Transect"
  ) 

total_rich_eff + ggview::canvas(4,4)

# calculate percentiles 
percentiles <- p_only_rich_rare %>%
  group_by(year) %>%
  summarize(meanrich = mean(rich_rarefied),
            rich_25 = quantile(rich_rarefied,.25),
            rich_75 = quantile(rich_rarefied,.75))

total_rich_eff_with_segments <- total_rich_eff +
  geom_segment(data = percentiles,
               aes(x = year, 
                   xend = year,
                   y = rich_25,
                   yend = rich_75),
               linewidth = .2) +
  geom_point(data = percentiles,
               aes(x = year, 
                   y = meanrich),
             size = .7) +
  coord_cartesian(xlim = c(1980, 2024),
                  ylim = c(7, 25.5),
                  expand = F)


total_rich_eff_with_segments + ggview::canvas(4,4)


total_rich_eff_with_segments_annotate <- total_rich_eff_with_segments +
  annotate(geom = "text",
           x = 1981.25,
           y = 24.75, 
           label = "Year effect: 0.05\np < 0.01",
           hjust = 0,
           vjust = 1,
           lineheight = .85)

total_rich_eff_with_segments_annotate + ggview::canvas(4,4)


ggsave(total_rich_eff_with_segments_annotate,
       filename = "outputs/richness_change/total_rich_change.png",
       width = 4,
       height = 4)

saveRDS(total_rich_eff_with_segments_annotate,
        file = "outputs/richness_change/overall_rich_plot.rds")

total_rich_eff_with_segments2 <- total_rich_eff +
  geom_point(data = p_only_rich_rare,
             aes(x = year, 
                 y = rich_rarefied),
             size = .7) 

