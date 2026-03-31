# find richness change at the transect level over time
# to get the overall trend of richness increase on Appledore Island


library(ggplot2)
library(dplyr)



# load presence only data -------------------------------------------------
# first, upload presence only data
p_only <- readr::read_csv(
  here::here("data-processed",
             "spp_pres_by_replicate.csv")
) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  arrange(tidalheight) %>%
  mutate(transect_label = paste0("Transect ",transect)) %>%
  arrange(transect) %>%
  mutate(transect_label = forcats::fct_inorder(transect_label))



# standardize replicates --------------------------------------------------
# first, cut to only replicate = 1 to standardize effort across quadrays
# this loses a lot of data, but I think necessary
p_only_rep1 <- p_only %>%
  filter(replicate == 1)
rm(p_only)



# cut to highly sampled tidal heights -------------------------------------
# next, identify highly-sampled tidal heights
p_only_rep1 %>% distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point() +
  facet_wrap(~transect_label) +
  geom_hline(yintercept = -.25) +
  geom_hline(yintercept = 3)
# this leaves 9 tidalheights between these bounds (0-3). 

p_only_rep1 %>% distinct(tidalheight) %>% 
  count(tidalheight >= 0 & tidalheight < 3)

p_only_bound <- p_only_rep1 %>%
  filter(tidalheight >= 0 & tidalheight < 3)


# filter out undersampled years ------------------------------------------------------------
# within these tidal height bounds, there are 9 tidal heights sampled. 
# cut out years in which there are 6 or fewer of these 9 sampled. 
p_only_bound %>%
  distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point() +
  facet_wrap(~transect_label) +
  geom_hline(yintercept = -.25) +
  geom_hline(yintercept = 3)

transect_years_to_keep <- p_only_bound %>% 
  group_by(year, transect_label) %>% 
  distinct(tidalheight) %>% 
  count(name = "levels_sampled") %>%
  mutate(keep = ifelse(levels_sampled > 6, "yes","no")) %>%
  ungroup()

p_only_bound_yearskeep <- p_only_bound %>%
  left_join(transect_years_to_keep) %>%
  filter(keep == "yes")

# plot
p_only_rep1 %>% distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point(color = "black",
             shape = 21,
             fill = "transparent",
             stroke = .2) +
  facet_wrap(~transect_label) +
  geom_hline(yintercept = -.25) +
  geom_hline(yintercept = 3) +
  geom_point(data = p_only_bound_yearskeep,
             color = "black",
             fill = "cyan4",
             shape = 21,
             stroke = .1)

# filter out short timeseries ---------------------------------------------
# some transects are measured for only ~3 years in the beginning of the 
# time series. These aren't going to form a good trend in such little time. 
# filter to only transects in which sampling occurs over 
transects_to_keep <- 
  p_only_bound_yearskeep %>%
  group_by(transect_label) %>%
  distinct(year) %>%
  count() %>%
  mutate(keep_transect = ifelse(n < 20,"no","yes")) %>% ungroup()

transects_to_keep %>% count(keep_transect)

p_only_bound_yearskeep_transectskeep <- 
  p_only_bound_yearskeep %>%
  left_join(transects_to_keep) %>%
  filter(keep_transect == "yes") 


# plot again
p_only_rep1 %>% distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point(color = "grey50",
             size = 1) +
  facet_wrap(~transect_label, nrow = 3) +
  geom_hline(yintercept = -.25) +
  geom_hline(yintercept = 3) +
  geom_point(data = p_only_bound_yearskeep_transectskeep,
             color = "cyan4",
             size = 1)


# calculate richness ------------------------------------------------------
rich_over_time <- p_only_bound_yearskeep_transectskeep %>% 
  group_by(year, transect_label) %>%
  distinct(organism) %>% 
  count(name = "richness")

rich_over_time %>% 
  ggplot(aes(x = year, 
             y = richness,
             group = transect_label)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~transect_label)


# model -------------------------------------------------------------------
richmod <- lm(data = rich_over_time,
              "richness ~ year + transect_label")
summary(richmod)

au <- broom::augment(richmod,
                     rich_over_time,
                     interval = "prediction")

# use emmeans or marginal effects to get only effect of year
p_faceted <- au %>% 
  ggplot(aes(x = year,
             y = richness,
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
year_effect <- effect("year", richmod)

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
    title = "Total Richness per Transect",
    x = "Year",
    y = "Predicted Species Richness"
  ) 

total_rich_eff + ggview::canvas(4,4)

ggsave(total_rich_eff,
       filename = "outputs/richness_change/total_rich_change.png",
       width = 4,
       height = 4)


# combine with richness over depth plot
library(patchwork)
full_rich_plot <- (total_rich_eff + p2.2.2) 

full_rich_plot + ggview::canvas(width = 8, height = 4)

ggsave(full_rich_plot, 
       filename = "outputs/richness_change/full_rich_plot.png",
       width = 8,
       height = 4)

sjPlot::plot_model(richmod, type = "pred", terms = c("year"),
                   color = "cyan4",
                   alpha = .6) +
  geom_line(data = au,
            inherit.aes=F,
            aes(x = year, y  = .fitted, group = transect_label),
            alpha = .5) +
  geom_segment(data = au,
               inherit.aes = F,
               aes(x = year,
                   xend = year,
                   y = .fitted,
                   yend = richness),
               linetype = "longdash",
               linewidth = .3,
               alpha = .4) +
  geom_point(data = au,
             inherit.aes = F,
             aes(x= year, 
                 y = richness),
             size = .5,
             alpha = .4)


sjPlot::plot_model(richmod, type = "pred", terms = c("year","transect_label"),
                   color = "cyan4",
                   alpha = .6, show.legend = F) 
  geom_point(data = au,
             inherit.aes = F,
             aes(x= year, 
                 y = richness,
                 color = transect_label),
             size = 1,
             alpha = .8) +
  scale_color_viridis_d() +
  guides(fill = guide_none(),
         color = guide_legend())



# scratch -----------------------------------------------------------------

# try with an interaction term
richmod2 <- lm(data = rich_over_time,
               "richness ~ year * transect_label")
summary(richmod2)

car::Anova(richmod2)

au2 <- broom::augment(richmod2,
                      rich_over_time,
                      interval = "confidence")

p_faceted2 <- au2 %>% 
  ggplot(aes(x = year,
             y = richness,
             group = transect_label)) +
  geom_point() +
  # gghighlight::gghighlight(use_direct_label = F,
  #                          unhighlighted_params = list(alpha = .2)) +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              alpha = .4,
              fill = "cyan4") +
  geom_line(aes(y = .fitted)) +
  facet_wrap(~transect_label, nrow = 2) +
  labs(y = "Full-Transect Species Richness",
       x = NULL) +
  scale_x_continuous(breaks = c(1980, 2000, 2020)) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1))

p_faceted2 + ggview::canvas(4,4)

performance::compare_performance(richmod, richmod2, rank=T)
# now do it a different way -----------------------------------------------
# find full island species richness
set.seed(1)
# here, let's find the richness on the full island, rarefying to the 
# minimum number of transects sampled, or something like that







