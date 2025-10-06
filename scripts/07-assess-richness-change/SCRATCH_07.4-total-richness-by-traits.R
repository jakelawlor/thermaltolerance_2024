# find richness change at the transect level over time
# to get the overall trend of richness increase on Appledore Island


library(ggplot2)
library(dplyr)
traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed.csv")



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
  mutate(transect_label = forcats::fct_inorder(transect_label)) %>%
  left_join(traits %>%
              rename(organism = gen_spp) %>%
              mutate(group = stringr::str_to_title(motility_adult)))



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
  group_by(year, transect_label, group) %>%
  distinct(organism) %>% 
  count(name = "richness")

rich_over_time %>% 
  ggplot(aes(x = year, 
             y = richness,
             color = group)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~transect_label)


# model -------------------------------------------------------------------
richmod <- lm(data = rich_over_time,
              "richness ~ year*group + transect_label")
summary(richmod)

au <- broom::augment(richmod,
                     rich_over_time,
                     interval = "prediction")

# use emmeans or marginal effects to get only effect of year
p_faceted <- au %>% 
  ggplot(aes(x = year,
             y = richness,
             color = group)) +
  geom_point() +
  # gghighlight::gghighlight(use_direct_label = F,
  #                          unhighlighted_params = list(alpha = .2)) +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  fill = group),
              alpha = .3,
              color = "transparent"
              ) +
  geom_line(aes(y = .fitted)) +
  facet_wrap(~transect_label, nrow = 2) +
  labs(y = "Full-Transect Species Richness",
       x = NULL,
       color = "Motility Group",
       fill = "Motility Group") +
  scale_x_continuous(breaks = c(1980, 2000, 2020)) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1,
                                   size = 7),
        panel.spacing.x = unit(.6,"lines"),
        strip.text = element_text(size = 8)) +
  theme(legend.position = "inside",
        legend.position.inside = c(.98,.02),
        legend.justification = c(1,0)) +
  scale_fill_manual(values = c("cyan2",
                               "darkcyan")) +
  scale_color_manual(values = c("cyan2",
                                "darkcyan"))

p_faceted + ggview::canvas(6,4)




# plot marginal effect ----------------------------------------------------
library(effects)

# Get the effect of year
year_effect <- margialavg_slopes("year:groupSessile", richmod)

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

