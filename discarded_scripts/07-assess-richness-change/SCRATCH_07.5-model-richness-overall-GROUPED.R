# find island-wide richness change  over time - GROUPED
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
  mutate(tidalheight = (13-level)*.348)  %>%
  arrange(tidalheight) %>%
  mutate(transect_label = paste0("Transect ",transect)) %>%
  arrange(transect) %>%
  mutate(transect_label = forcats::fct_inorder(transect_label))

# upload traits to see how many of each species are in each group ---------
traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed.csv") %>%
  mutate(trait_group = case_when(motility_adult == "sessile" & group == "Invertebrate" ~ "Sessile Invertebrate",
                                 motility_adult == "motile" & group == "Invertebrate" ~ "Motile Invertebrate",
                                 group == "Algae" ~ "Algae"
  )) %>%
  rename(organism = gen_spp)

p_only <- p_only %>% left_join(traits)


# find richness in replicates
p_only_rich <- p_only %>%
  group_by(year, transect, transect_label, replicate, trait_group) %>%
  # find richness
  summarize(richness = n_distinct(organism),
            n_tidal_level = n_distinct(tidalheight, level),
            .groups = "drop")  %>%
  filter(n_tidal_level >= 6) %>%
  select(-n_tidal_level)
p_only_rich %>% glimpse()


# plot
p_only_rich %>%
  ggplot(aes(x = year, 
             y = richness, 
             color = trait_group)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~transect_label)

# rarify replicates -------------------------------------------------------
p_only_rich_rare <- p_only_rich %>%
  group_by(year, transect, transect_label, trait_group) %>%
  summarise(rich_rarefied = mean(richness),
            rich_min = min(richness),
            rich_max = max(richness),
            n_reps = n_distinct(replicate))

p_only_rich_rare %>%
  ggplot(aes(x = year, 
             y = rich_rarefied,
             color = trait_group)) +
  geom_segment(aes(y = rich_min, 
                   yend = rich_max,
                   xend = year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~transect_label)





# model -------------------------------------------------------------------
sess_df <- p_only_rich_rare %>% filter(trait_group == "Sessile Invertebrate")
sess_mod <- lm(data = sess_df,
              "rich_rarefied ~ year + transect_label")
summary(sess_mod)

mot_df <- p_only_rich_rare %>% filter(trait_group == "Motile Invertebrate")
mot_mod <- lm(data = mot_df,
               "rich_rarefied ~ year + transect_label")
summary(mot_mod)

alg_df <- p_only_rich_rare %>% filter(trait_group == "Algae")
alg_mod <- lm(data = alg_df,
               "rich_rarefied ~ year + transect_label")
summary(alg_mod)

# augment datasets
sess_au <- broom::augment(sess_mod,
                          sess_df,
                     interval = "prediction")

mot_au <- broom::augment(mot_mod,
                          mot_df,
                          interval = "prediction")

alg_au <- broom::augment(alg_mod,
                         alg_df,
                         interval = "prediction")

p_faceted_mot <- mot_au%>% 
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

p_faceted_mot + ggview::canvas(4,4)

p_faceted_sess <- sess_au %>% 
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

p_faceted_sess + ggview::canvas(4,4)

p_faceted_al <- alg_au%>% 
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

p_faceted_al + ggview::canvas(4,4)




# plot marginal effect ----------------------------------------------------
library(effects)

# Get the effect of year
year_effect_sess <- effect("year", sess_mod, 
                      xlevels = max(sess_df$year)-min(sess_df$year)+1)
year_effect_mot <- effect("year", mot_mod, 
                           xlevels = max(mot_df$year)-min(mot_df$year)+1)
year_effect_alg <- effect("year", alg_mod, 
                           xlevels = max(alg_df$year)-min(alg_df$year)+1)

# Convert to data frame for ggplot
year_effect_df_sess <- as.data.frame(year_effect_sess)
year_effect_df_mot <- as.data.frame(year_effect_mot)
year_effect_df_alg <- as.data.frame(year_effect_alg)

year_effect_df <- 
  year_effect_df_sess %>% mutate(trait_group = "Sessile Invertebrates") %>%
  rbind(year_effect_df_mot %>% mutate(trait_group = "Motile Invertebrates")) %>%
  rbind(year_effect_df_alg %>% mutate(trait_group = "Algae"))

# Plot
total_rich_eff <- 
  ggplot(year_effect_df, aes(x = year, y = fit,
                             #color = trait_group,
                             fill = trait_group)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.4) +
  geom_line(color = "black") +
  ggthemes::theme_few() +
  theme(legend.position = "bottom",
        #legend.position.inside = c(.98,.02),
        #legend.justification = c(1,0)
        ) +
  labs(
    title = "Overall Richness by Trait Group",
    x = "Year",
    y = "Species Richness per Transect",
    fill = "Organism Group",
    color = "Organism Group"
  ) 

total_rich_eff + ggview::canvas(4,4)

# calculate percentiles 
percentiles_sess <- sess_df %>%
  group_by(year) %>%
  summarize(meanrich = mean(rich_rarefied),
            rich_25 = quantile(rich_rarefied,.25),
            rich_75 = quantile(rich_rarefied,.75)) %>%
  mutate(trait_group = "Sessile Invertebrates")
percentiles_mot <- mot_df %>%
  group_by(year) %>%
  summarize(meanrich = mean(rich_rarefied),
            rich_25 = quantile(rich_rarefied,.25),
            rich_75 = quantile(rich_rarefied,.75)) %>%
  mutate(trait_group = "Motile Invertebrates")
percentiles_alg <- alg_df %>%
  group_by(year) %>%
  summarize(meanrich = mean(rich_rarefied),
            rich_25 = quantile(rich_rarefied,.25),
            rich_75 = quantile(rich_rarefied,.75)) %>%
  mutate(trait_group = "Algae")

group_percentiles <- 
  percentiles_sess %>%
  rbind(percentiles_mot) %>%
  rbind(percentiles_alg)
  

total_rich_eff_with_segments <- total_rich_eff +
  geom_segment(data = group_percentiles,
               aes(x = year, 
                   xend = year,
                   y = rich_25,
                   yend = rich_75,
                   color = trait_group),
               linewidth = .2,
               show.legend = F) +
  geom_point(data = group_percentiles,
             aes(x = year, 
                 y = meanrich),
             size = .7,
             shape = 21,
             show.legend = F)  +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top"))) +
  theme(legend.box.margin = margin(t = -7, 0,0,0,'mm'),
        legend.background = element_blank())  +
  scale_fill_manual(breaks = c("Algae","Motile Invertebrates","Sessile Invertebrates"),
                    labels = c("Algae","Motile\nInvertebrates","Sessile\nInvertebrates"),
                    values = RColorBrewer::brewer.pal(3,"Dark2")) +
  scale_color_manual(breaks = c("Algae","Motile Invertebrates","Sessile Invertebrates"),
                    labels = c("Algae","Motile\nInvertebrates","Sessile\nInvertebrates"),
                    values =  RColorBrewer::brewer.pal(3,"Dark2"))


total_rich_eff_with_segments + ggview::canvas(4,4)

#ggsave(total_rich_eff_with_segments,
#       filename = "outputs/richness_change/overall_rich_plot_by_group.png",
#       width = 4,
#       height = 4)
#saveRDS(total_rich_eff_with_segments, 
#        "outputs/richness_change/overall_rich_plot_by_group.rds")
list <- list(year_effect = year_effect_df,
             percentiles = group_percentiles)
#saveRDS(list,
#        "outputs/richness_change/overall_rich_plot_by_group_data.rds")

total_rich_eff_with_segments_annotate + ggview::canvas(4,4)



