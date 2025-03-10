# track changes in species richness across tidal heights - BY GROUP



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(broom)
library(ggtext)
theme_set(ggthemes::theme_few())


# data --------------------------------------------------------------------
# first, upload presence only data
p_only <- readr::read_csv(
  here::here("data-processed",
             "appledore-survey-data",
             "pa-with-therm",
             "pa-with-therm-all.csv")
) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  arrange(tidalheight) 

traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed.csv") %>%
  rename(organism = gen_spp) %>%
  mutate(org_group = case_when(motility_adult == "sessile" & group == "Invertebrate" ~ "Sessile Invertebrate",
                   motility_adult == "motile" & group == "Invertebrate" ~ "Motile Invertebrate",
                   group == "Algae" ~ "Algae"
  ))

p_only <- p_only %>% left_join(traits)

p_only %>% distinct(organism, group, org_group) %>% count(org_group)

# find richness with presence only ----------------------------------------
# here, we omit quadrats that say they were sampled,
# but have a richness of zero. Instead, we'll use just quadrats that have
# some richness value. (below we'll include the zeros too)
richness_no_zeros <- p_only %>%
  group_by(year, transect, level, tidalheight, replicate, org_group) %>%
  # find richness
  summarize(richness = sum(pres),
            .groups = "drop") 

# rarify quadrats by replicate
# becuase some quadrats have up to 4 replicates, but many just have 1,
# we need to rarefy quadrats that have replicates to the average richness
# value in one replicate so we can fairly compare by quadrats with and without replicates
richness_no_zeros_rarified <- 
  richness_no_zeros %>%
  # some transect/tidalheights have multiple replicates
  # here, average those into one, because not all have replicates
  group_by(year, transect, tidalheight, level, org_group) %>%
  summarize(replicates = n(),
            rarified_richness = mean(richness),
            .groups = "drop")
rm(p_only)


# summarize by year and level ---------------------------------------------
# because both the number and selection of transects is so irregular,
# we're going to summarize to get a mean value of richness per tidal level
# around the island each year. these mean level richnesses will be 
# calculated from very different numbers of samples, but this will make it 
# more consistent. 
richness_no_zeros_by_level <- 
  richness_no_zeros_rarified %>%
  group_by(year,level,tidalheight, org_group) %>%
  summarize(mean_richness = mean(rarified_richness),
            n_transects_sampled = n(),
            .groups = "drop") %>%
  filter(level >= 0,
         level <= 15)

# quick plot to view 
#richness_no_zeros_by_level %>%# filter(n_samples >= 3) %>%
#  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) %>%
#  ggplot(aes(x = year,
#             y = mean_richness,
#             color = org_group)) +
#  geom_point() +
#  facet_wrap(~tidalheight_f) +
#  geom_smooth(method = "lm")
rm(richness_no_zeros_rarified)



# view richness across tidal heights
#richness_no_zeros_by_level %>%
#  ggplot(aes(x = tidalheight, 
#             y = mean_richness,
#             color = org_group)) +
#  geom_point() +
#  geom_smooth(method = "loess", se = TRUE)
#
#
#richness_no_zeros_by_level %>%
#  group_by(tidalheight, org_group) %>% 
#  summarize(meanrich = mean(mean_richness)) %>%
#  arrange(desc(meanrich))
## in both cases, tidal height 0.696 is the height with the highest 
## mean species richness. I think we should recenter around that value,
## at least for the quadratic model, to accurately capture the peak


# build models ------------------------------------------------------------
#richness_no_zeros_by_level %>%
#  ggplot(aes(x = mean_richness,
#             color = org_group)) +
#  geom_histogram(fill = "transparent",
#                 position = "identity") +
#  facet_wrap(~level)
#richness_no_zeros_by_level %>% 
#  ggplot(aes(x = mean_richness,
#             color = org_group)) +
#  geom_histogram(fill = "transparent",
#                 position = "identity")



# make function to test multiple models  -------------------------------------------
test_mods <- function(data){
  
  models <- list(
    mod1 = glm(mean_richness ~ year, data = data, family = Gamma(link = "log")),
    mod2 = glm(mean_richness ~ year + tidalheight, data = data, family = Gamma(link = "log")),
    mod3 = glm(mean_richness ~ year * tidalheight, data = data, family = Gamma(link = "log")),
    mod4 = glm(mean_richness ~ year * I(tidalheight^2), data = data, family = Gamma(link = "log")),
    mod5 = glm(mean_richness ~ year + I(tidalheight^2), data = data, family = Gamma(link = "log")),
    mod6 = glm(mean_richness ~ I(tidalheight^2), data = data, family = Gamma(link = "log")),
    mod7 = glm(mean_richness ~ tidalheight, data = data, family = Gamma(link = "log"))
  )
  
  perf <- performance::compare_performance(
    models,
    rank =T)
  print(perf)
  
  best_mod <- perf$Name[1]
  print(paste("best =", best_mod))
  out <- models[[best_mod]]
  return(out)
}


# test on spp groups ------------------------------------------------------
sess_inv <- richness_no_zeros_by_level %>% filter(org_group == "Sessile Invertebrate")
sess_mod <- test_mods(sess_inv)
# with full data, mod6: mean_richness ~ I(tidalheight^2)
# with HS data, mod6: mean_richness ~ I(tidalheight^2)
mot_inv <- richness_no_zeros_by_level %>% filter(org_group == "Motile Invertebrate")
mot_mod <- test_mods(mot_inv)
# with full data, mod3: mean_richness ~ year * tidalheight
# with HS data, mod5: year + I(tidalheight^2)
alg <- richness_no_zeros_by_level %>% filter(org_group == "Algae")
alg_mod <- test_mods(alg)
# with full data, mod4: mean_richness ~ year * I(tidalheight^2)
# with HS data, mod4: mean_richness ~ year * I(tidalheight^2)


predict_df_sess_inv <- data.frame(
  tidalheight = rep(unique(sess_inv$tidalheight),
                    times = length(unique(sess_inv$year))),
  year = rep(unique(sess_inv$year),
             each = length(unique(sess_inv$tidalheight)))
)
predict_df_mot_inv <- data.frame(
  tidalheight = rep(unique(mot_inv$tidalheight),
                    times = length(unique(mot_inv$year))),
  year = rep(unique(mot_inv$year),
             each = length(unique(mot_inv$tidalheight)))
)
predict_df_alg <- data.frame(
  tidalheight = rep(unique(alg$tidalheight),
                    times = length(unique(alg$year))),
  year = rep(unique(alg$year),
             each = length(unique(alg$tidalheight)))
)


sess_inv_au <- broom::augment(sess_mod,
                              newdata = predict_df_sess_inv,
                              type.predict = "response",
                              se_fit = T,
                              simulate_pi  = T) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 
mot_inv_au <- broom::augment(mot_mod,
                              newdata = predict_df_mot_inv,
                             type.predict = "response",
                             se_fit = T,
                              simulate_pi  = T) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 
alg_au <- broom::augment(alg_mod,
                              newdata = predict_df_alg,
                         type.predict = "response",
                         se_fit = T,
                              simulate_pi  = T) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 



# plot --------------------------------------------------------------------
sess_plot <- sess_inv_au %>%
  left_join(sess_inv) %>%
  ggplot(aes(x =.fitted, 
             y = tidalheight, 
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .9) +
  geom_path(data = . %>%  ungroup() %>% 
              filter(year == min(year))%>%
              arrange(tidalheight_f),
            aes(group = NULL),
            color = "black") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  {if(stringr::str_detect(deparse(sess_mod$formula),"year")){
    geom_path(arrow = arrow(length = unit(6,"pt"),
                            type = "open"),
              # position = position_nudge(y = .07),
              color = "black",
              linewidth = 1,
              linejoin = "mitre") 
  }} +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Sessile Invertebrates") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none") +
  annotate(x = 5.85, y = 4.85,
           geom = "richtext",
           hjust = 1, 
           vjust = 1,
           label.color = "transparent",
           label = glue::glue(
             deparse(sess_mod$formula) %>% 
               stringr::str_replace("mean_richness","Richness") %>%
               stringr::str_replace("I\\(tidalheight\\^2\\)","Tidal Level<sup>2</sup>") %>%
               stringr::str_to_title()
                             ),
           ) +
  ggview::canvas(4,4)

mot_plot <- mot_inv_au %>%
  left_join(mot_inv) %>%
  ggplot(aes(x =.fitted, 
             y = tidalheight, 
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .9) +
  geom_path(data = . %>%  ungroup() %>% 
              filter(year == min(year))%>%
              arrange(tidalheight_f),
            aes(group = NULL),
            color = "black") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  
  {if(stringr::str_detect(deparse(mot_mod$formula),"year")){
    geom_path(arrow = arrow(length = unit(6,"pt"),
                            type = "open"),
              # position = position_nudge(y = .07),
              color = "black",
              linewidth = 1,
              linejoin = "mitre") 
  }} +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Motile Invertebrates") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  annotate(x = 5.85, y = 4.85,
           geom = "richtext",
           hjust = 1, 
           vjust = 1,
           label.color = "transparent",
           label = glue::glue(
             deparse(mot_mod$formula) %>% 
               stringr::str_replace("mean_richness","Richness") %>%
               stringr::str_replace("I\\(tidalheight\\^2\\)","Tidal Level<sup>2</sup>") %>%
               stringr::str_to_title()
           ),
  ) +
  ggview::canvas(4,4)


alg_plot <- alg_au %>%
  left_join(alg) %>%
  ggplot(aes(x =.fitted, 
             y = tidalheight, 
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .9) +
  geom_path(data = . %>%  ungroup() %>% 
              filter(year == min(year))%>%
              arrange(tidalheight_f),
            aes(group = NULL),
            color = "black") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  
  {if(stringr::str_detect(deparse(alg_mod$formula),"year")){
    geom_path(arrow = arrow(length = unit(6,"pt"),
                            type = "open"),
              # position = position_nudge(y = .07),
              color = "black",
              linewidth = 1,
              linejoin = "mitre") 
  }} +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Algae") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank()) +
  annotate(x = 5.85, y = 3.6,
           geom = "richtext",
           hjust = 1, 
           vjust = 1,
           label.color = "transparent",
           label = glue::glue(
             deparse(alg_mod$formula) %>% 
               stringr::str_replace("mean_richness","Richness") %>%
               stringr::str_replace("I\\(tidalheight\\^2\\)","Tidal Level<sup>2</sup>") %>%
               stringr::str_to_title()
           ),
  ) +  ggview::canvas(4,4)


library(patchwork)
full_p <- alg_plot | sess_plot | mot_plot
full_p + ggview::canvas(12,4)
ggsave(full_p,
       filename = "outputs/richness_change/rich_per_level_by_group.png",
       width = 12, 
       height = 4)
saveRDS(full_p,
        "outputs/richness_change/rich_per_level_by_group.rds")

sess_inv_au %>% 
  left_join(sess_inv) %>%  
  mutate(tidalheight_f2 = paste0(round(as.numeric(as.character(tidalheight_f)),2),"m")) %>% 
  mutate(tidalheight_f2 = forcats::fct_reorder(tidalheight_f2, as.numeric(tidalheight_f))) %>%
  ggplot(aes(x = year,
             y = .fitted,
             color = tidalheight_f2,
             group = tidalheight_f2,
             fill = tidalheight_f2)) +
  geom_point(aes(y = mean_richness),
             shape = 21,
             alpha = .8, 
             color = "black") +
  gghighlight::gghighlight(
    unhighlighted_params = list(fill = "grey75",
                                alpha = .1),
    use_direct_label = F) +
  geom_ribbon(aes(ymin = .fitted - 2*.se.fit,
                  ymax = .fitted + 2*.se.fit),
              alpha = .75) +
  geom_line(linewidth = .5,
            color = "black") +
  
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(.~tidalheight_f2) +
  theme_bw() + 
  theme(legend.position = "none") +

  labs(x = NULL,
       y = "Rarified Species Richness",
       title = "Richness Changes Over Time"
  ) +
  geom_text(
    data = . %>%
      group_by(tidalheight, tidalheight_f2) %>%
      summarize(change = round(last(.fitted) - first(.fitted),2)),
    aes(x = 1982,
        y = 5,
        label = change),
    hjust = 0, 
    vjust = 1,
    size =3.5,
    color = "black",
    fontface = "bold",
    stat = "unique"
  ) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust  = 1),
        strip.background = element_rect(color = "black",
                                        fill = "white"),
        strip.text = element_text(size = 10,
                                  margin = margin( b = 0, t = 2)))




# repeat all with highly sampled data only --------------------------------


# data --------------------------------------------------------------------
# first, upload presence only data
p_only_hs <- readr::read_csv(
  here::here("data-processed",
             "appledore-survey-data",
             "pa-with-therm",
             "pa-with-therm-highly-sampled.csv")
) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  arrange(tidalheight) 

p_only_hs <- p_only_hs %>% left_join(traits)

p_only_hs %>% distinct(organism, group, org_group) %>% count(org_group)

# find richness with presence only ----------------------------------------
# here, we omit quadrats that say they were sampled,
# but have a richness of zero. Instead, we'll use just quadrats that have
# some richness value. (below we'll include the zeros too)
richness_no_zeros_hs <- p_only_hs %>%
  group_by(year, transect, level, tidalheight, replicate, org_group) %>%
  # find richness
  summarize(richness = sum(pres),
            .groups = "drop") 

# rarify quadrats by replicate
# becuase some quadrats have up to 4 replicates, but many just have 1,
# we need to rarefy quadrats that have replicates to the average richness
# value in one replicate so we can fairly compare by quadrats with and without replicates
richness_no_zeros_rarified_hs <- 
  richness_no_zeros_hs %>%
  # some transect/tidalheights have multiple replicates
  # here, average those into one, because not all have replicates
  group_by(year, transect, tidalheight, level, org_group) %>%
  summarize(replicates = n(),
            rarified_richness = mean(richness),
            .groups = "drop")
rm(p_only_hs)


# summarize by year and level ---------------------------------------------
# because both the number and selection of transects is so irregular,
# we're going to summarize to get a mean value of richness per tidal level
# around the island each year. these mean level richnesses will be 
# calculated from very different numbers of samples, but this will make it 
# more consistent. 
richness_no_zeros_by_level_hs <- 
  richness_no_zeros_rarified_hs %>%
  group_by(year,level,tidalheight, org_group) %>%
  summarize(mean_richness = mean(rarified_richness),
            n_transects_sampled = n(),
            .groups = "drop") %>%
  filter(level >= 0,
         level <= 15)

# quick plot to view 
#richness_no_zeros_by_level %>%# filter(n_samples >= 3) %>%
#  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) %>%
#  ggplot(aes(x = year,
#             y = mean_richness,
#             color = org_group)) +
#  geom_point() +
#  facet_wrap(~tidalheight_f) +
#  geom_smooth(method = "lm")
rm(richness_no_zeros_rarified_hs)


# test on spp groups ------------------------------------------------------
sess_inv_hs <- richness_no_zeros_by_level_hs %>% filter(org_group == "Sessile Invertebrate")
sess_mod_hs <- test_mods(sess_inv_hs)
# with HS data, mod6: mean_richness ~ I(tidalheight^2)
mot_inv_hs <- richness_no_zeros_by_level_hs %>% filter(org_group == "Motile Invertebrate")
mot_mod_hs <- test_mods(mot_inv_hs)
# with HS data, mod5: year + I(tidalheight^2)
alg_hs <- richness_no_zeros_by_level_hs %>% filter(org_group == "Algae")
alg_mod_hs <- test_mods(alg_hs)
# with HS data, mod4: mean_richness ~ year * I(tidalheight^2)


predict_df_sess_inv_hs <- data.frame(
  tidalheight = rep(unique(sess_inv_hs$tidalheight),
                    times = length(unique(sess_inv_hs$year))),
  year = rep(unique(sess_inv_hs$year),
             each = length(unique(sess_inv_hs$tidalheight)))
)
predict_df_mot_inv_hs <- data.frame(
  tidalheight = rep(unique(mot_inv_hs$tidalheight),
                    times = length(unique(mot_inv_hs$year))),
  year = rep(unique(mot_inv_hs$year),
             each = length(unique(mot_inv_hs$tidalheight)))
)
predict_df_alg_hs <- data.frame(
  tidalheight = rep(unique(alg_hs$tidalheight),
                    times = length(unique(alg_hs$year))),
  year = rep(unique(alg_hs$year),
             each = length(unique(alg_hs$tidalheight)))
)


sess_inv_au_hs <- broom::augment(sess_mod_hs,
                              newdata = predict_df_sess_inv_hs,
                              type.predict = "response",
                              se_fit = T,
                              simulate_pi  = T) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 
mot_inv_au_hs <- broom::augment(mot_mod_hs,
                             newdata = predict_df_mot_inv_hs,
                             type.predict = "response",
                             se_fit = T,
                             simulate_pi  = T) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 
alg_au_hs <- broom::augment(alg_mod_hs,
                         newdata = predict_df_alg_hs,
                         type.predict = "response",
                         se_fit = T,
                         simulate_pi  = T) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 



# plot --------------------------------------------------------------------
sess_inv_au_hs %>%
  left_join(sess_inv_hs) %>%
  ggplot(aes(x =.fitted, 
             y = tidalheight, 
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .9) +
  geom_path(data = . %>%  ungroup() %>% 
              filter(year == min(year))%>%
              arrange(tidalheight_f),
            aes(group = NULL),
            color = "black") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Sessile Invertebrates") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank()) +
  ggview::canvas(4,4)

mot_inv_au_hs %>%
  left_join(mot_inv_hs) %>%
  ggplot(aes(x =.fitted, 
             y = tidalheight, 
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .9) +
  geom_path(data = . %>%  ungroup() %>% 
              filter(year == min(year))%>%
              arrange(tidalheight_f),
            aes(group = NULL),
            color = "black") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Motile Invertebrates") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank()) +
  ggview::canvas(4,4)


alg_au_hs %>%
  left_join(alg_hs) %>%
  ggplot(aes(x =.fitted, 
             y = tidalheight, 
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .9) +
  geom_path(data = . %>%  ungroup() %>% 
              filter(year == min(year))%>%
              arrange(tidalheight_f),
            aes(group = NULL),
            color = "black") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre")  +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Algae") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank()) +
  ggview::canvas(4,4)



sess_inv_au_hs %>% 
  left_join(sess_inv_hs) %>%  
  mutate(tidalheight_f2 = paste0(round(as.numeric(as.character(tidalheight_f)),2),"m")) %>% 
  mutate(tidalheight_f2 = forcats::fct_reorder(tidalheight_f2, as.numeric(tidalheight_f))) %>%
  ggplot(aes(x = year,
             y = .fitted,
             color = tidalheight_f2,
             group = tidalheight_f2,
             fill = tidalheight_f2)) +
  geom_point(aes(y = mean_richness),
             shape = 21,
             alpha = .8, 
             color = "black") +
  gghighlight::gghighlight(
    unhighlighted_params = list(fill = "grey75",
                                alpha = .1),
    use_direct_label = F) +
  geom_ribbon(aes(ymin = .fitted - 2*.se.fit,
                  ymax = .fitted + 2*.se.fit),
              alpha = .75) +
  geom_line(linewidth = .5,
            color = "black") +
  
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(.~tidalheight_f2) +
  theme_bw() + 
  theme(legend.position = "none") +
  
  labs(x = NULL,
       y = "Rarified Species Richness",
       title = "Richness Changes Over Time"
  ) +
  geom_text(
    data = . %>%
      group_by(tidalheight, tidalheight_f2) %>%
      summarize(change = round(last(.fitted) - first(.fitted),2)),
    aes(x = 1982,
        y = 5,
        label = change),
    hjust = 0, 
    vjust = 1,
    size =3.5,
    color = "black",
    fontface = "bold",
    stat = "unique"
  ) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust  = 1),
        strip.background = element_rect(color = "black",
                                        fill = "white"),
        strip.text = element_text(size = 10,
                                  margin = margin( b = 0, t = 2)))

