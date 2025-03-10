# track changes in species richness across tidal heights


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(broom)
theme_set(theme_bw())


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


# find richness with presence only ----------------------------------------
# here, we omit quadrats that say they were sampled,
# but have a richness of zero. Instead, we'll use just quadrats that have
# some richness value. (below we'll include the zeros too)
richness_no_zeros <- p_only %>%
  group_by(year, transect, level, tidalheight, replicate) %>%
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
  group_by(year, transect, tidalheight, level) %>%
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
  group_by(year,level,tidalheight) %>%
  summarize(mean_richness = mean(rarified_richness),
            n_transects_sampled = n(),
            .groups = "drop") %>%
  filter(level >= 0,
         level <= 15)

# quick plot to view 
richness_no_zeros_by_level %>%# filter(n_samples >= 3) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) %>%
  ggplot(aes(x = year,
             y = mean_richness)) +
  geom_point() +
  facet_wrap(~tidalheight_f) +
  geom_smooth(method = "lm")
rm(richness_no_zeros_rarified)



# view richness across tidal heights
richness_no_zeros_by_level %>%
  ggplot(aes(x = tidalheight, 
             y = mean_richness)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)


richness_no_zeros_by_level %>%
  group_by(tidalheight) %>% 
  summarize(meanrich = mean(mean_richness)) %>%
  arrange(desc(meanrich))
# in both cases, tidal height 0.696 is the height with the highest 
# mean species richness. I think we should recenter around that value,
# at least for the quadratic model, to accurately capture the peak


# build models ------------------------------------------------------------
richness_no_zeros_by_level %>%
  ggplot(aes(x = mean_richness)) +
  geom_histogram() +
  facet_wrap(~level)
richness_no_zeros_by_level %>% 
  ggplot(aes(x = mean_richness)) +
  geom_histogram()

# | Mod 1. by year ---------------------------------------------------
# here, we're going to use a Gamma distribution
# even though out richness variable isn't an integer, 
# it shouldn't matter here, because both predictors are continuous
mod1 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year",
  family = Gamma(link = "log"))
summary(mod1)
AIC(mod1)
# no significance of predictors,
# AIC = 1385.353

# | Mod 2. by year + height ---------------------------------------------------
mod2 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year + tidalheight",
  family = Gamma(link = "log"))
summary(mod2)
AIC(mod2)
# significance of tidal height, not year
# AIC = 1212.117

# | Mod 3. by year * height ---------------------------------------------------
mod3 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year * tidalheight",
  family = Gamma(link = "log"))
summary(mod3)
AIC(mod3)
# significance year, tidal height, and interaction
# AIC = 1209.605

# | Mod 4. by year * height^2 ---------------------------------------------------
mod4 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ year * I(tidalheight^2)",
  family = Gamma(link = 'log'))
summary(mod4)
AIC(mod4)
# significance of year, tidal height and interaction
# AIC = 1102.042
# rm(mod4)
# I suppose best so far

# | Mod 5. by year + height^2 ---------------------------------------------------
mod5 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ year + I(tidalheight^2)",
  family = Gamma(link = "log"))
summary(mod5)
AIC(mod5)
# significance of tidal height, not year
# AIC =  1108.072




# | compare models ----------------------------------------------------------
AIC(mod1, mod2, mod3, mod4, mod5) %>% arrange(AIC)
MuMIn::AICc(mod1, mod2, mod3, mod4) %>% arrange(AICc)
performance::compare_performance(
  mod1, mod2, mod3, mod4, mod5,
  metrics = "all",
  rank =T)
# model 4 wins in all metrics 
#plot(mod4)
#check_model(mod4) # use this later
DHARMa::simulateResiduals(mod4) |> 
  DHARMa::plotQQunif() # Jarrett says this is ok


# predict and plot --------------------------------------------------------
predict_df <- data.frame(
  tidalheight = rep(unique(richness_no_zeros_by_level$tidalheight),
                    times = length(unique(richness_no_zeros_by_level$year))) #-peak_tidal_height
  ,
  year = rep(unique(richness_no_zeros_by_level$year),
             each = length(unique(richness_no_zeros_by_level$tidalheight)))
)

# package for getting prediction intervals from Gamma distributions
# https://search.r-project.org/CRAN/refmans/EnvStats/html/predIntGamma.html
# or Jarrett's simulation-based workflow: in combo with ggdisc
# https://github.com/jebyrnes/sinterval


no_zero_augment <- broom::augment(mod4,
                                  newdata = predict_df,
                                  type.predict = "response",
                                  se_fit = T,
                                  simulate_pi  = T) %>%
  # mutate(tidalheight = tidalheight_centered + peak_tidal_height) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 
# how / if to get the lower/upper bound

p1 <- no_zero_augment %>% 
  left_join(richness_no_zeros_by_level) %>% 
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
  coord_cartesian(ylim = c(0, NA )) +
  
  labs(x = NULL,
       y = "Rarified Species Richness",
       title = "Richness Changes Over Time"
  ) +
  geom_text(
    data = . %>%
      group_by(tidalheight, tidalheight_f2) %>%
      summarize(change = round(last(.fitted) - first(.fitted),2)),
    aes(x = 1982,
        y = 11,
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

p1 + ggview::canvas(6,5)
ggsave(p1, 
       filename = "outputs/richness_change/p1_no_zeros.png",
       width = 6,
       height = 5,
       unit = "in")



p2 <- no_zero_augment %>% 
  left_join(richness_no_zeros_by_level)  %>%
  arrange(tidalheight, year) %>%
  ggplot(aes(x = .fitted,
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
       title = "Richness Across Shore Levels") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(-.1,10.5),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank())

p2 + ggview::canvas(4,4)

ggsave(p2, 
       filename = "outputs/richness_change/p2.2.2_no_zeros.png",
       width = 4,
       height = 4,
       unit = "in")

# save RDS version so I can merge with other richness plot
saveRDS(p2,
        "outputs/richness_change/p2.2.2_no_zeros.rds")

p3 <- no_zero_augment %>%
  left_join(richness_no_zeros_by_level) %>%
  #filter(year %in% c(1982, 1990, 2000, 2010, 2020)) %>%
  arrange(tidalheight) %>%
  ggplot(aes(x = .fitted, y = tidalheight,
             color = year, 
             fill = year, 
             group = year)) +
  geom_point(aes(x = mean_richness),
             position = position_jitter(width = 0,
                                        height = .1),
             size = 1.5,
             stroke = .1,
             color = "black",
             shape = 21,
             alpha = .6, show.legend = F) +
  geom_path(linewidth = .5) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  scale_x_continuous(limits = c(0,10), 
                     breaks = scales::pretty_breaks()) + 
  scale_y_continuous(limits = c(-1,5), 
                     breaks = scales::pretty_breaks()) +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Richness Across Shore Levels") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  coord_cartesian(expand = F) +
  theme(panel.grid = element_blank())

p3 + ggview::canvas(4,4)
ggsave(p3, 
       filename = "outputs/richness_change/p3_no_zeros.png",
       width = 4,
       height = 4,
       unit = "in")
# remove all stuff we don't need
rm(mod1, mod2, mod3, mod5, predict_df)



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


# find richness with presence only ----------------------------------------
# here, we omit quadrats that say they were sampled,
# but have a richness of zero. Instead, we'll use just quadrats that have
# some richness value. (below we'll include the zeros too)
richness_no_zeros_hs <- p_only_hs %>%
  group_by(year, transect, level, tidalheight, replicate) %>%
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
  group_by(year, transect, tidalheight, level) %>%
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
  group_by(year,level,tidalheight) %>%
  summarize(mean_richness = mean(rarified_richness),
            n_transects_sampled = n(),
            .groups = "drop") %>%
  filter(level >= 0,
         level <= 15)

# quick plot to view 
richness_no_zeros_by_level_hs %>%# filter(n_samples >= 3) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) %>%
  ggplot(aes(x = year,
             y = mean_richness)) +
  geom_point() +
  facet_wrap(~tidalheight_f) +
  geom_smooth(method = "lm")
rm(richness_no_zeros_rarified_hs)



# view richness across tidal heights
richness_no_zeros_by_level_hs %>%
  ggplot(aes(x = tidalheight, 
             y = mean_richness)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)


richness_no_zeros_by_level_hs %>%
  group_by(tidalheight) %>% 
  summarize(meanrich = mean(mean_richness)) %>%
  arrange(desc(meanrich))
# in both cases, tidal height 0.696 is the height with the highest 
# mean species richness. I think we should recenter around that value,
# at least for the quadratic model, to accurately capture the peak


# build models ------------------------------------------------------------
richness_no_zeros_by_level_hs %>%
  ggplot(aes(x = mean_richness)) +
  geom_histogram() +
  facet_wrap(~level)
richness_no_zeros_by_level_hs %>% 
  ggplot(aes(x = mean_richness)) +
  geom_histogram()

# | Mod 1. by year ---------------------------------------------------
# here, we're going to use a Gamma distribution
# even though out richness variable isn't an integer, 
# it shouldn't matter here, because both predictors are continuous
mod1_hs <- glm(
  data = richness_no_zeros_by_level_hs,
  formula = "mean_richness ~ year",
  family = Gamma(link = "log"))
summary(mod1_hs)
AIC(mod1_hs)
# no significance of predictors,
# AIC = 1385.353

# | Mod 2. by year + height ---------------------------------------------------
mod2_hs <- glm(
  data = richness_no_zeros_by_level_hs,
  formula = "mean_richness ~ year + tidalheight",
  family = Gamma(link = "log"))
summary(mod2_hs)
AIC(mod2_hs)
# significance of tidal height, not year
# AIC = 1212.117

# | Mod 3. by year * height ---------------------------------------------------
mod3_hs <- glm(
  data = richness_no_zeros_by_level_hs,
  formula = "mean_richness ~ year * tidalheight",
  family = Gamma(link = "log"))
summary(mod3_hs)
AIC(mod3_hs)
# significance year, tidal height, and interaction
# AIC = 1209.605

# | Mod 4. by year * height^2 ---------------------------------------------------
mod4_hs <- glm(
  data = richness_no_zeros_by_level_hs ,
  formula = "mean_richness ~ year * I(tidalheight^2)",
  family = Gamma(link = 'log'))
summary(mod4_hs)
AIC(mod4_hs)
# significance of year, tidal height and interaction
# AIC = 1102.042
# rm(mod4)
# I suppose best so far

# | Mod 5. by year + height^2 ---------------------------------------------------
mod5_hs <- glm(
  data = richness_no_zeros_by_level_hs ,
  formula = "mean_richness ~ year + I(tidalheight^2)",
  family = Gamma(link = "log"))
summary(mod5_hs)
AIC(mod5_hs)
# significance of tidal height, not year
# AIC =  1108.072




# | compare models ----------------------------------------------------------
AIC(mod1_hs, mod2_hs, mod3_hs, mod4_hs, mod5_hs) %>% arrange(AIC)
MuMIn::AICc(mod1_hs, mod2_hs, mod3_hs, mod4_hs) %>% arrange(AICc)
performance::compare_performance(
  mod1_hs, mod2_hs, mod3_hs, mod4_hs, mod5_hs,
  metrics = "all",
  rank =T)
# model 4 wins in all metrics 
#plot(mod4)
#check_model(mod4) # use this later
DHARMa::simulateResiduals(mod4_hs) |> 
  DHARMa::plotQQunif() # Jarrett says this is ok


# predict and plot --------------------------------------------------------
predict_df_hs <- data.frame(
  tidalheight = rep(unique(richness_no_zeros_by_level_hs$tidalheight),
                    times = length(unique(richness_no_zeros_by_level_hs$year))) #-peak_tidal_height
  ,
  year = rep(unique(richness_no_zeros_by_level_hs$year),
             each = length(unique(richness_no_zeros_by_level_hs$tidalheight)))
)

# package for getting prediction intervals from Gamma distributions
# https://search.r-project.org/CRAN/refmans/EnvStats/html/predIntGamma.html
# or Jarrett's simulation-based workflow: in combo with ggdisc
# https://github.com/jebyrnes/sinterval


no_zero_augment_hs <- broom::augment(mod4_hs,
                                  newdata = predict_df_hs,
                                  type.predict = "response",
                                  se_fit = T,
                                  simulate_pi  = T) %>%
  # mutate(tidalheight = tidalheight_centered + peak_tidal_height) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 
# how / if to get the lower/upper bound

p1_hs <- no_zero_augment_hs %>% 
  left_join(richness_no_zeros_by_level_hs) %>% 
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
  coord_cartesian(ylim = c(0, NA )) +
  
  labs(x = NULL,
       y = "Rarified Species Richness",
       title = "Richness Changes Over Time"
  ) +
  geom_text(
    data = . %>%
      group_by(tidalheight, tidalheight_f2) %>%
      summarize(change = round(last(.fitted) - first(.fitted),2)),
    aes(x = 1982,
        y = 11,
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

p1_hs + ggview::canvas(6,5)
ggsave(p1, 
       filename = "outputs/richness_change/p1_no_zeros_HS.png",
       width = 6,
       height = 5,
       unit = "in")



p2_hs <- no_zero_augment_hs %>% 
  left_join(richness_no_zeros_by_level_hs)  %>%
  arrange(tidalheight, year) %>%
  ggplot(aes(x = .fitted,
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
       title = "Richness Across Shore Levels") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(-.1,10.5),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank())

p2_hs + ggview::canvas(4,4)

ggsave(p2_hs, 
       filename = "outputs/richness_change/p2.2.2_no_zeros_HS.png",
       width = 4,
       height = 4,
       unit = "in")

p3_hs <- no_zero_augment_hs %>%
  left_join(richness_no_zeros_by_level_hs) %>%
  #filter(year %in% c(1982, 1990, 2000, 2010, 2020)) %>%
  arrange(tidalheight) %>%
  ggplot(aes(x = .fitted, y = tidalheight,
             color = year, 
             fill = year, 
             group = year)) +
  geom_point(aes(x = mean_richness),
             position = position_jitter(width = 0,
                                        height = .1),
             size = 1.5,
             stroke = .1,
             color = "black",
             shape = 21,
             alpha = .6, show.legend = F) +
  geom_path(linewidth = .5) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  scale_x_continuous(limits = c(0,10), 
                     breaks = scales::pretty_breaks()) + 
  scale_y_continuous(limits = c(-1,5), 
                     breaks = scales::pretty_breaks()) +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Effort-Corrected Species Richness",
       y = "Tidal Height (m)",
       color = "Sample Year",
       title = "Richness Across Shore Levels") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  coord_cartesian(expand = F) +
  theme(panel.grid = element_blank())

p3_hs + ggview::canvas(4,4)
ggsave(p3_hs, 
       filename = "outputs/richness_change/p3_no_zeros_HS.png",
       width = 4,
       height = 4,
       unit = "in")
# remove all stuff we don't need
rm(mod1_hs, mod2_hs, mod3_hs, mod5_hs, predict_df_hs)



# add HS line to full plot ------------------------------------------------

all_rich_p <- (p2 + theme(legend.position = "none") +
   labs(title = "a) All Data")| p2_hs + labs(title = "b) Highly Sampled Data") )+
  plot_layout(#guides = "collect",
              axes = "collect")


all_rich_p + ggview::canvas(8,4)

ggsave(all_rich_p,
       file = "outputs/richness_change/rich_per_level_all_vs_hs.png",
       width = 8,
       height = 4,
       unit = "in")
