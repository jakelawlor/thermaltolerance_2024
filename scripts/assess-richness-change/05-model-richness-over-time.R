# track changes in species richness over time

# here, we will test species richness shifts over time within the 
# entire study area, as well as across depths 



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(broom)
theme_set(theme_bw())


# data --------------------------------------------------------------------
# first, upload presence only data
p_only <- readr::read_csv(
  here::here("data-processed",
             "spp_pres_by_replicate.csv")
) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  arrange(tidalheight) 

# then upload transects / years / tidalheights where
# both count and cover were marked as "data taken"
# we'll use this later to back-fill zeros when data was 
# taken but there is no richness value
both_taken <- readr::read_csv(
 here::here("data-processed",
            "quadrats_data_taken.csv")
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

# find richness with presence and absence ---------------------------------
# here, we're going to add in zeros when we know data was marked as taken,
# but no species were present

# Use the richness_no_zeros dataset above, but join with the data taken
# dataset so that the sites where data was taken but no species were present
# will be marked as zero richness
richness_with_zeros <- richness_no_zeros %>%
  full_join(both_taken) %>%
  tidyr::replace_na(replace = list(richness = 0)) 
richness_with_zeros %>% count(richness == 0)
richness_with_zeros %>%
  ggplot(aes(x = year, y = tidalheight,
             fill = richness > 0)) +
  geom_tile(alpha = .5) +
  facet_wrap(~transect)
# that adds about 2608 sites where richness is 0, of 9809 sites total
# about 26 percent of sites

# rarify quadrats by replicate
richness_with_zeros_rarified <- 
  richness_with_zeros %>%
  # find mean richness between replicates 
  group_by(year, transect, tidalheight, level) %>%
  summarize(replicates = n(),
            rarified_richness = mean(richness))
rm(richness_with_zeros, richness_no_zeros)


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


# do the same thing with the dataset including the zeros
richness_with_zeros_by_level <- 
  richness_with_zeros_rarified %>%
  group_by(year,level, tidalheight) %>%
  summarize(mean_richness = mean(rarified_richness),
            n_transects_sampled = n(),
            .groups = "drop") %>%
  filter(level >= 0,
         level <= 15)
rm(richness_with_zeros_rarified)

# again, view
richness_with_zeros_by_level %>% #filter(n_samples >= 3) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) %>%
  ggplot(aes(x = year,
             y = mean_richness)) +
  geom_point() +
  facet_wrap(~tidalheight_f) +
  geom_smooth(method = "lm")
# flipping back and forth between these plots shows that
# the trends aren't really that different, but a few outlier
# points get fixed a bit when the zeros are counted as well.

# view richness across tidal heights
richness_no_zeros_by_level %>%
  ggplot(aes(x = tidalheight, 
             y = mean_richness)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)

richness_with_zeros_by_level %>%
  ggplot( aes(x = tidalheight, 
              y = mean_richness)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)
# in both cases, there's a definite nonlinear trend of richness across tidal height,
# with the peak of richness being about .75m. 

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
# note: ask Jarrett what to do / think about data distributions
# weird distribution using all data, but more normal distribution
# when stratifying by tidal height


# | Mod 1. by year ---------------------------------------------------
# here, we're going to use a Gamma distribution
# even though out richness variable isn't an integer, 
# it shouldn't matter here, because both predictors are continuous
mod1 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year",
  family = Gamma(link = "log"))
summary(mod1)
# no significance of predictors,
# AIC = 2533

# | Mod 2. by year + height ---------------------------------------------------
mod2 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year + tidalheight",
  family = Gamma(link = "log"))
summary(mod2)
# significance of tidal height, not year
# AIC = 1886.9

# | Mod 3. by year * height ---------------------------------------------------
mod3 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year * tidalheight",
  family = Gamma(link = "log"))
summary(mod3)
# significance year, tidal height, and interaction
# AIC = 1878.6

# | Mod 4. by year * height^2 ---------------------------------------------------
mod4 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ year * I(tidalheight^2)",
  family = Gamma(link = 'log'))
summary(mod4)
# significance of tidal height and interaction
# AIC = 1591
# rm(mod4)
# I suppose best so far

# | Mod 5. by year + height^2 ---------------------------------------------------
mod5 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ year + I(tidalheight^2)",
  family = Gamma(link = "log"))
summary(mod5)
# significance of tidal height, not year
# AIC = 1597.1

# | Mod 6. by year + height_centered^2 ---------------------------------------------------
# rescale tidal height values
peak_tidal_height <- richness_no_zeros_by_level %>%
  group_by(tidalheight) %>%
  summarize(mean_rich = mean(mean_richness)) %>%
  arrange(desc(mean_rich)) %>% 
  slice(1) %>%
  pull(tidalheight)
richness_no_zeros_by_level$tidalheight_centered <- 
  richness_no_zeros_by_level$tidalheight - peak_tidal_height

richness_no_zeros_by_level %>%
  ggplot(aes(x = tidalheight_centered, 
             y = mean_richness)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)

mod6 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ year * I(tidalheight_centered^2)",
  family = Gamma(link = "log"))
summary(mod6)
# AIC = 1614.4

# | Mod 7. by year + height_centered^2 ---------------------------------------------------
mod7 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ year + I(tidalheight_centered^2)",
  family = Gamma(link = "log"))
summary(mod7)
# AIC = 1618.5

# | Mod 8. by tidal height only ---------------------------------------------------
mod8 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~  I(tidalheight^2)",
  family = Gamma(link = "log"))
summary(mod8)
# AIC = 1597.3



# | compare models ----------------------------------------------------------
AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8) %>% arrange(AIC)
MuMIn::AICc(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8) %>% arrange(AICc)
performance::compare_performance(
  mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8,
  metrics = "all",
  rank =T)
# model 4 wins in all metrics 
plot(mod4)
check_model(mod4) # use this later
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

# plot change in richness by tidalheight
p2 <- no_zero_augment %>% 
  left_join(richness_no_zeros_by_level)  %>%
  ggplot(aes(x = .fitted,
             y = tidalheight,
             group = tidalheight_f,
             fill = year)) +
  geom_point(aes(x = mean_richness),
             position = 
               position_jitter(width = 0,
                               height = .05),
             size = 1.5,
             stroke = .1,
             shape = 21,
             alpha = .75) +
  scale_fill_viridis_c(option = "viridis",
                        limits = c(1980,2020)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            position = position_nudge(y = .07),
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
       fill = "Sample Year",
       title = "Richness Change Across Tidalheights") +
  guides(fill = guide_colorbar(
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
       filename = "outputs/richness_change/p2_no_zeros.png",
       width = 4,
       height = 4,
       unit = "in")


p2.2 <- no_zero_augment %>% 
  left_join(richness_no_zeros_by_level)  %>%
  ggplot(aes(x = .fitted,
             y = tidalheight,
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .75) +
  scale_color_viridis_c(option = "viridis",
                        limits = c(1980,2020)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            position = position_nudge(y = .07),
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
       title = "Richness Change Across Tidalheights") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(-.1,10.5),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank())


p2.2 + ggview::canvas(4,4)
ggsave(p2.2, 
       filename = "outputs/richness_change/p2.2_no_zeros.png",
       width = 4,
       height = 4,
       unit = "in")

p2.2.2 <- no_zero_augment %>% 
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

p2.2.2 + ggview::canvas(4,4)

ggsave(p2.2.2, 
       filename = "outputs/richness_change/p2.2.2_no_zeros.png",
       width = 4,
       height = 4,
       unit = "in")

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
rm(mod1, mod2, mod3, mod5, mod6, mod7, mod8, predict_df)






# continue with zeros ------------------------------------------------------
# NOTE: Kylla says don't trust the "data taken", has a hard time
# imagining sites as 0 richness. Instead, stick with the presence only
# dataset, which is the more conservative way to do this, in this case. 

# change zeros to 0.1 so we can use gamma distribution
richness_with_zeros_by_level2 <- richness_with_zeros_by_level %>%
  mutate(mean_richness = case_when(mean_richness == 0 ~ 0.01,
                                   TRUE ~ mean_richness))

rm(richness_with_zeros_by_level)
# | mod 1. by year ---------------------------------------------------------
# here, we're going to use a poisson distribution
# even though out richness variable isn't an integer, 
# it shouldn't matter here, because both predictors are continuous
mod1_w_zero <- glm(
  data = richness_with_zeros_by_level2,
  formula = "mean_richness ~ year",
  family = Gamma(link = "log"))
summary(mod1_w_zero)
# no significance of predictors,
# AIC = 2718.6

# | mod 2. by year + tidalheight ---------------------------------------------------------
mod2_w_zero <- glm(
  data = richness_with_zeros_by_level2,
  formula = "mean_richness ~ year + tidalheight",
  family = Gamma(link = "log"))
summary(mod2_w_zero)
# significance of tidal height and year
# AIC = 2227.4

# | mod 3. by year * tidalheight ---------------------------------------------------------
mod3_w_zero <- glm(
  data = richness_with_zeros_by_level2,
  formula = "mean_richness ~ year * tidalheight",
  family = Gamma(link = "log"))
summary(mod3_w_zero)
# significance year, tidal height, and interaction
# AIC = 2211.3

# | mod 4. by year * tidalheight^2 ---------------------------------------------------------
mod4_w_zero <- glm(
  data = richness_with_zeros_by_level2,
  formula = "mean_richness ~ year * I(tidalheight^2)",
  family = Gamma)
summary(mod4_w_zero)
# significance of tidal height and interaction
# AIC = 2285.4

# | mod 4. by year + tidalheight^2 ---------------------------------------------------------
mod5_w_zero <- glm(
  data = richness_with_zeros_by_level2,
  formula = "mean_richness ~ year + I(tidalheight^2)",
  family = Gamma)
summary(mod5_w_zero)
# significance of tidal height and interaction
# AIC = 2288.7

# | Mod 6. by year + height_centered^2 ---------------------------------------------------
# rescale tidal height values
peak_tidal_height <- richness_with_zeros_by_level2 %>%
  group_by(tidalheight) %>%
  summarize(mean_rich = mean(mean_richness)) %>%
  arrange(desc(mean_rich)) %>% 
  slice(1) %>%
  pull(tidalheight)
richness_with_zeros_by_level2$tidalheight_centered <- 
  richness_with_zeros_by_level2$tidalheight - peak_tidal_height

richness_with_zeros_by_level2 %>%
  ggplot(aes(x = tidalheight_centered, 
             y = mean_richness)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)

mod6_w_zero <- glm(
  data = richness_with_zeros_by_level2 ,
  formula = "mean_richness ~ year * I(tidalheight_centered^2)",
  family = Gamma(link = "log"))
summary(mod6_w_zero)
# AIC = 1788.2
# NOTE: for jarrett, qqplot looks really bad

# | Mod 7. by year + height_centered^2 ---------------------------------------------------
mod7_w_zero <- glm(
  data = richness_with_zeros_by_level2 ,
  formula = "mean_richness ~ year + I(tidalheight_centered^2)",
  family = Gamma(link = "log"))
summary(mod7_w_zero)
# AIC = 1821.1

# | Mod 8. by tidal height only ---------------------------------------------------
mod8_w_zero <- glm(
  data = richness_with_zeros_by_level2 ,
  formula = "mean_richness ~  I(tidalheight^2)",
  family = Gamma(link = "log"))
summary(mod8_w_zero)
# AIC = 1897.5


# | compare models ----------------------------------------------------------
AIC(mod1_w_zero, mod2_w_zero, mod3_w_zero, mod4_w_zero, mod5_w_zero, 
    mod6_w_zero, mod7_w_zero, mod8_w_zero) %>% arrange(AIC)
MuMIn::AICc(mod1_w_zero, mod2_w_zero, mod3_w_zero, mod4_w_zero,
            mod5_w_zero, mod6_w_zero, mod7_w_zero, mod8_w_zero) %>% arrange(AICc)
performance::compare_performance(
  mod1_w_zero, mod2_w_zero, mod3_w_zero, mod4_w_zero, 
  mod5_w_zero, mod6_w_zero, mod7_w_zero, mod8_w_zero,
  rank = T)
# mod 6 wins by all metrics


# predict and plot --------------------------------------------------------
predict_df_zero <- data.frame(
  tidalheight_centered = rep(unique(richness_with_zeros_by_level2$tidalheight),
                    times = length(unique(richness_with_zeros_by_level2$year))) -peak_tidal_height
  ,
  year = rep(unique(richness_with_zeros_by_level2$year),
             each = length(unique(richness_with_zeros_by_level2$tidalheight)))
)

with_zero_augment <- broom::augment(mod6_w_zero,
                                  newdata = predict_df_zero,
                                  type.predict = "response",
                                  se_fit = T) %>%
  mutate(tidalheight = tidalheight_centered + peak_tidal_height) %>%
  mutate(tidalheight = round(tidalheight,3)) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) 


p1_w_zero <- with_zero_augment %>% 
  left_join(richness_with_zeros_by_level2 %>% mutate(tidalheight = round(tidalheight,3))) %>%
  ggplot(aes(x = year,
             y = .fitted,
             color = tidalheight_f,
             group = tidalheight_f,
             fill = tidalheight_f)) +
  geom_point(aes(y = mean_richness),
             shape = 21,
             alpha = .8, 
             color = "black") +
  gghighlight::gghighlight(
    unhighlighted_params = list(fill = "grey75",
                                alpha = .1),
    use_direct_label = F) +
  # geom_ribbon(aes(ymin = .lower,
  #                 ymax = .upper),
  #             alpha = .25) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(.~paste0(tidalheight_f,"m")) +
  theme_bw() + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, NA )) +
  
  labs(x = NULL,
       y = "Rarified Species Richness",
       title = "Richness Changes Over Time"
  ) +
  geom_text(
    data = with_zero_augment %>%
      group_by(tidalheight, tidalheight_f) %>%
      summarize(change = round(last(.fitted) - first(.fitted),2)),
    aes(x = 1982,
        y = 11,
        label = change),
    hjust = 0, vjust = 1,
    color = "black",
    fontface = "bold",
    stat = "unique"
  ) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust  = 1))

p1_w_zero + ggview::canvas(6,5)
ggsave(p1_w_zero, 
       filename = "outputs/richness_change/p1_with_zeros.png",
       width = 6,
       height = 5,
       unit = "in")


# plot change in richness by tidalheight
p2_w_zero <- with_zero_augment %>% 
  left_join(richness_with_zeros_by_level2 %>% mutate(tidalheight = round(tidalheight,3))) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight,
             group = tidalheight_f,
             fill = year)) +
  geom_point(aes(x = mean_richness),
             position = position_jitter(width = 0,
                                        height = .02),
             size = 1.5,
             stroke = .1,
             shape = 21,
             alpha = .6) +
  scale_fill_viridis_c(option = "viridis",
                       limits = c(1980,2020)) +
  geom_path(arrow = arrow(length = unit(5,"pt"),
                          type = "closed"),
            position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  theme(legend.position = c(.97,.98),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  labs(x = "Species Richness",
       y = "Tidal Height (m)",
       fill = "Sample Year",
       title = "Richness Change Across Tidalheights") +
  guides(fill = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(-.1,10.5),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank())


p2.2_w_zero <- with_zero_augment %>% 
  left_join(richness_with_zeros_by_level2 %>% mutate(tidalheight = round(tidalheight,3))) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight,
             group = tidalheight_f,
             color = year)) +
  geom_point(aes(x = mean_richness),
             size = 3,
             stroke = 2.5,
             shape = "|",
             alpha = .75) +
  scale_color_viridis_c(option = "viridis",
                        limits = c(1980,2020)) +
  geom_path(arrow = arrow(length = unit(3,"pt"),
                          type = "closed"),
            position = position_nudge(y = .07),
            color = "black",
            linewidth = .75,
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
       title = "Richness Change Across Tidalheights") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  coord_cartesian(xlim = c(-.1,10.5),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank())



p2_w_zero + ggview::canvas(4,4)
ggsave(p2_w_zero, 
       filename = "outputs/richness_change/p2_with_zeros.png",
       width = 4,
       height = 4,
       unit = "in")
p2.2_w_zero + ggview::canvas(4,4)
ggsave(p2.2_w_zero, 
       filename = "outputs/richness_change/p2.2_with_zeros.png",
       width = 4,
       height = 4,
       unit = "in")

p3_w_zero <- with_zero_augment %>%
  left_join(richness_with_zeros_by_level2 %>% mutate(tidalheight = round(tidalheight,3))) %>%
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
       title = "Richness Change Across Tidalheights") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  coord_cartesian(expand = F) +
  theme(panel.grid = element_blank())

p3_w_zero + ggview::canvas(4,4)
ggsave(p3_w_zero, 
       filename = "outputs/richness_change/p3_with_zeros.png",
       width = 4,
       height = 4,
       unit = "in")
# remove all stuff we don't need
rm(mod1_w_zero, mod2_w_zero, mod3_w_zero, mod4_w_zero,
   mod5_w_zero, mod7_w_zero, mod8_w_zero, predict_df_w_zero)




# scratch -----------------------------------------------------------------
library(glmmTMB)
mod5_w_zero <- glmmTMB(
  data = richness_with_zeros_by_level,
  formula = mean_richness ~ year * tidalheight, 
  family=ziGamma(link = "log"), 
  ziformula= ~.)
summary(mod5_w_zero)
AIC(mod4_w_zero)

