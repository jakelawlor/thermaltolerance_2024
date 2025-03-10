# track changes in species richness over time - by GROUP

# here, we will test species richness shifts over time within the 
# entire study area, as well as across depths 



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(broom)
theme_set(theme_bw())
traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed.csv")


# data --------------------------------------------------------------------
# first, upload presence only data
p_only <- readr::read_csv(
  here::here("data-processed",
             "spp_pres_by_replicate.csv")
) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  arrange(tidalheight) %>%
  left_join(traits %>% rename(organism = gen_spp) %>%
              mutate(group = stringr::str_to_title(motility_adult)))

#p_only <- readr::read_csv("data-processed/cti-data/cti-data-highly-sampled.csv") %>%
#  left_join(traits %>% rename(organism = gen_spp) %>%
#              mutate(group = stringr::str_to_title(motility_adult)))


# then upload transects / years / tidalheights where
# both count and cover were marked as "data taken"
# we'll use this later to back-fill zeros when data was 
# taken but there is no richness value
both_taken <- readr::read_csv(
  here::here("data-processed",
             "quadrats_data_taken.csv")
) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  arrange(tidalheight) %>%
  mutate(group1 = unique(p_only$group)[1],
         group2 = unique(p_only$group)[2]) %>%
  tidyr::pivot_longer(cols = c(group1, group2),
                      names_to = "delete",
                      values_to = "group") %>%
  select(-delete)


# find richness with presence only ----------------------------------------
# here, we omit quadrats that say they were sampled,
# but have a richness of zero. Instead, we'll use just quadrats that have
# some richness value. (below we'll include the zeros too)
richness_no_zeros <- p_only %>%
  group_by(year, transect, level, tidalheight, replicate, group) %>%
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
  group_by(year, transect, tidalheight, level, group) %>%
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
  full_join(both_taken ) %>%
  tidyr::replace_na(replace = list(richness = 0)) 
richness_with_zeros %>% count(richness == 0)
richness_with_zeros %>%
  ggplot(aes(x = year, y = tidalheight,
             fill = richness > 0)) +
  geom_tile(alpha = .5) +
  facet_wrap(group~transect)
# that adds about 2608 sites where richness is 0, of 9809 sites total
# about 26 percent of sites

# rarify quadrats by replicate
richness_with_zeros_rarified <- 
  richness_with_zeros %>%
  # find mean richness between replicates 
  group_by(year, transect, tidalheight, level, group) %>%
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
  group_by(year,level,tidalheight, group) %>%
  summarize(mean_richness = mean(rarified_richness),
            n_transects_sampled = n(),
            .groups = "drop") %>%
  filter(level >= 0,
         level <= 15)

# quick plot to view 
richness_no_zeros_by_level %>%# filter(n_samples >= 3) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight), -tidalheight)) %>%
  ggplot(aes(x = year,
             y = mean_richness, 
             color = group)) +
  geom_point() +
  facet_wrap(~tidalheight_f) +
  geom_smooth(method = "lm")
rm(richness_no_zeros_rarified)


# do the same thing with the dataset including the zeros
richness_with_zeros_by_level <- 
  richness_with_zeros_rarified %>%
  group_by(year,level, tidalheight,group) %>%
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
             y = mean_richness,
             color = group)) +
  geom_point() +
  facet_wrap(~tidalheight_f) +
  geom_smooth(method = "lm")
# flipping back and forth between these plots shows that
# the trends aren't really that different, but a few outlier
# points get fixed a bit when the zeros are counted as well.

# view richness across tidal heights
richness_no_zeros_by_level %>%
  ggplot(aes(x = tidalheight, 
             y = mean_richness,
             color = group)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)

richness_with_zeros_by_level %>%
  ggplot( aes(x = tidalheight, 
              y = mean_richness,
              color = group)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE)
# in both cases, there's a definite nonlinear trend of richness across tidal height,
# with the peak of richness being about .75m. 

richness_no_zeros_by_level %>%
  group_by(tidalheight, group) %>% 
  summarize(meanrich = mean(mean_richness)) %>%
  arrange(desc(meanrich))
# in both cases, tidal height 0.696 is the height with the highest 
# mean species richness. I think we should recenter around that value,
# at least for the quadratic model, to accurately capture the peak


# build models ------------------------------------------------------------
richness_no_zeros_by_level %>%
  ggplot(aes(x = mean_richness,
             fill = group)) +
  geom_histogram(alpha = .5,
                 position = 'identity') +
  facet_wrap(~level)
richness_no_zeros_by_level %>% 
  ggplot(aes(x = mean_richness,
             fill = group)) +
  geom_histogram(alpha = .5,
                 position = 'identity')
# note: ask Jarrett what to do / think about data distributions
# weird distribution using all data, but more normal distribution
# when stratifying by tidal height


# | Mod 1. by year ---------------------------------------------------
# here, we're going to use a Gamma distribution
# even though out richness variable isn't an integer, 
# it shouldn't matter here, because both predictors are continuous
mod1 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year * group",
  family = Gamma(link = "log"))
summary(mod1)
AIC(mod1)
# all significance of predictors,
# AIC = 3006.37

# | Mod 2. by year + height ---------------------------------------------------
mod2 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ (year + tidalheight) * group",
  family = Gamma(link = "log"))
summary(mod2)
AIC(mod2)
# all significant predictors
# AIC = 2141.548

# | Mod 3. by year * height ---------------------------------------------------
mod3 <- glm(
  data = richness_no_zeros_by_level,
  formula = "mean_richness ~ year * tidalheight * group",
  family = Gamma(link = "log"))
summary(mod3)
AIC(mod3)
# significance year, tidal height, and interaction
# AIC = 2132.761

# | Mod 4. by year * height^2 ---------------------------------------------------
mod4 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ (year * I(tidalheight^2))*group",
  family = Gamma(link = 'log'))
summary(mod4)
# significance of lots
AIC(mod4)
# AIC =  1842.096
# rm(mod4)
# I suppose best so far

# | Mod 5. by year + height^2 ---------------------------------------------------
mod5 <- glm(
  data = richness_no_zeros_by_level ,
  formula = "mean_richness ~ (year + I(tidalheight^2)) * group",
  family = Gamma(link = "log"))
summary(mod5)
# significance of lots
AIC(mod5)
# AIC =  1852.966

# | compare models ----------------------------------------------------------
#AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8) %>% arrange(AIC)
#MuMIn::AICc(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8) %>% arrange(AICc)
performance::compare_performance(
  mod1, mod2, mod3, mod4, mod5, #mod6, mod7, mod8,
  metrics = "all",
  rank =T)
# model 4 wins in all metrics 
plot(mod4)
#check_model(mod4) # use this later
#DHARMa::simulateResiduals(mod4) |> 
#  DHARMa::plotQQunif() # Jarrett says this is ok


# predict and plot --------------------------------------------------------
predict_df <- data.frame(
  tidalheight = rep(unique(richness_no_zeros_by_level$tidalheight),
                    times = length(unique(richness_no_zeros_by_level$year))) #-peak_tidal_height
  ,
  year = rep(unique(richness_no_zeros_by_level$year),
             each = length(unique(richness_no_zeros_by_level$tidalheight)))
)
predict_df <- rbind(predict_df %>% mutate(group = unique(richness_no_zeros_by_level$group)[1]),
                    predict_df %>% mutate(group = unique(richness_no_zeros_by_level$group)[2]))

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
             color = group,
             fill = group)) +
  geom_point(aes(y = mean_richness),
             shape = 21,
             alpha = .8, 
             stroke = .1,
             color = "black") +
 # gghighlight::gghighlight(
 #   unhighlighted_params = list(fill = "grey75",
 #                               alpha = .1),
 #   use_direct_label = F) +
  geom_ribbon(aes(ymin = .fitted - 2*.se.fit,
                  ymax = .fitted + 2*.se.fit),
              alpha = .75) +
  geom_line(linewidth = .5,
            color = "black") +
  
  facet_wrap(.~tidalheight_f2) +
  theme_bw() + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, NA )) +
  
  labs(x = NULL,
       y = "Rarified Species Richness",
       title = "Richness Changes Over Time"
  ) +
# geom_text(
#   data = . %>%
#     group_by(tidalheight, tidalheight_f2, group) %>%
#     summarize(change = round(last(.fitted) - first(.fitted),2)),
#   aes(x = 1982,
#       y = 11,
#       label = change),
#   hjust = 0, 
#   vjust = 1,
#   size =3.5,
#   color = "black",
#   fontface = "bold",
#   stat = "unique"
# ) +
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
#ggsave(p1, 
#       filename = "outputs/richness_change/p1_no_zeros.png",
#       width = 6,
#       height = 5,
#       unit = "in")
#
# plot across tidal heights
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
  coord_cartesian(xlim = c(-.1,8),
                  ylim = c(-1,5),
                  expand = F,
                  clip = "off") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~group)

p2.2.2 + ggview::canvas(8,4)

#ggsave(p2.2.2, 
#       filename = "outputs/richness_change/p2.2.2_no_zeros.png",
#       width = 4,
#       height = 4,
#       unit = "in")