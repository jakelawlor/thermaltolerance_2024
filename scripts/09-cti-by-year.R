# find mean thermal affinity over time in levels

# here, I use the presence/absence dataset to find the mean
# thermal affinity of communities within tidal heights to 
# ask whether the communities are becoming more warm-affinity
# over time within the study region


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(lme4)


# data --------------------------------------------------------------------
pa <- readr::read_csv(
  here::here(
    "data-processed",
    "spp_pres_by_replicate.csv"
  )
)

therm <- readr::read_csv(
  here::here(
    "data-processed",
    "spp_thermal_affinities.csv"
  )
)

pa <- pa %>% left_join(therm) %>%
  filter(pres == T) %>%
  # translate level into tidal height
  mutate(tidalheight = (13-level)*.348)  %>%
  # filter out species that we don't have thermal affinity for
  filter(!is.na(mean_monthly_mean))
# 5 organisms total
pa %>% distinct(data_taken)


# find CTI per level - TOTAL ----------------------------------------------
# first, find the community thermal index (mean of STIs) for all 
# unique species in each level/year. This is a bit weird as it 
# compares different numbers of species, but as long as we assume
# that thermal index is not correlated with rarity / detectability
# of species, it should be fine. 
cti_all <- pa %>%
  # group by tidal height per year
  group_by(level, tidalheight, year) %>% 
  # get rid of all sampling variables except the 3 above
  select(-replicate,
         -data_taken,
         -pres_cover,
         -pres,
         -both_taken,
         -transect,
         -pres_count) %>% 
  # find distinct species
  distinct() %>%
  summarize(cti_mean = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_spp = n()) %>%
  ungroup() %>%
  # filter to levels where there are more than 4 species
  filter(n_spp >= 4) %>%
  # keep only levels that have data 3 or more times so that
  # a trend can be formed 
  group_by(level) %>% 
  add_count(name = "n_times_level_represented") %>% 
  filter(n_times_level_represented > 2) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level))
  

cti_all %>%
  ggplot(aes(x = year, y = cti_mean,
             #color = level, 
             group = tidalheight)) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight_f) +
  geom_point(alpha = .6) +
  gghighlight::gghighlight(
    unhighlighted_params = list(color = "grey80",
                                alpha = .3),
    use_direct_label = F) 


# make a few models -------------------------------------------------------
cti_all_add <- lm(data = cti_all,
                  formula = "cti_mean ~ year + tidalheight")
summary(cti_all_add) 
AIC(cti_all_add) #-218.227

cti_all_mult <- lm(data = cti_all,
                   formula = "cti_mean ~ year * tidalheight")

summary(cti_all_mult)
AIC(cti_all_mult) #-236.4831

cti_all_year <- lm(data = cti_all,
                   formula = "cti_mean ~ year")

summary(cti_all_year)
AIC(cti_all_year) #-206.3783

cti_all_year_wt <- lm(data = cti_all,
                   formula = "cti_mean ~ year",
                   weights = n_spp)

summary(cti_all_year)
AIC(cti_all_year) #-206.3783

cti_all_add_wt <- lm(data = cti_all,
                      formula = "cti_mean ~ year + tidalheight",
                      weights = n_spp)
summary(cti_all_year_wt)
AIC(cti_all_year_wt) #-248.7909
 
cti_all_mult_wt <- lm(data = cti_all,
                     formula = "cti_mean ~ year * tidalheight",
                     weights = n_spp)
summary(cti_all_mult_wt)
AIC(cti_all_mult_wt) #-265.6291

cti_all_par <- lm(data = cti_all,
                  formula = "cti_mean ~ year + I(tidalheight^2)")
summary(cti_all_par)
AIC(cti_all_par) #-215.9753

marginaleffects::avg_slopes(cti_all_mult_wt)
performance::compare_performance(
  cti_all_add,
  cti_all_add_wt,
  cti_all_mult,
  cti_all_mult_wt,
  cti_all_par,
  cti_all_year,
  cti_all_year_wt,
  rank = T
)
# cti_all_mult_wt is the best
plot(cti_all_mult_wt)

pred_df <- data.frame(
  tidalheight = rep(unique(cti_all$tidalheight), times = length(unique(cti_all$year))),
  year = rep(unique(cti_all$year), each = length(unique(cti_all$tidalheight)))
)

# see if I need to add the weight variable (n_spp) in the "newdf" for prediction
au <- broom::augment(cti_all_mult_wt,
                     newdata = pred_df,
                     interval = "prediction"
                     )

p <- au %>% 
  left_join(cti_all) %>%
  ggplot(aes(x = year, 
             y = .fitted,
             group = tidalheight,
             color = tidalheight)) +
 # geom_ribbon(aes(ymin = .lower,
 #                 ymax = .upper, fill = tidalheight),
 #             alpha = .1) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(y = cti_mean,
                 fill = tidalheight),
             alpha = .7,
             color = "black",
             position = position_jitter(width = .3,
                                        height = 0),
             shape = 21,
             size = 2,
             stroke = .2) +
  scale_color_viridis_c(option = "viridis") +
  scale_fill_viridis_c(option = "viridis") +
  labs(x = NULL,
       y = "Community Thermal Index\n(Averaged Thermal Affinity of All Species",
       color = "Tidal Height",
       fill = "Tidal Height") +
  theme(legend.position = "inside",
        legend.position.inside = c(.02, .98),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  
  guides(fill  = guide_colorbar(theme = theme(legend.title.position = "top")),
         color = guide_colorbar(theme = theme(legend.title.position = "top"))) +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(labels = ~paste0(.,"°"))

p + ggview::canvas(4.5,4)
ggsave(p,
       filename = "outputs/cti/cti_p1.png",
       width = 4.5, height = 4, unit = "in")

p2 <- sjPlot::plot_model(cti_all_mult_wt, type = "pred", 
                   terms = "year", show.data = T,
                   jitter = .02,
                   color = "cyan4", 
                   shape=16,
                   dot.size = 1.2,
                   dot.shape = 16,
                   line.size = 1.75,
                   dot.alpha = .4,
                   se = T,
                   alpha = .4) +
  labs(x = NULL,
       y = "Community Thermal Index\n(Averaged Thermal Affinity of All Species",
       title = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(labels = ~paste0(.,"°"))

p2 + ggview::canvas(4,4)

ggsave(p2,
       filename = "outputs/cti/cti_p2.png",
       width = 4.5, height = 4, unit = "in")

marginaleffects::avg_slopes(cti_all_mult_wt, variables = "year")
marginaleffects::avg_slopes(cti_all_mult_wt, variables = c("year","tidalheight"))
marginaleffects::slopes(cti_all_mult_wt) %>%
  as_tibble() %>%
  group_by(year,tidalheight) %>% distinct(estimate)
  group_by(estimate) %>% distinct(year,tidalheight)
  filter(estimate == max(estimate)) %>% glimpse()


au %>% 
  left_join(cti_all) %>%
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
    geom_path(arrow = arrow()) +
  geom_point(data = cti_all,
             aes(x = cti_mean,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  scale_color_viridis_c()+
  ggview::canvas(4,4)



# average CTI per tidalheight ---------------------------------------------------

cti_per_replicate <- 
  pa %>%
  group_by(year, transect, tidalheight, replicate) %>%
  summarize(cti_mean = mean(mean_monthly_mean, na.rm=T),
            spp_in_replicate = n()) %>%
  filter(spp_in_replicate > 4) %>%
  group_by(tidalheight) %>%
  add_count(name = "n_in_tidalheight") %>%
  mutate(n_years = length(unique(year))) %>%
  filter(n_years > 3)

cti_per_replicate %>%
  ggplot(aes(x = year, 
             y = cti_mean,
             group = tidalheight)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight) +
  gghighlight::gghighlight(
    use_direct_label = F,
    unhighlighted_params = list(color = "grey80",
                                alpha = .3)
  )

range(cti_per_replicate$year)
cti_per_replicate <- cti_per_replicate %>% mutate(year = year - 1981)
cti_rep_add <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year + tidalheight + (1|transect)")
summary(cti_rep_add)

cti_rep_add_wt <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year + tidalheight + (1|transect)",
                  weights = spp_in_replicate)
summary(cti_rep_add_wt)


cti_rep_mult <- lmer(data = cti_per_replicate,
                     formula = "cti_mean ~ year * tidalheight + (1|transect)")
summary(cti_rep_mult)
marginaleffects::avg_slopes(cti_rep_mult)

cti_rep_mult_wt <- lmer(data = cti_per_replicate,
                   formula = "cti_mean ~ year * tidalheight + (1|transect)",
                   weights = spp_in_replicate)
summary(cti_rep_mult_wt)
marginaleffects::avg_slopes(cti_rep_mult_wt)

cti_rep_year <- lmer(data = cti_per_replicate,
                   formula = "cti_mean ~ year + (1|transect)")
summary(cti_rep_year)


cti_rep_year_wt <- lmer(data = cti_per_replicate,
                   formula = "cti_mean ~ year + (1|transect)",
                   weights = spp_in_replicate)
summary(cti_rep_year_wt)


cti_rep_year_par <- lmer(data = cti_per_replicate,
                      formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)")
summary(cti_rep_year_par)
marginaleffects::avg_slopes(cti_rep_year_par)

cti_rep_par_wt <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)",
                  weights = spp_in_replicate)
summary(cti_rep_par_wt)
marginaleffects::avg_slopes(cti_rep_par_wt)



cti_rep_multpar <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)")
summary(cti_rep_multpar)
marginaleffects::avg_slopes(cti_rep_multpar)

cti_rep_multpar_wt <- lmer(data = cti_per_replicate,
                     formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)",
                     weights = spp_in_replicate)
summary(cti_rep_multpar_wt)
marginaleffects::avg_slopes(cti_rep_multpar_wt)



performance::compare_performance(
  cti_rep_add,
  cti_rep_add_wt,
  cti_rep_mult,
  cti_rep_mult_wt, 
  cti_rep_year, 
  cti_rep_year_wt,
  #cti_rep_par,
  cti_rep_par_wt,
  cti_rep_multpar,
  cti_rep_multpar_wt,
  rank = T
)


au_rep <- broom.mixed::augment(cti_rep_mult,
                         cti_per_replicate)
marginaleffects::avg_slopes(cti_rep_mult)

au_rep %>%
  ggplot(aes(x = year, 
             y = .fitted,
             group = interaction(tidalheight, transect),
             color = tidalheight)) +
  geom_point(aes(fill = tidalheight,
                 y = cti_mean),
             size = 2,
             alpha = .3,
             position = position_jitter(width = .3,
                                        height = 0)) +
  geom_line(linewidth = 1.5) +
  scale_fill_viridis_c() +
  scale_color_viridis_c()
  

# NOTE: Jarrett says try modelbased
# https://easystats.github.io/modelbased/
install.packages("modelbased")
library(modelbased)

estimate_relation(cti_rep_mult) |> plot()


library(ggeffects)

by_rep_year <- ggpredict(cti_rep_mult, "year",
                         interval = "prediction")

by_rep_year %>%
  ggplot(aes(x = x,
             y = predicted)) +
  geom_point(data = cti_per_replicate,
             aes(x = year, 
                 y = cti_mean,
                 fill = tidalheight),
             shape = 21,
             size = 2,
             alpha = .5,
             position = position_jitter(width = .3, 
                                        height = 0)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = .5) +
  scale_x_continuous(breaks = c(-1,9,19,29,39),
                     labels = function(x) x + 1981) +
  geom_line(linewidth = 1.5)  +
  scale_fill_viridis_c() +
  labs(x = NULL,
       y = "Community Thermal Index",
       fill = "Tidal Height") +
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(theme = theme(legend.title.position = "top")))

predict_response(cti_rep_mult,
        terms = "year")

# find per transect -------------------------------------------------------


cti_per_transect <- 
  pa %>%
  group_by(year, transect, tidalheight) %>%
  distinct(organism, mean_monthly_mean) %>%
  summarize(cti_mean = mean(mean_monthly_mean, na.rm=T),
            spp_in_transect = n()) %>%
  filter(spp_in_transect > 4) %>%
  group_by(tidalheight) %>%
  add_count(name = "n_in_tidalheight") %>%
  mutate(n_years = length(unique(year))) %>%
  filter(n_years > 3)

cti_per_transect %>%
  ggplot(aes(x = year, 
             y = cti_mean,
             group = tidalheight)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight) +
  gghighlight::gghighlight(
    use_direct_label = F,
    unhighlighted_params = list(color = "grey80",
                                alpha = .3)
  )



cti_per_transect <- cti_per_transect %>% mutate(year = year - 1981)
cti_tran_add <- lmer(data = cti_per_transect,
                  formula = "cti_mean ~ year + tidalheight + (1|transect)")
summary(cti_tran_add)

cti_tran_add_wt <- lmer(data = cti_per_transect,
                     formula = "cti_mean ~ year + tidalheight + (1|transect)",
                     weights = spp_in_transect)
summary(cti_tran_add_wt)


cti_tran_mult <- lmer(data = cti_per_transect,
                   formula = "cti_mean ~ year * tidalheight + (1|transect)")
summary(cti_tran_mult)
marginaleffects::avg_slopes(cti_tran_mult)

cti_tran_mult_wt <- lmer(data = cti_per_transect,
                      formula = "cti_mean ~ year * tidalheight + (1|transect) ",
                      weights = spp_in_transect)
summary(cti_tran_mult_wt)
marginaleffects::avg_slopes(cti_tran_mult_wt)

cti_tran_year <- lmer(data = cti_per_transect,
                   formula = "cti_mean ~ year + (1|transect)")
summary(cti_tran_year)


cti_tran_year_wt <- lmer(data = cti_per_transect,
                      formula = "cti_mean ~ year + (1|transect)",
                      weights = spp_in_transect)
summary(cti_tran_year_wt)


cti_tran_par <- lmer(data = cti_per_transect,
                  formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)")
summary(cti_tran_par)
marginaleffects::avg_slopes(cti_tran_par)

cti_tran_par_wt <- lmer(data = cti_per_transect,
                     formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)",
                     weights = spp_in_transect)
summary(cti_tran_par_wt)
marginaleffects::avg_slopes(cti_tran_par_wt)


cti_tran_multpar <- lmer(data = cti_per_transect,
                      formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)")
summary(cti_tran_multpar)
marginaleffects::avg_slopes(cti_tran_multpar)

cti_tran_multpar_wt <- lmer(data = cti_per_transect,
                         formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)",
                         weights = spp_in_transect)
summary(cti_tran_multpar_wt)
marginaleffects::avg_slopes(cti_tran_multpar_wt)

performance::compare_performance(
  cti_tran_add,
  cti_tran_add_wt,
  cti_tran_mult,
  cti_tran_mult_wt, 
  cti_tran_year, 
  cti_tran_year_wt,
  cti_tran_par,
  cti_tran_par_wt,
  cti_tran_multpar,
  cti_tran_multpar_wt,
  rank = T
)


by_transect <- ggpredict(cti_tran_mult, "year",
                         interval = "prediction")

by_transect %>%
  ggplot(aes(x = x,
             y = predicted)) +
  geom_point(data = cti_per_transect,
             aes(x = year, 
                 y = cti_mean,
                 fill = tidalheight),
             shape = 21,
             size = 2,
             alpha = .5,
             position = position_jitter(width = .3, 
                                        height = 0)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = .5) +
  scale_x_continuous(breaks = c(-1,9,19,29,39),
                     labels = function(x) x + 1981) +
  geom_line(linewidth = 1.5)  +
  scale_fill_viridis_c() +
  labs(x = NULL,
       y = "Community Thermal Index",
       fill = "Tidal Height") +
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.key.width = unit(dev.size()[1] / 10, "inches"),
        legend.key.height = unit(dev.size()[1] / 20, "inches"),
        legend.background = element_blank()) +
  guides(fill = guide_colorbar(theme = theme(legend.title.position = "top")))


