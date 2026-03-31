# Find abundance change based for counts data -- HIGHLY SAMPLED ONLY

# here, we repeat abundance change models including only data from
# highly sampled transects and tidal heights. The motivation for this is that
# occurrences that are uncommon within a transect will count that transect
# into the height-averaged count values, but there are many non-highly-sampled
# transects that stopped being sampled throughout the sampling duration.
# therefore, if a species occurred once or a few times in a given transect,
# zeros will be averaged in for as long as that transect exists in sampling,
# then zeros will stop being averaged in when the transect stops being sampled.
# This is likely to have an effect on abundance change, likely inflating the 
# number of increasing species (if lots of zeros are averaged into abundance 
# values only at the front half of the time series).


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(brms)
library(purrr)
library(marginaleffects)
library(tidyr)

# upload prepped data -------------------------------------------------------------
# this dataset is from 06.1-prep-count-data-for-model.R
# which summarizes count species between tidal heights in each year
data <- readr::read_csv(
  here::here("data-processed",
             "abundance-data",
             "count-data-prepped-for-model_HS.csv")
) %>% mutate(tidalheight = as.character(tidalheight))
data %>% glimpse()




# source modeling functions -----------------------------------------------
source(here::here(
  "scripts",
  "06-abundance-modeling",
  "06.2-create-model-functions.R"
))




# quick test to see if the one high urchin count changes the slope
urchin1 <- data %>%
  filter(organism == "Strongylocentrotus droebachiensis")
urchin1mod <- find_regression_slopes(urchin1)
urchin1slope <- extract_slope(urchin1mod)

urchin2 <-  data %>%
  filter(organism == "Strongylocentrotus droebachiensis") %>%
  filter(density_nonzero < 500)
urchin2mod <- find_regression_slopes(urchin2)
urchin2slope <- extract_slope(urchin2mod)


# get count slopes overall
moddf <- readRDS( here::here("data-processed",
                             "abundance-data",
                             "abundance_change_slopes_counts_cattidalheight_HS.rds"))
moddf %>% glimpse()
moddf %>% filter(!is.na(mean_monthly_mean)) %>% distinct(organism)




# model slopes with therm -------------------------------------------------
slopemod <- lm(data = moddf,
               formula = "slope ~ mean_monthly_mean")
summary(slopemod)
slopemod_au <- broom::augment(slopemod,
                              moddf,
                              interval = "prediction")

slopemod_au %>%
  ggplot(aes(x = mean_monthly_mean,
             y = slope)) +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "cyan4",
              alpha = .2) +
  geom_line(aes(y = .fitted)) +
  geom_point() +
  # ggrepel::geom_text_repel(aes(label = organism)) +
  geom_hline(yintercept = 0)


# see if removing the high urchin count changes anything ------------------
moddf %>%
  ggplot(aes(x = mean_monthly_mean, 
             y = slope)) +
  geom_point(aes(color = organism %in% c("Strongylocentrotus droebachiensis")) )



mod <- lm(data = moddf,
          formula = slope ~ mean_monthly_mean)

moddf_edited <-  moddf %>% 
  mutate(slope = case_when(
    organism == "Strongylocentrotus droebachiensis" ~ urchin2slope,
    TRUE ~ slope,
  ))

mod2 <- lm(data = moddf_edited,
           formula = slope ~ mean_monthly_mean)

summary(mod)
summary(mod2)


au1 <- broom::augment(mod, moddf, interval = "prediction")
au2 <- broom::augment(mod2, moddf_edited,
                      interval = "prediction")

au_merge <- au1 %>% mutate(group = "a) Including High Urchin Count") %>%
  rbind(au2 %>% mutate(group = "b) Excluding High Urchin Count"))

theme_set(ggthemes::theme_few())

mod1_slope <- summary(mod)$coefficients[2,1]
mod2_slope <- summary(mod2)$coefficients[2,1]
mod1_p <- summary(mod)$coefficients[2,4]
mod2_p <- summary(mod2)$coefficients[2,4]

mod_anns <- data.frame(
  group = c("a) Including High Urchin Count", "b) Excluding High Urchin Count"),
  slope = c(mod1_slope, mod2_slope),
  p = c(mod1_p, mod2_p)
)
mod_anns <- mod_anns %>%
  mutate(p_ann = ifelse(p < 0.01,"< 0.01", ifelse(p < 0.05,"< 0.05",round(p,2))))

p_urchin <- au_merge %>% 
  ggplot(aes(color)) +
  geom_hline(yintercept = 0) +
  
  # add regression for all points
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper,
                  x = mean_monthly_mean,
                  fill = group,
                  color = group),
              alpha = .4) +
  geom_line(aes(y = .fitted,
                x = mean_monthly_mean),
            linewidth = 1,
            alpha = .5) +

  # add points for all organism
  geom_point(aes(x = mean_monthly_mean, 
                 y = slope,
                 shape = organism == "Strongylocentrotus droebachiensis"),
             show.legend = F,
             size = 2,
             alpha = .7) +
  scale_shape_manual(values = c(1,19))+
  facet_wrap(~group) +
  
  geom_text(data = mod_anns, 
            aes(x = 9, 
                y = .2,
                label = paste0("Slope = ",round(slope,3),
                               "\n",
                               "p ", p_ann)),
            hjust = 0) +
  scale_x_continuous(labels = ~paste0(.,"°")) +
  
  labs(x = "Mean Thermal Affinity (°C)",
       y = "Coefficient of Abundance Change (Density Spp)") +
  ggrepel::geom_label_repel(
    data = au_merge %>% 
      filter(organism == "Strongylocentrotus droebachiensis" ),
  aes(x = mean_monthly_mean, 
      y = slope,
      label = stringr::str_wrap(organism,10)),
  lineheight = .8,
  nudge_y = .05)  +
  theme(legend.position = "none",
        strip.text = element_text(size = 12))


p_urchin + ggview::canvas(8,4)

ggsave(p_urchin,
       file = file.path("review_response","figures","urchin_plot.png"),
       width = 8, 
       height = 4)
