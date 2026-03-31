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


# test worlflow on one species:
testdf <- data %>%
  filter(organism == unique(data$organism)[1])

mod <- find_regression_slopes(testdf)
summary(mod)
extract_slope(mod)



# model multi-tidal-height species ----------------------------------------
int_counts_models_fixed <- data %>%
  
  group_by(organism) %>%
  nest() %>%
  
  #  ungroup() %>%
  #  slice(1:3) %>%
  
  
  mutate(
    # mutate a column with full model details
    model = map(.x = data,
                .f = find_regression_slopes)
  )
# NOTE: one warning occurred for Pagurus acadianus


int_counts_models_fixed %>% distinct(organism) %>% pull
# models for 25 species

int_counts_coefs <- int_counts_models_fixed %>%
  
  mutate(
    
    # extract coef and intercept
    slope = map_dbl(model,extract_slope),
    slope_se = map_dbl(model,extract_slope_se),
    slope_q2.5 = map_dbl(model,extract_slope_q2.5),
    slope_q97.5 = map_dbl(model,extract_slope_q97.5),
    #intercept = map_dbl(model,extract_intercept),
    margeffs = map(model, extract_marginal_effects),
    pred = map(model, predict_model)
    
  ) %>%
  unnest(margeffs) %>%
  arrange(slope)  %>%
  ungroup() %>%
  mutate(organism = forcats::fct_reorder(organism, slope))




# plot slopes and color those that the conf.int doesn't pass zero
int_counts_coefs %>% glimpse()


int_counts_coefs %>%
  select(organism, slope, slope_se, estimate, slope_q2.5, slope_q97.5) %>%
  # mutate(organism = forcats::fct_reorder(name, slope)) %>%
  ggplot(aes(x=organism, y = slope)) +
  geom_hline(yintercept = 0) +
  # geom_segment(aes(xend = name, 
  #                   y = slope_q2.5, yend = slope_q97.5)) +
  geom_segment(aes(xend = organism, 
                   y = slope - slope_se,
                   yend = slope + slope_se,
                   color = abs(slope) > abs(slope_se))) +
  geom_point(aes(color = abs(slope) > abs(slope_se))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



# 8. merge datasets -------------------------------------------------------
int_counts_coefs_full <- 
  int_counts_coefs
rm(int_counts_coefs, int_counts_coefs_few)

int_counts_coefs_full[1,]$data[[1]] %>% distinct(tidalheight)
int_counts_coefs_full[1,]$model[[1]] %>% summary()


# merge with stis ---------------------------------------------------------
therm <- readr::read_csv(
  here::here(
    "data-processed",
    "species-thermal-affinities",
    "spp_thermal_affinities.csv"
  )
)


moddf <- int_counts_coefs_full %>% 
  left_join(therm)
moddf %>% distinct(organism)

# remove model from dataframe because they are very heavy
moddf <- moddf %>% select(-model)

# save --------------------------------------------------------------------
saveRDS(moddf,
        here::here("data-processed",
                   "abundance-data",
                   "abundance_change_slopes_counts_cattidalheight_HS.rds"))

rm(list = ls())

moddf <- readRDS( here::here("data-processed",
                             "abundance-data",
                             "abundance_change_slopes_counts_cattidalheight_HS.rds"))
moddf %>% glimpse()
moddf %>% filter(!is.na(mean_monthly_mean)) %>% distinct(organism)



# plot all for species ---------------------------------------------------------
moddf %>% glimpse()


int_count_coefs_full_plotdf <- moddf %>% 
  # put organisms in order by slope
  mutate(organism = forcats::fct_reorder(organism, slope)) %>%
  select(organism, slope, pred) %>%
  # unnest the predictions value
  unnest(pred) 

counts_points <- moddf %>% 
  select(organism, data) %>%
  unnest(data) %>%
  mutate(organism = factor(organism, 
                            levels = levels(int_count_coefs_full_plotdf$organism)))

# add in coordinate for the max point of each species
int_count_coefs_full_plotdf2 <- int_count_coefs_full_plotdf %>% 
  left_join(counts_points %>% 
              group_by(organism)%>% 
              summarize(max_y = max(density_nonzero)) %>%
              ungroup()) %>%
  group_by(organism) %>%
  mutate(max_y = pmax(max_y, max(upper__))) %>% ungroup() %>%
  mutate(organism = factor(organism, 
                           levels = levels(int_count_coefs_full_plotdf$organism)))

# make the max y coordinate whichever is bigger of the
# max point, or the max upper bound

count_trends <- int_count_coefs_full_plotdf2 %>%
  ggplot(aes(y = estimate__, x = year_zero+1981,
             group = tidalheight)) +
  geom_ribbon(aes(ymin=lower__,
                  ymax=upper__,
                  #fill = height_char
  ),
  fill = "cyan3",
  alpha=0.4) +
  geom_hline(yintercept = 0)+
  geom_text(data = . %>% 
              group_by(organism) %>% arrange(desc(upper__)) %>% slice(1) %>%
              mutate(annotate = round(slope, 3)),
            aes(
              label = annotate, 
              x = 1+1981, y = max_y,
              hjust = -0.15, vjust = 1),
            size = 4,
            stat="unique") +
  geom_point(data = counts_points %>%
               mutate(organism = factor(organism, levels = levels(int_count_coefs_full_plotdf2$organism))),
             aes(y = density_nonzero),
             size = .3, alpha = .3) +
  geom_line(linewidth = .2) +
  facet_wrap(~organism, 
             scales = "free_y",
             labeller = label_wrap_gen(width=20),
             nrow = 6)  +
  ggthemes::theme_few(base_size = 14) +
  theme(strip.text.x = element_text(size=13,
                                    face = "italic",
                                    margin = margin(b=0,t=1),
                                    lineheight=.85),
        panel.border = element_rect(color = "transparent",fill = "transparent"),
        panel.spacing.y = unit(10,"pt"),
        panel.spacing.x = unit(20,"pt"),
        plot.title = element_text(size=24),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = c(.5,1),
        legend.title = element_text(size=16),
        legend.text = element_text(size=13),
        legend.background  = element_rect(color = "black", size=.5),
        strip.clip = "off"
  ) +
  guides(fill=guide_legend(title.position="top", 
                           title.hjust =0)) +
  labs(y=NULL,
       x=NULL,
       fill = "Tidal Height",
       title = expression(paste("Density trends over time"))) +
  theme(axis.text.x = element_text(angle= 45, 
                                   hjust = 1, 
                                   vjust = 1),
        panel.spacing.x = unit(0.2,"lines"),
        panel.spacing.y = unit(.2, "lines"))

count_trends + ggview::canvas(10, 10)

# save
ggsave(count_trends,
       filename = "outputs/abundance-change/all_trends_counts_HS.png",
       width = 10,
       height = 10,
       dpi = 300)





# model slopes with therm -------------------------------------------------
slopemod <- lm(data = moddf[1:21,],
               formula = "slope ~ mean_monthly_mean")
summary(slopemod)
plot(slopemod)
slopemod_au <- broom::augment(slopemod,
                              moddf[1:21,],
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


rm(list = ls())
gc()


