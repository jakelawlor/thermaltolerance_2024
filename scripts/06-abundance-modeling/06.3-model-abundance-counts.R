# Find abundance change based for counts data


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
             "count-data-prepped-for-model.csv")
) %>% mutate(tidalheight = as.character(tidalheight))
data %>% glimpse()



# split data --------------------------------------------------------------
# most models will be run with tidalheight as a fixed categorical intercept,
# but a few species are found in only one tidal height, so errors
# are thrown, or species are found in 2 tidalheights at low densities,
# so models won't converge. 
# Separate dataset into species found in >2 tidal tidalheights or not.

## going to try to ignore this for now because I think models will 
# converge better now. 
# separate species only present in one or two levels
#counts_mult <- data %>%
#  group_by(organism,level) %>% nest() %>% 
#  group_by(organism) %>% add_count() %>%
#  filter(n>2) %>%
#  select(-n) %>%
#  unnest(data) %>%
#  # change tidal height to a categorical variable
#  mutate(tidalheight = as.character(tidalheight))
#counts_mult %>% distinct(organism)
## 25 species

#counts_sing <- data %>%
#  group_by(organism,level) %>% nest() %>% 
#  group_by(organism) %>% add_count() %>%
#  filter(n<=2) %>%
#  select(-n) %>%
#  unnest(data) %>%
#  # change tidal height to a categorical variable
#  mutate(tidalheight = as.character(tidalheight))
#counts_sing %>% distinct(organism)
## 1 species total was found in 2 or fewer levels



# make a df to test on ----------------------------------------------------
# i <- 2
# counts_mult %>%
#   filter(organism == unique(counts_mult$organism)[i]) %>%
#   ggplot(aes(x = year_zero, y = density_nonzero, group = tidalheight)) +
#   geom_point() +
#   facet_wrap(~tidalheight, scales = "free_y") +
#   geom_smooth(method = "lm", se = F)
# 
# testdf <- 
#   counts_mult %>%
#   filter(organism == unique(counts_mult$organism)[i])

# source modeling functions -----------------------------------------------
source(here::here(
  "scripts",
  "06-abundance-modeling",
  "06.2-create-model-functions.R"
))

testdf <- data %>%
  filter(organism == unique(data$organism)[1])

mod <- find_regression_slopes(testdf)
summary(mod)


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


# No longer need to do this with fixed effect model:
# 7. Repeat for level-constrained species ----------------------------------
# since 2 species were seen in 2 or fewer todai heights, we will model them separately
# without an effect of tidalheight becuase it wouldn't be doing anything.
#
## 7.2 apply new functions to level-constricted species --------------------
#counts_sing2 <- counts_sing #%>%
#
#  # add in the one that didn't converge above
# # rbind(counts_mult %>%
##          filter(organism == "Pagurus acadianus"))
## note, we didn't get model convergance issues with the fixed effect
## model, so no longer going to do this here
#
#
#int_counts_models_few <- counts_sing2 %>%
#  
#  group_by(organism) %>%
#  nest() %>%
#  
#  mutate(
#    # mutate a column with full model details
#    model = map(.x = data,
#                .f = find_regression_slopes_few)
#  )
#
#
#int_counts_coefs_few <- int_counts_models_few %>%
#  
#  mutate(
#    
#    # extract coef and intercept
#    slope = map_dbl(model,extract_slope),
#    slope_se = map_dbl(model,extract_slope_se),
#    slope_q2.5 = map_dbl(model,extract_slope_q2.5),
#    slope_q97.5 = map_dbl(model,extract_slope_q97.5),
#    #intercept = map_dbl(model,extract_intercept_few),
#    margeffs = map(model, extract_marginal_effects_few),
#    pred = map(model, predict_model_few)
#    
#  ) %>%
#  unnest(margeffs) %>%
#  arrange(slope)  %>%
#  ungroup() %>%
#  mutate(organism = forcats::fct_reorder(organism, slope))
#
#
#
## plot slopes and color those that the conf.int doesn't pass zero
#int_counts_coefs_few %>% glimpse()
#
#
#
#int_counts_coefs_few %>%
#  select(organism, slope, slope_se, estimate, slope_q2.5, slope_q97.5) %>%
#  mutate(organism = forcats::fct_reorder(organism, slope)) %>%
#  ggplot(aes(x=organism, y = slope)) +
#  geom_hline(yintercept = 0) +
#  # geom_segment(aes(xend = name, 
#  #                   y = slope_q2.5, yend = slope_q97.5)) +
#  geom_segment(aes(xend = organism, 
#                   y = slope - slope_se,
#                   yend = slope + slope_se,
#                   color = abs(slope) > abs(slope_se))) +
#  geom_point(aes(color = abs(slope) > abs(slope_se))) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#
#
#int_counts_coefs %>%
#  filter(organism =="Pagurus acadianus")
#

# 8. merge datasets -------------------------------------------------------
int_counts_coefs_full <- 
  int_counts_coefs #%>%
 # filter(organism != "Pagurus acadianus") %>%
 # rbind(int_counts_coefs_few  )
rm(int_counts_coefs, int_counts_coefs_few)

int_counts_coefs_full[1,]$data[[1]] %>% distinct(tidalheight)
int_counts_coefs_full[1,]$model[[1]] %>% summary()

# 9. plot FULL dataset ---glimpse()# 9. plot FULL dataset -------------------------------------
int_count_coefs_full_plotdf <- int_counts_coefs_full %>% 
  # put organisms in order by slope
  mutate(organism = forcats::fct_reorder(organism, slope)) %>%
  # unnest the predictions value
  unnest(pred) #%>%
#group_by(organism, tidalheight) %>% nest() %>% 
#group_by(organism) %>% arrange(desc(tidalheight)) %>%
#mutate(height_char = as.character(1:n())) %>%
#add_count() %>%
#ungroup() %>% 
#mutate(height_char = case_when(
#  height_char ==  "1" & n >= 3 ~ "High Intertidal",
#  height_char == "2" ~ "Mid Intertidal",
#  height_char == "3" ~ "Low Intertidal",
#  n < 3 ~ "Single Level")) %>%
#mutate(height_char = forcats::fct_reorder(height_char,tidalheight, .desc=T)) %>%

counts_points <- int_counts_coefs_full %>% 
  select(organism, data) %>%
  unnest(data)

# add in coordinate for the max point of each species
int_count_coefs_full_plotdf2 <- int_count_coefs_full_plotdf %>% 
  left_join(counts_points %>% 
              group_by(organism)%>% 
              summarize(max_y = max(density_nonzero)) %>%
              ungroup()) %>%
  group_by(organism) %>%
  mutate(max_y = pmax(max_y, max(upper__))) %>% ungroup()
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
              # previous annotation for percent change, but makes more sense as raw coefficient
              #mutate(perc_change = round((exp(slope)*100)-100),2),
              mutate(annotate = round(slope, 2)),
            aes(
              #label = paste0(perc_change,"%"),
              label = annotate, 
              x = 1+1981, y = max_y,
              hjust = -0.15, vjust = 1),
            size = 4.5,
            stat="unique") +
   geom_point(data = counts_points %>%
                mutate(organism = factor(organism, levels = levels(int_count_coefs_full_plotdf2$organism))),
              aes(y = density_nonzero),
              size = .3, alpha = .3) +
  geom_line(linewidth = .2) +
  facet_wrap(~organism, 
             scales = "free_y",
             labeller = label_wrap_gen(width=20))  +
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
 # scale_fill_manual(values = list("High Intertidal" = "#6ce1e6",
 #                                 "Mid Intertidal" = "#3499ad",
 #                                 "Low Intertidal" = "#156594",
 #                                 "Single Level" = "grey80")) +
  labs(y=NULL,
       x=NULL,
       fill = "Tidal Height",
       title = expression(paste("Density trends over time"))) +
  theme(axis.text.x = element_text(angle= 45, 
                                   hjust = 1, 
                                   vjust = 1),
        panel.spacing.x = unit(0.2,"lines"),
        panel.spacing.y = unit(.2, "lines"))

count_trends + ggview::canvas(13.5, 8)

# save
ggsave(count_trends,
       filename = "outputs/abundance-change/all_trends_counts.png",
       width = 13.5,
       height = 8,
       dpi = 300)



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



# save --------------------------------------------------------------------
saveRDS(moddf,
        here::here("data-processed",
                   "abundance-data",
                   "abundance_change_slopes_counts_cattidalheight.rds"))

moddf <- readRDS( here::here("data-processed",
                             "abundance-data",
                             "abundance_change_slopes_counts_cattidalheight.rds"))
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


rm(list = ls())
gc()
