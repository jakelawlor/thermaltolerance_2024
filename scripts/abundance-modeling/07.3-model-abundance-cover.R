# Find abundance change based for cover data



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(modelsummary)
library(ordbetareg)
library(tidyr)
library(purrr)
library(marginaleffects)


# upload prepped data -------------------------------------------------------------
# this dataset is from 06.1-prep-count-data-for-model.R
# which summarizes count species between tidal heights in each year
data <- readr::read_csv(
  here::here("data-processed",
             "abundance-data",
             "cover-data-prepped-for-model.csv")
) %>% mutate(tidalheight = as.character(tidalheight))
data %>% glimpse()

# make test set
i <- 4
testdf <- data %>%
  filter(organism == unique(data$organism)[i])
testdf %>% ggplot(aes(x = year_zero,
                      y = pc_num)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tidalheight)



# source modelling functions -----------------------------------------------
# source modeling functions -----------------------------------------------
source(here::here(
  "scripts",
  "abundance-modeling",
  "07.2-create-model-functions-cover.R"
))






# 6. Apply functions to full dataset -------------------------------------


# map a column to add a separate regression model to each species
int_cover_models <- data %>%
  
  group_by(organism) %>%
  nest() %>%
  
  mutate(intercept_prior = map(.x = data, .f=make_prior)) %>%
  
  mutate(
    # mutate a column with full model details
    model = map2(.x = data, 
                 .y = intercept_prior,
                 ~ find_regression_slopes(df = .x,
                                          prior = .y),
                 .progress = T)
  )
# NOTE: this threw warnings for 3 species with 
# divergent transitions when we used the random effect mode, 
# 1. Bonnemaisonia hamifera
# 2. Colpomenia peregrina
# 3. Palmaria palmata
# but that was not the case using this fixed effect model, 
# so we'll skip this step. 


# extract coefficients
int_cover_coeffs <- int_cover_models %>%
  
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


#
#redo_test <- data %>%
#  filter(organism %in% c("Bonnemaisonia hamifera",
#                         "Colpomenia peregrina",
#                         "Palmaria palmata")) %>%
#  group_by(organism) %>%
#  nest() %>%
#  mutate(intercept_prior = map(.x = data, .f=make_prior)) %>%
#  
#  mutate(
#    # mutate a column with full model details
#    model = map2(.x = data, 
#                 .y = intercept_prior,
#                 ~ find_regression_slopes_nolevel(df = .x,
#                                                  prior = .y),
#                 .progress = T)
#  )
#
## extract coefficients
#redo_coeffs <- redo_test %>%
#  
#  mutate(
#    
#    # extract coef and intercept
#    slope = map_dbl(model,extract_slope),
#    slope_se = map_dbl(model,extract_slope_se),
#    slope_q2.5 = map_dbl(model,extract_slope_q2.5),
#    slope_q97.5 = map_dbl(model,extract_slope_q97.5),
#    #intercept = map_dbl(model,extract_intercept),
#    margeffs = map(model, extract_marginal_effects),
#    pred = map(model, predict_model_nolevel)
#    
#  ) %>%
#  unnest(margeffs) %>%
#  arrange(slope)  %>%
#  ungroup() %>%
#  mutate(organism = forcats::fct_reorder(organism, slope))
#
# merge these redos into the original coefficients dataset
all_coeffs <- int_cover_coeffs #%>%
  #filter(!organism %in% c("Bonnemaisonia hamifera",
  #                        "Colpomenia peregrina",
  #                        "Palmaria palmata")) %>%
  #rbind(redo_coeffs)
rm(int_cover_coeffs, redo_coeffs)


# plot slopes and color those that the conf.int doesn't pass zero
all_coeffs %>% glimpse()

all_coeffs %>%
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



# 9. plot FULL dataset -------------------------------------
int_cover_coefs_plotdf <- all_coeffs %>% 
  mutate(organism = forcats::fct_reorder(organism, slope)) %>%
  unnest(pred)

# extract raw data points
cover_points <- all_coeffs %>% 
  select(organism, data) %>%
  unnest(data)

# add in coordinate for the max point of each species, 
# or the max upper limit, whichever is higher
int_cover_coefs_plotdf2 <- int_cover_coefs_plotdf %>% 
  left_join(cover_points %>% 
              group_by(organism)%>% 
              summarize(max_y = max(percent_cover_beta)) %>%
              ungroup()) %>%
  group_by(organism) %>%
  mutate(max_y = pmax(max_y, max(upper__))) %>% ungroup()
# make the max y coordinate whichever is bigger of the
# max point, or the max upper bound

cover_trends <- int_cover_coefs_plotdf2 %>%
  ggplot(aes(y = estimate__, 
             x = year_zero+1981,
             group = tidalheight
             )) +
  geom_ribbon(aes(ymin=lower__,
                  ymax=upper__),
              fill = "cyan3",
              alpha=0.4) +
  geom_hline(yintercept = 0)+
  geom_text(data = . %>% 
              group_by(organism) %>% arrange(desc(upper__)) %>% slice(1),
            aes(label = round(slope,2),
                #label = paste0(perc_change,"%"),
                x = 1+1981, y = max_y,
                hjust = -0.15, vjust = 1),
            size = 4.5,
            stat="unique") +
  geom_point(data = cover_points,
             aes(y = percent_cover_beta),
             size = .3, alpha = .3) +
  geom_line(linewidth = .1) +
  facet_wrap(~organism, 
             scales = "free_y",
             labeller = label_wrap_gen(width=20),
             ncol = 6)  +
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
  scale_y_continuous(labels = scales::percent_format()) +
  guides(fill=guide_legend(title.position="top", 
                           title.hjust =0)) +
 # scale_fill_manual(values = list("High Intertidal" = "#6ce1e6",
 #                                 "Mid Intertidal" = "#3499ad",
 #                                 "Low Intertidal" = "#156594",
 #                                 "Single Level" = "grey80")) +
  labs(y=NULL,
       x=NULL,
       fill = "Tidal Height",
       title = expression(paste("Percent cover trends over time"))) +
  theme(axis.text.x = element_text(angle= 45, 
                                   hjust = 1, 
                                   vjust = 1),
        panel.spacing.x = unit(0.2,"lines"),
        panel.spacing.y = unit(.2, "lines"))


cover_trends + ggview::canvas( width = 13.5, height = 14)
ggsave(cover_trends,
       filename = "outputs/abundance-change/all_trends_cover.png",
       width = 13.5,
       height = 14,
       dpi = 300)


# merge with stis ---------------------------------------------------------
therm <- readr::read_csv(
  here::here(
    "data-processed",
    "spp_thermal_affinities.csv"
  )
)


moddf <- all_coeffs %>% 
  left_join(therm) 

moddf %>% glimpse()
# save
saveRDS(moddf,
        here::here("data-processed",
                   "abundance-data",
                   "abundance_change_slopes_cover_cattidalheight.rds"))

moddf <- readRDS(
  here::here("data-processed",
             "abundance-data",
             "abundance_change_slopes_cover_cattidalheight.rds")
)

mod <- lm(data = moddf,
               formula = "slope ~ mean_monthly_mean")
summary(mod)
# NOTE: 48 species in the regression because Clathromo doesn't have therm

mod_au <- broom::augment(mod,
                         moddf %>% filter(!is.na(mean_monthly_mean)),
                         interval = "prediction")

mod_au %>%
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

# -------------------------------------------------------------------------
# scratch
# -------------------------------------------------------------------------



moddf %>% filter(is.na(mean_monthly_mean))
au_orig <- broom::augment(mod_orig,
                          interval = "prediction")

au_orig %>% 
  ggplot(aes(x = mean_monthly_mean,
             y = .fitted)) +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "red",
              alpha = .2) +
  geom_line() +
  geom_point(data = compare,
             aes(y = slope_orig))

mod_rand <- lm(data = compare,
               formula = "slope ~ mean_monthly_mean")
summary(mod_rand)
au_rand <- broom::augment(mod_rand,
                          interval = "prediction")

redo_coeffs <- redo_coeffs %>% left_join(therm)
au_rand %>% 
  ggplot(aes(x = mean_monthly_mean,
             y = .fitted)) +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              fill = "orange",
              alpha = .2) +
  geom_line() +
  geom_point(data = compare,
             aes(y = slope)) +
  geom_point(data = compare %>%
               filter(organism %in% c("Bonnemaisonia hamifera",
                                      "Colpomenia peregrina",
                                      "Palmaria palmata")),
             aes(y = slope),
             color = "red",
             shape = 21,
             size = 2,
             fill = "transparent") +
  geom_point(data = redo_coeffs,
             aes(y = slope),
             color = "purple")








# model slopes with therm -------------------------------------------------


slopemod <- lm(data = moddf,
               formula = "slope ~ mean_monthly_mean")
summary(slopemod)
slopemod_au <- broom::augment(slopemod,
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


moddf %>% filter(mean_monthly_mean < 10) %>% distinct(organism, slope, mean_monthly_mean)







  