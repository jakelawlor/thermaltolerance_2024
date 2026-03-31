# model counts with tweedie distribution as sensitivity test

# HIGHLY SAMPLED ONLY

library(glmmTMB)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# upload prepped data -------------------------------------------------------------
# this dataset is from 06.1-prep-count-data-for-model.R
# which summarizes count species between tidal heights in each year
data <- readr::read_csv(
  here::here("data-processed",
             "abundance-data",
             "count-data-prepped-for-model_HS.csv")
) %>% mutate(tidalheight = as.character(tidalheight))
data %>% glimpse()


# split data --------------------------------------------------------------
# most models will be run with tidalheight as a random intercept,
# but a few species are found in only one tidal height, so errors
# are thrown, or species are found in 2 tidalheights at low densities,
# so models won't converge. 
# Separate dataset into species found in >2 tidal tidalheights or not.


# separate species only present in one or two levels
counts_mult <- data %>%
  group_by(organism,level) %>% nest() %>% 
  group_by(organism) %>% add_count() %>%
  filter(n>2) %>%
  select(-n) %>%
  unnest(data)
counts_mult %>% distinct(organism)
# 22 species

counts_sing <- data %>%
  group_by(organism,level) %>% nest() %>% 
  group_by(organism) %>% add_count() %>%
  filter(n<=2) %>%
  select(-n) %>%
  unnest(data)
counts_sing %>% distinct(organism)
# 0 species total was found in 2 or fewer levels

rm(data)

i <- 2
glmmTMB(data = counts_mult %>% filter(organism == unique(counts_mult$organism)[i]),
        formula = density ~ year_zero + tidalheight,
        family = tweedie(link = "log"))

# create tweedie modeling functions ---------------------------------------
find_regression_slopes_tweedie <- function(df){
  glmmTMB(density ~ year_zero + tidalheight,
          data=df, family=tweedie(link="log"))
}

extract_slope_tweedie <- function(model){
  fixef(model)$cond["year_zero"]
}

find_regression_slopes_few_tweedie <-function(df){
  glmmTMB(density ~ year_zero,
          data=df, family=tweedie(link="log"))
}



# apply functtions to multispecies data -----------------------------------


int_counts_models_tweedie <- counts_mult %>%
  
  group_by(organism) %>%
  nest() %>%
  
  mutate(
    # mutate a column with full model details
    model = map(.x = data,
                .f = find_regression_slopes_tweedie,
                .progress = T)
  )
# ok, here, models wouldn't converge for 1 species
# 1. "Anurida maritima"
# actually in all cases, the specific warning was "singular convergence (7)"
# looks like these errors are caused by tidal heights where there are only zeros
# let's just ignore for now... maybe later we can remove the levels with 
# only zeros, but the issue is that these are cases in which the zero-only 
# levels are between levels that have values. It makes more sense to leave them for now
# than to lump these into the single-level models.
int_counts_models_tweedie[int_counts_models_tweedie$organism == "Mytilus edulis", ]$data[[1]] %>% 
  ggplot(aes(x = year, y = density )) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tidalheight)
int_counts_models_tweedie[int_counts_models_tweedie$organism == "Mytilus edulis", ]$data[[1]] %>% 
  group_by(tidalheight) %>% 
  summarize(max = max(density))
# find coefficients

int_counts_coefs_tweedie <- int_counts_models_tweedie %>%
  
  mutate(
    # extract coef 
    slope_tweedie = map_dbl(model,extract_slope_tweedie)
  )
rm(int_counts_models_tweedie)






# merge -------------------------------------------------------------------
tweedie_merge <-
  int_counts_coefs_tweedie 
rm(int_counts_coefs_tweedie,
   counts_mult, counts_sing)


tweedie_merge <- tweedie_merge %>% select(-model)

# save --------------------------------------------------------------------
saveRDS(tweedie_merge,
        here::here(
          "data-processed",
          "abundance-data",
          "abundance_change_slopes_counts_tweedie_cattidalheight_HS.rds"
        ))


