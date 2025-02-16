# model counts with tweedie distribution as sensitivity test

library(glmmTMB)
library(dplyr)
library(tidyr)
library(purrr)

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
# 25 species

counts_sing <- data %>%
  group_by(organism,level) %>% nest() %>% 
  group_by(organism) %>% add_count() %>%
  filter(n<=2) %>%
  select(-n) %>%
  unnest(data)
counts_sing %>% distinct(organism)
# 1 species total was found in 2 or fewer levels

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
# ok, here, models wouldn't converge for 3 species
# 1. "Anurida maritima"
# 2. "Asterias forbesi"
# 3. "Mytilus edulis"
# actually in all cases, the specific warning was "singular convergence (7)"
# looks like these errors are caused by tidal heights where there are only zeros
# let's just ignore for now... maybe later we can remove the levels with 
# only zeros, but the issue is that these are cases in which the zero-only 
# levels are between levels that have values. It makes more sense to leave them for now
# than to lump these into the single-level models.
int_counts_models_tweedie[int_counts_models_tweedie$organism =="Mytilus edulis", ]$model$convergence
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


# apply to single species -------------------------------------------------

int_counts_models_few_tweedie <- counts_sing %>%
  
  group_by(organism) %>%
  nest() %>%
  
  mutate(
    # mutate a column with full model details
    model = map(.x = data,
                .f = find_regression_slopes_few_tweedie)
  )


int_counts_coefs_few_tweedie <- int_counts_models_few_tweedie %>%
  
  mutate(
    
    # extract coef and intercept
    slope_tweedie = map_dbl(model,extract_slope_tweedie)
  )
rm(int_counts_models_few_tweedie)




# merge -------------------------------------------------------------------
tweedie_merge <-
  int_counts_coefs_tweedie %>%
  rbind(int_counts_coefs_few_tweedie)
rm(int_counts_coefs_tweedie,
   int_counts_coefs_few_tweedie,
   counts_mult, counts_sing)



# save --------------------------------------------------------------------
saveRDS(tweedie_merge,
        here::here(
          "data-processed",
          "abundance-data",
          "abundance_change_slopes_counts_tweedie_cattidalheight.rds"
        ))



# -------------------------------------------------------------------------
# scratch ---------------------------------------------------------------
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Repeat with ZIGamma ------------------------------------------------------------------
# -------------------------------------------------------------------------
library(glmmTMB)
library(brms)

i <- 8
testorg <- unique(counts6$organism)[i]
testdf <- counts6 %>%
  filter(organism == testorg) %>%
  filter(!is.na(density)) %>%
  mutate(density_nonzero = ifelse(density == 0, .01, density))

testdf %>%
  #group_by(level) %>% add_count() %>% 
  # filter(n>3) %>% filter(level <= 13) %>%
  mutate(height_label = paste(round(tidalheight,2), "m")) %>%
  mutate(height_label = forcats::fct_reorder(height_label, level)) %>%
  ggplot(aes(x=year, y=density)) +
  geom_hline(yintercept = 0)+
  geom_point(size=.75, alpha=.5) +
  geom_smooth(method = "lm", color = "black", fill = "#3499ad", size=.75) +
  theme(panel.border = element_blank()) +
  facet_wrap(~height_label) +
  labs(x=NULL,
       y=NULL,
       title = paste(testorg,"density/m^2 change over time")) +
  theme(plot.title = element_text(size=18),
        plot.title.position = "plot") 


gamma_mod <- brm(density_nonzero ~ year_zero + (1|tidalheight),
                 data=testdf, family=Gamma(link="log"),
                 chains=4,iter=4000, cores=4, init = "0", inits=NULL,
                 control = list(adapt_delta = 0.99) )
summary(gamma_mod)
fixef(gamma_mod)
conditional_effects(gamma_mod)

zi_gamma_mod <- glmmTMB(density ~ year_zero + (1|tidalheight),
                        ziformula = ~.,
                        data = testdf, family = ziGamma(link = "log"))
summary(zi_gamma_mod)
fixef(zi_gamma_mod)


tweedie_mod <- glmmTMB(density ~ year_zero + (1|tidalheight),
                       data = testdf, family = tweedie(link = "log"))
summary(tweedie_mod)

#pred <- data.frame(year_zero = rep(min(testdf$year_zero):max(testdf$year_zero),
#                                   each = length(unique(testdf$tidalheight))),
#                   tidalheight = rep(unique(testdf$tidalheight),times = length(min(testdf$year_zero):max(testdf$year_zero))))
pred <- data.frame(year_zero = c(min(testdf$year_zero):max(testdf$year_zero)))
pred$predtweedie <- predict(tweedie_mod,
                            pred,
                            re.form = NA,
                            type=  "response",
                            se.fit  =T)$fit
pred$predgamma <- predict(gamma_mod,
                          pred,
                          type = "response",
                          re_formula = NA)[,"Estimate"]
pred$predzi <- predict(zi_gamma_mod,
                       pred, type = "response",
                       re.form = NA)
pred %>% glimpse()
ggplot(pred,
       aes(x=year_zero)) +
  geom_line(aes(y = predtweedie),
            color = "blue") +
  geom_line(aes(y = predgamma),
            color = "green") +
  geom_line(aes(y = predzi),
            color = "purple")



find_regression_slopes <- function(df){
  brm(density_nonzero ~ year_zero + (1|tidalheight),
      data=df, family=Gamma(link="log"),
      #prior=c(prior(normal(0,2),class="Intercept"),
      #        prior(normal(0,2),class="b"),
      #        prior(gamma(0.01,0.01),class="shape")),
      chains=4,iter=4000, cores=4, init = "0", inits=NULL,
      control = list(adapt_delta = 0.99) )
}




extract_slope <- function(model){
  fixef(model)["year_zero","Estimate"]
}


remove.packages(c("rstan", "StanHeaders", "brms"))
install.packages("rstan", type = "source")
install.packages("StanHeaders", type = "source")
install.packages("brms", type = "source")



