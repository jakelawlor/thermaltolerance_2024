# functions for modeling abundance counts

library(brms)


# models with random intercept --------------------------------------------
# find density as a relationship between year and categorical tidalheight
find_regression_slopes <- function(df){
  brm(density_nonzero ~ year_zero + tidalheight,
      data=df, family=Gamma(link="log"),
      #prior=c(prior(normal(0,2),class="Intercept"),
      #        prior(normal(0,2),class="b"),
      #        prior(gamma(0.01,0.01),class="shape")),
      chains=4,iter=4000, cores=4, init = "0", inits=NULL,
      control = list(adapt_delta = 0.99)#,
      #backend = "cmdstanr"
      )
}


extract_slope <- function(model){
  fixef(model)["year_zero","Estimate"]
}

extract_slope_se <- function(model){
  fixef(model)["year_zero","Est.Error"]
}

extract_slope_q2.5 <- function(model){
  fixef(model)["year_zero","Q2.5"]
}

extract_slope_q97.5 <- function(model){
  fixef(model)["year_zero","Q97.5"]
}


extract_marginal_effects <- function(model){
  marg_effs <- avg_slopes(model,
                          variables = c("year_zero"))
  marg_effs %>%
    as.data.frame() %>%
    select(estimate,
           conf.low,
           conf.high) %>%
    return()
}


predict_model <- function(model){
  conditional_effects(model, 
                      effects = c("year_zero:tidalheight"),
                      resolution = length(unique(model$data$year_zero))#,
                      # here is where i would add tidal level if i want to! 
                      # as conditions_init
                      #int_conditions = list(tidalheight = c("4.524","2.436","0"))
  )[[1]][,c("tidalheight","estimate__","lower__","upper__","year_zero")]
}

# practice plotting model predictions
#predict_model(mod) %>% ggplot(aes(x = year_zero,
#                                  y = estimate__,
#                                  ymax = upper__,
#                                  ymin = lower__,
#                                  group = tidalheight)) +
#  geom_ribbon(alpha = .5) +
#  geom_line(linetype = "dashed", linewidth = .25) +
#  geom_point(data = testdf, 
#             inherit.aes=F,
#             aes(x=year_zero,
#                 y = density_nonzero),
#             alpha = .5,
#             size = .2)  + facet_wrap(~tidalheight)


# predict with tidybayes instead
# tidybayes::epred_draws(mod,
#                        newdata = data.frame(year_zero = rep(c(1:40), times = 4),
#                                             tidalheight = rep(c("-0.348","1.044","2.784","4.524"),
#                                                               each = length(c(1:40)))),
#                        ndraws = 100) %>%
#   ggplot(aes(x = year_zero,
#              y = .epred,
#              color = tidalheight,
#              group = paste0(tidalheight,.draw))) +
#   geom_line(linewidth = .25, alpha = .25)




# test on subset ----------------------------------------------------------

data <- readr::read_csv(
  here::here("data-processed",
             "abundance-data",
             "count-data-prepped-for-model.csv")
)

testdf <- data %>%
  filter(organism == "Cancer borealis") %>%
  mutate(tidalheight = as.character(tidalheight))

# 5.2 Test functions on subset --------------------------------------------
# 
# testmod <- find_regression_slopes(testdf)
# extract_slope(testmod)
# extract_slope_se(testmod) 
# extract_slope_q2.5(testmod) 
# extract_slope_q97.5(testmod) 
# #extract_intercept(testmod)
# extract_marginal_effects(testmod)
# summary(testmod)
# predict_model(testmod) 
# conditional_effects(testmod, effects = c("year_zero"))
# 
# # ok, functions all work, so we should be able to apply them to the full dataset




# models without intercept ------------------------------------------------

find_regression_slopes_few <- function(df){
  brm(density_nonzero ~ year_zero,
      data=df, family=Gamma(link="log"),
      #prior=c(prior(normal(0,2),class="Intercept"),
      #        prior(normal(0,2),class="b"),
      #        prior(gamma(0.01,0.01),class="shape")),
      chains=4,iter=4000, cores=4, init = "0", inits=NULL#,
      #backend = "cmdstanr"#
      )
}


extract_marginal_effects_few <- function(model){
  marg_effs <- avg_slopes(model,
                          variables = c("year_zero"))
  marg_effs %>%
    as.data.frame() %>%
    select(estimate,
           conf.low,
           conf.high) %>%
    return()
}


predict_model_few <- function(model){
  conditional_effects(model, 
                      effects = c("year_zero"),
                      resolution = length(unique(model$data$year_zero))
  )[[1]][,c("estimate__","lower__","upper__","year_zero")] %>%
    # add in tidal height as random factor
    mutate(tidalheight=as.factor("0"))
}

#testdf2 <- counts6 %>%
#  filter(organism == "Henricia sanguinolenta") %>%
#  filter(!is.na(density)) %>%
#  mutate(density_nonzero = ifelse(density == 0, .01, density))
#
#testmod2 <- find_regression_slopes_few(testdf2)
#extract_slope(testmod2)
#extract_slope_se(testmod2) 
#extract_slope_q2.5(testmod2) 
#extract_slope_q97.5(testmod2) 
##extract_intercept(testmod)
#extract_marginal_effects_few(testmod2)
#summary(testmod2)
#predict_model_few(testmod2) 


