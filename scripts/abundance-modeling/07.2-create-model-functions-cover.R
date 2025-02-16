# create functions for modeling cover species with ordbetareg

library(ordbetareg)
library(modelsummary)
library(dplyr)

# NOTE, changing tidal height here from random effect to fixed categorical
# effect to give them fixed random intercepts that don't matter which is
# the highest or lowest

# Create Functions to Map on data -------------------------------------------------
make_prior <- function(df){
  mean_transformed <- c(brms::logit_scaled(mean(df$percent_cover_beta)))
  intercept_prior <- tapply(mean_transformed,1,
                            FUN = function(x) paste0("prior(normal(",x,", 1), coef= 'Intercept')"))
  prior2 <- eval(parse(text = intercept_prior))
  return(prior2)
}

find_regression_slopes <- function(df,prior){
  ordbetareg(percent_cover_beta ~ year_zero + tidalheight,
             data=df, 
             true_bounds = c(0,1),
             chains=4,iter=4000,refresh=0,
             cores = 4,
             control = list(adapt_delta = 0.99),
            # backend = "cmdstanr"
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
                          newdata = model$data,
                          variables = c("year_zero"),
                         # allow_new_levels = T,
                          type = "response"
                       )
  marg_effs %>%
    as.data.frame() %>%
    select(term,
           estimate,
           conf.low,
           conf.high) %>%
    return()
}


predict_model <- function(model){
  conditional_effects(model, 
                      effects = c("year_zero:tidalheight"),
                      resolution = length(unique(model$data$year_zero))
  )[[1]][,c("tidalheight","estimate__","lower__","upper__","year_zero")]
}




# make versions with no tidal height factor -------------------------------

# manually redo the 3 species that didn't converge ------------------------
find_regression_slopes_nolevel <- function(df,prior){
  ordbetareg(percent_cover_beta ~ year_zero,
             data=df, 
             true_bounds = c(0,1),
             chains=4,iter=4000,refresh=0,
             cores = 4,
             control = list(adapt_delta = 0.99)#,
             #backend = "cdmstanr"
  )
}
predict_model_nolevel <- function(model){
  conditional_effects(model, 
                      effects = c("year_zero"),
                      resolution = length(unique(model$data$year_zero))
  )[[1]][,c("estimate__","lower__","upper__","year_zero")]
}


# # # test the functions on a subset ------------------------------------------
# # 
#  data <- readr::read_csv(
#    here::here("data-processed",
#               "abundance-data",
#               "cover-data-prepped-for-model.csv")
#  )
#  
#  testdata <- data %>%
#    filter(organism == "Chondrus crispus")
# # 
# # 
# intercept_prior <- make_prior(testdata)
# testmod <- find_regression_slopes(testdata,intercept_prior)
# extract_slope(testmod)
# extract_slope_se(testmod) 
# extract_slope_q2.5(testmod) 
# extract_slope_q97.5(testmod) 
# #extract_intercept(testmod)
# extract_marginal_effects(testmod)
# summary(testmod)
# predict_model(testmod) 
# conditional_effects(testmod, effects = c("year_zero"))
# # 
# # # ok, functions all work, so we should be able to apply them to the full dataset
# # 
# 


# avg_slopes(testmod)
# avg_slopes(testmod,
#        newdata = testmod$data)
# -------------------------------------------------------------------------

# 
# 
# # 4. Test moel on subset ----------------------------------------------------------
# mean_transformed <- c(logit_scaled(mean(testdf$percent_cover_beta)))
# intercept_prior <- prior(normal(paste0(mean_transformed), 1), coef= "Intercept")
# intercept_prior <- tapply(mean_transformed,1,
#                           FUN = function(x) paste0("prior(normal(",x,", 1), coef= 'Intercept')"))
# prior2 <- eval(parse(text = intercept_prior))
# 
# ### NOTE: ordbetareg() is giving us warnings because it is using the function
# # "inits" in brms under the hood instead of "init" - this shouldn't affect results, 
# # but does prevent us from seeing other potential errors, so we're going to use 
# # "trace()" to edit each instance of "inits" to "init" in the function so that 
# # we no longer receive this error
# 
# 
# #trace(ordbetareg, edit=TRUE) # change 6 occurrences in the popup window
# range(testdf$year_zero)
# 
# testmod <- ordbetareg(percent_cover_beta ~ year_zero + (1|tidalheight),
#                       data=testdata, 
#                       true_bounds = c(0,1),
#                       chains=4,iter=4000,refresh=0,
#                       control = list(adapt_delta = 0.99)
#                       #coef_prior_mean = 0,
#                       #coef_prior_sd = 1,
#                       #extra_prior = prior2,
#                       # init = "0",
#                       #  inits = NULL
# )
# ?ordbetareg()
# 
# testmod$prior
# summary(testmod)
# 
# all_draws <- prepare_predictions(testmod)
# cutzero <- plogis(all_draws$dpars$cutzero)
# cutone <- plogis(all_draws$dpars$cutzero + exp(all_draws$dpars$cutone))
# testdata %>%
#   ggplot(aes(x=percent_cover_beta)) +
#   geom_histogram(bins = 100) +
#   geom_vline(xintercept = mean(cutzero), linetype=2)+
#   geom_vline(xintercept = mean(cutone), linetype=2)
# 
# 
# pp_check_ordbeta(testmod)$discrete
# pp_check_ordbeta(testmod)$continuous
# 
# # view the model coefficients for effect of year and range
# modelsummary(testmod,
#              statistic = "conf.int",
#              metrics = "RMSE",
#              coef_map=c("b_year_zero" = "Year",
#                         "b_tidalheight" = "Tidal Height"
#              ))
# 
# # see marginal effect of year (e.g. actual value that percent cover changes in each year)
# marg_effs <- avg_slopes(testmod,
#                         variables = c("year_zero"),
#                         type = "link")
# marg_effs %>%
#   as.data.frame() %>%
#   select(term,
#          estimate,
#          conf.low,
#          conf.high) %>%
#   return()
# # so in this case, we see that that percent cover decreases by about .4% per year
# # across tidal heights, and that percent cover decreases by 14% per meter up the shoreline
# 
# ord_pred <- conditional_effects(testmod,
#                                 effects = c("year_zero"))[[1]]
# 
# range(ord_pred$year_zero)
# 
# 
# example <- ord_pred %>% 
#   group_by(tidalheight) %>% nest() %>% ungroup() %>% arrange(tidalheight) %>%
#   mutate(height_char = as.character(1:n())) %>% unnest(data) %>%
#   mutate(height_char = dplyr::recode(height_char,
#                                      "1" = "Low Intertidal",
#                                      "2" = "Mid Intertidal",
#                                      "3" = "High Intertidal")) %>%
#   mutate(height_char = forcats::fct_reorder(height_char,-tidalheight)) %>%
#   ggplot(aes(y=estimate__*100,x=year_zero+1981,group = height_char)) +
#   geom_hline(yintercept = 0) +
#   geom_ribbon(aes(ymin=lower__*100,
#                   ymax=upper__*100,
#                   fill = height_char),#fill="#3499ad",
#               alpha=0.5) +
#   scale_fill_manual(values = list("High Intertidal" = "#6ce1e6",
#                                   "Mid Intertidal" = "#3499ad",
#                                   "Low Intertidal" = "#156594")) +
#   geom_point(data = testdf,
#              inherit.aes = F,
#              aes(x=year_zero+1981, y=pc_num),
#              alpha = .4) +
#   # geom_hline(yintercept=1,linetype=2) +
#   geom_line() +
#   scale_y_continuous(labels = scales::percent_format(scale=1)) +
#   labs(y="Percent Cover",
#        x=NULL,
#        fill = NULL,
#        #title = glue::glue(unique(testdf$name)," Abundance Change")
#   )+
#   theme(panel.border = element_blank(),
#         legend.position = c(1,1),
#         legend.justification = c(1,1))+
#   coord_cartesian(ylim = c(0,100),
#                   #xlim = c(1 + 1981,42 + 1981),
#                   expand=F,
#                   clip="off")
# 
# example
# 
# ggview::ggview(example, width = 5, height = 5, units = "in")
# 
# 
# 
