# Script to prep count species for modeling


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(purrr)



# data --------------------------------------------------------------------
counts <- readr::read_csv(here::here(
  "data-processed",
  "appledore-survey-data",
  "counts_abundance_filtered.csv"
))
counts %>% distinct(organism)
# 36 species to start

# note that sometimes counts were listed as "p", or other,
# which are included in the presence/absence dataset, but 
# are changed to NA here.
counts %>%
  filter(is.na(countnum)) %>%
  distinct(count)




# find mean between replicates --------------------------------------------
# because not all quadrats had replicates, we'll calculate a mean count per
# quadrats (e.g., averaging the 1 to 4 replicates of one transect x tidalheight)
counts2 <- counts %>%
  group_by(year, transect, level, organism) %>% 
  summarize(n_rep = n(),
            countnum = mean(countnum,na.rm=T),
            pres = any(pres),
            in_sample = unique(in_sample),
            .groups = "drop") %>% 
  ungroup() 
counts2 %>% filter(is.na(countnum))
rm(counts)
range(counts2$n_rep)
mean(counts2$n_rep)



# find mean per level -----------------------------------------------------
# because transects and tidalheights were sampled inconsistently over time,
# we are going to summarize data as the mean count per tidal level in each year,
# regardless of how many transects were sampled in that year.
counts3 <- counts2 %>%
  group_by(organism, year, level) %>%
  summarize(countnum = mean(countnum, na.rm = T),
            n_transects = n()) %>%
  ungroup() %>%
  
  # add recentered year to make plotting easier
  mutate(year_zero = year - ( min(year)-1)) %>%
  
  # translate level into tidal height
  mutate(tidalheight = (13-level)*.348)  

counts3 %>% glimpse()
counts3 %>% filter(is.na(countnum))
range(counts3$n_transects)
mean(counts3$n_transects)
rm(counts2)


# filter to only species present in 3+ years of the data ------------------
# make a dataset of species and how many years they had a >0 count
multiyear_spp <- counts3 %>% 
  filter(countnum > 0) %>%
  group_by(organism) %>%
  distinct(year) %>%
  count(name = "years_pres")
multiyear_spp %>% glimpse()

counts4 <- counts3 %>%
  left_join(multiyear_spp) %>%
  filter(years_pres > 3)
counts4 %>% distinct(organism)
# NOTE: by subsetting the data to only organisms in 3+ years of sampling, 
# we cut down to 26 species
rm(counts3)
# view the species that were removed
multiyear_spp %>%
  filter(years_pres <= 3)
rm(multiyear_spp)




# remove tidal height singletons ------------------------------------------
# because we have level as a random factor in models, remove instances
# in which organisms are only >0 count in a given level once. 
counts5 <- counts4 %>%
  group_by(organism, level) %>% 
  add_count(name = "nrow_in_level") %>%
  filter(nrow_in_level > 1) %>%
  select(-nrow_in_level) %>%
  ungroup()
rm(counts4)
counts5 %>% distinct(organism)
# still 26 species



# remove tidalheights that end sampling -----------------------------------
# some tidal levels end sampling, or are just not sampled as much as the other
# levels. (particularly level -1 and 15). Remove the tidal levels that are 
# sampled in fewer than 20 years of the 40 year dataset. 

# make dataframe of how many years each tidal level has data
levels_sampled <- 
  counts5 %>%
  group_by(level) %>%
  distinct(year) %>% 
  count(name = "n_years_level_sampled")

# filter counts to tidal levels sampled in 20+ years
counts6 <- counts5 %>%
  left_join(levels_sampled) %>%
  filter(n_years_level_sampled > 20) %>%
  select(-n_years_level_sampled)
rm(counts5, levels_sampled)
counts6 %>% distinct(level,tidalheight) %>% arrange(level)
# we end with levels 0-14 (e.g., tidalheights  4.52 to -3.48ft)



# translate count to density per m sq ----------------------------------------
# counts were measured with 20cm x 20cm (400cm2) quadrats.
# to scale this up to density per m2, we will divide by 
counts7 <- counts6 %>%
  mutate(density = countnum / 400 * 10000) %>%
  select(-countnum) %>%
  mutate(density_nonzero = ifelse(density == 0, 0.01, density)) %>%
  filter(!is.na(density_nonzero))
rm(counts6)

readr::write_csv(counts7,
                 here::here("data-processed",
                            "abundance-data",
                            "count-data-prepped-for-model.csv"))

rm(list = ls())

#  # 3. Subset to test model workflow -------------------------------------------------------
#  # subset data to only one species to test modeling function
#  testdf <- counts6 %>%
#    filter(organism == "Pagurus acadianus") %>%
#    filter(!is.na(density)) 
#  hist(testdf$density)
#  hist(log10(testdf$density))
#  
#  testdf %>%
#    filter(!is.na(density)) %>%
#    group_by(tidalheight) %>%
#    count(density > 0) %>%
#    rename(pres = `density > 0` ) %>%
#    mutate(prop = n / sum(n) * 100) %>%
#    select(-n) %>%
#    tidyr::pivot_wider(names_from = pres,
#                       values_from = c(prop))
#  
#  testdf %>%
#    ggplot(aes(x=year, y=density,
#               group = tidalheight)) +
#    geom_point() +
#    geom_smooth(method="lm",
#                se = F)
#  
#  testdf %>%
#    ggplot(aes(x=year, y=density)) +
#    geom_point() +
#    geom_smooth(method="lm") +
#    facet_wrap(~tidalheight, scales = "free_y") 
#  # in this case, trends show fairly consistent declines in density over time
#  # in most tidal levels, but major differences in magnitude between levels



# # plot trends over time with just lm trend line
# example_unmodeled <- 
#   testdf %>%
#   #group_by(level) %>% add_count() %>% 
#   # filter(n>3) %>% filter(level <= 13) %>%
#   mutate(height_label = paste(round(tidalheight,2), "m")) %>%
#   mutate(height_label = forcats::fct_reorder(height_label, level)) %>%
#   ggplot(aes(x=year, y=density)) +
#   geom_hline(yintercept = 0)+
#   geom_point(size=.75, alpha=.5) +
#   geom_smooth(method = "lm", color = "black", fill = "#3499ad", size=.75) +
#   theme(panel.border = element_blank()) +
#   facet_wrap(~height_label) +
#   labs(x=NULL,
#        y=NULL,
#        title = "Semibalanus balanoides density/m^2 change over time") +
#   theme(plot.title = element_text(size=18),
#         plot.title.position = "plot") 
# 
# example_unmodeled


# 
# # 4. Test model on subset ----------------------------------------------------------
# testdf <- testdf %>%
#   filter(!is.na(density)) %>%
#   mutate(density_nonzero = ifelse(density == 0, .01, density))
# brm_glm_reg <- brm(density_nonzero ~ year_zero + (1|tidalheight),
#                    data=testdf, family=Gamma(link="log"),
#                    #prior=c(prior(normal(0,2),class="Intercept"),
#                    #        prior(normal(0,2),class="b"),
#                    #        prior(gamma(0.01,0.01),class="shape")),
#                    chains=4,iter=4000, cores=2,  init = "0", inits=NULL,
#                    control = list(adapt_delta = 0.99) )
# 
# summary(brm_glm_reg)
# 
# pairs(brm_glm_reg)
# 
# # view conditional effects
# # (the effect of year on density, given the condition of tidal height at 3 levels)
# conditional_effects(brm_glm_reg, 
#                     effects = c("year_zero"),
#                     #int_conditions = list(level = c(15,7,1))
# ) 
# # in all three cases of tidal height, we see a negative trend of density over time
# 
# 
# # posterior predictive check to make sure model is reflecting values
# pp_check(brm_glm_reg) + theme_minimal()

# 
# test <- glm(data = testdf,
#            formula = "countnum_nonzero ~ year_zero + tidalheight",
#            family = Gamma(link = "log"))
# summary(test)
# 
# library(pscl)
# 
# library(pscl)
# 
# 
# model <- hurdle(countnum ~ year_zero + tidalheight, data = testdf,
#                 dist = "negbin")
# 
# install.packages("tweedie")
# library(tweedie)
# library(lme4)
# library(statmod)
# 
# model <- glm(countnum ~ year_zero  + tidalheight, 
#              family = tweedie(var.power = 1.5,
#                               link.power = 0),
#              data = testdf)
# summary(model)
# 
# au <- broom::augment(model,
#                testdf,
#                type.predict = "response") %>% glimpse()
# 
# au %>%
#   mutate(height_label = paste(round(tidalheight,2), "m")) %>%
#   mutate(height_label = forcats::fct_reorder(height_label, level)) %>%
#   ggplot(aes(x = year_zero,
#              y = countnum)) +
#   geom_point() +
#   geom_line(aes(y = .fitted)) +
#   facet_wrap(~height_label) 
# 



# 
# 
# # separate species only present in one or two levels
# counts_mult <- counts6 %>%
#   group_by(organism,level) %>% nest() %>% 
#   group_by(organism) %>% add_count() %>%
#   filter(n>2) %>%
#   select(-n) %>%
#   unnest(data)
# counts_mult %>% distinct(organism)
# 
# counts_sing <- counts6 %>%
#   group_by(organism,level) %>% nest() %>% 
#   group_by(organism) %>% add_count() %>%
#   filter(n<=2) %>%
#   select(-n) %>%
#   unnest(data)
# counts_sing %>% distinct(organism)
# # 1 species total was found in only one level
# 


