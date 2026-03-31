# script to filter to only highly sampled transects 
# to more fairly assess thermal affinity. 


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
theme_set(ggthemes::theme_few())

# data --------------------------------------------------------------------
pa <- readr::read_csv(
  here::here(
    "data-processed",
    "appledore-survey-data",
    "spp_pres_by_replicate.csv"
  )
)
pa %>% filter(pres == T) %>% distinct(organism) %>% nrow()
pa %>% glimpse()
# 93 species have presences to start with

therm <- readr::read_csv(
  here::here(
    "data-processed",
    "species-thermal-affinities",
    "spp_thermal_affinities.csv"
  )
)
therm %>% filter(!is.na(mean_monthly_mean)) %>% distinct(organism) %>% nrow()
therm %>% glimpse()
# 89 species total have therm

pa_therm <- pa %>% left_join(therm) %>%
  filter(pres == T) %>%
  # translate level into tidal height
  mutate(tidalheight = (13.5-level)*.3048)  %>%
  # filter out species that we don't have thermal affinity for
 # filter(is.na(mean_monthly_mean)) %>% # 5 organisms total get removed (4 don't have therm)
  mutate(transect_label = paste0("Transect ",transect)) %>%
  arrange(transect) %>%
  mutate(transect_label = forcats::fct_inorder(transect_label))

pa_therm %>% distinct(organism) %>% nrow()
# 88 species have both...
# note that for 4 species that are present, we don't have therm,
# which makes sense, but for 1 species that we DO ahve therm for,
# it's never present... that species is Nemalion elminthoides. 
# maybe the species list got processed for therm before filtering
# had been done. Possibly check this out later.
unique(therm$organism) %in% unique(pa$organism) %>% table()
therm$organism[!therm$organism %in% pa$organism]
rm(pa, therm)


# upload emerical appledore temps
temps <- 
  readr::read_csv(here::here("data-processed",
                             "appledore-island-env-data",
                             "appledore_temps_1982_2023.csv")) %>%
  mutate(year = sample_year)
tempmod <- lm(data = temps,
              mean ~ year)
summary(tempmod) # 0.035722 increase in temp per year
#first_temp <- predict(tempmod, list(year= min(pa_therm$year)))
#last_temp <- predict(tempmod, list(year= max(pa_therm$year)))
temptrend <- broom::augment(tempmod, interval = "confidence")
rm(temps)

# note that here we don't have to rarify number of replicates,
# assuming that sampling effort has no correlation to CTI
# (as opposed to richness, where more sampling = higher richness)


# cut to highly sampled tidal heights -------------------------------------
# next, identify highly-sampled tidal heights
pa_therm %>% distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point() +
  facet_wrap(~transect_label) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2.75)
sort(unique(pa_therm$tidalheight))
# this leaves 9 tidal heights between these bounds (0-2.75). 


pa_therm %>% distinct(tidalheight) %>% 
  count(tidalheight >= 0 & tidalheight < 2.75)

pa_therm_bound <- pa_therm %>%
  filter(tidalheight >= 0 & tidalheight < 2.75)



# filter out undersampled years ------------------------------------------------------------
# within these tidal height bounds, there are 9 tidal heights sampled. 
# cut out years in which there are 6 or fewer of these 9 sampled. 
pa_therm_bound %>%
  distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point() +
  facet_wrap(~transect_label) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2.75)

transect_years_to_keep <- pa_therm_bound %>% 
  group_by(year, transect_label) %>% 
  distinct(tidalheight) %>% 
  count(name = "levels_sampled") %>%
  mutate(keep = ifelse(levels_sampled > 6, "yes","no")) %>%
  ungroup()

pa_therm_bound_yearskeep <- pa_therm_bound %>%
  left_join(transect_years_to_keep) %>%
  filter(keep == "yes")
rm(transect_years_to_keep)

# plot
pa_therm %>% distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point(color = "black",
             shape = 21,
             fill = "transparent",
             stroke = .2) +
  facet_wrap(~transect_label) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2.75) +
  geom_point(data = pa_therm_bound_yearskeep,
             color = "black",
             fill = "#aaeeffff",
             shape = 21,
             stroke = .1)



# filter out short timeseries ---------------------------------------------
# some transects are measured for only ~3 years in the beginning of the 
# time series. These aren't going to form a good trend in such little time. 
# filter to only transects in which sampling occurs over 
# cut out transects with fewer than 20 years of sampling
transects_to_keep <- 
  pa_therm_bound_yearskeep %>%
  group_by(transect_label) %>%
  distinct(year) %>%
  count() %>%
  mutate(keep_transect = ifelse(n < 20,"no","yes")) %>% ungroup()

transects_to_keep %>% count(keep_transect)

pa_therm_bound_yearskeep_transectkeep <- 
  pa_therm_bound_yearskeep %>%
  left_join(transects_to_keep) %>%
  filter(keep_transect == "yes") 
rm(transects_to_keep)

# remove unneeded objects
pa_therm_filtered <- pa_therm_bound_yearskeep_transectkeep
rm(pa_therm_bound, pa_therm_bound_yearskeep, pa_therm_bound_yearskeep_transectkeep)

# plot again
pa_therm %>% distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point(color = "black",
             shape = 21,
             fill = "transparent",
             stroke = .2) +
  facet_wrap(~transect_label, nrow = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2.75) +
  geom_point(data = pa_therm_filtered,
             color = "black",
             fill = "#aaeeffff",
             shape = 21,
             stroke = .1)


# cut out years when fewer than 3 transects are included ------------------
# since we're going to average cti for the whole island each year,
# lets cut to years when at least 3 of the 7 remaining transects 
# are sampled, since we don't want an island-wide cti to be from
# just one or two transects
years_to_keep2 <- pa_therm_filtered %>%
  group_by(year) %>% 
  distinct(transect_label) %>% 
  count() %>% 
  filter(n >= 3) %>% pull(year)


pa_therm_filtered2 <- pa_therm_filtered %>%
  filter(year %in% years_to_keep2)
rm(years_to_keep2)

# plot again
filtered_transect_p <- pa_therm %>% distinct(transect_label, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  geom_point(color = "grey70",
             shape = 16,
             size = .65,
             fill = "transparent",
             stroke = .1) +
  facet_wrap(~transect_label, nrow = 3) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = .25) +
  geom_hline(yintercept = 2.75,
             linetype = "dashed",
             linewidth = .25) +
  geom_point(data = pa_therm_filtered2,
             size = .65,
             color = "darkturquoise",
             shape = 16,
             stroke = .03) +
  scale_x_continuous(breaks = c(1980, 2000, 2020)) +
  
  labs(y = "Tidal Level (m)",
       x = NULL,
       title = "Retained Highly Sampled Transects") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, 
                                   vjust = 1),
        panel.spacing.y = unit(0,"mm"),
        strip.text = element_text(margin = ggplot2::margin(t=7,0,b=2,0,"pt")),
        plot.title = element_text(margin = ggplot2::margin(b=0)))

filtered_transect_p + ggview::canvas(8,3.5)

ggsave(filtered_transect_p,
       filename = "outputs/data-plots/highly_sampled_transects.png",
       width = 8,
       height = 3.5)

rm(pa_therm_filtered)
pa_therm_filtered <- pa_therm_filtered2
rm(pa_therm_filtered2)

pa_therm %>% distinct(year) %>% count()
pa_therm %>% glimpse()
pa_therm_filtered %>% glimpse()


readr::write_csv(pa_therm,
                 file = "data-processed/appledore-survey-data/pa-with-therm/pa-with-therm-all.csv")
readr::write_csv(pa_therm_filtered,
                 file = "data-processed/appledore-survey-data/pa-with-therm/pa-with-therm-highly-sampled.csv")
saveRDS(tempmod,
        file = "data-processed/appledore-island-env-data/temp-change-model.rds")
readr::write_csv(temptrend,
        file = "data-processed/appledore-island-env-data/temp-change-augment.csv")
rm(list = ls())
