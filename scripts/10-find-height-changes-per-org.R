# Assess when/where species enter and exit


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)


# data --------------------------------------------------------------------
pa <- readr::read_csv(
  here::here(
    "data-processed",
    "spp_pres_by_replicate.csv"
  )
)

therm <- readr::read_csv(
  here::here(
    "data-processed",
    "spp_thermal_affinities.csv"
  )
)

pa <- pa %>% left_join(therm) %>%
  filter(pres == T) %>%
  # translate level into tidal height
  mutate(tidalheight = (13-level)*.348)  
pa %>% distinct(data_taken)


# find cti of new spp -----------------------------------------------------
# first, let's break up the timeseries 
range(pa$year)
pa %>% 
  group_by(organism) %>%
  filter(year == min(year)) %>% 
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_group = case_when(
    year <= 1989 ~ "1980s",
    year >= 1990 & year <= 1999 ~ "1990s",
    year >= 2000 & year <= 2009 ~ "2000s",
    year >= 2010 & year <= 2019 ~ "2010s",
    year >= 2020 ~ "2020s"
  )) %>%
  ggplot(aes(x = year,
             y = mean_monthly_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year of First Occurrence",
       y = "Thermal Affinity") 



pa %>% 
  group_by(organism) %>%
  filter(year == max(year)) %>% 
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_group = case_when(
    year <= 1989 ~ "1980s",
    year >= 1990 & year <= 1999 ~ "1990s",
    year >= 2000 & year <= 2009 ~ "2000s",
    year >= 2010 & year <= 2019 ~ "2010s",
    year >= 2020 ~ "2020s"
  )) %>%
  filter(year != max(year)) %>%
  ggplot(aes(x = year,
             y = mean_monthly_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year of First Occurrence",
       y = "Thermal Affinity") 




# find how long species lasted in the dataset -----------------------------
first_year <- pa %>%
  group_by(organism) %>%
  filter(year == min(year)) %>%
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_type = "first_year")

last_year <- pa %>%
  group_by(organism) %>%
  filter(year == max(year)) %>%
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_type = "last_year")

year_breadth <-  
  first_year %>%
  rbind(last_year)


year_breadth %>% 
  tidyr::pivot_wider(names_from = year_type,
                     values_from = year) %>% 
  filter(last_year - first_year > 0) %>% 
  ggplot() +
  geom_segment(aes(x = first_year,
                   xend = last_year,
                   y = mean_monthly_mean,
                   yend = mean_monthly_mean),
               linewidth = .3,
               alpha = .4) +
  geom_point(data = . %>% filter(first_year > min(first_year)),
             aes(x = first_year,
                 y = mean_monthly_mean),
             color = "cyan4") +
  geom_point(data = . %>% filter(last_year < max(last_year)),
             aes(x = last_year,
                 y = mean_monthly_mean),
             color = "tomato2") +
  geom_smooth(data = . %>% filter(first_year > min(first_year)),
              aes(x = first_year,
                  y = mean_monthly_mean),
              color = "cyan4",
              method = "lm") +
  geom_smooth(data = . %>% filter(last_year < max(last_year)),
              aes(x = last_year,
                  y = mean_monthly_mean),
              color = "tomato2",
              method = "lm")

pa %>%
  distinct(transect, level, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = level)) +
  geom_tile() +
  facet_wrap(~transect)
# longest most consistent transects are:
# 5, 7, 15, 20, 22, 26, 28
# between levels 0 and 15
good_transects <- c(5, 7, 15, 20, 22, 26, 28)
min_level <- 0
max_level <- 15

# track individual species min/max depths ---------------------------------

pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>% 
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  ggplot(aes(x = year,
             group = organism)) +
  geom_line(aes(y = minheight),
             color = "grey80",
            linewidth = .2,
            alpha = .3) +
  geom_smooth(aes(y = minheight,
                  group = organism),
              method = "lm",
              se = F,
              color = "cyan4",
              linetype = "dashed",
              linewidth = .5) 

pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>% 
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight),
            mean_monthly_mean = unique(mean_monthly_mean)) %>%
  group_by(organism) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(mean_monthly_mean),
         n > 3) %>%
  mutate(organism_f = forcats::fct_reorder(organism, mean_monthly_mean)) %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = minheight),
            color = "tomato2",
            linewidth = .5,
            alpha = .3) +
  geom_line(aes(y = maxheight),
            color = "cyan4",
            linewidth = .5,
            alpha = .3) +
# geom_smooth(aes(y = maxheight),
#             method = "lm",
#             se = F,
#             color = "tomato2",
#             linetype = "dashed",
#             linewidth = .5) +
  facet_wrap(~organism_f)


pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>%   
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  ggplot(aes(x = year,
             group = organism)) +
  geom_ribbon(aes(ymin = minheight,
                  ymax = maxheight)) +
  facet_wrap(~organism)


pa %>%
  filter(!is.na(mean_monthly_mean)) %>%
  mutate(organism_f = forcats::fct_reorder(organism, mean_monthly_mean)) %>%
  filter(transect %in% good_transects,
        level >= min_level,
        level <= max_level) %>% 
  group_by(organism, organism_f, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  group_by(organism) %>%
  add_count() %>%
  ungroup() %>%
  filter(n > 3) %>%
  ggplot(aes(x = year)) +
  geom_smooth(aes(y = maxheight),
              method = "lm", 
              se = F,
              color = "tomato2") +
  geom_smooth(aes(y = minheight),
              method = "lm", 
              se = F,
              color = "cyan4") +
  facet_wrap(~organism_f)
