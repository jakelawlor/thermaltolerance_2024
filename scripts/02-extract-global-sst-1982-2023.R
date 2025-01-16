# workflow to find the actual temperature 
# at the actual time (year) and place (grid cell)
# of OBIS sampling locations. 


#using vignette at https://cran.r-project.org/web/packages/heatwaveR/vignettes/OISST_preparation.html

library(dplyr) # A staple for modern data management in R
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing
library(purrr)
setwd(here::here())


# get global yearly average temps -----------------------------------------

# The information for the NOAA OISST data
rerddap::info(datasetid = "erdHadISST", 
              url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# pull all erdHadiSST data for the entire world 
# from Jan 1 1982 to dec 31, 2023
# Average Sea Surface Temperature, 1°, Global, Monthly, 1870-present
# here, we're going to take mean, min, and max monthly average sst
temps <- griddap(datasetx = "erdHadISST", 
                 url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                 time = c("1981-01-01", "2023-12-31"), 
                 latitude = c(-89.5, 89.5),
                 longitude = c(-179.5, 179.5),
                 fields = "sst")

data <- temps$data[temps$data$sst > -999,]
data$date <- lubridate::as_date(data$time)
data$year <- lubridate::year(data$date)
data$month <- lubridate::month(data$date)
data <- data %>% select(-time)
data %>% glimpse()
data <- data %>% tidyr::drop_na(date)
range(data$date)
rm(temps)
gc()



data %>% 
  filter(year == 1981, 
         month == 1) %>% 
  ggplot(aes(x = longitude,
             y = latitude,
             fill = sst)) +
  geom_raster() +
  coord_equal()





# split by year -----------------------------------------------------------
data_split <- data %>%
  split(data$year)


# pivot so each month is a column
data_split_pivoted <- purrr::map(
  .x = data_split,
  .f = ~ .x %>%
    select(-date) %>%
    tidyr::pivot_wider(
      names_from = month,
      values_from = sst,
      names_prefix = "month_"
    ) ,
  .progress = T
)
rm(data_split)


# find min, max, mean of monthly temps
data_split_summarized <- purrr::map(
  .x = data_split_pivoted,
  .f = ~.x %>% 
    mutate(
      # Compute the minimum temperature
      min_monthly_sst = reduce(across(starts_with("month")), pmin, na.rm = TRUE),
      # Compute the maximum temperature
      max_monthly_sst = reduce(across(starts_with("month")), pmax, na.rm = TRUE),
      # Compute the mean temperature
      mean_monthly_sst = rowMeans(across(starts_with("month")), na.rm = TRUE)
    ),
  .progress = T
)
rm(data_split_pivoted)

# remove extra columns
data_split_summarized2 <- purrr::map(
  .x = data_split_summarized,
  .f = ~.x %>% 
    ungroup() %>%
    select(-starts_with("month_"))
)

# view why they have different cells - not sure still
purrr::map_vec(
  .x = data_split_summarized2,
  .f = ~nrow(.x)
)


library(viridis)
data_split_summarized2[[20]] %>%
  ggplot(aes(x = longitude,
             y = latitude,
             fill = mean_monthly_sst)) +
  geom_raster() +
  scale_fill_gradientn(colors = viridis(9), 
                       limits=c(-2, 35), 
                       na.value = "grey50") +
  coord_equal()


# save --------------------------------------------------------------------
saveRDS(data_split_summarized2,
        file = here::here("data-processed",
                          "global_temps_1982_2023.rds"
        ))

names(data_split_summarized2)

rm(list = ls())
