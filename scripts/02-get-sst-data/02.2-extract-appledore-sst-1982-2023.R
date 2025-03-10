# workflow to find the actual temperature on Appledore Island during the study period


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
appledore_temp <- griddap(datasetx = "erdHadISST", 
                 url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                 time = c("1981-01-01", "2023-12-31"), 
                 latitude = c(42.5, 42.5),
                 longitude = c(-70.75, -70.5),
                 fields = "sst")

data <- appledore_temp$data[appledore_temp$data$sst > -999,]
data$date <- lubridate::as_date(data$time)
data$year <- lubridate::year(data$date)
data$month <- lubridate::month(data$date)
data <- data %>% select(-time)
data %>% glimpse()
data <- data %>% tidyr::drop_na(date)
range(data$date)
rm(appledore_temp)
gc()

data %>% distinct(latitude)

data %>% filter(year == 1982,
                month == 1) %>%
  ggplot() +
  geom_raster(aes(x = longitude,
                  y = latitude,
                  fill = sst)) +
  coord_equal() +
  geom_sf(data = app,
          color = "red",
          fill = "red")



# add offset year ---------------------------------------------------------
# because appledore samples are taken in July or August, 
# we're going to make the "year" temperatures for each year of sampling
# represent the year preceding the sample, so septempber n-1 to august n
data <- data %>% 
  mutate(sample_year = case_when(month >= 9 ~ year+1,
                                 month <= 8 ~ year)) 

data_sum <- data %>%
  group_by(sample_year) %>%
  group_by(sample_year) %>%
  summarize(min = min(sst),
            max = max(sst),
            mean = mean(sst),
            median = median(sst)) %>%
  filter(sample_year > 1981,
         sample_year < 2024) 

# save as csv
readr::write_csv(
  data_sum,
  here::here(
  "data-processed",
  "appledore-island-env-data",
  "appledore_temps_1982_2023.csv"
))
