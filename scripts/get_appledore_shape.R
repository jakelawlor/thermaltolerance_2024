
# script to get appledore shape ------------------------------------------------

library(sf)
library(dplyr)
library(ggplot2)

# begin by uploading maine government shapefile for counties,
# downloaded from here: https://maine.hub.arcgis.com/datasets/b47e1e44ce70411f91fe4c1c3f4c7cd6/explore?location=43.028428%2C-70.608461%2C11.71
maine <- read_sf(here::here(
  "data-raw",
  "Maine_County_Boundary_Polygons_Feature "
))
maine %>%
  glimpse()

maine %>%
  filter(COUNTY == "York") %>%
  filter(ISLAND == "y") %>%
  #st_crop(xmin=-71, 
  #        xmax=-70,
  #        ymin=41, 
  #        ymax=43.0) %>%
  filter(OBJECTID == "6983") %>%
  st_transform(crs = st_crs("EPSG:4326")) %>%
  #distinct(OBJECTID) %>%
  ggplot() +
  geom_sf(aes(fill= as.character(OBJECTID)))

  
app <- maine %>%
  filter(COUNTY == "York") %>%
  filter(ISLAND == "y") %>%
  #st_crop(xmin=-71, 
  #        xmax=-70,
  #        ymin=41, 
  #        ymax=43.0) %>%
  filter(OBJECTID == "6983") %>%
  st_transform(crs = st_crs("EPSG:4326"))

# save appledore ----------------------------------------------------------
appledore <- maine %>%
  filter(OBJECTID == "6983")

saveRDS(appledore,
        file = here::here("data-processed",
                          "appledore_shp.rds"))

