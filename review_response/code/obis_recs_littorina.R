# make distribution map for Littorina littorea


# libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
#remotes::install_github("ropensci/taxize")
library(taxize)
library(worrms)
library(robis)
library(sf)
library(ggtext)
sf_use_s2(F)
theme_set(theme_light())


# load spp ----------------------------------------------------------------
spp_pres <- readr::read_csv(
  here::here(
    "data-processed",
    "appledore-survey-data",
    "pres_spp_list.csv"
  )
) %>%
  mutate(organism = recode(organism,
                           "Nemalion elmnithoides"=  "Nemalion elminthoides"))


world <- rnaturalearth::ne_countries(returnclass = "sf") %>%
  st_union()

# upload temp data --------------------------------------------------------
temp <- readRDS(here::here(
  "data-processed",
  "global-temp-data",
  "global_temps_1982_2023.rds"
))
# note that these are in 1x1 grid cells with the .5 as the
# coordinate in each direction (e.g., 0.5 N, 0.5 W)




# isolate names variable
names <- spp_pres %>% distinct(organism) %>% pull()

# use taxize package to confirm all names are correct
names_clean <- taxize::gna_verifier(names) %>% 
  filter(submittedName == matchedCanonicalSimple) %>%
  distinct(matchedCanonicalSimple) %>% pull(matchedCanonicalSimple) %>% sort()

setequal(names, names_clean) 
# ok, both vectors are the same, so all of our names are valid
rm(names)


# get APHIA IDs -----------------------------------------------------------
aphiaids <- list()

for(i in names_clean){
  
  id <- wm_records_name(i) %>%
    filter(status == "accepted") 
  
  if(nrow(id) > 1){
    id <- id %>%
      filter(scientificname == i)
  }
  
  id <- id %>%  pull(AphiaID)
  
  aphiaids[[i]] <- id
  
  print(paste(i, id))
  
}

# get aphiaID for littorina littorea
spp_pres_taxize_lookup <- purrr::map(
  .x = aphiaids[55],
  .f = ~as_tibble(.x)) %>%
  bind_rows(.id = "organism") %>%
  rename(AphiaID_taxize = value)


spp_pres <- spp_pres %>%
  left_join(spp_pres_taxize_lookup)
# here, 4 species don't have the same AphiaID from the 
# dataset as they do from taxize, but by manually searching
# those 4, I realize the taxize ones are correct
spp_pres <- spp_pres %>%
  select(-AphiaID)
rm(spp_pres_taxize_lookup)

# get obis records --------------------------------------------------------
obis_recs_raw <- purrr::map(
  .x = aphiaids[55],
  .f = ~ .x %>%
    occurrence(taxonid = .,
               enddate = as.Date("2023-12-31")) %>%
    as_tibble(),
  .progress = T
)
obis_recs_raw[[1]]

p1_raw_records_map <- ggplot() +
  geom_sf(data = world) +
  geom_point(data = obis_recs_raw[[1]],
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             alpha = .3, 
             size = .7,
             shape = 16) +
  labs(title = paste0("<i>",names(obis_recs_raw),"</i>", " raw records"),
       subtitle = paste0("n = ", scales::comma(nrow(obis_recs_raw[[1]]))),
       x = NULL,
       y = NULL) +
  coord_sf(expand = F) +
  theme(plot.title = element_markdown())

obis_recs_raw[[1]] %>% glimpse()


obis_recs_number <- purrr::map(
  .x = obis_recs_raw,
  .f = ~nrow(.x) %>%
    as_tibble()
) %>%
  bind_rows(.id = "organism") %>%
  rename(obis_recs_raw = value)
obis_recs_number


contains_depth <- purrr::map(
  .x = obis_recs_raw,
  .f = ~"depth" %in% colnames(.x)
) %>%
  unlist()
contains_depth %>% table()


obis_recs_filt <- purrr::map(
  .x = obis_recs_raw[contains_depth],
  .f = ~.x %>% 
    dplyr::select(c("decimalLongitude", 
                    "decimalLatitude", 
                    "depth", 
                    "year", 
                    "bathymetry",
                    "sst",
                    "month",
                    "scientificName", 
                    "aphiaID")) %>% 
    mutate(depth = as.numeric(depth),
           year = formatC(year),
           month = formatC(month, width = 2, flag = "0"),
           depth0 = case_when(
             is.na(depth) ~ 0,
             depth < 0 ~ 0,
             TRUE ~ depth))
)
obis_recs_filt[[1]] %>% glimpse()


obis_recs_filt2 <- purrr::map(
  .x = obis_recs_filt,
  .f = ~.x %>%
    mutate(year2 = as.numeric(year)) %>%
    filter(!is.na(year2))
)
# many years come out as NA

obis_recs_with_year <- purrr::map(
  .x = obis_recs_filt2,
  .f =~nrow(.x) %>%
    as_tibble
) %>%
  bind_rows(.id = "organism") %>%
  rename(obis_recs_with_year = value)

obis_recs_number <- obis_recs_number %>%
  left_join(obis_recs_with_year)
rm(obis_recs_with_year)

# filter to depths shallower than 10m
# and year after 1981
obis_recs_filt3 <- 
  purrr::map(
    .x = obis_recs_filt2,
    .f = ~.x %>%
      filter(depth0 < 10,
             year2 > 1981,
             year2 <= 2023)
  )

obis_recs_with_depth <- purrr::map(
  .x = obis_recs_filt3,
  .f =~nrow(.x) %>%
    as_tibble
) %>%
  bind_rows(.id = "organism") %>%
  rename(obis_recs_with_depth = value)



obis_recs_number <- obis_recs_number %>%
  left_join(obis_recs_with_depth)
rm(obis_recs_with_depth)
obis_recs_number %>%
  arrange(obis_recs_with_depth)

rm(obis_recs_filt2, obis_recs_filt, obis_recs_raw)



# match to temps ----------------------------------------------------------
obis_recs_matched <- vector(mode = "list",
                            length = length(obis_recs_filt3))
length(obis_recs_matched)
length(obis_recs_filt3)
names(obis_recs_matched) <- names(obis_recs_filt3)

# change lat and lon to have coordinates ending in .5 so that they can be 
# matched to the year and coordinate of SST data. 
for(i in 1:length(obis_recs_filt3)){
  
  df <- furrr::future_pmap(
    .l = list(year = as.character(obis_recs_filt3[[i]]$year2), 
              depth = obis_recs_filt3[[i]]$depth0,
              lat = stringr::str_replace(obis_recs_filt3[[i]]$decimalLatitude,"\\.[0-9]*",".5"),
              lon = stringr::str_replace(obis_recs_filt3[[i]]$decimalLongitude,"\\.[0-9]*",".5")
    ),
    .f = function(year, depth, lat, lon) {
      temp[[year]][temp[[year]]$longitude == lon &
                     temp[[year]]$latitude == lat , ]
      
    },
    .progress = T
  ) %>%
    bind_rows(.id = "rownum")
  
  obis_recs_matched[[i]] <- df
  print(i)
  
  rm(df)
  gc()
  
}

obis_recs_with_temp_match <- purrr::map(
  .x = obis_recs_matched,
  .f =~nrow(.x) %>%
    as_tibble
) %>%
  bind_rows(.id = "organism") %>%
  rename(obis_recs_with_temp_match = value)

obis_recs_number <- obis_recs_number %>%
  left_join(obis_recs_with_temp_match)


# remove any zero length matches
non_zero_matches <- obis_recs_with_temp_match %>% 
  filter(obis_recs_with_temp_match>0) %>%
  pull(organism)

obis_recs_matched_nonzero <- obis_recs_matched[non_zero_matches]




# match back into occurrence dataset --------------------------------------
# add rowname column to original dataset so we can merge them
# (we're doing this so we can keep the coordinates to plot)
obis_recs_filt4 <- 
  purrr::map(
    .x = obis_recs_filt3[non_zero_matches],
    .f = ~.x %>% 
      mutate(rownum = as.character(c(1:n()))) %>%
      select(-year) %>%
      rename(year = year2)
  )
rm(obis_recs_filt3)

# now join both datasets 
obis_recs_joined <-
  purrr::map(
    .x = names(obis_recs_filt4),
    .f = ~ obis_recs_filt4[[.x]] %>%
      left_join(obis_recs_matched_nonzero[[.x]],
                by = join_by(year, rownum))
  )
names(obis_recs_joined) <- names(obis_recs_filt4)
obis_recs_joined[[1]] %>% head()


obis_recs_joined[[1]] %>% 
  mutate(longitude = purrr::map(.x = longitude, .f = ~unlist(.x)),
         latitude = purrr::map(.x = latitude, .f = ~unlist(.x))) %>% 
  slice(5000) %>%
  glimpse()

obis_recs_joined[[1]]$longitude <- purrr::map(.x = obis_recs_joined[[1]]$longitude,
                                              .f = ~unlist(.x))
obis_recs_joined[[1]]$latitude <- purrr::map(.x = obis_recs_joined[[1]]$latitude,
                                              .f = ~unlist(.x))


p2_matched_records_map <- ggplot() +
  geom_sf(data = world) +
  geom_point(data = obis_recs_joined[[1]],
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             alpha = .3, 
             size = .7,
             shape = 16) +
  labs(title = paste0("<i>",names(obis_recs_joined),"</i>", " filtered & matched records"),
       subtitle = paste0("n = ", scales::comma(nrow(obis_recs_joined[[1]]))),
       x = NULL,
       y = NULL) +
  coord_sf(expand = F) +
  theme(plot.title = element_markdown())


good_recs <- obis_recs_joined[[1]] %>%
  filter((decimalLongitude > -85 & 
            decimalLongitude < -40 &
            decimalLatitude > 30 &
            decimalLatitude < 60) |
           (decimalLongitude > -15 & 
              decimalLongitude < 55 &
              decimalLatitude > 35 &
              decimalLatitude < 80)  )

bad_recs <- obis_recs_joined[[1]] %>%
  filter(!((decimalLongitude > -85 & 
              decimalLongitude < -40 &
              decimalLatitude > 30 &
              decimalLatitude < 60) |
             (decimalLongitude > -15 & 
                decimalLongitude < 55 &
                decimalLatitude > 35 &
                decimalLatitude < 80)  ))


p3_matched_records_map_good_bad <- ggplot() +
  geom_sf(data = rnaturalearth::ne_countries(returnclass = "sf") %>%
            sf::st_union()) +
  geom_point(data = bad_recs,
             aes(x = decimalLongitude,
                 y = decimalLatitude,
                 color = mean_monthly_sst),
             alpha = .3, 
             size = .7,
             shape = 16,
             color = "red") +
  geom_point(data = good_recs,
             aes(x = decimalLongitude,
                 y = decimalLatitude,
                 color = mean_monthly_sst),
             alpha = .3, 
             size = .7,
             shape = 16,
             color = "blue") +
  annotate(geom = "rect",
           xmin = -85,
           xmax = -40,
           ymin = 30, 
           ymax = 60,
           color = "black",
           fill = "transparent")+
  annotate(geom = "rect",
           xmin = -15,
           xmax = 55,
           ymin = 35, 
           ymax = 80,
           color = "black",
           fill = "transparent") +
  labs(x = NULL,
       y = NULL,
       title = "<i>Littorina Littorea</i> likely range records",
       subtitle = paste0("n in squares = ", scales::comma(nrow(good_recs)),
                         "; n outside = ", scales::comma(nrow(bad_recs)))) +
  theme(plot.title = element_markdown()) +
  coord_sf(expand = F)


p_joined <- cowplot::plot_grid(p1_raw_records_map,
                               p2_matched_records_map,
                               p3_matched_records_map_good_bad,
                               labels = c("a","b","c"),
                               nrow = 1)
p_joined + ggview::canvas(12,2.5)
ggsave(p_joined, 
       filename = "review_response/figures/littorina_map.png",
       width = 12, height = 2.5)


# calculate STIs using good and "bad" records -------------------------------
mean(obis_recs_joined[[1]]$mean_monthly_sst, na.rm=T)
mean(good_recs$mean_monthly_sst, na.rm=T)

quantile(
  obis_recs_joined[[1]]$max_monthly_sst, na.rm=T,
  .95)
quantile(
  good_recs$max_monthly_sst, na.rm=T,
  .95)

quantile(
  obis_recs_joined[[1]]$min_monthly_sst, na.rm=T,
  .05)
quantile(
  good_recs$min_monthly_sst, na.rm=T,
  .05)

mean(obis_recs_joined[[1]]$mean_monthly_sst, na.rm=T)
mean(good_recs$mean_monthly_sst, na.rm=T)

