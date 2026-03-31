# Script to match occurrence records from OBIS
# to the min, mean, max annual temperatures on the
# years that they were collected.


# libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
#remotes::install_github("ropensci/taxize")
library(taxize)
library(worrms)
library(robis)


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

taxize::gna_verifier(names) %>% select(submittedName, matchedName, matchType, matchedCanonicalSimple, matchedCanonicalFull)

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

spp_pres_taxize_lookup <- purrr::map(
  .x = aphiaids,
  .f = ~as_tibble(.x)) %>%
  bind_rows(.id = "organism") %>%
  rename( AphiaID_taxize = value)


spp_pres <- spp_pres %>%
  left_join(spp_pres_taxize_lookup)
# here, 4 species don't have the same AphiaID from the 
# dataset as they do from taxize, but by manually searching
# those 4, I realize the taxize ones are correct
spp_pres <- spp_pres %>%
  select(-AphiaID)
rm(spp_pres_taxize_lookup)

# set up parallel processing ----------------------------------------------
# pause to source script to set up parallel processing
source(here::here("quick-scripts","setup_parallel.R"))



# get obis records --------------------------------------------------------
obis_recs_raw <- furrr::future_map(
  .x = aphiaids,
  .f = ~ .x %>%
    occurrence(taxonid = .,
               enddate = as.Date("2023-12-31")) %>%
    as_tibble(),
  .progress = T
)
obis_recs_raw[[1]]

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
# 2 species don't contain depth values, 
# so we'll remove those for now

obis_recs_raw[contains_depth]

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
rm(aphiaids2)



# match to temps ----------------------------------------------------------
obis_recs_matched <- vector(mode = "list",
                            length = length(obis_recs_filt3))
length(obis_recs_matched)
length(obis_recs_filt3)
names(obis_recs_matched) <- names(obis_recs_filt3)

for(i in 1:length(obis_recs_filt3)){
  
  df <- furrr::future_pmap(
    .l = list(year = as.character(obis_recs_filt3[[i]]$year2), 
              depth = obis_recs_filt3[[i]]$depth0,
              lat = stringr::str_replace(obis_recs_filt3[[i]]$decimalLatitude,"\\.[0-9]*",".5"),
              lon = stringr::str_replace(obis_recs_filt3[[i]]$decimalLongitude,"\\.[0-9]*",".5")
    ),
    .f = function(year, depth, lat, lon) {
      # yeardat <-  temp[[year]]
      # yeardat[yeardat$longitude == stringr::str_replace(lon,
      #                                                   "\\.[0-9]*",".5") &
      #           yeardat$latitude == stringr::str_replace(lat,
      #                                                    "\\.[0-9]*",".5") , ]
      #
      temp[[year]][temp[[year]]$longitude == lon &
                     temp[[year]]$latitude == lat , ]
      
    },
    # .id = "rowname",
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
  filter(!is.na(min_monthly_sst)) %>%
  glimpse()

# explore some affinities -------------------------------------------------

obis_recs_joined[[1]] %>%
  filter(!is.na(mean_monthly_sst)) %>%
  ggplot(aes(x = mean_monthly_sst)) +
  geom_histogram(alpha = .7) +
  geom_histogram(aes(x = max_monthly_sst),
                 fill = "tomato2",
                 alpha = .5) +
  geom_histogram(aes(x = min_monthly_sst),
                 fill = "lightblue",
                 alpha = .8) +
  theme_minimal()


i <- 10
temp[[1]]%>%
  ggplot(aes(x = longitude,
             y = latitude,
             fill = mean_monthly_sst)) +
  geom_raster() +
  geom_point(inherit.aes = F,
             data = obis_recs_joined[[i]],
             aes(y = decimalLatitude,
                 x = decimalLongitude,
                 color = is.na(mean_monthly_sst)),
  ) 
# there are a lot of records here that are just outside
# of temperature raster cells. these are counting as NA in
# the matched dataset. In the future, switching to the 
# higher resolution sst map might help this. 


# save out list before I forget -------------------------------------------
saveRDS(
  obis_recs_joined,
  file = here::here("data-processed",
                    "species-thermal-affinities",
                    "obis_recs_matched_temp.rds")
)

obis_recs_joined <- readRDS( here::here("data-processed",
                                        "species-thermal-affinities",
                                        "obis_recs_matched_temp.rds"))

sapply(obis_recs_joined, nrow) %>% unname() %>% sort()

readr::write_csv(
  obis_recs_number,
  here::here("data-processed",
             "species-thermal-affinities",
             "number_obis_recs_per_spp.csv")
)

rm(list = ls())