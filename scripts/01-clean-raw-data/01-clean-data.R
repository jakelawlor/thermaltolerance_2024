
# Script to change clean data



# 1. libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)



# 2. load species list -------------------------------------------------------
spp <- readr::read_csv(
  here::here(
    "data-raw",
    "EDI Data thru 2023",
    "species_list.csv"
  )
)


# remove all "(canopy)" or "(primary)" 
spp2 <- spp %>%
  mutate(valid_name = str_replace(valid_name," \\(canopy\\)",""),
         valid_name = str_replace(valid_name," \\(primary\\)","")) %>%
  distinct()

spp2 %>% glimpse()
spp2 %>% distinct(valid_name)

spp2 %>% 
  group_by(valid_name) %>%
  count() %>%
  arrange(desc(n))

rm(spp)

# 3. process counts data ---------------------------------------------------------------
counts <- readr::read_csv(
  here::here(
    "data-raw",
    "EDI Data thru 2023",
    "counts_data.csv"
  )
) %>%
  janitor::clean_names() 

## | first, cut to real organisms --------------------------------------------
counts %>% distinct(organism) %>%
  arrange(organism) %>% pull()



counts2 <- counts %>%
  # cut out one-word organisms (Amphipoda, Cancer, etc)
  filter(str_detect(organism, " ")) %>%
  # filter out "eggs"
  filter(!str_detect(organism,"\\(eggs\\)")) %>%
  # filter out when year / transect / level are NA
  tidyr::drop_na(year, transect, level) %>%
  # filter to data_taken == yes, 
  # this column is a bit weird but I think this is the best solution
  filter(data_taken == "yes") %>%
  # filter out replicates with "a" and "b",
  # and change "NA" replicates to replicate 1
  # (these cases there was only one)
  mutate(replicate = case_when(is.na(replicate) ~ "1",
                               TRUE ~ replicate)) %>%
  filter(!str_detect(replicate,"a|b")) %>% # note this is only one instance
  mutate(replicate = as.numeric(replicate)) %>%
  # filter to only whole-number levels
  filter(level %% 1 == 0)
rm(counts)


## | deal with "count" column -------------------------------------------
counts3 <- counts2 %>% 
  # change count to numeric - this will create some NAs
  mutate(countnum = as.numeric(count)) %>%
  # this leaves a bunch of cases where count was NA, 
  # and some where count was "p".
  
  # add presence column
  mutate(pres = case_when(
    countnum > 0 ~ TRUE,
    count == "p" ~ TRUE,
    count == "sp100" ~ TRUE,
    countnum == 0 ~ FALSE,
    is.na(count) ~ NA,
    count == "ND" ~ NA
  ))
rm(counts2)

# sometimes, one organism is listed twice in one year/transect/
# level/organism. Usually, as count = 0 and count == 0.0,
# but sometimes count = 0 and count = p, etc. 
# in these cases, arrange by rev(pres) so TRUEs are on top, then 
# slice only the first row
counts4 <- counts3 %>%
  group_by(year, transect, level, replicate, organism) %>%
  arrange(-pres) %>%
  slice(1) %>% ungroup() %>%
  
  # filter to only species that 
  # are present at least once
  group_by(organism) %>%
  mutate(max_pres = sum(pres, na.rm = T)) %>%
  ungroup() %>%
  filter(max_pres > 0 ) %>%
  select(-max_pres)

rm(counts3)



## | make a pres abs dataset -------------------------------------------------
counts_pa <- counts4 %>%
  select(-count, -countnum)
# this will later be merged with cover_pa for one large
# presence absence dataset


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# 4. Begin cover -------------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
cover <- readr::read_csv(
  here::here("data-raw",
             "EDI Data thru 2023",
             "percent_cover_data.csv")
) %>%
  janitor::clean_names()


# view species names
cover %>% distinct(organism) %>%
  arrange(organism) %>% pull()


## | cut to real organisms ---------------------------------------------------
cover2 <- cover %>%
  # filter out single name organisms
  filter(str_detect(organism, " ")) %>%
  # often, species have "(canopy)" and "(primary)"
  # keep only canopy
  filter(!str_detect(organism, "\\(primary\\)")) %>%
  # then remove "canopy" from remaining names
  mutate(organism = str_remove(organism," \\(canopy\\)")) %>%
  # manually filter out some non-species
  filter(!organism %in% c(
    "Shell hash",
    "Bare Rock",
    "Red lichen",
    "Black Zone Spp",
    "Brown Lichen",
    "Coralline (crust)",
    "Little Round Green Things",
    "Algal parasite on Mastocarpus"
  ),
  !str_detect(organism,"Unknown"),
  !str_detect(organism, "\\(dead\\)"),
  !str_detect(organism, "\\(eggs\\)"),
  !organism == "Mastocarpus stellatus (crust)",
  ) %>%
  
  # Bonnemaisonia hamifera as normal and as "(Trailiella phase)"
  # but there are only 9 occurrences total so just change them into one
  # also, two names are invalid - change those here. 
  mutate(organism = recode(organism,
                           "Bonnemaisonia hamifera (Trailiella phase)" = "Bonnemaisonia hamifera",
                           "Leathesia difformis" = "Leathesia marina",
                           "Prasiola crispa" = "Prasiola stipitata")) %>%
  
  # filter out single name organisms again because some reappeared after removing "canopy"
  filter(str_detect(organism, " ")) %>%
  
  # filter to data_taken == yes
  filter(data_taken == "yes") %>%
  
  # filter out replicates with "a" and "b",
  # and change "NA" replicates to replicate 1
  # (these cases there was only one replicate)
  mutate(replicate = case_when(is.na(replicate) ~ "1",
                               TRUE ~ replicate)) %>%
  filter(!str_detect(replicate,"a|b")) %>% # again note this is only one instance
  mutate(replicate = as.numeric(replicate)) %>%
  # filter to only whole-number levels
  filter(level %% 1 == 0)

rm(cover)


## | make percent cover numeric ------------------------------------------
cover3 <- cover2 %>% 
  mutate(pc_num = as.numeric(percent_cover)) %>%
  # there are a few character values that weren't transcribed properly.
  # fix them manually here
  mutate(pc_num = case_when(
    percent_cover == "0. 5" ~ 0.5,
    percent_cover == "10% (est'd canopy)" ~ 10,
    percent_cover == "15 (80 sub canopy)" ~ 15,
    percent_cover == "15% (est'd canopy)" ~ 15,
    percent_cover == "25% (est'd canopy)" ~ 25,
    percent_cover == "upper 15.6%/substrate 1.2%" ~ 15.6,
    TRUE ~ pc_num
  )) %>%
  
  # add p/a column
  mutate(pres = case_when(
    pc_num > 0 ~ TRUE,
    pc_num == 0 ~ FALSE,
    percent_cover %in% c(
      "p",
      "1patch",
      "1stipe",
      "2patches",
      "6stipes",
      "<1"
    ) ~ TRUE,
    percent_cover %in% c("nd","ND") ~ NA,
    is.na(percent_cover) ~ NA,
    TRUE ~ NA
  ))

rm(cover2)

cover3 %>% glimpse()

# here, we have organisms showing up as repeats because a few steps ago,
# we changed Bonnemaisonia hamifera (Trailiella phase) to Bonnemaisonia hamifera,
# so many appear repeated now. group them by year/level/transect/level/org
# and take the top row, arranged by presence, so when they are present,
# that will be the row kept. (most cases they are absent for both options)
cover4 <- cover3 %>%
  group_by(year, transect, level, replicate, organism) %>%
  arrange(-pres) %>%
  # slice top row (the present row if there is one)
  slice(1) %>% ungroup() %>%
  
  # filter to only species that are 
  # present at least once
  group_by(organism) %>%
  mutate(max_pres = sum(pres, na.rm = T)) %>%
  ungroup() %>%
  filter(max_pres > 0 ) %>%
  select(-max_pres)

# here we've filtered out 3 species that appear never
# to be present: 1 Celleporella hyalina, 2 Spermothamnion repens, 3 Pilayella littoralis   
rm(cover3)



## | make a pres abs dataset -------------------------------------------------
cover_pa <- cover4 %>% 
  select(-pc_num,
         -percent_cover)



# 5. merge all P/A --------------------------------------------

## | find data taken quadrats --------------------------------------------
# first, find sites where data was taken in both years
# we will later use this to know when richness was not counted 
# (e.g., data not taken), or when data was taken and richness was zero
# and also to filter to only quadrats where both count and cover data
# were collected. 
count_taken <- counts_pa %>%
  distinct(year, transect, level, replicate) %>%
  mutate(count_taken = "yes")

cover_taken <- cover_pa %>%
  distinct(year, transect, level, replicate) %>%
  mutate(cover_taken = "yes")

both_taken <- count_taken %>%
  full_join(cover_taken) %>%
  mutate(both_taken = case_when(
    count_taken == "yes" & cover_taken == "yes" ~ "both",
    count_taken == "yes" & is.na(cover_taken) ~ "count only",
    is.na(count_taken) & cover_taken == "yes" ~ "cover only"
  )) %>%
  filter(both_taken == "both") %>%
  select(-count_taken,
         -cover_taken)

rm(count_taken,
   cover_taken)


## | save data taken ---------------------------------------------------------
readr::write_csv(
  both_taken,
  here::here("data-processed",
             "appledore-survey-data",
             "quadrats_data_taken.csv")
  
)

cover_pa %>% glimpse()
counts_pa %>% glimpse()

## | merge all PA --------------------------------------------------------
all_pa <- 
  cover_pa %>%
  mutate(df = "cover") %>%
  rbind(counts_pa %>% mutate(df = "count"))
rm(counts_pa,
   cover_pa)



## | get P only --------------------------------------------------------------
# make presence only dataset, since we don't really need absences
# as long as data are all taken. 
p_only <- all_pa %>%
  # filter to pres only
  filter(pres == TRUE) %>%
  # pivot to separate count from cover pres
  # (also this will fix species that are in both datasets)
  tidyr::pivot_wider(names_from = df,
                     names_prefix = "pres_",
                     values_from = pres,
                     values_fill = FALSE) %>%
  # add presence column when spp is pres in either dataset
  mutate(pres = case_when(pres_cover == TRUE | pres_count == TRUE ~ TRUE)) %>%
  # join to the "both_taken" dataframe to keep only quadrats
  # where both count and cover were taken
  # removes 157 rows of about 37,500
  inner_join(both_taken)
rm(all_pa)

# find all the organism names to make sure
# they're in the species dataset
p_orgs <- p_only %>%
  distinct(organism)
length(p_orgs$organism)
# 93 species total

(p_orgs$organism %in% spp2$valid_name) %>% table()
# all of them are in the species list valid names column

# save all species list ---------------------------------------------------
spp_pres <- spp2 %>%
  filter(valid_name %in% p_orgs$organism) %>%
  distinct(AphiaID, valid_name) %>%
  
  # note that one species (Idotea balthica) is present twice,
  # with 2 aphiaIDs. the correct one is 119039, so we'll slice out the other
  group_by(valid_name) %>%
  #filter(valid_name == "Idotea balthica") %>%
  arrange(AphiaID) %>%
  slice(1) %>%
  ungroup() %>%
  rename(organism = valid_name)

readr::write_csv(spp_pres,
                 here::here("data-processed",
                            "appledore-survey-data",
                            "pres_spp_list.csv"))

rm(p_orgs, spp_pres)



## | save P only -------------------------------------------------------------
readr::write_csv(p_only,
                 file = here::here("data-processed",
                                   "appledore-survey-data",
                                   "spp_pres_by_replicate.csv"))

p_only %>% filter(organism == "Crisia eburnea" )

# 6. return to count data ----------------------------------------------------
counts_sep <- counts4 %>%
  
  # keep only organisms that have a numeric count value > 0
  group_by(organism) %>%
  mutate(max_count = max(countnum, na.rm = T)) %>%
  filter(max_count > 0) %>% 
  # note, this filters out "Crisia eburnea", but it stays in the p_only
  # dataset twice becuase it has two instances where count == "p"
  select(-max_count) %>%
  # split into list by organism
  split(f = .$organism)

rm(counts4)
counts_sep[[7]] %>%
  # first, add a dummy column to signify this was in the dataset
  mutate(in_sample = "yes") %>%
  
  # join with the both_taken df using full_join to keep all rows
  full_join(both_taken,
            by = join_by(year, transect, level,
                         replicate)) %>%
  distinct(data_taken, both_taken) 

## | split and filter --------------------------------------------------------
# here, we'll split by species, filter each one to only transects and
# levels where it is present, then recombine.
counts_sep_filt <- purrr::map(
  .x = counts_sep,
  .f = ~.x %>%
    
    # first, add a dummy column to signify this was in the dataset
    mutate(in_sample = "yes") %>%
    
    # join with the both_taken df using full_join to keep all rows
    full_join(both_taken,
              by = join_by(year, transect, level,
                           replicate)) %>%
    
    # filter to only quadrats where both data were taken
    filter(both_taken == "both") %>%
    
    # this leaves some rows where data were taken, but the species wasn't seen.
    # so, where in_df == NA, we're going to fill in count_num as zero
    # because this means data WERE taken, but the species wasn't there
    mutate(countnum = 
             case_when(is.na(in_sample) & both_taken == "both" ~ 0,
                       TRUE ~ countnum),
           organism = unique(.$organism, na.rm = T)[1],
           data_taken = "yes") %>%
    
    tidyr::replace_na(replace = list(in_sample = "zero_filled")) %>%
    
    # note - it's ok if there are still some NAs - those are cases
    # where the sample was "p" or similar, but the zeros are filled in
    # filter(is.na(countnum))
    
    # Filter to only levels where the species is found at least once
    # identify min and max bounds of species
    mutate(min_level = min(level[countnum > 0], na.rm = T),
           max_level = max(level[countnum > 0], na.rm = T)) %>%
    # filter to only levels between min and max 
    filter(level >= min_level & level <= max_level) %>%
    select(-min_level, -max_level)   %>% 
    
    ## filter species to only transects where they are found
    group_by(transect) %>%
    mutate(max_in_transect = max(countnum, na.rm=T),) %>%
    filter(max_in_transect > 0) %>%
    dplyr::select(-max_in_transect) %>%
    ungroup()
) %>%
  # bind list back together
  bind_rows()

counts_sep_filt %>% count(in_sample)
# note that we backfilled zeros in 1727 rows of 148943, 
# so ~1% of the data is assumed zeros


# split the counts for a few species just to check them out
i <- 2
counts_sep_filt %>%
  filter(organism == unique(counts_sep_filt$organism)[i]) %>%
  filter(level <= 15,
         level >= 0) %>%
  group_by(year,level, transect) %>% 
  summarize(countnum = mean(countnum)) %>%

  group_by(year, level) %>%
  summarize(countnum = mean(countnum)) %>%
#
  ggplot(aes(x = year, y = countnum + 1)) +
  geom_point(aes(color = countnum > 0)) +
  geom_smooth(method = "lm") +
  facet_wrap(~level) +
  labs(title = unique(counts_sep_filt$organism)[i]) +
  scale_y_log10()

rm(counts_sep)


## | save counts filtered ------------------------------------------------
readr::write_csv(counts_sep_filt,
                 here::here("data-processed",
                            "appledore-survey-data",
                            "counts_abundance_filtered.csv"))
# note - 36 species total
#rm(counts_sep_filt)

# 7. Return to cover --------------------------------------------------------
cover_sep <- cover4 %>%
  # keep only organisms that have a numeric count value > 0
  group_by(organism) %>%
  mutate(max_pc = max(pc_num, na.rm = T)) %>% # here we actually don't lose any rows
  filter(max_pc > 0) %>%
  select(-max_pc) %>%
  # split into list
  split(f = .$organism)
rm(cover4)
cover_sep[[1]] %>% glimpse()


## | split and filter --------------------------------------------------------
cover_sep_filt <- purrr::map(
  .x = cover_sep,
  .f = ~.x %>%
    
    # first, add a dummy column to signify this was in the dataset
    mutate(in_sample = "yes") %>%
    
    # join with the both_taken df using full_join to keep all rows
    full_join(both_taken,
              by = join_by(year, transect, level,
                           replicate)) %>%
    
    # filter to only quadrats where both data were taken
    filter(both_taken == "both") %>%
    
    # this leaves some rows where data were taken, but the species wasn't seen.
    # so, where in_df == NA, we're going to fill in percent cover as zero
    # because this means data WERE taken, but the species isn't there
    mutate(pc_num = 
             case_when(is.na(in_sample) & both_taken == "both" ~ 0,
                       TRUE ~ pc_num),
           organism = unique(.$organism, na.rm = T)[1],
           data_taken = "yes") %>%
    
    tidyr::replace_na(replace = list(in_sample = "zero_filled")) %>%
    
    # it's ok if there are still some NAs - those are cases
    # where the sample was "p" or similar, but the zeros are filled in
    # filter(is.na(pc_num))
    
    # Filter to only levels where the species is found at least once
    # identify min and max bounds of species
    mutate(min_level = min(level[pc_num > 0], na.rm = T),
           max_level = max(level[pc_num > 0], na.rm = T)) %>%
    # filter to only levels between min and max 
    filter(level >= min_level & level <= max_level) %>%
    select(-min_level, -max_level)   %>% 
    
    ## filter species to only transects where they are found
    group_by(transect) %>%
    mutate(max_in_transect = max(pc_num, na.rm=T),) %>%
    filter(max_in_transect > 0) %>%
    dplyr::select(-max_in_transect) %>%
    ungroup()
) %>%
  bind_rows()

cover_sep_filt %>% count(in_sample)
rm(cover_sep)
# here, we backfilled zeros in 19625 of 221667 rows, 
# or, about 9% of the data

# check some out again
i <- 4
cover_sep_filt %>%
  filter(organism == unique(cover_sep_filt$organism)[i]) %>%
  filter(level <= 15,
         level >= 0) %>%
  group_by(year,level, transect) %>% 
  summarize(pc_num = mean(pc_num)) %>%
  
  group_by(year, level) %>%
  summarize(pc_num = mean(pc_num)) %>%
  #
  ggplot(aes(x = year, y = pc_num + 1)) +
  geom_point(aes(color = pc_num > 0)) +
  geom_smooth(method = "lm") +
  facet_wrap(~level) +
  labs(title = unique(cover_sep_filt$organism)[i]) 


## | save cover filtered -----------------------------------------------------
readr::write_csv(cover_sep_filt,
                 here::here("data-processed",
                            "appledore-survey-data",
                            "cover_abundance_filtered.csv"))
cover_sep_filt %>% distinct(organism)
# 65 species total

  
# find total number of species (bc some are in both datasets)
abundance_spp <- unique(c(unique(counts_sep_filt$organism),
                          unique(cover_sep_filt$organism)))
length(abundance_spp) # 92 species total across both datasets
length(unique(p_only$organism)) # 93 species in p_only dataset 
unique(p_only$organism)[!unique(p_only$organism)%in% abundance_spp]
# "Crisia eburnea" is the one in the p/a dataset but not in the abundance

rm(list = ls())
