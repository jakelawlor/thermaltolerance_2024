# Script to find species that are shared between datasets
# create one master list of all species with where they are found


library(dplyr)
library(here)
library(stringr)
library(tidyr)


# Step 1. Upload Data 
#--------------------------------------------------------------
int_counts <- read.csv(here::here(
  "data-raw",
  "EDI Data thru 2023",
  "counts_data.csv"
)) %>% janitor::clean_names()
int_counts %>% glimpse()

int_cover <- read.csv(here::here(
  "data-raw",
  "EDI Data thru 2023",
  "percent_cover_data.csv"
)) %>% janitor::clean_names()
int_cover %>% glimpse()

spp_list <- read.csv(here::here(
  "data-raw",
  "EDI Data thru 2023",
  "species_list.csv"
)) %>% janitor::clean_names()
spp_list %>% glimpse()



# filter to only species names --------------------------------------------
# many observations are either vague ("arthropod") or non-species (bare rock)
# first, filter these out to cut to only data containing species
int_cover_org <- int_cover %>%
  # filter out words without spaces ("Arthropod")
  filter(str_detect(organism, " ")) %>%
  # filter out eggs
  filter(!str_detect(organism, "(egg)")) %>%
  # filter out specific non-species that are in the database
  filter(!organism %in% c("Bare Rock",
                          "Shell hash",
                          "Little Round Green Things",
                          "Black Zone Spp",
                          "Red lichen",
                          "Brown Lichen",
                          "Algal parasite on Mastocarpus"),
         !str_detect(organism, "Unknown"),
         !str_detect(organism, "dead")) %>%
  # for those measured as both (canopy) and (primary),
  # keep only canopy, as the data is recorded longer
  # additionally Mastocarpus is measured with crust. remove that too.
  # for Bonnermaisonia hamifera, most counts are (Trailela phase). Keep only those. 
  filter(!str_detect(organism,"primary"),
         !str_detect(organism, "crust"),
         organism != "Bonnemaisonia hamifera") %>%
  mutate(organism = str_replace(organism," \\(canopy\\)",""),
         organism = str_replace(organism, " \\(Trailiella phase\\)","")) %>%
  # and again, filter to organisms with space ("Fucus" is left after above filtering)
  filter(str_detect(organism, " ")) 


# view differences in primary/canopy species
int_cover %>%
  filter(str_detect(organism, "canopy|primary")) %>%
  mutate(percent_cover = as.numeric(percent_cover)) %>%
  mutate(cat = case_when(str_detect(organism,"canopy") ~ "canopy",
                         TRUE ~ "primary")) %>%
  mutate(organism = str_replace(organism,"\\(canopy\\)|\\(primary\\)","")) %>%
  ggplot(aes(x = year, y = percent_cover,
             color = cat)) +
  geom_point(alpha = .3,
             size = .15,
             position = position_jitter(width = .3,
                                        height = 0)) +
  geom_smooth(method = "lm") +
  facet_wrap(~organism) +
  scale_y_sqrt()

# view differences in Bonnemaisonia types
int_cover %>%
  filter(str_detect(organism, "Bonnemaisonia")) %>%
  mutate(percent_cover = as.numeric(percent_cover)) %>%
  ggplot(aes(x = year, 
             y = percent_cover,
             color = organism)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~organism) +
  scale_y_sqrt()


# change percent cover to numeric -----------------------------------------
# and add presense absence column
int_cover_org %>%
  
  # manually change some abnormal percent cover entries
  mutate(percent_cover = case_when(
    str_detect(percent_cover, "est'd") ~ str_remove(percent_cover, "% \\(est'd canopy\\)"),
    str_detect(percent_cover,"sub canopy") ~ str_remove(percent_cover," \\(80 sub canopy\\)"),
    percent_cover == "upper 15.6%/substrate 1.2%" ~ "15.6",
    str_detect(percent_cover,"patch|stipe") ~ "p",
    TRUE ~ percent_cover
  )) %>%
  
  mutate(percent_cover_num = as.numeric(percent_cover)) %>%
  mutate(pres = case_when(percent_cover_num > 0 ~ TRUE,
                          stringr::str_detect(percent_cover,"[1-9]") ~ TRUE,
                          percent_cover == "p" ~ TRUE,
                          percent_cover_num == "0" ~ FALSE,
                          TRUE ~ NA)) %>%

  filter(is.na(pres)) %>%
  distinct( percent_cover)

int_cover_org %>%
  mutate(percent_cover_num = as.numeric(percent_cover)) %>%
  filter(is.na(percent_cover_num)) %>%
  distinct(organism, percent_cover, percent_cover_num)

  stringr::str_extract_all("upper 15.6%/substrate 1.2%", "[0-9]*")


# find presence absence ---------------------------------------------------
int_counts_org_pres <- int_counts_org %>%
  # change count to numeric - introduces some NAs
  mutate(count_num = as.numeric(count)) %>%
  # make presence column when count_num > 0 or 
  # count is "p" or "sp100"
  mutate(pres = case_when(count_num > 0 | count == "p" | count == "sp100" ~ TRUE,
                          TRUE ~ FALSE)) %>%
  filter(pres == TRUE) %>%
  select(-count, -count_num) %>%
  mutate(dataset = "count")
int_counts_org_pres %>% glimpse()

int_cover_org_pres <- int_cover_org %>%
  # change percent cover to numeric - introduces NAs
  mutate(percent_cover_num = as.numeric(percent_cover)) %>%
  # make present columnn for when percent_cover_num > 0,
  # percent_cover == "p". or percent_cover has any number in it e.g., "15% (est'd canopy)"
  mutate(pres = case_when(percent_cover_num > 0 ~ TRUE,
                          str_detect(percent_cover,"[1-9]") ~ TRUE,
                          percent_cover == "p" ~ TRUE,
                          TRUE ~ FALSE)) %>%
  filter(pres == TRUE) %>%
  select(-percent_cover, -percent_cover_num) %>%
  mutate(dataset = "cover")
  
pres_full <- int_counts_org_pres %>%
  rbind(int_cover_org_pres) %>%
  mutate(dup = duplicated(.)) %>%
  filter(dup == TRUE)
  duplicated()
  
  int_counts %>%
    group_by(year, level, transect, replicate, data_taken, organism) %>%
    add_count()%>%
    filter(n > 1) %>% 
    
    arrange(desc(n))
  
  int_counts %>%
    filter(year == 1982, 
           level == 1, 
           transect == 5,
           replicate == 1, 
           str_detect(organism, "Littorina"))
    mutate(dup = duplicated(year, transect, level, replica))
    filter(year == 2022,
           transect == 22,
           level %in% c(6:7),
           replicate == 1,
           organism == "Semibalanus balanoides")

# clean counts data -------------------------------------------------------
#remove NA years 
int_counts <- int_counts %>%
  filter(!is.na(year))

# find all transects / levels with data
int_counts_data_taken <- int_counts %>%
  filter(data_taken == "yes" |
           (!is.na(count) & count != 0)) %>%
  distinct(year, transect, level, replicate) %>%
  mutate(count_taken = "yes")

int_counts_data_taken %>%
  ggplot(aes(x = year, 
             y = level)) +
  geom_tile() +
  facet_wrap(~transect)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# find instances where nothing was found but data was still taken
# e.g. true zeros
true_zeros_cover <- int_cover %>%
  mutate(cover2 = as.numeric(percent_cover)) %>%
  group_by(year, transect, level,data_taken) %>%
  summarize(sum = sum(cover2, na.rm=  T),
            nrow = n()) %>%
  ungroup() %>%
  arrange(sum) %>%
  filter(sum == 0) %>%
  filter(data_taken== "yes")


true_zeros_counts <- int_counts %>%
  mutate(count2 = as.numeric(count)) %>%
  group_by(year, transect, level,data_taken) %>%
  summarize(sum = sum(count2, na.rm=  T),
            nrow = n()) %>%
  ungroup() %>%
  arrange(sum) %>%
  filter(data_taken== "yes")



int_cover %>%
  filter(year == 1982,
         transect == 1,
         level %in% c(5:10))



# clean cover data --------------------------------------------------------
# find all transects / levels with data
int_cover_data_taken <- int_cover %>%
  filter(data_taken == "yes" |
           (!is.na(percent_cover) & percent_cover != 0)) %>%
  distinct(year, transect, level, replicate) %>%
  mutate(cover_taken = "yes")

int_cover_data_taken %>%
  ggplot(aes(x = year, 
             y = level)) +
  geom_tile() +
  facet_wrap(~transect)

data_taken_all <- int_counts_data_taken %>%
  full_join(int_cover_data_taken) %>%
  mutate(data_taken = case_when(
    count_taken == "yes" & cover_taken == "yes" ~ "Both",
    count_taken == "yes" & is.na(cover_taken) ~ "Count only",
    cover_taken == "yes" & is.na(count_taken) ~ "Cover only"
  )) 
  
data_taken_all %>% 
  ggplot(aes(x = year,
             y = level,
             fill = data_taken)) +
  geom_raster(alpha = .7) +
  facet_wrap(~transect) +
  labs(fill = "Data Taken")

data_taken_all %>%
  filter(data_taken == "Cover only") %>%
  
  filter(#year == 2017, 
         transect == 7) %>%
  distinct(year, level) %>% 
  arrange(year)


int_cover %>%
  filter(year == 2016, 
         transect == 7,
         level %in% c(0:5))

# note that some transects are sampled considerably less than others.
# those with longest timeseries are:
# 5, 7, 15, 20, 22, 26, 28
# consider making a dataset with just these transects 
# to use primarily, or as a sensitivity test




# find data taken ---------------------------------------------------------

range(int_counts$year, na.rm = T)

int_counts %>%
  filter(data_taken == "no" & (!is.na(count) & count != 0 ))

int_counts %>%
  group_by(year, transect, level, replicate) %>%
  filter((data_taken == "no" | is.na(data_taken)) & !is.na(count) & count != 0) %>%
  ggplot(aes(x = year, 
             y = level)) +
  geom_tile() +
  facet_wrap(~transect)


int_cover %>% 
  group_by(year, transect, level, replicate) %>%
  count(data_taken)
filter(data_taken == "yes" | 
         (is.na(data_taken) & !is.na(count))) %>%
  distinct(year, transect, level, replicate) %>%
  ungroup()


# 1. Get Int Counts Species --------------------------------------------------
int_spp %>% glimpse()

# see unique species in counts dataset
int_counts %>% distinct(organism) 
# 76 distinct organisms initially, 
# but many are too general ("Cancer", "Annelid")

# start with full data
#int_counts_spp <- 
int_counts %>% 
  
  # first, remove all times when data_taken == no
  filter(data_taken == "yes" | is.na(data_taken),
         !is.na(count)) %>%
  
  
  # rename `organism` column to `name` for consistency
  rename(name = organism) %>%
  
  # filter to only multi-word entries (i.e. filter out "amphipoda")
  filter(str_detect(name, " ")) %>%
  # this already lowers us to 40 names 
  
  # remove parenthetical outliers (e.g. canopy, cover, egg case)
  separate(col = name, 
           into =c("name","parenthetical"), 
           sep=" \\(",
           fill = "right") %>%
  mutate(parenthetical = gsub(")","",parenthetical))  %>% 
  
  #distinct(parenthetical)
  # in this case, the only parentheticals are "egg case", which we'll remove
  filter(is.na(parenthetical)) %>%
  
  # manually change 2 names that I know are multiple names for the same 2 species,
  # where both versions of the name are accepted and messing up counts later:
  mutate(name = recode(name,
                       # fix typo
                       "Idotea baltica" = "Idotea balthica",
                       # resolve species
                       "Tectura testudinalis" = "Testudinalia testudinalis",
                       "Tectura testinalis" = "Testudinalia testudinalis")) %>%
  
  # remove parenthetical column
  select(-parenthetical) %>%
  
  # currently, count is a character column, because sometimes it 
  # has "p" or "sp100". first, separate numeric values, where possible
  mutate(count_numeric = as.numeric(count)) %>%
  
  #filter(is.na(count_numeric)) %>%
  #distinct(count)
  # just a note that distinct non-numeric values for count are "p" and "sp100"
  
  # add presence / absence column binary
  mutate(present = case_when(count == "p" | 
                               count == "sp100" | 
                               count_numeric > 0  ~ 1,
                             TRUE ~ 0))  
count(data_taken)

# now pull only spp that were present at some time (i.e. remove all present = 0 rows)
filter(!present == 0) %>%
  
  # find number of years in which the species was observed
  group_by(year, name) %>%
  slice(1) %>%
  group_by(name) %>%
  summarize(years_present = sum(present)) %>%
  
  
  # add column for which dataset we're talking about 
  mutate(dataset = "int_counts") %>%
  arrange(desc(years_present)) %>%
  
  # filter out singletons only present in one year of dataset
  filter(years_present > 1)


int_counts_spp
# this leaves us with 27 species, which were shown as present
# from between 2 and 34 years
hist(int_counts_spp$years_present)
# pretty bimodal distribution of years present, 
# with most species either found in all years or very few
int_counts_spp %>% pull(name) %>% sort()

#### NOTE: 27 species after filtering after 1985 and removing spp only present in one year

# 2. Get Int Cover Species ---------------------------------------------------
int_cover_spp <- 
  int_cover %>%  
  
  # cut to data after 1985
  filter(year > 1985) %>%
  
  # rename organism column to name for consistency
  rename(name = organism) %>% 
  # initially 111 entries but many are vague ("Red Lichen"), or 
  # not organisms at all ("Feather","Brown ground")
  
  # remove one-word names
  filter(str_detect(name, " ")) %>% 
  # this cuts down to 88 entries
  
  # remove parenthetical outliers (e.g. canopy, cover, dead)
  separate(col = name, 
           into =c("name","parenthetical"), 
           sep=" \\(",
           fill = "right") %>% 
  mutate(parenthetical = gsub(")","",parenthetical)) %>%
  
  # filter to multi-word names again 
  # (since some passed the first filter due to parentheticals )
  filter(str_detect(name, " ")) %>% 
  
  # keep only parentheticals of "canopy" or nothing,
  # thus removing "dead", "crust","Trailiella phase", and "primary"
  filter(is.na(parenthetical) | parenthetical == "canopy") %>%
  
  # remove parenthetical specificity for now
  select(-parenthetical) %>% 
  
  # many percent_cover values have notes in them. 
  # need to separate those from the numbers so like "75 phymato" becomes 75
  drop_na(percent_cover) %>% 
  mutate(percent_cover_numeric = 
           as.numeric(
             str_extract(percent_cover, "\\-*\\d+\\.*\\d*")
           )
  ) %>%
  
  
  # add presence / absence column binary - when cover is > 0 or marked as "present". 
  mutate(present = case_when(percent_cover_numeric > 0  | 
                               percent_cover == "p" ~ 1,
                             TRUE ~ 0))  %>%
  
  
  # now pull only spp that were present at some time
  filter(present == 1) %>%
  
  # find years present
  arrange(desc(percent_cover_numeric)) %>%
  group_by(name, year) %>% # make one row per spp/year
  slice(1) %>%
  group_by(name) %>% 
  # now count number of years each species is present (either "p" or number)
  # and number of years each species is counted (numeric value > 0)
  summarize(years_present = sum(present),
            #years_counted = sum(percent_cover_numeric > 0, na.rm=T)
  ) %>%
  
  
  # add dataset name
  mutate(dataset = "int_cover") %>%
  arrange(desc(years_present)) %>%
  
  # filter out species present in only one year
  filter(years_present > 1)


# this gives us a dataset of 70 species ranging from one 
# observation to 38, but with many species not listed as 
# proper species names (e.g. "red lichen")
hist(int_cover_spp$years_present)
# many are only observed a few times

# 3. merge species lists from all datasets -----------------------------------
full_intertidal_spp_list <- 
  rbind(int_counts_spp,
        int_cover_spp
  ) %>% 
  group_by(dataset) %>%
  arrange(dataset,desc(years_present),name)

rm(int_counts, int_cover,
   int_counts_spp, int_cover_spp)



# 4. identify species not in the int_spp dataset -----------------------------------
nonvalid_species <- full_intertidal_spp_list %>%
  left_join(int_spp %>% select(name, valid_name, aphia_id)) %>%
  filter(is.na(valid_name))  %>% arrange(name)

nonvalid_species %>% select(name)
# most of these are nonvalid species names (e.g. Shell hash), 
# but some are just typos or misentries. let's fix them below:

## 4.1 manually validate some nonvalid species ---------------------------------
corrected_nonvalid_names <- nonvalid_species %>%
  ungroup() %>%
  mutate(valid_name = case_when(name == "Carcinus maenus" ~ "Carcinus maenas",
                                name == "Crisia eburna" ~ "Crisia eburnea",
                                name == "Dumontia concortum" ~ "Dumontia contorta",
                                name == "Giffordia granulosa" ~ "Hincksia granulosa",
                                name == "Nemalion helminthoides" ~ "Nemalion elminthoides",
                                name %in% c("Tectura testudinalis","Tectura testinalis") ~ "Testudinalia testudinalis",
                                name == "Tricellaria inopinata" ~ "Tricellaria inopinata",
                                TRUE ~ name)) %>%  
  mutate(aphia_id = as.numeric(aphia_id)) %>%
  mutate(aphia_id = case_when(valid_name == "Carcinus maenas" ~ 107381,
                              valid_name =="Crisia eburnea" ~ 111696,
                              valid_name == "Dumontia contorta" ~ 145228,
                              valid_name == "Hincksia granulosa" ~ 145433,
                              valid_name == "Nemalion elminthoides" ~ 145765,
                              valid_name == "Testudinalia testudinalis" ~ 234208,
                              valid_name == "Tricellaria inopinata" ~ 111254,
                              TRUE ~ aphia_id
  )) %>%
  filter(!is.na(aphia_id)) %>% select(name, valid_name, aphia_id)

# add those to the int_spp list
valid_species <- full_intertidal_spp_list %>%
  ungroup() %>%
  left_join(int_spp %>% 
              select(name, valid_name, aphia_id)%>%
              rbind(corrected_nonvalid_names)
  ) %>%
  filter(!is.na(valid_name)) %>%
  arrange(dataset,name)  

valid_species %>% 
  distinct(valid_name) 
# 91 valid species names overall, split to 94 in-text names
valid_species %>%
  filter(name != valid_name)

#### NOTE: only 79 after removing the singletons and filtering to >1985

## 4.2 double check WORMS IDs --------------------------------------------------

# some species have multiple listed WORMS id's. 
# let's check them in WORMS and correct the data frame
u<-worms::wormsbyid(unique(valid_species$aphia_id))
#recursively retrive information on the taxa they refer to
v<-worms::wormsconsolidate(u)
# what are the currently correct "accepted" taxa? Answer: "accepted_id".
w<-worms::wormsaccepted(v)
name_corrections <- w[,c("scientificname","AphiaID","status","valid_AphiaID","valid_name","accepted_id")] 

# check that there are no valid species not in the name corrections df
valid_species[!valid_species$valid_name %in% name_corrections$scientificname,]

rm(u,v,w, nonvalid_species, corrected_nonvalid_names)

valid_species_fixed <- valid_species %>% 
  rename(AphiaID = aphia_id) %>%
  left_join(name_corrections, by = c("AphiaID", "valid_name")) %>% 
  mutate(accepted_id = as.numeric(accepted_id)) %>% 
  mutate(AphiaID = case_when(status == "unaccepted" ~ accepted_id,
                             TRUE ~ AphiaID))  %>%
  select(name, valid_name, years_present, dataset, AphiaID) %>%
  rename(aphia_id = AphiaID) %>%
  arrange(dataset, desc(years_present))


valid_species_fixed %>% View()
# so there are always a few times name does not match valid name. 
# we will have to change these in all datasets
valid_species_fixed[valid_species_fixed$name != valid_species_fixed$valid_name,]

# write out full list ------------------------------------------------------
write.csv(valid_species_fixed,
          row.names = F,
          file = here("final_materials","data_processed","01_full_intertidal_spp_list.csv"))


rm(full_intertidal_spp_list,
   int_spp, name_corrections, valid_species)


valid_species_fixed %>% distinct(valid_name)

# make plot ============================================================
library(ggplot2)
p <- valid_species_fixed %>% 
  select(-years_present) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = dataset,
              values_from = present,
              values_fill = 0 ) %>%
  mutate(aphia_id = as.character(aphia_id)) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  pivot_longer(cols = contains("int"),names_to = "dataset",values_to = "present") %>%
  group_by(name) %>%
  arrange(desc(sum)) %>%
  filter(sum > 0) %>%
  mutate(label = case_when(dataset == "int_counts" ~ "Intertidal\nCounts",
                           dataset == 'int_cover' ~ "Intertidal\n% Cover")) %>%
  # plot it as a table
  ggplot(aes(x=label,y=reorder(name, desc(name)),fill=as.character(present))) +
  geom_tile(show.legend = F,color="grey20",size=.25) +
  scale_fill_manual(values = c("transparent","lightblue"))+
  labs(x=NULL,y=NULL,
       title = "Species present in intertidal datasets") + 
  ggthemes::theme_few() +
  scale_x_discrete(position = "top")+
  coord_cartesian(xlim=c(.5,3.5),expand = F,clip = "off") +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(family="Open Sans"),
        axis.text.x=element_text(family="Open Sans Semibold"),
        plot.title = element_text(family="Open Sans",hjust=0,margin=margin(b=5), size=18),
        plot.title.position = "plot"
  )+
  geom_text(aes(x=3, label=aphia_id), 
            stat="unique",
            family = "Open Sans",
            size=2) +
  theme(plot.margin = margin(t=5,r=10,b=5,l=10,unit = "mm")) 

p

#ggsave(filename = here("outputs","species_in_surveys.png"),
#       height = 10,
#       width = 9,
#       dpi=200)


