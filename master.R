
# master script for thermal tolerance project -----------------------------


# first, clean the datasets
source(here::here("scripts",
                  "01-clean-data.R"))
# this creates 4 outputs:
# 1. data-processed/quadrats_data_taken.csv
#    dataframe of the quadrats sampled in both datasets
# 2. data-processed/spp_pres_by_replicate.csv
#    dataset of all species present in all quadrats
#    (merged count and cover species - 93 species total)
# 3. data-processed/counts_abundance_filtered.csv
#    abundance data for all count species, filtered to only levels and 
#    transects where they are present at least once
# 4. data-processed/cover_abundance_filtered.csv
#    abundance data for all cover species, filtered to only levels and 
#    transects where they are present at least once
# 5. data-processed/spp_pres_list.csv
#    list of the Aphia IDs and valid names of all species ever present
#    in both datasets (93 species total)


# second, find global sst to match to occurrence records
source(here::here("scripts",
                  "02-extract-global-sst-1982-2023.R"))
# this script uses the rerddap package to collect the 
# Hadi SST monthly mean temperatures at 1x1 grid cells
# from 1981 to 2023. These are then summarized to find the 
# hottest monthly temp, coldest monthly temp, and mean monthly temp
# for every grid cell in each year.
# outputs:
# 1. data-processed/global_temps_1982_2023.rds
#    this is a list of temperature values in each 1x1 grid cells, 
#    with each element of the list being one year. 



# third, find global occurrences for every species
source(here::here("scripts",
                  "03-metch-occurrences-to-global-temp.R"))
# this script collects OBIS records for every species in the dataset
# filters to occurrences after 1981, with depth values <10m.
# then, matches every geo-located occurrence record to the temperature
# values from that grid cell in that year, which we collected above.
# outputs:
# 1. data-processed/obis_recs_matched_temp.rds
#    list of every obis record (both temp-matched and not) per species,
#    with temperatures matched when possible
# 2. data-processed/number_obis_recs_per_spp.csv
#    a df keeping track of the number of records obtained originally,
#    and at every step of filtering. 




# fourth, calculate and plot Species Thermal Affinities
source(here::here("scripts",
                  "04-calculate-thermal-affinities.R"))
# this script summarizes the OBIS temp-matched occurrences from above
# to find the summary values of temperature at all occurrences.
# outputs:
# 1. data-processed/spp_thermal_affinities.csv
#    dataframe of all thermal affinities of 89 species, with 12 
#    summary values describing the temperature of occurrences
# 2. outputs/Fig1_thermal_affinities.png
#    plot of all 89 species' thermal affinities, showing 
#    the mean monthly mean temperature of all occurrences. 





