
# master script for thermal tolerance project -----------------------------



# clean raw survey data ---------------------------------------------------
# first, clean the datasets
source(here::here("scripts",
                  "01-clean-raw-data",
                  "01-clean-data.R"))
# this creates 4 outputs:
# all in data-processedd/appledore-survey-data/
# 1. quadrats_data_taken.csv
#    dataframe of the quadrats sampled in both datasets
# 2. spp_pres_by_replicate.csv
#    dataset of all species present in all quadrats
#    (merged count and cover species - 93 species total)
# 3. counts_abundance_filtered.csv
#    abundance data for all count species, filtered to only levels and 
#    transects where they are present at least once
# 4. cover_abundance_filtered.csv
#    abundance data for all cover species, filtered to only levels and 
#    transects where they are present at least once
# 5. spp_pres_list.csv
#    list of the Aphia IDs and valid names of all species ever present
#    in both datasets (93 species total)

# get appledore shapefile for mapping
source(here::here("scripts",
                  "01-clean-raw-data",
                  "01.2-get-appledore-shape.R"))
# outputs 
# 1. data-processed/appledore-island-env-data/appledore_shp.rds
#    shapefile of appledore island


# get HadiSST data --------------------------------------------------------
# second, find global sst to match to occurrence records
source(here::here("scripts",
                  "02-get-sst-data",
                  "02-extract-global-sst-1982-2023.R"))
# this script uses the rerddap package to collect the 
# Hadi SST monthly mean temperatures at 1x1 grid cells
# from 1981 to 2023. These are then summarized to find the 
# hottest monthly temp, coldest monthly temp, and mean monthly temp
# for every grid cell in each year.
# outputs:
# in data-processed/global-temp-data/
# 1. global_temps_1982_2023.rds
#    this is a list of temperature values in each 1x1 grid cells, 
#    with each element of the list being one year. 

# find appledore island sst 
source(here::here("scripts",
                  "02-get-sst-data",
                  "02.2-extract-appledore-sst-1982-2023.R"))
# data-processed/appledore-island-env-data/appledore_temps_1982_2023.csv
# yearly temperature values of appledore island throughouth the study period



# calculate species thermal affinities ------------------------------------
# third, find global occurrences for every species
source(here::here("scripts",
                  "03-calculate-STIs",
                  "03-match-occurrences-to-global-temp.R"))
# this script collects OBIS records for every species in the dataset
# filters to occurrences after 1981, with depth values <10m.
# then, matches every geo-located occurrence record to the temperature
# values from that grid cell in that year, which we collected above.
# outputs:
# 1. data-processed/species-thermal-affinities/obis_recs_matched_temp.rds
#    list of every obis record (both temp-matched and not) per species,
#    with temperatures matched when possible
# 2. data-processed/species-thermal-affinities/number_obis_recs_per_spp.csv
#    a df keeping track of the number of records obtained originally,
#    and at every step of filtering. 


# calculate and plot Species Thermal Affinities
source(here::here("scripts",
                  "03-calculate-STIs",
                  "03.2-calculate-thermal-affinities.R"))
# this script summarizes the OBIS temp-matched occurrences from above
# to find the summary values of temperature at all occurrences.
# outputs
# 1. data-processed/species-thermal-affinities/spp_thermal_affinities.csv
#    dataframe of all thermal affinities of 89 species, with 12 
#    summary values describing the temperature of occurrences
# 2. outputs/Fig1_thermal_affinities.png
#    plot of all 89 species' thermal affinities, showing 
#    the mean monthly mean temperature of all occurrences. 



# make highly sampled dataset ---------------------------------------------
# some analyses in this project are best suited to data with standardized
# sampling across transects and depths. Here, we calculate a highly sampled 
# dataset, which keeps only transects sampled in over 20 years, tidal heights
# that are consistently sampled through time (9 levels), and years within 
# these bounds where 6+ of 9 levels were sampled, and years when at least
# 3 of 7 transects fit these criteria. we retain 7 transects and 37 years 
# of sampling
source(here::here("scripts",
                  "04-make-highly-sampled-df.R"))
# outputs, in data-processed/appledore-island-survey-data/pa-with-therm
# 1. pa-with-therm-all.csv
#    presence absence data across all transects, years, levels, joined 
#    with thermal tolerance value for all species
# 2. pa-with-therm-highly-sampled.csv
#    same as above, but for only highly sampled quadrats


# find species distribution shifts ----------------------------------------

# first, try to assess depth shifts of individual species
source(here::here("scripts",
                  "05-spp-distribution-shifts",
                  "05-find-height-changes-per-org.R"))
# outputs:
# 1. outputs/depth_shifts/all_data_trends_p.png
#    boxplots of species distribtuions across tidal heights over 5-year bins
# 2. outputs/depth_shifts/all_data_bars.png
#    trends of sign and significance of models predicting depth shifts over time
#    arranged by distributional position (upper edge, lower edge, mean, median)
# 3-4. same as above, but with highly sampled dataset

# second, find species that appear or disappear throughout the sampling span
source(here::here("scripts",
                  "05-spp-distribution-shifts",
                  "05.2-find-spp-appear-disappear.R"))
# outputs
# 1. outputs/depth_shifts/spp_appear_disappear.png
#    multipanel figure showing species that appear and disappear throughout 
#    the sampling span, and thermal affinities (mean, min, max) of each group.



# model abundance change --------------------------------------------------
# prep abundance data for modeling
source(here::here("scripts",
                  "06-abindance-modeling",
                  "06.1-prep-count-data-for-model.R"))
# this script cleans the abundance data for modeling.
# basically, here we summarize data to mean abundance 
# per tidal height per year, using all data available
# we remove species that are present in <3 years, remove 
# tidal height singletons per species, and keep only tidal levels
# sampled in over 20 years of sampling. also reformat to density per m2
# outputs:
# 1. data-processed/abundance-data/count-data-prepped-for-model.csv
#    data prepped for adundance change model

# model density species trends over time
source(here::here("scripts",
                  "06-abindance-modeling",
                  "06.3-model-abundance-counts.R"))
# outputs:
# 1. data-processed/abundance-data/abundance_change_slopes_counts_cattidalheight.rds
#    tibble of all model dfs, abundance model outputs, predicted values, 
#    confidence intervals, etc., and also all thermal affinity values



# redo with tweedie distributions as a sensitivity test
source(here::here("scripts",
                  "06-abindance-modeling",
                  "06.4-mdoel-tweedie-counts.R"))
# outputs:
# 1. data-processed/abundance-data/abundance_change_slopes_counts_tweedie_cattidalheight.rds"
#    tibble of all model dfs, abundance model outputs, predicted values, 
#    confidence intervals, etc., and also all thermal affinity values


# clean cover data for modeling
source(here::here("scripts",
                  "06-abindance-modeling",
                  "06.5-prep-cover-data-for-model.R"))
# this script cleans the abundance data for modeling.
# basically, here we summarize data to mean abundance 
# per tidal height per year, using all data available
# we remove species that are present in <3 years, remove 
# tidal height singletons per species, and keep only tidal levels
# sampled in over 20 years of sampling. 
# outputs:
# 1. data-processed/abundance-data/cover-data-prepped-for-model.csv
#    data prepped for abundance change model

# model percent cover species trends over time
source(here::here("scripts",
                  "06-abindance-modeling",
                  "06.7-model-abundance-cover.R"))
# outputs:
# 1. data-processed/abundance-data/abundance_change_slopes_cover_cattidalheight.rds
#    tibble of all model dfs, abundance model outputs, predicted values, 
#    confidence intervals, etc., and also all thermal affinity values


# model abundance coefficients by thermal affinity
source(here::here("scripts",
                  "06-abindance-modeling",
                  "06.8-model-slopes-by-therm.R"))
# outputs:
# 1. outputs/abundance-change/Fig3_abundance_change_slopes_by_temp.png
#    model output of slopes by thermal affinity for counts and cover
# 2. outputs/abundance-change/Fig3_abundance_change_slopes_by_temp_SHARED_SPP.png
#    model output of slopes by thermal affinity, highlighting shared species
# 3. outputs/abundance-change/counts_gamma_tweedie_compare.png
#    model output of density species with gamma and tweedie distributions
#    as sensitivity test for treatments of zero data
# 4. outputs/abundance-change/counts_gamma_tweedie_compare_2.png
#    as above, but in one plot
# 5. outputs/abundance-change/abundance_by_t_all_metrics.png
#    thermal affinity metrics using mean, min, max thermal affinity




