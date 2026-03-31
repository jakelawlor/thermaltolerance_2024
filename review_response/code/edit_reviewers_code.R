
# Response to Reviewer 1's analysis
# note that this analysis is conducted on separate, uncleaned data 
# that differ from the data which we used and provided. 
# still, we can use it to assess the comparisons that the reviewer
# makes to our values. 

# reviewer's script, as submitted: ----------------------------------------

library(readr)
library(vroom)
library(readxl)
library(tidyverse)
library(EDIutils)

## Download data directly from Environmental Data Initative
## Find the data package using Appledore as search term
search_results <- search_data_packages(query = 'q="Appledore"&fl=packageid,title')
search_results
## Get package ID
package_id <- search_results[1,1]
## Look at what is in the package
data_entities <- read_data_entity_names(packageId = package_id)
print(data_entities)
## Download count data, which is the second entity
raw <- read_data_entity(packageId = package_id, entityId = data_entities$entityId[2])
## Use vroom so you can check for problems
df <- vroom(file=raw)
df
unique(df$Year)
df.problems <- problems(df)
df.problems
p.condition <- as.list(unique(df.problems$actual))
p.condition


problems.summary <- NULL
for(condition in p.condition){
  tmp1 <- filter(df.problems,actual == condition)
  tmp2 <- slice_head(tmp1,n=1)
  problems.summary <- rbind(problems.summary,tmp2)
}
problems.summary
## Check what "1a","2a",etc. mean for rows 5716 to 5725
df.1a <- df %>% rownames_to_column()
df.1a[5716:5725,] ## "1a", etc. means NA


## A second check
problem.tmp <- filter(df.problems, actual == "1a")
problem.tmp
df[5716,] ## Looks like "NA" for "1a",etc
print(df[c(5716:5736,11817:11837),],n=42)
## It appears there are >4 replicates in some cases!!


## Find years and transects for which replicates are not numbered
problems.NA <- df %>%
  filter(is.na(Replicate)) %>%
  distinct(Year,Transect,Level,Data_taken) %>%
  arrange(Year,Transect,Level)
print(problems.NA,n=300)
## 1982, 1983, 1984 and 1992 are missing Replicate numbers
## 1982 NA in 7 transects, 1983 in 8, 198ti in 5 and 1992 in 1
df.transects <- df %>%
  filter(Year == 1982 | Year == 1983 | Year == 1984 | Year == 1992) %>%
  distinct(Year,Transect) %>%
  arrange(Year,Transect)
print(df.transects, n=63)
df.trans.summary <- df.transects %>%
  group_by(Year) %>%
  summarise(N = n())
df.trans.summary
## Proportion of missing transect:
## 1982 = 7/16, 1983 = 8/17, 198ti = 5/14, 1992 = 1/16

## To keep things simple, change NA to 5 to make it replicate 5
unique(df$Replicate)
df1 <- df %>%
  mutate(Replicate = replace_na(Replicate,5))
unique(df1$Replicate)

## OK: drop Count = "p", "sp100", and "casings present"
df2 <- df1 %>%
  filter(!Count == "p" | !Count == "sp100" | !Count == "casings present") %>%
  as_tibble()
df2

unique(df2$Data_taken)

## Drop Data_taken = "no" or NA
df3 <- df1 %>%
  filter(!Data_taken == "no") %>%
  filter(!is.na(Data_taken))
unique(df3$Data_taken)
df3 <- df3 %>% select(!Data_taken)

## Drop Organisms that are not species (i.e., egg cases)
unique(df3$Organism)
df4 <- df3 %>%
  filter(!str_detect(Organism, "egg case"))
df4 %>% filter(str_detect(Organism, "egg case")) ## OK
df4
unique(df4$Organism)

## Correct spelling of Tectura testundinalis vs Tectura testinalis
## Correct spelling of Idotea balthica vs Idotea baltica
df5 <- df4 %>%
  mutate(Organism = if_else(Organism == "Tectura testinalis","Tectura testudinalis",Organism),
         Organism = if_else(Organism == "Idotea baltica","Idotea balthica",Organism))
df5
unique(df5$Organism)

## Fix tide levels; Level is in feet and relative to benchmark pins
## Benchmark pins are set at 13.5 feet MLLW
## To convert Level to MLLW, subtract from 13.5, e.g., 13.5 - 3 = 10.5 MLLW
## Convert to meters
df5$Level <- 0.3408*(13.5 - df5$Level)
unique(df5$Level)

## Summaries
unique(df5$Year) ## Data from 1982 to 2017
length(unique(df5$Year)) ## Data covers diﬀerent 35 years
length(unique(df5$Organism)) ## 72 species
## Check the number of replicates
df5.replicates <- df5 %>%
  group_by(Year,Transect,Level,Organism) %>%
  summarise(N.rep = n()) %>%
  filter(Organism == "Littorina littorea") %>%
  filter(N.rep == 6)
df5.replicates
unique(df5.replicates$N.rep)

## The inclusion of Mya arenaria is very weird since it occurs in so< sediments
## Look at its occurence in the data
df5 %>% filter(Organism == "Mya arenaria" & Count >= 1)
## Mya occurs 5 times; only in 1982 and only in 3 transects

## Check for unusually large counts
df.max <- df5 %>%
  group_by(Organism) %>%
  summarise(max_value = max(Count, na.rm = TRUE))
print(df.max, n=76)
## Note max for Strongylocentrotus droebachiensis is 1082 in a 20 x 20 cm area
df5 %>% filter(Organism == "Strongylocentrotus droebachiensis" & Count == 1082)

## Look at distribution of counts (most are zeros)
df.summary.tmp <- df5 %>%
  group_by(Organism) %>%
  summarise(Zeros = sum(Count == 0,na.rm = TRUE),
            Ones = sum(Count == 1,na.rm = TRUE),
            Two.up = sum(Count >= 2,na.rm = TRUE),
            Eleven.up = sum(Count >= 11,na.rm = TRUE),
            Twenty_one.up = sum(Count >= 21,na.rm = TRUE),
            Fifty_one.up = sum(Count >= 51,na.rm = TRUE),
            More.than.100 = sum(Count >= 101,na.rm = TRUE),
            Max = max(Count,na.rm = TRUE),
            N = n())
df.summary.tmp

df.summary <- df.summary.tmp %>%
  mutate(Two.to.10 = Two.up-Eleven.up,
         Eleven.to.20 = Eleven.up-Twenty_one.up,
         Twenty_one.to.50 = Twenty_one.up-Fifty_one.up,
         Fifty_one.to.100 = Fifty_one.up - More.than.100,
         Percent.zero = 100*Zeros/N) %>%
  select(Organism:Ones, Two.to.10:Fifty_one.to.100,More.than.100, Max,N,Percent.zero)
print(df.summary, n=76)
median(df.summary$N) ## 7143
median(df.summary$Percent.zero) ## 99.4%
min(df.summary$Percent.zero) ## 56.8% for Semibalanus balaniodes
write_csv(df.summary, "Summary.csv")

## Number of zeros suggests there are transects and/or levels that never have certain organisms
## These are non-informative zeros and should be excluded
## Look at two species as examples

Species <- c("Tectura testudinalis","Strongylocentrotus droebachiensis")
9## Include filter of Transect and Level if you want to focus on highly sampled subset
df.species <- df5 %>%
  filter(Transect %in% c(5,7,15,20,22,26,28)) %>%
  filter(Level <= 3 & Level >= 0) %>%
  filter(Organism == Species[1]) %>%
  mutate(Present = if_else(Count == 0,0,1))
df.species
unique(df.species$Level)

## OK: arrange the data to look for non-informative zeros
df.species.table <- df.species %>%
  group_by(Transect,Level) %>%
  summarise(Sum_count = sum(Count,na.rm = TRUE),
            N = n(),na.rm = TRUE,
            Sum_present = sum(Present,na.rm = TRUE)) %>%
  mutate(p.count = Sum_present/N,na.rm = TRUE) %>%
  arrange(Level)
df.species.table

## A couple ways to look at the data
## Number of quadrats in each cell across all years

df.species.table %>%
  select(Transect,Level,N) %>%
  pivot_wider(names_from = Transect, values_from = N) %>%
  arrange(Level) %>%
  print(n=24) ## Well sampled

## Total number of limpets in each cell across all years
## The zeros are non-informative and those Level x Transect combinations should be dropped
df.species.table %>%
  select(Transect,Level,Sum_count) %>%
  pivot_wider(names_from = Transect, values_from = Sum_count) %>%
  print()
## Should drop Transect 22 and tidal levels >1.87 MLLW

## Proportion of quadrats in each cell across all years for which limpets were present
df.species.table %>%
  select(Transect,Level,p.count) %>%
  pivot_wider(names_from = Transect, values_from = p.count) %>%
  arrange(Level) %>%
  print()

## Now let's make two datasets to run gamma regressions
## Dataset 1: keep everything and average as explained on lines 279-282
## And filtered as explained on lines 282-285
df5
df.gamma1 <- df5 %>% # NOTE: here is the reviewer's mistake
  mutate(Present = if_else(Count == 0,0,1)) %>%
  group_by(Year,Transect,Level) %>%
  summarise(quad_mean = mean(Count, na.rm = TRUE),
            quad_count = sum(Present, na.rm = TRUE)) %>%
  group_by(Year,Level) %>%
  summarize(level_mean = mean(quad_mean, na.rm = TRUE),
            level_count = sum(quad_count, na.rm = TRUE))

## Find tidal levels where 1 or fewer occurrances of limpets
df.check <- df.gamma1 %>%
  group_by(Level) %>%
  summarise(sum_count = sum(level_count, na.rm = TRUE),
            sum_mean = mean(level_mean, na.rm = TRUE))
print(df.check, n=30)

## Limpets not found at 5.28 and 5.62 m MLLW; delete these levels
## add 0.1 if level_mean = 0
df.gamma1.final <- df.gamma1 %>%
  filter(Level <= 5) %>%
  mutate(Y = if_else(level_mean == 0,level_mean+0.1,level_mean))
unique(df.gamma1.final$Level)
df.gamma1.final %>% filter(level_mean == 0)

## Dataset 2: Drop non-informative zeros
## Step 1: join df5 and df.species table
## Step 2: drop if Sum_present = 0
## Repeat steps for df.gamma1.final
df5
df.species.table
df.gamma2 <- full_join(df5,df.species.table, by = c("Transect","Level")) %>%
  filter(Sum_present >=1) %>%
  select(Year,Transect,Level,Replicate,Count) %>%
  mutate(Present = if_else(Count == 0,0,1)) %>%
  group_by(Year,Transect,Level) %>%
  summarise(quad_mean = mean(Count, na.rm = TRUE),
            quad_count = sum(Present, na.rm = TRUE)) %>%
  group_by(Year,Level) %>%
  summarize(level_mean = mean(quad_mean, na.rm = TRUE),
            level_count = sum(quad_count, na.rm = TRUE))
df.gamma2
min(df.gamma2$level_mean) ## There are zeros

## Find tidal levels where 1 or fewer occurances of limpets
## All OK
df.check2 <- df.gamma2 %>%
  group_by(Level) %>%
  summarise(sum_count = sum(level_count, na.rm = TRUE),
            sum_mean = mean(level_mean, na.rm = TRUE))
print(df.check2, n=30)

## add 0.1 if level_mean = 0
df.gamma2.final <- df.gamma2 %>%
  mutate(Y = if_else(level_mean == 0,level_mean+0.1,level_mean))
df.gamma2.final %>% filter(level_mean == 0)

library(INLA) ## Bayesian analysis using integrated nested Laplace approximations
## This approach is much faster

formula <- Y ~ Year + as.factor(Level)
fm1 <- inla(formula, family = "gamma", data = df.gamma1.final)
round(fm1$summary.fixed[1:2,],4) ## Year log-slope = -0.011
fm2 <- inla(formula, family = "gamma", data = df.gamma2.final)
round(fm2$summary.fixed[1:2,],4) ## Year log-slope = +0.004


## A quick check of percent cover data
raw.pc <- read_data_entity(packageId = package_id, entityId = data_entities$entityId[1])
## Use vroom so you can check for problems
df.pc <- vroom(file=raw.pc)
df.pc ## Lots of NAs
unique(df.pc$Percent_cover)
df.pc.2 <- df.pc %>%
  filter(!is.na(Percent_cover))
df.pc.2


# end reviewer's code -----------------------------------------------------




# response: ===============================================================
# the reviewer makes some errors in their code that impede comparisons with ours.
# first, on line 226, the reviewer attempts to make a dataset to conduct a gamma
# regression for a target species: "Tectura testudinalis", but does not filter
# to that organism when creating their dataframe, and thus, is summarizing
# all organisms in "df5" (72 organisms). Additionally, the reviewer conducts some 
# steps that we describe differently than we did, or not at all, and also does not
# conduct some steps which we completed but did not describe in our initial submission.
# Here, we show here how how applying the steps that we described in the main text
# of our initial submission (and properly filterint to the target species) achieve
# a much closer result to that which we found and presented. 


# upload our processed count data, prepped for abundance modeling: 
# (here using the Full dataset, although note that the second submission of our
# manuscript uses the Highly Sampled dataset for this task)
data_jl <- readr::read_csv("data-processed/abundance-data/count-data-prepped-for-model.csv") 
data_jl <- data_jl %>%
  filter(organism == "Testudinalia testudinalis") %>%
  mutate(tidalheight = as.character(tidalheight)) %>%
  #filter to year <= 2017 to compare with EDI data
  filter(year <= 2017)

# test with our provided modeling workflow:
source("scripts/06-abundance-modeling/06.2-create-model-functions.R")
mod_jl <- find_regression_slopes(df = data_jl)
summary(mod_jl)
extract_slope(mod_jl)
# with data until 2017, we get a slope of -0.058 using our cleaned data

# try with the reviewer's method (INLA)
mod_jl_inla_formula <- density_nonzero ~ year + tidalheight
mod_jl_inla <- inla(mod_jl_inla_formula, family = "gamma", data = data_jl)
round(mod_jl_inla$summary.fixed[1:2,],4) ## Year log-slop
# coefficient = -0.058 again. 

# plot counts over time
data_jl %>%
  ggplot(aes(x = year, 
             y = density_nonzero)) +
  geom_point() +
  facet_wrap(~level)

# fix reviewer's code to parallel ours ----------------------------------------------------
df5 <- df4 %>%
  mutate(Organism = if_else(Organism == "Tectura testinalis","Tectura testudinalis",Organism),
         Organism = if_else(Organism == "Idotea baltica","Idotea balthica",Organism))
df5$tidalheight <- 0.3408*(13.5 - df5$Level)
unique(df5$tidalheight)
unique(df5$Organism)

# first, the reviewer didn't filter to the proper species:
# here Tectura testudinalis
tect <- df5 %>%
  filter(Organism == "Tectura testudinalis") 
  

# complete only the steps that we describe in our main text initial submission:  
tect_mt <- tect %>% 
  # convert counts to density / m^2
  mutate(density = Count*25) %>%
  # summarize replicates of each quadrat (year/transect/level combination), 
  group_by(Year, Transect, Level, Organism, tidalheight) %>%
  summarize(mean_density = mean(density, na.rm=T)) %>%
  # then summarize quadrat means per level across all transects
  group_by(Year, Level, Organism, tidalheight) %>%
  summarize(mean_level_density = mean(mean_density, na.rm=T)) %>%
  # remove levels in which the organism was only present once:
  # (which we describe in our manuscript)
  group_by(Level) %>%
  mutate(n_pres = sum(mean_level_density > 0)) %>%
  filter(n_pres > 1) %>%
  # change zero density values to 0.01 so the models will run
  # (not 0.1 as the reviewer did)
  mutate(density_nonzero = ifelse(mean_level_density == 0, 0.01, mean_level_density))
  
# plot density over time:
tect_mt %>%
  ggplot(aes(x = Year,
             y = density_nonzero)) +
  geom_point() +
  facet_wrap(~Level)

# model with reviewer's method
formula <- density_nonzero ~ Year + as.factor(tidalheight)
tect_mod_mt <- inla(formula, family = "gamma", data = tect_mt)
round(tect_mod_mt$summary.fixed[1:2,],4) ## Year log-slop
# coefficient = -0.064

# model with our method:
tect_mod_jl <- find_regression_slopes(df = tect_mt %>% mutate(year_zero = Year, 
                                                              tidalheight = as.character(tidalheight)))
summary(tect_mod_jl)
extract_slope(tect_mod_jl)
# using our method, we also get -0.064



# now, do all the steps we did --------------------------------------------
# complete extra steps which we didn't detain in the first submission
# of our main text, including:
# - fill in zeros for transects that were measured but species wasn't found
# - filter out zeros from transects in which species was never present
# - filter out levels above and below vertical limits
# - those done above:

# first find all samples taken
all_taken <- df %>%
  filter(Data_taken == "yes") %>%
  distinct(Year, Transect, Level, Replicate)

tect_full <- tect %>% 
  full_join(all_taken) %>%
  tidyr::replace_na(list(Organism = "Tectura testudinalis",
                         Count = 0))

tect_full$tidalheight <- 0.3408*(13.5 - tect_full$Level)

tect_correct <- tect_full %>%
  
  #  filter to transects in which the species was seen 
  # (a detail which we completed, but did not describe in the main text)
  group_by(Transect) %>%
  mutate(transect_max = max(Count)) %>%
  filter(transect_max > 0) %>%
  select(-transect_max) %>%
  
  # filter to levels between the highest and lowest occurrence of the species
  # (a detail which we completed, but did not describe in the first submission)
  ungroup() %>%
  mutate(max_level = max(Level[Count > 0]),
         min_level = min(Level[Count > 0])) %>%
  filter(Level <= max_level,
         Level >= min_level) %>%
  select(-max_level, -min_level) %>%
  
  # convert count to density (as we did and as we described)
  mutate(density = Count*25) %>%
  
  #  summarize replicates of each quadrat, 
  group_by(Year, Transect, Level, Organism, tidalheight) %>%
  summarize(mean_density = mean(density, na.rm=T)) %>%
  
  # then summarize quadrat means per level across all transects
  group_by(Year, Level, Organism,tidalheight) %>%
  summarize(mean_level_density = mean(mean_density, na.rm=T)) %>%
  
  # remove levels in which the organism was only present once:
  # (which we describe in our manuscript)
  group_by(Level,tidalheight) %>%
  mutate(n_pres = sum(mean_level_density > 0)) %>%
  filter(n_pres > 1) %>%
  
  # change zero density values to 0.01 so the models will run
  # (not 0.1 as the reviewer did)
  ungroup() %>%
  mutate(density_nonzero = ifelse(mean_level_density == 0, 0.01, mean_level_density))

# plot:
tect_correct %>%
  ggplot(aes(x = Year, 
             y = density_nonzero)) +
  geom_point() +
  facet_wrap(~Level)

## Model with INLA
formula <- density_nonzero ~ Year + as.factor(Level)
tect_mod <- inla(formula, family = "gamma", data = tect_correct)
round(tect_mod$summary.fixed[1:2,],4) ## Year log-slope = -0.0675

# now check our dataframe with inla:
formula_jl <- density_nonzero ~ year + tidalheight
inla_tect_jl <- inla(formula_jl, family = "gamma", data = data_jl)
round(inla_tect_jl$summary.fixed[1:2,],4) ## Year log-slope = -0.011

# Just from that, we have achieved similar slopes:
# from reviewer: -0.0598
# from our data: -0.0580
# the remaining differences are: 
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------



# second, filter to transects in which the species was seen 
# (a detail which we completed, but did not describe in the main text)
tect2 <- tect %>%  
  group_by(Transect) %>%
  mutate(transect_max = max(Count)) %>%
  filter(transect_max > 0) %>%
  select(-transect_max)

# third, filter to levels between the highest and lowest occurrence of the species
# (a detail which we completed, but did not describe in the first submission)
good_levels <- tect2 %>%
  ungroup() %>%
  filter(Count > 0) %>%
  summarize(max_level = max(Level),
            min_level = min(Level))
tect3 <- tect2 %>%
  filter(Level <= good_levels$max_level,
         Level >= good_levels$min_level)

# fourth, convert count to density (as we did and as we described)
tect4 <- tect3 %>%
  mutate(density = Count*25)

# fifth, summarize replicates of each quadrat, 
# then summarize quadrat means per level across all transects
tect5 <- tect4 %>% 
  group_by(Year, Transect, Level, Organism) %>%
  summarize(mean_density = mean(density, na.rm=T)) %>%
  group_by(Year, Level, Organism) %>%
  summarize(mean_level_density = mean(mean_density, na.rm=T))


# sixth, remove levels in which the organism was only present once:
# (which we describe in our manuscript)
tect6 <- tect5 %>%
  group_by(Level) %>%
  mutate(n_pres = sum(mean_level_density > 0)) %>%
  filter(n_pres > 1)


# sixth, change zero density values to 0.01 so the models will run
# (not 0.1 as the reviewer did)
tect7 <- tect6 %>%
  ungroup() %>%
  mutate(density_nonzero = ifelse(mean_level_density == 0, 0.01, mean_level_density))

# plot:
tect7 %>%
  ggplot(aes(x = Year, 
             y = density_nonzero)) +
  geom_point() +
  facet_wrap(~Level)

## This approach is much faster

formula <- density_nonzero ~ Year + as.factor(Level)
tect_mod <- inla(formula, family = "gamma", data = tect7)
round(tect_mod$summary.fixed[1:2,],4) ## Year log-slope = -0.011
# now check our dataframe with inla:
formula_jl <- density_nonzero ~ year + tidalheight
inla_tect_jl <- inla(formula_jl, family = "gamma", data = data_jl)
round(inla_tect_jl$summary.fixed[1:2,],4) ## Year log-slope = -0.011

# Just from that, we have achieved similar slopes:
# from reviewer: -0.0598
# from our data: -0.0580
# the remaining differences are: 