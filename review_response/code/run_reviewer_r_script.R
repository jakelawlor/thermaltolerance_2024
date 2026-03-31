
# R SCRIPT

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
df5$Level_orig <- df5$Level
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
#write_csv(df.summary, "Summary.csv")


## Number of zeros suggests there are transects and/or levels that never have certain organisms
## These are non-informative zeros and should be excluded
## Look at two species as examples
Species <- c("Tectura testudinalis","Strongylocentrotus droebachiensis")


## Include filter of Transect and Level if you want to focus on highly sampled subset
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
df.gamma1 <- df.species %>%
  
  # NOTE: the following line was not included in the reviewer's script
  # but we assume the review meant to use this instead of averaging all species
  filter(str_detect(Organism, "Tectura")) %>%
  
  mutate(Present = if_else(Count == 0,0,1)) %>%
  group_by(Year,Transect,Level, Level_orig) %>%
  summarise(quad_mean = mean(Count, na.rm = TRUE),
            quad_count = sum(Present, na.rm = TRUE)) %>%
  group_by(Year,Level,Level_orig) %>%
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
  
  # NOTE: the following line was not included in the reviewer's original code:
  #filter(Organism == "Tectura testudinalis") %>%
  
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





# JL addition: ------------------------------------------------------------
# test our data with INLA package

# upload processed count data ready for abundance modeling:
data <- readr::read_csv("data-processed/abundance-data/count-data-prepped-for-model_HS.csv") 
data <- data %>%
  mutate(tidalheight = as.character(tidalheight))

testdf <- data %>%
  filter(organism == "Testudinalia testudinalis") 

testdf %>% 
  ggplot(aes(x = year, 
             y = density_nonzero))+
  geom_point()

df.gamma2.final %>%
  ggplot(aes(x = Year, 
             y = Y)) + geom_point()

# test with JL modeling workflow:
source("scripts/06-abundance-modeling/06.2-create-model-functions.R")
mod <- find_regression_slopes(df = testdf)
summary(mod)
extract_slope(mod)

# test with INLA
testdf_formula <- density_nonzero ~ year + as.factor(level)
testdfmod <- inla(testdf_formula, family = "gamma", data = testdf)
round(testdfmod$summary.fixed[1:2,],4) ## Year log-slope = -0.011

