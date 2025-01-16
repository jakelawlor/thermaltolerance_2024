# prep percent cover species for abundance modeling


library(dplyr)
library(ggplot2)



# data --------------------------------------------------------------------
cover <- readr::read_csv(here::here(
  "data-processed",
  "cover_abundance_filtered.csv"
))

cover %>% glimpse()
cover %>% filter(is.na(pc_num)) %>% distinct(percent_cover)
cover %>% distinct(organism)
# note that we start with 65 organisms total

# model abundance change --------------------------------------------------
# first, find mean by replicate 
# because most do not have replicates, but some have up to 4
# find the mean between replicates for each sample
cover2 <- cover %>%
  group_by(year, transect, level, organism) %>% 
  summarize(
    # cout how many replicates in each group
    n_rep = n(),
    # average percent cover between replicates 
    # removing NAs here, which are correct and should not be zero
    pc_num = mean(pc_num,na.rm=T),
    # find if the organism is present in any replicate
    pres = any(pres)) %>% 
  ungroup() 
range(cover2$n_rep)
mean(cover2$n_rep)
cover2 %>% glimpse()
rm(cover)


# now find mean per level across all transects
# again, we're doing this because the number and selection of 
# transects varies widely from year to year, but we just want 
# mean abundance per tidal level
cover3 <- cover2 %>%
  group_by(organism, year, level) %>%
  # find mean percent cover per tidal level
  summarize(pc_num = mean(pc_num,na.rm = T),
            n_transects = n()) %>%
  
  ungroup() %>%
  
  # add year column starting at 1
  # doing this so that the coefficient is easier to interpret when modeling
  mutate(year_zero = year - ( min(year)-1)) %>%
  
  # translate level into tidal height
  mutate(tidalheight = (13-level)*.348)  
cover3 %>% glimpse()
rm(cover2)
range(cover3$n_transects)
mean(cover3$n_transects)



# filter to only species present in 3+ years of the data ------------------
multiyear_spp <- cover3 %>% 
  filter(pc_num > 0) %>%
  group_by(organism) %>%
  distinct(year) %>%
  count(name = "years_pres")
cover3 %>% distinct(organism) # 65 organisms to start

cover4 <- cover3 %>%
  left_join(multiyear_spp) %>%
  filter(years_pres > 3)
cover4 %>% distinct(organism)
# by subsetting the data to only organisms in 3+ years of sampling, 
# we cut down to 49 species
rm(cover3)
multiyear_spp %>%
  filter(years_pres <= 3)
rm(multiyear_spp)


# remove instances where an organism only has data (including zeros)
# once in a given level
cover5 <- cover4 %>%
  group_by(organism, level) %>% 
  add_count(name = "nrow_in_level") %>%
  filter(nrow_in_level > 1) %>%
  select(-nrow_in_level) %>%
  ungroup()
rm(cover4)

# some tidal levels are never sampled again past 2004 - remove those
# technically, we'll remove levels that are sampled for < 20 years
levels_sampled <- 
  cover5 %>%
  group_by(level) %>%
  distinct(year) %>% 
  count(name = "n_years_level_sampled")

cover6 <- cover5 %>%
  left_join(levels_sampled) %>%
  filter(n_years_level_sampled > 20) %>%
  select(-n_years_level_sampled)
rm(cover5, levels_sampled)


# first, make dataset with percent cover ranging 0-1
cover7 <- cover6 %>%
  mutate(percent_cover_beta = pc_num/100) %>%
  filter(!is.na(pc_num)) 
rm(cover6)
range(cover7$percent_cover_beta)


readr::write_csv(cover7,
                 here::here("data-processed",
                            "abundance-data",
                            "cover-data-prepped-for-model.csv"))



# 3. Subset to test -------------------------------------------------------
# subset data to only one species to test modeling function
testdf <- cover7 %>%
  filter(organism == "Chondrus crispus") 
hist(testdf$pc_num)

testdf %>%
  ggplot(aes(x=year_zero, y=pc_num)) +
  geom_point() +
  geom_smooth(method="lm")

testdf %>%
  ggplot(aes(x=year_zero, y=percent_cover_beta)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~tidalheight, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format())
# in this case, trends show clear declines in percent cover 
# at basically every tidal height over time


# plot trends over time with just lm trend line
example_unmodeled <- 
  testdf %>%
  #group_by(level) %>% add_count() %>% 
  # filter(n>3) %>% filter(level <= 13) %>%
  mutate(height_label = paste(round(tidalheight,2), "m")) %>%
  mutate(height_label = forcats::fct_reorder(height_label, level)) %>%
  ggplot(aes(x=year_zero, y=pc_num)) +
  geom_hline(yintercept = 0)+
  geom_point(size=.75, alpha=.5) +
  geom_smooth(method = "lm", color = "black", fill = "#3499ad", size=.75) +
  theme(panel.border = element_blank()) +
  facet_wrap(~height_label) +
  coord_cartesian(ylim = c(0,80)) +
  scale_y_continuous(breaks = c(0,25,50,75),
                     labels = scales::percent_format(scale=1)) +
  labs(x=NULL,
       y=NULL,
       title = "Chondrus crispus percent cover change over time") +
  theme(plot.title = element_text(size=18),
        plot.title.position = "plot") 

example_unmodeled

ggview::ggview(example_unmodeled, height = 5, width = 7, units = "in")


rm(list = ls())
