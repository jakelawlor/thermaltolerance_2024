# find height changes weighted by abundance


# get highly sampled quadrats
hs <- readr::read_csv("data-processed/cti-data/cti-data-highly-sampled.csv")
hs <- hs %>% distinct(year, transect, level) %>%
  mutate(tidalheight = (13-level)*.348)  


cover <- readr::read_csv("data-processed/appledore-survey-data/cover_abundance_filtered.csv") %>%
  inner_join(hs)

cover %>% glimpse()


count <- readr::read_csv("data-processed/appledore-survey-data/counts_abundance_filtered.csv") %>%
  inner_join(hs) 




# start with count --------------------------------------------------------
# first, filter to species present in 10+ years
count2 <- count %>%
  group_by(organism) %>%
  
  # drop NAs for count and keep only >0
  tidyr::drop_na(countnum) %>%
  filter(count > 0) %>%
  
  # filter to 10+ years
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years > 10) %>%
  select(-n_years) %>%
  
  # filter to replicates 1-3
  filter(replicate %in% c(1:3))

# now find distributions of actual count by organism
count_dist <- count2 %>%
  group_by(organism, year) %>%
  
  # add vector of tidal heights of all occurrences
  summarize(occ_heights = list(rep(tidalheight, countnum))) %>%
  
  # find quantiles
  mutate(
    q5 = map_dbl(occ_heights, ~ quantile(.x, 0.05, names = FALSE)),
    q50 = map_dbl(occ_heights, ~ quantile(.x, 0.50, names = FALSE)),
    q95 = map_dbl(occ_heights, ~ quantile(.x, 0.95, names = FALSE)),
    mean = map_dbl(occ_heights, ~ mean(.x))
  ) %>%
  select(organism, year, q5, q50, q95, mean) 


count_dist %>%
  tidyr::pivot_longer(cols = c("q5","q50","q95","mean"),
                      names_to = "param",
                      values_to = "tidalheight") %>%
  ggplot(aes(x = year, 
             y = tidalheight,
             color = param)) +
  geom_point(size = .2,
             alpha = .4) +
  geom_line(linewidth = .2, alpha = .4) + 
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~organism) +
  geom_hline(yintercept =  c(min(hs$tidalheight), max(hs$tidalheight)),
             linewidth = .25, 
             linetype = "dashed")



# repeat for cover --------------------------------------------------------
# first, filter to species present in 10+ years
cover2 <- cover %>%
  group_by(organism) %>%
  
  # drop NAs for count and keep only >0
  tidyr::drop_na(pc_num) %>%
  filter(pc_num > 0) %>%
  
  # filter to 10+ years
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years > 10) %>%
  select(-n_years) %>%
  
  # filter to replicates 1-3
  filter(replicate %in% c(1:3))

# now find distributions of actual count by organism
cover_dist <- cover2 %>%
  group_by(organism, year) %>%
  
  # add vector of tidal heights of all occurrences
  summarize(occ_heights = list(rep(tidalheight, pc_num))) %>%
  
  # find quantiles
  mutate(
    q5 = map_dbl(occ_heights, ~ quantile(.x, 0.05, names = FALSE)),
    q50 = map_dbl(occ_heights, ~ quantile(.x, 0.50, names = FALSE)),
    q95 = map_dbl(occ_heights, ~ quantile(.x, 0.95, names = FALSE)),
    mean = map_dbl(occ_heights, ~ mean(.x))
  ) %>%
  select(organism, year, q5, q50, q95, mean) 


cover_dist %>%
  tidyr::pivot_longer(cols = c("q5","q50","q95","mean"),
                      names_to = "param",
                      values_to = "tidalheight") %>%
  ggplot(aes(x = year, 
             y = tidalheight,
             color = param)) +
  geom_point(size = .2,
             alpha = .4) +
  geom_line(linewidth = .2, alpha = .4) + 
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~organism) +
  geom_hline(yintercept =  c(min(hs$tidalheight), max(hs$tidalheight)),
             linewidth = .25, 
             linetype = "dashed")




# now build models --------------------------------------------------------
model <- function(df, param){
  
  mods <- purrr::map(.x = df, 
                     .f = ~lm(data = .x, 
                              formula = get(param) ~ year))
  
  mods_p <- purrr::map(.x = mods,
                       .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
  ) %>%
    bind_rows(.id = "organism") %>%
    rename(p_value = value)
  
  
  mods_trend <- purrr::map(
    .x = mods,
    .f = ~ data.frame(trend = summary(.x)$coefficients["year","Estimate"],
                      sign = ifelse(sign(summary(.x)$coefficients["year","Estimate"]) == 1,
                                    "(+)",
                                    "(-)"
                      )) )%>%
    bind_rows(.id = "organism") 
  
  mods_single_val <- purrr::map(
    .x = df,
    .f = ~ifelse(length(unique(.x[[param]])) == 1,"single value","no issue") %>% as_tibble()
  ) %>%
    bind_rows(.id = "organism") %>%
    rename(single_value = value)
  
#  # add excludions
#  if(param %in% c("mean","q50")){
#    mod_exclusion <- exclusions %>%
#      select(organism, exclude_mid) %>%
#      rename(exclude = exclude_mid)
#  } else if(param %in% c("q5")){
#    mod_exclusion <- exclusions %>%
#      select(organism, exclude_min) %>%
#      rename(exclude = exclude_min)
#  } else if(param %in% c("q95")){
#    mod_exclusion <- exclusions %>%
#      select(organism, exclude_max) %>%
#      rename(exclude = exclude_max)
#  }
  
  # merge all
  out <- mods_p %>%
    left_join(mods_trend,
              by = join_by(organism)) %>%
    left_join(mods_single_val,
              by = join_by(organism)) %>%
   # left_join(mod_exclusion, 
   #           by = join_by(organism)) %>%
    mutate(param = param) %>%
    relocate(param, .after = organism) %>%
    mutate(sign = case_when(single_value == "single value" ~ "(0)",
                            TRUE ~ sign))
  
  return(out)
  
}


#mods_min <- model(hs_unc_sum_split, "heightmin")
mods_05_count <- model(count_dist %>% split(f = .$organism), "q5")
mods_05_count %>% count(p_value < .05)
mods_05_count %>% count(sign)
mods_05_cover <- model(cover_dist %>% split(f = .$organism), "q5")
mods_05_cover %>% count(p_value < .05)
mods_05_cover %>% count(sign)
mods_05_cover %>% rbind(mods_05_count) %>%  count(sign)
mods_05_cover %>% rbind(mods_05_count) %>%  count(sign, p_value < .05)

mods_95_count <- model(count_dist %>% split(f = .$organism), "q95")
mods_95_count %>% count(p_value < .05)
mods_95_count %>% count(sign)
mods_95_cover <- model(cover_dist %>% split(f = .$organism), "q95")
mods_95_cover %>% count(p_value < .05)
mods_95_cover %>% count(sign)
mods_95_cover %>% rbind(mods_95_count) %>% count(sign)
mods_95_cover %>% rbind(mods_95_count) %>%  count(sign, p_value < .05)



mods_50_count <- model(count_dist %>% split(f = .$organism), "q50")
mods_50_count %>% count(p_value < .05)
mods_50_count %>% count(sign)
mods_50_cover <- model(cover_dist %>% split(f = .$organism), "q50")
mods_50_cover %>% count(p_value < .05)
mods_50_cover %>% count(sign)
mods_50_cover %>% rbind(mods_50_count) %>% count(sign)
mods_50_cover %>% rbind(mods_50_count) %>%  count(sign, p_value < .05)



mods_mean_count <- model(count_dist %>% split(f = .$organism), "mean")
mods_mean_count %>% count(p_value < .05)
mods_mean_count %>% count(sign)
mods_mean_cover <- model(cover_dist %>% split(f = .$organism), "mean")
mods_mean_cover %>% count(p_value < .05)
mods_mean_cover %>% count(sign)
mods_mean_cover %>% rbind(mods_mean_count) %>% count(sign)
mods_mean_cover %>% rbind(mods_mean_count) %>%  count(sign, p_value < .05)



mods_mean <- model(hs_quants_split, "heightmean")
mods_med <- model(hs_quants_split, "heightmed")
mods_95 <- model(hs_quants_split, "height95")
#mods_max <- model(hs_unc_sum_split, "heightmax")

