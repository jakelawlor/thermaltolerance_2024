# Assess intertidal height per species as a function of year
# ungrouped version: meaning no 5-year bins, as previously done. 

# SPLIT count vs cover species

# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)

# data --------------------------------------------------------------------
# upload PA with thermal tolerances from highly sampled dataset
pa_therm_filtered <- readr::read_csv(here::here("data-processed",
                                                "appledore-survey-data",
                                                "pa-with-therm",
                                                "pa-with-therm-highly-sampled.csv"))



# clean a few columns -----------------------------------------------
pa_hs <- pa_therm_filtered%>%
  tidyr::separate(organism, into = c("Genus","Species")) %>%
  mutate(G = paste0(stringr::str_sub(Genus,1,1),".")) %>%
  mutate(organism = paste(G, Species)) %>%
  select(-G, -Genus, -Species) %>%
  mutate(tidalheight = (13.5-level)*.3048)  
rm(pa_therm_filtered)

# find for only counts dataset
pa_hs_count <- pa_hs %>% filter(pres_count == T)
pa_hs_cover <- pa_hs %>% filter(pres_cover == T)
rm(pa_hs)

# find height quantiles
hs_quants_count <- pa_hs_count %>%
  
  # count total occurrences per species
  group_by(organism) %>%
  mutate(total_count = n()) %>%
  
  group_by(organism, year, total_count) %>%
  
  # find distributions of height occurrences 
  summarize(heightmax = max(tidalheight),
            height95 = quantile(tidalheight,.95),
            heightmed = median(tidalheight),
            heightmean = mean(tidalheight),
            height05 = quantile(tidalheight,.05),
            heightmin = min(tidalheight)) %>%
  
  # filter to organisms present in 10+ years
  group_by(organism) %>%
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years >= 10) %>%
  ungroup() %>%
  select(-n_years) 

hs_quants_cover <- pa_hs_cover %>%
  
  # count total occurrences per species
  group_by(organism) %>%
  mutate(total_count = n()) %>%
  
  group_by(organism, year, total_count) %>%
  
  # find distributions of height occurrences 
  summarize(heightmax = max(tidalheight),
            height95 = quantile(tidalheight,.95),
            heightmed = median(tidalheight),
            heightmean = mean(tidalheight),
            height05 = quantile(tidalheight,.05),
            heightmin = min(tidalheight)) %>%
  
  # filter to organisms present in 10+ years
  group_by(organism) %>%
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years >= 10) %>%
  ungroup() %>%
  select(-n_years) 


hs_quants_split_count <- hs_quants_count %>% split(f = .$organism)
hs_quants_split_cover <- hs_quants_cover %>% split(f = .$organism)
# find shared organisms
intersect(names(hs_quants_split_count),
          names(hs_quants_split_cover))

#rm(hs_quants)

# find exclusions ---------------------------------------------------------
# assessing qwuantile shifts won't work if the species is found all the
# way to the bounds of the domain in most years. 

# find domain
domain_count <-unique(pa_hs_count$tidalheight)[c(1,9)]
domain_cover <-unique(pa_hs_cover$tidalheight)[c(1,9)]
domain <- c(max(domain_count[1], domain_cover[1]),
            min(domain_count[2], domain_cover[2]))
rm(domain_count, domain_cover)

# identify species with 50% or more of their highest occurrences at domain edge
max_to_exclude_count <- purrr::map(
  .x = hs_quants_split_count,
  .f = ~.x %>%  summarise(n = n(),
                          at_domain = sum(heightmax == domain[1])) %>%
    mutate(at_domain_prop = at_domain/n) %>%
    mutate(exclude_max = at_domain_prop > .5) %>% 
    pull(exclude_max) %>% as_tibble
) %>% bind_rows(.id = "organism") %>%
  rename(exclude_max  = value)
table(max_to_exclude_count$exclude_max)

max_to_exclude_cover <- purrr::map(
  .x = hs_quants_split_cover,
  .f = ~.x %>%  summarise(n = n(),
                          at_domain = sum(heightmax == domain[1])) %>%
    mutate(at_domain_prop = at_domain/n) %>%
    mutate(exclude_max = at_domain_prop > .5) %>% 
    pull(exclude_max) %>% as_tibble
) %>% bind_rows(.id = "organism") %>%
  rename(exclude_max  = value)
table(max_to_exclude_cover$exclude_max)

# repeat for min
min_to_exclude_count <- purrr::map(
  .x = hs_quants_split_count,
  .f = ~.x %>%  summarise(n = n(),
                          at_domain = sum(heightmin == domain[2])) %>%
    mutate(at_domain_prop = at_domain/n) %>%
    mutate(exclude_min = at_domain_prop > .5) %>% 
    pull(exclude_min) %>% as_tibble
) %>% bind_rows(.id = "organism") %>%
  rename(exclude_min  = value)
table(min_to_exclude_count$exclude_min)

min_to_exclude_cover <- purrr::map(
  .x = hs_quants_split_cover,
  .f = ~.x %>%  summarise(n = n(),
                          at_domain = sum(heightmin == domain[2])) %>%
    mutate(at_domain_prop = at_domain/n) %>%
    mutate(exclude_min = at_domain_prop > .5) %>% 
    pull(exclude_min) %>% as_tibble
) %>% bind_rows(.id = "organism") %>%
  rename(exclude_min  = value)
table(min_to_exclude_cover$exclude_min)

exclusions_count <- max_to_exclude_count %>% left_join(min_to_exclude_count)
rm(min_to_exclude_count, max_to_exclude_count)

exclusions_cover <- max_to_exclude_cover %>% left_join(min_to_exclude_cover)
rm(min_to_exclude_cover, max_to_exclude_cover)

# add mid to exclude- when the species' max and min are both at domain bounds
exclusions_count <- exclusions_count %>%
  mutate(exclude_mid = case_when(exclude_max == T & exclude_min == T ~ T,
                                 TRUE ~ F))

table(exclusions_count$exclude_mid)

exclusions_cover <- exclusions_cover %>%
  mutate(exclude_mid = case_when(exclude_max == T & exclude_min == T ~ T,
                                 TRUE ~ F))

table(exclusions_cover$exclude_mid)


# model observation heights over time -------------------------------------
model <- function(df, param, group = "count"){
  
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
  
  # add excludions
  if(param %in% c("heightmean","heightmed")){
    mod_exclusion <- get(paste0("exclusions","_",group)) %>%
      select(organism, exclude_mid) %>%
      rename(exclude = exclude_mid)
  } else if(param %in% c("heightmin","height05")){
    mod_exclusion <-  get(paste0("exclusions","_",group)) %>%
      select(organism, exclude_min) %>%
      rename(exclude = exclude_min)
  } else if(param %in% c("heightmax","height95")){
    mod_exclusion <-  get(paste0("exclusions","_",group)) %>%
      select(organism, exclude_max) %>%
      rename(exclude = exclude_max)
  }
  
  # merge all
  out <- mods_p %>%
    left_join(mods_trend,
              by = join_by(organism)) %>%
    left_join(mods_single_val,
              by = join_by(organism)) %>%
    left_join(mod_exclusion, 
              by = join_by(organism)) %>%
    mutate(param = param) %>%
    relocate(param, .after = organism) %>%
    mutate(sign = case_when(single_value == "single value" ~ "(0)",
                            TRUE ~ sign))
  
  return(out)
  
}


#mods_min <- model(hs_unc_sum_split, "heightmin")
mods_05_count <- model(hs_quants_split_count, "height05", group = "count")
mods_05_cover <- model(hs_quants_split_cover, "height05", group = "cover")
mods_mean_count <- model(hs_quants_split_count, "heightmean", group = "count")
mods_mean_cover <- model(hs_quants_split_cover, "heightmean", group = "cover")
mods_med_count <- model(hs_quants_split_count, "heightmed", group = "count")
mods_med_cover <- model(hs_quants_split_cover, "heightmed", group = "cover")
mods_95_count <- model(hs_quants_split_count, "height95", group = "count")
mods_95_cover <- model(hs_quants_split_cover, "height95", group = "cover")
#mods_max <- model(hs_unc_sum_split, "heightmax")

# merge all ------------------------------------------------------

mods_all_count <- rbind( 
  mods_05_count,
  mods_mean_count,
  mods_med_count,
  mods_95_count
)

mods_all_cover <- rbind( 
  mods_05_cover,
  mods_mean_cover,
  mods_med_cover,
  mods_95_cover
)


mods_all_count %>% count(exclude)
mods_all_cover %>% count(exclude)

mods_all_annotate_count <- mods_all_count %>%
  mutate(sign = case_when(p_value < 0.05 ~ paste0(sign,"*"),
                          TRUE ~ sign))
mods_all_annotate_cover <- mods_all_cover %>%
  mutate(sign = case_when(p_value < 0.05 ~ paste0(sign,"*"),
                          TRUE ~ sign))




# plot --------------------------------------------------------------------

p1_count <- mods_all_annotate_count %>%
  # filter(exclude == F) %>% 
  group_by(param) %>%
  count(sign) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  mutate(text_color = if_else(sign %in% c("(-)*"), "white", "black")) %>%
  
  mutate(sign = factor(sign, 
                       levels = c("(+)*","(+)","(0)","(-)","(-)*"))) %>%
  mutate(param = factor(param,
                        levels = c("heightmin","height05",
                                   "heightmean","heightmed",
                                   "height95","heightmax"))) %>%
  filter(param %in% c("height95",'heightmed',"heightmean","height05")) %>%
  
  ggplot(aes(y = param, 
             x = prop,
             fill = sign)) +
  geom_col() +
  geom_vline(xintercept = .5, linetype = "dashed", linewidth = .25) +
  geom_vline(xintercept = c(0,1), linetype = "solid",
             linewidth = .25) +
  scale_fill_manual(breaks = c("(-)*","(-)","(0)","(+)","(+)*"),
                    values = c(c("darkmagenta",colorspace::lighten("magenta",.4),"grey93",
                                 colorspace::lighten("turquoise",.4), "darkturquoise")),
                    labels = c("Significant\nShift Down","Trend Down","No Change","Trend Up","Significant\nShift Up                   ")) +
  geom_text(data = . %>% group_by(param) %>% slice(1),
            aes(x = 1.03, 
                label = paste0("n = ",total)),
            stat = "unique",
            size = 3.5,
            hjust = 0) +
  theme_minimal() +
  geom_text(aes(label = n,
                color = text_color),
            position = position_stack(vjust = 0.5)) +
  scale_color_identity() +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  coord_cartesian(expand = F,
                  xlim = c(0, 1.2),
                  clip = "off") +
  scale_y_discrete(labels = c("5%\nOccurrence\nHeight",
                              "Mean\nOccurrence\nHeight",
                              "Median\nOccurrence\nHeight",
                              "95%\nOccurrence\nHeight")) +
  labs(x = "Proportion of Species",
       y = NULL,
       fill = "Shift Trend",
       title = "All Species") +
  # guides(fill = guide_legend(theme = theme(legend.title.position = "top"))) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = scales::percent_format())

p1_count


p1_cover <- mods_all_annotate_cover %>%
  # filter(exclude == F) %>% 
  group_by(param) %>%
  count(sign) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  mutate(text_color = if_else(sign %in% c("(-)*"), "white", "black")) %>%
  
  mutate(sign = factor(sign, 
                       levels = c("(+)*","(+)","(0)","(-)","(-)*"))) %>%
  mutate(param = factor(param,
                        levels = c("heightmin","height05",
                                   "heightmean","heightmed",
                                   "height95","heightmax"))) %>%
  filter(param %in% c("height95",'heightmed',"heightmean","height05")) %>%
  
  ggplot(aes(y = param, 
             x = prop,
             fill = sign)) +
  geom_col() +
  geom_vline(xintercept = .5, linetype = "dashed", linewidth = .25) +
  geom_vline(xintercept = c(0,1), linetype = "solid",
             linewidth = .25) +
  scale_fill_manual(breaks = c("(-)*","(-)","(0)","(+)","(+)*"),
                    values = c(c("darkmagenta",colorspace::lighten("magenta",.4),"grey93",
                                 colorspace::lighten("turquoise",.4), "darkturquoise")),
                    labels = c("Significant\nShift Down","Trend Down","No Change","Trend Up","Significant\nShift Up                   ")) +
  geom_text(data = . %>% group_by(param) %>% slice(1),
            aes(x = 1.03, 
                label = paste0("n = ",total)),
            stat = "unique",
            size = 3.5,
            hjust = 0) +
  theme_minimal() +
  geom_text(aes(label = n,
                color = text_color),
            position = position_stack(vjust = 0.5)) +
  scale_color_identity() +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  coord_cartesian(expand = F,
                  xlim = c(0, 1.2),
                  clip = "off") +
  scale_y_discrete(labels = c("5%\nOccurrence\nHeight",
                              "Mean\nOccurrence\nHeight",
                              "Median\nOccurrence\nHeight",
                              "95%\nOccurrence\nHeight")) +
  labs(x = "Proportion of Species",
       y = NULL,
       fill = "Shift Trend",
       title = "All Species") +
  # guides(fill = guide_legend(theme = theme(legend.title.position = "top"))) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = scales::percent_format())
p1_cover




p2_count <- mods_all_annotate_count %>%
  filter(exclude == F) %>% 
  group_by(param) %>%
  count(sign) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  mutate(text_color = if_else(sign %in% c("(-)*"), "white", "black")) %>%
  
  mutate(sign = factor(sign, 
                       levels = c("(+)*","(+)","(0)","(-)","(-)*"))) %>%
  mutate(param = factor(param,
                        levels = c("heightmin","height05",
                                   "heightmean","heightmed",
                                   "height95","heightmax"))) %>%
  filter(param %in% c("height95",'heightmed',"heightmean","height05")) %>%
  
  ggplot(aes(y = param, 
             x = prop,
             fill = sign)) +
  geom_col() +
  geom_vline(xintercept = .5, linetype = "dashed", linewidth = .25) +
  geom_vline(xintercept = c(0,1), linetype = "solid",
             linewidth = .25) +
  scale_fill_manual(breaks = c("(-)*","(-)","(0)","(+)","(+)*"),
                    values = c(c("darkmagenta",colorspace::lighten("magenta",.4),"grey93",
                                 colorspace::lighten("turquoise",.4), "darkturquoise")),
                    labels = c("Significant\nShift Down","Trend Down","No Change","Trend Up","Significant\nShift Up                   ")) +
  geom_text(aes(x = 1.03, 
                label = paste0("n = ",total)),
            stat = "unique",
            size = 3.5,
            hjust = 0) +
  theme_minimal() +
  geom_text(aes(label = n,
                color = text_color),
            position = position_stack(vjust = 0.5)) +
  scale_color_identity() +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  coord_cartesian(expand = F,
                  xlim = c(0, 1.2),
                  clip = "off") +
  scale_y_discrete(labels = c("5%\nOccurrence\nHeight",
                              "Mean\nOccurrence\nHeight",
                              "Median\nOccurrence\nHeight",
                              "95%\nOccurrence\nHeight")) +
  labs(x = "Proportion of Species",
       y = NULL,
       fill = "Shift Trend",
       title = "Species Inside Domain") +
  #guides(fill = guide_legend(theme = theme(legend.title.position = "top"))) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = scales::percent_format())
p2_count

p2_cover <- mods_all_annotate_cover %>%
  filter(exclude == F) %>% 
  group_by(param) %>%
  count(sign) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  mutate(text_color = if_else(sign %in% c("(-)*"), "white", "black")) %>%
  
  mutate(sign = factor(sign, 
                       levels = c("(+)*","(+)","(0)","(-)","(-)*"))) %>%
  mutate(param = factor(param,
                        levels = c("heightmin","height05",
                                   "heightmean","heightmed",
                                   "height95","heightmax"))) %>%
  filter(param %in% c("height95",'heightmed',"heightmean","height05")) %>%
  
  ggplot(aes(y = param, 
             x = prop,
             fill = sign)) +
  geom_col() +
  geom_vline(xintercept = .5, linetype = "dashed", linewidth = .25) +
  geom_vline(xintercept = c(0,1), linetype = "solid",
             linewidth = .25) +
  scale_fill_manual(breaks = c("(-)*","(-)","(0)","(+)","(+)*"),
                    values = c(c("darkmagenta",colorspace::lighten("magenta",.4),"grey93",
                                 colorspace::lighten("turquoise",.4), "darkturquoise")),
                    labels = c("Significant\nShift Down","Trend Down","No Change","Trend Up","Significant\nShift Up                   ")) +
  geom_text(aes(x = 1.03, 
                label = paste0("n = ",total)),
            stat = "unique",
            size = 3.5,
            hjust = 0) +
  theme_minimal() +
  geom_text(aes(label = n,
                color = text_color),
            position = position_stack(vjust = 0.5)) +
  scale_color_identity() +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  coord_cartesian(expand = F,
                  xlim = c(0, 1.2),
                  clip = "off") +
  scale_y_discrete(labels = c("5%\nOccurrence\nHeight",
                              "Mean\nOccurrence\nHeight",
                              "Median\nOccurrence\nHeight",
                              "95%\nOccurrence\nHeight")) +
  labs(x = "Proportion of Species",
       y = NULL,
       fill = "Shift Trend",
       title = "Species Inside Domain") +
  #guides(fill = guide_legend(theme = theme(legend.title.position = "top"))) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = scales::percent_format())
p2_cover


both_p1 <- (p1_count + labs(title = "Density Dataset Only") |
              p1_cover + labs(title = "Percent Cover Dataset Only")) + plot_layout(guides = "collect",
                                               axes = "collect") &
  theme(legend.position = "bottom",
        legend.margin = ggplot2::margin(t = -5,0,0,0,"pt"),
        plot.margin = ggplot2::margin(0,0,0,0,"pt"))

both_p1 + ggview::canvas(8,3.5)

both_p2 <- (p2_count + labs(title = "Density HS Dataset") |
              p2_cover + labs(title = "Percent Cover HS Dataset")) + plot_layout(guides = "collect",
                                                                                   axes = "collect") &
  theme(legend.position = "bottom",
        legend.margin = ggplot2::margin(t = -5,0,0,0,"pt"),
        plot.margin = ggplot2::margin(r = 5,0,0,0,"pt"))

both_p2 + ggview::canvas(8,3.5)




# save as rds -------------------------------------------------------------

out <- list(p1_count, p1_cover, p2_count, p2_cover)
names(out) <- c("all_count","all_cover","in_domain_count","in_domain_cover")

saveRDS(out,
        file = "outputs/depth_shifts/depth_proportions_by_group.rds")


rm(list = ls())