# Assess when/where species enter and exit


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)


# data --------------------------------------------------------------------
pa_therm_filtered <- readr::read_csv(here::here("data-processed",
                                                "appledore-survey-data",
                                                "pa-with-therm",
                                                "pa-with-therm-highly-sampled.csv"))


pa <- readr::read_csv(
#  here::here(
#    "data-processed",
#    "appledore-survey-data",
#    "spp_pres_by_replicate.csv"
#  )
  here::here("data-processed",
             "appledore-survey-data",
             "pa-with-therm",
             "pa-with-therm-all.csv")) %>%
  tidyr::separate(organism, into = c("Genus","Species")) %>%
  mutate(G = paste0(stringr::str_sub(Genus,1,1),".")) %>%
  mutate(organism = paste(G, Species)) %>%
  select(-G, -Genus, -Species)

unique <- unique(pa$year) %>% sort()
y1 <- unique[1:5]
y2 <- unique[6:10]
y3 <- unique[11:15]
y4 <- unique[16:20]
y5 <- unique[21:25]
y6 <- unique[26:30]
y7 <- unique[31:35]
y8 <- unique[36:40]

org_by_yeargroup <- pa %>% 
  #filter(replicate == 1) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  mutate(year_group = case_when(
    year %in% y1 ~ 1,
    year %in% y2 ~ 2,
    year %in% y3 ~ 3,
    year %in% y4 ~ 4,
    year %in% y5 ~ 5,
    year %in% y6 ~ 6,
    year %in% y7 ~ 7,
    year %in% y8 ~ 8
  ))

minmax_samples <- org_by_yeargroup %>%
  group_by(year_group) %>%
  summarize(max_sample = max(tidalheight),
            min_sample = min(tidalheight))

org_by_yeargroup_sum <- org_by_yeargroup %>%
  group_by(organism, year_group) %>%
  
  summarize(height95 = quantile(tidalheight,.95),
            heightmed = median(tidalheight),
            heightmean = mean(tidalheight),
            height05 = quantile(tidalheight,.05)) %>%
  
  # filter to organisms present in 3+ year groups
  group_by(organism) %>%
  mutate(n_yeargroups = n_distinct(year_group)) %>%
  filter(n_yeargroups > 3) %>%
  ungroup() %>%
  select(-n_yeargroups)
  
org_by_yeargroup_sum_split <- org_by_yeargroup_sum %>% split(f = .$organism)


# model and assess mean observation heights
mean_mods <- purrr::map(
  .x = org_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "heightmean ~ year_group")
)
mean_mods_p <- purrr::map(
  .x = mean_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(mean_p = value)
mean_mods_trend <- purrr::map(
  .x = mean_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                        "(+)",
                                        "(-)"
                )) )%>%
  bind_rows(.id = "organism") %>%
  rename(mean_trend = trend_value,
         mean_sign = trend_sign)


# model and assess .05 quantile observation heights
min_mods <- purrr::map(
  .x = org_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "height05 ~ year_group")
)
min_mods_p <- purrr::map(
  .x = min_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(min_p = value)
min_mods_trend <- purrr::map(
  .x = min_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                        "(+)",
                                        "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(min_trend = trend_value,
         min_sign = trend_sign)
min_perfect <- purrr::map(
  .x = org_by_yeargroup_sum_split,
  .f = ~ifelse(length(unique(.x$height05)) == 1,"single value","no issue")
)
min_mods_p[min_perfect == "single value","min_p"] <- 1
min_mods_trend[min_perfect == "single value","min_trend"] <- 0
min_mods_trend[min_perfect == "single value","min_sign"] <- "None"


# model and assess .95 quantile observation heights
max_mods <- purrr::map(
  .x = org_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "height95 ~ year_group")
)
max_mods_p <- purrr::map(
  .x = max_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(max_p = value)
max_mods_trend <- purrr::map(
  .x = max_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                   "(+)",
                                   "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(max_trend = trend_value,
           max_sign = trend_sign)
max_perfect <- purrr::map(
  .x = org_by_yeargroup_sum_split,
  .f = ~ifelse(length(unique(.x$height95)) == 1,"single value","no issue")
)
max_mods_p[max_perfect == "single value","max_p"] <- 1
max_mods_trend[max_perfect == "single value","max_trend"] <- 0
max_mods_trend[max_perfect == "single value","max_sign"] <- "None"

# model and assess median  observation heights
med_mods <- purrr::map(
  .x = org_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "heightmed ~ year_group")
)
med_mods_p <- purrr::map(
  .x = med_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(med_p = value)
med_mods_trend <- purrr::map(
  .x = med_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                        "(+)",
                                        "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(med_trend = trend_value,
         med_sign = trend_sign)
med_perfect <- purrr::map(
  .x = org_by_yeargroup_sum_split,
  .f = ~ifelse(length(unique(.x$heightmed)) == 1,"single value","no issue")
)
med_mods_p[med_perfect == "single value","med_p"] <- 1
med_mods_trend[med_perfect == "single value","med_trend"] <- 0
med_mods_trend[med_perfect == "single value","med_sign"] <- "None"


# merge all p values ------------------------------------------------------
all_p <- mean_mods_p %>%
  left_join(mean_mods_trend) %>%
  left_join(max_mods_p) %>%
  left_join(max_mods_trend) %>%
  left_join(min_mods_p) %>%
  left_join(min_mods_trend) %>%
  left_join(med_mods_p) %>%
  left_join(med_mods_trend)


all_p <- all_p %>%
  mutate(mean_p_sig = case_when(mean_p < .05 ~ paste0("<b>mean ",mean_sign,"</b><br>"),
                                TRUE ~ ""),
         max_p_sig = case_when(max_p < .05 ~ paste0("<b>max ",max_sign,"</b><br>"),
                                TRUE ~ ""),
         min_p_sig = case_when(min_p < .05 ~ paste0("<b>min ",min_sign,"</b><br>"),
                                TRUE ~ ""),
         med_p_sig = case_when(med_p < .05 ~ paste0("<b>med ",med_sign,"</b><br>"),
                                TRUE ~ "")) %>%
  mutate(any_p = paste0(max_p_sig, mean_p_sig, med_p_sig, min_p_sig))

# plot all ----------------------------------------------------------------
library(ggtext)

all_data_trends_p <- org_by_yeargroup %>%
  filter(organism %in% names(org_by_yeargroup_sum_split)) %>%
  ggplot(aes(x = as.character(year_group),
             y = tidalheight)) +
  
  # add min and max sampling 
  geom_segment(data = minmax_samples,
               aes(x = year_group -.45,
                   xend = year_group + .45,
                   y = min_sample,
                   yend = min_sample),
               linetype = "dashed",
               linewidth = .2) +
  geom_segment(data = minmax_samples,
               aes(x = year_group -.45,
                   xend = year_group + .45,
                   y = max_sample,
                   yend = max_sample),
               linetype = "dashed",
               linewidth = .2) +
  
  geom_boxplot(outliers = F,
               color = "grey50",
               alpha = .5,
               linewidth = .2,
               width = .3) +
  facet_wrap(~organism) +
  geom_point(data = org_by_yeargroup_sum,
             aes(y = height95), color = "darkturquoise",
             size = .25) +
  geom_smooth(data = org_by_yeargroup_sum,
              method= "lm",
             aes(x = year_group,y = height95), 
             color = "darkturquoise",
             se = F,
             size = .25) +
  geom_point(data = org_by_yeargroup_sum,
             aes(y = height05), color = "darkturquoise",
             size = .25) +
  geom_smooth(data = org_by_yeargroup_sum,
              method= "lm",
              aes(x = year_group,y = height05), 
              color = "darkturquoise",
              se = F,
              size = .25) +
  geom_point(data = org_by_yeargroup_sum,
             aes(y = heightmean), color = "darkturquoise",
             size = .25) +
  geom_smooth(data = org_by_yeargroup_sum,
              method= "lm",
              aes(x = year_group,y = heightmean), 
              color = "darkturquoise",
              se = F,
              size = .25) +
  
  # add max annotation
  geom_richtext(data = all_p,
            aes(x = .5, y = -2,
                label = max_p_sig),
            fill = NA,
            label.colour = NA,
            lineheight = .8,
            size  = 2.5,
            vjust =.5,
            hjust = 0) +
  
  geom_richtext(data = all_p,
                aes(x = 4.5, y = -2,
                    label = mean_p_sig),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =.5,
                hjust = .5) +
  
  geom_richtext(data = all_p,
                aes(x = 8.5, y = -2,
                    label = min_p_sig),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =.5,
                hjust = 1) +
  
  theme(strip.text = element_text(margin = margin(0,0,0,0, "pt")),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_x_discrete(labels = c("1982",#-86",
                              "1987",#-91",
                              "1992",#-96",
                              "1997",#-01",
                              "2002",#-06",
                              "2008",#-12",
                              "2013",#-17",
                              "2019"#-23"
                              )) +
  labs(x = "Year Group",
       y = "Tidal Level (m)",
       title = "Species Depth Shifts - All Data")

all_data_trends_p + ggview::canvas(12,8)
ggsave(all_data_trends_p, 
       file = "outputs/depth_shifts/all_data_trends_p.png",
       width = 12,
       height = 8)

# plot number of changes
all_p %>% select(-contains("sig"), -any_p) %>%
  tidyr::pivot_longer(cols = c(mean_p, mean_trend, 
                               max_p, max_trend,
                               min_p, min_trend,
                               med_p, med_trend),
                      names_to = c("stat", ".value"),
                      names_sep = "_") %>%
  mutate(sig = case_when(p < .05 & trend == "(+)" ~ "Shift Up",
                         p < .05 & trend == "(-)" ~ "Shift Down",
                         p > .05 ~ "No Change")) %>%
  mutate(sig = factor(sig, levels = c("Shift Up","No Change","Shift Down"))) %>%
  mutate(stat = factor(stat, levels = c("min","med","mean","max"))) %>%
  ggplot(aes(x = stat,
             fill = sig)) +
  geom_bar() +
  coord_flip() +
  scale_x_discrete(labels = ~paste(stringr::str_to_title(.),"\nObservation \nHeight")) +
  scale_fill_manual(values = c("tomato2","grey85","cyan4")) +
  labs(x = NULL,
       y = "Count",
       fill = "Modeled Shift") +
  theme(legend.position = "top")

all_p2 <- all_p %>% select(-ends_with("sig"), -any_p) %>%
  tidyr::pivot_longer(cols = c(mean_p, mean_sign, 
                               max_p, max_sign,
                               min_p, min_sign,
                               med_p, med_sign),
                      names_to = c("stat", ".value"),
                      names_sep = "_") %>%
  mutate(sig = case_when(p < .05 & sign == "(+)" ~ "Sig Shift Up",
                         p > .05 & sign == "(+)" ~ "Shift Up",
                         p < .05 & sign == "(-)" ~ "Sig Shift Down",
                         p > .05 & sign == "(-)" ~ "Shift Down",
                         sign == "None" ~ "No Change")) %>%
  mutate(sig = factor(sig, levels = c("Sig Shift Up","Shift Up","No Change","Shift Down","Sig Shift Down"))) %>%
  mutate(stat = factor(stat, levels = c("min","med","mean","max"))) 

#all_p2 %>% 
#  filter(stat == "min") %>%
#  count(trend)
#  filter(stringr::str_detect(sig,"Sig")) %>% distinct(organism)
#  group_by(stat) %>%
#  summarize(n_sig_total = sum(stringr::str_detect(sig, "Sig Shift Up")))
#  filter(stringr::str_detect(sig,"Sig"))

all_data_depth_bars <- all_p2 %>% 
  ggplot(aes(x = stat,
             fill = sig)) +
  geom_bar(color = "black", width = 1,
           linewidth = .25) +
  coord_flip(ylim = c(0,60),
             xlim = c(.5, 4.5),
             expand = F) +
  scale_x_discrete(labels = ~paste(stringr::str_to_title(.),"\nObservation \nHeight")) +
  scale_fill_manual(breaks = c("Sig Shift Down",
                               "Shift Down",
                               "No Change",
                               "Shift Up",
                               "Sig Shift Up"),
                    labels = c("Significant\nDownward Shift",
                               "Downward\nTrend",
                               "No\nChange",
                               "Upward\nTrend",
                               "Significant\nUpward Shift"),
                    values = c("darkmagenta",
                               colorspace::lighten("darkmagenta",.8),
                               "grey98",
                               colorspace::lighten("darkturquoise",.8),
                               "darkturquoise")) +
  labs(x = NULL,
       y = "Number of Species",
       fill = "Modeled Shift",
       title = "Species' depth shifts") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(theme= theme(legend.title.position = "top"))) +
  theme(legend.margin = margin(t = -15,0,0,0,"pt"),
        plot.title.position = "plot",
        legend.text = element_text(size = 7,
                                   margin = margin(l = 1,1,1,1,"mm")),
        legend.key.width  = unit(.45,"cm"),
        legend.key.spacing.x = unit(.1,"cm"),
        legend.justification = c(1,.5))


all_data_depth_bars + ggview::canvas(4,4)
ggsave(all_data_depth_bars,
       file = "outputs/depth_shifts/all_data_bars.png",
       width = 4, height = 4)


# find richness trends as functions of maximum height ---------------------
# first, find species by maximum upper tidal height
org_by_yeargroup_sum %>%
  group_by(organism) %>%
  summarize(max_max = mean(height95)) %>%
  left_join(all_p %>% select(organism, max_p, max_trend, max_sign, max_p_sig)) %>%
  ggplot(aes(x = max_max, 
             y = max_trend)) +
  geom_point(aes(color = max_p_sig)) +
  geom_smooth(method = "lm")

org_by_yeargroup_sum %>%
  group_by(organism) %>%
  summarize(mean_mean = max(heightmean)) %>%
  left_join(all_p %>% select(organism, mean_p, mean_trend, mean_sign, mean_p_sig)) %>%
  ggplot(aes(x = mean_mean, 
             y = mean_trend)) +
  geom_point(aes(color = mean_p_sig)) +
  geom_smooth(method = "lm") +
  theme(legend.text = element_markdown())

org_by_yeargroup_sum %>%
  group_by(organism) %>%
  summarize(mean_mean = max(heightmean)) %>%
  left_join(all_p %>% select(organism, mean_p, mean_trend, mean_sign, mean_p_sig)) %>%
  ggplot(aes(x = mean_mean, 
             y = mean_trend)) +
  geom_point(aes(color = mean_p_sig)) +
  geom_smooth(method = "lm") +
  theme(legend.text = element_markdown())


# find if all species higher limits shift downwards
# (starting with all species present in the initial cohort)
initial_cohord <- org_by_yeargroup_sum %>%
  filter(year_group == 1) %>% distinct(organism) %>% pull(organism)
org_by_yeargroup_sum %>%
  filter(organism %in% initial_cohord) %>%
  ggplot(aes(x = year_group,
             y = height95)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_point(aes(y = heightmean)) +
  geom_smooth(aes(y = heightmean), method = "lm") +
  geom_point(aes(y = height05)) +
  geom_smooth(aes(y = height05), method = "lm")


# repeat all highly sampled -----------------------------------------------
pa_hs <- pa_therm_filtered%>%
  tidyr::separate(organism, into = c("Genus","Species")) %>%
  mutate(G = paste0(stringr::str_sub(Genus,1,1),".")) %>%
  mutate(organism = paste(G, Species)) %>%
  select(-G, -Genus, -Species)


hs_by_yeargroup <- pa_hs %>% 
  #filter(replicate == 1) %>%
  mutate(tidalheight = (13-level)*.348)  %>%
  mutate(year_group = case_when(
    year %in% y1 ~ 1,
    year %in% y2 ~ 2,
    year %in% y3 ~ 3,
    year %in% y4 ~ 4,
    year %in% y5 ~ 5,
    year %in% y6 ~ 6,
    year %in% y7 ~ 7,
    year %in% y8 ~ 8
  ))

hs_by_yeargroup_sum <- hs_by_yeargroup %>%
  group_by(organism, year_group) %>%
  
  summarize(height95 = quantile(tidalheight,.95),
            heightmed = median(tidalheight),
            heightmean = mean(tidalheight),
            height05 = quantile(tidalheight,.05)) %>%
  
  # filter to organisms present in 3+ year groups
  group_by(organism) %>%
  mutate(n_yeargroups = n_distinct(year_group)) %>%
  filter(n_yeargroups > 3) %>%
  ungroup() %>%
  select(-n_yeargroups)

hs_by_yeargroup_sum_split <- hs_by_yeargroup_sum %>% split(f = .$organism)

# model and assess mean observation heights
hs_mean_mods <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "heightmean ~ year_group")
)
hs_mean_mods_p <- purrr::map(
  .x = hs_mean_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(mean_p = value)
hs_mean_mods_trend <- purrr::map(
  .x = hs_mean_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                        "(+)",
                                        "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(mean_trend = trend_value,
         mean_sign = trend_sign)
hs_mean_perfect <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~ifelse(length(unique(.x$heightmean)) == 1,"single value","no issue")
)
hs_mean_mods_p[hs_mean_perfect == "single value","mean_p"] <- 1
hs_mean_mods_trend[hs_mean_perfect == "single value","mean_value"] <- 0
hs_mean_mods_trend[hs_mean_perfect == "single value","mean_sign"] <- "None"


# model and assess .05 quantile observation heights
hs_min_mods <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "height05 ~ year_group")
)
hs_min_mods_p <- purrr::map(
  .x = hs_min_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(min_p = value)
hs_min_mods_trend <- purrr::map(
  .x = hs_min_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                        "(+)",
                                        "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(min_trend = trend_value,
         min_sign = trend_sign)
hs_min_perfect <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~ifelse(length(unique(.x$height05)) == 1,"single value","no issue")
)
hs_min_mods_p[hs_min_perfect == "single value","min_p"] <- 1
hs_min_mods_trend[hs_min_perfect == "single value","min_trend"] <- 0
hs_min_mods_trend[hs_min_perfect == "single value","min_sign"] <- "None"


# model and assess .95 quantile observation heights
hs_max_mods <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "height95 ~ year_group")
)
hs_max_mods_p <- purrr::map(
  .x = hs_max_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(max_p = value)
hs_max_mods_trend <- purrr::map(
  .x = hs_max_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                        "(+)",
                                        "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(max_trend = trend_value,
         max_sign = trend_sign)
hs_max_perfect <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~ifelse(length(unique(.x$height95)) == 1,"single value","no issue")
)
hs_max_mods_p[hs_max_perfect == "single value","max_p"] <- 1
hs_max_mods_trend[hs_max_perfect == "single value","max_trend"] <- 0
hs_max_mods_trend[hs_max_perfect == "single value","max_sign"] <- "None"

# model and assess median  observation heights
hs_med_mods <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~lm(data = .x, 
           formula = "heightmed ~ year_group")
)
hs_med_mods_p <- purrr::map(
  .x = hs_med_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(med_p = value)
hs_med_mods_trend <- purrr::map(
  .x = hs_med_mods,
  .f = ~ data.frame(trend_value = summary(.x)$coefficients["year_group","Estimate"],
                    trend_sign = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                        "(+)",
                                        "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(med_trend = trend_value,
         med_sign = trend_sign)
hs_med_perfect <- purrr::map(
  .x = hs_by_yeargroup_sum_split,
  .f = ~ifelse(length(unique(.x$heightmed)) == 1,"single value","no issue")
)
hs_med_mods_p[hs_med_perfect == "single value","med_p"] <- 1
hs_med_mods_trend[hs_med_perfect == "single value","med_trend"] <- 0
hs_med_mods_trend[hs_med_perfect == "single value","med_sign"] <- "None"

# merge all p values ------------------------------------------------------
hs_all_p <- hs_mean_mods_p %>%
  left_join(hs_mean_mods_trend) %>%
  left_join(hs_max_mods_p) %>%
  left_join(hs_max_mods_trend) %>%
  left_join(hs_min_mods_p) %>%
  left_join(hs_min_mods_trend) %>%
  left_join(hs_med_mods_p) %>%
  left_join(hs_med_mods_trend)


hs_all_p <- hs_all_p %>%
  mutate(mean_p_sig = case_when(mean_p < .05 ~ paste0("<b>mean ",mean_sign,"</b><br>"),
                                TRUE ~ ""),
         max_p_sig = case_when(max_p < .05 ~ paste0("<b>max ",max_sign,"</b><br>"),
                               TRUE ~ ""),
         min_p_sig = case_when(min_p < .05 ~ paste0("<b>min ",min_sign,"</b><br>"),
                               TRUE ~ ""),
         med_p_sig = case_when(med_p < .05 ~ paste0("<b>med ",med_sign,"</b><br>"),
                               TRUE ~ "")) %>%
  mutate(any_p = paste0(max_p_sig, mean_p_sig, med_p_sig, min_p_sig))

# plot all ----------------------------------------------------------------
library(ggtext)
hs_trends_p <- hs_by_yeargroup %>%
  filter(organism %in% names(hs_by_yeargroup_sum_split)) %>%
  ggplot(aes(x = as.character(year_group),
             y = tidalheight)) +
  geom_boxplot(outliers = F,
               color = "grey50",
               alpha = .5,
               linewidth = .2,
               width = .3) +
  facet_wrap(~organism) +
  geom_point(data = hs_by_yeargroup_sum,
             aes(y = height95), color = "darkturquoise",
             size = .25) +
  geom_smooth(data = hs_by_yeargroup_sum,
              method= "lm",
              aes(x = year_group,y = height95), 
              color = "darkturquoise",
              se = F,
              size = .25) +
  geom_point(data = hs_by_yeargroup_sum,
             aes(y = height05), color = "darkturquoise",
             size = .25) +
  geom_smooth(data = hs_by_yeargroup_sum,
              method= "lm",
              aes(x = year_group,y = height05), 
              color = "darkturquoise",
              se = F,
              size = .25) +
  geom_point(data = hs_by_yeargroup_sum,
             aes(y = heightmean), color = "darkturquoise",
             size = .25) +
  geom_smooth(data = hs_by_yeargroup_sum,
              method= "lm",
              aes(x = year_group,y = heightmean), 
              color = "darkturquoise",
              se = F,
              size = .25) +
  
  # add max annotation
  geom_richtext(data = hs_all_p,
                aes(x = .5, y = -.6,
                    label = max_p_sig),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =.5,
                hjust = 0) +
  
  geom_richtext(data = hs_all_p,
                aes(x = 4.5, y = -.6,
                    label = med_p_sig),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =.5,
                hjust = .5) +
  
  geom_richtext(data = hs_all_p,
                aes(x = 8.5, y = -.6,
                    label = min_p_sig),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =.5,
                hjust = 1) +
  
  # add bounds
  geom_hline(yintercept = max(hs_by_yeargroup$tidalheight),
             linewidth = .25, linetype = "dashed") +
  geom_hline(yintercept = min(hs_by_yeargroup$tidalheight),
             linewidth = .25, linetype = "dashed") +
  
  theme(strip.text = element_text(margin = margin(0,0,0,0, "pt")),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  scale_x_discrete(labels = c("1982",#-86",
                              "1987",#-91",
                              "1992",#-96",
                              "1997",#-01",
                              "2002",#-06",
                              "2008",#-12",
                              "2013",#-17",
                              "2019"#-23"
  )) +
  labs(x = "Year Group",
       y = "Tidal Level (m)",
       title = "Species Depth Shifts - Highly Sampled Data")
  
hs_trends_p + ggview::canvas(12,8)
ggsave(hs_trends_p, 
       file = "outputs/depth_shifts/hs_trends_p.png",
       width = 12,
       height = 8)


hs_all_p %>% select(-ends_with("sig"), -any_p) %>%
  tidyr::pivot_longer(cols = c(mean_p, mean_trend, 
                               max_p, max_trend,
                               min_p, min_trend,
                               med_p, med_trend),
                      names_to = c("stat", ".value"),
                      names_sep = "_") %>%
  mutate(sig = case_when(p < .05 & trend == "(+)" ~ "Shift Up",
                         p < .05 & trend == "(-)" ~ "Shift Down",
                         p > .05 ~ "No Change",
                         is.na(p) ~ "No Change")) %>%
  mutate(sig = factor(sig, levels = c("Shift Up","No Change","Shift Down"))) %>%
  mutate(stat = factor(stat, levels = c("min","med","mean","max"))) %>%
  ggplot(aes(x = stat,
             fill = sig)) +
  geom_bar() +
  coord_flip() +
  scale_x_discrete(labels = ~paste(stringr::str_to_title(.),"\nObservation \nHeight")) +
  scale_fill_manual(values = c("tomato2","grey85","cyan4")) +
  labs(x = NULL,
       y = "Count",
       fill = "Modeled Shift") +
  theme(legend.position = "top")

hs_all_p2 <- hs_all_p %>% select(-ends_with("sig"), -any_p) %>%
  tidyr::pivot_longer(cols = c(mean_p, mean_sign, 
                               max_p, max_sign,
                               min_p, min_sign,
                               med_p, med_sign),
                      names_to = c("stat", ".value"),
                      names_sep = "_") %>%
  mutate(sig = case_when(p < .05 & sign == "(+)" ~ "Sig Shift Up",
                         p > .05 & sign == "(+)" ~ "Shift Up",
                         p < .05 & sign == "(-)" ~ "Sig Shift Down",
                         p > .05 & sign == "(-)" ~ "Shift Down",
                         sign == "None" ~ "No Change")) %>%
  mutate(sig = factor(sig, levels = c("Sig Shift Up","Shift Up","No Change","Shift Down","Sig Shift Down"))) %>%
  mutate(stat = factor(stat, levels = c("min","med","mean","max"))) 

  

hs_depth_bars <- hs_all_p2 %>%
  ggplot(aes(x = stat,
             fill = sig)) +
  geom_bar(color = "black", width = 1,
           linewidth = .2) +
  coord_flip(ylim = c(0,49),
             xlim = c(.5, 4.5),
             expand = F) +
  scale_x_discrete(labels = ~paste(stringr::str_to_title(.),"\nObservation \nHeight")) +
  scale_fill_manual(breaks = c("Sig Shift Down",
                               "Shift Down",
                               "No Change",
                               "Shift Up",
                               "Sig Shift Up"),
                    values = c("darkmagenta",
                               colorspace::lighten("darkmagenta",.8),
                               "grey95",
                               colorspace::lighten("darkturquoise",.8),
                               "darkturquoise")) +
  labs(x = NULL,
       y = "Number of Species",
       fill = "Modeled Shift") +
  theme(legend.position = "bottom")  +
  geom_hline(yintercept = 49/2,
             linetype = "dashed") +
  guides(fill = guide_legend(theme= theme(legend.title.position = "top"))) +
  theme(legend.margin = margin(t = -15,0,0,0,"pt"))

hs_depth_bars + ggview::canvas(6.2,5)
ggsave(hs_depth_bars,
       file = "outputs/depth_shifts/hs_bars.png",
       width = 6.2, height = 5)


# find richness trends as functions of maximum height ---------------------
# first, find species by maximum upper tidal height
hs_by_yeargroup_sum %>%
  group_by(organism) %>%
  summarize(mean_max = mean(height95)) %>%
  left_join(hs_all_p %>% select(organism, contains("max"))) %>%
  ggplot(aes(x = mean_max, 
             y = max_trend)) +
  geom_point(aes(color = max_p_sig)) +
  geom_smooth(method = "lm")

hs_by_yeargroup_sum %>%
  group_by(organism) %>%
  summarize(mean_mean = mean(heightmean)) %>%
  left_join(hs_all_p %>% select(organism, contains("mean"))) %>%
  ggplot(aes(x = mean_mean, 
             y = mean_trend)) +
  geom_point(aes(color = mean_p_sig)) +
  geom_smooth(method = "lm") +
  theme(legend.text = element_markdown())

hs_by_yeargroup_sum %>%
  group_by(organism) %>%
  summarize(mean_min = max(height05)) %>%
  left_join(hs_all_p %>% select(organism,contains("min"))) %>%
  ggplot(aes(x = mean_min, 
             y = min_trend)) +
  geom_point(aes(color = min_p_sig)) +
  geom_smooth(method = "lm") +
  theme(legend.text = element_markdown())


# find if all species higher limits shift downwards
# (starting with all species present in the initial cohort)
hs_initial_cohort <- hs_by_yeargroup_sum %>%
  filter(year_group == 1) %>% distinct(organism) %>% pull(organism)
hs_by_yeargroup_sum %>%
  filter(organism %in% hs_initial_cohort) %>%
  ggplot(aes(x = year_group,
             y = height95)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_point(aes(y = heightmean)) +
  geom_smooth(aes(y = heightmean), method = "lm") +
  geom_point(aes(y = height05)) +
  geom_smooth(aes(y = height05), method = "lm")




# try without clustering --------------------------------------------------

unclustered <- hs_by_yeargroup %>%
  group_by(organism, year) %>%
  summarize(min = quantile(tidalheight, .05),
            med = quantile(tidalheight, .5),
            max = quantile(tidalheight,.95)) %>%
  group_by(organism) %>%
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years > 10)  %>%
  split(f = .$organism)

unclustered_max_mods <- purrr::map(
  .x = unclustered,
  .f = ~lm(data = .x,
           formula = "max ~ year")
)
unclustered_max_p <- purrr::map(
  .x = unclustered_max_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(med_p = value)
hs_med_mods_trend <- purrr::map(
  .x = hs_med_mods,
  .f = ~ data.frame(value = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                   "(+)",
                                   "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(med_trend = value)

# test for change in sampling over time -----------------------------------
samples_by_yeargroup <- org_by_yeargroup %>%
  distinct(year_group, tidalheight) %>%
  group_by(year_group) %>%
  
  summarize(height95 = quantile(tidalheight,.95),
            heightmed = median(tidalheight),
            heightmean = mean(tidalheight),
            height05 = quantile(tidalheight,.05)) 

year_max_mod <- lm(data = samples_by_yeargroup,
                   formula = "height95 ~year_group")
summary(year_max_mod)
year_mean_mod <- lm(data = samples_by_yeargroup,
                   formula = "heightmean ~year_group")
summary(year_mean_mod)
year_med_mod <- lm(data = samples_by_yeargroup,
                    formula = "heightmed ~year_group")
summary(year_med_mod)
year_min_mod <- lm(data = samples_by_yeargroup,
                   formula = "height05 ~year_group")
summary(year_min_mod)

org_by_yeargroup %>%
  ggplot(aes(x = as.character(year_group),
             y = tidalheight)) +
  geom_boxplot(outliers = F)

org_by_yeargroup %>% 
  group_by(year_group, tidalheight) %>%
  summarize(n_obs = n(),
            n_tran = n_distinct(transect)) %>% 
  ggplot(aes(x = year_group,
             y = tidalheight)) +
  geom_tile(aes(fill = n_tran ))

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------



mean_mods <- purrr::map(
  .x = org_by_height,
  .f = ~lm(data = .x, 
           formula = "heightmean ~ year_group")
)
mean_mods_p <- purrr::map(
  .x = mean_height_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
)

mean_height_mods$`Crepidula fornicata` %>% summary()

mod <- lm(data = org_by_height[[1]],
          formula = "heightmean ~ year_group")
summary(mod)

  ggplot(aes(x = year_group,
             y = heightmean)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "blue",
              se = F) +
  geom_point(aes(y = height95), color = "red") +
  geom_smooth(method = "lm", se = F, aes(y = height95),
              color = "red") +
  geom_point(aes(y = height05), color = "green") +
  geom_smooth(method = "lm", se = F, aes(y = height05),
              color = "green") +
  facet_wrap(~organism)
  





# -------------------------------------------------------------------------
# Scratch / old ---------------------------------------------------------------------
# -------------------------------------------------------------------------





therm <- readr::read_csv(
  here::here(
    "data-processed",
    "spp_thermal_affinities.csv"
  )
)

pa <- pa %>% left_join(therm) %>%
  filter(pres == T) %>%
  # translate level into tidal height
  mutate(tidalheight = (13-level)*.348)  
pa %>% distinct(data_taken)


# find cti of new spp -----------------------------------------------------
# first, let's break up the timeseries 
range(pa$year)
pa %>% 
  group_by(organism) %>%
  filter(year == min(year)) %>% 
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_group = case_when(
    year <= 1989 ~ "1980s",
    year >= 1990 & year <= 1999 ~ "1990s",
    year >= 2000 & year <= 2009 ~ "2000s",
    year >= 2010 & year <= 2019 ~ "2010s",
    year >= 2020 ~ "2020s"
  )) %>%
  ggplot(aes(x = year,
             y = mean_monthly_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year of First Occurrence",
       y = "Thermal Affinity") 



pa %>% 
  group_by(organism) %>%
  filter(year == max(year)) %>% 
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_group = case_when(
    year <= 1989 ~ "1980s",
    year >= 1990 & year <= 1999 ~ "1990s",
    year >= 2000 & year <= 2009 ~ "2000s",
    year >= 2010 & year <= 2019 ~ "2010s",
    year >= 2020 ~ "2020s"
  )) %>%
  filter(year != max(year)) %>%
  ggplot(aes(x = year,
             y = mean_monthly_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year of First Occurrence",
       y = "Thermal Affinity") 




# find how long species lasted in the dataset -----------------------------
first_year <- pa %>%
  group_by(organism) %>%
  filter(year == min(year)) %>%
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_type = "first_year")

last_year <- pa %>%
  group_by(organism) %>%
  filter(year == max(year)) %>%
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_type = "last_year")

year_breadth <-  
  first_year %>%
  rbind(last_year)


year_breadth %>% 
  tidyr::pivot_wider(names_from = year_type,
                     values_from = year) %>% 
  filter(last_year - first_year > 0) %>% 
  ggplot() +
  geom_segment(aes(x = first_year,
                   xend = last_year,
                   y = mean_monthly_mean,
                   yend = mean_monthly_mean),
               linewidth = .3,
               alpha = .4) +
  geom_point(data = . %>% filter(first_year > min(first_year)),
             aes(x = first_year,
                 y = mean_monthly_mean),
             color = "cyan4") +
  geom_point(data = . %>% filter(last_year < max(last_year)),
             aes(x = last_year,
                 y = mean_monthly_mean),
             color = "tomato2") +
  geom_smooth(data = . %>% filter(first_year > min(first_year)),
              aes(x = first_year,
                  y = mean_monthly_mean),
              color = "cyan4",
              method = "lm") +
  geom_smooth(data = . %>% filter(last_year < max(last_year)),
              aes(x = last_year,
                  y = mean_monthly_mean),
              color = "tomato2",
              method = "lm")

pa %>%
  distinct(transect, level, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = level)) +
  geom_tile() +
  facet_wrap(~transect)
# longest most consistent transects are:
# 5, 7, 15, 20, 22, 26, 28
# between levels 0 and 15
good_transects <- c(5, 7, 15, 20, 22, 26, 28)
min_level <- 0
max_level <- 15

# track individual species min/max depths ---------------------------------

pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>% 
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  ggplot(aes(x = year,
             group = organism)) +
  geom_line(aes(y = minheight),
             color = "grey80",
            linewidth = .2,
            alpha = .3) +
  geom_smooth(aes(y = minheight,
                  group = organism),
              method = "lm",
              se = F,
              color = "cyan4",
              linetype = "dashed",
              linewidth = .5) 

pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>% 
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight),
            mean_monthly_mean = unique(mean_monthly_mean)) %>%
  group_by(organism) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(mean_monthly_mean),
         n > 3) %>%
  mutate(organism_f = forcats::fct_reorder(organism, mean_monthly_mean)) %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = minheight),
            color = "tomato2",
            linewidth = .9,
            alpha = .3) +
  geom_line(aes(y = maxheight),
            color = "cyan4",
            linewidth = .9,
            alpha = .3) +
# geom_smooth(aes(y = maxheight),
#             method = "lm",
#             se = F,
#             color = "tomato2",
#             linetype = "dashed",
#             linewidth = .5) +
  facet_wrap(~organism_f)


pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>%   
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  ggplot(aes(x = year,
             group = organism)) +
  geom_ribbon(aes(ymin = minheight,
                  ymax = maxheight)) +
  facet_wrap(~organism)


pa %>%
  filter(!is.na(mean_monthly_mean)) %>%
  mutate(organism_f = forcats::fct_reorder(organism, mean_monthly_mean)) %>%
  filter(transect %in% good_transects,
        level >= min_level,
        level <= max_level) %>% 
  group_by(organism, organism_f, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  group_by(organism) %>%
  add_count() %>%
  ungroup() %>%
  filter(n > 3) %>%
  ggplot(aes(x = year)) +
  geom_smooth(aes(y = maxheight),
              method = "lm", 
              se = F,
              color = "tomato2") +
  geom_smooth(aes(y = minheight),
              method = "lm", 
              se = F,
              color = "cyan4") +
  facet_wrap(~organism_f)
