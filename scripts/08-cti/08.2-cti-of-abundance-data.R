# find mean thermal affinity over time overall
# by ABUNDANCE


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(lme4)
theme_set(ggthemes::theme_few())



# upload data --------------------------------------------------------------
# get raw cover data
cover <- readr::read_csv("data-processed/appledore-survey-data/cover_abundance_filtered.csv")
cover %>% glimpse()
cover <- cover %>% 
  select(organism, year, transect, level, replicate, percent_cover,pc_num) %>%
  filter(!is.na(pc_num))
# get STIs
therm <- readr::read_csv("data-processed/species-thermal-affinities/spp_thermal_affinities.csv")
therm %>% glimpse()
# get highly sampled quadrats
hs <- readr::read_csv("data-processed/cti-data/cti-data-highly-sampled.csv")
hs <- hs %>% distinct(year, transect, level)
temptrend <- readr::read_csv("data-processed/cti-data/temp-change-augment.csv")
tempmod <- readRDS("data-processed/cti-data/temp-change-model.rds")
tempmod %>% summary()

df <- cover %>%
  #filter(replicate == 1) %>%
  inner_join(hs) %>%
  left_join(therm %>%
              select(organism,
                     perc_monthly_min_05, 
                     mean_monthly_mean, 
                     perc_monthly_max_95)) %>%
  tidyr::drop_na(mean_monthly_mean) %>%
  filter(pc_num != 0)


# first, raw cti per year
df_sum <- df  %>%
  group_by(year) %>%
  summarise(cti = weighted.mean(mean_monthly_mean,
                                pc_num, na.rm=T),
            cti.max = weighted.mean(perc_monthly_max_95,
                                    pc_num, na.rm=T),
            cti.min = weighted.mean(perc_monthly_min_05,
                                    pc_num, na.rm=T)) 

df_sum %>%
  ggplot(aes(x = year, 
                 y = cti)) +
  geom_point() +
  geom_smooth(method = "lm")

df_sum %>%
  ggplot(aes(x = year, 
             y = cti.max)) +
  geom_point() +
  geom_smooth(method = "lm")

df_sum %>%
  ggplot(aes(x = year, 
             y = cti.min)) +
  geom_point() +
  geom_smooth(method = "lm")

df_mod <- lm(data = df_sum,
             cti ~ year)
summary(df_mod)


# repeat with counts ------------------------------------------------------
counts <- readr::read_csv("data-processed/appledore-survey-data/counts_abundance_filtered.csv")
counts %>% glimpse()
counts <- counts %>% 
  select(organism, year, transect, level, replicate, count, countnum) %>%
  filter(!is.na(countnum))

df_count <- counts %>%
#  filter(replicate == 1) %>%
  inner_join(hs) %>%
  left_join(therm %>%
              select(organism,
                     perc_monthly_min_05, 
                     mean_monthly_mean, 
                     perc_monthly_max_95)) %>%
  tidyr::drop_na(mean_monthly_mean) %>%
  filter(countnum != 0)

# first, raw cti per year
df_count_sum <- df_count  %>%
  group_by(year) %>%
  summarise(cti = weighted.mean(mean_monthly_mean,
                                countnum, na.rm=T),
            cti.max = weighted.mean(perc_monthly_max_95,
                                    countnum, na.rm=T),
            cti.min = weighted.mean(perc_monthly_min_05,
                                    countnum, na.rm=T)) 

df_count_sum %>%
  ggplot(aes(x = year, 
             y = cti)) +
  geom_point() +
  geom_smooth(method = "lm")

df_count_sum %>%
  ggplot(aes(x = year, 
             y = cti.max)) +
  geom_point() +
  geom_smooth(method = "lm")

df_count_sum %>%
  ggplot(aes(x = year, 
             y = cti.min)) +
  geom_point() +
  geom_smooth(method = "lm")

df_count_mod <- lm(data = df_count_sum,
             cti ~ year)
summary(df_count_mod)


# merge the two together --------------------------------------------------

df_merge <- df_count_sum %>% mutate(df = "count") %>%
  rbind(df_sum %>% mutate(df = "cover")) 

df_merge %>% 
  ggplot(aes(x = year,
             y = cti)) +
  geom_point(aes(color = df)) +
  geom_smooth(method = "lm")

df_merge_mod1 <- lm(data = df_merge,
                   cti ~ year)
#plot(df_merge_mod1)
summary(df_merge_mod1)
df_merge_mod2 <- lm(data = df_merge,
                    cti ~ year + df)
#plot(df_merge_mod1)
summary(df_merge_mod2)
df_merge_mod3 <- lm(data = df_merge,
                    cti ~ year * df)
#plot(df_merge_mod1)
summary(df_merge_mod3)

performance::compare_performance(df_merge_mod1,
                                 df_merge_mod2, 
                                 df_merge_mod3, rank = T)
summary(df_merge_mod2)
# 
library(performance)
perf <- plot(check_model(df_merge_mod2))
perf[[5]] <- perf[[5]] + coord_cartesian(xlim = c(.5, 2.5), expand = F)
library(patchwork)
perf <- perf & labs(subtitle = NULL)

perf + ggview::canvas(12,8)
ggsave(perf, 
       filename = "outputs/cti/total_cti_assumptions.png",
       width = 12,
       height = 8)

df_merge %>% filter(is.na(df))

au <- broom::augment(df_merge_mod2,
                     df_merge,
                     interval = "confidence")

cti_abund <- au %>%
  ggplot(aes(x = year,
             y = cti,
             ymin = .lower,
             ymax = .upper,
             group = df,
             color = df, 
             fill = df
             )) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .3) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed",
            alpha = .5) +
  
  # add CTI ribbons
  geom_ribbon(alpha = .5, 
              color = "transparent",
              show.legend = F) +
  geom_line(aes(y = .fitted),
            show.legend = F) +
  geom_point(aes(color = df)) +
  
  # scales
  scale_color_manual(values = c("cyan3","magenta3"),
                     breaks = c("cover","count"),
                     labels = c("Percent Cover Species","Density Species")) +
  scale_fill_manual(values = c("cyan","magenta3"),
                    breaks = c("cover","count"),
                    labels = c("Percent Cover Species","Density Species")) +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  
  # fix legend
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1),
        legend.title = element_blank()) +
  
  # add model outputs
  annotate(geom = "text",
           x = 2023,
           y = 10.3,
           hjust = 1, 
           vjust = 0,
           label = "Coefficient: 0.05°/Decade\np < 0.05") +
  coord_cartesian(xlim = c(1981, 2024),
                  ylim = c(10.25, 12.25),
                  expand = F) +
  labs(x = NULL,
       y = 'Community Thermal Index (°C)')


cti_abund + ggview::canvas(4.5,4)
saveRDS(cti_abund,
        file = "outputs/cti/cti_by_abundance_plot.rds")
ggsave(cti_abund,
       filename = "outputs/cti/cti_by_abundance_plot.png",
       width = 4.5, height = 4, unit = "in")


au %>% 
  filter(year %in% c(min(year),max(year))) 

# try CTImax --------------------------------------------------------------
raw_temps <- readr::read_csv("data-processed/appledore-island-env-data/appledore_temps_1982_2023.csv")
tempmod.max <- lm(data = raw_temps,
                  max ~ sample_year)
tempmod.max.au <- broom::augment(tempmod.max,
                                 interval = "confidence")

df_merge_mod1.max <- lm(data = df_merge,
                    cti.max ~ year)
#plot(df_merge_mod1)
summary(df_merge_mod1.max)
df_merge_mod2.max <- lm(data = df_merge,
                    cti.max ~ year + df)
#plot(df_merge_mod1)
summary(df_merge_mod2.max)
df_merge_mod3.max <- lm(data = df_merge,
                    cti.max ~ year * df)
#plot(df_merge_mod1)
summary(df_merge_mod3.max)

performance::compare_performance(df_merge_mod1.max,
                                 df_merge_mod2.max, 
                                 df_merge_mod3.max, rank = T)
summary(df_merge_mod2.max)



au.max <- broom::augment(df_merge_mod2.max,
                     df_merge,
                     interval = "confidence")

cti_abund.max <- au.max %>%
  ggplot(aes(x = year,
             y = cti.max,
             ymin = .lower,
             ymax = .upper,
             group = df,
             color = df, 
             fill = df
  )) +
  
  # first, add real temperature trend
  geom_ribbon(data = tempmod.max.au, inherit.aes =F,
              aes(x = sample_year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .3) +
  geom_line(data = tempmod.max.au,  inherit.aes = F,
            aes(x = sample_year, y = .fitted),
            linetype = "dashed",
            alpha = .5) +
  
  # add CTI ribbons
  geom_ribbon(alpha = .5, 
              color = "transparent",
              show.legend = F) +
  geom_line(aes(y = .fitted),
            show.legend = F) +
  geom_point(aes(color = df)) +
  
  # scales
  scale_color_manual(values = c("cyan3","magenta3"),
                     breaks = c("cover","count"),
                     labels = c("Percent Cover Species","Density Species")) +
  scale_fill_manual(values = c("cyan","magenta3"),
                    breaks = c("cover","count"),
                    labels = c("Percent Cover Species","Density Species")) +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  
  # fix legend
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1),
        legend.title = element_blank()) +
  
  # add model outputs
  annotate(geom = "text",
           x = 2023,
           y = 17.7,
           hjust = 1, 
           vjust = 0,
           label = "Coefficient: 0.005°/Decade\np = 0.86") +
  coord_cartesian(xlim = c(1981, 2024),
                  ylim = c(17.6, 21),
                  expand = F) +
  labs(x = NULL,
       y = 'Community Thermal Index (°C)')


cti_abund.max + ggview::canvas(4.5,4)
#saveRDS(cti_abund,
#        file = "outputs/cti/cti_by_abundance_plot.rds")
#ggsave(cti_abund,
#       filename = "outputs/cti/cti_by_abundance_plot.png",
#       width = 4.5, height = 4, unit = "in")
#



# try with CTI.min --------------------------------------------------------

tempmod.min <- lm(data = raw_temps,
                  min ~ sample_year)
summary(tempmod.min)
tempmod.min.au <- broom::augment(tempmod.min,
                                 interval = "confidence")

df_merge_mod1.min <- lm(data = df_merge,
                        cti.min ~ year)
#plot(df_merge_mod1)
summary(df_merge_mod1.min)
df_merge_mod2.min <- lm(data = df_merge,
                        cti.min ~ year + df)
#plot(df_merge_mod1)
summary(df_merge_mod2.min)
df_merge_mod3.min <- lm(data = df_merge,
                        cti.min ~ year * df)
#plot(df_merge_mod1)
summary(df_merge_mod3.min)

performance::compare_performance(df_merge_mod1.min,
                                 df_merge_mod2.min, 
                                 df_merge_mod3.min, rank = T)
summary(df_merge_mod2.min)



au.min <- broom::augment(df_merge_mod2.min,
                         df_merge,
                         interval = "confidence")

cti_abund.min <- au.min %>%
  ggplot(aes(x = year,
             y = cti.min,
             ymin = .lower,
             ymax = .upper,
             group = df,
             color = df, 
             fill = df
  )) +
  
  # first, add real temperature trend
  geom_ribbon(data = tempmod.min.au, inherit.aes =F,
              aes(x = sample_year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .3) +
  geom_line(data = tempmod.min.au,  inherit.aes = F,
            aes(x = sample_year, y = .fitted),
            linetype = "dashed",
            alpha = .5) +
  
  # add CTI ribbons
  geom_ribbon(alpha = .5, 
              color = "transparent",
              show.legend = F) +
  geom_line(aes(y = .fitted),
            show.legend = F) +
  geom_point(aes(color = df)) +
  
  # scales
  scale_color_manual(values = c("cyan3","magenta3"),
                     breaks = c("cover","count"),
                     labels = c("Percent Cover Species","Density Species")) +
  scale_fill_manual(values = c("cyan","magenta3"),
                    breaks = c("cover","count"),
                    labels = c("Percent Cover Species","Density Species")) +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  
  # fix legend
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1),
        legend.title = element_blank()) +
  
  # add model outputs
  annotate(geom = "text",
           x = 2023,
           y = 0.27,
           hjust = 1, 
           vjust = 0,
           label = "Coefficient: 0.19°/Decade\np < 0.01") +
  coord_cartesian(xlim = c(1981, 2024),
                  ylim = c(0.15, 5),
                  expand = F) +
  labs(x = NULL,
       y = 'Community Thermal Index (°C)')


cti_abund.min + ggview::canvas(4.5,4)


cti_all <- ((cti_abund.max + labs(title  =expression("a) "*CTI[max]*" and Max Temp"),
                      y = expression(CTI[max]))|
     cti_abund + labs(title  =expression("b) "*CTI[mean]*" and Mean Temp"),
                      y = expression(CTI[mean]))|
    cti_abund.min + labs(title  =expression("c) "*CTI[min]*" and Min Temp"),
                         y = expression(CTI[min]))) &
  theme(legend.background = element_blank())) 
  
cti_all + ggview::canvas(12,4)

#ggsave(cti_all,
#       filename = "outputs/cti/cti_trends_all_parameters.png",
#       width = 12, height = 4)



# average CTI between two groups ------------------------------------------

df_merge2 <- df_count_sum %>%
  rbind(df_sum) %>%
  group_by(year) %>%
  summarize(cti = mean(cti)) 

df_merge2 %>% 
  ggplot(aes(x = year,
             y = cti)) +
  geom_point() +
  geom_smooth(method = "lm")

df_merge2_mod2 <- lm(data = df_merge2,
                   cti ~ year)

summary(df_merge2_mod2)

performance::compare_performance(df_merge2_mod2,
                                 df_merge_mod2, 
                                 rank = T)


# weight groups by number of species --------------------------------------
count_n_spp <- df_count %>%
  group_by(year) %>%
  summarise(n_spp = n_distinct(organism),
            .groups = "drop") %>%
  mutate(df = "count")
cover_n_spp <- df %>%
  group_by(year) %>%
  summarise(n_spp = n_distinct(organism),
            .groups = "drop") %>%
  mutate(df = "cover")
  
df_merge3 <- df_merge %>%
  left_join(cover_n_spp %>% rbind(count_n_spp)) %>%
  group_by(year) %>%
  summarise(cti = weighted.mean(cti, n_spp))

df_merge3 %>% 
  ggplot(aes(x = year,
             y = cti)) +
  geom_point() +
  geom_smooth(method = "lm")

df_merge3_mod <- lm(data = df_merge3,
                    cti ~ year)

summary(df_merge3_mod)


