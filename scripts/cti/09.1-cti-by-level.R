# find mean thermal affinity over time in levels

# here, I use the presence/absence dataset to find the mean
# thermal affinity of communities within tidal heights to 
# ask whether the communities are becoming more warm-affinity
# over time within the study region


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(lme4)
theme_set(ggthemes::theme_few())



# upload data -------------------------------------------------------------
pa_therm <- readr::read_csv("data-processed/cti-data/cti-data-all.csv")
pa_therm_filtered <- readr::read_csv("data-processed/cti-data/cti-data-highly-sampled.csv")
temptrend <- readr::read_csv("data-processed/cti-data/temp-change-augment.csv")
tempmod <- readRDS("data-processed/cti-data/temp-change-model.rds")
first_temp <- predict(tempmod, list(year= min(pa_therm$year)))
last_temp <- predict(tempmod, list(year= max(pa_therm$year)))



# make function to test models --------------------------------------------
test_mods <- function(df){
  mod_add <- lm(data = df,
                formula = cti ~ year + tidalheight)
  mod_mult <- lm(data = df,
                 formula = cti ~ year * tidalheight)
  cti_add_par <- lm(data = df,
                    formula = cti ~ year + I(tidalheight^2)) 
  cti_mult_par <- lm(data = df,
                     formula = cti ~ year * I(tidalheight^2))
  perf <- performance::compare_performance(
    mod_add,
    mod_mult,
    cti_add_par,
    cti_mult_par,
    rank = T
  )
  
  best_mod <- get(perf$Name[[1]])
  print(paste0("best model:", perf$Name[[1]]))
  return(best_mod)
}

get_mod_p <- function(mod){
  
  full_p <- pf(summary(mod)$fstatistic[1], summary(mod)$fstatistic[2], summary(mod)$fstatistic[3], lower.tail = FALSE)

  mod_p_text <- NA
  if(full_p < 0.01){mod_p_text <- "Model p < 0.01"}
  else if(full_p < 0.05){mod_p_text <- "Model p < 0.05"}
  else if(full_p > 0.05){mod_p_text <- "Model p > 0.05"}
  
  term_p <- summary(mod)$coefficients
  if(any(stringr::str_detect(rownames(term_p),":"))){
   int_term <- term_p[which(stringr::str_detect(rownames(term_p),":")),"Pr(>|t|)"]
   int_term_text <- NA
   if(int_term < 0.01){int_term_text <- "Interaction p < 0.01"}
   else if(int_term < 0.05){int_term_text <- "Interaction p < 0.05"}
   else if(int_term > 0.05){int_term_text <- "Interaction p > 0.05"}
   
  }
  
  out <- list("modp" = mod_p_text,
              "term" = int_term_text)
  
  return(out)
}



# p1 highly sampled, per level ---------------------------------------------
df1 <- pa_therm_filtered %>%
  # find CTI per level/year
  group_by(level, tidalheight, year) %>%
  summarize(cti = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are 4+ species
  #filter(n_spp >= 4) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level))


df1 %>%   ggplot(aes(x = year, y = cti,
                     #color = level, 
                     group = tidalheight)) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight_f) +
  geom_point(alpha = .6) +
  gghighlight::gghighlight(
    unhighlighted_params = list(color = "grey80",
                                alpha = .3),
    use_direct_label = F) 

mod1 <- test_mods(df1)
summary(mod1)
mod1_p <- get_mod_p(mod1)
au1 <- broom::augment(mod1, 
                      newdata = tidyr::complete(df1,year,tidalheight), 
                      interval = "prediction")
p1 <- au1 %>% 
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.2,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
  annotate(geom = "richtext",
           label = paste0("<b>Best Model:</b><br>",
                          mod1$call[[2]], "<br>",
                          mod1_p[[1]],"<br>",
                          mod1_p[[2]]),
           label.color = NA,
           fill = NA,
           x = 10.2,
           y = 2.9,
           size = 3,
           hjust = 0,
           vjust = 1) +
  coord_cartesian(xlim = c(10.05,12.15),
                  ylim = c(-.4,3),
                  expand = F) +

  theme(legend.position = c(.08,.1),
        legend.justification = c(0,0),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       #title = "CTI Across Shore Levels"
       ) +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(), 
        plot.title.position = "plot")

p1 + ggview::canvas(4,4)

# save as rds object to combine later
saveRDS(p1,
        "outputs/cti/bylevel_cti_p1.rds")


# p2 repeat, first replicate only -----------------------------------------
df2 <- pa_therm_filtered %>%
  # find CTI per level/year
  filter(replicate == 1) %>%
  group_by(level, tidalheight, year) %>%
  summarize(cti = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are 4+ species
  #filter(n_spp >= 4) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level))


mod2 <- test_mods(df2)
summary(mod2)
mod2_p <- get_mod_p(mod2)
au2 <- broom::augment(mod2, 
                      newdata = tidyr::complete(df2,year,tidalheight), 
                      interval = "prediction")
p2 <- au2 %>% 
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.2,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
  annotate(geom = "richtext",
           label = paste0("<b>Best Model:</b><br>",
                          mod2$call[[2]], "<br>",
                          mod2_p[[1]],"<br>",
                          mod2_p[[2]]),
           label.color = NA,
           fill = NA,
           x = 10.2,
           y = 2.9,
           size = 3,
           hjust = 0,
           vjust = 1) +
  coord_cartesian(xlim = c(10.05,12.15),
                  ylim = c(-.4,3),
                  expand = F) +
  
  theme(legend.position = c(.08,.1),
        legend.justification = c(0,0),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "HS, first replicates only") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(),
        plot.title.position = "plot")

p2 + ggview::canvas(4,4)



# p3 HS, unsummarized -----------------------------------------------------
df3 <- pa_therm_filtered %>%
  rename(cti = mean_monthly_mean) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level))


mod3 <- test_mods(df3)
summary(mod3)
mod3_p <- get_mod_p(mod3)
au3 <- broom::augment(mod3, df3, interval = "prediction")
p3 <- au3 %>% 
 # arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = last_temp + 0.1, 
           y = -.25,
           label = "True temperature change",
           hjust = 0,
           vjust =0.5,
           size =3) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 1.75,
             alpha = .5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_line(data = . %>% distinct(year, tidalheight, .fitted),
              arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1) +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
  annotate(geom = "richtext",
           label = paste0("<b>Best Model:</b><br>",
                          mod3$call[[2]], "<br>",
                          mod3_p[[1]],"<br>",
                          mod3_p[[2]]),
           label.color = NA,
           fill = NA,
           x = 19,
           y = 2.9,
           size = 3,
           hjust = 1,
           vjust = 1) +
  theme(legend.position = c(.98,.75),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "HS, unsummarized") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(),
        plot.title.position = "plot")

p3 + ggview::canvas(4,4)





# p4 all data, per level ---------------------------------------------
df4 <- pa_therm %>%
  # find CTI per level/year
  group_by(level, tidalheight, year) %>%
  summarize(cti = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are 4+ species
  filter(n_spp >= 4) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level)) %>%
  # filter to levels where there are at least 3 years of sampling
  # (because we can't tell a trend when a level is only sampled once)
  group_by(level) %>%
  add_count() %>%
  filter(n >= 20) %>%
  select(-n) %>% ungroup()

df4 %>%   ggplot(aes(x = year, y = cti,
                     #color = level, 
                     group = tidalheight)) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight_f) +
  geom_point(alpha = .6) +
  gghighlight::gghighlight(
    unhighlighted_params = list(color = "grey80",
                                alpha = .3),
    use_direct_label = F) 

mod4 <- test_mods(df4)
summary(mod4)
mod4_p <- get_mod_p(mod4)
au4 <- broom::augment(mod4, 
                      newdata = tidyr::complete(df4,year,tidalheight),
                      #newdata = tidyr::complete(df4, year, tidalheight), 
                      interval = "prediction")
p4 <- au4 %>% 
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.2,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
  annotate(geom = "richtext",
           label = paste0("<b>Best Model:</b><br>",
                          mod4$call[[2]], "<br>",
                          mod4_p[[1]],"<br>",
                          mod4_p[[2]]),
           label.color = NA,
           fill = NA,
           x = 10.2,
           y = 3.6,
           size = 3,
           hjust = 0,
           vjust = 1) +
  coord_cartesian(xlim = c(10.05,12.15),
                  ylim = c(-.4,3.7),
                  expand = F) +
  
  theme(legend.position = c(.08,.1),
        legend.justification = c(0,0),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "All data, all replicates") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(), 
        plot.title.position = "plot")

p4 + ggview::canvas(4,4)

# p5 all data, first replicates ---------------------------------------------
df5 <- pa_therm %>%
  filter(replicate == 1) %>%
  # find CTI per level/year
  group_by(level, tidalheight, year) %>%
  summarize(cti = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are 4+ species
  filter(n_spp >= 4) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level)) %>%
  # filter to levels where there are at least 3 years of sampling
  # (because we can't tell a trend when a level is only sampled once)
  group_by(level) %>%
  add_count() %>%
  filter(n >= 20) %>%
  select(-n) %>% ungroup()

df5 %>%   ggplot(aes(x = year, y = cti,
                     #color = level, 
                     group = tidalheight)) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight_f) +
  geom_point(alpha = .6) +
  gghighlight::gghighlight(
    unhighlighted_params = list(color = "grey80",
                                alpha = .3),
    use_direct_label = F) 

mod5 <- test_mods(df5)
summary(mod5)
mod5_p <- get_mod_p(mod5)

au5 <- broom::augment(mod5, 
                      newdata = tidyr::complete(df5, year, tidalheight), 
                      interval = "prediction")
p5 <- au5 %>% 
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.2,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
  annotate(geom = "richtext",
           label = paste0("<b>Best Model:</b><br>",
                          mod5$call[[2]], "<br>",
                          mod5_p[[1]],"<br>",
                          mod5_p[[2]]),
           label.color = NA,
           fill = NA,
           x = 10.2,
           y = 3.6,
           size = 3,
           hjust = 0,
           vjust = 1) +
  coord_cartesian(xlim = c(10.05,12.15),
                  ylim = c(-.4,3.7),
                  expand = F) +
  
  theme(legend.position = c(.08,.1),
        legend.justification = c(0,0),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "All data, first replicates only") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(), 
        plot.title.position = "plot")

p5 + ggview::canvas(4,4)


# p6 all, unsummarized -----------------------------------------------------
df6 <- pa_therm %>%
  rename(cti = mean_monthly_mean) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level)) %>%
  # filter to levels where there are at least 3 years of sampling
  # (because we can't tell a trend when a level is only sampled once)
  group_by(level) %>%
  mutate(n_years = n_distinct(year),
         n_spp = n_distinct(organism)) %>%
  filter(n_years > 20) %>%
  filter(n_spp >= 4) %>%
  ungroup() %>%
  select(-n_years,-n_spp)

df6 %>% 
  ggplot(aes(x = year,
             y = cti)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tidalheight_f)

mod6 <- test_mods(df6)
summary(mod6)
mod6_p <- get_mod_p(mod6)
au6 <- broom::augment(mod6,
                      newdata = tidyr::complete(df6,year,tidalheight),
                      interval = "prediction")
p6 <- au6 %>% 
  # arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.75, yend = -.75,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = last_temp + 0.2, 
           y = -.75,
           label = "True temperature change",
           hjust = 0,
           vjust =0.5,
           size =3) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 1.75,
             alpha = .5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_line(data = . %>% distinct(year, tidalheight, .fitted),
            arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1) +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
  annotate(geom = "richtext",
           label = paste0("<b>Best Model:</b><br>",
                          mod6$call[[2]], "<br>",
                          mod6_p[[1]],"<br>",
                          mod6_p[[2]]),
           label.color = NA,
           fill = NA,
           x = 19,
           y = 4.6,
           size = 3,
           hjust = 1,
           vjust = 1) +
  theme(legend.position = c(.98,.75),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "All data, unsummarized") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(),
        plot.title.position = "plot")

p6 + ggview::canvas(4,4)




# combine -----------------------------------------------------------------
p1_sens <- p1 + 
  labs(title = "HS, all replicates") +
  theme(panel.border = element_rect(linewidth = 2))

library(patchwork)
cti_sens <- (
  (( p1_sens  | p4) + plot_layout(axis_titles = "collect")) /
    (( p2  | p5) + plot_layout(axis_titles = "collect")) /
    (( p3  | p6) + plot_layout(axis_titles = "collect")) 
)

cti_sens + ggview::canvas(8,12)

ggsave(cti_sens,
       filename = "outputs/cti/bylevel_cti_sens_test.png",
       height = 12,
       width = 8)


# -------------------------------------------------------------------------
# Scratch / old ----------------------------------------------------------------
# -------------------------------------------------------------------------



# p4 HS, transect effect --------------------------------------------------
df4 <- pa_therm_filtered %>%
  # find CTI per level/year
  group_by(level, tidalheight, year, transect) %>%
  summarize(cti = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are 4+ species
  filter(n_spp >= 4) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level))

df4 %>% ggplot(aes(x = year, y = cti, group = transect)) +
  geom_point() +
  geom_smooth(method = "lm", se=F) +
  facet_wrap(~tidalheight)

mod2 <- test_mods(df2)
summary(mod2)
mod2_p <- get_mod_p(mod2)
au2 <- broom::augment(mod2, df2, interval = "prediction")
p2 <- au2 %>% 
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.2,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
  annotate(geom = "richtext",
           label = paste0("<b>Best Model:</b><br>",
                          mod2$call[[2]], "<br>",
                          mod2_p[[1]],"<br>",
                          mod2_p[[2]]),
           label.color = NA,
           fill = NA,
           x = 10.2,
           y = 2.9,
           size = 3,
           hjust = 0,
           vjust = 1) +
  coord_cartesian(xlim = c(10.05,12.15),
                  ylim = c(-.4,3),
                  expand = F) +
  
  theme(legend.position = c(.08,.1),
        legend.justification = c(0,0),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "HS, First replicates only") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(),
        plot.title.position = "plot")

p2 + ggview::canvas(4,4)





# find CTI per level - TOTAL ----------------------------------------------
# first, find the community thermal index (mean of STIs) for all 
# unique species in each level/year. This is a bit weird as it 
# compares different numbers of species, but as long as we assume
# that thermal index is not correlated with rarity / detectability
# of species, it should be fine. 
cti_all <- pa %>%
  # group by tidal height per year
  group_by(level, tidalheight, year) %>% 
  # get rid of all sampling variables except the 3 above
  select(-replicate,
         -data_taken,
         -pres_cover,
         -pres,
         -both_taken,
         -transect,
         -pres_count) %>% 
  # find distinct species
  #distinct() %>%
  summarize(cti_mean = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are more than 4 species
  filter(n_spp >= 3) %>%
  # keep only levels that have data 3 or more times so that
  # a trend can be formed 
  group_by(level) %>% 
  add_count(name = "n_times_level_represented") %>% 
  filter(n_times_level_represented > 2) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level))
  

cti_all %>%
  ggplot(aes(x = year, y = cti_mean,
             #color = level, 
             group = tidalheight)) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight_f) +
  geom_point(alpha = .6) +
  gghighlight::gghighlight(
    unhighlighted_params = list(color = "grey80",
                                alpha = .3),
    use_direct_label = F) 


# make a few models -------------------------------------------------------
cti_all_add <- lm(data = cti_all,
                  formula = "cti_mean ~ year + tidalheight")
summary(cti_all_add) 
AIC(cti_all_add) #-206.8323

cti_all_mult <- lm(data = cti_all,
                   formula = "cti_mean ~ year * tidalheight")

summary(cti_all_mult)
AIC(cti_all_mult) #-224.7636

cti_all_year <- lm(data = cti_all,
                   formula = "cti_mean ~ year")

summary(cti_all_year)
AIC(cti_all_year) # -204.6188

#cti_all_year_wt <- lm(data = cti_all,
#                   formula = "cti_mean ~ year",
#                   weights = n_spp)
#
#summary(cti_all_year_wt)
#AIC(cti_all_year_wt) #-206.3783

#cti_all_add_wt <- lm(data = cti_all,
#                      formula = "cti_mean ~ year + tidalheight",
#                      weights = n_spp)
#summary(cti_all_year_wt)
#AIC(cti_all_year_wt) #-248.7909
 
#cti_all_mult_wt <- lm(data = cti_all,
#                     formula = "cti_mean ~ year * tidalheight",
#                     weights = n_spp)
#summary(cti_all_mult_wt)
#AIC(cti_all_mult_wt) #-265.6291

cti_all_par <- lm(data = cti_all,
                  formula = "cti_mean ~ year + I(tidalheight^2)")
summary(cti_all_par)
AIC(cti_all_par) #-203.6201

cti_all_mult_par <- lm(data = cti_all,
                  formula = "cti_mean ~ year * I(tidalheight^2)")
summary(cti_all_mult_par)
AIC(cti_all_mult_par) #-214.0151

performance::compare_performance(
  cti_all_add,
  #cti_all_add_wt,
  cti_all_mult,
  #cti_all_mult_wt,
  cti_all_par,
  cti_all_year,
  cti_all_mult_par,
  #cti_all_year_wt,
  rank = T
)
# cti_all_mult_wt is the best
plot(cti_all_mult)
rm(cti_all_add, cti_all_par, cti_all_year, cti_all_mult_par)

pred_df <- data.frame(
  tidalheight = rep(unique(cti_all$tidalheight), times = length(unique(cti_all$year))),
  year = rep(unique(cti_all$year), each = length(unique(cti_all$tidalheight)))
)

# see if I need to add the weight variable (n_spp) in the "newdf" for prediction
au <- broom::augment(cti_all_mult,
                     newdata = pred_df,
                     interval = "prediction"
                     )

p <- au %>% 
  left_join(cti_all) %>%
  ggplot(aes(x = year, 
             y = .fitted,
             group = tidalheight,
             color = tidalheight)) +
  
  # add true appledore temp trend - save both ways
  geom_ribbon(data = temptrend,
              inherit.aes =F,
              aes(x = year, 
                  ymin = .lower,
                  ymax = .upper),
              fill = "grey70",
              alpha = .5) +
  geom_line(data = temptrend,
            inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
 # geom_ribbon(aes(ymin = .lower,
 #                 ymax = .upper, fill = tidalheight),
 #             alpha = .1) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(y = cti_mean,
                 fill = tidalheight),
             alpha = .7,
             color = "black",
             position = position_jitter(width = .3,
                                        height = 0),
             shape = 21,
             size = 2,
             stroke = .2) +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  scale_fill_gradient(low = "grey40",
                       high = "cyan3") +
  labs(x = NULL,
       y = "Community Thermal Index\n(Averaged Thermal Affinity of All Species)",
       color = "Tidal Height",
       fill = "Tidal Height") +
  theme(legend.position = "inside",
        legend.position.inside = c(.02, .98),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  
  guides(fill  = guide_colorbar(theme = theme(legend.title.position = "top")),
         color = guide_colorbar(theme = theme(legend.title.position = "top"))) +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  coord_cartesian(ylim  = c(10.5, 12))

p + ggview::canvas(4.5,4)

ggsave(p,
       filename = "outputs/cti/cti_p1_with_temptrend.png",
       width = 4.5, height = 4, unit = "in")

# plot across tidal height transect
p2 <- au %>% 
  left_join(cti_all) %>%
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  geom_point(data = cti_all,
             aes(x = cti_mean,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  theme(legend.position = c(.02,1),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "CTI Across Shore Levels") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank())


p2 + ggview::canvas(4,4)
ggsave(p2,
       filename = "outputs/cti/cti_p2.png",
       width = 4.5, height = 4, unit = "in")

p2_with_temp <- p2 +
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -1.1, yend = -1.1,
           arrow = arrow(#ends = "both",
                         length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -1,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) 
p2_with_temp +  ggview::canvas(4,4)

ggsave(p2_with_temp,
       filename = "outputs/cti/cti_p2_with_temptrend.png",
       width = 4.5, height = 4, unit = "in")



# weighted CTI per tidal height -------------------------------------------
# since CTI is normally abundance-weighted, let's do the same thing again,
# (e.g., average STIs per tidal height), but this time we'll maintain the 
# number of transects in which each species was found for a proxy for abundance
# e.g., still an average CTI per tidalheight, but if one species was found 
# in one transect, and another in 5, species b will count towards the average
# 5 times more. 
cti_abundance_weighted <- pa %>%
  
  # first, cut to replicate 1 only
  #filter(replicate == 1) %>%
  
  # then, find distinct species within transects within years within tidal heights
  group_by(transect, year, tidalheight, replicate) %>%
  distinct(organism, mean_monthly_mean) %>%
  
  # now group by tidal height and year, but this will serve as a
  # pseudo-abundance-weighted CTI where some species count more than others
  group_by(tidalheight, year) %>%
  summarize(cti = mean(mean_monthly_mean),
            n_unique_spp = n_distinct(organism),
            n_transects = n_distinct(transect),
            n_obs = n()) %>% ungroup() %>%
  # filter out where the total number of speices is < 4
  filter(n_unique_spp >= 3) 

cti_abundance_weighted %>%
  ggplot(aes(x  = year, 
             y = cti,
             group = tidalheight)) +
  geom_point() +
  geom_smooth(method = "lm")

# make and test models
cti_ab_add <- lm(data = cti_abundance_weighted,
                  formula = "cti ~ year + tidalheight")
summary(cti_ab_add) 
AIC(cti_ab_add) #-145.319

cti_ab_mult <- lm(data = cti_abundance_weighted,
                   formula = "cti ~ year * tidalheight")

summary(cti_ab_mult)
AIC(cti_ab_mult) #-175.8601

cti_ab_year <- lm(data = cti_abundance_weighted,
                   formula = "cti ~ year")

summary(cti_ab_year)
AIC(cti_ab_year) #-113.0945


cti_ab_par <- lm(data = cti_abundance_weighted,
                  formula = "cti ~ year + I(tidalheight^2)")
summary(cti_ab_par)
AIC(cti_ab_par) #-119.9018

AIC(cti_ab_add, cti_ab_mult, cti_ab_par, cti_ab_year) %>% arrange((AIC))
performance::compare_performance(cti_ab_add, cti_ab_mult, cti_ab_par, cti_ab_year, rank = T)

# augment and plot

au_ab <- broom::augment(cti_ab_mult,
                     newdata = pred_df,
                     interval = "prediction"
)

p_ab <- au_ab %>% 
  left_join(cti_abundance_weighted) %>%
  ggplot(aes(x = year, 
             y = .fitted,
             group = tidalheight,
             color = tidalheight)) +
  
  # add true appledore temp trend - save both ways
  geom_ribbon(data = temptrend,
              inherit.aes =F,
              aes(x = year, 
                  ymin = .lower,
                  ymax = .upper),
              fill = "grey70",
              alpha = .5) +
  geom_line(data = temptrend,
            inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  # geom_ribbon(aes(ymin = .lower,
  #                 ymax = .upper, fill = tidalheight),
  #             alpha = .1) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(y = cti,
                 fill = tidalheight),
             alpha = .7,
             color = "black",
             position = position_jitter(width = .3,
                                        height = 0),
             shape = 21,
             size = 2,
             stroke = .2) +
  scale_color_gradient(low = "grey40",
                       high = "cyan3") +
  scale_fill_gradient(low = "grey40",
                      high = "cyan3") +
  labs(x = NULL,
       y = "Community Thermal Index\n(Averaged Thermal Affinity of All Species",
       color = "Tidal Height",
       fill = "Tidal Height") +
  theme(legend.position = "inside",
        legend.position.inside = c(.02, .98),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(30,"pt")) +
  
  guides(fill  = guide_colorbar(theme = theme(legend.title.position = "top")),
         color = guide_colorbar(theme = theme(legend.title.position = "top"))) +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  coord_cartesian(ylim = c(10.5, 12.2))

p_ab + ggview::canvas(4.5,4)
ggsave(p_ab,
       filename = "outputs/cti/cti_p1_abwt_with_temptrend.png",
       width = 4.5, height = 4, unit = "in")

cti_abundance_weighted %>% glimpse()
cti_all %>% glimpse()

# plot across tidal height transect
p2_ab <- au_ab %>% 
  left_join(cti_abundance_weighted) %>%
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  geom_point(#data = cti_all,
             aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  theme(legend.position = c(.98,1),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "CTI Across Shore Levels") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank())

p2_ab + ggview::canvas(4,4)
ggsave(p2_ab,
       filename = "outputs/cti/cti_p2_abwt.png",
       width = 4.5, height = 4, unit = "in")

p2_ab_with_temp <- p2_ab +
  geom_vline(xintercept = first_temp) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash") +
  theme(legend.position.inside = c(.05,1),
        legend.justification = c(0,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -1.1, yend = -1.1,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -1,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) 
p2_ab_with_temp +  ggview::canvas(4,4)
ggsave(p2_ab_with_temp,
       filename = "outputs/cti/cti_p2_abwt_with_temptrend.png",
       width = 4.5, height = 4, unit = "in")


library(patchwork)
p2_both <- 
  p2 + 
  labs(title = "CTI Across Shore Levels Unweighted",
            subtitle = "Mean of STIs for unique species per shore level") +
  theme(plot.title.position = "plot",
        plot.subtitle = element_text(size = 10)) +
  p2_ab +
  labs(title = "CTI Across Shore Levels Weighted",
       subtitle = "Mean of STIs per shore level, weighted by presences") +
  theme(plot.title.position = "plot",
        plot.subtitle = element_text(size = 10)) 

p2_both + ggview::canvas(8,4)
ggsave(p2_both, 
       filename = "outputs/cti/cti_p2_compare_mods.png",
       width = 8,
       height = 4
)

p2_both_with_temp <- 
  p2_with_temp + 
  labs(title = "CTI Across Shore Levels Unweighted",
       subtitle = "Mean of STIs for unique species per shore level") +
  theme(plot.title.position = "plot",
        plot.subtitle = element_text(size = 10)) +
  p2_ab_with_temp +
  labs(title = "CTI Across Shore Levels Weighted",
       subtitle = "Mean of STIs per shore level, weighted by presences") +
  theme(plot.title.position = "plot",
        plot.subtitle = element_text(size = 10)) 

p2_both_with_temp + ggview::canvas(8,4)
ggsave(p2_both_with_temp, 
       filename = "outputs/cti/cti_p2_compare_mods_with_temp.png",
       width = 8,
       height = 4
)

# cti per replicate per depth ---------------------------------------------
# here, we'll average the CTI in every individual replicate 
cti_per_replicate <- 
  pa %>%
  group_by(year, transect, tidalheight, replicate) %>%
  summarize(cti_mean = mean(mean_monthly_mean, na.rm=T),
            spp_in_replicate = n()) %>%
  filter(spp_in_replicate > 4) %>%
  group_by(tidalheight) %>%
  add_count(name = "n_in_tidalheight") %>%
  mutate(n_years = length(unique(year))) %>%
  filter(n_years > 3)

cti_per_replicate %>%
  ggplot(aes(x = year, 
             y = cti_mean,
             group = tidalheight)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight) +
  gghighlight::gghighlight(
    use_direct_label = F,
    unhighlighted_params = list(color = "grey80",
                                alpha = .3)
  )

range(cti_per_replicate$year)
cti_per_replicate <- cti_per_replicate %>% mutate(year = year - 1981)
cti_rep_add <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year + tidalheight + (1|transect)")
summary(cti_rep_add)

cti_rep_add_wt <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year + tidalheight + (1|transect)",
                  weights = spp_in_replicate)
summary(cti_rep_add_wt)


cti_rep_mult <- lmer(data = cti_per_replicate,
                     formula = "cti_mean ~ year * tidalheight + (1|transect)")
summary(cti_rep_mult)
marginaleffects::avg_slopes(cti_rep_mult)

cti_rep_mult_wt <- lmer(data = cti_per_replicate,
                   formula = "cti_mean ~ year * tidalheight + (1|transect)",
                   weights = spp_in_replicate)
summary(cti_rep_mult_wt)
marginaleffects::avg_slopes(cti_rep_mult_wt)

cti_rep_year <- lmer(data = cti_per_replicate,
                   formula = "cti_mean ~ year + (1|transect)")
summary(cti_rep_year)


cti_rep_year_wt <- lmer(data = cti_per_replicate,
                   formula = "cti_mean ~ year + (1|transect)",
                   weights = spp_in_replicate)
summary(cti_rep_year_wt)


cti_rep_year_par <- lmer(data = cti_per_replicate,
                      formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)")
summary(cti_rep_year_par)
marginaleffects::avg_slopes(cti_rep_year_par)

cti_rep_par_wt <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)",
                  weights = spp_in_replicate)
summary(cti_rep_par_wt)
marginaleffects::avg_slopes(cti_rep_par_wt)



cti_rep_multpar <- lmer(data = cti_per_replicate,
                  formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)")
summary(cti_rep_multpar)
marginaleffects::avg_slopes(cti_rep_multpar)

cti_rep_multpar_wt <- lmer(data = cti_per_replicate,
                     formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)",
                     weights = spp_in_replicate)
summary(cti_rep_multpar_wt)
marginaleffects::avg_slopes(cti_rep_multpar_wt)



performance::compare_performance(
  cti_rep_add,
  cti_rep_add_wt,
  cti_rep_mult,
  cti_rep_mult_wt, 
  cti_rep_year, 
  cti_rep_year_wt,
  #cti_rep_par,
  cti_rep_par_wt,
  cti_rep_multpar,
  cti_rep_multpar_wt,
  rank = T
)


au_rep <- broom.mixed::augment(cti_rep_mult,
                         cti_per_replicate)
marginaleffects::avg_slopes(cti_rep_mult)

au_rep %>%
  ggplot(aes(x = year, 
             y = .fitted,
             group = interaction(tidalheight, transect),
             color = tidalheight)) +
  geom_point(aes(fill = tidalheight,
                 y = cti_mean),
             size = 2,
             alpha = .3,
             position = position_jitter(width = .3,
                                        height = 0)) +
  geom_line(linewidth = 1.5) +
  scale_fill_viridis_c() +
  scale_color_viridis_c()
  

# NOTE: Jarrett says try modelbased
# https://easystats.github.io/modelbased/
#install.packages("modelbased")
library(modelbased)

estimate_relation(cti_rep_mult) |> plot()


library(ggeffects)

by_rep_year <- ggpredict(cti_rep_mult, "year",
                         interval = "prediction")

by_rep_year %>%
  ggplot(aes(x = x,
             y = predicted)) +
  geom_point(data = cti_per_replicate,
             aes(x = year, 
                 y = cti_mean,
                 fill = tidalheight),
             shape = 21,
             size = 2,
             alpha = .5,
             position = position_jitter(width = .3, 
                                        height = 0)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = .5) +
  scale_x_continuous(breaks = c(-1,9,19,29,39),
                     labels = function(x) x + 1981) +
  geom_line(linewidth = 1.5)  +
  scale_fill_viridis_c() +
  labs(x = NULL,
       y = "Community Thermal Index",
       fill = "Tidal Height") +
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(theme = theme(legend.title.position = "top")))

predict_response(cti_rep_mult,
        terms = "year")

# find per transect -------------------------------------------------------


cti_per_transect <- 
  pa %>%
  group_by(year, transect, tidalheight) %>%
  distinct(organism, mean_monthly_mean) %>%
  summarize(cti_mean = mean(mean_monthly_mean, na.rm=T),
            spp_in_transect = n()) %>%
  filter(spp_in_transect > 4) %>%
  group_by(tidalheight) %>%
  add_count(name = "n_in_tidalheight") %>%
  mutate(n_years = length(unique(year))) %>%
  filter(n_years > 3)

cti_per_transect %>%
  ggplot(aes(x = year, 
             y = cti_mean,
             group = tidalheight)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight) +
  gghighlight::gghighlight(
    use_direct_label = F,
    unhighlighted_params = list(color = "grey80",
                                alpha = .3)
  )



cti_per_transect <- cti_per_transect %>% mutate(year = year - 1981)
cti_tran_add <- lmer(data = cti_per_transect,
                  formula = "cti_mean ~ year + tidalheight + (1|transect)")
summary(cti_tran_add)

cti_tran_add_wt <- lmer(data = cti_per_transect,
                     formula = "cti_mean ~ year + tidalheight + (1|transect)",
                     weights = spp_in_transect)
summary(cti_tran_add_wt)


cti_tran_mult <- lmer(data = cti_per_transect,
                   formula = "cti_mean ~ year * tidalheight + (1|transect)")
summary(cti_tran_mult)
marginaleffects::avg_slopes(cti_tran_mult)

cti_tran_mult_wt <- lmer(data = cti_per_transect,
                      formula = "cti_mean ~ year * tidalheight + (1|transect) ",
                      weights = spp_in_transect)
summary(cti_tran_mult_wt)
marginaleffects::avg_slopes(cti_tran_mult_wt)

cti_tran_year <- lmer(data = cti_per_transect,
                   formula = "cti_mean ~ year + (1|transect)")
summary(cti_tran_year)


cti_tran_year_wt <- lmer(data = cti_per_transect,
                      formula = "cti_mean ~ year + (1|transect)",
                      weights = spp_in_transect)
summary(cti_tran_year_wt)


cti_tran_par <- lmer(data = cti_per_transect,
                  formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)")
summary(cti_tran_par)
marginaleffects::avg_slopes(cti_tran_par)

cti_tran_par_wt <- lmer(data = cti_per_transect,
                     formula = "cti_mean ~ year + I(tidalheight^2) + (1|transect)",
                     weights = spp_in_transect)
summary(cti_tran_par_wt)
marginaleffects::avg_slopes(cti_tran_par_wt)


cti_tran_multpar <- lmer(data = cti_per_transect,
                      formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)")
summary(cti_tran_multpar)
marginaleffects::avg_slopes(cti_tran_multpar)

cti_tran_multpar_wt <- lmer(data = cti_per_transect,
                         formula = "cti_mean ~ year * I(tidalheight^2) + (1|transect)",
                         weights = spp_in_transect)
summary(cti_tran_multpar_wt)
marginaleffects::avg_slopes(cti_tran_multpar_wt)

performance::compare_performance(
  cti_tran_add,
  cti_tran_add_wt,
  cti_tran_mult,
  cti_tran_mult_wt, 
  cti_tran_year, 
  cti_tran_year_wt,
  cti_tran_par,
  cti_tran_par_wt,
  cti_tran_multpar,
  cti_tran_multpar_wt,
  rank = T
)


by_transect <- ggpredict(cti_tran_mult, "year",
                         interval = "prediction")

by_transect %>%
  ggplot(aes(x = x,
             y = predicted)) +
  geom_point(data = cti_per_transect,
             aes(x = year, 
                 y = cti_mean,
                 fill = tidalheight),
             shape = 21,
             size = 2,
             alpha = .5,
             position = position_jitter(width = .3, 
                                        height = 0)) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = .5) +
  scale_x_continuous(breaks = c(-1,9,19,29,39),
                     labels = function(x) x + 1981) +
  geom_line(linewidth = 1.5)  +
  scale_fill_viridis_c() +
  labs(x = NULL,
       y = "Community Thermal Index",
       fill = "Tidal Height") +
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.key.width = unit(dev.size()[1] / 10, "inches"),
        legend.key.height = unit(dev.size()[1] / 20, "inches"),
        legend.background = element_blank()) +
  guides(fill = guide_colorbar(theme = theme(legend.title.position = "top")))







# -------------------------------------------------------------------------


# scratch -----------------------------------------------------------------


# -------------------------------------------------------------------------


p2 <- sjPlot::plot_model(cti_all_mult_wt, type = "pred", 
                         terms = "year", show.data = T,
                         jitter = .02,
                         color = "cyan3", 
                         shape=16,
                         dot.size = 1.2,
                         dot.shape = 16,
                         line.size = 1.75,
                         dot.alpha = .4,
                         se = T,
                         alpha = .4) +
  labs(x = NULL,
       y = "Community Thermal Index\n(Averaged Thermal Affinity of All Species",
       title = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(labels = ~paste0(.,"°"))

p2 + ggview::canvas(4,4)

ggsave(p2,
       filename = "outputs/cti/cti_p2.png",
       width = 4.5, height = 4, unit = "in")

marginaleffects::avg_slopes(cti_all_mult_wt, variables = "year")
marginaleffects::avg_slopes(cti_all_mult_wt, variables = c("year","tidalheight"))
marginaleffects::slopes(cti_all_mult_wt) %>%
  as_tibble() %>%
  group_by(year,tidalheight) %>% distinct(estimate)
group_by(estimate) %>% distinct(year,tidalheight)
filter(estimate == max(estimate)) %>% glimpse()

