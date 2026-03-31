# model total CTI per transect


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(lme4)
library(effects)
library(ggtext)
theme_set(ggthemes::theme_few())



# upload data -------------------------------------------------------------
pa_therm <- readr::read_csv("data-processed/cti-data/cti-data-all.csv")
pa_therm_filtered <- readr::read_csv("data-processed/cti-data/cti-data-highly-sampled.csv")
temptrend <- readr::read_csv("data-processed/cti-data/temp-change-augment.csv")
tempmod <- readRDS("data-processed/cti-data/temp-change-model.rds")


# calculate CTI -----------------------------------------------------------
# first, we'll do the most logical thing:
# CTI as mean STIs of all species on the island per year. 
# Here, we retain the most information, and summarize per year 
# so that all years count the same in the regression
df1 <- pa_therm_filtered %>%
  group_by(year) %>% 
  summarise(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean,.25),
            cti_75 = quantile(mean_monthly_mean,.75))
mod1 <- lm(data = df1,
           "cti ~ year")
summary(mod1)
au1 <- broom::augment(mod1, df1, interval = "prediction")

p1 <- au1 %>%
  ggplot(aes(x = year, y = cti)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_line(aes(y = .fitted))+
  
  # now add points and 25 - 75 quantiles
  geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15) +
  geom_point(aes(y = cti), 
             shape = 16, size = 1, alpha = .7) +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod1$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod1)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                              "<br>",
                              "<b>Temp. Trend:</b> ",
                              round(tempmod$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(tempmod)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                              ),
           hjust = 1,vjust = .5,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index") +
  scale_y_continuous(labels = ~paste0(.,"°"))

p1 + ggview::canvas(4,4)

ggsave(p1,
       file = "outputs/cti/total_cti_p1.png",
       width = 4,
       height = 4)
# save as rds so we can access it to combine with other plot
saveRDS(p1, file = "outputs/cti/total_cti_p1.rds")



# calculate CTI first rep -----------------------------------------------
df2 <- pa_therm_filtered %>%
  filter(replicate == 1) %>%
  group_by(year) %>% 
  summarise(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean,.25),
            cti_75 = quantile(mean_monthly_mean,.75))
mod2 <- lm(data = df2,
           "cti ~ year")
summary(mod2)
au2 <- broom::augment(mod2, df2, interval = "prediction")

p2 <- au2 %>%
  ggplot(aes(x = year, y = cti)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_line(aes(y = .fitted))+
  
  # now add points and 25 - 75 quantiles
  geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15) +
  geom_point(aes(y = cti), 
             shape = 16, size = 1, alpha = .7) +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod2$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod2)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                              ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "HS, first replicates averaged per year") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(plot.title.position = "plot")

p2 + ggview::canvas(4,4)


# p3 all species in filtered quadrats -------------------------------------
mod3 <- lm(pa_therm_filtered, 
           formula = "mean_monthly_mean ~ year")
summary(mod3)
au3 <- broom::augment(mod3, pa_therm_filtered, interval = "prediction")

p3 <- au3 %>%
  ggplot(aes(x = year, y = mean_monthly_mean)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey40", alpha = .8) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_point(aes(y = mean_monthly_mean),
             size = .25, alpha = .25, shape = 16,
             position = position_jitter(width = .3,height = .05)) +
  geom_line(aes(y = .fitted))+
  
  # add annotations
  annotate(geom = "richtext",
           x = 2024, y = 8,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod3$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod3)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                              ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "HS, all species unsummarized") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(plot.title.position = "plot")  +
  coord_cartesian(xlim = c())

p3 + ggview::canvas(4,4)



# p4: all species summarized per year -------------------------------------
df4 <- pa_therm %>%
  group_by(year) %>% 
  summarise(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean,.25),
            cti_75 = quantile(mean_monthly_mean,.75))
mod4 <- lm(data = df4,
           "cti ~ year")
summary(mod4)
au4 <- broom::augment(mod4, df4, interval = "prediction")

p4 <- au4 %>%
  ggplot(aes(x = year, y = cti)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_line(aes(y = .fitted))+
  
  # now add points and 25 - 75 quantiles
  geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15) +
  geom_point(aes(y = cti), 
             shape = 16, size = 1, alpha = .7) +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod4$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod4)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"}
                              else(if(summary(mod4)$coefficients["year","Pr(>|t|)"] < 0.05){"< 0.05"}),
                              ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "All, summarized per year") +
  theme(plot.title.position = "plot") +
  scale_y_continuous(labels = ~paste0(.,"°"))

p4 + ggview::canvas(4,4)



# p5. all data unsummarized -----------------------------------------------
mod5 <- lm(pa_therm, 
           formula = "mean_monthly_mean ~ year")
summary(mod5)
au5 <- broom::augment(mod5, pa_therm, interval = "prediction")

p5 <- au5 %>%
  ggplot(aes(x = year, y = mean_monthly_mean)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey40", alpha = .8) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_point(aes(y = mean_monthly_mean),
             size = .25, alpha = .25, shape = 16,
             position = position_jitter(width = .3,height = .05)) +
  geom_line(aes(y = .fitted))+
  
  # add annotations
  annotate(geom = "richtext",
           x = 2024, y = 7,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod5$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod5)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                               ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "All, all species unsummarized") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(plot.title.position = "plot")  +
  coord_cartesian(xlim = c())

p5 + ggview::canvas(4,4)



# p6. filtered transects, transect effect ---------------------------------
df6 <- pa_therm_filtered %>%
  group_by(year, transect_label) %>%
  summarize(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean, .25),
            cti_75 = quantile(mean_monthly_mean, .75))
mod6 <- lm(data = df6,
           "cti ~ year + transect_label")
summary(mod6)
au6 <- broom::augment(mod6, df6, interval = "prediction")

p6 <- au6 %>%
  ggplot(aes(x = year, y = cti, group = transect_label)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .3) +
  
  # now add points and 25 - 75 quantiles
 # geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
 #              linewidth = .15) +
 # geom_point(aes(y = cti), 
 #            shape = 16, size = 1, alpha = .7) +
  
  geom_pointrange(aes(y = cti, ymin = cti_25, ymax = cti_75),
                  position = position_dodge(width = 1),
                  linewidth = .15,
                  shape = 16, size = .25,
                  alpha = .5) +
  
  geom_line(aes(y = .fitted))  +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod6$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod6)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"}
           ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "HS, summarized per year, transect effect") +
  theme(plot.title.position = "plot") +
  scale_y_continuous(labels = ~paste0(.,"°"))


p6 + ggview::canvas(4,4)



# p7. all data, transect effect -------------------------------------------
df7 <- pa_therm %>%
  group_by(year, transect_label) %>%
  summarize(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean, .25),
            cti_75 = quantile(mean_monthly_mean, .75))
mod7 <- lm(data = df7,
           "cti ~ year + transect_label")
summary(mod7)
au7 <- broom::augment(mod7, df7, interval = "prediction")

p7 <- au7 %>%
  ggplot(aes(x = year, y = cti, group = transect_label)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .3) +
  
  # now add points and 25 - 75 quantiles
  # geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
  #              linewidth = .15) +
  # geom_point(aes(y = cti), 
  #            shape = 16, size = 1, alpha = .7) +
  
  geom_pointrange(aes(y = cti, ymin = cti_25, ymax = cti_75),
                  position = position_dodge(width = 1),
                  linewidth = .15,
                  shape = 16, size = .25,
                  alpha = .5) +
  
  geom_line(aes(y = .fitted))  +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod7$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod7)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                              ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "All, summarized per year, transect effect") +
  theme(plot.title.position = "plot") +
  scale_y_continuous(labels = ~paste0(.,"°"))


p7 + ggview::canvas(4,4)



# p8 unique spp -----------------------------------------------------------
df8 <- pa_therm_filtered %>% 
  group_by(year) %>%
  distinct(organism, mean_monthly_mean) %>%
  summarize(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean, .25),
            cti_75 = quantile(mean_monthly_mean, .75))
mod8 <- lm(data = df8,
           "cti ~ year")
summary(mod8)
au8 <- broom::augment(mod8, df8, interval = "prediction")


p8 <- au8 %>%
  ggplot(aes(x = year, y = cti)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_line(aes(y = .fitted))+
  
  # now add points and 25 - 75 quantiles
  geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15) +
  geom_point(aes(y = cti), 
             shape = 16, size = 1, alpha = .7) +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod8$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod8)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                              ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "HS, unique species per year") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(plot.title.position = "plot")

p8 + ggview::canvas(4,4)


# p9 first replicates, all data -------------------------------------------
df9 <- pa_therm %>%
  filter(replicate == 1) %>%
  group_by(year) %>% 
  summarise(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean,.25),
            cti_75 = quantile(mean_monthly_mean,.75))
mod9 <- lm(data = df9,
           "cti ~ year")
summary(mod9)
au9 <- broom::augment(mod9, df9, interval = "prediction")

p9 <- au9 %>%
  ggplot(aes(x = year, y = cti)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_line(aes(y = .fitted))+
  
  # now add points and 25 - 75 quantiles
  geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15) +
  geom_point(aes(y = cti), 
             shape = 16, size = 1, alpha = .7) +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod9$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod9)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"}
                              else if(summary(mod9)$coefficients["year","Pr(>|t|)"] < 0.05){"< 0.05"}
                              else {"= 0.12"},
                              ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "All, first replicates averaged per year") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(plot.title.position = "plot")

p9 + ggview::canvas(4,4)



# p10, all transects unique spp -------------------------------------------

df10 <- pa_therm %>% 
  group_by(year) %>%
  distinct(organism, mean_monthly_mean) %>%
  summarize(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean, .25),
            cti_75 = quantile(mean_monthly_mean, .75))
mod10 <- lm(data = df10,
           "cti ~ year")
summary(mod10)
au10 <- broom::augment(mod10, df10, interval = "prediction")


p10 <- au10 %>%
  ggplot(aes(x = year, y = cti)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              fill = "darkturquoise",
              alpha = .5) +
  geom_line(aes(y = .fitted))+
  
  # now add points and 25 - 75 quantiles
  geom_segment(aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15) +
  geom_point(aes(y = cti), 
             shape = 16, size = 1, alpha = .7) +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod10$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod10)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
           ),
           hjust = 1,vjust = 1,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       title = "All, unique species per year") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(plot.title.position = "plot")

p10 + ggview::canvas(4,4)



# combine plots -----------------------------------------------------------

library(patchwork)

p1_sens <- p1 + 
  labs(title = "HS, summarized per year") +
  theme(plot.title.position = "plot",
        panel.border = element_rect(linewidth = 2))

cti_sens <- ( p1_sens | p2 | p3 | p4 | p5) /
  ( p6 | p7 | p8 | p9 | p10)

cti_sens + ggview::canvas(16,8)

# rearrange in a way that makes sense:
# p1  p4
# p3  p5
# p6  p7
# p2  p9
# p8  p10

cti_sens_2 <- (
  (( p1_sens  | p4) + plot_layout(axis_titles = "collect")) /
    (( p3  | p5) + plot_layout(axis_titles = "collect")) /
    (( p6  | p7) + plot_layout(axis_titles = "collect")) /
    (( p2  | p9) + plot_layout(axis_titles = "collect")) /
    (( p8  | p10) + plot_layout(axis_titles = "collect")) 
)

cti_sens_2 + ggview::canvas(8,14)

ggsave(cti_sens_2,
       filename = "outputs/cti/total_cti_p2_sens_test.png",
       height = 14,
       width = 8)


## scratch -----------------------------------------------------------------
## try to plot first and last appearances
#whole_island_ind %>%
#  ggplot(aes(x = year,
#             y = mean_monthly_mean)) +
#  geom_point(#position = position_jitter(width = .4),
#             alpha = .2, shape = 16) +
#  geom_point(data = . %>%
#               #filter(year < (max(year) -5)) %>%
#               group_by(organism, mean_monthly_mean) %>%
#               summarize(last_year = max(year),
#                         n_years = n()) %>%
#               filter(n_years > 1),
#             aes(x = last_year),
#             color = "red", size = 2
#             ) +
#  geom_point(data = . %>%
#               #filter(year < (max(year) -5)) %>%
#               group_by(organism, mean_monthly_mean) %>%
#               summarize(first_year = min(year),
#                         n_years = n()) %>%
#               filter(n_years > 1),
#             aes(x = first_year),
#             color = "blue", size = 2
#  )
#
#
#whole_island_ind %>%
#  group_by(organism, mean_monthly_mean) %>%
#  summarize(first_year = min(year),
#            last_year = max(year),
#            n_years = n()) %>%
#  filter(n_years > 1) %>%
#  ungroup() %>%
#  mutate(first_year_exclusive = case_when(first_year < 1982 + 2 ~ NA,
#                                          TRUE ~ first_year),
#         last_year_exclusive = case_when(last_year > 2023 - 2 ~ NA,
#                                         TRUE ~ last_year)) %>%
#  ggplot(aes(x = first_year,
#             xend = last_year,
#             y = mean_monthly_mean,
#             yend = mean_monthly_mean)) +
#  geom_segment(linewidth = .1) +
#  geom_point(aes(x = first_year_exclusive),
#             color = "black",
#             fill = "cyan3",
#             stroke = .3,
#             shape = 21) +
#  geom_smooth(aes(x = first_year_exclusive),
#              method = "lm",
#              fill = "cyan3",
#              color = "black") +
#  geom_point(aes(x = last_year_exclusive),
#             color = "black",
#             fill = "tomato2",
#             stroke = .3,
#             shape = 21) +
#  geom_smooth(aes(x = last_year_exclusive),
#              method = "lm",
#              fill = "tomato2",
#              color = "black") +
#  labs(x = NULL,
#       y = "Species Thermal Affinity",
#       title = "First and Last Appearances by STI")
#
#
#whole_island_ind %>%
#  ungroup() %>%
#  filter(year < (max(year) -5)) %>%
#  group_by(organism, mean_monthly_mean) %>%
#  summarize(last_year = max(year))
#  #