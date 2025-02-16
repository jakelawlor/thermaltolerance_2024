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
traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed.csv")



pa_therm <- pa_therm %>%
  left_join(traits %>% rename(organism = gen_spp) %>% 
              mutate(group = motility_adult))
pa_therm_filtered <- pa_therm_filtered %>%
  left_join(traits %>% rename(organism = gen_spp)%>% 
              mutate(group = motility_adult))


# calculate CTI by group-----------------------------------------------------------
# first, we'll do the most logical thing:
# CTI as mean STIs of all species on the island per year. 
# Here, we retain the most information, and summarize per year 
# so that all years count the same in the regression
df1 <- pa_therm_filtered %>%
  group_by(year,group) %>% 
  summarise(cti = mean(mean_monthly_mean),
            cti_25 = quantile(mean_monthly_mean,.25),
            cti_75 = quantile(mean_monthly_mean,.75))

df1 %>%
  ggplot(aes(x = year,
             y = cti,
             color = group)) +
  geom_segment(aes(y = cti_25, 
                   yend = cti_75,
                   xend = year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~group)

mod1 <- lm(data = df1,
           "cti ~ year * group")
summary(mod1)
mod1.2 <- lm(data = df1,
             "cti ~ year + group")
summary(mod1.2)
au1 <- broom::augment(mod1.2, df1, interval = "confidence")

p1 <- au1 %>%
  ggplot(aes(x = year, y = cti,color = group)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = .5) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = group),
            #  fill = "darkturquoise",
              alpha = .5,
            color = "transparent") +
  geom_line(aes(y = .fitted))+
  
  # now add points and 25 - 75 quantiles
  geom_segment(data = . %>% filter(group == unique(df1$group)[1]),
               aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15,
               position = position_nudge(x = .25)) +
  geom_segment(data = . %>% filter(group ==  unique(df1$group)[2]),
               aes(xend = year, y = cti_25, yend = cti_75),
               linewidth = .15,
               position = position_nudge(x = -.25)) +
  geom_point(data = . %>% filter(group ==  unique(df1$group)[1]),
             aes(y = cti), 
             shape = 16, size = 1, alpha = .7,
             position = position_nudge(x = .25)) +
  geom_point(data = . %>% filter(group ==  unique(df1$group)[2]),
             aes(y = cti), 
             shape = 16, size = 1, alpha = .7,
             position = position_nudge(x = -.25)) +
  
  # add annotations
  annotate(geom = "richtext",
           x = 2023, y = 10,
           label = glue::glue("<b>CTI Trend:</b> ", 
                              round(mod1.2$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(mod1.2)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
                              "<br>",
                              "<b>CTI Offset:</b> ", 
                              unique(df1$group)[[2]], " +", round(mod1.2$coefficients[["groupsessile"]],2),"°C, p ",
                              if(summary(mod1.2)$coefficients["groupsessile","Pr(>|t|)"] < 0.01){"< 0.01"},
                              "<br>",
                              "<b>Temp. Trend:</b> ",
                              round(tempmod$coefficients[["year"]],3)*10,
                              " °C / decade, p ",
                              if(summary(tempmod)$coefficients["year","Pr(>|t|)"] < 0.01){"< 0.01"},
           ),
           hjust = 1,vjust = .5,size = 3.25,label.color = NA,fill = NA) +
  
  # labels
  labs(x = NULL,
       y = "Community Thermal Index",
       color = "Motility Group",
       fill = "Motility Group") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1)) +
  scale_color_manual(values = c("cyan",
                                "darkcyan"),
                     labels = ~stringr::str_to_title(.)) +
  
  scale_fill_manual(values = c("cyan",
                                "darkcyan"),
                    labels = ~stringr::str_to_title(.))

p1 + ggview::canvas(4,4)






# -------------------------------------------------------------------------
# p3 all species in filtered quadrats -------------------------------------
mod3 <- lm(pa_therm_filtered, 
           formula = "mean_monthly_mean ~ year * group")
summary(mod3)
au3 <- broom::augment(mod3, pa_therm_filtered, interval = "confidence")

p3 <- au3 %>%
  ggplot(aes(x = year, y = mean_monthly_mean, color = group,
             fill = group)) +
  
  # first, add real temperature trend
  geom_ribbon(data = temptrend, inherit.aes =F,
              aes(x = year,  ymin = .lower, ymax = .upper),
              fill = "grey40", alpha = .8) +
  geom_line(data = temptrend,  inherit.aes = F,
            aes(x = year, y = .fitted),
            linetype = "dashed") +
  
  # then, add model output
  geom_point(aes(y = mean_monthly_mean),
             size = .25, alpha = .1, shape = 16,
             position = position_jitter(width = .3,height = .05)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = .8, color = "black") +
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
  coord_cartesian(xlim = c()) +
  theme(legend.position = "inside",
        legend.position.inside = c(.02,.98),
        legend.justification = c(0,1))

p3 + ggview::canvas(4,4)


