
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# upload raw temperature values in the study area
temp <- read.csv(here::here(
  "data-processed",
  "appledore-island-env-data",
  "appledore_temps_1982_2023.csv"
))


## plot min, mean, max   ------------------------------
color_max <- "#F21A00"
color_mean <- "#EBCC2A"
color_min <- "#3B9AB2"


minmod <- lm(data = temp, min ~ sample_year)
summary(minmod)
minmod_au <- broom::augment(minmod, 
                            interval = "confidence")

meanmod <- lm(data = temp, mean ~ sample_year)
summary(meanmod)
meanmod_au <- broom::augment(meanmod, 
                            interval = "confidence")

maxmod <- lm(data = temp, max ~ sample_year)
summary(maxmod)
maxmod_au <- broom::augment(maxmod, 
                             interval = "confidence")




# format data -------------------------------------------------------------
temp_pivot <- temp %>%
  select(-median) %>%
  tidyr::pivot_longer(cols = c(mean,min, max),
                      names_to = "metric",
                      values_to = "temp")


au_pivot <- 
  meanmod_au %>% rename("temp" = "mean") %>%
  mutate(metric = "mean") %>%
  rbind(
    minmod_au %>% rename("temp" = "min") %>%
      mutate(metric = "min")
  ) %>%
  rbind(
    maxmod_au %>% rename("temp" = "max") %>%
      mutate(metric = "max")
  )


summary(maxmod)
summary(meanmod)
summary(minmod)

labels_df <- data.frame(
  metric = c("max","mean","min"),
  trend = c(round(maxmod$coefficients[["sample_year"]]*10,2),
            round(meanmod$coefficients[["sample_year"]]*10,2),
            round(minmod$coefficients[["sample_year"]]*10,2)
  ),
  std_err = c(round(summary(maxmod)$coefficients["sample_year","Std. Error"]*10,2),
              round(summary(meanmod)$coefficients["sample_year","Std. Error"]*10,2),
              round(summary(minmod)$coefficients["sample_year","Std. Error"]*10,2)
  ),
  p = c("p < 0.01", "p < 0.01", "p < 0.05"),
  last_y = c(au_pivot[au_pivot$sample_year == 2023 & au_pivot$metric == "max",]$.fitted,
             au_pivot[au_pivot$sample_year == 2023 & au_pivot$metric == "mean",]$.fitted,
             au_pivot[au_pivot$sample_year == 2023 & au_pivot$metric == "min",]$.fitted)
) %>%
  mutate(label = paste0(trend," ± ",std_err,"°","\nper decade",
                        "\n", p))

# plot --------------------------------------------------------------------
temp_p <- temp_pivot %>%
  ggplot(aes(x = sample_year,
              y = temp,
              color = metric,
              fill = metric)) +
  geom_line() +
  geom_point() +
  geom_ribbon(data = au_pivot,
              aes(ymin = .lower,
                  ymax = .upper),
              alpha  = .25,
              color = "transparent") +
  scale_alpha_identity() +
  geom_line(data = au_pivot,
            aes(y = .fitted),
            linetype = "dashed") +
  
  # add labels
  geom_label(data = labels_df, 
             aes(x = 2024.5,
                 y = last_y,
                 label = label),
            hjust = 0,
            #fill = "white",
            alpha = .3,
            fontface = "bold",
            size = 3,
            lineheight = .88,
            label.size = .1,
            show.legend = F,
            color = "black"
            ) +
  
  scale_color_manual(values = c(color_max, color_mean, color_min),
                     labels = c("Warmest Month", "Annual Mean", "Coldest Month")) +
  scale_fill_manual(values = c(color_max, color_mean, color_min),
                    labels = c("Warmest Month", "Annual Mean", "Coldest Month")) +
  scale_x_continuous(breaks = c(1980, 1990, 2000, 2010, 2020)) +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  coord_cartesian(xlim = c(1982, 2035),
                  ylim = c(0, 22)) +
  labs(x = NULL,
       y = "Sea Surface Temperature (°C)",
       fill = NULL,#"Temperature Metric",
       color = NULL,#"Temperature Metric"
       ) +
  
  theme(legend.position = "inside",
        legend.position.inside = c(.5, .01),
        legend.justification = c(.5, 0),
        legend.direction = "horizontal",
        panel.grid = element_blank(),
        legend.box.margin = margin(t = -5,b = -10,0,0,"pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top")),
         color = guide_legend(theme = theme(legend.title.position = "top"))) 

temp_p + ggview::canvas(4,3)  
ggsave("outputs/env_plots/appledore_temp_plot.png",
       temp_p,
       width = 4, 
       height = 3)
