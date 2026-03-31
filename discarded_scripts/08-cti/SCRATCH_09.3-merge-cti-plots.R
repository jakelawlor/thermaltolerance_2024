# merge total and by level CTI plots

# libraries
library(patchwork)
library(ggplot2)
library(ggtext)
library(dplyr)
theme_set(ggthemes::theme_few())


p1 <- readRDS("outputs/cti/total_cti_p1.rds") +
  labs(title = "Overall CTI Change",
       x = "Year")
p2 <- readRDS("outputs/cti/bylevel_cti_p1.rds") +
  labs(title = "CTI Change by Level") +
  theme(plot.title.position = "panel")

fullp <- p1 | p2
fullp + ggview::canvas(8,4)

ggsave(fullp,
       file = "outputs/cti/cti_plot_full.png",
       width = 8,
       height = 4)
