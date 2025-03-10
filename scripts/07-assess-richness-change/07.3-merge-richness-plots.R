# script to merge and save richness plots



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)

# upload rich by level plot and join --------------------------------------
rich_by_level <- readRDS("outputs/richness_change/p2.2.2_no_zeros.rds")
overall_rich <- readRDS("outputs/richness_change/overall_rich_plot.rds")

theme_set(ggthemes::theme_few())

full_rich_plot <- (overall_rich +
                     theme(panel.border  = element_rect(linewidth = .7)) +
                     labs(title = "a) Overall richness")
                   | rich_by_level + 
                     theme(panel.border  = element_rect(linewidth = .4)) +
                     labs(title = "b) Richness across shore levels")) 

full_rich_plot + ggview::canvas(width = 8, height = 4)

ggsave(full_rich_plot, 
       filename = "outputs/richness_change/full_rich_plot.png",
       width = 8,
       height = 4)

