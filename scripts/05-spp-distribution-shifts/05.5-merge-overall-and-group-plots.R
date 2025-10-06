# upload and merge species proportion depth shifts figures

library(ggplot2)
library(patchwork)
library(dplyr)

all <- readRDS("outputs/depth_shifts/depth_proportions.rds")
gr <- readRDS("outputs/depth_shifts/depth_proportions_by_group.rds")


# the legend themes differ slightly, preventing guides from collecting in patchwork
# extract legend theme here
legend_theme <- list(
  theme_minimal(),
  guides(fill = 
           guide_legend(theme = theme(legend.title.position = "top",
                                      legend.margin = ggplot2::margin(r = 20,t = -15,0,0, unit= "pt"))))

)



# merge all unfiltered ----------------------------------------------------
p_unfilt <- ((
  all$all_sp_unfiltered + legend_theme +
    labs(title = "a) Both datasets")|
    gr$all_count + legend_theme +
    labs(title = "b) Density only")|
    gr$all_cover + legend_theme + 
    labs(title = "c) Percent cover only")
  
) #&
 # scale_y_discrete(labels = c("Min (5%)\nOccurrence Height","Mean\nOccurrence Height","Median\nOccurrence Height","Max (95%)\nOccurrence Height"))
  )+
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom",
        axis.text.y = element_text(lineheight = .8,
                                   size = 7),
        plot.title = element_text(size = 10),
        panel.grid = element_blank(),
        plot.margin = margin(r = 20,0,0,0,"pt"))


p_unfilt + ggview::canvas(8,2.9)

ggsave(p_unfilt,
       filename = "outputs/depth_shifts/depth_proportions_unfiltered.pdf",
       width = 8, 
       height = 2.9,
       device = cairo_pdf())




# merge all filtered ------------------------------------------------------
# merge plots with species at domain margins excluded
p_inside_domain <- ((
  all$all_sp_inside_domain + legend_theme +
    labs(title = "a) Both datasets")|
    gr$in_domain_count + legend_theme +
    labs(title = "b) Density only")|
    gr$in_domain_cover + legend_theme + 
    labs(title = "c) Percent cover only")
  
) #&
# scale_y_discrete(labels = c("Min (5%)\nOccurrence Height","Mean\nOccurrence Height","Median\nOccurrence Height","Max (95%)\nOccurrence Height"))
)+
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom",
        axis.text.y = element_text(lineheight = .8,
                                   size = 7),
        plot.title = element_text(size = 10),
        panel.grid = element_blank(),
        plot.margin = margin(r = 20,0,0,0,"pt"))



p_inside_domain + ggview::canvas(8,2.5)

ggsave(p_inside_domain,
       filename = "outputs/depth_shifts/depth_proportions_inside_domain.pdf",
       width = 8, 
       height = 2.5,
       device = cairo_pdf())



