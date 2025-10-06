library(ggplot2)
p1 <- readRDS("outputs/cti/cti_by_abundance_plot.rds") +
  labs(title = "Total CTI with abundance data")
p2 <- readRDS("outputs/cti/p_by_level_by_abundance.rds") 


library(patchwork)


pfull <- (p1 | p2) + plot_layout(widths = c(1,1))
pfull +
  ggview::canvas(8,5)

ggsave(pfull,
       file = "outputs/cti/p_by_abundance_full.png",
       width = 8,
       height = 5
       )

ggsave(pfull,
       file = "outputs/cti/p_by_abundance_full.pdf",
       width = 8,
       height = 5,
       device = cairo_pdf()
)
