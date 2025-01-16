# make obis recs match example fig
library(terra)
library(tidyterra)
library(ggplot2)

# test_spp = 
# chondrus crispus
# fucus distichus
# carcenus maenas
# nucella lapilus 

test <- readRDS("data-processed/obis_recs_matched_temp.rds")
test[["Chondrus crispus"]] %>% distinct(year) %>% arrange(year)

test %>% bind_rows(.id = "organism") %>% count(year) %>% arrange(desc(year))

temps <- readRDS("data-processed/global_temps_1982_2023.rds")


rast_85 <- temps[["1985"]] %>% 
  select(longitude, latitude,mean_monthly_sst) %>%
  terra::rast()
world <- rnaturalearth::ne_countries(returnclass = "sf")
crs(rast_85) <- crs(world)

theme_set(theme_bw(base_size = 6) +
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.ticks.length = unit(0,"mm"),
                  plot.background = element_blank(),
                  plot.margin = margin(0,0,0,0),
                  panel.spacing = margin(0,0,0,0),
                  panel.border = element_rect(linewidth = .5,
                                              color = "grey90")))

baseplot <- ggplot() +
  geom_spatraster(data = rast_85,
                  #alpha = .8
                  ) +
  scale_fill_viridis_c(#option = "inferno",
                       na.value = "transparent") +
  geom_sf(data = world, fill = "grey98",
          color = "grey50", linewidth = .1) +
  theme(legend.position = "none") +
  coord_sf(expand = F) +
  labs(x = NULL,
       y = NULL)
  



chondrus_plot <- baseplot + 
  geom_point(data = test[["Chondrus crispus"]] %>%filter(year == 1985)
             ,
             shape = 21,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             size = 1.5, 
             alpha = .5,
             #shape = 16,
             fill = "#FF3131",
             stroke = .05) 



fucus_plot <- baseplot + 
  geom_point(data = test[["Fucus distichus"]] %>%filter(year == 1985)
             ,
             shape = 21,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             size = 1.5, 
             alpha = .5,
             #shape = 16,
             fill = "#FF3131",
             stroke = .05) 


carcinus_plot <- baseplot + 
  geom_point(data = test[["Carcinus maenas"]] %>%filter(year == 1985)
             ,
             shape = 21,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             size = 1.5, 
             alpha = .5,
             #shape = 16,
             fill = "#FF3131",
             stroke = .05) 



nucella_plot <- baseplot + 
  geom_point(data = test[["Nucella lapillus"]] %>%filter(year == 1985)
             ,
             shape = 21,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             size = 1.5, 
             alpha = .5,
             #shape = 16,
             fill = "#FF3131",
             stroke = .05) 



# save plots --------------------------------------------------------------
carcinus_plot + ggview::canvas(2,1)
chondrus_plot + ggview::canvas(2,1)
fucus_plot + ggview::canvas(2,1)
nucella_plot + ggview::canvas(2,1)

ggsave(carcinus_plot, 
       filename = "outputs/methods_fig_examples/carcinus_plot.svg",
       width = 2,
       height = 1,
       unit = "in",
       device = svglite::svglite
       )

ggsave(chondrus_plot, 
       filename = "outputs/methods_fig_examples/chondrus_plot.svg",
       width = 2,
       height = 1,
       unit = "in",
       device = svglite::svglite
)

ggsave(fucus_plot, 
       filename = "outputs/methods_fig_examples/fucus_plot.svg",
       width = 2,
       height = 1,
       unit = "in",
       device = svglite::svglite
)


ggsave(nucella_plot, 
       filename = "outputs/methods_fig_examples/nucella_plot.svg",
       width = 2,
       height = 1,
       unit = "in",
       device = svglite::svglite
)


