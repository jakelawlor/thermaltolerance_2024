
library(ggplot2)

library(sf)
world <- rnaturalearth::ne_countries(country = c("United States of America",
                                                 "Canada","Mexico"),
                                     returnclass = "sf",
                                    # scale = "large"
                                    )

p <- world %>%
  ggplot() +
  geom_sf(fill = "lightgreen",
          alpha = .5) +
  coord_sf(xlim = c(-100, -50),
           ylim =  c(25, 70)) + 
  theme_void() 


p + ggview::canvas(
               width = 3.75,
               height = 5, unit = "in")
ggsave(p,
       filename = here::here("shoals-pres-figs",
                             "map.png"),
       width = 3.75,
       height = 5,
       unit = "in")
