# Script to use obis records matched to temperatures
# to calculate thermal affinity of the species



# libraries ---------------------------------------------------------------
library(dplyr)
library(purrr)
library(ggplot2)
library(ggtext)


# upload data -------------------------------------------------------------
obis_recs_joined <- readRDS(
  here::here("data-processed",
             "obis_recs_matched_temp.rds")
)
obis_recs_number <- readr::read_csv(
  here::here("data-processed",
             "number_obis_recs_per_spp.csv")
)
obis_recs_number %>% filter(obis_recs_with_temp_match == 0)
# note that of 93 total species, 2 have no obis records
# that we could match temperatures to (Clathromorphum circumscriptum & Ptilota serrata)
obis_recs_number %>% filter(is.na(obis_recs_with_temp_match))
# and 2 had no obis recs that we could match a year to 
# (Micrura affinis &  Protectocarpus speciosus), so our total 
# thermal affinities are 89 species instead of 93
obis_recs_number$obis_recs_with_temp_match[obis_recs_number$obis_recs_with_temp_match!= 0] %>%
  quantile(c(.05,.5,.95), na.rm=T)

# find affinities ---------------------------------------------------------
# we'll find a range of values here:
# 5th percentile, mean, median, 95th percentile
# of coldest month, warmest month, and averaged month SSTs
obis_derived_affinities <- purrr::map(
  .x = obis_recs_joined,
  .f = ~.x %>%
    filter(!is.na(mean_monthly_sst)) %>%
    summarize(perc_monthly_min_05 = quantile(min_monthly_sst,.05),
              median_monthly_min = median(min_monthly_sst),
              mean_monthly_min = mean(min_monthly_sst),
              perc_monthly_min_95 = quantile(min_monthly_sst, .95),#min_monthly_min = min(min_monthly_sst),
              #mean_monthly_min = mean(min_monthly_sst),
              #median_monthly_min = median(min_monthly_sst),
              #min_monthly_mean = min(mean_monthly_sst),
              perc_monthly_mean_05 = quantile(mean_monthly_sst,.05),
              median_monthly_mean = median(mean_monthly_sst),
              mean_monthly_mean = mean(mean_monthly_sst),
              perc_monthly_mean_95 = quantile(mean_monthly_sst, .95),
              #max_monthly_mean = max(mean_monthly_sst),
              #mean_monthly_max = mean(max_monthly_sst),
              #median_monthly_max = median(max_monthly_sst),
              #max_monthly_max = max(max_monthly_sst),
              perc_monthly_max_05 = quantile(max_monthly_sst,.05),
              median_monthly_max = median(max_monthly_sst),
              mean_monthly_max = mean(max_monthly_sst),
              perc_monthly_max_95 = quantile(max_monthly_sst, .95),
              n_records = n())
) %>%
  bind_rows(.id = "organism") %>%
  # rank by monthly median
  mutate(organism = forcats::fct_reorder(organism, median_monthly_mean)) %>%
  arrange(organism) %>%
  # add a numeric rank value
  mutate(spp_rank = c(1:n())) 

obis_derived_affinities2 <- obis_derived_affinities %>%
  # separate genus and species
  tidyr::separate(organism, into=c("genus","species")) %>%
  
  # keep only first letter of genus
  mutate(genus = paste0(substr(genus,1,1),".")) %>%
  
  # make a G. species column, with asterisks for markdown code
  mutate(organism = paste0("*",genus," ",species,"*")) %>%
  
  # delete genus and species separate columns
  dplyr::select(-genus, -species) %>%
  
  # add rank, as the axis will be numeric,
  # and label, numbering the abbreviated spp names
  mutate(label = paste0(spp_rank,"\\. ",organism)) %>% 
  
  # turn label to ordered factor
  mutate(label = forcats::fct_reorder(label, spp_rank))




## plot ------------------------------
color1 <- "#782391"

p <- obis_derived_affinities2 %>% 
  ggplot(aes(x=spp_rank, fill = label)) +
  
  geom_point(aes(y=perc_monthly_mean_05), 
             color = color1, alpha=.5, size = .5)+
  geom_point(aes(y=perc_monthly_mean_95),  
             color = color1, alpha=.5, size = .5) +
  geom_point(aes(y=median_monthly_mean),  
             shape=17, color = color1, alpha=1, size=1 ) +
  geom_point(aes(y=mean_monthly_mean), 
             shape=8, color = color1, alpha=1, size=1, stroke=.4 ) +
 # theme(legend.position = "none") +
  geom_segment(aes(xend=spp_rank,
                   y=perc_monthly_mean_05,
                   yend=perc_monthly_mean_95),  
               color = color1, alpha=.5, linewidth = .5) +
  ggthemes::theme_few(base_size = 11) +
 # theme(legend.position = "none") +
  
  # add a fake legend
  annotate(geom = "rect",
           xmin = 4,
           xmax = 39,
           ymin = 17,
           ymax = 25,
           color = "white",
           linewidth=.5,
           fill = "white"
  )+
  annotate(geom = "point",
           x = 5,
           y = 24, 
           color = color1,
           alpha = .5) +
  annotate(geom = "richtext",
           x = 7,
           y = 24, 
           label = "95<sup>th</sup> Percentile<br>Occurrence Temperature",
           lineheight = .7, 
           vjust = .5, 
           hjust = 0,
           size = 2.75,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  
  
  annotate(geom = "point",
           x = 5,
           y = 18, 
           color = color1,
           alpha = .5) +
  annotate(geom = "richtext",
           x = 7,
           y = 18, 
           label = "5<sup>th</sup> Percentile<br>Occurrence Temperature",
           lineheight = .7, 
           vjust = .5, 
           hjust = 0,
           size = 2.75,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt") ) +
  
  annotate(geom = "segment",
           x = 5,
           xend = 5,
           y = 18, 
           yend = 24,
           color = color1,
           alpha = .5) +
  
  annotate(geom = "point",
           shape=17,
           x = 5,
           y = 20.5, 
           color = color1, 
           size = 2) +
  annotate(geom = "richtext",
           x = 7,
           y = 20.5, 
           label = "Median Occurrence Temperature",
           lineheight = .7, 
           vjust = .5, 
           hjust = 0,
           size = 2.75,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  
  annotate(geom = "point",
           shape=8,
           x = 5,
           y = 22, 
           color = color1, 
           size = 2) +
  annotate(geom = "richtext",
           x = 7,
           y = 22, 
           label = "Mean Occurrence Temperature",
           lineheight = .7, 
           vjust = .5, 
           hjust = 0,
           size = 2.75,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  
  scale_y_continuous(breaks = c(10,20,30),
                     labels = c("10°","20°","30°")
  ) +
  coord_cartesian(xlim = c(-.5, 89.5), expand = F, clip = "off") +
  labs(y= "Thermal Affinity (°C)",
       x = "Species",
       title = "Occupancy-Derived Thermal Affinities",
       subtitle = "Summary values of mean monthly sea surface temperature for OBIS occurrence records"
  ) +
  theme(legend.text = ggtext::element_markdown(size=7),
        legend.key = element_blank(),
        # Change legend key size and key width
        legend.key.size = unit(0, "cm"),
        legend.key.width = unit(0,"cm"),
        legend.key.height = unit(.25,"cm"),
        legend.background = element_blank(),
        legend.box.margin = margin(0,0,0,-5, unit = "pt"),
        legend.margin = margin(5,0,0,0, unit="pt")) +
  guides(fill=guide_legend(ncol=3,
                           title = NULL,
                           override.aes=list(size = 0,
                                             color = "transparent"))) +
  
  theme(plot.margin = margin(l=0,r=0,b=0,t=5,unit="pt"),
        panel.border = element_blank(),
        plot.title.position = "plot",
        axis.line.x = element_line(linewidth  = .5),
        axis.ticks = element_line(linewidth = .5),
        panel.grid.major.y = element_line(linewidth = .25, 
                                          color = "grey20", 
                                          linetype = "longdash")) 

  
p + ggview::canvas(height = 4, width = 8.5)


ggsave(p,
       file = here::here("outputs",
                         "Fig1_thermal_affinties.png"),
       height = 4,
       width = 8.5,
       unit = "in"
       )


readr::write_csv(
  obis_derived_affinities,
  here::here("data-processed",
             "spp_thermal_affinities.csv")
)

rm(list  = ls())
