# Script to use obis records matched to temperatures
# to calculate thermal affinity of the species



# libraries ---------------------------------------------------------------
library(dplyr)
library(purrr)
library(ggplot2)
library(ggtext)
library(patchwork)


# upload data -------------------------------------------------------------
obis_recs_joined <- readRDS(
  here::here("data-processed",
             "species-thermal-affinities",
             "obis_recs_matched_temp.rds")
)
obis_recs_number <- readr::read_csv(
  here::here("data-processed",
             "species-thermal-affinities",
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
i <- 20
bw = .5
color_max <- "#F21A00"
color_mean <- "#EBCC2A"
color_min <- "#3B9AB2"

example_p <- obis_recs_joined[[i]] %>% 
  ggplot() +
  geom_density(aes(x = mean_monthly_sst), fill = color_mean, alpha = .5,
               bw = bw) +
  geom_density(aes(x = min_monthly_sst), fill = color_min, alpha = .5,
               bw = bw) +
  geom_density(aes(x = max_monthly_sst), fill = color_max, alpha = .5,
               bw = bw) +
  # add mean point and label
  geom_point(data = . %>% summarize(mean_mean = mean(mean_monthly_sst, na.rm=T)),
             aes(x = mean_mean, y = 0), fill = color_mean, size = 4, shape = 21) +
  geom_label(data = . %>% summarize(mean_mean = mean(mean_monthly_sst, na.rm=T)),
             aes(x = mean_mean, y = 0.08, label = "Mean of\nannual means"),
             fill = colorspace::lighten(color_mean,.2), size = 2.5, lineheight = .8) +
  # add min point and label
  geom_point(data = . %>% summarize(min_05 = quantile(min_monthly_sst,.05, na.rm=T)),
             aes(x = min_05, y = 0), fill = color_min, size = 4, shape = 24)  +
  geom_label(data = . %>% summarize(min_05 = quantile(min_monthly_sst,.05, na.rm=T)),
             aes(x = min_05, y = 0.08, label = "5th % of\ncoldest months"), 
             fill = colorspace::lighten(color_min,.2), size = 2.5, lineheight = .8)  +
  # add max point and label
  geom_point(data = . %>% summarize(max_95 = quantile(max_monthly_sst,.95, na.rm=T)),
             aes(x = max_95, y = 0), fill = color_max, size = 4, shape = 22) +
  geom_label(data = . %>% summarize(max_95 = quantile(max_monthly_sst,.95, na.rm=T)),
             aes(x = max_95, y = 0.08, label = "95th % of\nwarmest months"), 
             fill = colorspace::lighten(color_max,.2), size = 2.5,  lineheight = .8) +
  labs(title = paste0("Example species: ",names(obis_recs_joined)[i]),
       subtitle = paste0("n records: ", scales::comma(nrow(obis_recs_joined[[i]]))),
       x = "Sea Surface Temperature (°C)",
       y = "Density") +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(labels = ~paste0(.,"°")) +
  geom_density(data = data.frame(x = 0, y = 0, 
                               group = c("Warmest Month of Occurrences",
                                         "Annual Mean of Occurrences",
                                         "Coldest Month of Occurrences")),
             aes(x=x,y=y, fill = group), alpha = 0) +
  scale_fill_manual(
    breaks = c(
      "Warmest Month of Occurrences",
      "Annual Mean of Occurrences",
      "Coldest Month of Occurrences"),
    labels = c("Warmest Month SST of Occurrences",
               "Annual Mean SST of Occurrences",
               "Coldest Month SST of Occurrences"),
    values = c(color_max,
               color_mean,
               color_min)) +
  labs(fill = NULL) +
  coord_cartesian(clip = "off", expand = F) +
  theme(panel.border = element_blank(),
        plot.title.position = "plot") +
  guides(fill = guide_legend(override.aes = list(alpha = .65))) 
  
example_p + ggview::canvas(8.5, 2)

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
  mutate(organism = forcats::fct_reorder(organism, mean_monthly_mean)) %>%
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





## plot min, mean, max   ------------------------------
color_max <- "#F21A00"
color_mean <- "#EBCC2A"
color_min <- "#3B9AB2"
bump <- 8
bump2 <- -4

p_min_mean_max <- obis_derived_affinities2 %>%  
  
  ggplot(aes(x=spp_rank, fill = label)) +
  
  ggthemes::theme_few(base_size = 11) +
  
  # add segment from min to max
  geom_segment(aes(xend = spp_rank,
                   y = perc_monthly_min_05, 
                   yend = perc_monthly_max_95),
               linewidth = .1) + 


  # add mean thermal affinity: mean of monthly means
  geom_point(aes(y=mean_monthly_mean), 
             shape=19, color = color_mean, alpha=1, size=1, stroke=.4 ) +

  # add thermal minimum: 5th % of monthly min
  geom_point(aes(y=perc_monthly_min_05), 
             shape=17, color = color_min, alpha=1, size=1, stroke=.4 ) +
  
  # add thermal maxmimum: 95th % of monthly max
  geom_point(aes(y=perc_monthly_max_95), 
             shape=15, color = color_max, alpha=1, size=1, stroke=.4 ) +

  
  # add a fake legend
  annotate(geom = "rect",
           xmin = 4 + bump2,
           xmax = 47 + bump2,
           ymin = 17+bump,
           ymax = 25 + bump,
           color = "white",
           linewidth=.5,
           fill = "white",
           alpha = .8
  )+
  
  annotate(geom = "segment",
           x = 5 + bump2,
           xend = 5 + bump2,
           y = 18 + bump, 
           yend = 24 + bump,
           color = "black",
           alpha = .5) +
  
  annotate(geom = "point",
           x = 5 + bump2,
           y = 24 + bump, 
           color = color_max,
           alpha = .9,
           shape = 15) +
  annotate(geom = "richtext",
           x = 7 + bump2,
           y = 24 + bump, 
           label = "<b>Thermal Maximum</b><br><span style = 'font-size:7pt;'>95<sup>th</sup> Percentile Warmest Month of Occurrences</span>",
           lineheight = .7, 
           vjust = .5, 
           hjust = 0,
           size = 2.75,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  
  
  annotate(geom = "point",
           x = 5 + bump2,
           y = 18 + bump, 
           color = color_min,
           alpha = .9,
           shape = 17) +
  annotate(geom = "richtext",
           x = 7 + bump2,
           y = 18+bump, 
           label =  "<b>Thermal Minimum</b><br><span style = 'font-size:7pt;'>5<sup>th</sup> Percentile Coldest Month of Occurrences</span>",
           lineheight = .7, 
           vjust = .5, 
           hjust = 0,
           size = 2.75,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt") ) +
  
  
  annotate(geom = "point",
           shape=19,
           x = 5 + bump2,
           y = 21 + bump, 
           color = color_mean, 
           size = 2) +
  annotate(geom = "richtext",
           x = 7 + bump2,
           y = 21 + bump, 
           label =  "<b>Mean Thermal Affinity</b><br><span style = 'font-size:7pt;'>Mean of Annual Mean Temp. of Occurrences</span>",
           lineheight = .7, 
           vjust = .5, 
           hjust = 0,
           size = 2.75,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt")) +
  
  
  # add summary violins and boxplots
  geom_violin(aes(x = 93, 
                  y = perc_monthly_min_05),
              fill = color_min, 
              color = "transparent",
              width = 5,
              alpha = .8) +
    geom_violin(aes(x = 93, 
                  y = perc_monthly_max_95),
              fill = color_max, 
              color = "transparent",
              width = 5, alpha = .8) +
  geom_violin(aes(x = 93, 
                  y = mean_monthly_mean),
              fill = color_mean, 
              color = "transparent",
              width = 5, alpha = .8) +
  
  geom_boxplot(aes(x = 93, 
                   y = perc_monthly_min_05),
               fill = "transparent", 
               color = "black",
               linewidth = .2,
               outliers = F,
               width = 1) +
  geom_boxplot(aes(x = 93, 
                   y = perc_monthly_max_95),
               fill = "transparent", 
               color = "black",
               linewidth = .2,
               outliers = F,
               width = 1) +
  geom_boxplot(aes(x = 93, 
                   y = mean_monthly_mean),
               fill = "transparent", 
               color = "black",
               linewidth = .2,
               outliers = F,
               width = 1) +
  

  
  scale_y_continuous(breaks = c(0,10,20,30),
                     labels = c("0°","10°","20°","30°")
  ) +
  coord_cartesian(xlim = c(-.5, 94), expand = F, clip = "off") +
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
                                             alpha = 0,
                                             color = "transparent"))) +
  
  theme(plot.margin = margin(l=0,r=0,b=0,t=5,unit="pt"),
        panel.border = element_blank(),
        plot.title.position = "plot",
        axis.line.x = element_line(linewidth  = .5),
        axis.ticks = element_line(linewidth = .5),
        panel.grid.major.y = element_line(linewidth = .25, 
                                          color = "grey20", 
                                          linetype = "longdash")) 


p_min_mean_max + ggview::canvas(height = 4, width = 8.5)


ggsave(p_min_mean_max,
       file = here::here("outputs",
                         "Fig1_thermal_affinties_min_mean_max.png"),
       height = 4,
       width = 8.5,
       unit = "in"
)




# merge with example ------------------------------------------------------
p_therms_and_example <- ((p_min_mean_max +
   labs(title = "a) Occupancy-derived Thermal Affinities",
        subtitle = NULL))/
   (example_p + 
      theme(legend.justification = "left") +
   labs(title = paste0("b) Example Species: ",names(obis_recs_joined)[[i]]," (n records = ", scales::comma(nrow(obis_recs_joined[[i]])),")"),
        subtitle = NULL))) +
  plot_layout(heights = c(1,.2)) 

p_therms_and_example +  ggview::canvas(8.5,5.5)
ggsave(p_therms_and_example,
       file = here::here("outputs",
                         "Fig1_thermal_affinties_min_mean_max_and_example.png"),
       height = 5.5,
       width = 8.5,
       unit = "in"
)



readr::write_csv(
  obis_derived_affinities,
  here::here("data-processed",
             "species-thermal-affinities",
             "spp_thermal_affinities.csv")
)


rm(list  = ls())
