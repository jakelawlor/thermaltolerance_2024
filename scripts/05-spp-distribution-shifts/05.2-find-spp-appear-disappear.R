# find if species that enter are warm-affinity and 
# species that leave are cold-affinity as a whole



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(multcompView)
library(ggtext)

# data --------------------------------------------------------------------
# upload full sampling to get year groups
pa <- readr::read_csv(
  here::here(
    "data-processed",
    "appledore-survey-data",
    "pa-with-therm",
    "pa-with-therm-all.csv"
  )
) %>%  filter(!is.na(mean_monthly_mean))


# use the highly-sampled quadrats only to standardize sampling
pa_therm_filtered <- readr::read_csv(
  here::here(
    "data-processed",
    "appledore-survey-data",
    "pa-with-therm",
    "pa-with-therm-highly-sampled.csv"
  )
) %>% filter(!is.na(mean_monthly_mean))


# upload traits to see how many of each species are in each group ---------
traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed_2.csv") %>%
  mutate(trait_group = case_when(motility_adult == "sessile" & group == "Invertebrate" ~ "Sessile Invertebrate",
                               motility_adult == "motile" & group == "Invertebrate" ~ "Motile Invertebrate",
                               group == "Algae" ~ "Algae"
  )) %>%
  rename(organism = gen_spp)

pa <- pa %>% left_join(traits)
pa_therm_filtered <- pa_therm_filtered %>% left_join(traits)


# separate by year group
unique <- unique(pa$year) %>% sort()
y1 <- unique[1:10]#1980:1989
y2 <- unique[11:20]#1990:1999
y3 <- unique[21:30]#2000:2009
y4 <- unique[31:40]#2010:2023
#y1 <- 1982:1991#1980:1989
#y2 <- 1992:2001#1990:1999
#y3 <- 2002:2012#2000:2009
#y4 <- 2013:2023#2010:2023

pa_therm_filtered <- pa_therm_filtered %>% 
  mutate(tidalheight = (13-level)*.348)  %>%
  mutate(year_group = case_when(
    year %in% y1 ~ 1,
    year %in% y2 ~ 2,
    year %in% y3 ~ 3,
    year %in% y4 ~ 4
  ))
#rm(unique, y1,y2,y3,y4, pa)



# find entering species ---------------------------------------------------
pa_across_years <- pa_therm_filtered %>% 
  filter(replicate == 1) %>%
  group_by(organism, mean_monthly_mean, perc_monthly_min_05, mean_monthly_min, perc_monthly_max_95, mean_monthly_max, trait_group, invasive_gom) %>% 
  distinct(year_group) %>%
  mutate(present = TRUE) %>%
  ungroup() %>%
  mutate(organism = forcats::fct_reorder(organism, mean_monthly_mean)) 

pa_across_years %>%  
  ggplot(aes(x = year_group,
             y = organism)) +
  geom_tile() +
  coord_equal()

pa_pivot <- pa_across_years %>%
  tidyr::pivot_wider(names_from = year_group,
                     names_prefix = "yg_",
                     values_from = present)  %>%
  mutate(org_group = case_when(
    # first, if species is present in first 2 groups, but not last 2, consider disappearing
    (yg_1 == TRUE ) & (is.na(yg_4)) ~ "Disappearing",
    # second, if species is absent in first 2 groups, but present in last 2, added
    (is.na(yg_1)) & (yg_4 == TRUE) ~ "Arriving",
    yg_1 == T & yg_2 == T & yg_3 == T & yg_4 == T ~ "Persistent",
    TRUE ~ "Rare / Other"
    
  ))

pa_pivot2 <- pa_pivot %>%  
  tidyr::pivot_longer(cols = c("yg_1":"yg_4"),
                      names_to = "year_group",
                      values_to = "present") %>%
  tidyr::drop_na(present) %>%
  # add asterisk to organisms that are invasive
  mutate(organism_label = paste0("<i>",organism,"</i>")) %>%
  mutate(organism_label = case_when(invasive_gom == "yes" ~ paste0(organism_label,"<b> *</b>"),
                                    TRUE ~ organism_label)) %>%
  mutate(organism_label = forcats::fct_reorder(organism_label, mean_monthly_mean)) %>%
  mutate(org_group = factor(org_group,
                            levels = c("Disappearing","Persistent","Rare / Other","Arriving")))


# make fill scale that will apply to all plots
fill_scale <- {
  scale_fill_manual(values = c(Arriving = "magenta3",
                               Disappearing = "darkturquoise",
                               Persistent = "grey70",
                               `Rare / Other` = "grey85")) 
}

p_time <- pa_pivot2 %>%
  ggplot(aes(x = year_group,
             y = organism_label, 
             fill = org_group)) +
  geom_tile() +
  geom_segment(aes(x = "yg_4",
                   xend = "yg_4",
                   y = 1,
                   yend = length(unique(organism))),
               position = position_nudge(x=  1.5),
               arrow = arrow(ends = "both")) +
  annotate(geom = "text",
           x = 5.8,
           y = 15,
           angle = -90,
           label = "Cooler Affinity",
           vjust = 0,
           hjust = .5,
           size = 5) +
  annotate(geom = "text",
           x = 5.8,
           y = 65,
           angle = -90,
           label = "Warmer Affinity",
           vjust = 0,
           hjust = .5,
           size = 5) +
  fill_scale +
  scale_x_discrete(breaks = c("yg_1","yg_2","yg_3","yg_4"),
                   labels = c("1982-91","1992-01","2002-12","2013-23")
                   ) +
  ggthemes::theme_few() +
  annotate(geom = "rect", 
           xmin = .5, 
           xmax = 4.5,
           ymin = .5,
           ymax = length(unique(pa_pivot2$organism))+.5,
           color = "black",
           fill = "transparent") +
  coord_cartesian(xlim = c(0.5,7),
                  expand = F) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1),
        legend.position = "none",
        panel.border = element_blank(),
        axis.text.y = element_markdown())  +
  labs(x = NULL,
       y = NULL)

p_time + ggview::canvas(width = 4, height = 12, unit = "in") 
  

pa_pivot2_summary <- pa_pivot2 %>% 
  distinct(organism, across(contains("monthly_")), org_group, trait_group, invasive_gom)


pa_pivot2_summary %>% 
  count(org_group)


pa_pivot2_summary %>% 
  count(org_group, invasive_gom) %>%
  tidyr::pivot_wider(names_from = invasive_gom,
                     values_from = n) %>% 
  mutate(total = no  + yes,
         percent_invasive = yes / total * 100)

# bar plot

n_plot <- pa_pivot2_summary %>% 
  count(org_group, invasive_gom) %>%
  tidyr::pivot_wider(names_from = invasive_gom,
                     values_from = n) %>% 
  mutate(total = no  + yes,
         percent = paste0(round(yes/ total * 100,0),"%")) %>%
  ggplot(aes(x = total,
             y = org_group,
             fill = org_group)) +
  theme_void()+
  theme(#axis.line.y = element_line(),
        #axis.text.y = element_text(hjust = 1,
        #                           margin = margin(r = 1,0,0,0,"mm")),
        legend.position = "none") +
  coord_cartesian(expand = F,
                  ylim = c(.5,5),
                  xlim = c(-10, 37),
                  clip = "off") +
  geom_col(width = .85) +
  geom_text(aes(label = org_group,
                x = -1),
            hjust = 1,
            lineheight = .85) +
  geom_text(aes(label = total),
            hjust = -.15,
            vjust = .4,
            stat = "unique",
            size = 3,
            fontface = "bold") +
  fill_scale +
  geom_col(aes(x = yes),
           width = .85,
           fill = "black",
           alpha = .4) +
  geom_text(aes(label = yes,
                x = yes),
            hjust = -.35,
            vjust = .4,
            stat = "unique",
            size = 2.75) +
  annotate(geom = "segment",
           x = 0,
           xend = 0,
           y = .5,
           yend = 4.5) +
  annotate(geom = "text",
           x = 10.1,
           y = 5,
           label = "n Non-native*",
           hjust = 0,
           size = 2.75) +
  annotate(geom = "curve",
           x = 9.9,
           xend = 8.2,
           y = 5,
           yend = 4.4,
           curvature = .3,
           linewidth = .3) +
  annotate(geom = "text",
           x = 27,
           y = 3.25,
           label = "n Total",
           hjust = 0.5,
           size = 2.75) +
  annotate(geom = "curve",
           x = 27,
           xend = 22,
           y = 3.6,
           yend = 4.05,
           curvature = .2,
           linewidth = .3)

n_plot +ggview::canvas(4,2)

pa_pivot2_summary %>% 
  count(org_group, trait_group) %>%
  group_by(org_group) %>%
  mutate(perc = n/sum(n)*100) %>%
  group_split(org_group)

org_types <- pa_pivot2_summary %>%
  select(org_group) %>%
  mutate(trait_group = "All")
org_types <- org_types %>% rbind(pa_pivot2 %>% select(org_group, trait_group))
org_types %>%  
  mutate(trait_group = factor(trait_group,
                              levels = c("All","Algae","Motile Invertebrate","Sessile Invertebrate")))%>%
  ggplot(aes(x = trait_group, 
             fill = org_group)) +
  geom_bar(position = "fill") +
  fill_scale

org_types2 <- pa_pivot2_summary %>%
  select(trait_group) %>%
  mutate(org_group = "All")
org_types2 <- org_types2 %>% rbind(pa_pivot2 %>% select(org_group, trait_group))
org_types2 %>%  
  mutate(org_group = factor(org_group,
                             levels = c("All",
                                        "Arriving",
                                        "Persistent",
                                        "Rare / Other",
                                        "Disappearing")))%>%
  ggplot(aes(x = org_group, 
             fill = trait_group)) +
  geom_bar(position = "fill") 


pa_therm_filtered %>%
  distinct(organism, trait_group) %>%
  count(trait_group) %>%
  mutate(perc = n/sum(n)*100) 


# upload appledore temps --------------------------------------------------

# upload raw temperature values in the study area
temp <- read.csv(here::here(
  "data-processed",
  "appledore-island-env-data",
  "appledore_temps_1982_2023.csv"
))

temp$sample_year[which.max(temp$mean)]
temp %>% arrange(desc(mean))
temp$min[temp$sample_year == 1982]

tempmod_mean <- lm(data = temp,
                   mean ~ sample_year
              )
tempmod_mean_au <- broom::augment(tempmod_mean)

tempmod_max <- lm(data = temp,
                   max ~ sample_year
)
tempmod_max_au <- broom::augment(tempmod_max)

tempmod_min <- lm(data = temp,
                   min ~ sample_year
)
tempmod_min_au <- broom::augment(tempmod_min)
rm(tempmod_mean, tempmod_min, tempmod_max)



# model responses ---------------------------------------------------------
# run anova
mod <- lm(data = pa_pivot2_summary,
          formula = "mean_monthly_mean ~ org_group")
summary(mod)
# significant differences between groups

# run tukey for pairwise comparisons
tukey <- TukeyHSD(aov(mod))

cld <- multcompLetters4(mod, tukey, reversed = T)
cld_df <- as.data.frame.list(cld$org_group) %>%
  mutate(org_group = row.names(.)) %>%
  select(org_group, Letters) %>%
  mutate(org_group = factor(org_group,
                            levels = c("Disappearing","Persistent","Rare / Other","Arriving")))

p <- pa_pivot2_summary %>% 
  left_join(cld_df) %>% 
  group_by(org_group) %>%
  ggplot(aes(x = org_group,
             y = mean_monthly_mean,
             fill = org_group)) +
  # add true temp change
  geom_hline(yintercept = tempmod_mean_au[tempmod_mean_au$sample_year %in% range(temp$sample_year),]$.fitted,
             linewidth = .25,
             linetype = c("solid","dashed")) +  
  geom_text(data = tempmod_mean_au %>% select(.fitted,sample_year) %>% filter(sample_year %in% range(sample_year)),
            inherit.aes=F,
            aes(x = "Arriving",
                y = .fitted,
                label = sample_year),
            hjust = 0,
            vjust = -.15,
            nudge_x = .5,
            size = 3.5) +
  coord_cartesian(xlim = c(1,4.7)) +
  
  geom_rug(data = temp,
           inherit.aes=F,
           linewidth = .1,
           aes(y = mean)) +
  
  geom_violin(alpha = .5) +
  ggbeeswarm::geom_beeswarm(alpha = .4, size = .6, shape = 21, cex = 1.5, show.legend = F) +
  #geom_boxplot(width = .1, outliers = F, show.legend = F) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, color = "black") +
  stat_summary(fun = mean,
               geom = 'point',
               show.legend = F) +
  geom_text(aes(y = 20,
                label = Letters),
            stat = "unique",
            fontface = "bold",
            size = 5) +
  scale_fill_manual(values = c(Arriving = "magenta3",
                               Disappearing = "darkturquoise",
                               Persistent = "grey70",
                               `Rare / Other` = "grey95"),
                    labels = ~stringr::str_wrap(paste(.,"Species"),12)) +
  theme(panel.grid = element_blank()) +
  labs(#title = "Thermal Affinity of Organismal Groups",
       x = NULL,
       y = "Mean Thermal Affinity",
       fill = "Organism Group") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(#legend.position = "none",
        axis.text.x = element_text(face = "bold",
                                   angle = 45, 
                                   hjust = 1,
                                   vjust = 1,
                                   #size = 16,
                                   #lineheight = .8)
        ),
        legend.key.spacing.y   = unit(6,"pt"),
        legend.background = element_rect(color = "black")) +
  guides(fill = guide_legend(nrow = 2))


p + theme(legend.position = "none" ) + ggview::canvas(4,4)

# count species who the temp line crosses
meantemp_start <- tempmod_mean_au[tempmod_mean_au$sample_year == 1982,]$.fitted
meantemp_end <- tempmod_mean_au[tempmod_mean_au$sample_year == 2023,]$.fitted



pa_pivot2_summary %>% 
  summarize(n_below_start = sum(mean_monthly_mean < meantemp_start),
            n_below_end = sum(mean_monthly_mean < meantemp_end),
            n_total = n()) %>%
  mutate(n_below_start_prop = n_below_start/n_total*100,
         n_below_end_prop = n_below_end/n_total*100) %>%
  glimpse()


pa_pivot2_summary %>% 
  group_by(org_group) %>%
  summarize(n_below_start = sum(mean_monthly_mean < meantemp_start),
            n_below_end = sum(mean_monthly_mean < meantemp_end),
            n_total = n()) %>%
  mutate(n_below_start_prop = n_below_start/n_total*100,
         n_below_end_prop = n_below_end/n_total*100) %>%
  glimpse()

# repeat with max ---------------------------------------------------------
# run anova
mod_max <- lm(data = pa_pivot2_summary,
          formula = "perc_monthly_max_95 ~ org_group")
summary(mod_max)
# significant differences between groups

# run tukey for pairwise comparisons
tukey_max <- TukeyHSD(aov(mod_max))

cld_max <- multcompLetters4(mod_max, tukey_max, reversed = T)
cld_df_max <- as.data.frame.list(cld_max$org_group) %>%
  mutate(org_group = row.names(.)) %>%
  select(org_group, Letters) %>%
  mutate(org_group = factor(org_group,
                            levels = c("Disappearing","Persistent","Rare / Other","Arriving")))

p_max <- pa_pivot2_summary %>% 
  left_join(cld_df_max) %>% 
  group_by(org_group) %>%
  ggplot(aes(x = org_group,
             y = perc_monthly_max_95,
             fill = org_group)) +
  
  geom_hline(yintercept = tempmod_max_au[tempmod_max_au$sample_year %in% range(temp$sample_year),]$.fitted,
             linewidth = .25,
             linetype = c("solid","dashed")) +  
  geom_text(data = tempmod_max_au %>% select(.fitted,sample_year) %>% filter(sample_year %in% range(sample_year)),
            inherit.aes=F,
            aes(x = "Arriving",
                y = .fitted,
                label = sample_year),
            hjust = 0,
            vjust = -.15,
            nudge_x = .5,
            size = 3.5) +
  coord_cartesian(xlim = c(1,4.7)) +
  
  geom_rug(data = temp,
           inherit.aes=F,
           linewidth = .1,
           aes(y = max)) +
  
  geom_violin(alpha = .5) +
  ggbeeswarm::geom_beeswarm(alpha = .4, size = .6, shape = 21, cex = 1.5, show.legend = F) +
  
  annotate(geom = "text",
           x = 4.5, 
           y = 17,
           label = "Heat\nStress",
           size = 2.75,
           lineheight = .8,
           hjust = 0) +
  annotate(geom = "text",
           x = 4.5, 
           y = 24,
           label = "No\nHeat\nStress",
           size = 2.75,
           lineheight = .8,
           hjust = 0) +
  
  #geom_boxplot(width = .1, outliers = F, show.legend = F) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, color = "black") +
  stat_summary(fun = mean,
               geom = 'point',
               show.legend = F) +
  geom_text(aes(y = 31,
                label = Letters),
            stat = "unique",
            fontface = "bold",
            size = 5) +
  fill_scale +
  theme(panel.grid = element_blank()) +
  labs(#title = "Thermal Affinity of Organismal Groups",
       x = NULL,
       y = "Thermal Maximum") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold",
                                   angle = 45, 
                                   hjust = 1,
                                   vjust = 1,
                                   #size = 16,
                                   #lineheight = .8)
        ))


p_max + ggview::canvas(4,4)

# if we think more about our hypothesis, it should really be that maximum temps drive disappearances.
# lets' try again but change groups to "disappearing" and "not disappearing"
# run anova

maxtemp_start <- tempmod_max_au[tempmod_max_au$sample_year == 1982,]$.fitted
maxtemp_end <- tempmod_max_au[tempmod_max_au$sample_year == 2023,]$.fitted



pa_pivot2_summary %>% 
  summarize(n_below_start = sum(perc_monthly_max_95 < maxtemp_start),
            n_below_end = sum(perc_monthly_max_95 < maxtemp_end),
            n_total = n()) %>%
  mutate(n_below_start_prop = n_below_start/n_total*100,
         n_below_end_prop = n_below_end/n_total*100) %>%
  glimpse()

pa_pivot2_summary %>% 
  group_by(org_group) %>%
  summarize(n_below_start = sum(perc_monthly_max_95 < maxtemp_start),
            n_below_end = sum(perc_monthly_max_95 < maxtemp_end),
            n_total = n()) %>%
  mutate(n_below_start_prop = n_below_start/n_total*100,
         n_below_end_prop = n_below_end/n_total*100) %>%
  glimpse()


# repeat with min ---------------------------------------------------------
# run anova
mod_min <- lm(data = pa_pivot2_summary,
              formula = "perc_monthly_min_05 ~ org_group")
summary(mod_min)
# significant differences between groups

# run tukey for pairwise comparisons
tukey_min <- TukeyHSD(aov(mod_min))

cld_min <- multcompLetters4(mod_min, tukey_min, reversed = T)
cld_df_min <- as.data.frame.list(cld_min$org_group) %>%
  mutate(org_group = row.names(.)) %>%
  select(org_group, Letters) %>%
  mutate(org_group = factor(org_group,
                            levels = c("Disappearing","Persistent","Rare / Other","Arriving")))

p_min <- pa_pivot2_summary %>% 
  left_join(cld_df_min) %>% 
  group_by(org_group) %>%
  ggplot(aes(x = org_group,
             y = perc_monthly_min_05,
             fill = org_group)) +
  
  geom_hline(yintercept = tempmod_min_au[tempmod_min_au$sample_year %in% range(temp$sample_year),]$.fitted,
             linewidth = .25,
             linetype = c("solid","dashed")) +
  
  geom_text(data = tempmod_min_au %>% select(.fitted,sample_year) %>% filter(sample_year %in% range(sample_year)),
            inherit.aes=F,
            aes(x = "Arriving",
                y = .fitted,
                label = sample_year),
            hjust = 0,
            vjust = -.15,
            nudge_x = .5,
            size = 3.5) +
  coord_cartesian(xlim = c(1,4.7)) +
  
  geom_rug(data = temp,
           inherit.aes=F,
           linewidth = .1,
           aes(y = min)) +
  
  annotate(geom = "text",
           x = 4.5, 
           y = 7.5,
           label = "Cold\nStress",
           size = 2.75,
           lineheight = .8,
           hjust = 0) +
  annotate(geom = "text",
           x = 4.5, 
           y = 1,
           label = "No \nCold\nStress",
           size = 2.75,
           lineheight = .8,
           hjust = 0) +
  
  
  geom_violin(alpha = .5) +
  ggbeeswarm::geom_beeswarm(alpha = .4, size = .6, shape = 21, cex = 1.5, show.legend = F) +
  #geom_boxplot(width = .1, outliers = F, show.legend = F) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.15, color = "black") +
  stat_summary(fun = mean,
               geom = 'point',
               show.legend = F) +
  geom_text(aes(y = 9,
                label = Letters),
            stat = "unique",
            fontface = "bold",
            size = 5) +
  fill_scale +
  
  theme(panel.grid = element_blank()) +
  labs(#title = "Thermal Affinity of Organismal Groups",
       x = NULL,
       y = "Thermal Minimum") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold",
                                   angle = 45, 
                                   hjust = 1,
                                   vjust = 1,
                                   #size = 16,
                                   #lineheight = .8)
                                   ))


p_min + ggview::canvas(4,4)

mintemp_start <- tempmod_min_au[tempmod_min_au$sample_year == 1982,]$.fitted
mintemp_end <- tempmod_min_au[tempmod_min_au$sample_year == 2023,]$.fitted



pa_pivot2_summary %>% 
  summarize(n_below_start = sum(perc_monthly_min_05 < mintemp_start),
            n_below_end = sum(perc_monthly_min_05 < mintemp_end),
            n_total = n()) %>%
  mutate(n_below_start_prop = n_below_start/n_total*100,
         n_below_end_prop = n_below_end/n_total*100) %>%
  glimpse()

pa_pivot2_summary %>% 
  group_by(org_group) %>%
  summarize(n_below_start = sum(perc_monthly_min_05 < mintemp_start),
            n_below_end = sum(perc_monthly_min_05 < mintemp_end),
            n_total = n()) %>%
  mutate(n_below_start_prop = n_below_start/n_total*100,
         n_below_end_prop = n_below_end/n_total*100) %>%
  glimpse()




# merge all plots --------------------------------------------------------
library(patchwork)

# get legend
leg <- ggpubr::get_legend(p)

# Convert to a ggplot and print
leg_grob <- ggpubr::as_ggplot(leg)

p_noleg <- p + theme(legend.position = "none")
p_therms <- ( guide_area() / p_max / p / p_min )  + 
  plot_layout(axes = "collect_x",
              widths = c(1,1,1,1),
              heights = c(.5,1,1,1),
              guides = "collect")
p_therms  + ggview::canvas(4,10) 



p_full <- (p_time | p_therms) + 
  plot_layout(widths = c(.35,.65)) +
  plot_annotation(tag_levels = list("a","","b","c","d")) &
  theme(plot.margin = margin(0,0,0,0))

p_full + ggview::canvas(8*.8,
                        12*.8)


ggsave(p_full,
       filename = "outputs/depth_shifts/spp_appear_disappear.png",
       width = 8*.8,
       height = 12*.8,
       unit = "in")



# make alternate version with n plot --------------------------------------


p_therms2 <- (n_plot/ p_max / (p +theme(legend.position = "none"))/ p_min )  + 
  plot_layout(axes = "collect_x",
              widths = c(1,1,1,1),
              heights = c(.4,1,1,1)) 
p_therms2  + ggview::canvas(4,10) 



p_full2 <- (p_time | p_therms2) + 
  plot_layout(widths = c(.35,.65)) +
  plot_annotation(tag_levels = list("a","b","c","d","e")) &
  theme(plot.margin = margin(0,0,0,0))

p_full2 + ggview::canvas(8*.8,
                        12*.8)

ggsave(p_full2,
       filename = "outputs/depth_shifts/spp_appear_disappear_with_nonnative.png",
       width = 8*.8,
       height = 12*.8,
       unit = "in")



pa_pivot2_summary %>% filter(org_group == "Disappearing",
                             perc_monthly_max_95 > 25) %>% glimpse()

pa_therm_filtered %>% filter(organism == "Asterias forbesi")




# make plot of disappearing species vs temperature through time -----------
pa_pivot_lag <- 
  pa_pivot %>%  
  tidyr::pivot_longer(cols = c("yg_1":"yg_4"),
                      names_to = "year_group",
                      values_to = "present") %>% 
  mutate(present = case_when(is.na(present) ~ FALSE,
                             TRUE ~ TRUE)) %>% 
  group_by(organism) %>%
  mutate(change = case_when(
    present == TRUE & lag(present, default = FALSE) == FALSE ~ "appeared",  # Appeared in this decade
    present == FALSE & lag(present, default = TRUE) == TRUE ~ "disappeared", # Disappeared in this decade
    TRUE ~ "persisted"  # No change
  )) %>% 
  arrange(organism, year_group) %>%
  group_by(organism) %>%
  mutate(change = case_when(
    # For decade 1, if present, mark as 'initially present'
    year_group == "yg_1" & present == TRUE ~ "initially present",
    # For decade 1, if present, mark as 'initially absent'
    year_group == "yg_1" & present == FALSE ~ "initially absent",
    # Appeared in a subsequent decade (present, but was absent in the prior decade)
    present == TRUE & lag(present, default = FALSE) == FALSE & year_group != "decade1" ~ "appeared",  
    # Disappeared (absent, but was present in the prior decade)
    present == FALSE & lag(present, default = TRUE) == TRUE ~ "disappeared", 
    # If both time points are false, consider it as 'no change'
    present == FALSE & lag(present, default = TRUE) == FALSE ~ NA_character_, 
    TRUE ~ NA_character_  # In other cases, we leave it blank (NA)
  )) %>%
  ungroup() %>%
    mutate(change = case_when(is.na(change) ~ "no special status",
                              change %in% c("initially present","initially absent") ~ "no special status",
                              TRUE ~ change)) %>%
    filter(present == T | change == "disappeared") %>%
  group_by(organism) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(jitter_width = runif(nrow(.),-.23, .23)) %>%
  #mutate(jitter_width = rnorm(nrow(.),0, .1)) %>%
  tidyr::unnest(data) %>%
  mutate(year_group = factor(year_group,
         levels = c("yg_1","yg_2","yg_3","yg_4")))
  
pa_pivot_lag %>% glimpse()
as.numeric(pa_pivot_lag$year_group)



# make mean temp plot -----------------------------------------------------
tempmod_means  <- tempmod_mean_au %>%
  mutate(year_group = case_when(sample_year %in% y1 ~ "yg_1",
                                sample_year %in% y2 ~ "yg_2",
                                sample_year %in% y3 ~ "yg_3",
                                sample_year %in% y4 ~ "yg_4")) %>%
  mutate(year_group = factor(year_group, levels = c("yg_1","yg_2","yg_3","yg_4")))

p_loop_mean <- 
  pa_pivot_lag %>%
  ggplot(aes(x = as.numeric(year_group),
             y = mean_monthly_mean)) +
  # first add temp ribbons
  geom_rect(data = tempmod_means %>% group_by(year_group) %>%
              summarize(min = min(.fitted),
                        max = max(.fitted)), 
            inherit.aes=F,
            aes(xmin = as.numeric(year_group) - .4,
                xmax = as.numeric(year_group) + .4,
                ymin = min,
                ymax  = max),
            alpha = .15,
            fill= "grey50",
            stat = "unique") +
  geom_point(data = tempmod_means, 
             inherit.aes=F,
             aes(x = as.numeric(year_group) - .4,
                 y = mean),
             alpha = 1,
             stroke = .1,
             shape = "-", 
             size = 5)  +
  scale_shape_manual(breaks = c("appeared","no special status",
                                "disappeared"),
                     labels = c("Appeared this decade",
                                "No Change",
                                "Disappeared this decade"),
                     values = c("disappeared" = 1,
                                "appeared" = 16,
                                "no special status" = 16)) +
  
  scale_color_manual(values = c(Arriving = "magenta3",
                                Disappearing = "darkturquoise",
                                Persistent = "grey70",
                                `Rare / Other` = "grey85")) 

# now add points in a loop
# Loop to add two points for each row, one at a time
for(i in 1:nrow(pa_pivot_lag)) {
  p_loop_mean <- p_loop_mean + 
    
    # add first layer of points, shaped by change category
    geom_point(
      data = pa_pivot_lag[i,],
      size = 1.5, stroke = 1,
      aes(x = as.numeric(year_group) + jitter_width,
          shape = change, color = org_group),
      alpha = .8,
      # position = position_jitter(width = .35, height = 0, seed = 1)
    ) +
    # add second layer of points, which outlines the appeared cells
    geom_point(data = pa_pivot_lag %>% 
                 mutate(mean_monthly_mean = ifelse(change == 'appeared',mean_monthly_mean, NA)) %>%
                 slice(i),
               aes(x = as.numeric(year_group) + jitter_width),
               size = 1.5, shape = 1, stroke = .5, color = "black",
               # position = position_jitter(width = .35, height = 0, seed=1)
    ) +
    # add third layer of points, exing out lost points:
    geom_point(data = pa_pivot_lag %>% 
                 mutate(mean_monthly_mean = ifelse(change == 'disappeared',mean_monthly_mean, NA)) %>%
                 slice(i),
               aes(x = as.numeric(year_group) + jitter_width),
               size = 1.5, shape = 4, stroke = .5, color = "black",
               # position = position_jitter(width = .35, height = 0, seed=1)
    ) 
  print(i)
}

p_loop_mean <- p_loop_mean +
  guides(shape = guide_legend(
    override.aes = list(
      shape = c(21, 16, 1),
      color = c("black", "grey80","grey80"),
      fill = c("grey80", NA, NA)
    ))
    ) +
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c("Decade 1","Decade 2","Decade 3","Decade 4")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  theme(panel.grid.major.x = element_blank()) +
  labs(x = NULL,
       y = "Mean Thermal Affinity (°C)",
       shape = "Decadal Change",
       color = "Long-Term Change") +
  scale_y_continuous(labels = ~paste0(.,"°"))

#p_loop_mean + ggview::canvas(7,4)

# make max temp plot -----------------------------------------------------
tempmod_maxes  <- tempmod_max_au %>%
  mutate(year_group = case_when(sample_year %in% y1 ~ "yg_1",
                                sample_year %in% y2 ~ "yg_2",
                                sample_year %in% y3 ~ "yg_3",
                                sample_year %in% y4 ~ "yg_4")) %>%
  mutate(year_group = factor(year_group, levels = c("yg_1","yg_2","yg_3","yg_4")))


p_loop_max <- 
  pa_pivot_lag %>%
  ggplot(aes(x = as.numeric(year_group),
             y = perc_monthly_max_95)) +
  # first add temp ribbons
  geom_rect(data = tempmod_maxes %>% group_by(year_group) %>%
              summarize(min = min(.fitted),
                        max = max(.fitted)), 
            inherit.aes=F,
            aes(xmin = as.numeric(year_group) - .4,
                xmax = as.numeric(year_group) + .4,
                ymin = min,
                ymax  = max),
            alpha = .15,
            fill= "grey50",
            stat = "unique") +
  geom_point(data = tempmod_maxes, 
             inherit.aes=F,
             aes(x = as.numeric(year_group) - .4,
                 y = max),
             alpha = 1,
             stroke = .1,
             shape = "-", 
             size = 5)  +
  scale_shape_manual(breaks = c("appeared","no special status",
                                "disappeared"),
                     labels = c("Appeared this decade",
                                "No Change",
                                "Disappeared this decade"),
                     values = c("disappeared" = 1,
                                "appeared" = 16,
                                "no special status" = 16)) +
  
  scale_color_manual(values = c(Arriving = "magenta3",
                                Disappearing = "darkturquoise",
                                Persistent = "grey70",
                                `Rare / Other` = "grey85")) 

# now add points in a loop
# Loop to add two points for each row, one at a time
for(i in 1:nrow(pa_pivot_lag)) {
  p_loop_max <- p_loop_max + 
    
    # add first layer of points, shaped by change category
    geom_point(
      data = pa_pivot_lag[i,],
      size = 1.5, stroke = 1,
      aes(x = as.numeric(year_group) + jitter_width,
          shape = change, color = org_group),
      alpha = .8,
      # position = position_jitter(width = .35, height = 0, seed = 1)
    ) +
    # add second layer of points, which outlines the appeared cells
    geom_point(data = pa_pivot_lag %>% 
                 mutate(perc_monthly_max_95 = ifelse(change == 'appeared',perc_monthly_max_95, NA)) %>%
                 slice(i),
               aes(x = as.numeric(year_group) + jitter_width),
               size = 1.5, shape = 1, stroke = .5, color = "black",
               # position = position_jitter(width = .35, height = 0, seed=1)
    ) +
    # add third layer of points, exing out lost points:
    geom_point(data = pa_pivot_lag %>% 
                 mutate(perc_monthly_max_95 = ifelse(change == 'disappeared',perc_monthly_max_95, NA)) %>%
                 slice(i),
               aes(x = as.numeric(year_group) + jitter_width),
               size = 1.5, shape = 4, stroke = .5, color = "black",
               # position = position_jitter(width = .35, height = 0, seed=1)
    ) 
  print(i)
}

p_loop_max <- p_loop_max +
  guides(shape = guide_legend(
    override.aes = list(
      shape = c(21, 16, 1),
      color = c("black", "grey80","grey80"),
      fill = c("grey80", NA, NA)
    ))
  ) +
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c("Decade 1","Decade 2","Decade 3","Decade 4")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  theme(panel.grid.major.x = element_blank()) +
  labs(x = NULL,
       y = "Thermal Maximum (°C)",
       shape = "Decadal Change",
       color = "Long-Term Change") +
  scale_y_continuous(labels = ~paste0(.,"°"))

# make min temp plot -----------------------------------------------------
tempmod_mins  <- tempmod_min_au %>%
  mutate(year_group = case_when(sample_year %in% y1 ~ "yg_1",
                                sample_year %in% y2 ~ "yg_2",
                                sample_year %in% y3 ~ "yg_3",
                                sample_year %in% y4 ~ "yg_4")) %>%
  mutate(year_group = factor(year_group, levels = c("yg_1","yg_2","yg_3","yg_4")))


p_loop_min <- 
  pa_pivot_lag %>%
  ggplot(aes(x = as.numeric(year_group),
             y = perc_monthly_min_05)) +
  # first add temp ribbons
  geom_rect(data = tempmod_mins %>% group_by(year_group) %>%
              summarize(min = min(.fitted),
                        max = max(.fitted)), 
            inherit.aes=F,
            aes(xmin = as.numeric(year_group) - .4,
                xmax = as.numeric(year_group) + .4,
                ymin = min,
                ymax  = max),
            alpha = .15,
            fill= "grey50",
            stat = "unique") +
  geom_point(data = tempmod_mins, 
             inherit.aes=F,
             aes(x = as.numeric(year_group) - .4,
                 y = min),
             alpha = 1,
             stroke = .1,
             shape = "-", 
             size = 5)  +
  scale_shape_manual(breaks = c("appeared","no special status",
                                "disappeared"),
                     labels = c("Appeared this decade",
                                "No Change",
                                "Disappeared this decade"),
                     values = c("disappeared" = 1,
                                "appeared" = 16,
                                "no special status" = 16)) +
  
  scale_color_manual(values = c(Arriving = "magenta3",
                                Disappearing = "darkturquoise",
                                Persistent = "grey70",
                                `Rare / Other` = "grey85")) 

# now add points in a loop
# Loop to add two points for each row, one at a time
for(i in 1:nrow(pa_pivot_lag)) {
  p_loop_min <- p_loop_min + 
    
    # add first layer of points, shaped by change category
    geom_point(
      data = pa_pivot_lag[i,],
      size = 1.5, stroke = 1,
      aes(x = as.numeric(year_group) + jitter_width,
          shape = change, color = org_group),
      alpha = .8,
      # position = position_jitter(width = .35, height = 0, seed = 1)
    ) +
    # add second layer of points, which outlines the appeared cells
    geom_point(data = pa_pivot_lag %>% 
                 mutate(perc_monthly_min_05 = ifelse(change == 'appeared',perc_monthly_min_05, NA)) %>%
                 slice(i),
               aes(x = as.numeric(year_group) + jitter_width),
               size = 1.5, shape = 1, stroke = .5, color = "black",
               # position = position_jitter(width = .35, height = 0, seed=1)
    ) +
    # add third layer of points, exing out lost points:
    geom_point(data = pa_pivot_lag %>% 
                 mutate(perc_monthly_min_05 = ifelse(change == 'disappeared',perc_monthly_min_05, NA)) %>%
                 slice(i),
               aes(x = as.numeric(year_group) + jitter_width),
               size = 1.5, shape = 4, stroke = .5, color = "black",
               # position = position_jitter(width = .35, height = 0, seed=1)
    ) 
  print(i)
}

p_loop_min <- p_loop_min +
  guides(shape = guide_legend(
    override.aes = list(
      shape = c(21, 16, 1),
      color = c("black", "grey80","grey80"),
      fill = c("grey80", NA, NA)
    ))
  ) +
  scale_x_continuous(breaks = c(1,2,3,4),
                     labels = c("Decade 1","Decade 2","Decade 3","Decade 4")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +
  theme(panel.grid.major.x = element_blank()) +
  labs(x = NULL,
       y = "Thermal Minimum (°C)",
       shape = "Decadal Change",
       color = "Long-Term Change") +
  scale_y_continuous(labels = ~paste0(.,"°"))



# merge -------------------------------------------------------------------

p_loop_full <- egg::ggarrange(p_loop_max + theme(legend.position = "none"),
                              p_loop_mean + theme(legend.position = "none"),
                              p_loop_min + theme(legend.position = "bottom",
                                                 legend.box = "vertical"),
                              nrow = 3)

ggsave(p_loop_full, 
       filename = 'outputs/depth_shifts/spp_appear_disappear_decadal.png',
       width = 6, 
       height = 8)


p_loop_full + ggview::canvas(4,12)
pa_pivot_lag %>%
ggplot(aes(x = as.numeric(year_group) ,
             y = mean_monthly_mean)) +
  
  # first add temp ribbons
  geom_rect(data = tempmod_means %>% group_by(year_group) %>%
              summarize(min = min(.fitted),
                        max = max(.fitted)), 
            inherit.aes=F,
            aes(xmin = as.numeric(year_group) - .4,
                xmax = as.numeric(year_group) + .4,
                ymin = min,
                ymax  = max),
            alpha = .15,
            fill= "grey50",
            stat = "unique") +
  geom_point(data = tempmod_means, 
            inherit.aes=F,
            aes(x = as.numeric(year_group) - .4,
                y = mean),
            alpha = 1,
            stroke = .1,
            shape = "-", 
            size = 5) +
  
  geom_point(
    aes(shape = change,
        color = case_when(change == "appeared" ~ "black",))
  )


  geom_point(
    data = . %>% arrange(change),
    size = 2.5, stroke = 1,
    aes(x = as.numeric(year_group) + jitter_width),
    alpha = .8,
   # position = position_jitter(width = .35, height = 0, seed = 1)
   ) +
  geom_point(data = . %>% mutate(mean_monthly_mean = ifelse(change == 'appeared',mean_monthly_mean, NA)),
             aes(x = as.numeric(year_group) + jitter_width),
             size = 2.5, shape = 1, stroke = .5, color = "black",
            # position = position_jitter(width = .35, height = 0, seed=1)
            ) +
 # geom_point(data = . %>% mutate(mean_monthly_mean = ifelse(change == 'disappeared',mean_monthly_mean, NA)),
 #            aes(x = as.numeric(year_group) + jitter_width),
 #            size = 3, shape = 4, stroke = .5, color = "black",
 #            # position = position_jitter(width = .35, height = 0, seed=1)
 # ) +
  scale_x_continuous(breaks = c(1:4),
                     labels = c("Year Group 1",
                                "Year Group 2",
                                "Year Group 3",
                                "Year Group 4")) +
  scale_shape_manual(values = c("appeared" = 16,
                                "no special status" = 16,
                                "disappeared" = 1)) +
  scale_color_manual(values = c(Arriving = "magenta3",
                               Disappearing = "darkturquoise",
                               Persistent = "grey70",
                               `Rare / Other` = "grey85")) +
  guides(shape = guide_legend(override.aes = list(shape = c(21, 1, 19),
                                                  color = c("black","grey80","grey80"),
                                                  fill = c("grey80")))) +
  ggview::canvas(7,4)



  # Original data
  df <- data.frame(
    species = c('a', 'b', 'c'),
    value = c(1, 5, 3),
    decade1 = c(TRUE, FALSE, TRUE),
    decade2 = c(FALSE, TRUE, FALSE),
    decade3 = c(FALSE, TRUE, TRUE),
    decade4 = c(FALSE, TRUE, TRUE)
  )
  
  # Pivoting data to a longer format
  df_long <- df %>%
    tidyr::pivot_longer(cols = starts_with("decade"), 
                 names_to = "decade", 
                 values_to = "presence") %>%
    arrange(species, decade)
  
  # Create the 'change' column to indicate appearance or disappearance
  df_long2 <- df_long %>%
    group_by(species) %>%
    mutate(change = case_when(
      # For decade 1, if present, mark as 'initially present'
      decade == "decade1" & presence == TRUE ~ "initially present",
      # For decade 1, if present, mark as 'initially absent'
      decade == "decade1" & presence == FALSE ~ "initially absent",
      # Appeared in a subsequent decade (present, but was absent in the prior decade)
      presence == TRUE & lag(presence, default = FALSE) == FALSE & decade != "decade1" ~ "appeared",  
      # Disappeared (absent, but was present in the prior decade)
      presence == FALSE & lag(presence, default = TRUE) == TRUE ~ "disappeared", 
      # If both time points are false, consider it as 'no change'
      presence == FALSE & lag(presence, default = TRUE) == FALSE ~ NA_character_, 
      TRUE ~ NA_character_  # In other cases, we leave it blank (NA)
    )) %>%
    ungroup()
  
  # View the transformed data
  print(df_long2)
  
  df_long2 %>%
    mutate(change = case_when(is.na(change) ~ "no special status",
                              change %in% c("initially present","initially absent") ~ "no special status",
                              TRUE ~ change)) %>%
    filter(presence == T | change == "disappeared") %>%
    ggplot(aes(x = decade, y = value,
               shape = change, color = species)) +
    geom_point(#color = "lightblue", 
               size = 6, stroke = 2) +
    geom_point(data = . %>% filter(change == "appeared"),
               size = 6, shape = 1, stroke = 2, color = "black") +
    scale_shape_manual(values = c("disappeared" = 13,
                                  "no special status" = 19,
                                  "appeared" = 16)) +
    guides(shape = guide_legend(override.aes = list(
      color = c("black", "grey80","grey80"),
      shape = c(21, 13, 19),
      fill = c("grey80",NA,NA)
    )))
  
  df_long2 %>%
    mutate(change = case_when(is.na(change) ~ "no special status",
                              change %in% c("initially present","initially absent") ~ "no special status",
                              TRUE ~ change)) %>%
    filter(presence == T | change == "disappeared") %>%
  ggplot( aes(x = decade, y = value, color = species)) +
    geom_point(aes(
      fill = case_when(
        change == "appeared" ~ species,
        change == "no special status" ~ NA_character_,
        change == "disappeared" ~ NA_character_
      ),
      shape = change,
      color = case_when(
        change == "appeared" ~ change,
        change == "no special status" ~ species,
        change == "disappeared" ~ species
      )
    ), size = 5, stroke = 2) +
    scale_shape_manual(values = c("disappeared" = 1,
                                  "appeared" = 21,
                                  "no special status" = 19)) +
    scale_color_manual(values = c("appeared" = "black",
                                  "a" = "tomato",
                                  "b" = "lightblue",
                                  "c" = "lightgreen")) +
    scale_fill_manual(values = c(#"no fill" = "transparent",
                                 b = "lightblue",
                                 c = "lightgreen"),
                      na.translate = FALSE)
    scale_fill_identity() +  # Ensure fill is based on species color
    scale_color_identity() + # Ensure color is based on species color
    theme_minimal()

  # add asterisk to organisms that are invasive
  mutate(organism_label = paste0("<i>",organism,"</i>")) %>%
  mutate(organism_label = case_when(invasive_gom == "yes" ~ paste0(organism_label,"<b> *</b>"),
                                    TRUE ~ organism_label)) %>%
  mutate(organism_label = forcats::fct_reorder(organism_label, mean_monthly_mean)) %>%
  mutate(org_group = factor(org_group,
                            levels = c("Disappearing","Persistent","Rare / Other","Arriving")))



pa_pivot2 %>%
  filter(year_group == "yg_2") %>%
  mutate(status = "Present") %>%
  rbind(pa_pivot2 %>%
          filter(organism %in% c(pa_pivot$organism[is.na(pa_pivot$yg_2) & !is.na(pa_pivot$yg_1)])) %>%
          mutate(status = "Disappeared")) %>%
  ggplot(aes(x = 1, y = mean_monthly_mean, shape = status)) +
  geom_point(position = position_jitter(width = 1)) +
  scale_shape_manual(values = c(1,))
pa_pivot %>% glimpse()
  filter(year)


pa_pivot2 %>%
  ggplot(aes(x = year_group, y = mean_monthly_mean)) +
  geom_point(position = position_jitter(width = .2)) +
  facet_grid(~year_group, scales = "free_x")

# test nonnative vs native sti --------------------------------------------
pa_pivot2_summary %>% glimpse()
# run anova
mod_nat <- lm(data = pa_pivot2_summary,
          formula = "mean_monthly_mean ~ invasive_gom")
summary(mod_nat)
# significant differences between groups

# run tukey for pairwise comparisons
tukey_nat <- TukeyHSD(aov(mod_nat))

cld_nat <- multcompLetters4(mod_nat, tukey_nat, reversed = T)
cld_df_nat <- as.data.frame.list(cld_nat$invasive_gom) %>%
  mutate(invasive_gom = row.names(.)) %>%
  select(invasive_gom, Letters) %>%
  mutate(invasive_gom = factor(invasive_gom,
                            levels = c("no","yes")))

p_nat <- pa_pivot2_summary %>% 
  left_join(cld_df_nat) %>%  
  group_by(invasive_gom) %>%
  ggplot(aes(x = invasive_gom,
             y = mean_monthly_mean,
             fill = invasive_gom)) +
  geom_violin(alpha = .5) +
  geom_boxplot(width = .1, outliers = F, show.legend = F) +
  geom_text(aes(y = 20,
                label = Letters),
            stat = "unique",
            fontface = "bold",
            size = 5) +
  theme(panel.grid = element_blank()) +
  labs(#title = "Thermal Affinity of Organismal Groups",
    x = NULL,
    y = "Thermal Affinity (Mean)",
    fill = "Native Group") +
  scale_y_continuous(labels = ~paste0(.,"°")) +
  scale_x_discrete(labels = c("Native","Non-Native"))+
  scale_fill_manual(values = c("grey90","grey40")) +
  theme(legend.position = "none",
   # axis.text.x = element_text(face = "bold",
   #                            angle = 45, 
   #                            hjust = 1,
   #                            vjust = 1,
   #                            #size = 16,
   #                            #lineheight = .8)
   # ),
    legend.key.spacing.y   = unit(6,"pt"),
    legend.background = element_rect(color = "black"))

p_nat



