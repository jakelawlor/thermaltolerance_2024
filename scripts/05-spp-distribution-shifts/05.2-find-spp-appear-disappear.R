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
  group_by(organism, mean_monthly_mean, mean_monthly_min, mean_monthly_max, trait_group, invasive_gom) %>% 
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
    (is.na(yg_1)) & (yg_4 == TRUE) ~ "Added",
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
                            levels = c("Disappearing","Persistent","Rare / Other","Added")))


# make fill scale that will apply to all plots
fill_scale <- {
  scale_fill_manual(values = c(Added = "magenta3",
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
  distinct(organism, across(contains("mean_")), org_group, trait_group, invasive_gom)


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
                                        "Added",
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
                            levels = c("Disappearing","Persistent","Rare / Other","Added")))

p <- pa_pivot2_summary %>% 
  left_join(cld_df) %>% 
  group_by(org_group) %>%
  ggplot(aes(x = org_group,
             y = mean_monthly_mean,
             fill = org_group)) +
  geom_violin(alpha = .5) +
  geom_boxplot(width = .1, outliers = F, show.legend = F) +
  geom_text(aes(y = 20,
                label = Letters),
            stat = "unique",
            fontface = "bold",
            size = 5) +
  scale_fill_manual(values = c(Added = "magenta3",
                               Disappearing = "darkturquoise",
                               Persistent = "grey70",
                               `Rare / Other` = "grey95"),
                    labels = ~stringr::str_wrap(paste(.,"Species"),12)) +
  theme(panel.grid = element_blank()) +
  labs(#title = "Thermal Affinity of Organismal Groups",
       x = NULL,
       y = "Thermal Affinity (Mean)",
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


p + ggview::canvas(4,4)



# repeat with max ---------------------------------------------------------
# run anova
mod_max <- lm(data = pa_pivot2_summary,
          formula = "mean_monthly_max ~ org_group")
summary(mod_max)
# significant differences between groups

# run tukey for pairwise comparisons
tukey_max <- TukeyHSD(aov(mod_max))

cld_max <- multcompLetters4(mod_max, tukey_max, reversed = T)
cld_df_max <- as.data.frame.list(cld_max$org_group) %>%
  mutate(org_group = row.names(.)) %>%
  select(org_group, Letters) %>%
  mutate(org_group = factor(org_group,
                            levels = c("Disappearing","Persistent","Rare / Other","Added")))

p_max <- pa_pivot2_summary %>% 
  left_join(cld_df_max) %>% 
  group_by(org_group) %>%
  ggplot(aes(x = org_group,
             y = mean_monthly_max,
             fill = org_group)) +
  geom_violin(alpha = .5) +
  geom_boxplot(width = .1, outliers = F, show.legend = F) +
  geom_text(aes(y = 26.5,
                label = Letters),
            stat = "unique",
            fontface = "bold",
            size = 5) +
  fill_scale +
  theme(panel.grid = element_blank()) +
  labs(#title = "Thermal Affinity of Organismal Groups",
       x = NULL,
       y = "Thermal Affinity (Max)") +
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


# repeat with min ---------------------------------------------------------
# run anova
mod_min <- lm(data = pa_pivot2_summary,
              formula = "mean_monthly_min ~ org_group")
summary(mod_min)
# significant differences between groups

# run tukey for pairwise comparisons
tukey_min <- TukeyHSD(aov(mod_min))

cld_min <- multcompLetters4(mod_min, tukey_min, reversed = T)
cld_df_min <- as.data.frame.list(cld_min$org_group) %>%
  mutate(org_group = row.names(.)) %>%
  select(org_group, Letters) %>%
  mutate(org_group = factor(org_group,
                            levels = c("Disappearing","Persistent","Rare / Other","Added")))

p_min <- pa_pivot2_summary %>% 
  left_join(cld_df_min) %>% 
  group_by(org_group) %>%
  ggplot(aes(x = org_group,
             y = mean_monthly_min,
             fill = org_group)) +
  geom_violin(alpha = .5) +
  geom_boxplot(width = .1, outliers = F, show.legend = F) +
  geom_text(aes(y = 18,
                label = Letters),
            stat = "unique",
            fontface = "bold",
            size = 5) +
  fill_scale +
  
  theme(panel.grid = element_blank()) +
  labs(#title = "Thermal Affinity of Organismal Groups",
       x = NULL,
       y = "Thermal Affinity (Min)") +
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
  plot_annotation(tag_levels = list("a","","b","c","d"))

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
  plot_annotation(tag_levels = list("a","b","c","d","e"))

p_full2 + ggview::canvas(8*.8,
                        12*.8)

ggsave(p_full2,
       filename = "outputs/depth_shifts/spp_appear_disappear_with_nonnative.png",
       width = 8*.8,
       height = 12*.8,
       unit = "in")







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



