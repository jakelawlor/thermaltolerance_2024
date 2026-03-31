# Assess intertidal height per species as a function of year
# ungrouped version: meaning no 5-year bins, as previously done. 

# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)

# data --------------------------------------------------------------------
# upload PA with thermal tolerances from highly sampled dataset
pa_therm_filtered <- readr::read_csv(here::here("data-processed",
                                                "appledore-survey-data",
                                                "pa-with-therm",
                                                "pa-with-therm-highly-sampled.csv"))



# clean a few columns -----------------------------------------------
pa_hs <- pa_therm_filtered%>%
  tidyr::separate(organism, into = c("Genus","Species")) %>%
  mutate(G = paste0(stringr::str_sub(Genus,1,1),".")) %>%
  mutate(organism = paste(G, Species)) %>%
  select(-G, -Genus, -Species) %>%
  mutate(tidalheight = (13.5-level)*.3048)  
rm(pa_therm_filtered)

pa_hs %>% glimpse()


# find height quantiles
hs_quants <- pa_hs %>%
  
  # count total occurrences per species
  group_by(organism) %>%
  mutate(total_count = n()) %>%
  
  group_by(organism, year, total_count) %>%
  
  # find distributions of height occurrences 
  summarize(heightmax = max(tidalheight),
            height95 = quantile(tidalheight,.95),
            heightmed = median(tidalheight),
            heightmean = mean(tidalheight),
            height05 = quantile(tidalheight,.05),
            heightmin = min(tidalheight)) %>%
  
  # filter to organisms present in 10+ years
  group_by(organism) %>%
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years >= 10) %>%
  ungroup() %>%
  select(-n_years) #%>%
  
# don't need if filtering to 10 years:
#  # filter to first obs and last obs at least 10 years apart
#  group_by(organism) %>%
#  mutate(first_year = min(year),
#         last_year = max(year)) %>%
#  mutate(year_breadth = last_year - first_year) %>%
#  filter(year_breadth >= 10) %>%
#  ungroup() %>%
#  select(-first_year, last_year, year_breadth) %>%
#  filter(total_count >= 10)

hs_quants_split <- hs_quants %>% split(f = .$organism)
#rm(hs_quants)

# find exclusions ---------------------------------------------------------
# assessing qwuantile shifts won't work if the species is found all the
# way to the bounds of the domain in most years. 

# view all occurrences of species
# if LE or TE is right at domain edges, we won't be ableto detect a shift
hs_quants_split[[1]] %>%
  ggplot(aes(x = year, y = heightmax)) +
  geom_point() +
  geom_point(aes(y = height95)) +
  geom_point(aes(y = heightmed)) +
  geom_point(aes(y = heightmean)) +
  geom_point(aes(y = height05)) +
  geom_point(aes(y = heightmin)) +
  geom_smooth(method = "lm",se=F, ) +
  geom_smooth(method = "lm",se=F, aes(y = height95)) +
  geom_smooth(method = "lm",se=F, aes(y = heightmed)) +
  geom_smooth(method = "lm",se=F, aes(y = heightmean)) +
  geom_smooth(method = "lm",se=F, aes(y = height05)) +
  geom_smooth(method = "lm",se=F, aes(y = heightmin)) +
  geom_hline(yintercept = unique(pa_hs$tidalheight)[c(1,9)],
             linetype = "dashed",
             color = "red")

# find domain
domain <-unique(pa_hs$tidalheight)[c(1,9)]

# identify species with 50% or more of their highest occurrences at domain edge
max_to_exclude <- purrr::map(
  .x = hs_quants_split,
  .f = ~.x %>%  summarise(n = n(),
                          at_domain = sum(heightmax == domain[1])) %>%
    mutate(at_domain_prop = at_domain/n) %>%
    mutate(exclude_max = at_domain_prop > .5) %>% 
    pull(exclude_max) %>% as_tibble
) %>% bind_rows(.id = "organism") %>%
  rename(exclude_max  = value)
table(max_to_exclude$exclude_max)

min_to_exclude <- purrr::map(
  .x = hs_quants_split,
  .f = ~.x %>%  summarise(n = n(),
                          at_domain = sum(heightmin == domain[2])) %>%
    mutate(at_domain_prop = at_domain/n) %>%
    mutate(exclude_min = at_domain_prop > .5) %>% 
    pull(exclude_min) %>% as_tibble
) %>% bind_rows(.id = "organism") %>%
  rename(exclude_min  = value)
table(min_to_exclude$exclude_min)

exclusions <- max_to_exclude %>% left_join(min_to_exclude)
rm(min_to_exclude, max_to_exclude)

# add mid to exclude- when the species' max and min are both at domain bounds
exclusions <- exclusions %>%
  mutate(exclude_mid = case_when(exclude_max == T & exclude_min == T ~ T,
                                 TRUE ~ F))

table(exclusions$exclude_mid)

# model observation heights over time -------------------------------------
model <- function(df, param){
  
  mods <- purrr::map(.x = df, 
                     .f = ~lm(data = .x, 
                            formula = get(param) ~ year))
  
  mods_p <- purrr::map(.x = mods,
                       .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
  ) %>%
    bind_rows(.id = "organism") %>%
    rename(p_value = value)
  
  
  mods_trend <- purrr::map(
    .x = mods,
    .f = ~ data.frame(trend = summary(.x)$coefficients["year","Estimate"],
                      sign = ifelse(sign(summary(.x)$coefficients["year","Estimate"]) == 1,
                                          "(+)",
                                          "(-)"
                      )) )%>%
    bind_rows(.id = "organism") 
  
  mods_single_val <- purrr::map(
    .x = df,
    .f = ~ifelse(length(unique(.x[[param]])) == 1,"single value","no issue") %>% as_tibble()
  ) %>%
    bind_rows(.id = "organism") %>%
    rename(single_value = value)
  
  # add excludions
  if(param %in% c("heightmean","heightmed")){
    mod_exclusion <- exclusions %>%
      select(organism, exclude_mid) %>%
      rename(exclude = exclude_mid)
  } else if(param %in% c("heightmin","height05")){
    mod_exclusion <- exclusions %>%
      select(organism, exclude_min) %>%
      rename(exclude = exclude_min)
  } else if(param %in% c("heightmax","height95")){
    mod_exclusion <- exclusions %>%
      select(organism, exclude_max) %>%
      rename(exclude = exclude_max)
  }
  
  # merge all
  out <- mods_p %>%
    left_join(mods_trend,
              by = join_by(organism)) %>%
    left_join(mods_single_val,
              by = join_by(organism)) %>%
    left_join(mod_exclusion, 
              by = join_by(organism)) %>%
    mutate(param = param) %>%
    relocate(param, .after = organism) %>%
    mutate(sign = case_when(single_value == "single value" ~ "(0)",
                            TRUE ~ sign))
  
  return(out)

}


#mods_min <- model(hs_unc_sum_split, "heightmin")
mods_05 <- model(hs_quants_split, "height05")
mods_mean <- model(hs_quants_split, "heightmean")
mods_med <- model(hs_quants_split, "heightmed")
mods_95 <- model(hs_quants_split, "height95")
#mods_max <- model(hs_unc_sum_split, "heightmax")

# merge all ------------------------------------------------------

mods_all <- rbind(#mods_min, 
                  mods_05,
                  mods_mean,
                  mods_med,
                  mods_95#,
                  #mods_max
                  )


mods_all %>% count(exclude)
mods_all %>% filter(single_value == "single value") %>% glimpse()

mods_all_annotate <- mods_all %>%
  mutate(sign = case_when(p_value < 0.05 ~ paste0(sign,"*"),
                          TRUE ~ sign))
mods_all_annotate %>% 
  filter(exclude == F) %>%
  count(param)


# plot --------------------------------------------------------------------

p1 <- mods_all_annotate %>%
 # filter(exclude == F) %>% 
  group_by(param) %>%
  count(sign) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
    mutate(text_color = if_else(sign %in% c("(-)*"), "white", "black")) %>%

  mutate(sign = factor(sign, 
                       levels = c("(+)*","(+)","(0)","(-)","(-)*"))) %>%
  mutate(param = factor(param,
                        levels = c("heightmin","height05",
                                   "heightmean","heightmed",
                                   "height95","heightmax"))) %>%
  filter(param %in% c("height95",'heightmed',"heightmean","height05")) %>%

  ggplot(aes(y = param, 
             x = prop,
             fill = sign)) +
  geom_col() +
  geom_vline(xintercept = .5, linetype = "dashed", linewidth = .25) +
  geom_vline(xintercept = c(0,1), linetype = "solid",
             linewidth = .25) +
  scale_fill_manual(breaks = c("(-)*","(-)","(0)","(+)","(+)*"),
                    values = c(c("darkmagenta",colorspace::lighten("magenta",.4),"grey93",
                                 colorspace::lighten("turquoise",.4), "darkturquoise")),
                    labels = c("Significant\nShift Down","Trend Down","No Change","Trend Up","Significant\nShift Up                   ")) +
  geom_text(aes(x = 1.03, 
                label = paste0("n = ",total)),
            stat = "unique",
            size = 3.5,
            hjust = 0) +
  theme_minimal() +
  geom_text(aes(label = n,
                color = text_color),
            position = position_stack(vjust = 0.5)) +
  scale_color_identity() +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  coord_cartesian(expand = F,
                  xlim = c(0, 1.2),
                  clip = "off") +
  scale_y_discrete(labels = c("5%\nOccurrence\nHeight",
                              "Mean\nOccurrence\nHeight",
                              "Median\nOccurrence\nHeight",
                              "95%\nOccurrence\nHeight")) +
  labs(x = "Proportion of Species",
       y = NULL,
       fill = "Shift Trend",
       title = "All Species") +
 # guides(fill = guide_legend(theme = theme(legend.title.position = "top"))) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = scales::percent_format())

p2 <- mods_all_annotate %>%
  filter(exclude == F) %>% 
  group_by(param) %>%
  count(sign) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  mutate(text_color = if_else(sign %in% c("(-)*"), "white", "black")) %>%
  
  mutate(sign = factor(sign, 
                       levels = c("(+)*","(+)","(0)","(-)","(-)*"))) %>%
  mutate(param = factor(param,
                        levels = c("heightmin","height05",
                                   "heightmean","heightmed",
                                   "height95","heightmax"))) %>%
  filter(param %in% c("height95",'heightmed',"heightmean","height05")) %>%
  
  ggplot(aes(y = param, 
             x = prop,
             fill = sign)) +
  geom_col() +
  geom_vline(xintercept = .5, linetype = "dashed", linewidth = .25) +
  geom_vline(xintercept = c(0,1), linetype = "solid",
             linewidth = .25) +
  scale_fill_manual(breaks = c("(-)*","(-)","(0)","(+)","(+)*"),
                    values = c(c("darkmagenta",colorspace::lighten("magenta",.4),"grey93",
                                 colorspace::lighten("turquoise",.4), "darkturquoise")),
                    labels = c("Significant\nShift Down","Trend Down","No Change","Trend Up","Significant\nShift Up                   ")) +
  geom_text(aes(x = 1.03, 
                label = paste0("n = ",total)),
            stat = "unique",
            size = 3.5,
            hjust = 0) +
  theme_minimal() +
  geom_text(aes(label = n,
                color = text_color),
            position = position_stack(vjust = 0.5)) +
  scale_color_identity() +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  coord_cartesian(expand = F,
                  xlim = c(0, 1.2),
                  clip = "off") +
  scale_y_discrete(labels = c("5%\nOccurrence\nHeight",
                              "Mean\nOccurrence\nHeight",
                              "Median\nOccurrence\nHeight",
                              "95%\nOccurrence\nHeight")) +
  labs(x = "Proportion of Species",
       y = NULL,
       fill = "Shift Trend",
       title = "Species Inside Domain") +
  #guides(fill = guide_legend(theme = theme(legend.title.position = "top"))) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = scales::percent_format())
p2
 mods_all_annotate %>% 
  group_by(param) %>% 
  filter(param %in% c("height95","heightmed","heightmean","height05")) %>%
  count(sign)


full_p <- ((p1 | p2) + plot_layout(guides = "collect",
                         axes = "collect") &
   theme(legend.position = "bottom",
         legend.margin = ggplot2::margin(t = -5,0,0,0,"pt"),
         plot.margin = ggplot2::margin(0,0,0,0,"pt"))) 
full_p + ggview::canvas(8,3.5)

# save
p1_edit <- p1 + 
  theme(legend.text = element_text(size = 7),
        legend.key.width= unit(13.5,"pt"),
        legend.key.height = unit(5,"pt")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top",
                                           legend.margin = ggplot2::margin(r = 20,t = -15,0,0, unit= "pt")))) +
  theme(legend.position = "bottom") +
  labs(title = NULL) 
p1_edit +   ggview::canvas(4,3.5)
# ggsave(p1_edit,
#        filename = "outputs/depth_shifts/ungrouped_all_spp.png",
#        width = 4,
#        height = 3.5)
# ggsave(p1_edit,
#        filename = "outputs/depth_shifts/ungrouped_all_spp.pdf",
#        width = 4,
#        height = 3.5,
#        device = cairo_pdf())

     
# save p2
p2_edit <- p2 + 
  theme(legend.text = element_text(size = 7),
        legend.key.width= unit(13.5,"pt"),
        legend.key.height = unit(5,"pt")) +
  guides(fill = guide_legend(theme = theme(legend.title.position = "top",
                                           legend.margin = ggplot2::margin(r = 20,t = -15,0,0, unit= "pt")))) +
  theme(legend.position = "bottom") +
  labs(title = NULL) 
p2_edit +   ggview::canvas(4,3.5)
# ggsave(p2_edit,
#        filename = "outputs/depth_shifts/ungrouped_nondomain_spp.png",
#        width = 4,
#        height = 3.5)
# ggsave(p2_edit,
#        filename = "outputs/depth_shifts/ungrouped_nondomain_spp.pdf",
#        width = 4,
#        height = 3.5,
#        device = cairo_pdf())


# save as rds
out <- list(p1_edit, p2_edit)
names(out) <- c("all_sp_unfiltered","all_sp_inside_domain")
saveRDS(out,
        file = "outputs/depth_shifts/depth_proportions.rds")


# plot all trends ---------------------------------------------------------

cols <- RColorBrewer::brewer.pal("Set1",n=4)

trends_p <- pa_hs %>%
  filter(organism %in% names(hs_quants_split)) %>%
  ggplot(aes(x = year,
             y = tidalheight)) +
  ggthemes::theme_few() +
  geom_point(size = .01, alpha = .2,
             position = position_jitter(width = .4, height = .15),
             color = "grey10",
             shape = 16) +
  
  geom_hline(yintercept = domain,
             linetype = "dashed",
             linewidth = .25) +
  facet_wrap(~organism,
             nrow = 5) +
  # add 95 percentile
 #geom_point(data = hs_quants,
 #           aes(y = height95), color = "darkmagenta",
 #           size = .5,
 #           stroke = .35,
 #           shape = 4) +
  geom_smooth(data = hs_quants,
              method= "lm",
              aes(x = year,y = height95,
                  group = NULL,
                  linetype = organism %in% mods_all_annotate$organism[mods_all_annotate$param == "height95" & mods_all_annotate$p_value <.05]), 
              color = cols[1],
              se = F,
              linewidth = .5) +
 # geom_point(data = hs_quants,
 #            aes(y = heightmed), color = "blue",
 #            size = .5,
 #            stroke = .35,
 #            shape = 3) +
  geom_smooth(data = hs_quants,
              method= "lm",
              aes(x = year,y = heightmed,
                  group = NULL,
                  linetype = organism %in% mods_all_annotate$organism[mods_all_annotate$param == "heightmed" & mods_all_annotate$p_value <.05]), 
              color = cols[2],
              se = F,
              linewidth = .5) +
 # geom_point(data = hs_quants,
 #            aes(y = heightmean), color = "green3",
 #            size = .5,
 #            stroke = .35,
 #            shape = 1) +
  geom_smooth(data = hs_quants,
              method= "lm",
              aes(x = year,y = heightmean,
                  group = NULL,
                  linetype = organism %in% mods_all_annotate$organism[mods_all_annotate$param == "heightmean" & mods_all_annotate$p_value <.05]), 
              color = cols[3],
              se = F,
              linewidth = .5) +
 # geom_point(data = hs_quants,
 #            aes(y = height05), color = "orange3",
 #            size = .5,
 #            stroke = .35,
 #            shape = 3) +
  geom_smooth(data = hs_quants,
              method= "lm",
              aes(x = year,y = height05,
                  group = NULL,
                  linetype = organism %in% mods_all_annotate$organism[mods_all_annotate$param == "height05" & mods_all_annotate$p_value <.05]), 
              color = cols[4],
              se = F,
              linewidth = .5) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(margin = ggplot2::margin(b = 0,t = 3, 0, 0, unit = "pt"),
                                  size = 7)) +
  labs(linetype = "Significant") +
  scale_linetype_manual(breaks = c(TRUE,FALSE),
                        labels = c("Yes","No"),
                        values = c("solid","11")) +
  
  # add significance
  geom_text(data = mods_all_annotate %>%
              filter(param == "height95",
                     p_value < .05),
            aes(x = 1981,
                y = -.15,
                label = sign),
            color = cols[1],
            vjust = 1,
            hjust = 0, 
            size = 2.5) +
  geom_text(data = mods_all_annotate %>%
              filter(param == "heightmed",
                     p_value < .05),
            aes(x = 1991.5,
                y = -.15,
                label = sign),
            color = cols[2],
            vjust = 1,
            hjust = 0, 
            size = 2.5) +
  geom_text(data = mods_all_annotate %>%
              filter(param == "heightmean",
                     p_value < .05),
            aes(x = 2002,
                y = -.15,
                label = sign),
            color = cols[3],
            vjust = 1,
            hjust = 0, 
            size = 2.5) +
  geom_text(data = mods_all_annotate %>%
              filter(param == "height05",
                     p_value < .05),
            aes(x = 2012.5,
                y = -.15,
                label = sign),
            color = cols[4],
            vjust = 1,
            hjust = 0, 
            size = 2.5) +
  coord_cartesian(xlim = c(1980, 2025),
                  ylim = c(-.8, 3.6),
                  expand = F) +
  scale_x_continuous(breaks = c(1980, 2000, 2020)) +
  
  geom_text(data = exclusions %>%
              filter(exclude_max == T & exclude_min == F),
            aes(x = 2002.5,
                y = 3.15,
                label = "Max at Domain"),
            size = 2) +
  
  geom_text(data = exclusions %>%
              filter(exclude_min == T & exclude_max == F),
            aes(x = 2002.5,
                y = 3.15,
                label = "Min at Domain"),
            size = 2) +
  
  geom_text(data = exclusions %>%
              filter(exclude_min == T & exclude_max == T),
            aes(x = 2002.5,
                y = 3.15,
                label = "Both at Domain"),
            size = 2) +
  theme(axis.text.x = element_text(size = 7,
                                   angle = -45,
                                   hjust = 0, 
                                   vjust = 0),
        axis.text.y = element_text(size = 7),
        legend.position = "bottom",
        #legend.position.inside = c(.52, 0),
        #legend.justification = c(0, 0),
        legend.box = "horizontal",
        legend.box.margin = ggplot2::margin(t = -15,0,0,0,"pt")) +
  
  geom_line(data = data.frame(
    param = c("height95","heightmed","heightmean","height05"),
    year = 1982,
    tidalheight = 1
  ),
  aes(color = param)) +
  scale_color_manual(breaks = c("height95",
                                "heightmed",
                                "heightmean",
                                "height05"),
                     labels = c("Max (95%)",
                                "Median",
                                "Mean",
                                "Min (5%)"),
                     values = c(cols)) +
  
  guides(linetype = guide_legend(override.aes = list(color = "grey20"),
                                 theme = theme(legend.title.position = "top")),
         color = guide_legend(nrow = 1,
                              override.aes = list(linewidth = 1),
                              theme = theme(legend.title.position = "top"))) +
  
  labs(x = NULL,
       y = "Intertidal Height (m)",
       color = "Occurrence Parameter") +
  theme(strip.clip = "off",
        strip.text = element_text(face = "italic"),
        plot.margin = ggplot2::margin(0,0,0,r=10)) 

  
  
trends_p + ggview::canvas(7,6)

ggsave(trends_p,
       filename = "outputs/depth_shifts/ungrouped_all_trends.png",
       width =7,
       height = 6)


# make table --------------------------------------------------------------
library(reactable)
library(reactablefmtr)
mods_all_annotate %>% filter(p_value <= .05) %>%
  distinct(organism)

mods_all_annotate %>%
  mutate(p_value = round(p_value,3),
         trend = round(trend * 10,3)) %>%
  select(-single_value, -sign) %>%
  mutate(exclude = case_when(exclude == T ~ "yes",
                                  exclude == F ~ "")) %>%

  tidyr::pivot_wider(
    names_from = param,
    values_from = c( trend, p_value, exclude),
    names_glue = "{param}_{.value}"
  ) %>%
  relocate(organism, contains("height95"), contains("heightmed"),
           contains("heightmean"), contains("height05")) %>%
  reactable(
    
    columns = list(
      organism = colDef(name = "Organism",
                        style = cell_style(font_style = "italic"),
                        width = 150),
      height95_trend = colDef(name = "Shift\n(m/Decade)",
                              style = function(value) {
                                color <- if (value > 0) "green" else if (value < 0) "red" else "black"
                                list(color = color)
                              }),
      height95_p_value = colDef(name = "p value",
                                cell = function(value) {
                                  if (value < 0.05) {
                                    paste0(value, " *")
                                  } else {
                                    value
                                  }
                                },
                                style = function(value) {
                                  fontWeight <- if (value < 0.05) "bold" else "normal"
                                  list(fontWeight = fontWeight)
                                }),
      height95_exclude = colDef(name = "Exclude;\ndomain edge"),
      heightmed_trend = colDef(name = "Shift\n(m/Decade)",
                               style = function(value) {
                                 color <- if (value > 0) "green" else if (value < 0) "red" else "black"
                                 list(color = color)
                               }),
      heightmed_p_value = colDef(name = "p value",
                                 cell = function(value) {
                                   if (value < 0.05) {
                                     paste0(value, " *")
                                   } else {
                                     value
                                   }
                                 },
                                 style = function(value) {
                                   fontWeight <- if (value < 0.05) "bold" else "normal"
                                   list(fontWeight = fontWeight)
                                 }),
      heightmed_exclude = colDef(name = "Exclude;\ndomain edge"),
      heightmean_trend = colDef(name ="Shift\n(m/Decade)",
                                style = function(value) {
                                  color <- if (value > 0) "green" else if (value < 0) "red" else "black"
                                  list(color = color)
                                }),
      heightmean_p_value = colDef(name = "p value",
                                  cell = function(value) {
                                    if (value < 0.05) {
                                      paste0(value, " *")
                                    } else {
                                      value
                                    }
                                  },
                                  style = function(value) {
                                    fontWeight <- if (value < 0.05) "bold" else "normal"
                                    list(fontWeight = fontWeight)
                                  }),
      heightmean_exclude = colDef(name = "Exclude;\ndomain edge"),
      height05_trend = colDef(name = "Shift\n(m/Decade)",
                              style = function(value) {
                                color <- if (value > 0) "green" else if (value < 0) "red" else "black"
                                list(color = color)
                              }),
      height05_p_value = colDef(name = "p value",
                                cell = function(value) {
                                  if (value < 0.05) {
                                    paste0(value, " *")
                                  } else {
                                    value
                                  }
                                },
                                style = function(value) {
                                  fontWeight <- if (value < 0.05) "bold" else "normal"
                                  list(fontWeight = fontWeight)
                                }),
      height05_exclude = colDef(name = "Exclude;\ndomain edge")
    ),
    
    columnGroups = list(
      colGroup(name = "Max Height (95%)", columns = c("height95_trend", "height95_p_value", "height95_exclude")),
      colGroup(name = "Median Height", columns = c("heightmed_trend", "heightmed_p_value", "heightmed_exclude")),
      colGroup(name = "Mean Height", columns = c("heightmean_trend", "heightmean_p_value", "heightmean_exclude")),
      colGroup(name = "Min Height (5%)", columns = c("height05_trend", "height05_p_value", "height05_exclude"))
      
    ),
    
    pagination = F
  )
  

  

library(dplyr)
library(tidyr)
library(gt)

mods_table <- mods_all_annotate %>%
  mutate(
    p_value = round(p_value, 3),
    trend = round(trend * 10, 3)
  ) %>%
  select(-single_value, -sign) %>%
  mutate(
    exclude = case_when(exclude == TRUE ~ "yes", TRUE ~ "")
  ) %>%
  pivot_wider(
    names_from = param,
    values_from = c(trend, p_value, exclude),
    names_glue = "{param}_{.value}"
  ) %>%
  relocate(
    organism,
    contains("height95"), contains("heightmed"),
    contains("heightmean"), contains("height05")
  ) %>%
  gt(rowname_col = "organism") %>%
  
  # Italicize organism names
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = 1)
  ) %>%
  
  # color text 
  tab_style(
    style = cell_text(color = "green3"),
    locations = cells_body(
      columns = c(height95_trend),
      rows = height95_trend > 0
    )
  ) %>%
  tab_style(
    style = cell_text(color = "green3"),
    locations = cells_body(
      columns = c(heightmed_trend),
      rows = heightmed_trend > 0
    )
  ) %>%
  tab_style(
    style = cell_text(color = "green3"),
    locations = cells_body(
      columns = c(heightmean_trend),
      rows = heightmean_trend > 0
    )
  ) %>%
  tab_style(
    style = cell_text(color = "green3"),
    locations = cells_body(
      columns = c(height05_trend),
      rows = height05_trend > 0
    )
  ) %>%
  
  # Red for negative trends
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(
      columns = c(height95_trend),
      rows = height95_trend < 0
    )
  ) %>%
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(
      columns = c(heightmed_trend),
      rows = heightmed_trend < 0
    )
  ) %>%
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(
      columns = c(heightmean_trend),
      rows = heightmean_trend < 0
    )
  ) %>%
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(
      columns = c(height05_trend),
      rows = height05_trend < 0
    )
  ) %>%
  
  # Add asterisks to significant p-values + bold
  text_transform(
    locations = cells_body(
      columns = c(
        height95_p_value, heightmed_p_value,
        heightmean_p_value, height05_p_value
      )
    ),
    fn = function(x) {
      ifelse(as.numeric(x) < 0.05,
             paste0(x, " *"),
             x)
    }
  ) %>%
  # height95_p_value
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(height95_p_value),
      rows = height95_p_value < 0.05
    )
  ) %>%
  # heightmed_p_value
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(heightmed_p_value),
      rows = heightmed_p_value < 0.05
    )
  ) %>%
  # heightmean_p_value
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(heightmean_p_value),
      rows = heightmean_p_value < 0.05
    )
  ) %>%
  # height05_p_value
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(height05_p_value),
      rows = height05_p_value < 0.05
    )
  ) %>%
  
  # Add spanner column groups
  tab_spanner(
    label = "Max Height (95%)",
    columns = c(
      height95_trend, height95_p_value, height95_exclude
    )
  ) %>%
  tab_spanner(
    label = "Median Height",
    columns = c(
      heightmed_trend, heightmed_p_value, heightmed_exclude
    )
  ) %>%
  tab_spanner(
    label = "Mean Height",
    columns = c(
      heightmean_trend, heightmean_p_value, heightmean_exclude
    )
  ) %>%
  tab_spanner(
    label = "Min Height (5%)",
    columns = c(
      height05_trend, height05_p_value, height05_exclude
    )
  ) %>%
  
  # Set nice column titles
  cols_label(
    organism = "Organism",
    height95_trend = html("Shift<br>(m/Decade)"),
    height95_p_value = "p value",
    height95_exclude = html("Exclude;<br>domain edge"),
    heightmed_trend = html("Shift<br>(m/Decade)"),
    heightmed_p_value = "p value",
    heightmed_exclude = html("Exclude;<br>domain edge"),
    heightmean_trend = html("Shift<br>(m/Decade)"),
    heightmean_p_value = "p value",
    heightmean_exclude = html("Exclude;<br>domain edge"),
    height05_trend = html("Shift<br>(m/Decade)"),
    height05_p_value = "p value",
    height05_exclude = html("Exclude;<br>domain edge")
  ) %>%
  
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_stub()
  )

mods_table
gtsave(mods_table, "outputs/depth_shifts/trends_table.rtf")


mods_all_annotate %>%
  filter(p_value < .05) %>%
  distinct(organism)





