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
  mutate(tidalheight = (13-level)*.348)  
rm(pa_therm_filtered)


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
ggsave(p1_edit,
       filename = "outputs/depth_shifts/ungrouped_all_spp.png",
       width = 4,
       height = 3.5)
ggsave(p1_edit,
       filename = "outputs/depth_shifts/ungrouped_all_spp.pdf",
       width = 4,
       height = 3.5,
       device = cairo_pdf())

     
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
ggsave(p2_edit,
       filename = "outputs/depth_shifts/ungrouped_nondomain_spp.png",
       width = 4,
       height = 3.5)
ggsave(p2_edit,
       filename = "outputs/depth_shifts/ungrouped_nondomain_spp.pdf",
       width = 4,
       height = 3.5,
       device = cairo_pdf())


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
        strip.text = element_text(margin = margin(b = 0,t = 3, 0, 0, unit = "pt"),
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
        legend.box.margin = margin(t = -15,0,0,0,"pt")) +
  
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
       y = "Tidal Height (m)",
       color = "Occurrence Parameter") +
  theme(strip.clip = "off",
        strip.text = element_text(face = "italic"),
        plot.margin = margin(0,0,0,0)) 

  
  
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






# -------------------------------------------------------------------------
# SCRATCH ---------------------------------------------------------------------
# -------------------------------------------------------------------------


# plot all ----------------------------------------------------------------
library(ggtext)
hs_unc_trends_p <- hs_unc %>%
  filter(organism %in% names(hs_unc_sum_split)) %>%
  ggplot(aes(x = (year),
             group = year,
             y = tidalheight)) +
  # geom_boxplot(outliers = F,
  #              color = "grey50",
  #              alpha = .5,
  #              linewidth = .2,
  #              width = .3) +
  facet_wrap(~organism) +
  geom_point(data = hs_unc_sum,
             aes(y = height95), color = "black",
             size = .5,
             stroke = .2,
             shape = 4) +
  geom_smooth(data = hs_unc_sum,
              method= "lm",
              aes(x = year,y = height95,
                  group = NULL,
                  color = organism %in% hs_unc_all_p$organism[hs_unc_all_p$max_trend > 0]), 
              #color = "darkmagenta",
              se = F,
              size = .25) +
  geom_point(data = hs_unc_sum,
             aes(y = height05), color = "black",
             size = .5,
             shape = 3,
             stroke = .2) +
  geom_smooth(data = hs_unc_sum,
              method= "lm",
              aes(x = year,y = height05,
                  group = NULL,
                  color = organism %in% hs_unc_all_p$organism[hs_unc_all_p$min_trend > 0]), 
              # color = "darkturquoise",
              # color = ifelse(organism %in% hs_unc_all_p$organism[hs_unc_all_p$min_trend > 0],
              #                "blue","red"),
              se = F,
              size = .25) +
  geom_point(data = hs_unc_sum,
             aes(y = heightmed), color = "black",
             size = .5,
             shape = 10,
             stroke = .2) +
  geom_smooth(data = hs_unc_sum,
              method= "lm",
              aes(x = year,y = heightmean,
                  group = NULL,
                  color = organism %in% hs_unc_all_p$organism[hs_unc_all_p$med_trend > 0]), 
              #color = "black",
              se = F,
              size = .25) +
  
  # add max annotation
  geom_richtext(data = hs_unc_all_p,
                aes(x = 1980, y = 0,
                    label = max_p_sig,
                    group = NULL),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =1,
                hjust = 0) +
  
  geom_richtext(data = hs_unc_all_p,
                aes(x = 2002.5, y = 0,
                    label = med_p_sig,
                    group = NULL),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =1,
                hjust = .5) +
  
  geom_richtext(data = hs_unc_all_p,
                aes(x = 2025, y = 0,
                    label = min_p_sig,
                    group = NULL),
                fill = NA,
                label.colour = NA,
                lineheight = .8,
                size  = 2.5,
                vjust =1,
                hjust = 1) +
  
  # add bounds
  geom_hline(yintercept = max(hs_unc$tidalheight),
             linewidth = .25, linetype = "dashed") +
  geom_hline(yintercept = min(hs_unc$tidalheight),
             linewidth = .25, linetype = "dashed") +
  
  theme(strip.text = element_text(margin = margin(0,0,0,0, "pt")),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_rect(color = "transparent"),
        panel.border = element_rect(linewidth = .25))+
  coord_cartesian(xlim = c(1980, 2025),
                  ylim = c(-1,3),
                  expand = F)+
  labs(x = "Year Group",
       y = "Tidal Level (m)",
       title = "Species Depth Shifts - Highly Sampled Data",
       color = "Sign",
       shape = "Parameter") +
  scale_color_manual(breaks = c(TRUE,FALSE),
                     labels = c("Upward Shift","Downward Shift"),
                     values = c("darkcyan","darkmagenta")) +
  geom_point(data = data.frame(cat = c("Median","5th %","95th %"),
                               year = NA_integer_,
                               y = NA_integer_),
             inherit.aes=F,
             aes(x = year, 
                 y = y,
                 shape = cat)) +
  scale_shape_manual(breaks = c("95th %","Median","5th %"),
                     values = c(4,10,3))

hs_unc_trends_p + ggview::canvas(12,8)
# ggsave(hs_unc_trends_p, 
#        file = "outputs/depth_shifts/hs_ungrouped_trends_p.png",
#        width = 12,
#        height = 8)



hs_unc_all_p2 <- hs_unc_all_p %>% select(-ends_with("sig"), -any_p) %>%
  tidyr::pivot_longer(cols = c(mean_p, mean_sign, 
                               max_p, max_sign,
                               min_p, min_sign,
                               med_p, med_sign),
                      names_to = c("stat", ".value"),
                      names_sep = "_") %>%
  mutate(sig = case_when(p < .05 & sign == "(+)" ~ "Sig Shift Up",
                         p > .05 & sign == "(+)" ~ "Shift Up",
                         p < .05 & sign == "(-)" ~ "Sig Shift Down",
                         p > .05 & sign == "(-)" ~ "Shift Down",
                         sign == "None" ~ "No Change")) %>%
  mutate(sig = factor(sig, levels = c("Sig Shift Up","Shift Up","No Change","Shift Down","Sig Shift Down"))) %>%
  mutate(stat = factor(stat, levels = c("min","med","mean","max"))) 



hs_unc_depth_bars <- hs_unc_all_p2 %>%
  filter(stat != "mean") %>%
  ggplot(aes(x = stat,
             fill = sig)) +
  geom_bar(color = "black", width = 1,
           linewidth = .2) +
  coord_flip(ylim = c(0,length(unique(hs_unc_all_p2$organism))),
             xlim = c(.5, 3.5),
             expand = F) +
  scale_x_discrete(labels = ~paste(stringr::str_to_title(.),"\nObservation \nHeight")) +
  scale_fill_manual(breaks = c("Sig Shift Down",
                               "Shift Down",
                               "No Change",
                               "Shift Up",
                               "Sig Shift Up"),
                    values = c("darkmagenta",
                               colorspace::lighten("darkmagenta",.8),
                               "grey95",
                               colorspace::lighten("darkturquoise",.8),
                               "darkturquoise")) +
  labs(x = NULL,
       y = "Number of Species",
       fill = "Modeled Shift") +
  theme(legend.position = "bottom")  +
  geom_hline(yintercept = length(unique(hs_unc_all_p2$organism))/2,
             linetype = "dashed") +
  guides(fill = guide_legend(theme= theme(legend.title.position = "top"))) +
  theme(legend.margin = margin(t = -15,0,0,0,"pt"))

hs_unc_depth_bars + ggview::canvas(6.2,5)
#ggsave(hs_depth_bars,
#       file = "outputs/depth_shifts/hs_bars.png",
#       width = 6.2, height = 5)
#


# make table of trends and significance
library(reactable)
hs_unc_all_p %>%
  select(organism,
         min_trend, #min_p,
         med_trend, #med_p,
         max_trend, #max_p
  ) %>%
  mutate(min_trend = min_trend * 10,
         med_trend = med_trend * 10,
         max_trend = max_trend * 10) %>%
  mutate(organism = paste0("<i>",organism,"</i>")) %>%
  reactable(
    pagination = F,
    defaultColDef = colDef(format = colFormat(digits = 2),
                           html  = T),
    columns = list(
      organism = colDef(name = "Organism",
                        style = list(font_style = "italic")),
      min_trend = colDef(name = "Lower Limit Shift (m/decade)",
                         style = function(value, index){
                           color <- ifelse(value < 0,"magenta","cyan")
                           bold <- ifelse(hs_unc_all_p$min_p[index] <= 0.05,"bold","normal")
                           list(background = alpha(color,.4), fontWeight = bold)
                         },
                         format = colFormat(digits = 2)
      ),
      med_trend = colDef(name = "Distribution Median Shift (m/decade)",
                         style = function(value, index){
                           color <- ifelse(value < 0,"magenta","cyan")
                           bold <- ifelse(hs_unc_all_p$med_p[index] <= 0.05,"bold","normal")
                           list(background = alpha(color,.4), fontWeight = bold)
                         },
                         format = colFormat(digits = 2)
      ),
      max_trend = colDef(name = "Upper Limit Shift (m/decade)",
                         style = function(value, index){
                           color <- ifelse(value < 0,"magenta","cyan")
                           bold <- ifelse(hs_unc_all_p$max_p[index] <= 0.05,"bold","normal")
                           list(background = alpha(color,.4), fontWeight = bold)
                         },
                         format = colFormat(digits = 2)
      )
      # med_trend = colDef(name = "")
    )
  )



# Example dataset
df <- data.frame(
  species = c("Species A", "Species B", "Species C", "Species D"),
  coeff = c(1.5, -2.3, 0.7, -1.8),
  p_value = c(0.04, 0.08, 0.03, 0.001)
)

# Create reactable table
# Create reactable table
reactable(df, 
          columns = list(
            species = colDef(name = "Species", width = 200),
            coeff = colDef(
              name = "Model Coefficient",
              style = function(value, index) {
                # Color based on the sign of the coefficient
                color <- ifelse(value > 0, "blue", "red")
                
                # Make the coefficient bold if the p-value is significant
                bold <- ifelse(df$p_value[index] < 0.05, "bold", "normal")
                
                list(color = color, fontWeight = bold)
              }
            )
          )
)


# -------------------------------------------------------------------------
# scratch -----------------------------------------------------------------
# -------------------------------------------------------------------------



unclustered <- hs_by_yeargroup %>%
  group_by(organism, year) %>%
  summarize(min = quantile(tidalheight, .05),
            med = quantile(tidalheight, .5),
            max = quantile(tidalheight,.95)) %>%
  group_by(organism) %>%
  mutate(n_years = n_distinct(year)) %>%
  filter(n_years > 10)  %>%
  split(f = .$organism)

i <- 23
unclustered[[i]] %>%
  ggplot(aes(x = year, 
             y = min)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", 
              se=F) +
  geom_point(aes(y = med), 
             color = "black") +
  geom_smooth(aes(y = med),
              color = "black",
              method = "lm",
              se=F) +
  geom_point(aes(y = max), 
             color = "red",) +
  geom_smooth(aes(y = max),
              color = "red",
              method = "lm", 
              se = F) +
  labs(title = unique(unclustered[[i]]$organism))

summary(lm(data = unclustered[[i]],
           formula = max ~ year))

unclustered_max_mods <- purrr::map(
  .x = unclustered,
  .f = ~lm(data = .x,
           formula = "max ~ year")
)
unclustered_max_p <- purrr::map(
  .x = unclustered_max_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
) %>%
  bind_rows(.id = "organism") %>%
  rename(med_p = value)
hs_med_mods_trend <- purrr::map(
  .x = hs_med_mods,
  .f = ~ data.frame(value = ifelse(sign(summary(.x)$coefficients["year_group","Estimate"]) == 1,
                                   "(+)",
                                   "(-)"
  )) )%>%
  bind_rows(.id = "organism") %>%
  rename(med_trend = value)

# test for change in sampling over time -----------------------------------
samples_by_yeargroup <- org_by_yeargroup %>%
  distinct(year_group, tidalheight) %>%
  group_by(year_group) %>%
  
  summarize(height95 = quantile(tidalheight,.95),
            heightmed = median(tidalheight),
            heightmean = mean(tidalheight),
            height05 = quantile(tidalheight,.05)) 

year_max_mod <- lm(data = samples_by_yeargroup,
                   formula = "height95 ~year_group")
summary(year_max_mod)
year_mean_mod <- lm(data = samples_by_yeargroup,
                    formula = "heightmean ~year_group")
summary(year_mean_mod)
year_med_mod <- lm(data = samples_by_yeargroup,
                   formula = "heightmed ~year_group")
summary(year_med_mod)
year_min_mod <- lm(data = samples_by_yeargroup,
                   formula = "height05 ~year_group")
summary(year_min_mod)

org_by_yeargroup %>%
  ggplot(aes(x = as.character(year_group),
             y = tidalheight)) +
  geom_boxplot(outliers = F)

org_by_yeargroup %>% 
  group_by(year_group, tidalheight) %>%
  summarize(n_obs = n(),
            n_tran = n_distinct(transect)) %>% 
  ggplot(aes(x = year_group,
             y = tidalheight)) +
  geom_tile(aes(fill = n_tran ))

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------



mean_mods <- purrr::map(
  .x = org_by_height,
  .f = ~lm(data = .x, 
           formula = "heightmean ~ year_group")
)
mean_mods_p <- purrr::map(
  .x = mean_height_mods,
  .f = ~ pf(summary(.x)$fstatistic[1], summary(.x)$fstatistic[2], summary(.x)$fstatistic[3], lower.tail = FALSE)
)

mean_height_mods$`Crepidula fornicata` %>% summary()

mod <- lm(data = org_by_height[[1]],
          formula = "heightmean ~ year_group")
summary(mod)

ggplot(aes(x = year_group,
           y = heightmean)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "blue",
              se = F) +
  geom_point(aes(y = height95), color = "red") +
  geom_smooth(method = "lm", se = F, aes(y = height95),
              color = "red") +
  geom_point(aes(y = height05), color = "green") +
  geom_smooth(method = "lm", se = F, aes(y = height05),
              color = "green") +
  facet_wrap(~organism)






# -------------------------------------------------------------------------
# Scratch / old ---------------------------------------------------------------------
# -------------------------------------------------------------------------





therm <- readr::read_csv(
  here::here(
    "data-processed",
    "spp_thermal_affinities.csv"
  )
)

pa <- pa %>% left_join(therm) %>%
  filter(pres == T) %>%
  # translate level into tidal height
  mutate(tidalheight = (13-level)*.348)  
pa %>% distinct(data_taken)


# find cti of new spp -----------------------------------------------------
# first, let's break up the timeseries 
range(pa$year)
pa %>% 
  group_by(organism) %>%
  filter(year == min(year)) %>% 
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_group = case_when(
    year <= 1989 ~ "1980s",
    year >= 1990 & year <= 1999 ~ "1990s",
    year >= 2000 & year <= 2009 ~ "2000s",
    year >= 2010 & year <= 2019 ~ "2010s",
    year >= 2020 ~ "2020s"
  )) %>%
  ggplot(aes(x = year,
             y = mean_monthly_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year of First Occurrence",
       y = "Thermal Affinity") 



pa %>% 
  group_by(organism) %>%
  filter(year == max(year)) %>% 
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_group = case_when(
    year <= 1989 ~ "1980s",
    year >= 1990 & year <= 1999 ~ "1990s",
    year >= 2000 & year <= 2009 ~ "2000s",
    year >= 2010 & year <= 2019 ~ "2010s",
    year >= 2020 ~ "2020s"
  )) %>%
  filter(year != max(year)) %>%
  ggplot(aes(x = year,
             y = mean_monthly_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year of First Occurrence",
       y = "Thermal Affinity") 




# find how long species lasted in the dataset -----------------------------
first_year <- pa %>%
  group_by(organism) %>%
  filter(year == min(year)) %>%
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_type = "first_year")

last_year <- pa %>%
  group_by(organism) %>%
  filter(year == max(year)) %>%
  select(organism, 
         year, 
         pres,
         contains("monthly")) %>%
  distinct() %>%
  ungroup() %>%
  mutate(year_type = "last_year")

year_breadth <-  
  first_year %>%
  rbind(last_year)


year_breadth %>% 
  tidyr::pivot_wider(names_from = year_type,
                     values_from = year) %>% 
  filter(last_year - first_year > 0) %>% 
  ggplot() +
  geom_segment(aes(x = first_year,
                   xend = last_year,
                   y = mean_monthly_mean,
                   yend = mean_monthly_mean),
               linewidth = .3,
               alpha = .4) +
  geom_point(data = . %>% filter(first_year > min(first_year)),
             aes(x = first_year,
                 y = mean_monthly_mean),
             color = "cyan4") +
  geom_point(data = . %>% filter(last_year < max(last_year)),
             aes(x = last_year,
                 y = mean_monthly_mean),
             color = "tomato2") +
  geom_smooth(data = . %>% filter(first_year > min(first_year)),
              aes(x = first_year,
                  y = mean_monthly_mean),
              color = "cyan4",
              method = "lm") +
  geom_smooth(data = . %>% filter(last_year < max(last_year)),
              aes(x = last_year,
                  y = mean_monthly_mean),
              color = "tomato2",
              method = "lm")

pa %>%
  distinct(transect, level, tidalheight, year) %>%
  ggplot(aes(x = year,
             y = level)) +
  geom_tile() +
  facet_wrap(~transect)
# longest most consistent transects are:
# 5, 7, 15, 20, 22, 26, 28
# between levels 0 and 15
good_transects <- c(5, 7, 15, 20, 22, 26, 28)
min_level <- 0
max_level <- 15

# track individual species min/max depths ---------------------------------

pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>% 
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  ggplot(aes(x = year,
             group = organism)) +
  geom_line(aes(y = minheight),
            color = "grey80",
            linewidth = .2,
            alpha = .3) +
  geom_smooth(aes(y = minheight,
                  group = organism),
              method = "lm",
              se = F,
              color = "cyan4",
              linetype = "dashed",
              linewidth = .5) 

pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>% 
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight),
            mean_monthly_mean = unique(mean_monthly_mean)) %>%
  group_by(organism) %>%
  add_count() %>%
  ungroup() %>%
  filter(!is.na(mean_monthly_mean),
         n > 3) %>%
  mutate(organism_f = forcats::fct_reorder(organism, mean_monthly_mean)) %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = minheight),
            color = "tomato2",
            linewidth = .9,
            alpha = .3) +
  geom_line(aes(y = maxheight),
            color = "cyan4",
            linewidth = .9,
            alpha = .3) +
  # geom_smooth(aes(y = maxheight),
  #             method = "lm",
  #             se = F,
  #             color = "tomato2",
  #             linetype = "dashed",
  #             linewidth = .5) +
  facet_wrap(~organism_f)


pa %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>%   
  group_by(organism, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  ggplot(aes(x = year,
             group = organism)) +
  geom_ribbon(aes(ymin = minheight,
                  ymax = maxheight)) +
  facet_wrap(~organism)


pa %>%
  filter(!is.na(mean_monthly_mean)) %>%
  mutate(organism_f = forcats::fct_reorder(organism, mean_monthly_mean)) %>%
  filter(transect %in% good_transects,
         level >= min_level,
         level <= max_level) %>% 
  group_by(organism, organism_f, year) %>%
  summarize(minheight = min(tidalheight),
            maxheight = max(tidalheight)) %>%
  group_by(organism) %>%
  add_count() %>%
  ungroup() %>%
  filter(n > 3) %>%
  ggplot(aes(x = year)) +
  geom_smooth(aes(y = maxheight),
              method = "lm", 
              se = F,
              color = "tomato2") +
  geom_smooth(aes(y = minheight),
              method = "lm", 
              se = F,
              color = "cyan4") +
  facet_wrap(~organism_f)
