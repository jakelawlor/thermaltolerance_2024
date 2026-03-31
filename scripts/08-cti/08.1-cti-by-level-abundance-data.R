# find mean thermal affinity over across levels
# by ABUNDANCE


# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(lme4)
theme_set(ggthemes::theme_few())



# upload data --------------------------------------------------------------
# get raw cover data
temptrend <- readr::read_csv("data-processed/cti-data/temp-change-augment.csv")
tempmod <- readRDS("data-processed/cti-data/temp-change-model.rds")
first_temp <- predict(tempmod, list(year= 1982))
last_temp <- predict(tempmod, list(year= 2023))

cover <- readr::read_csv("data-processed/appledore-survey-data/cover_abundance_filtered.csv")
cover %>% glimpse()
cover <- cover %>% 
  select(organism, year, transect, level, replicate, percent_cover,pc_num) %>%
  filter(!is.na(pc_num)) %>%
  mutate(tidalheight = (13.5-level)*.3048) 
# get STIs
therm <- readr::read_csv("data-processed/species-thermal-affinities/spp_thermal_affinities.csv")
therm %>% glimpse()
# get highly sampled quadrats
hs <- readr::read_csv("data-processed/cti-data/cti-data-highly-sampled.csv")
hs <- hs %>% distinct(year, transect, level)


# make function to test models --------------------------------------------
test_mods <- function(df){
  mod_add <- lm(data = df,
                formula = cti ~ year + tidalheight)
  mod_mult <- lm(data = df,
                 formula = cti ~ year * tidalheight)
  cti_add_par <- lm(data = df,
                    formula = cti ~ year + I(tidalheight^2)) 
  cti_mult_par <- lm(data = df,
                     formula = cti ~ year * I(tidalheight^2))
  perf <- performance::compare_performance(
    mod_add,
    mod_mult,
    cti_add_par,
    cti_mult_par,
    rank = T
  )
  
  best_mod <- get(perf$Name[[1]])
  print(paste0("best model:", perf$Name[[1]]))
  return(best_mod)
}




# process -----------------------------------------------------------------
df <- cover %>%
  #filter(replicate == 1) %>%
  inner_join(hs) %>%
  left_join(therm %>%
              select(organism,
                     perc_monthly_min_05, 
                     mean_monthly_mean, 
                     perc_monthly_max_95)) %>%
  tidyr::drop_na(mean_monthly_mean) %>%
  filter(pc_num != 0)
df %>% glimpse()



# summarize within levels -------------------------------------------------
df_sum <- df %>%
  group_by(year, tidalheight) %>%
  summarize(cti = weighted.mean(mean_monthly_mean, pc_num, na.rm=T),
            
            cti.min = weighted.mean(perc_monthly_min_05, 
                                    pc_num, na.rm=T),
            cti.max = weighted.mean(perc_monthly_max_95, 
                                    pc_num, na.rm=T),
            n_spp = n_distinct(organism)) #%>%
 # filter(n_spp >= 4)

df_sum %>%
  ggplot(aes(x = year,
             y = cti)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tidalheight)

df_sum %>%
  ggplot(aes(x = year,
             y = cti.min)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tidalheight)

df_sum %>%
  ggplot(aes(x = year,
             y = cti.max)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tidalheight)

df_sum %>%
  arrange(year) %>%
  ggplot(aes(x = cti,
             y = tidalheight)) +
  geom_point(shape = "|",
             aes(color = year),
             size = 4) +
  scale_y_reverse() 


mod1 <- lm(data = df_sum,
           cti ~ year + tidalheight)
summary(mod1)
mod2 <- lm(data = df_sum,
           cti ~ year * tidalheight)
summary(mod2)
mod3 <- lm(data = df_sum,
           cti ~ year + I(tidalheight^2))
summary(mod3)
mod4 <- lm(data = df_sum,
           cti ~ year * I(tidalheight^2))
summary(mod4)

performance::compare_performance(mod1, mod2, mod3, mod4,
                                 rank = T)
test_mods(df_sum)







# repeat with counts ------------------------------------------------------
counts <- readr::read_csv("data-processed/appledore-survey-data/counts_abundance_filtered.csv")
counts %>% glimpse()
counts <- counts %>% 
  select(organism, year, transect, level, replicate, count, countnum) %>%
  filter(!is.na(countnum)) %>%
  mutate(tidalheight = (13.5-level)*.3048) 


df_count <- counts %>%
 # filter(replicate == 1) %>%
  inner_join(hs) %>%
  left_join(therm %>%
              select(organism,
                     perc_monthly_min_05, 
                     mean_monthly_mean, 
                     perc_monthly_max_95)) %>%
  tidyr::drop_na(mean_monthly_mean) %>%
  filter(countnum != 0)



# first, raw cti per year
df_count_sum <- df_count  %>%
  group_by(year,tidalheight) %>%
  summarise(cti = weighted.mean(mean_monthly_mean,
                                countnum, na.rm=T),
            cti.min = weighted.mean(perc_monthly_min_05, 
                                    countnum, na.rm=T),
            cti.max = weighted.mean(perc_monthly_max_95, 
                                    countnum, na.rm=T),
            n_spp = n_distinct(organism)) #%>%
 # filter(n_spp >= 4)


# merge the two together --------------------------------------------------

df_merge <- df_count_sum %>% mutate(df = "count") %>%
  rbind(df_sum %>% mutate(df = "cover")) %>%
  group_by(year)  %>%
  filter(n_spp >= 4)

df_merge %>% 
  ggplot(aes(x = year,
             y = cti)) +
  geom_point(aes(color = df)) +
  geom_smooth(method = "lm") +
  facet_wrap(~tidalheight) 


formulas_list <- list(
  # first, test relationships with no accounting for dataset
  cti ~ year + tidalheight,
  cti ~ year * tidalheight,
  cti ~ year + I(tidalheight^2),
  cti ~ year * I(tidalheight^2),
  # then, with dataset added as a categorical factor
  cti ~ (year + tidalheight) + df,
  cti ~ (year * tidalheight) + df,
  cti ~ (year + I(tidalheight^2)) + df,
  cti ~ (year * I(tidalheight^2)) + df,
  # then, with dataset as an interacting categorical factor
  cti ~ (year + tidalheight) * df,
  cti ~ (year * tidalheight) * df,
  cti ~ (year + I(tidalheight^2)) * df,
  cti ~ (year * I(tidalheight^2)) * df,
  # with dataset interacting with only one additional factor
  cti ~ year + (tidalheight*df),
  cti ~ year + (I(tidalheight^2)*df),
  cti ~ (year*df) + tidalheight,
  cti ~ (year*df) + I(tidalheight^2)

) %>%
  purrr::set_names(paste0("mod",c(1:length(.))))


mods <- purrr::map(
  .x = formulas_list,
  .f = ~lm(data = df_merge, 
           formula = .x)
) 

performance::compare_performance(mods, rank = T)


formulas_list$mod13
mergemod <- mods$mod13
summary(mergemod)
mergemod
# 
library(performance)
perf <- plot(check_model(mergemod))
library(patchwork)
perf <- perf & labs(subtitle = NULL)
perf + ggview::canvas(12,8)

ggsave(perf, 
      filename = "outputs/cti/cti_by_depth_assumptions.png",
      width = 12,
      height = 8)


pred_df <- data.frame(
  tidalheight = rep(unique(df_merge$tidalheight), times = length(unique(df_merge$year))),
  year = rep(unique(df_merge$year), each = length(unique(df_merge$tidalheight)))
) 
pred_df <- pred_df %>% mutate(df = "cover") %>%
  rbind(pred_df %>% mutate(df = "count")) %>%
  left_join(df_merge)

au_merge <- broom::augment(mergemod,
                           newdata = pred_df)

# make palette manually = gradient for two datasets individually 
# starting from different cyans
pal.df <- data.frame(year = c(1982:2023))
pal.df$pal1 <- colorRampPalette(c("grey10","magenta"))(length(1982:2023))
pal.df$pal2 <- colorRampPalette(c("grey10","cyan"))(length(1982:2023))
pal.df <- pal.df %>%
  tidyr::pivot_longer(cols = c("pal1","pal2"),
                      values_to = "pal") %>%
  mutate(df = case_when(name == "pal2" ~ "Percent Cover Species",
                        name == "pal1" ~ "Density Species"))


p <- 
  au_merge %>% 
  mutate(df = recode(df,
                     "cover" = "Percent Cover Species",
                     "count" = "Density Species")) %>%
    left_join(pal.df) %>%
 #left_join(df_merge) %>%
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, 
             group = interaction(tidalheight,df))) +
    # add background rectangle
    annotate(geom = "rect",
             xmin = first_temp,
             xmax = last_temp,
             ymin = -Inf,
             ymax = Inf,
             fill = "grey95") +
  geom_vline(xintercept = c(first_temp, last_temp),
             linetype = "twodash",
             alpha = .7,
             linewidth = .3) +
  # add percent cover dots
  geom_point(data = . %>% filter(df == "Percent Cover Species"),
             aes(x = cti,
                 color = pal),
             shape= "\\",
            #shape = 2,
             size= 1,
            stroke = 2.5,
             #stroke = 1.5,
            # alpha= .7,
             position = position_nudge(y = .04)) +
    # add density cover dots
    geom_point(data = . %>% filter(df == "Density Species") ,
               aes(x = cti,
                   color = pal),
               shape= "/",
               size= 1,
               stroke = 2.5,
               #shape = 1,
              # stroke =1.5,
              # alpha = .7,
               position = position_nudge(y = -.04)) +
    scale_color_identity() +
    # add min path percent cover
  geom_path(data = . %>% filter(year == min(year),
                                df == "Percent Cover Species"),
            aes(group = df),
            color = colorspace::darken("cyan3",.2),
            linewidth = .7,
           # position = position_nudge(y = .05)
           ) +
    # add min path percent cover
    geom_path(data = . %>% filter(year == min(year),
                                  df == "Density Species"),
              aes(group = df),
              color = "magenta3",
              linewidth = .7,
             # position = position_nudge(y = -.05)
             ) +
    
    # add percent cover arrows
    geom_path(data = . %>% filter(df == "Percent Cover Species"),
              arrow = arrow(length = unit(6,"pt"),
                            type = "open"),
              # position = position_nudge(y = .07),
              color = colorspace::darken("cyan3",.2),
              linewidth = .7,
              linejoin = "mitre",
              #position = position_nudge(y = .05)
              ) +
    # add percent cover arrows
  geom_path(data = . %>% filter(df == "Density Species"),
            arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "magenta3",
            linewidth = .7,
            linejoin = "mitre",
           # position = position_nudge(y = -.05)
           ) +
 # theme(legend.position = c(.47,.01),
 #       legend.justification = c(1,0),
 #       legend.direction = "horizontal",
 #       legend.key.height = unit(5,"pt"),
 #       legend.key.width = unit(15,"pt"),
 #       legend.text = element_text(size = 7),
 #       legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index (°C)",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "CTI Across Shore Levels") +
 # guides(color = guide_colorbar(
 #   theme = theme(legend.title.position = "top"))
 # ) +
  # add real temp change
  #geom_vline(xintercept = first_temp) +
  #geom_vline(xintercept = last_temp,
  #           linetype = "longdash") +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           linewidth = .3,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.2,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =3) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank()) +
  # add fake legend
    annotate(geom = "text",
             label = 
             "Percent Cover\nSpecies\n\nDensity Species",
             x = 12.5,
             y = 2.6,
             lineheight = .85,
             size = 2.5,
             hjust = 0,
             vjust = 1
             ) +
    annotate(geom = "point",
             x = 12.45,
             y = 2.52,
             shape = "\\",
             size = 2,
             stroke = 2.5,
    ) + 
    annotate(geom = "point",
             x = 12.45,
             y = 2.23,
             shape = "/",
             size = 2,
             stroke = 2.5,
    ) +
    
    ggview::canvas(4,4)
  
  p
# make blank plot with color scales
  # Now: create dummy data for each colorbar legend
  legend_df_a <- data.frame(year = unique(df_merge$year), 
                            value = 0, group = "Density Species")
  legend_df_b <- data.frame(year = unique(df_merge$year), 
                            value = 0, group = "Percent Cover Species")

  p_legend <- ggplot() +
    # These are transparent points used to generate the gradient legends
    geom_point(data = legend_df_a, aes(x = year, y = value, color = year), alpha = 0) +
    geom_point(data = legend_df_b, aes(x = year, y = value, fill = year), shape = 21, alpha = 0) +
    
    # Add color and fill scales that look like your gradients
    scale_color_gradientn(
      colors = colorRampPalette(c("grey10", "cyan"))(100),
      name = "Year\n(Percent Cover)"
    ) +
    scale_fill_gradientn(
      colors = colorRampPalette(c("grey10", "magenta"))(100),
      name = "Year\n(Density)"
    ) +
    theme_void() +
    guides(
      color = guide_colorbar(order = 1,
                             theme = theme(legend.title.position = "top")),
      fill = guide_colorbar(order = 2,
                            theme = theme(legend.title.position = "top"))
    ) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom")
  
  p_legend
  
# merge 
  library(patchwork)
  
p_full <-  (p | p_legend) + plot_layout(guides = "collect",
                             widths = c(1,0)) & theme(legend.position = 'bottom')
  
p_full2 <- p_full + 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.width = unit(5,"mm"),
        legend.key.height = unit(2.5,"mm"),
        legend.position = "bottom") 
p_full2 +ggview::canvas(4,5)

saveRDS(p_full2,
        file = "outputs/cti/p_by_level_by_abundance.rds")


