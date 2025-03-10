
# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(lme4)
theme_set(ggthemes::theme_few())



# upload data -------------------------------------------------------------
pa_therm <- readr::read_csv("data-processed/cti-data/cti-data-all.csv")
pa_therm_filtered <- readr::read_csv("data-processed/cti-data/cti-data-highly-sampled.csv")
temptrend <- readr::read_csv("data-processed/cti-data/temp-change-augment.csv")
tempmod <- readRDS("data-processed/cti-data/temp-change-model.rds")
first_temp <- predict(tempmod, list(year= min(pa_therm$year)))
last_temp <- predict(tempmod, list(year= max(pa_therm$year)))
traits <- readr::read_csv("data-raw/spp_traits/functional_traits_completed.csv") %>%
  mutate(trait_group = case_when(motility_adult == "sessile" & group == "Invertebrate" ~ "Sessile Invertebrate",
                                 motility_adult == "motile" & group == "Invertebrate" ~ "Motile Invertebrate",
                                 group == "Algae" ~ "Algae"
  )) %>%
  rename(organism = gen_spp)

pa_therm <- pa_therm %>%
  left_join(traits %>% #rename(organism = gen_spp) %>%
              mutate(group = trait_group))
pa_therm_filtered <- pa_therm_filtered %>%
  left_join(traits %>% #rename(organism = gen_spp) %>%
              mutate(group = trait_group))


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

get_mod_p <- function(mod){
  
  full_p <- pf(summary(mod)$fstatistic[1], summary(mod)$fstatistic[2], summary(mod)$fstatistic[3], lower.tail = FALSE)
  
  mod_p_text <- NA
  if(full_p < 0.01){mod_p_text <- "Model p < 0.01"}
  else if(full_p < 0.05){mod_p_text <- "Model p < 0.05"}
  else if(full_p > 0.05){mod_p_text <- "Model p > 0.05"}
  
  term_p <- summary(mod)$coefficients
  if(any(stringr::str_detect(rownames(term_p),":"))){
    int_term <- term_p[which(stringr::str_detect(rownames(term_p),":")),"Pr(>|t|)"]
    int_term_text <- NA
    if(int_term < 0.01){int_term_text <- "Interaction p < 0.01"}
    else if(int_term < 0.05){int_term_text <- "Interaction p < 0.05"}
    else if(int_term > 0.05){int_term_text <- "Interaction p > 0.05"}
  }
  
  out <- list("modp" = mod_p_text,
              "term" = int_term_text)
  
  return(out)
}



# p1 highly sampled, per level ---------------------------------------------
df1 <- pa_therm_filtered %>%
  # find CTI per level/year
  group_by(level, tidalheight, year, group) %>%
  summarize(cti = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are 4+ species
  #filter(n_spp >= 4) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level))


df1 %>%   ggplot(aes(x = year, y = cti,
                     color = group, 
                     group = group)) +
  geom_smooth(method = "lm",
              se = F,
              color = "grey20") +
  facet_wrap(~tidalheight_f) +
  geom_point(alpha = .6) 

mod_g1 <- test_mods(df1 %>% filter(group == unique(df1$group)[1]))
summary(mod_g1)
mod_g2 <- test_mods(df1 %>% filter(group ==  unique(df1$group)[2]))
summary(mod_g2)
mod_g3 <- test_mods(df1 %>% filter(group ==  unique(df1$group)[3]))
summary(mod_g3)

mod_g1_p <- get_mod_p(mod_g1)
mod_g2_p <- get_mod_p(mod_g2)
mod_g3_p <- get_mod_p(mod_g3)
au_g1 <- broom::augment(mod_g1, 
                      newdata = tidyr::complete(df1 %>% filter(group == unique(df1$group)[1]),year,tidalheight, group), 
                      interval = "prediction")
au_g2 <- broom::augment(mod_g2, 
                        newdata = tidyr::complete(df1 %>% filter(group == unique(df1$group)[2]),year,tidalheight, group), 
                        interval = "prediction")
au_g3 <- broom::augment(mod_g3, 
                        newdata = tidyr::complete(df1 %>% filter(group == unique(df1$group)[3]),year,tidalheight, group), 
                        interval = "prediction")
au1 <- rbind(au_g1, au_g2, au_g3)

# plot
p1 <-  au1 %>% 
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, 
             group = interaction(tidalheight, group),
             color = group)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp,
             alpha = .5, 
             linewidth = .25) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash", 
             alpha = .5, 
             linewidth = .25) +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -.25, yend = -.25,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.2,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =2) +
  facet_wrap(~group) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% group_by(group) %>% filter(year == min(year)),
            aes(color = NULL,
                group = NULL,
                y  = tidalheight,
                x = .fitted)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
# # add model annotation
# annotate(geom = "richtext",
#          label = paste0("<b>Best Model:</b><br>",
#                         mod1$call[[2]], "<br>",
#                         mod1_p[[1]],"<br>",
#                         mod1_p[[2]]),
#          label.color = NA,
#          fill = NA,
#          x = 10.2,
#          y = 2.9,
#          size = 3,
#          hjust = 0,
#          vjust = 1) +
  coord_cartesian(xlim = c(10.05,12.15),
                  ylim = c(-.4,3),
                  expand = F) +
  
  theme(legend.position = c(.03,.98),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "CTI Change Across Depths for Organism Groups"
  ) +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(), 
        plot.title.position = "plot")

p1 + ggview::canvas(8,4)

# save as rds object to combine later
#saveRDS(p1,
#        "outputs/cti/bylevel_cti_p1.rds")




# p4 all data, per level ---------------------------------------------
df4 <- pa_therm %>%
  # find CTI per level/year
  group_by(level, tidalheight, year, group) %>%
  summarize(cti = mean(mean_monthly_mean,na.rm=T),
            cti_med = mean(median_monthly_mean, na.rm=T),
            n_obs = n(),
            n_spp = n_distinct(organism)) %>%
  ungroup() %>%
  # filter to levels where there are 4+ species
  filter(n_spp >= 4) %>%
  mutate(tidalheight_f = forcats::fct_reorder(as.character(tidalheight),
                                              level)) %>%
  # filter to levels where there are at least 3 years of sampling
  # (because we can't tell a trend when a level is only sampled once)
  group_by(level) %>%
  add_count() %>%
  filter(n >= 20) %>%
  select(-n) %>% ungroup()

df4 %>%   ggplot(aes(x = year, y = cti,
                     color = group)) +
  geom_smooth(method = "lm",
              se = F) +
  facet_wrap(~tidalheight_f) +
  geom_point(alpha = .6) 


mod4_g1 <- test_mods(df4 %>% filter(group == unique(df4$group)[1]))
summary(mod4_g1)
mod4_g2 <- test_mods(df4 %>% filter(group ==  unique(df4$group)[2]))
summary(mod4_g2)

mod4_g1_p <- get_mod_p(mod4_g1)
mod4_g2_p <- get_mod_p(mod4_g2)
au4_g1 <- broom::augment(mod4_g1, 
                        newdata = tidyr::complete(df4 %>% filter(group == unique(df4$group)[1]),year,tidalheight, group), 
                        interval = "prediction")
au4_g2 <- broom::augment(mod4_g2, 
                        newdata = tidyr::complete(df4 %>% filter(group == unique(df4$group)[2]),year,tidalheight, group), 
                        interval = "prediction")
au4 <- rbind(au4_g1, au4_g2)


p4 <- au4 %>% 
  arrange(year) %>%
  ggplot(aes(x = .fitted,
             y = tidalheight, group = tidalheight)) +
  
  # add real temp change
  geom_vline(xintercept = first_temp,
             alpha = .5, 
             linewidth = .25) +
  geom_vline(xintercept = last_temp,
             linetype = "longdash",
             alpha = .5, 
             linewidth = .25) +
  theme(legend.position.inside = c(.08,1)) +
  annotate(geom = "segment",
           x = first_temp + .05,
           xend = last_temp - .05,
           y = -1, yend = -1,
           arrow = arrow(#ends = "both",
             length = unit(5,"pt"))) +
  annotate(geom = "text",
           x = 10.85, 
           y = -.9,
           label = "True temperature change",
           hjust = .5,
           vjust =0,
           size =2) +
  
  # add model data
  geom_point(aes(x = cti,
                 color = year),
             shape= "|",
             size= 2.5,
             stroke = 2.5) +
  geom_path(data = . %>% filter(year == min(year)),
            aes(group = NULL)) +
  geom_path(arrow = arrow(length = unit(6,"pt"),
                          type = "open"),
            # position = position_nudge(y = .07),
            color = "black",
            linewidth = 1,
            linejoin = "mitre") +
  scale_color_gradient(low = "grey40",
                       high = "cyan3",
                       limits = c(1980, 2023),
                       breaks = c(1980, 2000, 2020)) +
  
  # add model annotation
 # annotate(geom = "richtext",
 #          label = paste0("<b>Best Model:</b><br>",
 #                         mod4$call[[2]], "<br>",
 #                         mod4_p[[1]],"<br>",
 #                         mod4_p[[2]]),
 #          label.color = NA,
 #          fill = NA,
 #          x = 10.2,
 #          y = 3.6,
 #          size = 3,
 #          hjust = 0,
 #          vjust = 1) +
 # coord_cartesian(xlim = c(10.05,12.15),
 #                 ylim = c(-.4,3.7),
 #                 expand = F) +
  
  theme(legend.position = c(.03,.98),
        legend.justification = c(0,1),
        legend.direction = "horizontal",
        legend.key.height = unit(5,"pt"),
        legend.key.width = unit(15,"pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 10)) +
  labs(x = "Community Thermal Index",
       y = "Shore Level (m)",
       color = "Sample Year",
       title = "All data, all replicates") +
  guides(color = guide_colorbar(
    theme = theme(legend.title.position = "top"))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(),
                     labels = ~paste0(.,"°")) +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(), 
        plot.title.position = "plot") +
  facet_wrap(~group)

p4 + ggview::canvas(8,4)
