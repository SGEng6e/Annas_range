# Visualization code (Figures 1, S9, and S10) for:
# Supplemental feeding as a driver of population expansion and morphological change in Anna’s Hummingbirds

# NM Alexandre, FG Romero, SG English, E Grames, F Garzón-Agudelo, K Epperly, T Barnes, DR Powers, 
## AE Smith, Z Migicovsky, L Stein, S Akalu, H Sridhar, G Montross, E Collins, A Rico-Guevara

##### LIBRARIES #############################################################################################

library(tidyverse)
library(terra)
library(tidyterra)
library(MCMCvis)
library(coda)
library(nimble)
library(ggpubr)
library(magick)
library(gridExtra)

thm_obj <- theme(
  axis.text = element_text(size = 11),
  panel.background = element_blank(),
  strip.text = element_text(size = 12),
  axis.line = element_line(colour = "black"),
  legend.title = element_text(vjust = 2.5, size = 13, hjust = 1),
  legend.text = element_text(face = "bold"),
  legend.background = element_blank()
)

##### LOAD DATA #############################################################################################

load("./nimble_data.rda")
load("./nimble_constants.rda")

load("./anna_glmnuts_chain1.rda")
load("./anna_glmnuts_chain2.rda")
load("./anna_glmnuts_chain3.rda")
load("./anna_glmnuts_chain4.rda")
annas_chains <- mcmc.list(mget(ls(pattern = "chain\\d"))) ; rm(list = ls(pattern = "chain\\d"))

##### DATA SUMMARY ##########################################################################################

quants <- c(0.05,0.95)

marginals <- set_names(c("baseline","marginal_fa","marginal_ea","marginal_hp","marginal_yr"),
                       c("baseline","fa","^ea","hp","year"))

raw_dat <- imap(marginals,
                \(x,y) {
                  if (x == "baseline") {
                    tibble(
                      Term = paste0("baseline[",nimble.data$latlon_raw[,1],"]"),
                      route = nimble.data$latlon_raw[,1],
                      lat = nimble.data$latlon_raw[,2],
                      lon = nimble.data$latlon_raw[,3],
                      raw = NA)
                  } else {
                    dat <- unlist(nimble.data[which(str_detect(names(nimble.data), y))])
                    tibble(
                      Term = paste0(x,"[",1:length(dat), "]"),
                      route = nimble.constants$route,
                      year = nimble.data$year_raw,
                      lat = nimble.data$Latitude,
                      lon = nimble.data$Longitude,
                      raw = dat)
                  }
                }) %>% 
  bind_rows()

main_result_summary <- MCMCsummary(
  map(annas_chains, \(x) x[,which(!str_detect(colnames(x), "(N)|(marginal)|(baseline)|(ltln)"))]),
  probs = quants, pg0 = T, Rhat = T) %>% 
  rownames_to_column("Term") %>% 
  mutate(Term_name = factor(
    Term, levels = rev(c("beta_ea","beta_fa","beta_hp","beta_yr","beta_lt","beta_in","beta_ef","thet_ef",
                         "sd_ovds","Rn1","Rn2","bayes.p")),
    labels = rev(c("Eucalyptus\navailability","Feeder\navailability","Human\npopulation","Year",
                   "Year-by-\nlatitude","Abundance\nintercept","Effort\nslope","Effort\nsaturation",
                   "Overdispersion","Observed\nError","Simulated\nError","Bayesian\nP-value"))),
    signif = factor(case_when(`p>0` > 0.95 ~ 1,
                              `p>0` < 0.05 ~ -1,
                              `p>0` > 0.9 ~ 0.5,
                              `p>0` < 0.1 ~ -0.5,
                              between(`p>0`, 0.1, 0.9) ~ 0),
                    levels = paste(seq(-1,1, by = 0.5))))

GAM_component_summary <- MCMCsummary(
  map(annas_chains, \(x) x[,which(str_detect(colnames(x),"ltln"))]),
  probs = quants, pg0 = T, Rhat = T) %>% 
  rownames_to_column("Term")

marginal_effects_summary <- MCMCsummary(
  map(annas_chains, \(x) x[,which(str_detect(colnames(x), "(marginal)|(baseline)"))]),
  probs = quants, pg0 = T, Rhat = T) %>% 
  rownames_to_column("Term") %>% 
  left_join(raw_dat, by = "Term") %>% 
  mutate(Term_id = str_extract(Term, "(?<=\\[)\\d{1,4}(?=\\])"),
         Term = str_extract(Term, ".*(?=\\[)")) %>% 
  distinct(Term, lon, lat, route, raw, year, .keep_all = T) %>% 
  mutate(Term_name = factor(
    Term, levels = rev(c("baseline","marginal_ea","marginal_fa","marginal_hp","marginal_yr")),
    labels = rev(c("Baseline","Eucalyptus","Feeders","Human\npopulation","Year")))) %>% 
  relocate(lon, lat, .before = 1) %>% 
  pivot_wider(id_cols = c(lat,lon,route,year), names_from = Term_name, values_from = c(raw,mean,sd)) %>% 
  relocate(lat, lon, route, year,
           contains("Baseline"), contains("Feeders"), contains("Eucalyptus"), contains("Human\npopulation"))

marginal_change_map <- marginal_effects_summary %>% 
  filter(!is.na(year)) %>% 
  group_by(route, lon, lat) %>% 
  filter(year == min(year) | year == max(year)) %>% 
  filter(n() > 1) %>% 
  summarise(Min_year = min(year),
            Max_year = max(year),
            across(matches("(mean_)"), \(x) { ((x[which.max(year)]/x[which.min(year)])-1)*100 }, #|(raw_)
                   .names = "pct_change_{col}"),
            .groups = "drop") %>% 
  rename_with(\(x) str_remove(x, "mean_")) %>% 
  select(-contains("Baseline"))

##### MARGINAL MAPS #########################################################################################

calif <- vect(
  map_data("state") %>% 
    filter(region == "california") %>% 
    select(lon = long, lat) %>% 
    as.matrix(), 
  type= "polygon", crs="+proj=longlat")

base_df <- marginal_effects_summary %>% 
  filter(!is.na(mean_Baseline)) %>% 
  select(1:3,6,7)

baseline_map <- rast(nrow = length(unique(base_df$lat)),
                     ncol = length(unique(base_df$lon)),
                     xmin = min(base_df$lon), xmax = max(base_df$lon),
                     ymin = min(base_df$lat), ymax = max(base_df$lat),
                     crs = "EPSG:4326", resolution = 0.2) %>% 
  rasterize(vect(base_df), ., field = names(base_df)[-c(1:2)])

change_map <- rast(nrow = length(unique(marginal_change_map$lat)),
                   ncol = length(unique(marginal_change_map$lon)),
                   xmin = min(marginal_change_map$lon), xmax = max(marginal_change_map$lon),
                   ymin = min(marginal_change_map$lat), ymax = max(marginal_change_map$lat),
                   crs = "EPSG:4326", resolution = 0.2) %>% 
  rasterize(vect(marginal_change_map), ., field = names(marginal_change_map)[-c(2,3)]) %>% 
  select(-contains("Baseline"), -contains("raw")) %>% 
  rename_with(\(x) str_replace_all(x, c("Feeders"="fa","Eucalyptus"="ea","Human\npopulation"="pop","Year"="yr")))

counties <- map_data("county") %>% 
  filter(region == "california") %>% 
  group_by(subregion) %>% 
  group_map(\(x,y) {
    cnty <- x %>% 
      select(lon = long, lat) %>% 
      as.matrix()
    vect(cnty, type= "polygon", crs="+proj=longlat", atts = y)
  }) %>% 
  bind_spat_rows()

##### VISUALIZATION #########################################################################################

##### PERCENT CHANGE IN MARGINAL EFFECTS
######### FIGURE 1A
lt_map <- ggplot()+
  geom_spatvector(data = counties, fill = "#FFFFF2")+
  geom_spatraster(data = filter(change_map, Min_year < 1980, Max_year > 2015),
                  aes(fill = pct_change_yr), na.rm = T)+
  scale_fill_viridis_c(na.value = "transparent", 
                       name = "Change in abundance from\nyear-by-latitude effect",
                       breaks = seq(0,150,50), labels = paste0(seq(0,150,50),"%"), limits = c(0,175))+
  thm_obj+
  # annotation_raster(as.raster(image_read("./yr.png")), -123.25, -121, 33.5, 34.9)+
  annotate("text", x = -124, y = 33, label = "bold(A)", parse = T, size = 8)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.725, 0.7),
        axis.text.y = element_text(size = 13),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(t=0, r=0, b=0, l=0, unit = "mm"))

######### FIGURE 1B
hp_map <- ggplot()+
  geom_spatvector(data = counties, fill = "#FFFFF2")+
  geom_spatraster(data = filter(change_map, Min_year < 1980, Max_year > 2015),
                  aes(fill = pct_change_pop), na.rm = T)+
  scale_fill_viridis_c(na.value = "transparent", 
                       name = "Change in abundance from\nchange in human population",
                       breaks = seq(0,60,15), labels = paste0(seq(0,60,15),"%"), limits = c(0,60))+
  thm_obj+
  annotate("text", x = -124, y = 33, label = "bold(B)", parse = T, size = 8)+
  # annotation_raster(as.raster(image_read("./hp.png")), -123.25, -121.25, 33.35, 35.25)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.725, 0.7),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(t=0, r=0, b=0, l=0, unit = "mm"))

######### FIGURE 1C
fa_map <- ggplot()+
  geom_spatvector(data = counties, fill = "#FFFFF2")+
  geom_spatraster(data = filter(change_map, Min_year < 1980, Max_year > 2015),
                  aes(fill = pct_change_fa), na.rm = T)+
  scale_fill_viridis_c(na.value = "transparent", 
                       name = "Change in abundance from\nchange in feeder availability",
                       breaks = seq(0,125,25), labels = paste0(seq(0,125,25),"%"), limits = c(0,125))+
  thm_obj+
  annotate("text", x = -124, y = 33, label = "bold(C)", parse = T, size = 8)+
  # annotation_raster(as.raster(image_read("./fa.png")), -123.25, -121.25, 33.25, 35.25)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.725, 0.7),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13, hjust = 2/3),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(t=0, r=2, b=0, l=-10, unit = "mm"))

######### FIGURE 1D
ea_map <- ggplot()+
  geom_spatvector(data = counties, fill = "#FFFFF2")+
  geom_spatraster(data = filter(change_map, Min_year < 1980, Max_year > 2015),
                  aes(fill = pct_change_ea), na.rm = T)+
  scale_fill_viridis_c(na.value = "transparent", 
                       name = "\nChange in abundance from\nchange in Eucalytpus\navailability",
                       breaks = seq(0,1.5,0.3), labels = paste0(seq(0,1.5,0.3),"%"), limits = c(0,1.5))+
  thm_obj+
  # annotation_raster(as.raster(image_read("./ea.png")), -123.25, -121.25, 33.25, 34.75)+
  annotate("text", x = -124, y = 33, label = "bold(D)", parse = T, size = 8)+
  theme(legend.position = "inside",
        legend.position.inside = c(0.725, 0.7),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 13, hjust = 2/3),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(t=0, r=2, b=0, l=-10, unit = "mm"))

######### FIGURE 1
ggarrange(plotlist = list(lt_map,hp_map,fa_map,ea_map), ncol = 2, nrow = 2, align = "hv")

######### FIGURE S9
# check main effects, excluding intercept, dispersion, smoothing parameters
ggplot(filter(main_result_summary, !str_detect(Term, "(sd)|(_in)|(_ef)|(Rn)|(bayes)")),
       aes(x = Term_name, y = mean, colour = signif))+
  scale_y_continuous(expand = c(0,0.1))+
  geom_hline(yintercept = 0, lty = 2, colour = "darkgrey")+
  geom_point(size = 2)+
  # crossbars correspond to 90% credible intervals
  geom_linerange(aes(ymin = `5%`, ymax = `95%`), linewidth = 1)+
  scale_colour_manual(values = c("#a50026","#d73027","black","#4575b4","#313695"), drop = F)+
  guides(colour = "none")+
  labs(y = "Effect sizes on expected counts\nof Anna's hummingbirds", x = "", fill = "")+
  coord_flip()+
  thm_obj+
  theme(axis.text = element_text(face = "bold"),
        axis.title = element_text(size = 13, face = "bold"))

######### FIGURE S10
ggplot()+
  geom_spatvector(data = counties, fill = "#FFFFF2")+
  geom_spatraster(data = baseline_map, aes(fill = mean_Baseline), na.rm = T)+
  scale_fill_viridis_c(na.value = "transparent", name = "Expected counts of\nAnna's hummingbirds",
                       breaks = seq(0,125,25), labels = seq(0,125,25), limits = c(0,125))+
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.7),
        legend.title = element_text(vjust = 2.5, hjust = 1))+
  thm_obj
