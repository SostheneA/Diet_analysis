## =============================================================================
## DATA CLEANING AND ADDING SOME VARIABLES WORKFLOW
## =============================================================================

# 1. Read data and correct size info
# --- a) Load libraries and data -----------------------------------------------
rm(list = ls())

library(data.table)
library(gulf)
library(tidyverse)
library(sf)
library(readr)

load("data/prey_groups.RData")
#load("data/prey_groups_fixed.RData")
load("data/diet_corr.RData")

# --- b) Data Integration & Variable Engineering -------------------------------

# Merge predator diet data with prey functional groups
prey_groups_unique <- unique(prey_groups, by = "prey_species_common_name")
#prey_groups_unique <- unique(prey_groups_fixed, by = "prey_species_common_name")
diet_dt <- merge(
  diet_corr[, prey_species_latin_name := NULL],
  prey_groups_unique,
  by = "prey_species_common_name",
  all.x = TRUE
)

Database <- as.data.table(diet_dt)

# Standardize formats and handle missing values
Database <- Database %>%
  mutate(
    somatic_length_cm = as.numeric(somatic_length_cm),
    number_of_prey    = replace_na(as.numeric(number_of_prey), 0),
    somatic_wt_g      = replace_na(as.numeric(somatic_wt_g), 0),
    year              = cruise_year
  ) %>%
  # Define study periods for temporal comparative analysis
  mutate(period = case_when(
    cruise_year %in% 2004:2006 ~ "2004-2006",
    cruise_year %in% 2018:2019 ~ "2018-2019",
    TRUE ~ NA_character_
  ))

Database[, period := factor(period, levels = c("2004-2006", "2018-2019"))]

# --- c) Biological Constraint Auditing ----------------------------------------

# Load biological reference table (Lmax per species)
#lmat_predators <- read_csv("data/lmat_predators.csv") # LL20260715 fishbase in source column I think it is not enough informative, maybe a geographic information column
lmat_predators <- fread("data/lmat_predators.csv", encoding = "Latin-1")
setDT(lmat_predators)
lmat_predators[lmat == "Unknown" | is.na(lmat), lmat := as.character(lmax)]
names(lmat_predators)[1]="predator_species_code"

refs <- lmat_predators %>% select(predator_species_code, lmax, lmat)

# Identify biological anomalies (Observations > Theoretical Max)
errors <- Database %>%
  left_join(refs, by = "predator_species_code") %>%
  filter(somatic_length_cm > lmax) %>%
  select(stomach_id, predator_species_code, predator_species_common_name,
         somatic_length_cm, lmax)

# --- d) Surgical Data Correction (Double-Threshold Logic) --------------------

# Create a clean working copy and join reference lengths
Database_clean <- Database %>% left_join(refs, by = "predator_species_code")
setDT(Database_clean)

# STRATEGY 1: Unit Correction (Millimeters to Centimeters)
# If length is > 2x the biological max, it is likely a unit entry error (mm) # LL20260715 it has a problem with smooth skate lmax
Database_clean[somatic_length_cm > (2 * lmax),
               somatic_length_cm_p := somatic_length_cm / 10]

# STRATEGY 2: Capping Protocol
# If length still exceeds Lmax, cap the value at Lmax to handle rounding/extremes
Database_clean[somatic_length_cm > lmax,
               somatic_length_cm := lmax]

# Update the main Database with corrected values
Database$somatic_length_cm <- Database_clean$somatic_length_cm

# Verification check: Ensure zero errors remain
errors_test <- Database_clean %>%
  filter(somatic_length_cm > lmax) %>%
  select(stomach_id, predator_species_code,predator_species_common_name, somatic_length_cm, lmax)


Database = merge(Database, refs, by="predator_species_code")

##########


# 2. Access Fall sGLS trawl survey 'set' data via gulf package
sets <- read.card(
  survey = "rv",
  sampling = "research",
  card.type = "set",
  year = c(2004:2006, 2018:2019),
  as.data.frame = TRUE
) %>% as.data.table()

# 2. Variable Engineering for spatial analysis
sets[, ":="(
  samp.date = ymd(paste(year, month, day, sep = "-")),
  latitude  = (latitude.start + latitude.end) / 2,
  longitude = (longitude.start + longitude.end) / 2,
  depth     = (depth.start + depth.end) / 2
)]

# 3. Merge metadata with biological data
Database <- merge(
  x = Database,
  y = sets[, .(stratum, year, month, day, samp.date, set = set.number,
               vessel.code, hour = start.hour, latitude, longitude, depth,
               bottom_temp = bottom.temperature, salinity = bottom.salinity)],
  by = c("year", "set", "vessel.code"),
  all.x = TRUE
) #%>% filter(!(is.na(stratum) & is.na(latitude)))

Database$latitude[Database$samp.date=="2006-09-05" & Database$hour==19]=mean(Database$latitude[Database$samp.date=="2006-09-05" & Database$hour==19], na.rm = T)
Database$longitude[Database$samp.date=="2006-09-05" & Database$hour==19]=mean(Database$longitude[Database$samp.date=="2006-09-05" & Database$hour==19], na.rm = T)

# ==============================================================================
# SCRIPT: Spatial Assignment of Stomach Content Data
# ==============================================================================

library(sf)
library(dplyr)
library(tidyr)

# --- 1. Data Preparation ---
# Convert raw database to a spatial (sf) object
# We use CRS 4326 (WGS84) for initial coordinate handling
Database <- Database %>%
  mutate(Lat = latitude, Lon = longitude) %>%
  drop_na(Lon, Lat) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)


# Ensure geometries are valid before processing
EAR_sf <- st_read("data/Spatial_data/EAR_map/EAR_map.shp")
NAFO_4T_sf <- st_read("data/Spatial_data/NAFO/nafo_2014_02.shp")
NAFO_4A_sf <- st_read("data/Spatial_data/nafo_4T/nafo_2014_02.shp")

new_ecoregions_clean <- st_read("data/Spatial_data/new_ecoregions_final.shp")
new_regions_v2_clean <- st_read("data/Spatial_data/new_regions_v2_final.shp")
new_gulf_clean <- st_read("data/Spatial_data/new_gulf_final.shp")

# --- 2. CRS Harmonization ---
# Ensure all spatial layers match the specimen data projection
target_crs           <- st_crs(Database)
new_ecoregions_clean <- st_transform(new_ecoregions_clean, target_crs)
new_regions_v2_clean <- st_transform(new_regions_v2_clean, target_crs)
NAFO_4T_sf           <- st_transform(NAFO_4T_sf, target_crs)
EAR_sf               <- st_transform(EAR_sf, target_crs)
EAR_clean_sf         <- st_transform(new_gulf_clean, target_crs)

# --- 3. Spatial Attribution (Nearest Neighbor Method) ---
# We use st_nearest_feature to avoid row duplication and handle
# points located slightly offshore or on polygon boundaries.
message("Assigning spatial attributes...")

Database$EcoZone <- new_ecoregions_clean$EcoZone[st_nearest_feature(Database, new_ecoregions_clean)]
Database$Zone    <- new_regions_v2_clean$Zone[st_nearest_feature(Database, new_regions_v2_clean)]
Database$NAFO    <- NAFO_4T_sf$level_2[st_nearest_feature(Database, NAFO_4T_sf)]
Database$Region  <- EAR_sf$regn_nm[st_nearest_feature(Database, EAR_sf)]
Database$Area  <- EAR_clean_sf$Area[st_nearest_feature(Database, EAR_clean_sf)]

# --- 4. Quality Control & Handling Edge Cases ---
# Fill missing spatial values (NAs) for points outside the specific polygon extents
Database <- Database %>%
  mutate(
    EcoZone = replace_na(EcoZone, "4TF Only"),
    Zone    = replace_na(Zone, "4TF Only"),
    NAFO    = replace_na(NAFO, "4TF") # Specific label for outside NAFO 4T polygons
  )

Database <- Database %>%
  mutate(
    is_nutritional_prey = if_else(is_empty == "FALSE", 1, 0)
  )

save(Database, file = "data/Database.rda")

if (!dir.exists("output/Figures")) {
  dir.create("output/Figures", recursive = TRUE)
}

#######################

### Visualisation

# Prepare "All Years" summary data
diet_all_years <- Database %>% mutate(year = "All Years")
diet_plot_data <- bind_rows(Database %>% mutate(year = as.character(year)),
                            diet_all_years)

# Order facets so "All Years" is at the end
diet_plot_data$year <- factor(diet_plot_data$year,
                              levels = c(sort(unique(as.character(diet_dt$year))), "All Years"))

library(rnaturalearth)
library(rnaturalearthdata)

# Load country boundaries
canada <- ne_states(country = "canada", returnclass = "sf")
usa    <- ne_states(country = "united states of america", returnclass = "sf")

EcoZone <- ggplot() +
  geom_sf(data = canada, fill = "gray90", color = "white") +
  geom_sf(data = usa, fill = "gray95", color = "gray80") +
  geom_point(data = diet_plot_data,
             aes(x = longitude, y = latitude, color = EcoZone),
             size = 0.8, alpha = 0.6) +
  coord_sf(xlim = c(-67, -59.5), ylim = c(45, 49.5), expand = FALSE) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "EcoRegions Map: Gulf of St. Lawrence",
       subtitle = "Time series and Global Distribution of Diet Samples",
       x = "Longitude", y = "Latitude") +
  facet_wrap(~year, ncol = 3)

ggsave("output/Figures/EcoZone_diet.png", EcoZone, width = 12, height = 10)



Zone <- ggplot() +
  geom_sf(data = canada, fill = "gray90", color = "white") +
  geom_sf(data = usa, fill = "gray95", color = "gray80") +
  geom_point(data = diet_plot_data,
             aes(x = longitude, y = latitude, color = Zone),
             size = 0.8, alpha = 0.6) +
  coord_sf(xlim = c(-67, -59.5), ylim = c(45, 49.5), expand = FALSE) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Regions Map: Gulf of St. Lawrence",
       subtitle = "Time series and Global Distribution of Diet Samples",
       x = "Longitude", y = "Latitude") +
  facet_wrap(~year, ncol = 3)

ggsave("output/Figures/Zone_diet.png", Zone, width = 12, height = 10)


NAFO <- ggplot() +
  geom_sf(data = canada, fill = "gray90", color = "white") +
  geom_sf(data = usa, fill = "gray95", color = "gray80") +
  geom_point(data = diet_plot_data,
             aes(x = longitude, y = latitude, color = NAFO),
             size = 0.8, alpha = 0.6) +
  coord_sf(xlim = c(-67, -59.5), ylim = c(45, 49.5), expand = FALSE) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "EcoRegions Map: Gulf of St. Lawrence",
       subtitle = "Time series and Global Distribution of Diet Samples",
       x = "Longitude", y = "Latitude") +
  facet_wrap(~year, ncol = 3)

ggsave("output/Figures/NAFO_diet.png",NAFO, width = 12, height = 10)


# --- 4. Visualization: Spatiotemporal Distribution ---

# Re-enable S2 for plotting
sf_use_s2(TRUE)

# Define base map aesthetic to avoid repetition
base_map <- ggplot() +
  geom_sf(data = canada, fill = "gray90", color = "white") +
  geom_sf(data = new_gulf_clean, fill = NA, color = "black", size = 0.5) +
  geom_sf(data = usa, fill = "gray95", color = "gray80") +
  coord_sf(xlim = c(-67, -59.5), ylim = c(45, 49.5), expand = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Longitude", y = "Latitude", color = "Sampling Area")

# Plot A: Annual Distribution
map_yearly <- base_map +
  geom_point(data = diet_plot_data,
             aes(x = longitude, y = latitude, color = Area),
             size = 0.6, alpha = 0.5) +
  facet_wrap(~year, ncol = 3) +
  labs(title = "Fish Diet Sampling Distribution: Annual",
       subtitle = "Spatial coverage within defined Gulf ecoregions")

# Plot B: Period-based Distribution (Comparison across eras)
map_period <- base_map +
  geom_point(data = diet_plot_data,
             aes(x = longitude, y = latitude, color = Area),
             size = 0.8, alpha = 0.6) +
  facet_wrap(~period, ncol = 2) +
  labs(title = "Fish Diet Sampling Distribution: Historical vs Modern",
       subtitle = "Comparative spatial coverage across sampling periods")

# Display maps
print(map_yearly)
print(map_period)

ggsave("output/Figures/map_yearly_diet.png",NAFO, width = 12, height = 10)
ggsave("output/Figures/map_period_diet.png",NAFO, width = 12, height = 10)


######################

library(ggplot2)
library(patchwork) # Pour mettre les cartes côte à côte
setDT(Database)
# 1. Filtrage des données pour nos deux groupes de "Vers"
map_data <- Database

# 1. Préparation de la liste des taxons
unique_taxa <- sort(unique(map_data$prey_category_new))
# On divise les 36 taxons en 4 listes de 9
taxa_groups <- split(unique_taxa, ceiling(seq_along(unique_taxa) / 9))

# 2. Création de la fonction de plot pour éviter la répétition
make_diet_map <- function(taxa_list, plot_number) {

  ggplot(map_data[prey_category_new %in% taxa_list], aes(x = longitude, y = latitude)) +
    stat_summary_hex(aes(z = n_prey_count_per_prey), fun = sum, bins = 25) +
    scale_fill_viridis_c(option = "viridis", name = "Abondance") +
    # Contrainte 3x3 ici
    facet_wrap(~prey_category_new, ncol = 3, nrow = 3) +
    coord_quickmap() +
    theme_minimal() +
    labs(
      title = paste0("Distribution Spatiale - Groupe ", plot_number, "/4"),
      subtitle = "Grille 3x3 (Somme de l'abondance)",
      x = "Longitude", y = "Latitude"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 7), # Texte plus petit pour que ça rentre
      legend.position = "right",
      panel.spacing = unit(0.5, "lines")
    )
}

# 3. Génération et affichage/sauvegarde
# On utilise une boucle pour créer les 4 objets
for (i in 1:4) {
  p <- make_diet_map(taxa_groups[[i]], i)
  print(p)
  ggsave(paste0("output/Figures/map_taxa_part_", i, ".png"), p, width = 12, height = 10)
}



# Draw strata for the rv survey:
file_path <- "output/Figures/map_strata_sGSL.png"
png(filename = file_path, width = 2000, height = 1800, res = 300)

library(gulf)
gulf.map(xlim = c(-66.2, -60), ylim = c(45.5, 49.2), sea = TRUE)
map.strata(survey = "rv", labels = TRUE,region = "gulf")
dev.off()


#####################################  CREATED CLEAN DATA WITH EXISTING PREYS vs Vacuity index dataset

invalid_reasons <- c("content weight larger than stomach", "No set number", "Prey weights too large")
Database$ID= 1:nrow(Database)

diet_clean <-  Database%>%
  group_by(predator_species_common_name, predator_species_latin_name) %>%
  mutate(
    Total_N = n_distinct(stomach_id),
    # Calcul du N par période
    N_2004_2006 = n_distinct(stomach_id[period == "2004-2006"]),
    N_2018_2019 = n_distinct(stomach_id[period == "2018-2019"]),

    # Application de la double règle
    Status = if_else(
      Total_N >= 90 & N_2004_2006 >= 10 & N_2018_2019 >= 10,
      "Retained",
      "Excluded (low N)"
    )
  ) %>%
  filter(Status == "Retained") %>%
  filter(!predator_invalid_reason %in% invalid_reasons) %>%
  filter(!(prey_category_keep  %in% c("parasite", "empty",  "digested")))


diet_sub <-  subset(Database, !(ID %in% c(diet_clean$ID)))


############ Fullness data
diet_C <-  Database %>%
  group_by(predator_species_common_name, predator_species_latin_name) %>%
  mutate(
    Total_N = n_distinct(stomach_id),
    # Calcul du N par période
    N_2004_2006 = n_distinct(stomach_id[period == "2004-2006"]),
    N_2018_2019 = n_distinct(stomach_id[period == "2018-2019"]),

    # Application de la double règle
    Status = if_else(
      Total_N >= 90 & N_2004_2006 >= 10 & N_2018_2019 >= 10,
      "Retained",
      "Excluded (low N)"
    )
  ) %>%
  filter(Status == "Retained") %>%
  filter(!predator_invalid_reason %in% invalid_reasons)


stomach_fullness <- diet_C %>%
  dplyr::group_by(stomach_id,period,  year,  predator_species_common_name) %>%
  dplyr::summarise(
    has_valid_food = if_else(sum(is_nutritional_prey, na.rm = TRUE) > 0, 1, 0),
    .groups = "drop"
  )

save(diet_C, stomach_fullness, file="data/fullness.rda")

########################################

library(data.table)
library(gulf)

# 1. Passer en data.table
setDT(diet_clean)

# 2. Créer un tableau des combinaisons UNIQUES (Espèce, Taille, Année)
# C'est ce tableau qui servira de dictionnaire pour interroger Oracle.
# On filtre les NA pour ne pas envoyer de requêtes inutiles.
unique_samples <- diet_clean[!is.na(predator_species_code) &
                               !is.na(somatic_length_cm) &
                               !is.na(year),
                             .(n_fish = .N),
                             by = .(predator_species_code, somatic_length_cm, year)]

# 3. Calculer le poids sur ce petit tableau (Très efficace)
# On utilise une boucle par ligne sur ce tableau RÉDUIT
unique_samples[, pred_weight_est := {
  res <- tryCatch({
    # Appel direct à la fonction gulf
    w <- weight(somatic_length_cm, species = predator_species_code, year = year)*1000
    if(length(w) > 0) w[1] else NA_real_
  }, error = function(e) NA_real_)
  res
}, by = 1:nrow(unique_samples)]


### uSE 2017 info for herring in 2018 to fill its NA values

unique_samples[year == 2018 & predator_species_code == 60,
               pred_weight_est := weight(somatic_length_cm,
                                         species = predator_species_code,
                                         year = 2017) * 1000]

unique_samples[year == 2019 & predator_species_code == 60,
               pred_weight_est := weight(somatic_length_cm,
                                         species = predator_species_code,
                                         year = 2017) * 1000]


save(unique_samples, file="data/unique_samples.rda")

# 4. Joindre le résultat au gros dataset (Update-on-join)
# R va "mapper" instantanément les poids calculés sur vos milliers de lignes.
diet_clean[unique_samples, on = .(predator_species_code, somatic_length_cm, year),
           pred_weight_est := i.pred_weight_est]

# 1. Calculer le minimum des valeurs strictement positives
min_val <- diet_clean[somatic_wt_g > 0, min(somatic_wt_g, na.rm = TRUE)]

# 2. Remplacer les 0 par cette valeur
diet_clean[somatic_wt_g == 0, somatic_wt_g := min_val]

# 5. Calcul final de ton indice de réplétion (PFI)
diet_clean[, pfi := (somatic_wt_g / pred_weight_est) * 100]


save(diet_clean, file="data/diet_clean.rda")
