# =============================================================================
# 4_Database_with_others_factors.R - DIET RECORDS x PREY GROUPS x SURVEY SETS
# -----------------------------------------------------------------------------
# Input  : data/diet_corr.RData, data/prey_groups.RData,
#          data/prey_groups_PP.RData (optional, from 3PP_taxonomic_groups.R),
#          data/lmat_predators.csv, data/Spatial_data/*.shp,
#          survey set cards (gulf::read.card), data/unique_samples.rda (cache)
# Output : data/Database.rda        all prey records with covariates and Area
#          data/correction_log_lengths.csv   fish whose length was corrected
#          data/fullness.rda        retained predators, all stomachs (vacuity)
#          data/diet_clean.rda      retained predators, nutritional prey only,
#                                   with estimated predator weight and PFI
# =============================================================================
rm(list = ls())

library(data.table)
library(gulf)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)

load("data/prey_groups.RData")
load("data/diet_corr.RData")

# 1) Prey groups onto diet records -----------------------------------------------
prey_groups_unique <- unique(prey_groups, by = "prey_species_common_name")
Database <- merge(diet_corr[, prey_species_latin_name := NULL], prey_groups_unique,
                  by = "prey_species_common_name", all.x = TRUE)

# Per-predator family (_PP), keyed on (predator, prey); an all.x merge on a
# unique key must leave the row count unchanged.
if (file.exists("data/prey_groups_PP.RData")) {
  load("data/prey_groups_PP.RData")
  pp_cols <- grep("^(prey_category|tax_level)_\\d+_PP$", names(prey_groups_PP), value = TRUE)
  pp_map  <- unique(prey_groups_PP[, c("predator", "prey_species_common_name", pp_cols), with = FALSE])
  if (anyDuplicated(pp_map, by = c("predator", "prey_species_common_name")))
    stop("prey_groups_PP: (predator, prey) key is not unique.")
  n_before <- nrow(Database)
  Database <- merge(Database, pp_map,
                    by.x = c("predator_species_common_name", "prey_species_common_name"),
                    by.y = c("predator", "prey_species_common_name"), all.x = TRUE)
  stopifnot(nrow(Database) == n_before)
  first_pp <- grep("^prey_category_\\d+_PP$", names(Database), value = TRUE)[1]
  n_na <- sum(is.na(Database[[first_pp]]) & !is.na(Database$prey_species_common_name))
  cat("PP family merged:", length(pp_cols) %/% 2, "resolution columns;",
      n_na, "prey rows without PP classification\n")
}

Database[, `:=`(
  somatic_length_cm = as.numeric(somatic_length_cm),
  number_of_prey    = fcoalesce(as.numeric(number_of_prey), 0),
  somatic_wt_g      = fcoalesce(as.numeric(somatic_wt_g), 0),
  year              = cruise_year,
  period            = factor(fcase(cruise_year %in% 2004:2006, "2004-2006",
                                   cruise_year %in% 2018:2019, "2018-2019"),
                             levels = c("2004-2006", "2018-2019"))
)]

# 2) Length correction against Lmax ---------------------------------------------
# Lengths above 2 x Lmax are millimetres entered as centimetres (divided by 10);
# lengths still above Lmax are capped at Lmax.
lmat_predators <- fread("data/lmat_predators.csv", encoding = "Latin-1")
lmat_predators[lmat == "Unknown" | is.na(lmat), lmat := as.character(lmax)]
setnames(lmat_predators, 1, "predator_species_code")
refs <- lmat_predators[, .(predator_species_code, lmax, lmat)]

Database <- merge(Database, refs, by = "predator_species_code", all.x = TRUE)
Database[, somatic_length_cm_orig := somatic_length_cm]
Database[somatic_length_cm > 2 * lmax, somatic_length_cm := somatic_length_cm / 10]
Database[somatic_length_cm > lmax,     somatic_length_cm := lmax]

correction_log <- unique(
  Database[somatic_length_cm != somatic_length_cm_orig,
           .(stomach_id, predator_species_code, predator_species_common_name,
             original = somatic_length_cm_orig, corrected = somatic_length_cm, lmax,
             action = fcase(
               somatic_length_cm_orig > 2 * lmax & somatic_length_cm == somatic_length_cm_orig / 10, "unit mm->cm",
               somatic_length_cm_orig > 2 * lmax & somatic_length_cm == lmax, "unit mm->cm then cap",
               somatic_length_cm == lmax, "cap at Lmax",
               default = "other"))],
  by = "stomach_id")
print(correction_log[, .N, by = action])
fwrite(correction_log, "data/correction_log_lengths.csv")
Database[, somatic_length_cm_orig := NULL]

# 3) Survey set covariates ----------------------------------------------------------
sets <- as.data.table(read.card(survey = "rv", sampling = "research", card.type = "set",
                                year = c(2004:2006, 2018:2019), as.data.frame = TRUE))
sets[, `:=`(
  samp.date = ymd(paste(year, month, day, sep = "-")),
  latitude  = (latitude.start + latitude.end) / 2,
  longitude = (longitude.start + longitude.end) / 2,
  depth     = (depth.start + depth.end) / 2
)]

Database <- merge(
  Database,
  sets[, .(stratum, year, month, day, samp.date, set = set.number, vessel.code,
           hour = start.hour, latitude, longitude, depth,
           bottom_temp = bottom.temperature, salinity = bottom.salinity)],
  by = c("year", "set", "vessel.code"), all.x = TRUE)

# One set (2006-09-05, 19 h) has a missing position: use the mean of its records.
fix <- Database$samp.date == "2006-09-05" & Database$hour == 19
Database$latitude[fix]  <- mean(Database$latitude[fix],  na.rm = TRUE)
Database$longitude[fix] <- mean(Database$longitude[fix], na.rm = TRUE)

# 4) Spatial attribution (nearest polygon, planar) -----------------------------------
Database <- Database %>%
  mutate(Lat = latitude, Lon = longitude) %>%
  drop_na(Lon, Lat) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)

read_layer <- function(path) st_make_valid(st_transform(st_read(path, quiet = TRUE), st_crs(Database)))
new_ecoregions_clean <- read_layer("data/Spatial_data/new_ecoregions_final.shp")
new_regions_v2_clean <- read_layer("data/Spatial_data/new_regions_v2_final.shp")
new_gulf_clean       <- read_layer("data/Spatial_data/new_gulf_final.shp")
NAFO_4T_sf           <- read_layer("data/Spatial_data/NAFO/nafo_2014_02.shp")
EAR_sf               <- read_layer("data/Spatial_data/EAR_map/EAR_map.shp")

sf_use_s2(FALSE)
Database$EcoZone <- new_ecoregions_clean$EcoZone[st_nearest_feature(Database, new_ecoregions_clean)]
Database$Zone    <- new_regions_v2_clean$Zone[st_nearest_feature(Database, new_regions_v2_clean)]
Database$NAFO    <- NAFO_4T_sf$level_2[st_nearest_feature(Database, NAFO_4T_sf)]
Database$Region  <- EAR_sf$regn_nm[st_nearest_feature(Database, EAR_sf)]
Database$Area    <- new_gulf_clean$Area[st_nearest_feature(Database, new_gulf_clean)]
sf_use_s2(TRUE)

Database <- Database %>%
  mutate(EcoZone = replace_na(EcoZone, "4TF Only"),
         Zone    = replace_na(Zone, "4TF Only"),
         NAFO    = replace_na(NAFO, "4TF"),
         is_nutritional_prey = if_else(is_empty == FALSE, 1, 0))

save(Database, file = "data/Database.rda")

# 5) Retained predators ------------------------------------------------------------
# >= 90 stomachs over the study and >= 10 in each period; invalid records out.
invalid_reasons <- c("content weight larger than stomach", "No set number", "Prey weights too large")

retained <- Database %>%
  group_by(predator_species_common_name, predator_species_latin_name) %>%
  mutate(Total_N     = n_distinct(stomach_id),
         N_2004_2006 = n_distinct(stomach_id[period == "2004-2006"]),
         N_2018_2019 = n_distinct(stomach_id[period == "2018-2019"])) %>%
  ungroup() %>%
  filter(Total_N >= 90, N_2004_2006 >= 10, N_2018_2019 >= 10,
         !predator_invalid_reason %in% invalid_reasons)

# Vacuity: every stomach of the retained predators.
diet_C <- retained
stomach_fullness <- diet_C %>%
  st_drop_geometry() %>%
  group_by(stomach_id, period, year, predator_species_common_name) %>%
  summarise(has_valid_food = if_else(sum(is_nutritional_prey, na.rm = TRUE) > 0, 1, 0),
            .groups = "drop")
save(diet_C, stomach_fullness, file = "data/fullness.rda")

# Diet: nutritional prey only.
diet_clean <- retained %>%
  filter(!prey_category_keep %in% c("parasite", "empty", "digested"))
setDT(diet_clean)

# 6) Predator weight (length-weight, gulf::weight) and fullness index --------------
# Cached per (species, length, year); R_helpers/Update_unique_samples.R completes
# the cache after a length correction adds new combinations.
combos <- diet_clean[!is.na(predator_species_code) & !is.na(somatic_length_cm) & !is.na(year),
                     .(n_fish = .N), by = .(predator_species_code, somatic_length_cm, year)]

if (file.exists("data/unique_samples.rda")) {
  load("data/unique_samples.rda")
  setDT(unique_samples)
  n_missing <- nrow(combos[!unique_samples, on = .(predator_species_code, somatic_length_cm, year)])
  cat("unique_samples loaded from cache:", nrow(unique_samples), "combinations\n")
  if (n_missing > 0)
    warning(n_missing, " combinations missing from the cache: run R_helpers/Update_unique_samples.R",
            call. = FALSE)
} else {
  unique_samples <- copy(combos)
  unique_samples[, pred_weight_est := tryCatch({
    w <- weight(somatic_length_cm, species = predator_species_code, year = year) * 1000
    if (length(w) > 0) w[1] else NA_real_
  }, error = function(e) NA_real_), by = seq_len(nrow(unique_samples))]
  # Herring (code 60) has no relation for 2018/2019: use 2017.
  unique_samples[year %in% c(2018, 2019) & predator_species_code == 60,
                 pred_weight_est := weight(somatic_length_cm, species = predator_species_code,
                                           year = 2017) * 1000]
  save(unique_samples, file = "data/unique_samples.rda")
  cat("unique_samples computed:", nrow(unique_samples), "combinations saved\n")
}

diet_clean[unique_samples, on = .(predator_species_code, somatic_length_cm, year),
           pred_weight_est := i.pred_weight_est]

min_wt <- diet_clean[somatic_wt_g > 0, min(somatic_wt_g, na.rm = TRUE)]
diet_clean[somatic_wt_g == 0, somatic_wt_g := min_wt]
diet_clean[, pfi := (somatic_wt_g / pred_weight_est) * 100]

save(diet_clean, file = "data/diet_clean.rda")
