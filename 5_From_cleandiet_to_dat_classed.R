# =============================================================================
# STANDARDIZED TROPHIC ANALYSIS BY SIZE CLASS × AREA × PERIOD
# Design: Global Effect → Size Effect → Final Stratified Test
#
# Key Features:
#   • Rigorous statistical tests for ΔH' (Hutcheson + Bootstrap CI)
#   • Bootstrap CI for ΔBs (Levins Niche Breadth)
#   • Bray-Curtis + PERMANOVA (adonis2) by stratum
#   • SIMPER analysis to identify drivers of community shifts
#   • Decision matrix based on p-values rather than subjective thresholds
#   • Corrected N_MIN logic for strata inventory
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# 0. PACKAGE INITIALIZATION
# ─────────────────────────────────────────────────────────────────────────────
rm(list = ls()) # Clear workspace
packages <- c("tidyverse", "vegan", "DT", "plotly", "kableExtra",
              "RColorBrewer", "htmltools", "data.table")

# Install missing packages if necessary
new_pkg <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(new_pkg)) install.packages(new_pkg, dependencies = TRUE)

library(tidyverse)
library(vegan)
library(DT)
library(plotly)
library(kableExtra)
library(RColorBrewer)
library(data.table)

# Detach igraph if present to avoid conflicts with vegan::diversity
if ("igraph" %in% (.packages())) detach("package:igraph", unload = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 1. DATA LOADING & DYNAMIC COLUMN SELECTION
# ─────────────────────────────────────────────────────────────────────────────
load("data/diet_clean.rda")

# Define core metadata columns to retain
cols_to_keep <- c("predator_species_common_name", "period", "Area", "stomach_id",
                  "somatic_length_cm", "somatic_wt_g", "longitude", "latitude",
                  "depth", "bottom_temp", "salinity")

# Dynamically identify all prey category columns (Sensitivity X variables)
cols_prey_x <- names(diet_clean)[grep("prey_category_", names(diet_clean))]

# Create working dataset using data.table for performance
data_doc <- diet_clean[, c(cols_to_keep, cols_prey_x), with = FALSE]

# Adjust prey counts: Ensure 0 values are treated as 1 for presence-based logic
if("number_of_prey" %in% names(data_doc)) {
  data_doc$number_of_prey[data_doc$number_of_prey == 0] <- 1
}

cat("Dimensions:", nrow(data_doc), "rows ×", ncol(data_doc), "columns\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2. GLOBAL PARAMETERS & PRE-PROCESSING
# ─────────────────────────────────────────────────────────────────────────────
N_STRICT <- 25   # Threshold for full analysis without reliability warnings
ALPHA    <- 0.05 # Statistical significance level

# Visualization Colors
col_period <- c("2004-2006" = "#2c7bb6", "2018-2019" = "#d7191c")
col_Area   <- c("Northumberland"    = "#1b7837",
                "Magdalen Shallows" = "#762a83",
                "Central"           = "#e08214",
                "Chaleur Bay" = "#2166ac")

# Data Cleaning: Filter invalid lengths and factorize key variables
dat <- data_doc %>%
  filter(!is.na(somatic_length_cm), somatic_length_cm > 0) %>%
  mutate(
    period = factor(period, levels = c("2004-2006", "2018-2019")),
    Area   = factor(Area)
  )


# ─────────────────────────────────────────────────────────────────────────────
# BLOC 1 — SIZE CLASS ATTRIBUTION
# ─────────────────────────────────────────────────────────────────────────────

library(data.table)
setDT(diet_clean)

# 1. S'assurer que les colonnes sont numériques
diet_clean[, `:=`(lmat = as.numeric(lmat), lmax = as.numeric(lmax))]

diet_clean[, size_class := fcase(
  lmat == lmax, "one_class",
  somatic_length_cm < lmat, "juvenile",
  somatic_length_cm >= lmat, "adult"
)]

table(diet_clean$predator_species_common_name, diet_clean$size_class)
table(diet_clean$predator_species_common_name, diet_clean$size_class, diet_clean$period)

# Le margin c(1, 3) signifie : calcule les proportions sur l'index 2 (size_class)
# pour chaque combinaison de l'index 1 (espèce) et 3 (période).
tab_pct <- prop.table(table(diet_clean$predator_species_common_name,
                            diet_clean$size_class,
                            diet_clean$period), margin = c(1, 3)) * 100

# Vérification pour la période
P1=data.frame(tab_pct[,, "2004-2006"])
P2=data.frame(tab_pct[,, "2018-2019"])


library(dplyr)
library(sf)
stats_size <- diet_clean %>%
  st_drop_geometry() %>%
  group_by(predator_species_common_name, period) %>%
  count(size_class) %>%
  mutate(pct = (n / sum(n)) * 100) %>%
  ungroup()

# Si tu veux retrouver le format "Wide" (colonnes adult, juvenile, one_class)
stats_wide <- stats_size %>%
  select(-n) %>%
  pivot_wider(names_from = size_class, values_from = pct, values_fill = 0)

print(stats_wide)



dat_classed= diet_clean
save(dat_classed, file='data/dat_classed.rda')
