# =============================================================================
# 8_Tables.R - MANUSCRIPT TABLES 1-3 AND S1-S6 (csv)
# -----------------------------------------------------------------------------
# Run after 6b, 6c and 6d. Same aggregation rule as the figures
# (freq_by_resolution, Inconclusive excluded), so text, tables and figures
# cannot disagree.
#
#   Table1_typology      the nine states, their signal combination, family
#   Table2_headline      families by contrast and currency, Gulf-wide
#   Table3_crossscale    families by spatial level
#   TableS2_sampling     stomachs and sets per predator, size class and period
#   TableS3_states_gulf  the nine states, Gulf-wide
#   TableS4_by_ecoregion / TableS5_by_stratum   families per unit
#   TableS6_coverage     cells, testable, reliable, confounded dispersion
# Percentages in Tables 2 to S4 are computed on testable cells only; read
# Table S6 alongside them.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
})

PREY_FAMILY <- "_1"
source("R_helpers/Config_Mappings.R")

emit <- function(df, stem, caption, digits = 1) {
  write_csv(df, file.path(DIR_TABLES, paste0(stem, ".csv")))
  cat("\n\n===== ", caption, " =====\n", sep = "")
  print(as.data.frame(df), row.names = FALSE, digits = 4)
  cat("  -> ", stem, ".csv\n", sep = "")
  invisible(df)
}

# =============================================================================
# 1. READ
# =============================================================================
cat("\nReading pipeline results...\n")
res_all <- read_all_runs()

fam_gulf <- res_all %>%
  filter(level == "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast")) %>%
  add_families(keys = c("currency", "contrast"))

fam_unit <- res_all %>%
  filter(level != "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  add_families(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  left_join(
    n_units_by(filter(res_all, level != "all_gulf"),
               keys = c("currency", "contrast", "level", "spatial_unit")),
    by = c("currency", "contrast", "level", "spatial_unit")
  )

# =============================================================================
# TABLE 1 - THE TYPOLOGY
# =============================================================================
readings <- c(
  "Stable Diet"                 = "Same prey, same proportions, same breadth.",
  "Emerging Shift"              = "Early compositional signal, not yet significant.",
  "Ghost Shift"                 = "Prey identities turn over; diversity and breadth hold.",
  "Niche Compression/Expansion" = "Range of prey used changes on a stable prey list.",
  "Internal Rebalancing"        = "Proportions shift within the same prey set.",
  "Niche Restructuring"         = "Diversity and breadth both change without turnover.",
  "Partial Diet Shift"          = "Turnover accompanied by a change in diversity.",
  "Structural Shift"            = "Turnover accompanied by a change in breadth.",
  "Major Shift"                 = "Turnover with change in both diversity and breadth."
)

table1 <- STATE_LOGIC %>%
  mutate(
    family = as.character(family_of(diagnostic)),
    reading = unname(readings[diagnostic])
  ) %>%
  select(Family = family, State = diagnostic,
         Composition = composition, `Shannon H'` = H_prime, `Niche breadth Bs` = Bs,
         `Ecological reading` = reading) %>%
  arrange(factor(Family, levels = FAMILY_LEVELS),
          factor(State, levels = DIAG_LEVELS))

emit(table1, "Table1_typology",
     "Table 1. The nine diagnostic states, their defining signal combination, and the four families they collapse onto.")

# =============================================================================
# TABLE 2 - HEADLINE, GULF-WIDE
# =============================================================================
table2 <- fam_gulf %>%
  filter(!is.na(contrast)) %>%
  arrange(currency, contrast) %>%
  mutate(across(all_of(FAMILY_LEVELS), ~round(.x, 1))) %>%
  transmute(Currency = CURRENCY_LAB[as.character(currency)],
            Contrast = CONTRAST_LAB_1L[as.character(contrast)],
            Stability, Substitution,
            `Functional change`, Reorganisation,
            `Turnover total` = round(Substitution + Reorganisation, 1))

emit(table2, "Table2_headline",
     "Table 2. Frequency (%) of the four families by contrast and currency, Gulf-wide level. Turnover total = Substitution + Reorganisation.")

# =============================================================================
# TABLE 3 - CROSS-SCALE GRADIENT
# =============================================================================
# Gulf-wide is a single pooled test; the other two are means across units. The
# note is part of the table, not a footnote to be lost in copy-editing.
fam_level <- bind_rows(
  fam_gulf %>% mutate(level = "all_gulf", .before = 1),
  fam_unit %>%
    group_by(currency, contrast, level) %>%
    summarise(across(all_of(FAMILY_LEVELS), ~mean(.x, na.rm = TRUE)),
              n_units_total = dplyr::n(), .groups = "drop")
) %>%
  mutate(level = factor(level, levels = SPATIAL_LEVELS_ORD))

table3 <- fam_level %>%
  filter(!is.na(contrast)) %>%
  arrange(currency, contrast, level) %>%
  mutate(across(all_of(FAMILY_LEVELS), ~round(.x, 1))) %>%
  transmute(Currency = CURRENCY_LAB[as.character(currency)],
            Contrast = contrast_label(contrast),
            `Spatial level` = LEVEL_SHORT[as.character(level)],
            `Units` = ifelse(is.na(n_units_total), 1L, n_units_total),
            Stability, Substitution,
            `Functional change`, Reorganisation)

emit(table3, "Table3_crossscale",
     paste0("Table 3. Family frequency (%) by spatial level. Gulf-wide is a single ",
            "pooled test on all trawl sets; ecoregion and stratum values are means ",
            "across units, each unit tested on its own data."))

# =============================================================================
# TABLE S1 - THE PREDATORS
# =============================================================================
# One row per retained predator: common and scientific names, length range and
# number of stomachs with prey in each period, Lmat and Lmax with their source.
# Built from dat_classed (stomachs actually analysed) and data/lmat_predators.csv.
DATA_PATH <- "data/dat_classed.rda"
LMAT_PATH <- "data/lmat_predators.csv"

if (file.exists(DATA_PATH) && file.exists(LMAT_PATH)) {
  e <- new.env(); load(DATA_PATH, envir = e)
  raw <- get(ls(e)[1], envir = e)
  if (inherits(raw, "sf")) raw <- sf::st_drop_geometry(raw)
  raw <- as.data.frame(raw)

  has_code <- "predator_species_code" %in% names(raw)
  if (!has_code) raw$predator_species_code <- NA_integer_

  lmat_ref <- read_csv(LMAT_PATH, show_col_types = FALSE,
                       locale = locale(encoding = "Latin1")) %>%
    transmute(predator_species_code = as.integer(code),
              predator_species_common_name = Predator,
              `Scientific name` = latin_name,
              `Lmat (cm)` = ifelse(tolower(lmat) %in% c("unknown", "na", ""), NA_character_, lmat),
              `Lmax (cm)` = as.character(lmax),
              Source = coalesce(`Information source`, "")) %>%
    distinct(predator_species_code, .keep_all = TRUE)
  join_key <- if (has_code) "predator_species_code" else "predator_species_common_name"
  lmat_ref <- lmat_ref %>% select(-all_of(setdiff(c("predator_species_code", "predator_species_common_name"), join_key)))

  fmt_range <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_character_)
    paste0(format(round(min(x), 1), nsmall = 0), "-", format(round(max(x), 1), nsmall = 0))
  }

  per_period <- raw %>%
    filter(period %in% c("2004-2006", "2018-2019")) %>%
    distinct(predator_species_code, predator_species_common_name, period,
             stomach_id, somatic_length_cm) %>%
    group_by(predator_species_code, predator_species_common_name, period) %>%
    summarise(n_sto = n_distinct(stomach_id),
              range = fmt_range(somatic_length_cm), .groups = "drop")

  tableS1 <- per_period %>%
    pivot_wider(names_from = period, values_from = c(range, n_sto)) %>%
    left_join(lmat_ref, by = join_key) %>%
    transmute(Predator = predator_species_common_name,
              `Scientific name`,
              `Length range 2004-2006 (cm)` = `range_2004-2006`,
              `Length range 2018-2019 (cm)` = `range_2018-2019`,
              `Stomachs 2004-2006` = coalesce(`n_sto_2004-2006`, 0L),
              `Stomachs 2018-2019` = coalesce(`n_sto_2018-2019`, 0L),
              `Lmat (cm)`, `Lmax (cm)`, Source) %>%
    arrange(desc(`Stomachs 2004-2006` + `Stomachs 2018-2019`))

  emit(tableS1, "TableS1_predators",
       paste0("Table S1. The ", nrow(tableS1), " predators retained: length range and number ",
              "of stomachs with prey per period, length at maturity (Lmat) and maximum ",
              "length (Lmax) with their source. Predators without a documented Lmat form ",
              "a single size class."))
} else {
  message("TableS1 skipped: ", DATA_PATH, " or ", LMAT_PATH, " not found.")
}

# =============================================================================
# TABLE S2 - SAMPLING BY CELL
# =============================================================================
# Built from one representative run per level so that the counts are not
# multiplied by the resolution sweep. Stomach counts are identical across
# resolutions within a cell, so taking the first x_threshold is exact.
one_res <- res_all %>%
  filter(level == "all_gulf", currency == "biomass", contrast == "PTa") %>%
  filter(x_threshold == min(x_threshold, na.rm = TRUE))

tableS2 <- one_res %>%
  transmute(Predator = species, `Size class` = size_class,
            `Stomachs 2004-2006` = n_sto_P1, `Stomachs 2018-2019` = n_sto_P2,
            `Sets 2004-2006` = n_set_P1, `Sets 2018-2019` = n_set_P2,
            Testable = testable, Reliable = reliable) %>%
  arrange(Predator, `Size class`)

if (nrow(tableS2)) {
  emit(tableS2, "TableS2_sampling",
       paste0("Table S2. Stomachs and trawl sets per predator and size class, ", CONTRAST_LAB_1L[["PTa"]], " contrast, biomass currency, finest taxonomic resolution."))
} else {
  message("TableS2 skipped: no 2004-2006 vs 2018-2019 biomass run at the Gulf-wide level.")
}

# =============================================================================
# TABLE S3 - THE NINE STATES
# =============================================================================
# What the collapse into families conceals. Reorganisation in particular is
# three quite different states, and their relative weight is worth showing.
tableS3 <- res_all %>%
  filter(level == "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast")) %>%
  filter(!is.na(contrast)) %>%
  mutate(Family = as.character(family_of(diagnostic)),
         pct = round(pct, 1)) %>%
  mutate(contrast = contrast_label(contrast)) %>%
  select(Currency = currency, Contrast = contrast, Family,
         State = diagnostic, `Frequency (%)` = pct) %>%
  arrange(Currency, Contrast,
          factor(Family, levels = FAMILY_LEVELS),
          factor(State, levels = DIAG_LEVELS))

emit(tableS3, "TableS3_states_gulf",
     "Table S3. Frequency (%) of each of the nine diagnostic states, Gulf-wide level.")

# =============================================================================
# TABLES S3 AND S4 - PER UNIT
# =============================================================================
unit_table <- function(lvl, stem, caption) {
  d <- fam_unit %>%
    filter(level == lvl, !is.na(contrast)) %>%
    arrange(currency, contrast, spatial_unit) %>%
    mutate(across(all_of(FAMILY_LEVELS), ~round(.x, 1)),
           flag = ifelse(!is.na(n_units) & n_units < 5, "*", "")) %>%
    transmute(Currency = CURRENCY_LAB[as.character(currency)],
              Contrast = contrast_label(contrast),
              Unit = paste0(area_label(spatial_unit), flag),
              `n units` = n_units,
              Stability, Substitution,
              `Functional change`, Reorganisation)
  if (!nrow(d)) { message(stem, " skipped: level '", lvl, "' not present."); return(invisible(NULL)) }
  emit(d, stem, caption)
}

unit_table("ecoregion", "TableS4_by_ecoregion",
           "Table S4. Family frequency (%) per ecoregion. * = fewer than five predator x size-class units.")

unit_table("stratum", "TableS5_by_stratum",
           "Table S5. Family frequency (%) per survey stratum. * = fewer than five predator x size-class units.")

# =============================================================================
# TABLE S5 - COVERAGE AND CONFOUNDING
# =============================================================================
# The table that qualifies every other one. pct_testable is the share of cells
# where all three signals could be computed; pct_confounded is the share of
# testable cells where the composition and dispersion tests were both
# significant, i.e. where a location shift cannot be separated from a spread
# difference. A high value at the Gulf-wide level is the expected cost of
# pooling across heterogeneous prey fields.
tableS6 <- coverage_by(res_all, keys = c("level", "currency", "contrast")) %>%
  filter(!is.na(contrast)) %>%
  arrange(level, currency, contrast) %>%
  transmute(`Spatial level` = LEVEL_SHORT[as.character(level)],
            Currency = CURRENCY_LAB[as.character(currency)],
            Contrast = contrast_label(contrast),
            `Cells (all resolutions)` = n_cells,
            `Testable (%)` = pct_testable,
            `Reliable (%)` = pct_reliable,
            `Location/dispersion confounded (%)` = pct_confounded)

emit(tableS6, "TableS6_coverage",
     paste0("Table S6. Coverage and confounding by level and contrast. Percentages ",
            "in Tables 2 to S4 are computed on testable cells only."))

# =============================================================================
# CONSISTENCY CHECKS (warn, never stop)
# =============================================================================
cat("\n\n===== Consistency checks =====\n")

chk_sum <- fam_gulf %>%
  mutate(total = rowSums(across(all_of(FAMILY_LEVELS)))) %>%
  filter(abs(total - 100) > 0.01)
if (nrow(chk_sum)) {
  warning("Family percentages do not sum to 100 in ", nrow(chk_sum),
          " Gulf-wide row(s). Inconclusive cells may not have been excluded.",
          call. = FALSE)
  print(chk_sum)
} else {
  cat("  Family percentages sum to 100 at the Gulf-wide level.\n")
}

n_lv <- dplyr::n_distinct(res_all$level)
cat("  Spatial levels present: ", n_lv, " of 3",
    if (n_lv < 3) "  <- run the missing 6b/6c/6d before quoting Table 3" else "", "\n", sep = "")

miss_ct <- setdiff(CONTRAST_LEVELS, unique(as.character(res_all$contrast)))
if (length(miss_ct)) {
  cat("  Contrasts absent from the results: ", paste(miss_ct, collapse = ", "), "\n", sep = "")
}

low_cov <- tableS6 %>% filter(`Testable (%)` < 60)
if (nrow(low_cov)) {
  cat("  Level/contrast combinations below 60 % testable cells:\n")
  print(as.data.frame(low_cov), row.names = FALSE)
}

cat("\nTables written to ", normalizePath(DIR_TABLES), "\n", sep = "")
