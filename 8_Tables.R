# =============================================================================
# 8_Tables.R - MANUSCRIPT TABLES
# -----------------------------------------------------------------------------
# Every numbered and supplementary table, built from the same aggregation rule
# as the figures so that a number quoted in the text, a number in a table and a
# point on a figure cannot disagree.
#
# Each table is written three ways:
#   .csv   the machine-readable version, for the data archive
#   .html  a formatted version to paste into the manuscript
#   the console, so the numbers can be read without leaving the session
#
# WHAT THIS SCRIPT PRODUCES
# -----------------------------------------------------------------------------
#   Table1_typology          The nine diagnostic states, the signal combination
#                            that defines each, the family it belongs to, and
#                            its ecological reading. This is the key to every
#                            other table; it comes from the configuration
#                            module, not from the data.
#
#   Table2_headline          Family frequency by contrast and currency at the
#                            Gulf-wide level. The table behind the Results
#                            paragraph on the inter-decade inversion.
#
#   Table3_crossscale        Family frequency by spatial level for the
#                            between-decade contrast. The §3.3 gradient.
#
#   TableS1_sampling         Stomachs and cells per predator, size class and
#                            period. The sample-size table reviewers ask for.
#
#   TableS2_states_gulf      The nine states, not just the four families, at
#                            the Gulf-wide level. Shows what the collapse into
#                            families hides.
#
#   TableS3_by_ecoregion     Family frequency per ecoregion.
#   TableS4_by_stratum       Family frequency per stratum, with the low-sample
#                            flag.
#
#   TableS5_coverage         Per level and contrast: how many cells existed, how
#                            many were testable and reliable, and how often
#                            composition and dispersion were both significant.
#                            A family composition means little without this.
#
# READ TableS5 BEFORE THE OTHERS. Percentages in Tables 2 to S4 are computed on
# testable cells only. If coverage is 45 % at one level and 90 % at another, the
# two family compositions are not describing the same thing.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(knitr)
})

source("R_helpers/Config_Mappings.R")
PREY_FAMILY <- "_1"
has_kable_extra <- requireNamespace("kableExtra", quietly = TRUE)

# -----------------------------------------------------------------------------
# Writer: csv + html + console, under one stem
# -----------------------------------------------------------------------------
emit <- function(df, stem, caption, digits = 1) {
  write_csv(df, file.path(DIR_TABLES, paste0(stem, ".csv")))

  kt <- knitr::kable(df, format = "html", digits = digits, caption = caption)
  if (has_kable_extra) {
    kt <- kableExtra::kable_styling(
      kt, bootstrap_options = c("striped", "condensed"), full_width = FALSE)
  }
  writeLines(as.character(kt), file.path(DIR_TABLES, paste0(stem, ".html")))

  cat("\n\n===== ", caption, " =====\n", sep = "")
  print(as.data.frame(df), row.names = FALSE, digits = 4)
  cat("  -> ", stem, ".csv / .html\n", sep = "")
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
# Definitional, not empirical: it states what each label means before any
# frequency is reported. The ecological readings are the ones carried in the
# engine's diet_diagnostics_table, restated here so this script does not need
# the engine.
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
            Contrast = as.character(contrast),
            `Spatial level` = LEVEL_SHORT[as.character(level)],
            `Units` = ifelse(is.na(n_units_total), 1L, n_units_total),
            Stability, Substitution,
            `Functional change`, Reorganisation)

emit(table3, "Table3_crossscale",
     paste0("Table 3. Family frequency (%) by spatial level. Gulf-wide is a single ",
            "pooled test on all trawl sets; ecoregion and stratum values are means ",
            "across units, each unit tested on its own data."))

# =============================================================================
# TABLE S1 - SAMPLING
# =============================================================================
# Built from one representative run per level so that the counts are not
# multiplied by the resolution sweep. Stomach counts are identical across
# resolutions within a cell, so taking the first x_threshold is exact.
one_res <- res_all %>%
  filter(level == "all_gulf", currency == "biomass", contrast == "PTa") %>%
  filter(x_threshold == min(x_threshold, na.rm = TRUE))

tableS1 <- one_res %>%
  transmute(Predator = species, `Size class` = size_class,
            `Stomachs P1` = n_sto_P1, `Stomachs P2` = n_sto_P2,
            `Sets P1` = n_set_P1, `Sets P2` = n_set_P2,
            Testable = testable, Reliable = reliable) %>%
  arrange(Predator, `Size class`)

if (nrow(tableS1)) {
  emit(tableS1, "TableS1_sampling",
       "Table S1. Stomachs and trawl sets per predator and size class, PTa contrast, biomass currency, finest taxonomic resolution.")
} else {
  message("TableS1 skipped: no PTa biomass run at the Gulf-wide level.")
}

# =============================================================================
# TABLE S2 - THE NINE STATES
# =============================================================================
# What the collapse into families conceals. Reorganisation in particular is
# three quite different states, and their relative weight is worth showing.
tableS2 <- res_all %>%
  filter(level == "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast")) %>%
  filter(!is.na(contrast)) %>%
  mutate(Family = as.character(family_of(diagnostic)),
         pct = round(pct, 1)) %>%
  select(Currency = currency, Contrast = contrast, Family,
         State = diagnostic, `Frequency (%)` = pct) %>%
  arrange(Currency, Contrast,
          factor(Family, levels = FAMILY_LEVELS),
          factor(State, levels = DIAG_LEVELS))

emit(tableS2, "TableS2_states_gulf",
     "Table S2. Frequency (%) of each of the nine diagnostic states, Gulf-wide level.")

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
              Contrast = as.character(contrast),
              Unit = paste0(spatial_unit, flag),
              `n units` = n_units,
              Stability, Substitution,
              `Functional change`, Reorganisation)
  if (!nrow(d)) { message(stem, " skipped: level '", lvl, "' not present."); return(invisible(NULL)) }
  emit(d, stem, caption)
}

unit_table("ecoregion", "TableS3_by_ecoregion",
           "Table S3. Family frequency (%) per ecoregion. * = fewer than five predator x size-class units.")

unit_table("stratum", "TableS4_by_stratum",
           "Table S4. Family frequency (%) per survey stratum. * = fewer than five predator x size-class units.")

# =============================================================================
# TABLE S5 - COVERAGE AND CONFOUNDING
# =============================================================================
# The table that qualifies every other one. pct_testable is the share of cells
# where all three signals could be computed; pct_confounded is the share of
# testable cells where the composition and dispersion tests were both
# significant, i.e. where a location shift cannot be separated from a spread
# difference. A high value at the Gulf-wide level is the expected cost of
# pooling across heterogeneous prey fields.
tableS5 <- coverage_by(res_all, keys = c("level", "currency", "contrast")) %>%
  filter(!is.na(contrast)) %>%
  arrange(level, currency, contrast) %>%
  transmute(`Spatial level` = LEVEL_SHORT[as.character(level)],
            Currency = CURRENCY_LAB[as.character(currency)],
            Contrast = as.character(contrast),
            `Cells (all resolutions)` = n_cells,
            `Testable (%)` = pct_testable,
            `Reliable (%)` = pct_reliable,
            `Location/dispersion confounded (%)` = pct_confounded)

emit(tableS5, "TableS5_coverage",
     paste0("Table S5. Coverage and confounding by level and contrast. Percentages ",
            "in Tables 2 to S4 are computed on testable cells only."))

# =============================================================================
# CONSISTENCY CHECKS
# =============================================================================
# Cheap assertions that would have caught the aggregation errors this pipeline
# has already been through. They warn rather than stop, so a partial run still
# produces its tables.
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

low_cov <- tableS5 %>% filter(`Testable (%)` < 60)
if (nrow(low_cov)) {
  cat("  Level/contrast combinations below 60 % testable cells:\n")
  print(as.data.frame(low_cov), row.names = FALSE)
}

cat("\nTables written to ", normalizePath(DIR_TABLES), "\n", sep = "")
