# =============================================================================
# 6f_Compare_Levels.R - CROSS-LEVEL CHECKS, ALL PREY-GROUPING FAMILIES
# -----------------------------------------------------------------------------
# Run after 6b, 6c, 6d (and 6e for the alternative families). For each family
# present ("_1", "_2", "_PP") it reads the three all_runs files and writes to
# Sensitivity_comparison/ :
#   crosswalk_strata_ecoregions.csv         strata split across ecoregions (data)
#   <fam>_levels_sample.csv                 cells, units, stomachs, coverage per level
#   <fam>_levels_composition.csv            diagnostic composition per level
#   <fam>_levels_effect_vs_power.csv        effect sizes vs detection rates
#   <fam>_levels_confounded.csv             confounded location/dispersion cells
#   <fam>_eco_vs_strata.csv                 ecoregional verdict vs stratum verdicts
#   <fam>_eco_vs_strata_match.csv           share of ecoregional cells matching the
#                                           modal verdict of their strata
#   <fam>_diagnostic_composition_by_level.png
#   families_diagnostic_composition_by_level.png   all families side by side
# Ecoregion and stratum are alternative partitions, not a nested hierarchy; the
# ecoregion-stratum agreement uses the strata lying in one ecoregion only.
# The manuscript tables and figures come from 7, 8 and 9; this script is a check.
# =============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(sf); library(readr)
})

FAMILIES <- c("_1", "_2", "_PP")
ALPHA    <- 0.05
OBJ      <- "res_biomass_PT"          # reference comparison: between periods, biomass
OUT_DIR  <- "Sensitivity_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FAM_LAB <- c("_1" = "Pooled q = 1 (main)", "_2" = "Pooled q = 2", "_PP" = "Per predator")

# =============================================================================
# 0. STRATUM-TO-ECOREGION CROSSWALK (prey records, shared by all families)
# =============================================================================
load("data/dat_classed.rda")
xwalk_raw <- as.data.frame(dat_classed) %>%
  st_drop_geometry() %>%
  mutate(stratum = as.character(stratum)) %>%
  count(Area, stratum, name = "n_rec") %>%
  filter(n_rec > 0)
rm(dat_classed)

xwalk <- xwalk_raw %>%
  group_by(stratum) %>%
  mutate(n_stratum = sum(n_rec), share = n_rec / n_stratum, n_ecoregions = dplyr::n()) %>%
  slice_max(n_rec, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(stratum, eco_majority = Area, share_majority = round(share, 3),
            n_ecoregions, pure = n_ecoregions == 1)

split_strata <- xwalk_raw %>%
  group_by(stratum) %>%
  filter(dplyr::n() > 1) %>%
  mutate(share = round(100 * n_rec / sum(n_rec), 1)) %>%
  ungroup() %>%
  arrange(stratum, desc(n_rec))
write_csv(split_strata, file.path(OUT_DIR, "crosswalk_strata_ecoregions.csv"))

cat("\n--- Strata split across ecoregions ---\n"); print(split_strata, n = Inf)
cat("\n--- Crosswalk coverage ---\n")
print(xwalk %>% summarise(
  n_strata = dplyr::n(), n_pure = sum(pure), n_split = sum(!pure),
  pct_records_pure = round(100 * sum(xwalk_raw$n_rec[xwalk_raw$stratum %in% stratum[pure]]) /
                             sum(xwalk_raw$n_rec), 1)))

# =============================================================================
# 1. ONE FAMILY
# =============================================================================
tidy_level <- function(res) {
  d <- res$results
  sp_col <- intersect(c("Area", "str"), names(d))[1]
  d$spatial_unit <- if (is.na(sp_col)) NA_character_ else as.character(d[[sp_col]])
  dplyr::select(d, -dplyr::any_of(c("Area", "str")))
}

load_level <- function(level, fam) {
  f <- paste0("data/Sensitivity/all_runs_", level, fam, ".rda")
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e)
  if (!exists(OBJ, envir = e)) return(NULL)
  tidy_level(get(OBJ, envir = e))
}

emit <- function(df, stem, title) {
  write_csv(df, file.path(OUT_DIR, paste0(stem, ".csv")))
  cat("\n--- ", title, " ---\n", sep = ""); print(as.data.frame(df), row.names = FALSE)
  invisible(df)
}

compare_family <- function(fam) {
  L <- list(all_gulf = load_level("all_gulf", fam),
            ecoregion = load_level("ecoregion", fam),
            stratum   = load_level("stratum", fam))
  L <- L[!vapply(L, is.null, logical(1))]
  if (length(L) < 2) { message("Family ", fam, ": fewer than two levels found, skipped."); return(NULL) }
  cat("\n\n===================== FAMILY ", fam, " (", FAM_LAB[[fam]], ") =====================\n", sep = "")
  all3 <- bind_rows(L)

  emit(all3 %>% group_by(spatial_level) %>%
         summarise(n_cells = dplyr::n(), n_units = dplyr::n_distinct(spatial_unit),
                   n_resolutions = dplyr::n_distinct(x_threshold),
                   stomachs = sum(n_sto_P1 + n_sto_P2) / dplyr::n_distinct(x_threshold),
                   pct_testable = round(100 * mean(testable), 1),
                   pct_reliable = round(100 * mean(reliable), 1), .groups = "drop"),
       paste0(fam, "_levels_sample"), "Cells, units and stomach totals per level")

  comp <- all3 %>% filter(testable) %>% count(spatial_level, diagnostic) %>%
    group_by(spatial_level) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
  emit(comp %>% select(-n) %>% pivot_wider(names_from = spatial_level, values_from = pct, values_fill = 0),
       paste0(fam, "_levels_composition"), "Diagnostic composition per level (testable cells)")

  emit(all3 %>% filter(testable) %>% group_by(spatial_level) %>%
         summarise(n = dplyr::n(),
                   k_med = median(pmin(n_set_P1, n_set_P2), na.rm = TRUE),
                   BC_med = round(median(BC, na.rm = TRUE), 3),
                   R2_med = round(median(R2_comp, na.rm = TRUE), 3),
                   pct_H_sig = round(100 * mean(p_H < ALPHA, na.rm = TRUE), 1),
                   pct_Bs_sig = round(100 * mean(p_Bs < ALPHA, na.rm = TRUE), 1),
                   pct_comp_sig = round(100 * mean(p_comp < ALPHA, na.rm = TRUE), 1), .groups = "drop"),
       paste0(fam, "_levels_effect_vs_power"), "Effect size vs detection rate, by level")

  emit(all3 %>% filter(testable) %>% group_by(spatial_level) %>%
         summarise(pct_comp_sig = round(100 * mean(p_comp < ALPHA, na.rm = TRUE), 1),
                   pct_disp_sig = round(100 * mean(p_disp < ALPHA, na.rm = TRUE), 1),
                   pct_confounded = round(100 * mean(comp_dispersion, na.rm = TRUE), 1), .groups = "drop"),
       paste0(fam, "_levels_confounded"), "Confounded location/dispersion cells, by level")

  # Ecoregion vs stratum partitions, on strata lying entirely in one ecoregion.
  if (all(c("ecoregion", "stratum") %in% names(L))) {
    key <- c("x_threshold", "species", "size_class")
    L3 <- L$stratum %>% filter(testable) %>% left_join(xwalk, by = c("spatial_unit" = "stratum"))
    cat("\n--- Coverage of the ecoregion / stratum comparison ---\n")
    print(L3 %>% summarise(n_stratum_cells = dplyr::n(), n_kept = sum(pure, na.rm = TRUE),
                           pct_kept = round(100 * mean(pure, na.rm = TRUE), 1)))
    L3_pure <- L3 %>% filter(pure)
    L2 <- L$ecoregion %>% filter(testable) %>%
      select(all_of(key), eco_unit = spatial_unit, d_eco = diagnostic)

    emit(inner_join(L2, L3_pure %>% select(all_of(key), eco_unit = eco_majority, d_str = diagnostic),
                    by = c(key, "eco_unit"), relationship = "many-to-many") %>%
           count(d_eco, d_str) %>% group_by(d_eco) %>%
           mutate(pct_within_eco = round(100 * n / sum(n), 1)) %>% ungroup() %>% arrange(d_eco, desc(n)),
         paste0(fam, "_eco_vs_strata"), "Ecoregional verdict vs stratum verdicts (unambiguous strata)")

    emit(inner_join(L2, L3_pure %>% select(all_of(key), eco_unit = eco_majority, d_str = diagnostic),
                    by = c(key, "eco_unit"), relationship = "many-to-many") %>%
           count(across(all_of(c(key, "eco_unit"))), d_eco, d_str) %>%
           group_by(across(all_of(c(key, "eco_unit")))) %>%
           slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
           summarise(n_eco_cells = dplyr::n(), pct_match = round(100 * mean(d_eco == d_str), 1)),
         paste0(fam, "_eco_vs_strata_match"), "Ecoregion vs modal stratum verdict")
  }

  p <- comp %>%
    mutate(spatial_level = factor(spatial_level, levels = c("all_gulf", "ecoregion", "stratum"))) %>%
    ggplot(aes(x = reorder(diagnostic, pct), y = pct, fill = spatial_level)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2", name = "Spatial level") +
    labs(title = paste0("Diagnostic composition across spatial levels - ", FAM_LAB[[fam]]),
         subtitle = "Biomass, 2004-2006 vs 2018-2019, testable cells only",
         x = NULL, y = "Within-level frequency (%)") +
    theme_minimal() +
    theme(axis.text.y = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(OUT_DIR, paste0(fam, "_diagnostic_composition_by_level.png")), p,
         width = 10, height = 6, dpi = 300)

  comp %>% mutate(family = fam)
}

# =============================================================================
# 2. ALL FAMILIES
# =============================================================================
comp_all <- bind_rows(lapply(FAMILIES, compare_family))

if (nrow(comp_all)) {
  write_csv(comp_all, file.path(OUT_DIR, "families_levels_composition.csv"))
  p_all <- comp_all %>%
    mutate(spatial_level = factor(spatial_level, levels = c("all_gulf", "ecoregion", "stratum")),
           family = factor(FAM_LAB[family], levels = unname(FAM_LAB))) %>%
    ggplot(aes(x = reorder(diagnostic, pct), y = pct, fill = spatial_level)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
    coord_flip() +
    facet_wrap(~family, nrow = 1) +
    scale_fill_brewer(palette = "Set2", name = "Spatial level") +
    labs(title = "Diagnostic composition across spatial levels, by prey-grouping family",
         subtitle = "Biomass, 2004-2006 vs 2018-2019, testable cells only",
         x = NULL, y = "Within-level frequency (%)") +
    theme_minimal() +
    theme(axis.text.y = element_text(face = "bold"), legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT_DIR, "families_diagnostic_composition_by_level.png"), p_all,
         width = 4.5 * dplyr::n_distinct(comp_all$family) + 2, height = 6, dpi = 300)
}

cat("\nCross-level comparison written to ", normalizePath(OUT_DIR), "\n", sep = "")
