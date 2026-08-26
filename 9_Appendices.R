# =============================================================================
# 9_Appendices.R - APPENDICES, AUDITS AND REPRODUCIBILITY
# -----------------------------------------------------------------------------
# Everything a reviewer will ask for that is not a body figure or a numbered
# table: the sampling map, the taxonomic-resolution sweep, the power and
# coverage audits, the partition crosswalk, and a manifest recording what
# produced the outputs.
#
# WHAT THIS SCRIPT PRODUCES
# -----------------------------------------------------------------------------
#   Fig2_sampling_map            Manuscript Figure 2. Stomachs analysed per
#                                stratum and survey year, with ecoregion
#                                outlines. Answers "where and when were the
#                                data collected" before any result is shown.
#
#   Fig3_resolution_sweep        Manuscript Figure 3. Number of distinct prey
#                                categories retained as the aggregation
#                                threshold x rises from 10 to 1000, plus the
#                                family frequencies across the same sweep. The
#                                second panel is the robustness claim: the
#                                signal varies smoothly, so no single resolution
#                                drives the result.
#
#   FigS7_flagged_units          Which spatial units fall below the reliability
#                                thresholds, and what the composites look like
#                                with and without them. The sensitivity check
#                                behind "restricting to better-sampled strata
#                                left the composites essentially unchanged".
#
#   TableA1_set_depth            Trawl sets per cell by level. The binding
#                                constraint on every test in the framework.
#
#   TableA2_power_vs_effect      Effect sizes (BC, R2) next to detection rates,
#                                by level. Separates "smaller effects at finer
#                                grain" from "less power at finer grain".
#
#   TableA3_partition_crosswalk  Which strata straddle an ecoregion boundary.
#                                Ecoregion and stratum are alternative
#                                partitions, not a nested hierarchy, and this
#                                is the table that documents it.
#
#   TableA4_untestable           Why cells could not be classified, by level.
#
#   manifest.txt                 Session info, package versions, file
#                                inventory with sizes and timestamps, and the
#                                analytical parameters read back from the
#                                results. Attach to the data archive.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(sf); library(forcats)
})

PREY_FAMILY <- "_1"          # prey-grouping family used for the manuscript
source("R_helpers/Config_Mappings.R")   # builds RDA_DIRS from PREY_FAMILY

DATA_PATH   <- "data/dat_classed.rda"
STRATA_PATH <- "strata_rv_gulf.rds"
MAP_XLIM    <- c(-66.2, -60.0)
MAP_YLIM    <- c(45.5, 49.2)
N_FLAG      <- 5

save_app <- function(p, stem, w, h, dpi = 300) save_fig(p, stem, w, h, dpi, dir = DIR_APPEND)

emit_app <- function(df, stem, caption) {
  # Contrast codes never reach the appendix files: relabel to periods.
  if ("contrast" %in% names(df)) df <- dplyr::mutate(df, contrast = contrast_label(contrast))
  write_csv(df, file.path(DIR_APPEND, paste0(stem, ".csv")))
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

have_raw <- file.exists(DATA_PATH)
if (have_raw) {
  e <- new.env(); load(DATA_PATH, envir = e)
  raw <- get(ls(e)[1], envir = e)
  # dat_classed may still carry an sf geometry column; dropping it keeps the
  # distinct() and count() calls below from doing spatial work they don't need.
  if (inherits(raw, "sf")) raw <- sf::st_drop_geometry(raw)
  raw <- as.data.frame(raw)
  cat("Raw data: ", nrow(raw), " prey records.\n", sep = "")
} else {
  raw <- NULL
  message("Raw data not found at ", DATA_PATH,
          "; Figures 2 and 3 and the crosswalk are skipped.")
}

# =============================================================================
# 2. FIGURE 2 - SAMPLING
# =============================================================================
# Stomachs, not prey records: one record per prey item would over-weight
# stomachs holding many taxa.
if (have_raw) {
  cat("\nFigure 2: sampling\n")

  samp <- raw %>%
    filter(year %in% c(2004, 2005, 2006, 2018, 2019)) %>%
    distinct(stomach_id, year, Area, stratum) %>%
    count(year, Area, stratum, name = "n_stomachs")

  # Panel A: stomachs per ecoregion and year.
  pA <- samp %>%
    count(year, Area, wt = n_stomachs, name = "n") %>%
    ggplot(aes(factor(year), n, fill = Area)) +
    geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              size = 2.6, colour = "white", fontface = "bold") +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(title = "A. Stomachs analysed per year and ecoregion",
         x = NULL, y = "Stomachs") +
    theme_diag()

  save_app(pA, "Fig2a_sampling_by_year", 8, 4.6)

  # Panel B: the same on the map, if the geometry is available.
  strata_sf <- load_strata(STRATA_PATH)   # gulf.spatial shapefile, else the .rds
  if (!is.null(strata_sf)) {
    strata_sf <- st_make_valid(strata_sf)
    if (is.na(st_crs(strata_sf))) st_crs(strata_sf) <- 4326
    strata_sf$str <- as.character(strata_sf$str)

    d <- strata_sf %>%
      left_join(samp %>% mutate(str = as.character(stratum)) %>%
                  count(str, year, wt = n_stomachs, name = "n"),
                by = "str") %>%
      filter(!is.na(n))

    pB <- ggplot() +
      geom_sf(data = strata_sf, fill = "grey97", colour = "grey75", linewidth = 0.12) +
      geom_sf(data = d, aes(fill = n), colour = "grey45", linewidth = 0.12) +
      facet_wrap(~year, nrow = 1) +
      scale_fill_viridis_c(option = "rocket", direction = -1, trans = "sqrt",
                           name = "Stomachs") +
      coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
      labs(title = "B. Stomachs analysed per stratum and survey year") +
      theme_diag(base_size = 10) +
      theme(axis.title = element_blank(), axis.text = element_text(size = 5.5),
            strip.text = element_text(face = "bold"),
            panel.grid = element_line(colour = "grey93"))

    save_app(pB, "Fig2b_sampling_map", 14, 4.2, dpi = 200)
  }
}

# =============================================================================
# 3. FIGURE 3 - TAXONOMIC RESOLUTION SWEEP
# =============================================================================
# Panel A counts the prey categories that survive each aggregation threshold.
# Panel B carries the robustness claim: the family frequencies move smoothly
# with x, so the reported means are not an artefact of one resolution.
cat("\nFigure 3: resolution sweep\n")

if (have_raw) {
  prey_cols <- grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"), names(raw), value = TRUE)
  if (length(prey_cols)) {
    ncat <- tibble(
      col = prey_cols,
      x_threshold = as.numeric(str_extract(prey_cols, "\\d+")),
      n_cat = vapply(prey_cols, function(cc) dplyr::n_distinct(raw[[cc]], na.rm = TRUE),
                     numeric(1))
    ) %>% arrange(x_threshold)

    p3a <- ggplot(ncat, aes(x_threshold, n_cat)) +
      geom_line(linewidth = 0.8, colour = "#2F4A5A") +
      geom_point(size = 1.3, colour = "#2F4A5A") +
      labs(title = "A. Prey categories retained by aggregation threshold",
           subtitle = "Each record is kept at its finest level, or walked up the taxonomic tree until the node reaches x stomachs",
           x = "Threshold x (stomachs)", y = "Distinct prey categories") +
      theme_diag()

    save_app(p3a, "Fig3a_prey_categories", 8, 4.4)
    emit_app(ncat %>% select(-col), "TableA0_prey_categories",
             "Prey categories retained at each aggregation threshold.")
  } else {
    message("No prey_category_*", PREY_FAMILY, " columns found; Figure 3A skipped.")
  }
}

# Panel B: family frequency across the sweep, Gulf-wide, both currencies.
sweep_fam <- res_all %>%
  filter(level == "all_gulf", diagnostic %in% DIAG_LEVELS) %>%
  count(currency, contrast, x_threshold, diagnostic, name = "n") %>%
  group_by(currency, contrast, x_threshold) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(family = family_of(diagnostic)) %>%
  filter(!is.na(family), !is.na(contrast)) %>%
  group_by(currency, contrast, x_threshold, family) %>%
  summarise(pct = sum(pct), .groups = "drop")

p3b <- ggplot(sweep_fam, aes(x_threshold, pct, colour = family)) +
  geom_line(linewidth = 0.6) +
  facet_grid(currency ~ contrast,
             labeller = labeller(currency = CURRENCY_LAB,
                                 contrast = CONTRAST_LAB_1L)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "B. Sensitivity of the diagnostic composition to prey taxonomic resolution",
       subtitle = "Family frequencies across the full sweep; the reported values are the means of these curves",
       x = "Threshold x (stomachs)", y = "Frequency (%)") +
  theme_diag() +
  theme(strip.text = element_text(size = 8))

save_app(p3b, "Fig3b_sweep_families", 13, 5.6)

# Spread across the sweep: the uncertainty quoted as error bars in the body.
sweep_spread <- sweep_fam %>%
  group_by(currency, contrast, family) %>%
  summarise(mean_pct = round(mean(pct), 1),
            sd_pct   = round(sd(pct), 1),
            min_pct  = round(min(pct), 1),
            max_pct  = round(max(pct), 1),
            n_resolutions = dplyr::n(), .groups = "drop")

emit_app(sweep_spread, "TableA5_resolution_spread",
         "Family frequency across the taxonomic-resolution sweep: mean, spread and range.")

# =============================================================================
# 4. AUDITS
# =============================================================================
cat("\nAudits\n")

# Set depth: the number of tows per cell, the constraint everything rests on.
tableA1 <- res_all %>%
  mutate(k = pmin(n_set_P1, n_set_P2)) %>%
  group_by(level, currency) %>%
  summarise(
    n_cells   = dplyr::n(),
    k_min     = min(k, na.rm = TRUE),
    k_q25     = quantile(k, 0.25, na.rm = TRUE),
    k_median  = median(k, na.rm = TRUE),
    k_q75     = quantile(k, 0.75, na.rm = TRUE),
    pct_k_lt4 = round(100 * mean(k < 4, na.rm = TRUE), 1),
    pct_k_lt6 = round(100 * mean(k < 6, na.rm = TRUE), 1),
    stomachs_per_set = round(median((n_sto_P1 + n_sto_P2) /
                                      pmax(n_set_P1 + n_set_P2, 1), na.rm = TRUE), 1),
    .groups = "drop"
  )

emit_app(tableA1, "TableA1_set_depth",
         "Trawl sets per cell by spatial level. k = min(sets in 2004-2006, sets in 2018-2019).")

# Effect size next to detection rate.
tableA2 <- res_all %>%
  filter(testable) %>%
  group_by(level, currency, contrast) %>%
  summarise(
    n            = dplyr::n(),
    k_median     = median(pmin(n_set_P1, n_set_P2), na.rm = TRUE),
    BC_median    = round(median(BC, na.rm = TRUE), 3),
    R2_median    = round(median(R2_comp, na.rm = TRUE), 3),
    pct_H_sig    = round(100 * mean(p_H    < 0.05, na.rm = TRUE), 1),
    pct_Bs_sig   = round(100 * mean(p_Bs   < 0.05, na.rm = TRUE), 1),
    pct_comp_sig = round(100 * mean(p_comp < 0.05, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  filter(!is.na(contrast)) %>%
  arrange(currency, contrast, level)

emit_app(tableA2, "TableA2_power_vs_effect",
         "Effect sizes and detection rates by level. BC and R2 do not depend on significance; the pct_*_sig columns do.")

# Why cells were not classifiable.
tableA4 <- res_all %>%
  group_by(level, currency, contrast) %>%
  summarise(
    n_cells          = dplyr::n(),
    n_inconclusive   = sum(!testable),
    n_no_H           = sum(is.na(p_H)),
    n_no_Bs          = sum(is.na(p_Bs)),
    n_no_composition = sum(is.na(p_comp)),
    n_below_min_sets = sum(n_set_P1 < 3 | n_set_P2 < 3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(contrast), n_inconclusive > 0) %>%
  arrange(desc(n_inconclusive))

emit_app(tableA4, "TableA4_untestable",
         "Cells that could not be classified, and which signal was missing.")

# =============================================================================
# 5. PARTITION CROSSWALK
# =============================================================================
# Ecoregions are assigned per stomach by a spatial join; strata are the survey
# polygons. The two geometries do not coincide, so some strata are split across
# ecoregions and the two levels are alternative partitions rather than a nested
# hierarchy. This table is what the Methods should cite on that point.
if (have_raw && all(c("Area", "stratum") %in% names(raw))) {
  cat("\nPartition crosswalk\n")

  xw <- raw %>%
    mutate(stratum = as.character(stratum)) %>%
    count(Area, stratum, name = "n_rec") %>%
    group_by(stratum) %>%
    mutate(n_stratum = sum(n_rec),
           share = round(100 * n_rec / n_stratum, 1),
           n_ecoregions = dplyr::n()) %>%
    ungroup() %>%
    arrange(desc(n_ecoregions), stratum, desc(n_rec))

  tableA3 <- xw %>%
    filter(n_ecoregions > 1) %>%
    select(Stratum = stratum, Ecoregion = Area,
           `Records` = n_rec, `Share of stratum (%)` = share)

  emit_app(tableA3, "TableA3_partition_crosswalk",
           "Strata split across more than one ecoregion. Ecoregion is assigned per stomach by spatial join; strata are survey polygons.")

  summ <- xw %>%
    distinct(stratum, n_stratum, n_ecoregions) %>%
    summarise(
      n_strata     = dplyr::n(),
      n_split      = sum(n_ecoregions > 1),
      pct_records_in_split = round(
        100 * sum(n_stratum[n_ecoregions > 1]) / sum(n_stratum), 1)
    )
  cat("\nCrosswalk summary:\n"); print(as.data.frame(summ), row.names = FALSE)
}

# =============================================================================
# 6. FIGURE S7 - SENSITIVITY TO LOW-SAMPLE UNITS
# =============================================================================
# The composites are recomputed with and without the units below N_FLAG. If the
# two versions agree, the spatial signal does not rest on the sparse units.
cat("\nFigure S7: low-sample sensitivity\n")

fam_unit <- res_all %>%
  filter(level != "all_gulf") %>%
  freq_by_resolution(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  add_families(keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  left_join(
    n_units_by(filter(res_all, level != "all_gulf"),
               keys = c("currency", "contrast", "level", "spatial_unit")),
    by = c("currency", "contrast", "level", "spatial_unit")
  )

sens <- bind_rows(
  fam_unit %>% mutate(subset = "All units"),
  fam_unit %>% filter(!is.na(n_units), n_units >= N_FLAG) %>%
    mutate(subset = paste0("n units >= ", N_FLAG))
) %>%
  group_by(level, currency, contrast, subset) %>%
  summarise(across(all_of(FAMILY_LEVELS), ~mean(.x, na.rm = TRUE)),
            n_units_kept = dplyr::n(), .groups = "drop") %>%
  filter(!is.na(contrast)) %>%
  pivot_longer(all_of(FAMILY_LEVELS), names_to = "family", values_to = "pct") %>%
  mutate(family = factor(family, levels = FAMILY_LEVELS))

pS7 <- ggplot(sens, aes(contrast, pct, colour = family,
                        shape = subset, group = interaction(family, subset))) +
  geom_line(aes(linetype = subset), linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_x_discrete(labels = CONTRAST_LAB) +
  facet_grid(currency ~ level,
             labeller = labeller(currency = CURRENCY_LAB, level = LEVEL_SHORT)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity of the across-unit composites to sparsely sampled units",
       subtitle = paste0("Solid = every unit; dashed = units built on at least ",
                         N_FLAG, " predator x size-class cells"),
       x = NULL, y = "Frequency (%)") +
  theme_diag()

save_app(pS7, "FigS7_flagged_units", 10, 6)

emit_app(
  sens %>%
    pivot_wider(names_from = family, values_from = pct) %>%
    mutate(across(all_of(FAMILY_LEVELS), ~round(.x, 1))),
  "TableA6_lowN_sensitivity",
  "Across-unit composites with and without sparsely sampled units.")

# =============================================================================
# 7. MANIFEST
# =============================================================================
# Written last so it records everything the run produced. Attach to the archive.
cat("\nManifest\n")

manifest <- file.path(DIR_APPEND, "manifest.txt")
con <- file(manifest, "w")

wl <- function(...) writeLines(paste0(...), con)

wl("TROPHIC TRANSITION ANALYSIS - RUN MANIFEST")
wl("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
wl(strrep("-", 78))
wl("")
wl("SPATIAL LEVELS PRESENT")
for (lv in levels(res_all$level)) {
  n <- sum(res_all$level == lv, na.rm = TRUE)
  if (n > 0) wl("  ", lv, ": ", n, " cell-rows across all resolutions and contrasts")
}
wl("")
wl("CONTRASTS PRESENT")
for (ct in CONTRAST_LEVELS) {
  n <- sum(as.character(res_all$contrast) == ct, na.rm = TRUE)
  wl("  ", CONTRAST_LAB_1L[[ct]], ": ", n, " cell-rows")
}
wl("")
wl("TAXONOMIC RESOLUTIONS")
wl("  ", dplyr::n_distinct(res_all$x_threshold), " thresholds, from ",
   min(res_all$x_threshold, na.rm = TRUE), " to ",
   max(res_all$x_threshold, na.rm = TRUE))
wl("")
wl("TEST SETTINGS READ BACK FROM THE RESULTS")
wl("  Niche breadth test: ", paste(unique(res_all$Bs_test), collapse = ", "))
wl("  Prey-grouping family: ",
   if ("prey_family" %in% names(res_all)) paste(unique(res_all$prey_family), collapse = ", ") else PREY_FAMILY)
wl("  Result folders: ", paste(RDA_DIRS, collapse = ", "))
wl("")
wl("OUTPUT FILES")
for (d in c(DIR_FIGURES, DIR_TABLES, DIR_APPEND, DIR_DRIVERS)) {
  if (!dir.exists(d)) next
  fs <- list.files(d, full.names = TRUE)
  if (!length(fs)) next
  wl("  [", d, "]")
  info <- file.info(fs)
  for (i in seq_along(fs)) {
    wl(sprintf("    %-52s %8.1f KB  %s", basename(fs[i]),
               info$size[i] / 1024, format(info$mtime[i], "%Y-%m-%d %H:%M")))
  }
}
wl("")
wl("SESSION")
writeLines(capture.output(sessionInfo()), con)
close(con)

cat("  -> manifest.txt\n")
cat("\nAppendices written to ", normalizePath(DIR_APPEND), "\n", sep = "")
