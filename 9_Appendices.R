# =============================================================================
# 9_Appendices.R - FIGURES 2, 3 AND 4, SUPPLEMENTARY FIGURES S1-S2 AND S7,
#                  APPENDIX A TABLES AND THE RUN MANIFEST
# -----------------------------------------------------------------------------
# Run after 6b, 6c and 6d. Reads the saved results (read_all_runs) and
# dat_classed. File names follow the manuscript numbering.
#
#   Output_Figures/
#     Fig2_sets_ecoregions         trawl sets by year and pooled, four ecoregions
#     Fig3a_prey_categories        prey categories by aggregation threshold
#     Fig3b_sweep_families         family frequencies across the sweep
#     Fig4_typology_tree           decision tree of the typology (svg + png)
#   Output_Appendices/
#     FigS1_sampling_by_year       stomachs per year and ecoregion
#     FigS2_sampling_by_stratum    stomachs per stratum and year
#     FigS7_lowN_sensitivity       composites with / without sparse units
#     TableA0 prey categories per threshold   TableA1 sets per cell
#     TableA2 effect size vs detection        TableA3 stratum / ecoregion crosswalk
#     TableA4 unclassifiable cells            TableA5 spread across resolutions
#     TableA6 low-n sensitivity               TableA7 identification depth by period
#     manifest.txt                            session, parameters, file inventory
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(sf); library(forcats)
})

PREY_FAMILY <- "_1"
source("R_helpers/Config_Mappings.R")

DATA_PATH   <- "data/dat_classed.rda"
STRATA_PATH <- "strata_rv_gulf.rds"
ECO_PATH    <- "data/Spatial_data/new_gulf_final.shp"
MAP_XLIM    <- c(-66.2, -60.0)
MAP_YLIM    <- c(45.5, 49.2)
CRS_MAP     <- "+proj=lcc +lat_1=46 +lat_2=49 +lat_0=47 +lon_0=-63 +datum=WGS84 +units=m"
N_FLAG      <- 5

save_app <- function(p, stem, w, h, dpi = 300) save_fig(p, stem, w, h, dpi, dir = DIR_APPEND)

emit_app <- function(df, stem, caption) {
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
  if (inherits(raw, "sf")) raw <- sf::st_drop_geometry(raw)
  raw <- as.data.frame(raw)
  cat("Raw data: ", nrow(raw), " prey records.\n", sep = "")
} else {
  raw <- NULL
  message("Raw data not found at ", DATA_PATH,
          "; Figures 2 and 3 and the crosswalk are skipped.")
}

# =============================================================================
# 2. FIGURE 2 - TRAWL SETS BY YEAR IN THE FOUR ECOREGIONS (metric projection)
# =============================================================================
if (have_raw && file.exists(ECO_PATH) && all(c("longitude", "latitude") %in% names(raw))) {
  cat("\nFigure 2: trawl sets and ecoregions\n")

  eco <- st_make_valid(st_read(ECO_PATH, quiet = TRUE)) %>% mutate(Area = area_label(Area))
  sets_xy <- raw %>%
    filter(year %in% c(2004:2006, 2018:2019), !is.na(longitude), !is.na(latitude)) %>%
    distinct(year, vessel.code, set, longitude, latitude, Area) %>%
    mutate(Area = area_label(Area))
  panels <- c(sort(unique(as.character(sets_xy$year))), "All years")
  pts <- bind_rows(sets_xy %>% mutate(panel = as.character(year)),
                   sets_xy %>% mutate(panel = "All years")) %>%
    mutate(panel = factor(panel, levels = panels)) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

  land <- NULL
  if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    land <- tryCatch(rnaturalearth::ne_states(country = c("canada", "united states of america"),
                                              returnclass = "sf"),
                     error = function(e) tryCatch(
                       rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"),
                       error = function(e2) NULL))
  }
  bb <- st_bbox(st_transform(
    st_as_sfc(st_bbox(c(xmin = MAP_XLIM[1], xmax = MAP_XLIM[2],
                        ymin = MAP_YLIM[1], ymax = MAP_YLIM[2]), crs = st_crs(4326))),
    CRS_MAP))

  p2 <- ggplot()
  if (!is.null(land)) p2 <- p2 + geom_sf(data = land, fill = "grey90", colour = "grey70", linewidth = 0.2)
  p2 <- p2 +
    geom_sf(data = eco, aes(fill = Area), alpha = 0.25, colour = "grey30", linewidth = 0.3) +
    geom_sf(data = pts, aes(colour = Area), size = 0.6, alpha = 0.75, show.legend = FALSE) +
    facet_wrap(~panel, ncol = 3) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    scale_colour_brewer(palette = "Set2") +
    scale_x_continuous(breaks = seq(-65, -61, by = 1)) +
    scale_y_continuous(breaks = seq(46, 49, by = 1)) +
    coord_sf(crs = CRS_MAP, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
             expand = FALSE) +
    labs(title = "Trawl sets with analysed stomachs, by survey year and pooled",
         x = NULL, y = NULL) +
    theme_diag(base_size = 10) +
    theme(axis.text = element_text(size = 6), panel.grid = element_line(colour = "grey93"),
          panel.spacing.x = unit(1.2, "lines"), legend.position = "bottom")
  if (requireNamespace("ggspatial", quietly = TRUE))
    p2 <- p2 + ggspatial::annotation_scale(
      data = data.frame(panel = factor("All years", levels = panels)),
      location = "br", width_hint = 0.25, text_cex = 0.6)

  save_fig(p2, "Fig2_sets_ecoregions", width = 11, height = 7.6, dir = DIR_FIGURES)
}

# =============================================================================
# 2a. FIGURE 4 - DECISION TREE OF THE TYPOLOGY (SVG + PNG, R_helpers/Fig4_typology_tree.R)
# =============================================================================
cat("\nFigure 4: decision tree\n")
source("R_helpers/Fig4_typology_tree.R")
draw_typology_tree(out_prefix = file.path(DIR_FIGURES, "Fig4_typology_tree"))

# =============================================================================
# 2b. FIGURES S1-S2 - SAMPLING (stomachs, not prey records)
# =============================================================================
if (have_raw) {
  cat("\nFigures S1-S2: sampling\n")

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
    labs(title = "Stomachs analysed per year and ecoregion",
         x = NULL, y = "Stomachs") +
    theme_diag()

  save_app(pA, "FigS1_sampling_by_year", 8, 4.6)

  strata_sf <- load_strata(STRATA_PATH)
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
      labs(title = "Stomachs analysed per stratum and survey year") +
      theme_diag(base_size = 10) +
      theme(axis.title = element_blank(), axis.text = element_text(size = 5.5),
            strip.text = element_text(face = "bold"),
            panel.grid = element_line(colour = "grey93"))

    save_app(pB, "FigS2_sampling_by_stratum", 14, 4.2, dpi = 200)
  }
}

# =============================================================================
# 3. FIGURE 3 - TAXONOMIC RESOLUTION SWEEP (A: categories, B: families)
# =============================================================================
# =============================================================================
# 3. FIGURE 3A - PREY CATEGORIES RETAINED (Bar plot & refined axes)
# =============================================================================
cat("\nFigure 3A: prey categories bar plot\n")

if (have_raw) {
  prey_cols <- grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"), names(raw), value = TRUE)
  if (length(prey_cols)) {
    ncat <- tibble(
      col = prey_cols,
      x_threshold = as.numeric(str_extract(prey_cols, "\\d+")),
      n_cat = vapply(prey_cols, function(cc) dplyr::n_distinct(raw[[cc]], na.rm = TRUE),
                     numeric(1))
    ) %>% arrange(x_threshold)

    # Graphique en barres avec des échelles plus lisibles (pas de 25 pour y, 100 pour x)
    p3a <- ggplot(ncat, aes(x = x_threshold, y = n_cat)) +
      geom_col(fill = "#2F4A5A", width = 7, alpha = 0.85) + # geom_col pour un diagramme en barres
      scale_x_continuous(
        breaks = seq(0, 1000, by = 50),   # Graduations tous les 100 sur l'axe X
        limits = c(-10, 1020)
      ) +
      scale_y_continuous(
        breaks = seq(0, 200, by = 10),     # Graduations tous les 25 sur l'axe Y (plus lisible que 50)
        limits = c(0, 190)
      ) +
      labs(title = "A. Prey categories retained by aggregation threshold",
           subtitle = "Each record is kept at its finest level, or walked up the taxonomic tree until the node reaches x stomachs",
           x = "Threshold x (stomachs)", y = "Distinct prey categories") +
      theme_diag() +
      theme(
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid.minor = element_blank()
      )

    save_fig(p3a, "Fig3a_prey_categories", 9, 5, dir = DIR_FIGURES)
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

save_fig(p3b, "Fig3b_sweep_families", 13, 5.6, dir = DIR_FIGURES)



# =============================================================================
# FIGURE 3b - SWEEP FAMILIES: ECOREGIONS
# =============================================================================
cat("\nFigure 3b: taxonomic resolution sweep - Ecoregions\n")

sweep_eco <- res_all %>%
  filter(level == "ecoregion", diagnostic %in% DIAG_LEVELS, contrast != "2006 vs 2018") %>%
  count(spatial_unit, currency, contrast, x_threshold, diagnostic, name = "n") %>%
  group_by(spatial_unit, currency, contrast, x_threshold) %>%
  mutate(pct_unit = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(family = family_of(diagnostic)) %>%
  filter(!is.na(family)) %>%
  group_by(spatial_unit, currency, contrast, x_threshold, family) %>%
  summarise(pct_unit = sum(pct_unit), .groups = "drop") %>%
  group_by(currency, contrast, x_threshold, family) %>%
  summarise(pct = mean(pct_unit), .groups = "drop")

p3b_eco <- ggplot(sweep_eco, aes(x_threshold, pct, colour = family)) +
  geom_line(linewidth = 0.6) +
  facet_grid(currency ~ contrast,
             labeller = labeller(currency = CURRENCY_LAB, contrast = CONTRAST_LAB_1L)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity of diagnostic composition to taxonomic resolution — Ecoregions",
       subtitle = "Family frequencies across the full sweep (mean across ecoregions)",
       x = "Threshold x (stomachs)", y = "Frequency (%)") +
  theme_diag() +
  theme(strip.text = element_text(size = 8), legend.position = "top")

save_fig(p3b_eco, "Fig3b_sweep_families_ecoregion", 13, 5.6, dir = DIR_FIGURES)



# =============================================================================
# FIGURE 3b - SWEEP FAMILIES: SURVEY STRATA
# =============================================================================
cat("\nFigure 3b: taxonomic resolution sweep - Survey Strata\n")

sweep_strata <- res_all %>%
  filter(level == "stratum", diagnostic %in% DIAG_LEVELS, contrast != "2006 vs 2018") %>%
  count(spatial_unit, currency, contrast, x_threshold, diagnostic, name = "n") %>%
  group_by(spatial_unit, currency, contrast, x_threshold) %>%
  mutate(pct_unit = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(family = family_of(diagnostic)) %>%
  filter(!is.na(family)) %>%
  group_by(spatial_unit, currency, contrast, x_threshold, family) %>%
  summarise(pct_unit = sum(pct_unit), .groups = "drop") %>%
  group_by(currency, contrast, x_threshold, family) %>%
  summarise(pct = mean(pct_unit), .groups = "drop")

p3b_strata <- ggplot(sweep_strata, aes(x_threshold, pct, colour = family)) +
  geom_line(linewidth = 0.6) +
  facet_grid(currency ~ contrast,
             labeller = labeller(currency = CURRENCY_LAB, contrast = CONTRAST_LAB_1L)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity of diagnostic composition to taxonomic resolution — Survey Strata",
       subtitle = "Family frequencies across the full sweep (mean across strata)",
       x = "Threshold x (stomachs)", y = "Frequency (%)") +
  theme_diag() +
  theme(strip.text = element_text(size = 8), legend.position = "top")

save_fig(p3b_strata, "Fig3b_sweep_families_stratum", 13, 5.6, dir = DIR_FIGURES)



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
# 5. PARTITION CROSSWALK (strata split across ecoregions)
# =============================================================================
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
# 6. FIGURE S7 - SENSITIVITY TO LOW-SAMPLE UNITS (three scales, with intervals)
# =============================================================================
# Across-unit composites with and without the units built on fewer than N_FLAG
# predator x size-class cells; Gulf-wide (single unit) is the reference column.
# Composite computed within each resolution, then summarised across resolutions:
#   S5_BAND = "q95"  empirical 2.5-97.5 % interval across resolutions
#   S5_BAND = "sd"   mean +/- SD across resolutions
cat("\nFigure S7: low-sample sensitivity\n")

S5_BAND <- "q95"

res_s7 <- res_all %>%
  mutate(spatial_unit = if_else(level == "all_gulf", "Gulf-wide", spatial_unit))

unit_keys <- c("currency", "contrast", "level", "spatial_unit")

# cells (predator x size-class) behind each unit, for the >= N_FLAG filter
cells_per_unit <- n_units_by(res_s7, keys = unit_keys) %>%
  rename(n_cells = n_units)

# family percentages per unit AND per resolution (absent family = 0)
fam_ux <- family_by_resolution(res_s7, keys = unit_keys) %>%
  left_join(cells_per_unit, by = unit_keys) %>%
  filter(!is.na(contrast))

sub_lab <- paste0("n cells >= ", N_FLAG)

comp_by_x <- bind_rows(
  fam_ux %>% mutate(subset = "All units"),
  fam_ux %>%
    filter(level != "all_gulf", !is.na(n_cells), n_cells >= N_FLAG) %>%
    mutate(subset = sub_lab)
) %>%
  # composite for one resolution = mean over the units retained
  group_by(level, currency, contrast, subset, family, x_threshold) %>%
  summarise(pct = mean(pct), n_units_kept = n_distinct(spatial_unit),
            .groups = "drop")

sens <- comp_by_x %>%
  group_by(level, currency, contrast, subset, family) %>%
  summarise(pct_mean = mean(pct), pct_sd = sd(pct),
            q_lo = unname(quantile(pct, 0.025)), q_hi = unname(quantile(pct, 0.975)),
            n_res = n_distinct(x_threshold), n_units_kept = min(n_units_kept),
            .groups = "drop") %>%
  mutate(lo = if (S5_BAND == "sd") pmax(pct_mean - pct_sd, 0) else q_lo,
         hi = if (S5_BAND == "sd") pmin(pct_mean + pct_sd, 100) else q_hi,
         level    = factor(level, levels = SPATIAL_LEVELS_ORD),
         subset   = factor(subset, levels = c("All units", sub_lab)),
         family   = factor(family, levels = FAMILY_LEVELS),
         contrast = factor(contrast, levels = CONTRAST_LEVELS),
         # numeric x with a small offset by subset so the intervals do not overlap
         x = as.integer(contrast) + if_else(subset == "All units", -0.07, 0.07))

band_lab <- if (S5_BAND == "sd") "mean +/- SD across taxonomic resolutions" else
  "95 % interval (2.5-97.5 %) across taxonomic resolutions"

pS7 <- ggplot(sens, aes(x = x, y = pct_mean, colour = family,
                        shape = subset, group = interaction(family, subset))) +
  geom_errorbar(aes(ymin = lo, ymax = hi, linetype = subset),
                width = 0.1, linewidth = 0.35, alpha = 0.85) +
  geom_line(aes(linetype = subset), linewidth = 0.6) +
  geom_point(size = 2.2, fill = "white") +
  scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS), labels = unname(CONTRAST_LAB),
                     limits = c(0.5, length(CONTRAST_LEVELS) + 0.5), expand = c(0, 0)) +
  facet_grid(currency ~ level,
             labeller = labeller(currency = CURRENCY_LAB, level = LEVEL_SHORT)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_shape_manual(values = c(16, 21), name = NULL, drop = FALSE) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Sensitivity of the across-unit composites to sparsely sampled units",
       subtitle = paste0("Solid = every unit; dashed = units built on at least ", N_FLAG,
                         " predator x size-class cells. Bars: ", band_lab,
                         ".\nGulf-wide is a single unit and is shown as the reference."),
       x = NULL, y = "Frequency (%)") +
  theme_diag() +
  theme(legend.position = "top")

save_app(pS7, "FigS7_lowN_sensitivity", 13, 6.5)

emit_app(
  sens %>%
    select(level, currency, contrast, subset, family, pct_mean, lo, hi,
           n_res, n_units_kept) %>%
    mutate(across(c(pct_mean, lo, hi), ~round(.x, 1))) %>%
    pivot_wider(names_from = family, values_from = c(pct_mean, lo, hi),
                names_glue = "{family}_{.value}"),
  "TableA6_lowN_sensitivity",
  paste0("Across-unit composites with and without sparsely sampled units; ",
         band_lab, "."))

# =============================================================================
# 7. TABLE A7 - DEPTH OF PREY IDENTIFICATION BY PERIOD
# =============================================================================
# From tax_level_lowest (lowest recorded rank of each prey item): % of items
# whose lowest rank is each rank (all items), and cumulative % identified to
# that rank or finer, computed over the items carrying a taxonomic rank
# (life-stage prey categories excluded from that denominator).
if (!is.null(raw) && "tax_level_lowest" %in% names(raw)) {
  cat("\nTable A7: identification depth by period\n")

  e2 <- new.env(); sys.source("R_helpers/taxonomic_rank_order.R", envir = e2)
  rank_order <- get("ranknfile", envir = e2)

  # Sub-, super- and infra-ranks are counted with their parent rank.
  parent_rank <- function(r) {
    r <- tolower(as.character(r))
    out <- rep(NA_character_, length(r))
    for (p in c("phylum", "class", "order", "family", "genus", "species")) {
      hit <- grepl(paste0("^(super|sub|infra|parv|mega|giga|magn|grand)?", p, "$"), r)
      out[hit & is.na(out)] <- p
    }
    out
  }

  id <- raw %>%
    mutate(period = case_when(year %in% 2004:2006 ~ "2004-2006",
                              year %in% 2018:2019 ~ "2018-2019",
                              TRUE ~ NA_character_)) %>%
    filter(!is.na(period), !is.na(tax_level_lowest)) %>%
    transmute(period,
              rank_raw = tolower(as.character(tax_level_lowest)),
              rank     = coalesce(parent_rank(tax_level_lowest), "life-stage prey category"),
              depth    = match(rank_raw, rank_order))    # position in the full rank list

  rank_lv <- c("phylum", "class", "order", "family", "genus", "species",
               "life-stage prey category")
  i_of <- function(p) match(p, rank_order)

  at_rank <- id %>%
    count(period, rank, name = "n") %>%
    group_by(period) %>%
    mutate(pct_at = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    select(-n)

  reach <- id %>%
    filter(!is.na(depth)) %>%
    group_by(period) %>%
    summarise(n_items = n(),
              class   = round(100 * mean(depth >= i_of("class"),   na.rm = TRUE), 1),
              order   = round(100 * mean(depth >= i_of("order"),   na.rm = TRUE), 1),
              family  = round(100 * mean(depth >= i_of("family"),  na.rm = TRUE), 1),
              genus   = round(100 * mean(depth >= i_of("genus"),   na.rm = TRUE), 1),
              species = round(100 * mean(depth >= i_of("species"), na.rm = TRUE), 1),
              .groups = "drop") %>%
    pivot_longer(c(class, order, family, genus, species),
                 names_to = "rank", values_to = "pct_at_or_finer")

  tabA7 <- expand_grid(rank = rank_lv, period = unique(id$period)) %>%
    left_join(at_rank, by = c("period", "rank")) %>%
    left_join(reach %>% select(period, rank, pct_at_or_finer), by = c("period", "rank")) %>%
    mutate(pct_at = coalesce(pct_at, 0)) %>%
    pivot_wider(names_from = period, values_from = c(pct_at, pct_at_or_finer),
                names_glue = "{.value}_{period}") %>%
    mutate(rank = factor(rank, levels = rank_lv)) %>%
    arrange(rank)

  n_items <- id %>% count(period, name = "n_items")
  cat("Prey items: ", paste(n_items$period, n_items$n_items, collapse = " | "), "\n", sep = "")

  emit_app(tabA7, "TableA7_identification_depth",
           paste0("Depth of prey identification by period (% of prey items whose lowest ",
                  "recorded rank is each rank, and % identified to that rank or finer); n = ",
                  paste(n_items$n_items, collapse = " / "), " items"))
} else {
  message("tax_level_lowest not in dat_classed: Table A7 skipped.")
}

# =============================================================================
# 8. MANIFEST
# =============================================================================
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
for (d in c(DIR_FIGURES, DIR_TABLES, DIR_APPEND)) {
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

