# =============================================================================
# 7_Figures.R - MAIN FIGURES 5-10 AND SUPPLEMENTARY S3-S4
# -----------------------------------------------------------------------------
# Run after 6b, 6c and 6d. Nothing is recomputed: the saved results hold the
# per-cell diagnostics; this script aggregates (freq_by_resolution: percentages
# within each resolution, then the mean across resolutions; Inconclusive cells
# excluded) and draws. File names follow the manuscript numbering.
#
#   Fig5_regime_break_gulf        families by contrast, Gulf-wide
#   Fig6_regime_break_ecoregion   composition within each ecoregion
#   Fig7_regime_break_ecoregion   families by contrast, mean across ecoregions
#   Fig8_regime_break_stratum     families by contrast, mean across strata
#   Fig9_map_strata_biomass       family frequency mapped on strata, biomass
#   Fig10_crossscale_gradient     stratum -> ecoregion -> Gulf, between periods
#   FigS3_regime_break_by_ecoregion   as Fig 5, one panel per ecoregion (2 x 2)
#   FigS4_map_strata_occurrence   as Fig 9, occurrence
#   FigS5/FigS6_map_strata_*_masked  as Fig 9 / S4, strata below N_FLAG units greyed
#   data_Fig*.csv                 the numbers behind each figure
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(forcats)
  library(sf); library(readr); library(stringr)
})

PREY_FAMILY <- "_1"
source("R_helpers/Config_Mappings.R")

N_FLAG        <- 5                              # flag units built on fewer cells
STRATA_PATH   <- "strata_rv_gulf.rds"           # sf object, or the gulf.spatial package
COAST_PATH    <- "gulf.coast.intermediate.csv"  # optional coastline (x, y, pid)
MAP_XLIM      <- c(-66.2, -60.0)
MAP_YLIM      <- c(45.5, 49.2)
MAP_CONTRASTS <- c("P1a", "P1b", "P2", "PTa")   # 2006 vs 2018 is not estimable by stratum

# =============================================================================
# 1. READ AND AGGREGATE
# =============================================================================
cat("\nReading pipeline results...\n")
res_all <- read_all_runs()

cat("Runs read:\n")
print(res_all %>% count(level, currency, contrast) %>% arrange(level, currency, contrast),
      n = Inf)

# Gulf-wide and per-unit family tables.
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

# Mean family composition at each level (mean across units for ecoregion and
# stratum).
fam_level <- bind_rows(
  fam_gulf %>% mutate(level = "all_gulf", .before = 1),
  fam_unit %>%
    group_by(currency, contrast, level) %>%
    summarise(across(all_of(FAMILY_LEVELS), ~mean(.x, na.rm = TRUE)),
              .groups = "drop")
) %>%
  mutate(level = factor(level, levels = SPATIAL_LEVELS_ORD))

long_of <- function(d) {
  d %>%
    pivot_longer(all_of(FAMILY_LEVELS), names_to = "family", values_to = "pct") %>%
    mutate(family   = factor(family, levels = FAMILY_LEVELS),
           currency = factor(currency, levels = names(CURRENCY_LAB)))
}

# =============================================================================
# 2. FIGURE 5 - GULF-WIDE REGIME BREAK
# =============================================================================
# Lines per family, linetype per currency; error bars = SD across resolutions.
cat("\nFigure 5: Gulf-wide regime break\n")

fam_gulf_sd <- spread_by_resolution(filter(res_all, level == "all_gulf"),
                                    keys = c("currency", "contrast")) %>%
  select(currency, contrast, family, sd = pct_sd) %>%
  mutate(currency = factor(currency, levels = names(CURRENCY_LAB)))

d5 <- long_of(fam_gulf) %>%
  filter(!is.na(contrast)) %>%
  left_join(fam_gulf_sd, by = c("currency", "contrast", "family")) %>%
  mutate(x = as.integer(contrast) + ifelse(currency == "biomass", -0.06, 0.06))

fig5 <- ggplot(d5, aes(x, pct, colour = family, shape = family,
                       linetype = currency, group = interaction(family, currency))) +
  decade_bands() +
  annotate("text", x = 2,   y = Inf, label = "WITHIN PERIOD",
           colour = "#2F4A5A", fontface = 2, size = 3.6, vjust = 1.8) +
  annotate("text", x = 4.5, y = Inf, label = "BETWEEN PERIODS",
           colour = "#9B2226", fontface = 2, size = 3.6, vjust = 1.8) +
  geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                width = 0.08, linewidth = 0.4, linetype = "solid",
                alpha = 0.7, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.6, fill = "white") +
  scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS),
                     labels = unname(CONTRAST_LAB[CONTRAST_LEVELS]),
                     expand = expansion(add = 0.25)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_shape_manual(values = FAMILY_SHAPE, name = NULL) +
  scale_linetype_manual(values = c(biomass = "solid", occurrence = "22"),
                        labels = CURRENCY_LAB, name = NULL) +
  labs(
    title = "Diet regime break at the Gulf-wide scale",
    subtitle = paste0("Frequency of the four diagnostic families among predator x ",
                      "size-class cells: mean and\nstandard deviation across ~100 prey taxonomic resolutions"),
    x = NULL, y = "Frequency (%)"
  ) +
  theme_diag() +
  theme(legend.box = "horizontal")

save_fig(fig5, "Fig5_regime_break_gulf", width = 9, height = 5.6)

# =============================================================================
# 3. FIGURE 6 - COMPOSITION BY ECOREGION
# =============================================================================
cat("\nFigure 6: by ecoregion\n")

d6 <- long_of(filter(fam_unit, level == "ecoregion")) %>%
  filter(!is.na(contrast)) %>%
  mutate(unit_lab = area_label(spatial_unit))
d6_n <- d6 %>% distinct(currency, unit_lab, contrast, n_units)

fig6 <- ggplot(d6, aes(contrast, pct, fill = family)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
  geom_text(data = d6_n, aes(contrast, 101, label = paste0("n=", n_units)),
            inherit.aes = FALSE, size = 2.3, vjust = 0, colour = "grey30") +
  scale_x_discrete(labels = CONTRAST_LAB) +
  geom_hline(yintercept = c(25, 50, 75), colour = "white",
             linewidth = 0.3, linetype = "dotted") +
  facet_grid(currency ~ unit_lab,
             labeller = labeller(currency = CURRENCY_LAB)) +
  scale_fill_manual(values = FAMILY_PAL, name = NULL) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Diet regime break by ecoregion",
    subtitle = paste0("Composition of the four families within each ecoregion. ",
                      "n = predator x size-class units contributing to each bar."),
    x = NULL, y = "Within-ecoregion frequency (%)"
  ) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        panel.grid.major.x = element_blank())

save_fig(fig6, "Fig6_regime_break_ecoregion", width = 11, height = 6.4)

# =============================================================================
# 4. FIGURES 7-8 - REGIME BREAK AT THE ECOREGION AND STRATUM SCALES
# =============================================================================
# Same grammar as Figure 5 on the across-unit mean; error bars = SD across
# resolutions of that mean.
cat("\nFigures 7-8: ecoregion and stratum scales\n")

rb_sd <- spread_by_resolution_units(filter(res_all, level != "all_gulf"),
                                    keys = c("currency", "contrast", "level")) %>%
  transmute(currency = factor(currency, levels = names(CURRENCY_LAB)),
            contrast, level = as.character(level), family, sd = pct_sd)

make_regime_break <- function(lvl, stem, title_txt, subtitle_txt) {
  d <- long_of(filter(fam_level, level == lvl)) %>%
    filter(!is.na(contrast)) %>%
    mutate(level = as.character(level)) %>%
    left_join(rb_sd, by = c("currency", "contrast", "level", "family")) %>%
    mutate(x = as.integer(contrast) + ifelse(currency == "biomass", -0.06, 0.06))

  g <- ggplot(d, aes(x, pct, colour = family, shape = family,
                     linetype = currency,
                     group = interaction(family, currency))) +
    decade_bands() +
    annotate("text", x = 2,   y = Inf, label = "WITHIN PERIOD",
             colour = "#2F4A5A", fontface = 2, size = 3.6, vjust = 1.8) +
    annotate("text", x = 4.5, y = Inf, label = "BETWEEN PERIODS",
             colour = "#9B2226", fontface = 2, size = 3.6, vjust = 1.8) +
    geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                  width = 0.08, linewidth = 0.4, linetype = "solid",
                  alpha = 0.7, show.legend = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.6, fill = "white") +
    scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS),
                       labels = unname(CONTRAST_LAB[CONTRAST_LEVELS]),
                       expand = expansion(add = 0.25)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_colour_manual(values = FAMILY_PAL, name = NULL) +
    scale_shape_manual(values = FAMILY_SHAPE, name = NULL) +
    scale_linetype_manual(values = c(biomass = "solid", occurrence = "22"),
                          labels = CURRENCY_LAB, name = NULL) +
    labs(title = title_txt, subtitle = subtitle_txt, x = NULL, y = "Frequency (%)") +
    theme_diag() +
    theme(legend.box = "horizontal")

  save_fig(g, stem, width = 9, height = 5.6)
}

make_regime_break(
  "ecoregion", "Fig7_regime_break_ecoregion",
  "Diet regime break at the ecoregion scale",
  paste0("Family frequencies averaged across the four ecoregions: mean and\n",
         "standard deviation across ~100 prey taxonomic resolutions"))

make_regime_break(
  "stratum", "Fig8_regime_break_stratum",
  "Diet regime break at the stratum scale",
  paste0("Family frequencies averaged across survey strata: mean and\n",
         "standard deviation across ~100 prey taxonomic resolutions"))

# =============================================================================
# 5. FIGURE S3 - REGIME BREAK WITHIN EACH ECOREGION
# =============================================================================
# Same grammar as Figure 5, one panel per ecoregion; error bars = SD across
# resolutions within the ecoregion.
cat("\nFigure S3: regime break by ecoregion\n")

unit_sd <- spread_by_resolution(filter(res_all, level != "all_gulf"),
                                keys = c("currency", "contrast", "level", "spatial_unit")) %>%
  transmute(currency = factor(currency, levels = names(CURRENCY_LAB)),
            contrast, level = as.character(level), spatial_unit, family, sd = pct_sd)

make_unit_regime_break <- function(lvl, stem, title_txt, ncol_panels, width, height,
                                   base_size = 10, point_size = 1.8) {
  d <- long_of(filter(fam_unit, level == lvl)) %>%
    filter(!is.na(contrast)) %>%
    mutate(level = as.character(level)) %>%
    left_join(unit_sd, by = c("currency", "contrast", "level", "spatial_unit", "family")) %>%
    mutate(unit_lab = if (lvl == "ecoregion") area_label(spatial_unit) else paste("Stratum", spatial_unit),
           x = as.integer(contrast) + ifelse(currency == "biomass", -0.08, 0.08))
  n_cells <- d %>% group_by(unit_lab) %>% summarise(n = max(n_units, na.rm = TRUE), .groups = "drop")
  d <- d %>% left_join(n_cells, by = "unit_lab") %>%
    mutate(unit_lab = paste0(unit_lab, " (n = ", n, ")"))
  if (lvl == "stratum") d <- d %>% mutate(unit_lab = fct_reorder(unit_lab, as.numeric(spatial_unit)))

  g <- ggplot(d, aes(x, pct, colour = family, shape = family,
                     linetype = currency, group = interaction(family, currency))) +
    decade_bands() +
    geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                  width = 0.1, linewidth = 0.3, linetype = "solid",
                  alpha = 0.7, show.legend = FALSE) +
    geom_line(linewidth = 0.6) +
    geom_point(size = point_size, fill = "white") +
    facet_wrap(~unit_lab, ncol = ncol_panels) +
    scale_x_continuous(breaks = seq_along(CONTRAST_LEVELS),
                       labels = unname(CONTRAST_LAB[CONTRAST_LEVELS]),
                       expand = expansion(add = 0.25)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_colour_manual(values = FAMILY_PAL, name = NULL) +
    scale_shape_manual(values = FAMILY_SHAPE, name = NULL) +
    scale_linetype_manual(values = c(biomass = "solid", occurrence = "22"),
                          labels = CURRENCY_LAB, name = NULL) +
    labs(title = title_txt,
         subtitle = paste0("Family frequencies within each unit: mean and SD across ~100 prey resolutions.\n",
                           "Bands: within-period (left) and between-period (right) contrasts; ",
                           "n = predator x size-class cells."),
         x = NULL, y = "Frequency (%)") +
    theme_diag(base_size = base_size) +
    theme(legend.position = "bottom", legend.box = "horizontal",
          axis.text.x = element_text(size = base_size - 3),
          strip.text = element_text(size = base_size - 1))
  save_fig(g, stem, width = width, height = height)
}

if (any(fam_unit$level == "ecoregion"))
  make_unit_regime_break("ecoregion", "FigS3_regime_break_by_ecoregion",
                         "Diet regime break within each ecoregion",
                         ncol_panels = 2, width = 11, height = 8, base_size = 11, point_size = 2.2)

# =============================================================================
# 6. FIGURE 10 - CROSS-SCALE GRADIENT (between-period contrast)
# =============================================================================
cat("\nFigure 10: cross-scale gradient\n")

grad_sd <- bind_rows(
  spread_by_resolution(filter(res_all, level == "all_gulf"),
                       keys = c("currency", "contrast")) %>%
    mutate(level = "all_gulf"),
  spread_by_resolution_units(filter(res_all, level != "all_gulf"),
                             keys = c("currency", "contrast", "level"))
) %>%
  transmute(currency = factor(currency, levels = names(CURRENCY_LAB)),
            contrast, level = as.character(level), family, sd = pct_sd)

d_grad <- long_of(fam_level) %>%
  filter(contrast == "PTa") %>%
  mutate(level = as.character(level)) %>%
  left_join(grad_sd, by = c("currency", "contrast", "level", "family")) %>%
  mutate(level = factor(level, levels = rev(SPATIAL_LEVELS_ORD)))

fig_grad <- ggplot(d_grad, aes(level, pct, colour = family, group = family)) +
  geom_errorbar(aes(ymin = pmax(pct - sd, 0), ymax = pmin(pct + sd, 100)),
                width = 0.1, linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = family), size = 3.2) +
  geom_text(aes(label = sprintf("%.1f", pct)), vjust = -1.1,
            size = 2.9, show.legend = FALSE) +
  facet_wrap(~currency, labeller = labeller(currency = CURRENCY_LAB)) +
  scale_colour_manual(values = FAMILY_PAL, name = NULL) +
  scale_shape_manual(values = FAMILY_SHAPE, name = NULL) +
  scale_x_discrete(labels = c(stratum = "Stratum", ecoregion = "Ecoregion",
                              all_gulf = "Gulf-wide")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = paste0("Effect of spatial aggregation on the diagnostic composition (",
                   CONTRAST_LAB_1L[["PTa"]], ")"),
    subtitle = paste0("Stratum and ecoregion values are means across units; the ",
                      "Gulf-wide value is a single pooled test.\nBars = SD across resolutions. Aggregation ",
                      "raises apparent Reorganisation and lowers Stability, ",
                      "while Substitution is scale-invariant."),
    x = "Spatial unit, from finest to coarsest", y = "Frequency (%)"
  ) +
  theme_diag()

save_fig(fig_grad, "Fig10_crossscale_gradient", width = 9.5, height = 5.4)

# =============================================================================
# 7. FIGURE 9 AND S4 - STRATUM MAPS
# =============================================================================
cat("\nFigures 9 and S4: stratum maps\n")

strata <- load_strata(STRATA_PATH)

if (is.null(strata)) {
  message("Stratum geometry not found (neither gulf.spatial nor ", STRATA_PATH,
          "). Maps skipped; every other figure is written.")
} else {
  strata <- sf::st_make_valid(strata)
  if (is.na(sf::st_crs(strata))) sf::st_crs(strata) <- 4326

  lab_xy <- if (all(c("label_x", "label_y") %in% names(strata))) {
    strata %>% sf::st_drop_geometry() %>% transmute(spatial_unit = as.character(str),
                                                    X = label_x, Y = label_y)
  } else {
    cc <- suppressWarnings(sf::st_coordinates(sf::st_centroid(strata)))
    tibble(spatial_unit = as.character(strata$str), X = cc[, 1], Y = cc[, 2])
  }

  coast <- if (file.exists(COAST_PATH)) {
    read_csv(COAST_PATH, show_col_types = FALSE) %>% filter(!is.na(x), !is.na(y))
  } else NULL

  map_dat <- long_of(filter(fam_unit, level == "stratum")) %>%
    filter(!is.na(contrast), contrast %in% MAP_CONTRASTS) %>%
    mutate(low_n = !is.na(n_units) & n_units < N_FLAG)

  strata$spatial_unit <- as.character(strata$str)

  # mask_low_n = FALSE : every stratum coloured, low-n values starred (Fig 9, S4)
  # mask_low_n = TRUE  : strata below N_FLAG units greyed, no value (Fig S5, S6)
  render_map <- function(currencies, facet_formula, title, stem,
                         ncol_panels, nrow_panels, mask_low_n = FALSE) {
    d <- strata %>%
      left_join(filter(map_dat, currency %in% currencies), by = "spatial_unit") %>%
      filter(!is.na(pct)) %>%
      mutate(contrast = droplevels(factor(contrast, levels = CONTRAST_LEVELS)),
             currency = droplevels(factor(currency, levels = names(CURRENCY_LAB))))
    d_col  <- if (mask_low_n) filter(d, !low_n) else d
    d_grey <- if (mask_low_n) filter(d, low_n) else d[0, ]

    lab <- d_col %>%
      sf::st_drop_geometry() %>%
      left_join(lab_xy, by = "spatial_unit") %>%
      transmute(family, contrast, currency, X, Y,
                lab = ifelse(low_n & !mask_low_n, paste0(sprintf("%.0f", pct), "*"),
                             sprintf("%.0f", pct)),
                txt = ifelse(pct < 50, "grey95", "grey10"))

    subtitle <- if (mask_low_n) {
      paste0("Grey strata: fewer than ", N_FLAG, " predator x size-class units ",
             "(value not shown). 2006 vs 2018 is not estimable at this scale.")
    } else {
      paste0("Every stratum shown; * = built on < ", N_FLAG,
             " predator x size-class units (interpret with caution). ",
             "2006 vs 2018 is not estimable at this scale.")
    }

    g <- ggplot()
    if (!is.null(coast)) {
      g <- g + geom_polygon(data = coast, aes(x, y, group = pid),
                            fill = "grey88", colour = NA)
    }
    g <- g +
      geom_sf(data = strata, fill = "grey97", colour = "grey70", linewidth = 0.15) +
      geom_sf(data = d_col, aes(fill = pct), colour = "grey45", linewidth = 0.15)
    if (nrow(d_grey)) {
      g <- g + geom_sf(data = d_grey, fill = "grey78", colour = "grey45", linewidth = 0.15)
    }
    g <- g +
      geom_text(data = lab, aes(X, Y, label = lab, colour = txt),
                size = 2.1, fontface = "bold") +
      scale_colour_identity() +
      facet_grid(facet_formula,
                 labeller = labeller(contrast = CONTRAST_LAB_1L,
                                     currency = CURRENCY_LAB)) +
      scale_fill_viridis_c(option = "viridis", limits = c(0, 100),
                           breaks = seq(0, 100, 25),
                           name = "Within-stratum\nfrequency (%)") +
      scale_x_continuous(breaks = seq(-65, -61, by = 1)) +
      scale_y_continuous(breaks = seq(46, 49, by = 1)) +
      coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
      labs(title = title, subtitle = subtitle) +
      theme_diag(base_size = 11) +
      theme(panel.grid   = element_line(colour = "grey92"),
            panel.spacing.x = unit(0.8, "lines"),
            axis.title   = element_blank(),
            axis.text    = element_text(size = 6),
            strip.text   = element_text(face = "bold"),
            strip.text.y = element_text(angle = 0),
            plot.title   = element_text(face = "bold"))

    save_fig(g, stem,
             width  = 2.55 * ncol_panels + 1.8,
             height = 2.45 * nrow_panels + 1.2,
             dpi = 200)
  }

  n_ct <- length(intersect(MAP_CONTRASTS, levels(droplevels(map_dat$contrast))))
  n_fm <- length(FAMILY_LEVELS)

  render_map("biomass", family ~ contrast,
             "Diagnostic families across survey strata - biomass",
             "Fig9_map_strata_biomass", n_ct, n_fm)

  render_map("occurrence", family ~ contrast,
             "Diagnostic families across survey strata - occurrence",
             "FigS4_map_strata_occurrence", n_ct, n_fm)

  render_map("biomass", family ~ contrast,
             "Diagnostic families across survey strata - biomass, low-n strata masked",
             "FigS5_map_strata_biomass_masked", n_ct, n_fm, mask_low_n = TRUE)

  render_map("occurrence", family ~ contrast,
             "Diagnostic families across survey strata - occurrence, low-n strata masked",
             "FigS6_map_strata_occurrence_masked", n_ct, n_fm, mask_low_n = TRUE)
}

# =============================================================================
# 8. THE NUMBERS BEHIND THE FIGURES
# =============================================================================
export_csv <- function(x, name) {
  x %>%
    mutate(contrast = contrast_label(contrast)) %>%
    write_csv(file.path(DIR_FIGURES, name))
}
export_csv(fam_gulf,  "data_Fig5_gulf_families.csv")
export_csv(fam_unit,  "data_Fig6_9_unit_families.csv")
export_csv(fam_level, "data_Fig10_crossscale.csv")

cat("\nFigures and their source tables written to ", normalizePath(DIR_FIGURES), "\n", sep = "")
