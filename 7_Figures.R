# =============================================================================
# 7_Figures.R - MANUSCRIPT FIGURES
# -----------------------------------------------------------------------------
# Builds the body figures from the saved pipeline results. Nothing is
# recomputed: every .rda already holds the per-cell diagnostics, so this script
# only aggregates, maps to the four families, and draws.
#
# Run after 6b, 6c and 6d have finished. Independent of which spatial level was
# last sourced, because it reads all three from disk.
#
# WHAT THIS SCRIPT PRODUCES
# -----------------------------------------------------------------------------
#   Fig5_regime_break_gulf        Manuscript Figure 5.
#                                 Family frequencies by contrast and currency at
#                                 the Gulf-wide level. The headline figure: the
#                                 within-decade contrasts set the baseline, PTa
#                                 and PTb sit outside it.
#
#   Fig6_regime_break_ecoregion   Manuscript Figure 6.
#                                 Same quantities per ecoregion, showing where
#                                 the inter-decade break is expressed.
#
#   Fig7_map_families_biomass     Manuscript Figure 7.
#   FigS6_map_families_occurrence Supplementary S6.
#                                 Within-stratum family frequency mapped on the
#                                 survey strata, one row per family, one column
#                                 per contrast. Every stratum is shown; values
#                                 built on fewer than N_FLAG predator x size
#                                 units carry an asterisk.
#
#   Fig7b_map_families_combined   Both currencies in one figure, for review.
#
#   Fig_crossscale_gradient       The §3.3 result: as the spatial unit widens
#                                 from stratum to ecoregion to the whole Gulf,
#                                 Stability falls and Reorganisation rises while
#                                 Substitution holds. This is the figure that
#                                 justifies reporting three levels.
#
#   Fig_heatmap_units             Dense alternative to Figure 6: family by
#                                 spatial unit, faceted by currency x contrast.
#                                 Useful when the stratum panel count is too
#                                 high for a bar chart.
#
# HOW THE NUMBERS ARE BUILT
#   Percentages are computed within each taxonomic resolution and then averaged
#   across the ~100 resolutions (freq_by_resolution() in the config module).
#   Inconclusive cells are excluded before the percentage, so the four families
#   sum to 100. Coverage — how many cells were classifiable — is reported
#   separately by 8_Tables.R rather than hidden in these figures.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(forcats)
  library(sf); library(readr); library(stringr)
})

PREY_FAMILY <- "_1"          # prey-grouping family used for the manuscript
source("R_helpers/Config_Mappings.R")   # builds RDA_DIRS from PREY_FAMILY

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
N_FLAG      <- 5                        # below this many units: flag the value
STRATA_PATH <- "strata_rv_gulf.rds"     # sf object, or leave to the package
COAST_PATH  <- "gulf.coast.intermediate.csv"   # optional coastline (x, y, pid)
MAP_XLIM    <- c(-66.2, -60.0)
MAP_YLIM    <- c(45.5, 49.2)
# Contrasts to show on the maps. 2006 vs 2018 is not estimable at the stratum scale, so
# including it would produce a column of empty panels.
MAP_CONTRASTS <- c("P1a", "P1b", "P2", "PTa")

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

# Cross-scale summary: the mean family composition at each level. For the two
# disaggregated levels this is the mean across units, i.e. a pool of results,
# not a pooled test. The distinction matters and is stated on the figure.
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
# One line per family, one linetype per currency, contrasts on the x axis in
# manuscript order. The shaded bands separate the three within-decade contrasts
# from the two between-decade ones: the visual claim is that PTa and PTb leave
# the envelope traced by P1a, P1b and P2.
cat("\nFigure 5: Gulf-wide regime break\n")

# Error bars = SD of the family frequencies across the ~100 resolutions
# (spread_by_resolution() in Config_Mappings.R; same keys as fam_gulf).
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
# 3. FIGURE 6 - BY ECOREGION
# =============================================================================
# The same reading, one panel per ecoregion. Stacked bars rather than lines:
# with four families summing to 100 the stack shows composition directly, and
# the contrast axis stays readable across four panels.
cat("\nFigure 6: by ecoregion\n")

# n_units varies by contrast, so it must not be part of the facet label
# (that would split each ecoregion into one facet per contrast); it is
# printed above each bar instead.
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
# 3b. REGIME BREAK AT THE ECOREGION AND STRATUM SCALES (one figure per scale)
# =============================================================================
# Same grammar as Figure 5 (lines per family, linetype per currency, shaded
# within/between-period bands, SD error bars), applied to the across-unit mean
# of each disaggregated scale. Values are those of Table 3; error bars are the
# SD across the ~100 resolutions of the across-unit mean
# (spread_by_resolution_units in Config_Mappings.R).
cat("\nRegime break at the ecoregion and stratum scales\n")

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
  "ecoregion", "Fig7_regime_break_ecoregion_mean",
  "Diet regime break at the ecoregion scale",
  paste0("Family frequencies averaged across the four ecoregions: mean and\n",
         "standard deviation across ~100 prey taxonomic resolutions"))

make_regime_break(
  "stratum", "Fig8_regime_break_stratum_mean",
  "Diet regime break at the stratum scale",
  paste0("Family frequencies averaged across survey strata: mean and\n",
         "standard deviation across ~100 prey taxonomic resolutions"))

# =============================================================================
# 4. HEATMAP - DENSE ALTERNATIVE FOR THE STRATUM SCALE
# =============================================================================
# With ~25 strata a bar panel per unit is unreadable. The heatmap keeps every
# unit on one page: family on x, unit on y, faceted by currency x contrast.
cat("\nHeatmap by spatial unit\n")

make_heatmap <- function(lvl, stem, height) {
  d <- long_of(filter(fam_unit, level == lvl)) %>%
    filter(!is.na(contrast), contrast %in% MAP_CONTRASTS) %>%
    mutate(unit_lab = paste0(spatial_unit, " (n=", n_units, ")"),
           unit_lab = fct_reorder(unit_lab, pct * (family == "Stability"),
                                  .fun = sum, .desc = FALSE))

  g <- ggplot(d, aes(family, unit_lab, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f", pct),
                  colour = ifelse(pct > 50, "white", "grey15")),
              fontface = "bold", size = 2.7) +
    scale_colour_identity() +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 100),
                         name = "Frequency\n(%)") +
    facet_grid(currency ~ contrast,
               labeller = labeller(currency = CURRENCY_LAB,
                                   contrast = CONTRAST_LAB_1L)) +
    labs(title = paste0("Family composition by ", LEVEL_SHORT[[lvl]]),
         subtitle = "Rows ordered by Stability; n = predator x size-class units",
         x = NULL, y = NULL) +
    theme_diag(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(size = 6.5),
          strip.text  = element_text(face = "bold", size = 8),
          panel.grid  = element_blank())
  save_fig(g, stem, width = 12, height = height)
}

if (any(fam_unit$level == "ecoregion")) make_heatmap("ecoregion", "Fig_heatmap_ecoregion", 5.5)
if (any(fam_unit$level == "stratum"))   make_heatmap("stratum",   "Fig_heatmap_stratum",   10)

# =============================================================================
# 5. CROSS-SCALE GRADIENT
# =============================================================================
# The methodological result of §3.3 and §4.3. Pooling the diet matrices across
# regionally distinct prey fields inflates apparent Reorganisation at the
# expense of Stability, while Substitution is scale-invariant. Restricted to
# the between-decade contrast, where the effect is largest and the claim is
# made.
cat("\nCross-scale gradient\n")

# Error bars: SD across resolutions of the Gulf-wide percentage (all_gulf)
# and of the across-unit mean (ecoregion, stratum).
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

save_fig(fig_grad, "Fig_crossscale_gradient", width = 9.5, height = 5.4)

# =============================================================================
# 6. FIGURE 7 - STRATUM MAPS
# =============================================================================
# Geometry comes from the gulf.spatial shapefile when available, otherwise from
# a pre-extracted sf object. gulf.spatial depends on the retired rgdal and
# often fails to load, hence the fallback.
cat("\nFigure 7: stratum maps\n")

# load_strata() now lives in R_helpers/Config_Mappings.R (shared with 9).
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

  render_map <- function(currencies, facet_formula, title, stem,
                         ncol_panels, nrow_panels) {
    d <- strata %>%
      left_join(filter(map_dat, currency %in% currencies), by = "spatial_unit") %>%
      filter(!is.na(pct)) %>%
      mutate(contrast = droplevels(factor(contrast, levels = CONTRAST_LEVELS)),
             currency = droplevels(factor(currency, levels = names(CURRENCY_LAB))))

    lab <- d %>%
      sf::st_drop_geometry() %>%
      left_join(lab_xy, by = "spatial_unit") %>%
      transmute(family, contrast, currency, X, Y,
                lab = ifelse(low_n, paste0(sprintf("%.0f", pct), "*"),
                             sprintf("%.0f", pct)),
                txt = ifelse(pct < 50, "grey95", "grey10"))

    g <- ggplot()
    if (!is.null(coast)) {
      g <- g + geom_polygon(data = coast, aes(x, y, group = pid),
                            fill = "grey88", colour = NA)
    }
    g <- g +
      geom_sf(data = strata, fill = "grey97", colour = "grey70", linewidth = 0.15) +
      geom_sf(data = d, aes(fill = pct), colour = "grey45", linewidth = 0.15) +
      geom_text(data = lab, aes(X, Y, label = lab, colour = txt),
                size = 2.1, fontface = "bold") +
      scale_colour_identity() +
      facet_grid(facet_formula,
                 labeller = labeller(contrast = CONTRAST_LAB_1L,
                                     currency = CURRENCY_LAB)) +
      scale_fill_viridis_c(option = "viridis", limits = c(0, 100),
                           breaks = seq(0, 100, 25),
                           name = "Within-stratum\nfrequency (%)") +
      coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +
      labs(title = title,
           subtitle = paste0("Every stratum shown; * = built on < ", N_FLAG,
                             " predator x size-class units (interpret with ",
                             "caution). 2006 vs 2018 is not estimable at this scale.")) +
      theme_diag(base_size = 11) +
      theme(panel.grid   = element_line(colour = "grey92"),
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
             "Diet-regime families across survey strata - biomass currency",
             "Fig7_map_families_biomass", n_ct, n_fm)

  render_map("occurrence", family ~ contrast,
             "Diet-regime families across survey strata - occurrence currency",
             "FigS6_map_families_occurrence", n_ct, n_fm)

  render_map(c("biomass", "occurrence"), family ~ currency + contrast,
             "Diet-regime families across survey strata - both currencies",
             "Fig7b_map_families_combined", 2 * n_ct, n_fm)
}

# =============================================================================
# 7. SAVE THE AGGREGATED TABLES BEHIND THE FIGURES
# =============================================================================
# 8_Tables.R rebuilds these independently; writing them here means a figure can
# always be traced back to the exact numbers it was drawn from.
# Contrast codes are replaced by explicit period labels in the exported CSVs.
export_csv <- function(x, name) {
  x %>%
    mutate(contrast = contrast_label(contrast)) %>%
    write_csv(file.path(DIR_FIGURES, name))
}
export_csv(fam_gulf,  "data_Fig5_gulf_families.csv")
export_csv(fam_unit,  "data_Fig6_7_unit_families.csv")
export_csv(fam_level, "data_crossscale_gradient.csv")

cat("\nFigures and their source tables written to ", normalizePath(DIR_FIGURES), "\n", sep = "")
