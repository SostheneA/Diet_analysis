# =============================================================================
# 0_Config_Mappings.R - SHARED CONFIGURATION AND READERS
# -----------------------------------------------------------------------------
# Single source of truth for everything scripts 7 to 10 have in common: the
# diagnostic typology, the contrast naming, the spatial levels, the colour
# palettes, and the readers that turn the pipeline's .rda files into tidy
# tables. Change a label here and every table, figure and map follows.
#
# This file defines objects only. It loads no data and writes no output, so it
# is safe to source at the top of any script and safe to source twice.
#
# It deliberately does NOT source the engine (6a_Engine_Trophic.R). The engine
# fixes SPATIAL_SOURCE and SPATIAL_OUT to one level at source() time; scripts
# 7 to 10 read all three levels at once and must not be tied to any of them.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(stringr)
})

# =============================================================================
# 1. THE TYPOLOGY
# =============================================================================
# The engine writes one of nine diagnostic states per cell, plus "Inconclusive"
# where a signal could not be computed. The manuscript reports four families.
# The mapping answers two questions: did the prey set turn over (composition),
# and did the functional structure change (H' and/or Bs)?
#
#   Stability          C = 0, no functional change   (+ Emerging Shift, folded in)
#   Substitution       C = 1, no functional change
#   Functional change  C = 0, H' and/or Bs significant
#   Reorganisation     C = 1, H' and/or Bs significant
#
# Inconclusive sits outside the four families: it is a coverage flag, not an
# ecological outcome, and is excluded before any percentage is computed.

DIAG_LEVELS <- c(
  "Stable Diet", "Emerging Shift",                                  # Stability
  "Ghost Shift",                                                    # Substitution
  "Niche Compression/Expansion", "Internal Rebalancing",
  "Niche Restructuring",                                            # Functional change
  "Partial Diet Shift", "Structural Shift", "Major Shift"           # Reorganisation
)

INCONCLUSIVE_LAB <- "Inconclusive"

FAMILY_LEVELS <- c("Stability", "Substitution",
                   "Functional change", "Reorganisation")

STATE_TO_FAMILY <- c(
  "Stable Diet"                 = "Stability",
  "Emerging Shift"              = "Stability",
  "Ghost Shift"                 = "Substitution",
  "Niche Compression/Expansion" = "Functional change",
  "Internal Rebalancing"        = "Functional change",
  "Niche Restructuring"         = "Functional change",
  "Partial Diet Shift"          = "Reorganisation",
  "Structural Shift"            = "Reorganisation",
  "Major Shift"                 = "Reorganisation"
)

# The signal combination behind each state, for the typology table.
STATE_LOGIC <- tibble::tribble(
  ~diagnostic,                    ~composition,        ~H_prime,          ~Bs,
  "Stable Diet",                  "ns",                "ns",              "ns",
  "Emerging Shift",               "trend (0.05-0.10)", "ns",              "ns",
  "Ghost Shift",                  "significant",       "ns",              "ns",
  "Niche Compression/Expansion",  "ns",                "ns",              "significant",
  "Internal Rebalancing",         "ns",                "significant",     "ns",
  "Niche Restructuring",          "ns",                "significant",     "significant",
  "Partial Diet Shift",           "significant",       "significant",     "ns",
  "Structural Shift",             "significant",       "ns",              "significant",
  "Major Shift",                  "significant",       "significant",     "significant"
)

family_of <- function(diag) {
  factor(unname(STATE_TO_FAMILY[as.character(diag)]), levels = FAMILY_LEVELS)
}

# Collapses a long table carrying one percentage per diagnostic into one row
# per key with the four family columns. The nine states partition 100 %, so the
# four families do too; no residual category is needed.
add_families <- function(df, keys = character(0),
                         diag_col = "diagnostic", val_col = "pct") {
  long <- df %>%
    mutate(.family = family_of(.data[[diag_col]])) %>%
    filter(!is.na(.family)) %>%
    group_by(across(all_of(c(keys, ".family")))) %>%
    summarise(.v = sum(.data[[val_col]], na.rm = TRUE), .groups = "drop")

  wide <- pivot_wider(long, names_from = .family, values_from = .v, values_fill = 0)
  for (m in setdiff(FAMILY_LEVELS, names(wide))) wide[[m]] <- 0
  select(wide, all_of(c(keys, FAMILY_LEVELS)))
}

# =============================================================================
# 2. CONTRASTS
# =============================================================================
# The run scripts label scenarios P1 to PT; the manuscript uses P1a to PTa.
# Rather than recode from the filename, the contrast is derived from the period
# labels stored inside the results, which cannot drift out of sync with the
# data. CONTRAST_RECODE is kept as a fallback for older .rda files.

CONTRAST_LEVELS <- c("P1a", "P1b", "P2", "PTb", "PTa")

CONTRAST_FROM_PERIODS <- c(
  "2004|2006"                = "P1a",
  "2004-2005|2006"           = "P1b",
  "2018|2019"                = "P2",
  "2006|2018"                = "PTb",
  "2004-2006|2018-2019"      = "PTa"
)

CONTRAST_RECODE <- c(P1 = "P1a", P3 = "P1b", P2 = "P2", P4 = "PTb", PT = "PTa")

# Manuscript labels: explicit periods only. The short codes (P1a ... PTa) stay
# internal (filters, ordering) and never appear in figures, tables or CSVs.
CONTRAST_LAB_1L <- c(
  P1a = "2004 vs 2006",
  P1b = "2004-2005 vs 2006",
  P2  = "2018 vs 2019",
  PTb = "2006 vs 2018",
  PTa = "2004-2006 vs 2018-2019"
)

CONTRAST_LAB <- c(
  P1a = "2004\nvs 2006",
  P1b = "2004-2005\nvs 2006",
  P2  = "2018\nvs 2019",
  PTb = "2006\nvs 2018",
  PTa = "2004-2006\nvs 2018-2019"
)

# One-line label for a contrast code (or vector of codes); NA-safe. Returned as
# a factor in manuscript order so that arrange()/facets keep the P1a..PTa order
# (alphabetical order of the period strings would not).
contrast_label <- function(x) {
  factor(unname(CONTRAST_LAB_1L[as.character(x)]), levels = unname(CONTRAST_LAB_1L))
}

INTRA_CONTRASTS <- c("P1a", "P1b", "P2")
INTER_CONTRASTS <- c("PTb", "PTa")

contrast_from_periods <- function(p1, p2) {
  key <- paste(p1, p2, sep = "|")
  factor(unname(CONTRAST_FROM_PERIODS[key]), levels = CONTRAST_LEVELS)
}

relabel_contrast <- function(x) {
  factor(unname(CONTRAST_RECODE[as.character(x)]), levels = CONTRAST_LEVELS)
}

# =============================================================================
# 3. SPATIAL LEVELS
# =============================================================================
# The three levels produced by 6b, 6c and 6d. Ecoregion and stratum are two
# alternative partitions of the same domain, not a nested hierarchy: ecoregions
# are assigned per stomach by a spatial join, strata are the survey's polygons,
# and eleven strata straddle an ecoregion boundary. all_gulf contains both.

SPATIAL_LEVELS_ORD <- c("all_gulf", "ecoregion", "stratum")

LEVEL_LAB <- c(
  all_gulf  = "Gulf-wide (no spatial split)",
  ecoregion = "By ecoregion",
  stratum   = "By survey stratum"
)

LEVEL_SHORT <- c(all_gulf = "Gulf", ecoregion = "Ecoregion", stratum = "Stratum")

# Output column each level writes its spatial unit into.
LEVEL_SPATIAL_COL <- c(all_gulf = "Area", ecoregion = "Area", stratum = "str")

# Prey-grouping family used for the manuscript. The run scripts (6b-6d) write
# to Sensitivity_<level><PREY_FAMILY>/ ; "_1" is the main analysis, "_2" the
# original pooled grouping (kept for the appendix comparison, see 6f).
# Define PREY_FAMILY before sourcing this file to override.
if (!exists("PREY_FAMILY")) PREY_FAMILY <- "_1"

# Folders read by 7/8/9, derived from PREY_FAMILY. Override before sourcing if
# needed (e.g. RDA_DIRS for a different family).
if (!exists("RDA_DIRS")) {
  RDA_DIRS <- c(all_gulf  = paste0("Sensitivity_all_gulf",  PREY_FAMILY),
                ecoregion = paste0("Sensitivity_ecoregion", PREY_FAMILY),
                stratum   = paste0("Sensitivity_stratum",   PREY_FAMILY))
}

# =============================================================================
# 4. AESTHETICS
# =============================================================================

FAMILY_PAL <- c(
  "Stability"         = "#2F4A5A",   # dark slate
  "Substitution"      = "#D96C4A",   # terracotta
  "Functional change" = "#7F9CAB",   # muted steel
  "Reorganisation"    = "#9B2226"    # deep red
)

FAMILY_SHAPE <- c("Stability" = 16, "Substitution" = 15,
                  "Functional change" = 18, "Reorganisation" = 17)

CURRENCY_LAB <- c(biomass = "Biomass", occurrence = "Occurrence")
CURRENCY_PAL <- c(biomass = "#2c7bb6", occurrence = "#d7191c")

LEVEL_PAL <- c(all_gulf = "#4C6A78", ecoregion = "#C08552", stratum = "#5F7A4F")

# Shaded bands separating intra- from inter-decade contrasts on a 5-point axis.
decade_bands <- function() {
  list(
    ggplot2::annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
                      fill = "#D7E6EC", alpha = 0.55),
    ggplot2::annotate("rect", xmin = 3.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
                      fill = "#F6E0D6", alpha = 0.65),
    ggplot2::geom_vline(xintercept = 3.5, linetype = "dashed",
                        colour = "grey60", linewidth = 0.3)
  )
}

# Period colours, shared by every script that colours by period (10, 6f).
PERIOD_PAL <- c("2004-2006" = "#2c7bb6", "2018-2019" = "#d7191c")

# Ecoregion display names (co-author request: English names on all outputs).
AREA_LAB <- c("Baie des Chaleurs" = "Chaleur Bay",
              "Chaleur Bay"       = "Chaleur Bay")
area_label <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(AREA_LAB), unname(AREA_LAB[x]), x)
}

# One theme for every figure (7, 9, 10, 6f): same font family and sizes.
# Set FIG_FONT before sourcing to change the family; "sans" resolves to Arial /
# Helvetica on Windows and macOS PDF/PNG devices.
if (!exists("FIG_FONT")) FIG_FONT <- "sans"

theme_diag <- function(base_size = 11, base_family = FIG_FONT) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top",
      legend.title     = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(size = 8),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.caption     = ggplot2::element_text(hjust = 0, size = base_size - 2)
    )
}
# Make it the default for every ggplot built after sourcing this file.
if (requireNamespace("ggplot2", quietly = TRUE)) ggplot2::theme_set(theme_diag())

# =============================================================================
# 5. READING THE PIPELINE OUTPUT
# =============================================================================
# Filenames are {mode}_{scenario}_{level}_{period1}_vs_{period2}.rda. The level
# "all_gulf" contains an underscore, so positional splitting on "_" is wrong;
# the level is matched against the known set instead. Where the results carry
# the columns written by the engine (spatial_level, period_1, period_2, mode),
# those take precedence over anything parsed from the filename.

parse_run_name <- function(path) {
  b  <- sub("\\.rda$", "", basename(path))
  # Scenario token may carry the family suffix ("P1" or "P1_1" / "PT_PP").
  rx <- paste0("^(biomass|occurrence)_([A-Za-z0-9_]+)_(",
               paste(SPATIAL_LEVELS_ORD, collapse = "|"), ")_(.+)_vs_(.+)$")
  m  <- regmatches(b, regexec(rx, b))[[1]]
  if (!length(m)) {
    return(list(currency = NA_character_, scenario = NA_character_,
                family = NA_character_, level = NA_character_,
                period_1 = NA_character_, period_2 = NA_character_))
  }
  sc  <- m[3]
  fam <- if (grepl("_", sc)) sub("^[^_]+", "", sc) else ""
  sc  <- sub("_.*$", "", sc)
  list(currency = m[2], scenario = sc, family = fam, level = m[4],
       period_1 = gsub("_", "-", m[5]), period_2 = gsub("_", "-", m[6]))
}

# Returns the results tibble of one run, with currency / contrast / level
# attached and the spatial unit renamed to a common `spatial_unit` column.
read_run <- function(path) {
  e <- new.env(); load(path, envir = e)
  if (!exists("res", envir = e)) stop("No object 'res' in ", path)
  r <- tibble::as_tibble(e$res$results)
  if (!nrow(r)) return(NULL)

  meta <- parse_run_name(path)

  lvl <- if ("spatial_level" %in% names(r)) as.character(r$spatial_level[1]) else meta$level
  cur <- if ("mode" %in% names(r))          as.character(r$mode[1])          else meta$currency
  p1  <- if ("period_1" %in% names(r))      as.character(r$period_1[1])      else meta$period_1
  p2  <- if ("period_2" %in% names(r))      as.character(r$period_2[1])      else meta$period_2

  ct <- contrast_from_periods(p1, p2)
  if (is.na(ct)) ct <- relabel_contrast(meta$scenario)   # fallback for old files

  sp_col <- intersect(c("Area", "str"), names(r))[1]
  r$spatial_unit <- if (is.na(sp_col)) NA_character_ else as.character(r[[sp_col]])

  r %>%
    select(-any_of(c("Area", "str"))) %>%
    mutate(currency = cur, contrast = ct, level = lvl,
           period_1 = p1, period_2 = p2, .before = 1)
}

# Reads every run of one or more levels into a single long table.
read_all_runs <- function(levels = SPATIAL_LEVELS_ORD, dirs = RDA_DIRS) {
  out <- list()
  for (lv in levels) {
    d <- dirs[[lv]]
    if (is.null(d) || !dir.exists(d)) {
      warning("Folder not found for level '", lv, "': ", d, call. = FALSE)
      next
    }
    fs <- list.files(d, pattern = "\\.rda$", full.names = TRUE)
    if (!length(fs)) warning("No .rda in ", d, call. = FALSE)
    for (f in fs) out[[length(out) + 1]] <- read_run(f)
  }
  res <- bind_rows(out)
  if (!nrow(res)) stop("No results read. Check RDA_DIRS and that 6b-6d have run.")
  res %>% mutate(level = factor(level, levels = SPATIAL_LEVELS_ORD))
}

# =============================================================================
# 6. THE HEADLINE AGGREGATION RULE
# =============================================================================
# The manuscript rule, applied identically everywhere: within each taxonomic
# resolution (x_threshold) compute the percentage of cells in each diagnostic,
# then average those percentages across the ~100 resolutions. Averaging the
# percentages rather than pooling the cells keeps every resolution equally
# weighted, which is the point of the sweep.
#
# `keys` names the grouping beyond the resolution — typically currency,
# contrast, level and, for the disaggregated levels, spatial_unit.
# Inconclusive cells are dropped first: they are a coverage flag, and leaving
# them in would make the four families sum to less than 100.

freq_by_resolution <- function(results, keys = c("currency", "contrast", "level")) {
  base <- results %>% filter(diagnostic %in% DIAG_LEVELS)
  if (!nrow(base)) return(base[0, ])

  cnt <- base %>% count(across(all_of(c(keys, "x_threshold"))), diagnostic, name = "n")
  tot <- cnt %>%
    group_by(across(all_of(c(keys, "x_threshold")))) %>%
    summarise(n_tot = sum(n), .groups = "drop")

  # Explicit full grid rather than complete(): a diagnostic absent from a
  # resolution is a zero, not a missing value, and dropping it would inflate
  # the mean of the diagnostics that are present.
  grid <- tot %>%
    select(all_of(c(keys, "x_threshold"))) %>%
    distinct() %>%
    tidyr::crossing(diagnostic = DIAG_LEVELS)

  grid %>%
    left_join(cnt, by = c(keys, "x_threshold", "diagnostic")) %>%
    left_join(tot, by = c(keys, "x_threshold")) %>%
    mutate(n = dplyr::coalesce(n, 0L), pct = 100 * n / n_tot) %>%
    group_by(across(all_of(c(keys, "diagnostic")))) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    mutate(diagnostic = factor(diagnostic, levels = DIAG_LEVELS))
}

# Spread of the family frequencies across the resolution sweep, for error
# bars. Same construction as freq_by_resolution(): the four family percentages
# are computed WITHIN each resolution on the zero-filled grid, then summarised
# across resolutions. Returns one row per key x family with pct_mean, pct_sd,
# pct_min, pct_max and n_res. `keys` must be the same keys used for the
# corresponding mean table so that the two can be joined.
#   fam_gulf_sd <- spread_by_resolution(filter(res_all, level == "all_gulf"),
#                                       keys = c("currency", "contrast"))
family_by_resolution <- function(results, keys = c("currency", "contrast")) {
  results %>%
    filter(diagnostic %in% DIAG_LEVELS) %>%
    mutate(family = family_of(diagnostic)) %>%
    count(across(all_of(c(keys, "x_threshold"))), family, name = "n") %>%
    group_by(across(all_of(c(keys, "x_threshold")))) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    tidyr::complete(tidyr::nesting(!!!rlang::syms(c(keys, "x_threshold"))),
                    family = factor(FAMILY_LEVELS, levels = FAMILY_LEVELS),
                    fill = list(pct = 0))
}

spread_by_resolution <- function(results, keys = c("currency", "contrast")) {
  family_by_resolution(results, keys) %>%
    group_by(across(all_of(c(keys, "family")))) %>%
    summarise(pct_mean = mean(pct), pct_sd = sd(pct),
              pct_min = min(pct), pct_max = max(pct),
              n_res = dplyr::n_distinct(x_threshold), .groups = "drop")
}

# For the disaggregated levels the reported value is a mean across spatial
# units; its spread across resolutions is the SD of that across-unit mean,
# i.e. compute the per-unit percentages within each resolution, average them
# over units, then take the SD over resolutions.
spread_by_resolution_units <- function(results, keys = c("currency", "contrast", "level"),
                                       unit_col = "spatial_unit") {
  family_by_resolution(results, keys = c(keys, unit_col)) %>%
    group_by(across(all_of(c(keys, "x_threshold", "family")))) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(across(all_of(c(keys, "family")))) %>%
    summarise(pct_mean = mean(pct), pct_sd = sd(pct),
              pct_min = min(pct), pct_max = max(pct),
              n_res = dplyr::n_distinct(x_threshold), .groups = "drop")
}

# Number of independent predator x size-class units behind each group. This is
# the n used for the low-sample flag and for the cluster-corrected error bars;
# it is not the number of rows, which counts resolutions as well.
n_units_by <- function(results, keys = c("currency", "contrast", "level")) {
  results %>%
    distinct(across(all_of(c(keys, "species", "size_class")))) %>%
    count(across(all_of(keys)), name = "n_units")
}

# Coverage: the share of cells that could not be classified. Reported alongside
# every family table, because a family composition computed on 40 % of the
# cells means something different from one computed on 95 %.
coverage_by <- function(results, keys = c("currency", "contrast", "level")) {
  results %>%
    group_by(across(all_of(keys))) %>%
    summarise(
      n_cells        = dplyr::n(),
      pct_testable   = round(100 * mean(testable), 1),
      pct_reliable   = round(100 * mean(reliable), 1),
      pct_confounded = round(100 * mean(comp_dispersion, na.rm = TRUE), 1),
      .groups = "drop"
    )
}

# =============================================================================
# 7. SHARED OUTPUT FOLDERS
# =============================================================================
if (!exists("DIR_TABLES"))  DIR_TABLES  <- "Output_Tables"
if (!exists("DIR_FIGURES")) DIR_FIGURES <- "Output_Figures"
if (!exists("DIR_APPEND"))  DIR_APPEND  <- "Output_Appendices"
if (!exists("DIR_DRIVERS")) DIR_DRIVERS <- "Output_Drivers"

for (d in c(DIR_TABLES, DIR_FIGURES, DIR_APPEND, DIR_DRIVERS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Survey-stratum polygons: the gulf.spatial shapefile when the package is
# installed, otherwise a local sf object saved as .rds; NULL if neither exists.
load_strata <- function(rds_path = "strata_rv_gulf.rds") {
  shp <- suppressWarnings(
    system.file("extdata/shapefiles/survey.stratum.polygons.shp",
                package = "gulf.spatial"))
  if (nzchar(shp) && file.exists(shp)) {
    s <- sf::read_sf(shp)
    s <- s[s$survey == "rv" & s$region == "gulf" & s$type == "polygon", ]
    s$str <- as.character(s$stratum)
    return(s)
  }
  if (file.exists(rds_path)) return(readRDS(rds_path))
  NULL
}

# Writes a figure in both raster and vector form under the same stem.
save_fig <- function(plot, stem, width, height, dpi = 300, dir = DIR_FIGURES) {
  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), plot,
                  width = width, height = height, dpi = dpi,
                  bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(dir, paste0(stem, ".pdf")), plot,
                  width = width, height = height, limitsize = FALSE)
  cat("  wrote ", stem, ".png / .pdf\n", sep = "")
  invisible(plot)
}

message("0_Config_Mappings.R loaded: ",
        length(DIAG_LEVELS), " states -> ", length(FAMILY_LEVELS), " families, ",
        length(SPATIAL_LEVELS_ORD), " spatial levels.")
