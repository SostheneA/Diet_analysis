# =============================================================================
# CROSS-LEVEL COMPARISON
# -----------------------------------------------------------------------------
# Run after 6b, 6c and 6d have all finished. Reads the three saved runs and
# builds the tables that justify reporting more than one spatial level.
#
# This script sources no engine: the three levels write different spatial column
# names ("Area" for L1 and L2, "str" for L3), and the helpers below reconcile
# them without needing SPATIAL_OUT to be set. Sourcing an engine here would fix
# the session to one level and is not wanted.
#
# The three levels run on the same stomachs — the source data carries no missing
# Area and no missing stratum — so every difference below is attributable to the
# spatial treatment, not to a change of sample. Section 1 verifies that rather
# than assuming it.
#
# THE THREE LEVELS ARE NOT A NESTED HIERARCHY
#   Ecoregions are assigned per stomach by a nearest-feature spatial join on the
#   capture position; strata are the survey's own polygons. The two geometries
#   do not coincide, so a stratum that crosses an ecoregion boundary is split
#   between ecoregions. all_gulf contains both of the finer levels, but stratum
#   is not nested inside ecoregion: they are two alternative partitions of the
#   same domain at different resolutions.
#
#   Section 0 rebuilds the stratum-to-ecoregion crosswalk from the data and
#   reports which strata are split. Section 5 then compares the two partitions
#   only on the strata that sit in a single ecoregion, where the mapping is
#   unambiguous, and states how much of the data that covers.
#
#   The cross-level gradient in sections 2 to 4 does not require nesting — only
#   that resolution decreases from stratum to ecoregion to Gulf — so it is
#   unaffected.
# =============================================================================

rm(list = ls())

PREY_FAMILY <- "_1"

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

ALPHA <- 0.05

# Reserved names so that the three environments do not collide.
env_L1 <- new.env(); load(paste0("data/Sensitivity/all_runs_all_gulf",  PREY_FAMILY, ".rda"), envir = env_L1)
env_L2 <- new.env(); load(paste0("data/Sensitivity/all_runs_ecoregion", PREY_FAMILY, ".rda"), envir = env_L2)
env_L3 <- new.env(); load(paste0("data/Sensitivity/all_runs_stratum",   PREY_FAMILY, ".rda"), envir = env_L3)

# Maps the level-specific spatial column onto a common `spatial_unit` and keeps
# the `spatial_level` tag the engine wrote into every row.
tidy_level <- function(res) {
  d <- res$results
  sp_col <- intersect(c("Area", "str"), names(d))[1]
  d$spatial_unit <- if (is.na(sp_col)) NA_character_ else as.character(d[[sp_col]])
  dplyr::select(d, -dplyr::any_of(c("Area", "str")))
}

pick <- function(env, obj) tidy_level(get(obj, envir = env))

# The PT scenario in biomass is the reference comparison; change these two
# lines to compare any other scenario or mode.
OBJ  <- "res_biomass_PT"
all3 <- bind_rows(pick(env_L1, OBJ), pick(env_L2, OBJ), pick(env_L3, OBJ))

# =============================================================================
# 0. STRATUM-TO-ECOREGION CROSSWALK
# =============================================================================
# Built from the source data rather than assumed. Counts are prey records, not
# stomachs, so the percentages are indicative of where the two partitions
# disagree, not of a stomach tally.
load("data/dat_classed.rda")
# stratum is numeric in the source data but character in the results, where it
# was written through factor(); casting here keeps the join types compatible.
xwalk_raw <- as.data.frame(dat_classed) %>%
  mutate(stratum = as.character(stratum)) %>%
  count(Area, stratum, name = "n_rec") %>%
  filter(n_rec > 0)

xwalk <- xwalk_raw %>%
  group_by(stratum) %>%
  mutate(
    n_stratum   = sum(n_rec),
    share       = n_rec / n_stratum,
    n_ecoregions = dplyr::n()
  ) %>%
  slice_max(n_rec, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(stratum,
            eco_majority = Area,
            share_majority = round(share, 3),
            n_ecoregions,
            pure = n_ecoregions == 1)

cat("\n--- Strata split across ecoregions ---\n")
print(
  xwalk_raw %>%
    group_by(stratum) %>%
    filter(dplyr::n() > 1) %>%
    mutate(share = round(100 * n_rec / sum(n_rec), 1)) %>%
    ungroup() %>%
    arrange(stratum, desc(n_rec)),
  n = Inf
)

cat("\n--- Crosswalk coverage ---\n")
print(
  xwalk %>%
    summarise(
      n_strata        = dplyr::n(),
      n_pure          = sum(pure),
      n_split         = sum(!pure),
      pct_records_pure = round(
        100 * sum(xwalk_raw$n_rec[xwalk_raw$stratum %in% stratum[pure]]) /
          sum(xwalk_raw$n_rec), 1)
    )
)

# =============================================================================
# 1. SAMPLE CHECK
# =============================================================================
# Stomach totals must agree across levels. They are summed over cells, so they
# are not identical by construction — a cell exists only where both periods pass
# N_MIN, and that filter bites differently at each grain. What must hold is that
# no level is missing whole regions of the data.
cat("\n--- Cells, units and stomach totals per level ---\n")
print(
  all3 %>%
    group_by(spatial_level) %>%
    summarise(
      n_cells      = dplyr::n(),
      n_units      = dplyr::n_distinct(spatial_unit),
      n_resolutions = dplyr::n_distinct(x_threshold),
      stomachs     = sum(n_sto_P1 + n_sto_P2) / dplyr::n_distinct(x_threshold),
      pct_testable = round(100 * mean(testable), 1),
      pct_reliable = round(100 * mean(reliable), 1),
      .groups = "drop"
    )
)

# =============================================================================
# 2. DIAGNOSTIC COMPOSITION PER LEVEL
# =============================================================================
# The headline table: what changing the spatial grain does to the
# classification. Percentages are within level, over testable cells only.
cat("\n--- Diagnostic composition per level (testable cells) ---\n")
comp_table <- all3 %>%
  filter(testable) %>%
  count(spatial_level, diagnostic) %>%
  group_by(spatial_level) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

print(
  comp_table %>%
    select(-n) %>%
    pivot_wider(names_from = spatial_level, values_from = pct, values_fill = 0),
  n = Inf
)

# =============================================================================
# 3. EFFECT SIZE VS POWER
# =============================================================================
# Separates the two reasons a diagnostic class can empty out at a finer grain.
# BC and R2 are effect sizes and do not depend on significance; the pct_*_sig
# columns do. Effect sizes holding steady while detection rates fall is a power
# problem, driven by k_med — the number of tows per cell.
cat("\n--- Effect size vs detection rate, by level ---\n")
print(
  all3 %>%
    filter(testable) %>%
    group_by(spatial_level) %>%
    summarise(
      n            = dplyr::n(),
      k_med        = median(pmin(n_set_P1, n_set_P2), na.rm = TRUE),
      BC_med       = round(median(BC, na.rm = TRUE), 3),
      R2_med       = round(median(R2_comp, na.rm = TRUE), 3),
      pct_H_sig    = round(100 * mean(p_H    < ALPHA, na.rm = TRUE), 1),
      pct_Bs_sig   = round(100 * mean(p_Bs   < ALPHA, na.rm = TRUE), 1),
      pct_comp_sig = round(100 * mean(p_comp < ALPHA, na.rm = TRUE), 1),
      .groups = "drop"
    )
)

# =============================================================================
# 4. THE DISPERSION ARGUMENT FOR REPORTING MORE THAN ONE LEVEL
# =============================================================================
# At all_gulf, spatial heterogeneity is not modelled and lands in within-group
# variance. comp_dispersion flags the cells where the composition and dispersion
# tests are both significant, i.e. where location and spread effects cannot be
# told apart. A rate that falls as the grain gets finer is the quantitative
# argument that the disaggregated levels are doing work.
cat("\n--- Confounded location/dispersion cells, by level ---\n")
print(
  all3 %>%
    filter(testable) %>%
    group_by(spatial_level) %>%
    summarise(
      pct_comp_sig   = round(100 * mean(p_comp < ALPHA, na.rm = TRUE), 1),
      pct_disp_sig   = round(100 * mean(p_disp < ALPHA, na.rm = TRUE), 1),
      pct_confounded = round(100 * mean(comp_dispersion, na.rm = TRUE), 1),
      .groups = "drop"
    )
)

# =============================================================================
# 5. AGREEMENT BETWEEN THE TWO DISAGGREGATED PARTITIONS
# =============================================================================
# L2 and L3 share species, size class and resolution but not the spatial unit.
# Because the two partitions are not nested, the comparison is restricted to the
# strata that fall entirely within one ecoregion; a split stratum has no single
# ecoregional counterpart and would make the join arbitrary.
#
# Within that restriction one ecoregional cell still matches several stratum
# cells, so the table reads: given an ecoregional verdict, how were the strata
# inside that ecoregion classified.
key <- c("x_threshold", "species", "size_class")

L3_pure <- pick(env_L3, OBJ) %>%
  filter(testable) %>%
  left_join(xwalk, by = c("spatial_unit" = "stratum")) %>%
  filter(pure)

cat("\n--- Coverage of the L2 / L3 comparison ---\n")
print(
  pick(env_L3, OBJ) %>%
    filter(testable) %>%
    left_join(xwalk, by = c("spatial_unit" = "stratum")) %>%
    summarise(
      n_stratum_cells = dplyr::n(),
      n_kept          = sum(pure, na.rm = TRUE),
      pct_kept        = round(100 * mean(pure, na.rm = TRUE), 1)
    )
)

cat("\n--- Ecoregional verdict vs stratum verdicts (unambiguous strata only) ---\n")
print(
  inner_join(
    pick(env_L2, OBJ) %>% filter(testable) %>%
      select(all_of(key), eco_unit = spatial_unit, d_eco = diagnostic),
    L3_pure %>%
      select(all_of(key), str_unit = spatial_unit,
             eco_unit = eco_majority, d_str = diagnostic),
    by = c(key, "eco_unit"), relationship = "many-to-many"
  ) %>%
    count(d_eco, d_str) %>%
    group_by(d_eco) %>%
    mutate(pct_within_eco = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    arrange(d_eco, desc(n)),
  n = Inf
)

# Headline number: how often the ecoregional verdict matches the modal verdict
# of the strata it contains. A low value means the ecoregional reading averages
# over genuinely heterogeneous strata.
cat("\n--- Ecoregion vs modal stratum verdict ---\n")
print(
  inner_join(
    pick(env_L2, OBJ) %>% filter(testable) %>%
      select(all_of(key), eco_unit = spatial_unit, d_eco = diagnostic),
    L3_pure %>%
      select(all_of(key), eco_unit = eco_majority, d_str = diagnostic),
    by = c(key, "eco_unit"), relationship = "many-to-many"
  ) %>%
    count(across(all_of(c(key, "eco_unit"))), d_eco, d_str) %>%
    group_by(across(all_of(c(key, "eco_unit")))) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    summarise(
      n_eco_cells = dplyr::n(),
      pct_match   = round(100 * mean(d_eco == d_str), 1)
    )
)

# =============================================================================
# 6. FIGURE
# =============================================================================
p_levels <- comp_table %>%
  mutate(spatial_level = factor(spatial_level,
                                levels = c("all_gulf", "ecoregion", "stratum"))) %>%
  ggplot(aes(x = reorder(diagnostic, pct), y = pct, fill = spatial_level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", name = "Spatial level") +
  labs(
    title = "Diagnostic composition across spatial levels",
    subtitle = "Biomass, 2004-2006 vs 2018-2019, testable cells only",
    x = NULL, y = "Within-level frequency (%)",
    caption = paste0(
      "all_gulf pools the diet matrices before testing; ecoregion and stratum ",
      "test within units.\nSame stomachs at every level. Ecoregion and stratum ",
      "are alternative partitions, not a nested hierarchy."
    )
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"),
        legend.position = "bottom")

dir.create(paste0("Sensitivity_Plot_comparison", PREY_FAMILY), recursive = TRUE, showWarnings = FALSE)
ggsave("Sensitivity_Plot_comparison/diagnostic_composition_by_level.png",
       p_levels, width = 10, height = 6, dpi = 300)
print(p_levels)

cat("\nCross-level comparison complete.\n")
