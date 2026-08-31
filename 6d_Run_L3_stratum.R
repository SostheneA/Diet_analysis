# =============================================================================
# RUN SCRIPT L3 - stratum
# -----------------------------------------------------------------------------
# One of three run scripts driving 6a_Engine_Trophic.R. They differ only in
# SPATIAL_LEVEL and in their output folders; the analytical parameters below are
# identical across all three on purpose, so that any difference between levels
# is attributable to the spatial treatment alone.
#
#   6b_Run_L1_all_gulf.R    SPATIAL_LEVEL <- "all_gulf"
#   6c_Run_L2_ecoregion.R   SPATIAL_LEVEL <- "ecoregion"
#   6d_Run_L3_stratum.R     SPATIAL_LEVEL <- "stratum"
#
# The three are independent and write to separate folders, so they can be
# launched in three parallel R sessions. Each one must run in its OWN session:
# the engine derives SPATIAL_SOURCE, SPATIAL_OUT and SPATIAL_POOL from
# SPATIAL_LEVEL at source() time, so sourcing it twice with different levels in
# one session would leave the second level in force for everything.
#
# WHAT THIS LEVEL ANSWERS
#   The same question as L2 at the finest grain available. A stratum is far
#   smaller than an ecoregion, so cells rest on 3-5 tows instead of 20-25.
#   Two consequences to check before interpreting anything:
#     * more cells fall below MIN_SETS and come back "Inconclusive";
#     * the three tests lose power, which mechanically empties the diagnostic
#       classes that require two significant signals.
#   audit_set_depth() in section 6b is the first thing to read, and the effect
#   size table that follows it is what separates "smaller effects" from
#   "less power".
# =============================================================================

rm(list = setdiff(ls(), "PREY_FAMILY"))   # PREY_FAMILY may be set by a driver
if (!exists("PREY_FAMILY")) PREY_FAMILY <- "_1"   # "_1" (main), "_2" or "_PP"
# =============================================================================
# 1. GLOBAL ANALYSIS PARAMETERS & SWITCHES
# =============================================================================
# Set before sourcing the engine. The engine declares each of these with
# if (!exists(...)), so the values here are the ones that apply.

# Sample size and permutation constraints
N_MIN    <- 5       # Min stomachs per cell for individual metrics
N_STRICT <- 25      # Min total stomachs per period for reliable classification
MIN_SETS <- 3       # Min sets per period to compute jackknife / Welch
R_PERM   <- 999     # Permutations for PERMANOVA, dispersion and IndVal
ALPHA    <- 0.05    # Significance threshold

# Methodological switches
BS_TEST         <- "welch"    # Welch's t-test for niche breadth (Bs)
COMP_RELATIVE   <- TRUE       # Standardize by set totals before Bray-Curtis
OCC_BINARY      <- FALSE      # Bray-Curtis on proportional occurrence, not Jaccard
ADD_NICHE_PART  <- TRUE       # TNW / WIC / BIC at the set level
ADD_IND_NICHE   <- FALSE      # Stomach-level WIC/TNW; too heavy for the full sweep
DRIVER_PUNI     <- "adjusted" # Step-down adjusted p-values across taxa
DRIVER_RELATIVE <- TRUE       # Standardize abundance for IndVal

# Spatial level for this session
SPATIAL_LEVEL <- "stratum"

# Source the engine (defines objects only, runs nothing)
if (!exists("PREY_FAMILY")) PREY_FAMILY <- "_1"   # "_1" (main), "_2" or "_PP"
source("6a_Engine_Trophic.R")

# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("\nLoading data...\n")
load("data/dat_classed.rda")
diet_clean <- dat_classed %>%
  mutate(processing_date = as.Date(processing_date)) %>%
  as.data.frame()

# The engine stops on a missing spatial column; this reports the available
# names instead. The NA count is printed because every level assumes it is
# zero — a non-zero value means the levels no longer share one sample.
if (!SPATIAL_SOURCE %in% names(diet_clean)) {
  stop("Column '", SPATIAL_SOURCE, "' not found in diet_clean. Available: ",
       paste(names(diet_clean), collapse = ", "))
}
diet_clean[[SPATIAL_SOURCE]] <- as.factor(diet_clean[[SPATIAL_SOURCE]])
cat("Rows with a missing ", SPATIAL_SOURCE, ": ",
    sum(is.na(diet_clean[[SPATIAL_SOURCE]])), "\n", sep = "")


# =============================================================================
# 2b. SAMPLE-SIZE CHECK BEFORE COMMITTING TO THE SWEEP
# =============================================================================
# Stomachs per spatial unit per year. A unit below N_MIN in either period
# yields no cell at all, so this table sets the ceiling on what the sweep can
# return. Worth reading before launching a run that takes hours.
cat("\n--- Stomachs per stratum x year (scenario window) ---\n")
print(
  diet_clean %>%
    filter(year %in% c(2004, 2005, 2006, 2018, 2019)) %>%
    count(.data[[SPATIAL_SOURCE]], year) %>%
    tidyr::pivot_wider(names_from = year, values_from = n, values_fill = 0),
  n = Inf
)


# Regression guard on the biomass currency. `somatic_wt_g` holds the PREY
# weight; this confirms the upstream build still delivers it that way. Expect a
# low constant-within-stomach share and a flat slope against predator length.
check_prey_weight_column(diet_clean)

# =============================================================================
# 3. SCENARIOS (Temporal Period Maps)
# =============================================================================
period_map_P1 <- list("2004"      = c(2004),             "2006"      = c(2006))
period_map_P2 <- list("2018"      = c(2018),             "2019"      = c(2019))
period_map_P3 <- list("2004-2005" = c(2004, 2005),       "2006"      = c(2006))
period_map_P4 <- list("2006"      = c(2006),             "2018"      = c(2018))
period_map_PT <- list("2004-2006" = c(2004, 2005, 2006), "2018-2019" = c(2018, 2019))

# =============================================================================
# 4. BUILD DATASETS
# =============================================================================
# One dataset per scenario. The spatial treatment is applied inside
# run_pipeline(), from SPATIAL_LEVEL, so no second pooled copy is needed.
cat("Building classed datasets...\n")
dat_P1 <- make_dat_classed(diet_clean, period_map = period_map_P1)
dat_P2 <- make_dat_classed(diet_clean, period_map = period_map_P2)
dat_P3 <- make_dat_classed(diet_clean, period_map = period_map_P3)
dat_P4 <- make_dat_classed(diet_clean, period_map = period_map_P4)
dat_PT <- make_dat_classed(diet_clean, period_map = period_map_PT)

# =============================================================================
# 5. OUTPUT FOLDERS
# =============================================================================
# Named after the level, so the three parallel sessions cannot overwrite each
# other. Filenames also carry the level, so merging the folders later is safe.
out_dir_data  <- paste0("Sensitivity_stratum", PREY_FAMILY)
out_dir_plots <- paste0("Sensitivity_Plot_stratum", PREY_FAMILY)

dir.create(out_dir_data,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 6. RUNS (Engine Execution)
# =============================================================================
# do_driver_prey = FALSE for the sweep: the manylm bootstrap dominates runtime
# across ~100 taxonomic resolutions. It is turned on once, at the single
# retained resolution, in section 9.
cat("\nRunning Analysis Pipeline (stratum)...\n")

res_biomass_P1 <- run_save_plot(dat_P1, "P1", "biomass",    "2004", "2006", out_dir = out_dir_data, plot_dir = out_dir_plots)
res_occ_P1     <- run_save_plot(dat_P1, "P1", "occurrence", "2004", "2006", out_dir = out_dir_data, plot_dir = out_dir_plots)

res_biomass_P2 <- run_save_plot(dat_P2, "P2", "biomass",    "2018", "2019", out_dir = out_dir_data, plot_dir = out_dir_plots)
res_occ_P2     <- run_save_plot(dat_P2, "P2", "occurrence", "2018", "2019", out_dir = out_dir_data, plot_dir = out_dir_plots)

res_biomass_P3 <- run_save_plot(dat_P3, "P3", "biomass",    "2004-2005", "2006", out_dir = out_dir_data, plot_dir = out_dir_plots)
res_occ_P3     <- run_save_plot(dat_P3, "P3", "occurrence", "2004-2005", "2006", out_dir = out_dir_data, plot_dir = out_dir_plots)

res_biomass_P4 <- run_save_plot(dat_P4, "P4", "biomass",    "2006", "2018", out_dir = out_dir_data, plot_dir = out_dir_plots)
res_occ_P4     <- run_save_plot(dat_P4, "P4", "occurrence", "2006", "2018", out_dir = out_dir_data, plot_dir = out_dir_plots)

res_biomass_PT <- run_save_plot(dat_PT, "PT", "biomass",    "2004-2006", "2018-2019", out_dir = out_dir_data, plot_dir = out_dir_plots)
res_occ_PT     <- run_save_plot(dat_PT, "PT", "occurrence", "2004-2006", "2018-2019", out_dir = out_dir_data, plot_dir = out_dir_plots)

# =============================================================================
# 6b. DIAGNOSTICS
# =============================================================================
# audit_set_depth() comes first. With the trawl set as the replicate unit, the
# number of tows per cell is the binding constraint on every test downstream.
cat("\n--- Set depth (the binding constraint) ---\n")
print(audit_set_depth(res_biomass_PT$results))

cat("\n--- Untestable cells ---\n")
print(audit_untestable(res_biomass_PT$results))

cat("\n--- H' / Bs overlap, within resolution ---\n")
print(audit_signal_overlap(res_biomass_PT$results))

cat("\n--- Rank-test significance floor ---\n")
print(bs_test_floor())

# Effect sizes alongside detection rates. BC and R2 do not depend on
# significance, so comparing these across the three levels separates a
# difference in effect size from a difference in power.
cat("\n--- Effect size vs detection rate, PT ---\n")
print(
  res_biomass_PT$results %>%
    filter(testable) %>%
    summarise(
      n            = dplyr::n(),
      k_med        = median(pmin(n_set_P1, n_set_P2), na.rm = TRUE),
      BC_med       = median(BC, na.rm = TRUE),
      R2_med       = median(R2_comp, na.rm = TRUE),
      pct_H_sig    = round(100 * mean(p_H    < ALPHA, na.rm = TRUE), 1),
      pct_Bs_sig   = round(100 * mean(p_Bs   < ALPHA, na.rm = TRUE), 1),
      pct_comp_sig = round(100 * mean(p_comp < ALPHA, na.rm = TRUE), 1)
    )
)

cat("\n--- Cells lost to insufficient set depth ---\n")
# The count that decides whether the stratum grain is usable at all.
print(
  res_biomass_PT$results %>%
    summarise(
      n_cells        = dplyr::n(),
      pct_testable   = round(100 * mean(testable), 1),
      pct_reliable   = round(100 * mean(reliable), 1),
      pct_below_sets = round(100 * mean(n_set_P1 < MIN_SETS | n_set_P2 < MIN_SETS), 1)
    )
)

# =============================================================================
# 7. FINAL SAVE
# =============================================================================
cat("\nSaving all runs to a single file...\n")
dir.create("data/Sensitivity", recursive = TRUE, showWarnings = FALSE)
save(
  res_biomass_P1, res_occ_P1,
  res_biomass_P2, res_occ_P2,
  res_biomass_P3, res_occ_P3,
  res_biomass_P4, res_occ_P4,
  res_biomass_PT, res_occ_PT,
  file = paste0("data/Sensitivity/all_runs_stratum", PREY_FAMILY, ".rda")
)

# =============================================================================
# 8. QUICK TEST PLOT
# =============================================================================
# Uses the object already in memory rather than re-running the pipeline, so the
# figure shown is the one whose .rda was just written.
cat("\nGenerating test plot for P2 biomass...\n")

p <- make_diag_plot2(
  df           = res_biomass_P2$results,
  inventory_df = res_biomass_P2$inventory,
  mode_label   = "biomass",
  period_1     = "2018",
  period_2     = "2019",
  by_area      = TRUE,
  title_main   = "Biomass dietary transitions by stratum: 2018 vs 2019",
  y_lab        = "Within-unit frequency (%)"
)
print(p)

# =============================================================================
# 9. OPTIONAL PASS AT THE FINAL RESOLUTION ONLY
# =============================================================================
# Both switches below are too expensive for the full sweep. They are meant to be
# run once, on a dataset reduced to the single taxonomic resolution retained for
# the manuscript, to produce the driver-prey table and the individual-level
# WIC/TNW figures.
#
#   FINAL_X <- "prey_category_25_2"   # replace with the retained resolution
#   drop <- setdiff(grep("prey_category_\\d+_2", names(dat_PT), value = TRUE), FINAL_X)
#   dat_final <- dplyr::select(dat_PT, -dplyr::all_of(drop))
#
#   ADD_IND_NICHE <- TRUE
#   res_final <- run_pipeline(dat_final, "biomass", "2004-2006", "2018-2019",
#                             do_driver_prey = TRUE)
#
#   res_final$results %>%
#     dplyr::filter(testable) %>%
#     dplyr::select(species, size_class, diagnostic,
#                   WIC_TNW_ind_P1, WIC_TNW_ind_P2,
#                   BIC_within_set_P1, BIC_among_set_P1, driver_prey)

# Once all three levels have finished, 6e_Compare_Levels.R reads the three
# .rda files and builds the cross-level comparison.

cat("\nPipeline execution complete (stratum)!\n")
