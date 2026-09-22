# =============================================================================
# 6c_Run_L2_ecoregion.R - L2 ecoregion (cells split by ecoregion)
# -----------------------------------------------------------------------------
# One of three run scripts for 6a_Engine_Trophic.R (6b all_gulf, 6c ecoregion,
# 6d stratum). They differ only in SPATIAL_LEVEL and must each run in their own
# R session. Results: Sensitivity_ecoregion<family>/ and
# data/Sensitivity/all_runs_ecoregion<family>.rda. 6e_Run_families.R reruns this
# script for the "_2" and "_PP" families.
# =============================================================================
rm(list = setdiff(ls(), "PREY_FAMILY"))

if (!exists("PREY_FAMILY")) PREY_FAMILY <- "_1"   # "_1" (manuscript), "_2" or "_PP"
SPATIAL_LEVEL <- "ecoregion"

N_MIN    <- 5
N_STRICT <- 25
MIN_SETS <- 3
R_PERM   <- 999
ALPHA    <- 0.05
RUN_TAG         <- "_T"

source("6a_Engine_Trophic.R")
run_all_scenarios(data_path = "data/dat_classed.rda", tag = RUN_TAG)
