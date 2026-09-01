# =============================================================================
# 6e_Run_families.R - SENSITIVITY RUNS FOR THE ALTERNATIVE GROUPING FAMILIES
# -----------------------------------------------------------------------------
# The manuscript uses the "_1" family (default of 6b/6c/6d: just run them).
# This driver reruns the same three run scripts for the alternative families,
# "_2" (pooled, q = 2) and "_PP" (per-predator). Outputs land beside the _1 ones:
#   Sensitivity_<level>_2/ , Sensitivity_<level>_PP/ (+ _Plot_ folders)
#   data/Sensitivity/all_runs_<level>_2.rda , ..._PP.rda
# 6g_Compare_families.R then reads the three families and builds Appendix B.
#
# "_PP" requires 3PP_taxonomic_groups.R (writes data/prey_groups_PP.RData) and
# a rerun of 4 -> 5 so that dat_classed carries the prey_category_<x>_PP
# columns. "_2" needs nothing: its columns are already in dat_classed.
#
# Each run is launched in a fresh Rscript process, so no state leaks between
# families and a crash in one run does not kill the others.
# =============================================================================

FAMILIES <- c("_2", "_PP")          # "_1" is the main analysis, run directly
RUN_SCRIPTS <- c("6b_Run_L1_all_gulf.R",
                 "6c_Run_L2_ecoregion.R",
                 "6d_Run_L3_stratum.R")

rscript <- file.path(R.home("bin"), "Rscript")

for (fam in FAMILIES) {
  # sanity: does dat_classed carry this family's columns?
  e <- new.env(); load("data/dat_classed.rda", envir = e)
  obj <- get(ls(e)[1], envir = e)
  n_cols <- sum(grepl(paste0("^prey_category_\\d+", fam, "$"), names(obj)))
  rm(e, obj); gc()
  if (n_cols == 0) {
    warning("dat_classed has no prey_category_*", fam, " columns - family ",
            fam, " skipped. (For _PP: run 3PP, then 4 -> 5.)", call. = FALSE)
    next
  }
  cat("\n=========", fam, ":", n_cols, "resolution columns =========\n")

  for (rs in RUN_SCRIPTS) {
    cat("--", fam, rs, "\n")
    status <- system2(
      rscript,
      c("-e", shQuote(sprintf('PREY_FAMILY <- "%s"; source("%s")', fam, rs)))
    )
    if (status != 0) warning("FAILED: ", rs, " for family ", fam, call. = FALSE)
  }
}

