# =============================================================================
# 10_Prey_Drivers.R - PREY ASSOCIATED WITH THE INTER-DECADE TURNOVER
#                     (cell-consistent design)
# =============================================================================
# PURPOSE
#   Scripts 6b-6d establish THAT the diets reorganised between 2004-2006 and
#   2018-2019 (Substitution is the modal diagnostic at every scale, in both
#   currencies). This script asks WHICH prey categories carry that turnover -
#   which lost ground, which gained - and whether the answer survives the
#   choices that could manufacture it: the predator mix of the samples, body
#   size, ecoregion, taxonomic resolution and diet currency.
#
#   Everything here is ASSOCIATION, not mechanism. Nothing holds prey
#   availability or the environment constant, so a prey named below is a
#   candidate to be read against independent evidence on the Gulf prey field
#   (Savenkoff et al. 2007; Benoit & Swain 2008), not a demonstrated cause.
#
# WHY THE DESIGN CHANGED
#   The earlier draft pooled every stomach of a trawl set into one "set diet",
#   whatever predator had eaten the prey. That quantity depends on WHICH
#   predators were sampled in the set, and the predator mix of the survey
#   changed between decades; a prey "gained" at the set level can therefore be
#   an artefact of sampling rather than a change of diet. The manuscript's own
#   framework compares each predator x size x ecoregion CELL with itself for
#   exactly that reason, and this script now does the same.
#
# DESIGN
#   Row       = one predator species, one size class, in one trawl set: the
#               prey of its stomachs summed (biomass) or counted (occurrence).
#               The set remains the replicate (Hurlbert 1984); the predator is
#               no longer mixed.
#   Cell      = predator x size class x ecoregion, as in scripts 6b-6d. Only
#               cells sampled in BOTH periods with >= MIN_ROWS rows per period
#               enter the tests (the framework's testable-cell rule).
#
#   TEST A - covariate-adjusted multivariate linear model (mvabund; Warton et
#            al. 2012; Wang et al. 2012). Hellinger-transformed profiles
#            (Legendre & Gallagher 2001) ~ cell + period, the period term tested
#            AFTER the cell term, by permutation of residuals RESTRICTED WITHIN
#            CELLS (rows of one cell are exchangeable under H0). Univariate
#            p-values adjusted step-down across prey (Westfall & Young 1993).
#            Answers: at equal predator, size and ecoregion, which prey differ
#            between decades? The same model without the cell term is fitted
#            for comparison, so the size of the sampling-mix artefact is visible
#            (column mglm_p_unadjusted).
#
#   TEST B - cell-level consistency. For every cell, the mean relative profile
#            in each period and its change; across cells, a Wilcoxon signed-rank
#            test per prey (each cell one paired observation; Benjamini &
#            Hochberg 1995 across prey), the share of cells moving in the modal
#            direction, and an Indicator Value analysis within each cell
#            (Dufrene & Legendre 1997) summarised as the number of cells in
#            which the prey indicates one decade. Answers: is the change the
#            same thing in most cells, or the doing of a few?
#
#   ROBUSTNESS LAYERS
#     by axis         Test B summarised within each ecoregion, size class and
#                     predator (no new model: a reading of the cell table).
#     by resolution   Tests A (fewer permutations) and B at >= 10 taxonomic
#                     resolutions (taxonomic sufficiency: Ellis 1985; Warwick
#                     1988). A prey significant at one resolution only, or whose
#                     name exists in one decade only (identification practice
#                     changed: Gammaridae -> Amphipoda, Crustacea ->
#                     Arthropoda), is exposed here.
#     cross-currency  Biomass and occurrence tables joined (Hyslop 1980; Baker
#                     et al. 2014): a prey moving the same way, significantly,
#                     in both currencies is the strongest evidence produced.
#
#   CONTEXT (secondary to the cell framework of Methods 2.4)
#     PERMANOVA on the rows, D ~ predator + size_class + Area + period, period
#     last, permutations within cells (Anderson 2001), PERMDISP on period
#     (Anderson 2006; Anderson & Walsh 2013); NMDS of the cell centroids with
#     an arrow per cell across decades and a blocked PERMANOVA (one block per
#     cell); the robust prey as fitted vectors. Stress read as in Clarke (1993).
#
# WHAT WAS DROPPED, AND WHY
#   The pooled set-level test (predator mix confound); the set-level NMDS
#   (same confound, stress ~0.3); SIMPER (confounds mean and variance: Warton
#   et al. 2012).
#
# CRASH PROTECTION
#   mvabund's resampling runs in compiled code and can take the R session down
#   (a segfault, not an R error) on small or degenerate matrices. Every model
#   fit runs in a child R process (package callr); a crash comes back as an
#   error, is logged in driver_test_failures.txt, and the script carries on.
#
# OUTPUTS  (Output_Drivers/, DIR_DRIVERS from Config_Mappings.R)
#   drivers_<cur>.csv                    THE driver table (Table 4). One row per
#                                        prey: adjusted and unadjusted MGLM p,
#                                        cell-mean share per decade and change,
#                                        n cells gained / lost, % in the modal
#                                        direction, Wilcoxon BH p, IndVal cell
#                                        counts, criteria met.
#   drivers_cells_<cur>.csv              Per cell x prey: share per period and
#                                        change, IndVal within the cell.
#   drivers_by_axis_<cur>.csv            Test B within ecoregion / size / predator.
#   drivers_resolution_sweep_<cur>.csv   Driver table at every sweep resolution.
#   drivers_resolution_stability_<cur>.csv  Per prey across resolutions.
#   drivers_cross_currency.csv           Biomass x occurrence join, robustness class.
#   Fig_prey_turnover_<cur>              Cell-mean share per decade (dumbbell).
#   Fig_prey_cell_consistency_<cur>      % of cells where each prey gained / lost.
#   Fig_resolution_stability_<cur>       Prey x resolution tiles.
#   Fig_nmds_cells_<cur>                 Cell-centroid NMDS with decade arrows.
#   permanova_rows_<cur>.txt, permdisp_period_<cur>.txt,
#   permanova_cellpaired_<cur>.txt, driver_test_failures.txt, run_summary.txt
#
# REFERENCES
#   Anderson MJ (2001) Austral Ecol 26:32-46.  Anderson MJ (2006) Biometrics
#   62:245-253.  Anderson MJ, Walsh DCI (2013) Ecol Monogr 83:557-574.
#   Baker R, Buckland A, Sheaves M (2014) Fish Fish 15:170-177.
#   Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
#   Benoit HP, Swain DP (2008) Can J Fish Aquat Sci 65:2088-2104.
#   Clarke KR (1993) Aust J Ecol 18:117-143.
#   Dufrene M, Legendre P (1997) Ecol Monogr 67:345-366.
#   Ellis D (1985) Mar Pollut Bull 16:459.  Hurlbert SH (1984) Ecol Monogr 54:187-211.
#   Hyslop EJ (1980) J Fish Biol 17:411-429.
#   Legendre P, Gallagher ED (2001) Oecologia 129:271-280.
#   Oksanen J et al. (2022) vegan: Community Ecology Package.
#   Savenkoff C et al. (2007) Estuar Coast Shelf Sci 73:711-724.
#   Wang Y, Naumann U, Wright ST, Warton DI (2012) Methods Ecol Evol 3:471-474.
#   Warton DI, Wright ST, Wang Y (2012) Methods Ecol Evol 3:89-101.
#   Warwick RM (1988) Mar Pollut Bull 19:259-268.
#   Westfall PH, Young SS (1993) Resampling-based multiple testing. Wiley.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ggplot2)
  library(dplyr); library(tidyr); library(readr)
  library(mvabund); library(labdsv)
})

PREY_FAMILY <- "_1"          # prey-grouping family used for the manuscript
source("R_helpers/Config_Mappings.R")

# =============================================================================
# CONFIGURATION
# =============================================================================
CFG <- list(
  rda_path  = "data/dat_classed.rda",
  obj_name  = "dat_classed",

  col_pred   = "predator_species_common_name",
  col_size   = "size_class",
  col_area   = "Area",
  col_stom   = "stomach_id",
  col_year   = "year",
  col_vessel = "vessel.code",
  col_set    = "set",
  col_wt     = "somatic_wt_g",     # PREY item weight (one row per prey item)

  prey_family = PREY_FAMILY,
  res_col     = paste0("prey_category_440", PREY_FAMILY),   # manuscript resolution
  sweep_x     = c(50, 100, 150, 200, 250, 300, 350, 440, 500, 600, 750, 900, 1000),
  sweep_min   = 10,

  yr_early  = 2004:2006,  yr_late = 2018:2019,
  lab_early = "2004-2006", lab_late = "2018-2019",

  currencies = c("biomass", "occurrence"),
  drop_prey  = c("Empty", "empty", "Unidentified", "unidentified",
                 "Digested", "digested", "NA"),

  # Cells and rows
  min_rows_per_period = 3,   # testable cell: >= 3 rows (sets) in EACH period (= MIN_SETS of the engine)
  min_stom_per_row    = 1,   # rows (set x predator x size) built on fewer stomachs are dropped;
  # 1 keeps everything (as the engine), 2-3 trades rows for less noise
  min_rows_per_prey   = 3,   # a prey must occur in at least this many rows
  min_cells_wilcox    = 5,   # Wilcoxon across cells needs at least this many cells

  # Resampling
  seed        = 42,
  perms       = 999,   # Test A at the manuscript resolution
  perms_cell  = 499,   # IndVal within cells
  perms_sweep = 199,   # Test A at the sweep resolutions (Test B is instant)
  sweep_run_mglm = TRUE,     # FALSE: sweep on Test B only (minutes instead of hours)
  alpha       = 0.05,
  modal_share = 0.60,  # Test B: >= 60 % of cells in the same direction

  top_n_figs = 20,
  top_n_vectors = 10,

  # Crash protection (see header)
  isolate_tests = TRUE,
  test_timeout  = 3600
)

set.seed(CFG$seed)
PER_COL  <- PERIOD_PAL
DIR_COL  <- c(gained = "#d7191c", lost = "#2c7bb6")
save_drv <- function(p, stem, w, h, dpi = 320) save_fig(p, stem, w, h, dpi, dir = DIR_DRIVERS)

FAIL_LOG <- file.path(DIR_DRIVERS, "driver_test_failures.txt")
if (file.exists(FAIL_LOG)) file.remove(FAIL_LOG)
log_failure <- function(label, msg) {
  message("  FAILED '", label, "': ", msg)
  cat(format(Sys.time()), " | ", label, " | ", msg, "\n", sep = "", file = FAIL_LOG, append = TRUE)
}
SUMMARY <- list()
note <- function(cur, key, value) SUMMARY[[cur]][[key]] <<- value

# =============================================================================
# 1. LOAD
# =============================================================================
cat("\nLoading ", CFG$rda_path, "\n", sep = "")
e <- new.env(); loaded <- load(CFG$rda_path, envir = e)
obj <- if (CFG$obj_name %in% loaded) CFG$obj_name else loaded[1]
DT  <- get(obj, envir = e); rm(e)
if (inherits(DT, "sf")) DT <- sf::st_drop_geometry(DT)
DT  <- as.data.table(DT)

for (nm in c(CFG$col_pred, CFG$col_size, CFG$col_area, CFG$col_stom,
             CFG$col_year, CFG$col_vessel, CFG$col_set, CFG$col_wt))
  if (!nm %in% names(DT)) stop("Column '", nm, "' declared in CFG is not in the data.")

DT[, period := fifelse(get(CFG$col_year) %in% CFG$yr_early, CFG$lab_early,
                       fifelse(get(CFG$col_year) %in% CFG$yr_late,  CFG$lab_late, NA_character_))]
DT <- DT[!is.na(period) & !is.na(get(CFG$col_pred)) & !is.na(get(CFG$col_size)) &
           !is.na(get(CFG$col_area))]
DT[, period := factor(period, levels = c(CFG$lab_early, CFG$lab_late))]
DT[, set_uid := paste(get(CFG$col_year), get(CFG$col_vessel), get(CFG$col_set), sep = "_")]
DT[, `:=`(predator = as.character(get(CFG$col_pred)),
          size_cl  = as.character(get(CFG$col_size)),
          Area     = as.character(get(CFG$col_area)))]
DT[, cell := paste(predator, size_cl, Area, sep = "|")]
DT[, row_id := paste(set_uid, predator, size_cl, sep = "#")]

# Resolutions available and sweep actually run
avail_x <- sort(as.integer(sub(paste0("^prey_category_(\\d+)", CFG$prey_family, "$"), "\\1",
                               grep(paste0("^prey_category_\\d+", CFG$prey_family, "$"),
                                    names(DT), value = TRUE))))
if (!CFG$res_col %in% names(DT)) stop("CFG$res_col '", CFG$res_col, "' not in the data.")
sweep_x <- intersect(CFG$sweep_x, avail_x)
if (length(sweep_x) < CFG$sweep_min && length(avail_x) >= CFG$sweep_min)
  sweep_x <- avail_x[unique(round(seq(1, length(avail_x), length.out = CFG$sweep_min)))]
sweep_x <- sort(union(sweep_x, 440L))
cat("Resolution sweep: x = ", paste(sweep_x, collapse = ", "), " (", length(sweep_x), ")\n", sep = "")
cat("Prey records: ", nrow(DT), " | sets: ", uniqueN(DT$set_uid),
    " | rows (set x predator x size): ", uniqueN(DT$row_id),
    " | cells: ", uniqueN(DT$cell), "\n", sep = "")

# -----------------------------------------------------------------------------
# Row matrix builder: rows = set x predator x size class, columns = prey.
# Returns the matrix, its metadata (aligned) and the testable-cell filter.
# -----------------------------------------------------------------------------
build_rows <- function(dat, currency, res_col) {
  d <- dat[, c("row_id", "cell", "set_uid", "predator", "size_cl", "Area", "period",
               CFG$col_stom, CFG$col_wt, res_col), with = FALSE]
  d[, prey := as.character(get(res_col))]
  d <- d[!is.na(prey) & !(prey %chin% CFG$drop_prey)]
  if (!nrow(d)) return(NULL)

  agg <- if (currency == "biomass") {
    d[, .(val = sum(get(CFG$col_wt), na.rm = TRUE)), by = .(row_id, prey)]
  } else {
    d[, .(val = as.numeric(uniqueN(get(CFG$col_stom)))), by = .(row_id, prey)]
  }
  m   <- dcast(agg, row_id ~ prey, value.var = "val", fill = 0)
  mat <- as.matrix(m[, -1]); rownames(mat) <- m$row_id
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  meta <- unique(d[, .(row_id, cell, set_uid, predator, size_cl, Area, period)])
  meta <- meta[match(rownames(mat), row_id)]
  meta[, n_stom := d[, uniqueN(get(CFG$col_stom)), by = row_id][match(meta$row_id, row_id), V1]]
  keep0 <- meta$n_stom >= CFG$min_stom_per_row
  mat <- mat[keep0, , drop = FALSE]; meta <- meta[keep0]

  # Testable cells: both periods, >= min rows in each.
  ok <- meta[, .(n1 = sum(period == CFG$lab_early), n2 = sum(period == CFG$lab_late)), by = cell][
    n1 >= CFG$min_rows_per_period & n2 >= CFG$min_rows_per_period, cell]
  keep <- meta$cell %chin% ok
  mat <- mat[keep, , drop = FALSE]; meta <- meta[keep]
  mat <- mat[, colSums(mat > 0) >= CFG$min_rows_per_prey, drop = FALSE]
  mat <- mat[, apply(mat, 2, function(v) length(unique(v)) > 1), drop = FALSE]
  keep2 <- rowSums(mat) > 0
  mat <- mat[keep2, , drop = FALSE]; meta <- meta[keep2]
  meta[, cell := factor(cell)]
  meta[, period := factor(as.character(period), levels = c(CFG$lab_early, CFG$lab_late))]
  list(mat = mat, meta = meta, n_cells = uniqueN(meta$cell))
}

# Permutation ids restricted within cells (rows of one cell exchangeable under
# H0 of no period effect). Base R, 1-based, one row per resample.
block_perm_ids <- function(cell, n_perm, seed) {
  set.seed(seed)
  idx <- seq_along(cell); spl <- split(idx, cell)
  t(vapply(seq_len(n_perm), function(i) {
    out <- idx
    for (s in spl) if (length(s) > 1) out[s] <- s[sample.int(length(s))]
    out
  }, integer(length(cell))))
}

# =============================================================================
# 2. TEST A - covariate-adjusted MGLM (child process)
# =============================================================================
core_mglm_adjusted <- function(comm, cell, period, lev, n_perm, seed, with_unadjusted) {
  suppressPackageStartupMessages({ library(vegan); library(mvabund) })
  cell   <- factor(cell)
  period <- factor(as.character(period), levels = lev)
  mv  <- mvabund::mvabund(as.matrix(vegan::decostand(comm, method = "hellinger")))

  # ids: within-cell permutations (1-based, one row per resample)
  set.seed(seed)
  idx <- seq_along(cell); spl <- split(idx, cell)
  ids <- t(vapply(seq_len(n_perm), function(i) {
    out <- idx
    for (s in spl) if (length(s) > 1) out[s] <- s[sample.int(length(s))]
    out }, integer(length(cell))))

  fit <- mvabund::manylm(mv ~ cell + period)
  restricted <- TRUE
  an <- tryCatch(
    anova(fit, resamp = "perm.resid", test = "F", p.uni = "adjusted", bootID = ids),
    error = function(e) {
      # Older mvabund without bootID: unrestricted residual permutation. The
      # cell term is still in the model, so the period effect remains adjusted;
      # only the exchangeability set is wider. Flagged in the output.
      restricted <<- FALSE
      anova(fit, resamp = "perm.resid", test = "F", p.uni = "adjusted", nBoot = n_perm)
    })
  uni_p  <- an$uni.p["period", ]; uni_p <- uni_p[is.finite(uni_p)]
  mult_p <- tryCatch(an$table["period", "Pr(>F)"], error = function(e) NA_real_)
  # Effect direction and size: model coefficient of the period term per prey
  # (Hellinger scale), sign only is used downstream.
  coef_period <- tryCatch(coef(fit)[grep("^period", rownames(coef(fit)))[1], ],
                          error = function(e) rep(NA_real_, ncol(comm)))

  out <- list(uni_p = uni_p, mult_p = mult_p, coef_period = coef_period,
              restricted_permutations = restricted)

  if (isTRUE(with_unadjusted)) {
    set.seed(seed)
    fit0 <- mvabund::manylm(mv ~ period)
    an0  <- anova(fit0, resamp = "perm.resid", test = "F", p.uni = "adjusted", nBoot = n_perm)
    u0 <- an0$uni.p["period", ]; out$uni_p_unadj <- u0[is.finite(u0)]
    out$mult_p_unadj <- tryCatch(an0$table["period", "Pr(>F)"], error = function(e) NA_real_)
  }
  out
}

# IndVal within each cell (child process). Returns one row per cell x prey.
core_cell_indval <- function(comm, cell, period, lev, n_perm, seed, min_rows) {
  suppressPackageStartupMessages({ library(vegan); library(labdsv) })
  set.seed(seed)
  cell <- as.character(cell); period <- factor(as.character(period), levels = lev)
  rel  <- vegan::decostand(comm, "total")
  out  <- list()
  for (cl in unique(cell)) {
    i <- which(cell == cl); g <- droplevels(period[i])
    if (nlevels(g) != 2 || any(table(g) < min_rows)) next
    m <- rel[i, , drop = FALSE]
    m <- m[, colSums(m > 0) >= 2, drop = FALSE]
    if (ncol(m) < 2) next
    iv <- tryCatch(labdsv::indval(as.data.frame(m), as.integer(g), numitr = n_perm),
                   error = function(e) NULL)
    if (is.null(iv)) next
    out[[cl]] <- data.frame(cell = cl, prey = names(iv$indcls),
                            indval = as.numeric(iv$indcls), indval_p = as.numeric(iv$pval),
                            indval_grp = levels(g)[iv$maxcls], stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

run_isolated <- function(core, args, label) {
  t0 <- Sys.time()
  res <- tryCatch({
    if (isTRUE(CFG$isolate_tests)) {
      if (!requireNamespace("callr", quietly = TRUE))
        stop("CFG$isolate_tests = TRUE needs package 'callr': install.packages('callr')")
      callr::r(function(core, args) do.call(core, args), args = list(core = core, args = args),
               timeout = CFG$test_timeout, show = FALSE)
    } else do.call(core, args)
  }, error = function(err) { log_failure(label, conditionMessage(err)); NULL })
  cat(sprintf("    %s: %.1f min\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}

# =============================================================================
# 3. TEST B - cell-level consistency (in-session, cheap)
# =============================================================================
cell_shares <- function(comm, meta) {
  rel <- decostand(comm, "total")
  dt  <- as.data.table(rel); dt[, `:=`(cell = as.character(meta$cell), period = as.character(meta$period))]
  long <- melt(dt, id.vars = c("cell", "period"), variable.name = "prey", value.name = "share",
               variable.factor = FALSE)
  # mean relative share per cell x period, then the paired change per cell
  cp <- long[, .(share = mean(share), n_rows = .N), by = .(cell, period, prey)]
  w  <- dcast(cp, cell + prey ~ period, value.var = c("share", "n_rows"))
  setnames(w, paste0("share_", CFG$lab_early),  "share_early")
  setnames(w, paste0("share_", CFG$lab_late),   "share_late")
  setnames(w, paste0("n_rows_", CFG$lab_early), "n_rows_early")
  setnames(w, paste0("n_rows_", CFG$lab_late),  "n_rows_late")
  w[, d_share := share_late - share_early]
  w[]
}

cell_consistency <- function(cs, iv_cells, alpha = CFG$alpha) {
  # Unweighted: every cell counts once (the consistency question). Weighted:
  # cells weighted by their sampling effort, min(rows early, rows late), so
  # that a cell built on 3 + 3 sets does not move the mean share as much as one
  # built on 40 + 40 (the descriptive question).
  x <- cs[, .(
    n_cells   = .N,
    n_gain    = sum(d_share > 0), n_loss = sum(d_share < 0), n_flat = sum(d_share == 0),
    share_early = mean(share_early), share_late = mean(share_late),
    d_share   = mean(d_share),
    share_early_w = weighted.mean(share_early, pmin(n_rows_early, n_rows_late)),
    share_late_w  = weighted.mean(share_late,  pmin(n_rows_early, n_rows_late)),
    d_share_w     = weighted.mean(d_share,     pmin(n_rows_early, n_rows_late)),
    d_share_median = median(d_share),
    wilcox_p  = if (.N >= CFG$min_cells_wilcox && any(d_share != 0))
      suppressWarnings(wilcox.test(d_share, mu = 0, exact = FALSE)$p.value) else NA_real_
  ), by = prey]
  x[, direction := fifelse(d_share > 0, "gained", "lost")]
  x[, pct_modal := 100 * pmax(n_gain, n_loss) / n_cells]
  x[, wilcox_p_bh := p.adjust(wilcox_p, "BH")]
  if (!is.null(iv_cells) && nrow(iv_cells)) {
    ivs <- as.data.table(iv_cells)[indval_p <= alpha, .(
      iv_cells_early = sum(indval_grp == CFG$lab_early),
      iv_cells_late  = sum(indval_grp == CFG$lab_late)), by = prey]
    x <- merge(x, ivs, by = "prey", all.x = TRUE)
    x[is.na(iv_cells_early), iv_cells_early := 0L][is.na(iv_cells_late), iv_cells_late := 0L]
  } else x[, `:=`(iv_cells_early = NA_integer_, iv_cells_late = NA_integer_)]
  x[]
}

# Assemble the driver table from Test A and Test B.
assemble_drivers <- function(A, B, comm, label) {
  d <- copy(B)
  d[, mglm_p := if (!is.null(A)) unname(A$uni_p[prey]) else NA_real_]
  d[, mglm_coef_sign := if (!is.null(A)) sign(unname(A$coef_period[prey])) else NA_real_]
  d[, mglm_p_unadjusted := if (!is.null(A$uni_p_unadj)) unname(A$uni_p_unadj[prey]) else NA_real_]
  d[, mglm_multivariate_p := if (!is.null(A)) A$mult_p else NA_real_]
  d[, A_sig := !is.na(mglm_p) & mglm_p <= CFG$alpha]
  d[, B_sig := !is.na(wilcox_p_bh) & wilcox_p_bh <= CFG$alpha & pct_modal >= 100 * CFG$modal_share]
  d[, n_criteria := as.integer(A_sig) + as.integer(B_sig)]
  d[, level := label]
  d[, abs_d := abs(d_share)]; setorder(d, -n_criteria, -abs_d); d[, abs_d := NULL]
  setcolorder(d, c("prey", "direction", "n_criteria", "A_sig", "B_sig", "mglm_p",
                   "mglm_p_unadjusted", "mglm_multivariate_p", "share_early", "share_late",
                   "d_share", "d_share_median", "share_early_w", "share_late_w", "d_share_w",
                   "n_cells", "n_gain", "n_loss", "n_flat",
                   "pct_modal", "wilcox_p", "wilcox_p_bh", "iv_cells_early", "iv_cells_late"))
  d[]
}

check_stress <- function(nm, label, n_points) {
  s <- nm$stress
  if (is.finite(s) && s < 0.001 && n_points > 10)
    warning("Degenerate NMDS for ", label, ": stress = ", signif(s, 3), call. = FALSE)
  else if (is.finite(s) && s > 0.2)
    warning("High NMDS stress for ", label, ": ", round(s, 3), " (Clarke 1993).", call. = FALSE)
  invisible(s)
}

# =============================================================================
# 4. RUN, PER CURRENCY
# =============================================================================
DRV <- list()

for (cur in CFG$currencies) {
  cat("\n", strrep("=", 70), "\n", cur, "\n", strrep("=", 70), "\n", sep = "")
  SUMMARY[[cur]] <- list()

  R <- build_rows(DT, cur, CFG$res_col)
  if (is.null(R) || !ncol(R$mat)) { message("No matrix for ", cur); next }
  comm <- R$mat; meta <- R$meta
  cat("Rows: ", nrow(comm), " | prey: ", ncol(comm), " | testable cells: ", R$n_cells,
      " | rows early/late: ", sum(meta$period == CFG$lab_early), "/",
      sum(meta$period == CFG$lab_late), "\n", sep = "")
  note(cur, "rows", nrow(comm)); note(cur, "prey", ncol(comm)); note(cur, "cells", R$n_cells)

  # ---- 4.1 Test A --------------------------------------------------------------
  cat("\nTest A: adjusted MGLM (", CFG$perms, " within-cell permutations)...\n", sep = "")
  A <- run_isolated(core_mglm_adjusted,
                    list(comm = comm, cell = as.character(meta$cell), period = as.character(meta$period),
                         lev = levels(meta$period), n_perm = CFG$perms, seed = CFG$seed,
                         with_unadjusted = TRUE), "Test A (overall)")
  if (!is.null(A)) {
    cat(sprintf("  multivariate period effect: p = %s (adjusted) | p = %s (unadjusted, pooled)\n",
                format(A$mult_p, digits = 3), format(A$mult_p_unadj, digits = 3)))
    note(cur, "A_multivariate_p", A$mult_p); note(cur, "A_multivariate_p_unadjusted", A$mult_p_unadj)
    note(cur, "A_within_cell_permutations", A$restricted_permutations)
    if (!isTRUE(A$restricted_permutations))
      message("  NOTE: this mvabund has no bootID argument; Test A used unrestricted residual permutations.")
  }

  # ---- 4.2 Test B --------------------------------------------------------------
  cat("\nTest B: cell-level consistency (IndVal in ", R$n_cells, " cells)...\n", sep = "")
  cs <- cell_shares(comm, meta)
  iv_cells <- run_isolated(core_cell_indval,
                           list(comm = comm, cell = as.character(meta$cell), period = as.character(meta$period),
                                lev = levels(meta$period), n_perm = CFG$perms_cell, seed = CFG$seed,
                                min_rows = CFG$min_rows_per_period), "IndVal within cells")
  B <- cell_consistency(cs, iv_cells)

  cells_out <- merge(cs, if (!is.null(iv_cells)) as.data.table(iv_cells) else
    data.table(cell = character(), prey = character()), by = c("cell", "prey"), all.x = TRUE)
  cells_out <- merge(cells_out, unique(meta[, .(cell = as.character(cell), predator, size_cl, Area)]),
                     by = "cell", all.x = TRUE)
  write_csv(cells_out, file.path(DIR_DRIVERS, sprintf("drivers_cells_%s.csv", cur)))

  # ---- 4.3 The driver table -----------------------------------------------------
  drv <- assemble_drivers(A, B, comm, "overall")
  DRV[[cur]] <- drv
  write_csv(drv, file.path(DIR_DRIVERS, sprintf("drivers_%s.csv", cur)))
  note(cur, "n_both_criteria", sum(drv$n_criteria == 2)); note(cur, "n_one_criterion", sum(drv$n_criteria == 1))

  cat("\nDrivers (>= 1 criterion). A = adjusted MGLM, B = cell consistency:\n")
  print(as.data.frame(head(drv[n_criteria >= 1], 25)[, .(
    prey, direction, A = A_sig, B = B_sig, mglm_p = signif(mglm_p, 2),
    unadj_p = signif(mglm_p_unadjusted, 2),
    share_early = round(100 * share_early, 2), share_late = round(100 * share_late, 2),
    d_pct = round(100 * d_share, 2), cells = n_cells, gain = n_gain, loss = n_loss,
    modal = round(pct_modal), wilcox_bh = signif(wilcox_p_bh, 2))]), row.names = FALSE)

  # How much did adjusting for the cell change the answer? (the sampling-mix artefact)
  if (!is.null(A$uni_p_unadj)) {
    both <- drv[!is.na(mglm_p) & !is.na(mglm_p_unadjusted)]
    cat(sprintf("\nPrey significant unadjusted only (sampling-mix artefact candidates): %d | adjusted only: %d | both: %d\n",
                sum(both$mglm_p_unadjusted <= CFG$alpha & both$mglm_p > CFG$alpha),
                sum(both$mglm_p <= CFG$alpha & both$mglm_p_unadjusted > CFG$alpha),
                sum(both$mglm_p <= CFG$alpha & both$mglm_p_unadjusted <= CFG$alpha)))
  }

  # ---- 4.4 by axis (a reading of the cell table) --------------------------------
  cell_meta <- unique(meta[, .(cell = as.character(cell), predator, size_cl, Area)])
  cs2 <- merge(cs, cell_meta, by = "cell")
  by_axis <- rbindlist(lapply(list(c("Area", "Ecoregion"), c("size_cl", "Size class"),
                                   c("predator", "Predator")), function(ax) {
                                     cs2[, .(axis = ax[2], n_cells = .N, n_gain = sum(d_share > 0), n_loss = sum(d_share < 0),
                                             d_share = mean(d_share)), by = c("prey", ax[1])] %>%
                                       setnames(ax[1], "axis_level")
                                   }))
  by_axis[, direction := fifelse(d_share > 0, "gained", "lost")]
  write_csv(by_axis, file.path(DIR_DRIVERS, sprintf("drivers_by_axis_%s.csv", cur)))

  # ---- 4.5 figures: turnover and cell consistency -------------------------------
  top <- head(drv[n_criteria >= 1], CFG$top_n_figs)
  if (nrow(top)) {
    top[, prey_f := factor(prey, levels = rev(prey))]
    tl <- melt(top[, .(prey_f, share_early, share_late)], id.vars = "prey_f",
               variable.name = "period", value.name = "share")
    tl[, period := factor(fifelse(period == "share_early", CFG$lab_early, CFG$lab_late),
                          levels = c(CFG$lab_early, CFG$lab_late))]
    p_tv <- ggplot(top) +
      geom_segment(aes(x = 100 * share_early, xend = 100 * share_late, y = prey_f, yend = prey_f,
                       colour = direction), linewidth = 1.1, alpha = 0.6, show.legend = FALSE) +
      geom_point(data = tl, aes(x = 100 * share, y = prey_f, fill = period),
                 shape = 21, size = 2.8, colour = "grey20") +
      geom_text(aes(x = 100 * pmax(share_early, share_late), y = prey_f,
                    label = ifelse(n_criteria == 2, "A+B", ifelse(A_sig, "A", "B"))),
                hjust = -0.35, size = 2.8, colour = "grey30") +
      scale_colour_manual(values = DIR_COL) +
      scale_fill_manual(values = PER_COL, name = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
      labs(title = paste0("Prey associated with the inter-decade turnover, ", CURRENCY_LAB[[cur]], " currency"),
           subtitle = paste0("Mean share of the diet per cell (predator x size x ecoregion),\naveraged over ",
                             R$n_cells, " cells sampled in both decades.\nA = adjusted MGLM significant; ",
                             "B = consistent across cells\n(Wilcoxon BH p <= ", CFG$alpha, " and >= ",
                             100 * CFG$modal_share, " % of cells in the same direction)."),
           x = "Mean share of cell diet (%)", y = NULL) +
      theme_diag(base_size = 10) + theme(legend.position = "top")
    save_drv(p_tv, sprintf("Fig_prey_turnover_%s", cur), 7.8, 0.28 * nrow(top) + 2.8)

    cc <- top[, .(prey_f, gained = 100 * n_gain / n_cells, lost = -100 * n_loss / n_cells, n_cells)]
    ccl <- melt(cc, id.vars = c("prey_f", "n_cells"), variable.name = "dir", value.name = "pct")
    p_cc <- ggplot(ccl, aes(x = pct, y = prey_f, fill = dir)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = c(-100 * CFG$modal_share, 100 * CFG$modal_share),
                 linetype = "dashed", colour = "grey55", linewidth = 0.3) +
      geom_text(data = cc, aes(x = 104, y = prey_f, label = paste0("n=", n_cells)),
                inherit.aes = FALSE, size = 2.6, colour = "grey35", hjust = 0) +
      scale_fill_manual(values = DIR_COL, name = NULL) +
      scale_x_continuous(limits = c(-100, 125), breaks = seq(-100, 100, 50),
                         labels = function(x) abs(x)) +
      labs(title = paste0("Consistency across cells, ", CURRENCY_LAB[[cur]], " currency"),
           subtitle = paste0("Share of the cells in which each prey lost (left) or gained (right)\n",
                             "between decades. Dashed lines: the ", 100 * CFG$modal_share,
                             " % consistency threshold of criterion B."),
           x = "Cells (%)", y = NULL) +
      theme_diag(base_size = 10) + theme(legend.position = "top")
    save_drv(p_cc, sprintf("Fig_prey_cell_consistency_%s", cur), 7.2, 0.28 * nrow(top) + 2.6)
  }

  # ---- 4.6 resolution sweep ------------------------------------------------------
  cat("\nResolution sweep (", length(sweep_x), " resolutions; Test A with ",
      if (CFG$sweep_run_mglm) CFG$perms_sweep else "no", " permutations)...\n", sep = "")
  sweep <- rbindlist(lapply(sweep_x, function(x) {
    rc <- paste0("prey_category_", x, CFG$prey_family)
    Rx <- build_rows(DT, cur, rc)
    if (is.null(Rx) || ncol(Rx$mat) < 2) return(NULL)
    Ax <- if (isTRUE(CFG$sweep_run_mglm)) run_isolated(
      core_mglm_adjusted,
      list(comm = Rx$mat, cell = as.character(Rx$meta$cell), period = as.character(Rx$meta$period),
           lev = levels(Rx$meta$period), n_perm = CFG$perms_sweep, seed = CFG$seed,
           with_unadjusted = FALSE), paste0("Test A (x = ", x, ")")) else NULL
    Bx <- cell_consistency(cell_shares(Rx$mat, Rx$meta), NULL)
    dx <- assemble_drivers(Ax, Bx, Rx$mat, paste0("x = ", x))
    dx[, `:=`(x_threshold = x, n_categories = ncol(Rx$mat), n_cells = Rx$n_cells)]
    dx
  }), fill = TRUE)

  if (nrow(sweep)) {
    write_csv(sweep, file.path(DIR_DRIVERS, sprintf("drivers_resolution_sweep_%s.csv", cur)))
    n_res_run <- uniqueN(sweep$x_threshold)
    stability <- sweep[, .(
      n_res_present = uniqueN(x_threshold),
      n_res_sig     = sum(n_criteria >= 1),
      n_res_both    = sum(n_criteria == 2),
      direction_agree = uniqueN(direction[n_criteria >= 1]) <= 1,
      direction     = if (any(n_criteria >= 1)) names(which.max(table(direction[n_criteria >= 1]))) else NA_character_,
      d_share_mean  = mean(d_share), d_share_min = min(d_share), d_share_max = max(d_share)
    ), by = prey]
    stability[, pct_res_sig := round(100 * n_res_sig / n_res_present)]
    stability[, robust := n_res_present >= 3 & pct_res_sig >= 60 & direction_agree]
    stability[, n_res_run := n_res_run]
    stability[, abs_d := abs(d_share_mean)]; setorder(stability, -robust, -n_res_sig, -abs_d); stability[, abs_d := NULL]
    write_csv(stability, file.path(DIR_DRIVERS, sprintf("drivers_resolution_stability_%s.csv", cur)))
    cat("Resolutions completed: ", n_res_run, "/", length(sweep_x), " | prey robust across resolutions: ",
        sum(stability$robust), "\n", sep = "")
    print(as.data.frame(head(stability[robust == TRUE], 20)[, .(
      prey, direction, n_res_present, n_res_sig, n_res_both, d_pct = round(100 * d_share_mean, 2))]),
      row.names = FALSE)
    note(cur, "resolutions_completed", paste0(n_res_run, "/", length(sweep_x)))
    note(cur, "n_robust_across_resolutions", sum(stability$robust))

    show_prey <- union(head(drv[n_criteria >= 1], CFG$top_n_figs)$prey,
                       head(stability[robust == TRUE], CFG$top_n_figs)$prey)
    show_prey <- show_prey[!is.na(show_prey)]
    tiles <- sweep[prey %in% show_prey]
    if (nrow(tiles)) {
      tiles[, prey := factor(prey, levels = rev(show_prey))]
      tiles[, xf := factor(x_threshold, levels = sort(unique(sweep$x_threshold)))]
      tiles[, mark := fifelse(n_criteria == 2, "A+B", fifelse(A_sig, "A", fifelse(B_sig, "B", "")))]
      lim <- max(abs(100 * tiles$d_share))
      p_st <- ggplot(tiles, aes(x = xf, y = prey, fill = 100 * d_share)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = mark), size = 2.4, colour = "grey15") +
        scale_fill_gradient2(low = "#2c7bb6", mid = "grey96", high = "#d7191c",
                             limits = c(-lim, lim), name = "Change in\ncell-mean share\n(% of diet)") +
        scale_x_discrete(labels = function(x) paste0("x=", x)) +
        labs(title = paste0("Stability across taxonomic resolutions, ", CURRENCY_LAB[[cur]], " currency"),
             subtitle = paste0("Each column is the full cell-adjusted analysis at one grouping threshold (",
                               n_res_run, " resolutions). A / B / A+B = criteria met.\n",
                               "Blank cell = no category of that name at that resolution (merged)."),
             x = "Prey-grouping threshold", y = NULL) +
        theme_diag(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_drv(p_st, sprintf("Fig_resolution_stability_%s", cur), 9, 0.26 * nlevels(tiles$prey) + 2.8)
    }
  }

  # ---- 4.7 context: PERMANOVA on rows, PERMDISP ----------------------------------
  cat("\nPERMANOVA on rows (period after predator, size, ecoregion; permutations within cells)...\n")
  D <- vegdist(decostand(comm, "total"), method = "bray")
  pdf_ <- data.frame(predator = factor(meta$predator), size_class = factor(meta$size_cl),
                     Area = factor(meta$Area), period = meta$period, cell = meta$cell)
  ado <- try(adonis2(D ~ predator + size_class + Area + period, data = pdf_, by = "terms",
                     permutations = how(blocks = pdf_$cell, nperm = CFG$perms)), silent = TRUE)
  if (!inherits(ado, "try-error")) {
    print(ado)
    capture.output(ado, file = file.path(DIR_DRIVERS, sprintf("permanova_rows_%s.txt", cur)))
    note(cur, "permanova_period_R2", round(ado["period", "R2"], 4)); note(cur, "permanova_period_p", ado["period", "Pr(>F)"])
  } else message("  PERMANOVA failed: ", conditionMessage(attr(ado, "condition")))
  bd_a <- anova(betadisper(D, meta$period)); print(bd_a)
  capture.output(bd_a, file = file.path(DIR_DRIVERS, sprintf("permdisp_period_%s.txt", cur)))
  note(cur, "permdisp_p", bd_a[["Pr(>F)"]][1])

  # ---- 4.8 NMDS of cell centroids with decade arrows -----------------------------
  cat("\nNMDS (cell centroids)...\n")
  rel <- decostand(comm, "total")
  cdt <- as.data.table(rel); cdt[, `:=`(cell = as.character(meta$cell), period = as.character(meta$period))]
  cm  <- cdt[, lapply(.SD, mean), by = .(cell, period), .SDcols = colnames(rel)]
  cmeta <- cm[, .(cell, period = factor(period, levels = c(CFG$lab_early, CFG$lab_late)))]
  cmat  <- as.matrix(cm[, colnames(rel), with = FALSE]); rownames(cmat) <- paste(cm$cell, cm$period, sep = "@@")
  cmat  <- cmat[, colSums(cmat) > 0, drop = FALSE]

  if (nrow(cmat) >= 8) {
    nmc <- try(metaMDS(cmat, distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = 0), silent = TRUE)
    if (inherits(nmc, "try-error")) {
      message("  NMDS failed: ", conditionMessage(attr(nmc, "condition")))
    } else {
      check_stress(nmc, paste0("cell centroids, ", cur), nrow(cmat)); note(cur, "nmds_cells_stress", round(nmc$stress, 3))
      sc <- as.data.table(vegan::scores(nmc, display = "sites")); sc[, `:=`(cell = cmeta$cell, period = cmeta$period)]
      arr <- dcast(sc, cell ~ period, value.var = c("NMDS1", "NMDS2")); setnames(arr, c("cell", "x0", "x1", "y0", "y1"))
      p_cell <- ggplot(sc, aes(NMDS1, NMDS2, colour = period)) +
        geom_segment(data = arr, inherit.aes = FALSE, aes(x0, y0, xend = x1, yend = y1),
                     arrow = arrow(length = unit(0.12, "cm")), colour = "grey70", linewidth = 0.3, alpha = 0.6) +
        geom_point(size = 1.8, alpha = 0.85) +
        stat_ellipse(aes(group = period), type = "norm", level = 0.95, linewidth = 1)
      dv <- head(drv[n_criteria == 2], CFG$top_n_vectors)$prey; dv <- intersect(dv, colnames(cmat))
      if (length(dv) >= 2) {
        ef <- envfit(nmc, cmat[, dv, drop = FALSE], permutations = CFG$perms)
        ev <- as.data.frame(vegan::scores(ef, "vectors")); ev$prey <- dv
        ev$NMDS1 <- ev$NMDS1 * 0.8 * max(abs(sc$NMDS1)); ev$NMDS2 <- ev$NMDS2 * 0.8 * max(abs(sc$NMDS2))
        p_cell <- p_cell +
          geom_segment(data = ev, inherit.aes = FALSE, aes(0, 0, xend = NMDS1, yend = NMDS2),
                       arrow = arrow(length = unit(0.18, "cm")), colour = "grey25") +
          geom_text(data = ev, inherit.aes = FALSE, aes(NMDS1, NMDS2, label = prey), size = 2.5, colour = "grey15")
      }
      p_cell <- p_cell +
        scale_colour_manual(values = PER_COL, name = NULL) +
        labs(title = paste0("Cell centroids, ", CURRENCY_LAB[[cur]], " currency"),
             subtitle = sprintf(paste0("One point per predator x size class x ecoregion and decade; ",
                                       "Bray-Curtis, stress = %.2f.\nGrey arrows link the same cell across ",
                                       "decades; dark arrows are the prey meeting both criteria."), nmc$stress)) +
        theme_diag(base_size = 11) + theme(panel.grid = element_blank(), legend.position = "top")
      save_drv(p_cell, sprintf("Fig_nmds_cells_%s", cur), 7.4, 6.2)

      Dc <- vegdist(cmat, method = "bray")
      ado_c <- try(adonis2(Dc ~ period, data = data.frame(period = cmeta$period),
                           permutations = how(blocks = factor(cmeta$cell), nperm = CFG$perms)), silent = TRUE)
      if (!inherits(ado_c, "try-error")) {
        print(ado_c)
        capture.output(ado_c, file = file.path(DIR_DRIVERS, sprintf("permanova_cellpaired_%s.txt", cur)))
        note(cur, "permanova_cellpaired_R2", round(ado_c["period", "R2"], 3)); note(cur, "permanova_cellpaired_p", ado_c["period", "Pr(>F)"])
      }
    }
  }
}

# =============================================================================
# 5. CROSS-CURRENCY CONVERGENCE
# =============================================================================
cat("\n", strrep("=", 70), "\nCross-currency convergence\n", strrep("=", 70), "\n", sep = "")
if (!is.null(DRV$biomass) && !is.null(DRV$occurrence)) {
  pick <- function(d, sfx) {
    z <- d[, .(prey, mglm_p, wilcox_p_bh, d_share, direction, n_criteria, pct_modal)]
    old <- setdiff(names(z), "prey"); setnames(z, old, paste0(old, "_", sfx)); z
  }
  xc <- merge(pick(DRV$biomass, "B"), pick(DRV$occurrence, "O"), by = "prey", all = TRUE)
  xc[, `:=`(sig_B = !is.na(n_criteria_B) & n_criteria_B >= 1, sig_O = !is.na(n_criteria_O) & n_criteria_O >= 1)]
  xc[, class := fcase(
    sig_B & sig_O & direction_B == direction_O & n_criteria_B == 2 & n_criteria_O == 2, "A+B in both currencies",
    sig_B & sig_O & direction_B == direction_O, "robust (both currencies)",
    sig_B & sig_O, "discordant",
    sig_B | sig_O, "one currency",
    default = "not significant")]
  xc[, d_share_mean := rowMeans(cbind(d_share_B, d_share_O), na.rm = TRUE)]
  xc[, class := factor(class, levels = c("A+B in both currencies", "robust (both currencies)",
                                         "discordant", "one currency", "not significant"))]
  xc[, abs_d := abs(d_share_mean)]; setorder(xc, class, -abs_d); xc[, abs_d := NULL]
  write_csv(xc, file.path(DIR_DRIVERS, "drivers_cross_currency.csv"))
  cat("\nPrey by robustness class:\n"); print(table(xc$class))
  cat("\nRobust in both currencies:\n")
  print(as.data.frame(xc[class %in% levels(class)[1:2], .(prey, class, direction = direction_B,
                                                          d_pct_B = round(100 * d_share_B, 2),
                                                          d_pct_O = round(100 * d_share_O, 2))]), row.names = FALSE)
  SUMMARY$cross_currency <- as.list(table(xc$class))
} else message("Cross-currency join skipped: one currency has no driver table.")

# =============================================================================
# 6. RUN SUMMARY AND REPORTING GUIDANCE
# =============================================================================
con <- file(file.path(DIR_DRIVERS, "run_summary.txt"), "w")
writeLines(c("10_Prey_Drivers.R - run summary (cell-consistent design)", format(Sys.time()),
             paste0("Prey family: ", CFG$prey_family, " | manuscript resolution: ", CFG$res_col),
             paste0("Sweep resolutions: ", paste(sweep_x, collapse = ", ")),
             paste0("Permutations: ", CFG$perms, " (Test A) / ", CFG$perms_cell, " (IndVal in cells) / ",
                    CFG$perms_sweep, " (sweep)"),
             paste0("Testable cell: >= ", CFG$min_rows_per_period, " rows per period | criterion B: BH p <= ",
                    CFG$alpha, " and >= ", 100 * CFG$modal_share, " % of cells in the same direction"), ""), con)
for (nm in names(SUMMARY)) {
  writeLines(paste0("[", nm, "]"), con)
  for (k in names(SUMMARY[[nm]]))
    writeLines(paste0("  ", k, ": ", paste(unlist(SUMMARY[[nm]][[k]]), collapse = ", ")), con)
  writeLines("", con)
}
if (file.exists(FAIL_LOG)) writeLines(c("Failed / crashed tests:", paste0("  ", readLines(FAIL_LOG))), con) else
  writeLines("No test failed.", con)
close(con)

cat("\n", strrep("=", 70), "\nOutputs in ", normalizePath(DIR_DRIVERS), "\n\n", sep = "")
cat("Reporting guidance:\n")
cat("  * Association, not causation: 'prey associated with the transition'.\n")
cat("  * A prey earns a place in Results when it meets criteria A AND B at x = 440,\n")
cat("    is robust across the resolution sweep, and is 'robust' or 'A+B' in\n")
cat("    drivers_cross_currency.csv. One criterion only -> mention as tentative.\n")
cat("  * Compare mglm_p and mglm_p_unadjusted: prey significant only unadjusted are\n")
cat("    the sampling-mix artefacts the earlier pooled design would have reported.\n")
cat("  * A name present in one decade only (Gammaridae, Crustacea...) is an\n")
cat("    identification-level change: Methods caveat, not a Result.\n")
cat("  * Table 4 = drivers_<cur>.csv; Table S7 = drivers_by_axis; Fig = turnover +\n")
cat("    cell consistency; NMDS of cell centroids for the ordination.\n")
