# =============================================================================
# 10_Prey_Drivers.R - PREY ASSOCIATED WITH THE INTER-DECADE TURNOVER
#                     (cell-consistent design, self-contained engine)
# =============================================================================
# PURPOSE
#   Scripts 6b-6d establish THAT the diets reorganised between 2004-2006 and
#   2018-2019 (Substitution is the modal diagnostic at every scale, in both
#   currencies). This script asks WHICH prey categories carry that turnover -
#   which lost ground, which gained - and whether the answer survives the
#   choices that could manufacture it: the predator mix of the samples, body
#   size, ecoregion, taxonomic resolution, identification practice and diet
#   currency.
#
#   Everything here is ASSOCIATION, not mechanism. Nothing holds prey
#   availability or the environment constant, so a prey named below is a
#   candidate to be read against independent evidence on the Gulf prey field
#   (Savenkoff et al. 2007; Benoit & Swain 2008), not a demonstrated cause.
#
# WHY THE DESIGN IS CELL-BASED
#   Pooling every stomach of a trawl set into one "set diet" makes the result
#   depend on WHICH predators were sampled in the set, and the predator mix of
#   the survey changed between decades. The manuscript's framework compares
#   each predator x size x ecoregion CELL with itself for exactly that reason;
#   this script does the same.
#
# DESIGN
#   Row   = one predator species, one size class, in one trawl set: the prey of
#           its stomachs summed (biomass) or counted (occurrence), then turned
#           into a relative profile. The set remains the replicate (Hurlbert
#           1984); the predator is no longer mixed; the number of stomachs
#           behind a row cancels out in the profile.
#   Cell  = predator x size class x ecoregion, as in scripts 6b-6d. Only cells
#           sampled in BOTH periods with >= MIN_ROWS rows per period enter
#           (the framework's testable-cell rule).
#
#   TEST A - covariate-adjusted multivariate linear model on Hellinger
#            profiles (Legendre & Gallagher 2001): profile ~ cell + period, the
#            period term tested AFTER the cell term by permutation of the
#            reduced-model residuals RESTRICTED WITHIN CELLS (Freedman & Lane
#            1983; Anderson & Legendre 1999). Univariate p-values are adjusted
#            step-down across prey (maxT, Westfall & Young 1993). This is the
#            manylm test of mvabund (Warton et al. 2012; Wang et al. 2012),
#            used when the package is installed; otherwise the same test is
#            computed by the base-R engine below (identical resampling scheme,
#            sum-of-F multivariate statistic). The model without the cell term
#            is fitted for comparison (mglm_p_unadjusted), so the size of the
#            sampling-mix artefact is visible.
#
#   TEST B - cell-level consistency. Mean relative profile per cell and
#            decade, paired change per cell; across cells a Wilcoxon signed-rank
#            test per prey (Benjamini & Hochberg 1995 across prey), the share of
#            cells moving in the modal direction, and an Indicator Value
#            analysis within each cell (Dufrene & Legendre 1997; labdsv or the
#            base engine). Answers: is the change the same thing in most cells?
#
#   ROBUSTNESS LAYERS
#     harmonisation   CFG$harmonise merges labels that changed with
#                     identification practice (Gammaridae -> Amphipoda,
#                     Crustacea -> Arthropoda) BEFORE any test; the
#                     identification-reach tables document why.
#     by axis         Test B summarised within ecoregion / size class / predator.
#     by resolution   Tests A (fewer permutations) and B at two fixed-rank
#                     anchors (class, order) and >= 10 count-based resolutions
#                     (taxonomic sufficiency: Ellis 1985; Warwick 1988).
#     cross-currency  Biomass and occurrence joined (Hyslop 1980; Baker et al.
#                     2014): same direction, significant in both = strongest
#                     class of evidence produced here.
#
#   CONTEXT (secondary to Methods 2.4)
#     PERMANOVA on rows, D ~ predator + size_class + Area + period, period last,
#     permutations within cells (Anderson 2001; McArdle & Anderson 2001);
#     PERMDISP on period (Anderson 2006; Anderson & Walsh 2013); NMDS of cell
#     centroids with an arrow per cell across decades and a blocked PERMANOVA.
#
# ENGINE
#   CFG$engine = "auto": vegan / mvabund / labdsv when installed, base R
#   otherwise (MASS::isoMDS for the NMDS). The base engine implements the same
#   permutation tests and never segfaults; package fits run in a child R
#   process (callr) so that a crash cannot kill the session. run_summary.txt
#   records which engine produced each result.
#
# OUTPUTS  (Output_Drivers/)
#   drivers_<cur>.csv, drivers_cells_<cur>.csv, drivers_by_axis_<cur>.csv,
#   drivers_resolution_sweep_<cur>.csv, drivers_resolution_stability_<cur>.csv,
#   drivers_cross_currency.csv, identification_depth_by_period.csv,
#   identification_reach_by_period.csv, identification_reach_by_predator_period.csv,
#   Fig_prey_turnover_<cur>, Fig_prey_cell_consistency_<cur>,
#   Fig_resolution_stability_<cur>, Fig_nmds_cells_<cur>,
#   permanova_rows_<cur>.txt, permdisp_period_<cur>.txt,
#   permanova_cellpaired_<cur>.txt, driver_test_failures.txt, run_summary.txt
#
# REFERENCES
#   Anderson MJ (2001) Austral Ecol 26:32-46.  Anderson MJ (2006) Biometrics 62:245-253.
#   Anderson MJ, Legendre P (1999) J Stat Comput Simul 62:271-303.
#   Anderson MJ, Walsh DCI (2013) Ecol Monogr 83:557-574.
#   Baker R, Buckland A, Sheaves M (2014) Fish Fish 15:170-177.
#   Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
#   Benoit HP, Swain DP (2008) Can J Fish Aquat Sci 65:2088-2104.
#   Clarke KR (1993) Aust J Ecol 18:117-143.  Dufrene M, Legendre P (1997) Ecol Monogr 67:345-366.
#   Ellis D (1985) Mar Pollut Bull 16:459.  Freedman D, Lane D (1983) J Bus Econ Stat 1:292-298.
#   Hurlbert SH (1984) Ecol Monogr 54:187-211.  Hyslop EJ (1980) J Fish Biol 17:411-429.
#   Legendre P, Gallagher ED (2001) Oecologia 129:271-280.
#   McArdle BH, Anderson MJ (2001) Ecology 82:290-297.
#   Oksanen J et al. (2022) vegan.  Savenkoff C et al. (2007) Estuar Coast Shelf Sci 73:711-724.
#   Wang Y et al. (2012) Methods Ecol Evol 3:471-474.  Warton DI et al. (2012) Methods Ecol Evol 3:89-101.
#   Warwick RM (1988) Mar Pollut Bull 19:259-268.  Westfall PH, Young SS (1993) Wiley.
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr); library(readr)
})
HAVE <- list(vegan   = requireNamespace("vegan",   quietly = TRUE),
             mvabund = requireNamespace("mvabund", quietly = TRUE),
             labdsv  = requireNamespace("labdsv",  quietly = TRUE),
             MASS    = requireNamespace("MASS",    quietly = TRUE),
             callr   = requireNamespace("callr",   quietly = TRUE))

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
  # Fixed taxonomic ranks added to the sweep. A run at rank r is insensitive
  # to changes of identification depth BELOW r, but items stopping ABOVE r keep
  # a coarser label (finest coarser rank reached), so read the anchors with the
  # identification-reach tables. "class" (~88 % of items) is the safest anchor
  # but too coarse to name prey; "order" (~63 %) is the informative one.
  rank_levels    = c("class", "order"),
  rank_hierarchy = c("phylum", "subphylum", "class", "subclass", "order",
                     "infraorder", "family", "genus", "species"),

  yr_early  = 2004:2006,  yr_late = 2018:2019,
  lab_early = "2004-2006", lab_late = "2018-2019",

  currencies = c("biomass", "occurrence"),
  drop_prey  = c("Empty", "empty", "Unidentified", "unidentified",
                 "Digested", "digested", "NA"),

  # Taxonomic harmonisation BEFORE the tests (identification practice changed
  # between decades). Set list() for the raw picture; the merged run is the one
  # to report, the raw one goes to the Methods as the demonstration.
  harmonise = list(
    Amphipoda  = c("Gammaridae", "Amphipoda"),
    Arthropoda = c("Crustacea", "Arthropoda")
  ),

  # Cells and rows
  min_rows_per_period = 3,   # testable cell: >= 3 rows in EACH period (engine MIN_SETS)
  min_stom_per_row    = 1,   # rows built on fewer stomachs are dropped (1 = keep all)
  min_rows_per_prey   = 3,   # a prey must occur in at least this many rows
  min_cells_wilcox    = 5,

  # Resampling
  seed        = 42,
  perms       = 999,   # Test A at the manuscript resolution
  perms_cell  = 499,   # IndVal within cells
  perms_sweep = 199,   # Test A at the sweep resolutions
  perms_rows_permanova = 199,   # base engine only: row-level PERMANOVA is O(n^2) per permutation
  sweep_run_mglm = TRUE,
  alpha       = 0.05,
  modal_share = 0.60,

  top_n_figs = 20,
  top_n_vectors = 10,

  engine        = "auto",   # "auto" | "packages" | "base"
  isolate_tests = TRUE,     # package fits in a child R process (callr)
  test_timeout  = 3600
)

ENGINE <- if (CFG$engine == "base") "base" else if (CFG$engine == "packages") "packages" else
  if (HAVE$vegan && HAVE$mvabund && HAVE$labdsv) "packages" else "base"
if (CFG$engine == "packages" && !(HAVE$vegan && HAVE$mvabund && HAVE$labdsv))
  stop("CFG$engine = 'packages' but vegan, mvabund or labdsv is missing.")
cat("Engine:", ENGINE, if (ENGINE == "base") "(vegan/mvabund/labdsv not all installed - base-R implementations)", "\n")

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
SUMMARY <- list(engine = ENGINE)
note <- function(cur, key, value) SUMMARY[[cur]][[key]] <<- value

# =============================================================================
# 0. BASE-R ENGINE (used when the packages are absent; also the reference
#    implementation of every test in this script)
# =============================================================================
BASE <- list()
BASE$rel_profile <- function(m) { m <- as.matrix(m); rs <- rowSums(m); rs[rs == 0] <- 1; m / rs }
BASE$hellinger   <- function(m) sqrt(BASE$rel_profile(m))
BASE$bray_dist   <- function(m) {           # Bray-Curtis on rows (Bray & Curtis 1957)
  m <- as.matrix(m); n <- nrow(m); rs <- rowSums(m)
  D <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    j <- (i + 1):n
    num <- colSums(abs(t(m[j, , drop = FALSE]) - m[i, ]))
    D[i, j] <- D[j, i] <- num / (rs[i] + rs[j])
  }
  D[!is.finite(D)] <- 0
  stats::as.dist(D)
}
BASE$block_perm_ids <- function(block, n_perm, seed) {   # 1-based, one row per resample
  set.seed(seed)
  idx <- seq_along(block); spl <- split(idx, block)
  t(vapply(seq_len(n_perm), function(i) {
    out <- idx
    for (s in spl) if (length(s) > 1) out[s] <- s[sample.int(length(s))]
    out }, integer(length(block))))
}
# Test A, base implementation. Y = Hellinger matrix (rows x prey).
# Reduced model ~ cell (or ~ 1 if not adjusted); full model adds period.
# Freedman-Lane residual permutation within blocks; sum-of-F multivariate
# statistic; step-down maxT adjusted univariate p (Westfall & Young 1993).
BASE$mglm <- function(Y, cell, period, lev, n_perm, seed, adjusted = TRUE) {
  Y <- as.matrix(Y); n <- nrow(Y); p <- ncol(Y)
  period <- factor(as.character(period), levels = lev); cell <- factor(cell)
  X_red  <- if (adjusted) stats::model.matrix(~ cell) else matrix(1, n, 1)
  X_full <- cbind(X_red, period = as.numeric(period == lev[2]))
  q_red <- qr(X_red); q_full <- qr(X_full)
  df1 <- q_full$rank - q_red$rank; df2 <- n - q_full$rank
  if (df1 < 1 || df2 < 2) stop("Test A: period is confounded with the cells or too few rows.")
  fit_F <- function(Yb) {
    rss_red  <- colSums(qr.resid(q_red,  Yb)^2)
    rss_full <- colSums(qr.resid(q_full, Yb)^2)
    Fj <- ((rss_red - rss_full) / df1) / (rss_full / df2)
    Fj[!is.finite(Fj)] <- 0; Fj
  }
  F_obs <- fit_F(Y)
  fitted_red <- qr.fitted(q_red, Y); E_red <- Y - fitted_red
  ids <- BASE$block_perm_ids(if (adjusted) cell else rep(1L, n), n_perm, seed)
  F_perm <- matrix(NA_real_, n_perm, p)
  for (b in seq_len(n_perm)) F_perm[b, ] <- fit_F(fitted_red + E_red[ids[b, ], , drop = FALSE])
  S_obs <- sum(F_obs); p_mult <- (1 + sum(rowSums(F_perm) >= S_obs)) / (n_perm + 1)
  ord <- order(F_obs, decreasing = TRUE); p_adj <- numeric(p)
  for (k in seq_len(p)) {
    idx  <- ord[k:p]
    maxT <- if (length(idx) == 1) F_perm[, idx] else apply(F_perm[, idx, drop = FALSE], 1, max)
    p_adj[ord[k]] <- (1 + sum(maxT >= F_obs[ord[k]])) / (n_perm + 1)
  }
  p_adj[ord] <- cummax(p_adj[ord])
  coef_period <- qr.coef(q_full, Y)[ncol(X_full), ]
  list(uni_p = setNames(p_adj, colnames(Y)), mult_p = p_mult,
       coef_period = setNames(as.numeric(coef_period), colnames(Y)),
       F_obs = setNames(F_obs, colnames(Y)), restricted_permutations = adjusted, engine = "base")
}
# IndVal (Dufrene & Legendre 1997): specificity x fidelity, max over groups,
# p by permutation of group labels.
BASE$indval <- function(m, g, n_perm, seed) {
  set.seed(seed); m <- as.matrix(m); g <- as.integer(g); K <- sort(unique(g))
  stat <- function(gg) {
    A <- sapply(K, function(k) colMeans(m[gg == k, , drop = FALSE]))
    A <- A / pmax(rowSums(A), 1e-12)
    B <- sapply(K, function(k) colMeans(m[gg == k, , drop = FALSE] > 0))
    IV <- A * B
    list(indcls = apply(IV, 1, max), maxcls = apply(IV, 1, which.max))
  }
  obs <- stat(g); cnt <- integer(ncol(m))
  for (b in seq_len(n_perm)) cnt <- cnt + (stat(sample(g))$indcls >= obs$indcls)
  list(indcls = as.numeric(obs$indcls), maxcls = as.integer(obs$maxcls),
       pval = (cnt + 1) / (n_perm + 1), names = colnames(m))
}
# PERMANOVA with sequential terms (McArdle & Anderson 2001). p of the LAST
# term by permutation restricted within `blocks`; earlier terms free.
BASE$permanova <- function(D, terms, blocks = NULL, n_perm, seed) {
  Dm <- as.matrix(D); n <- nrow(Dm)
  A <- -0.5 * Dm^2
  G <- A - outer(rowMeans(A), rep(1, n)) - outer(rep(1, n), colMeans(A)) + mean(A)
  tot <- sum(diag(G))
  Qs <- list(); ranks <- integer(0); X <- matrix(1, n, 1)
  for (t in names(terms)) {
    X <- cbind(X, stats::model.matrix(~ f, data.frame(f = factor(terms[[t]])))[, -1, drop = FALSE])
    q <- qr(X); Qs[[t]] <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]; ranks[t] <- q$rank
  }
  ss_model <- function(Q, perm = NULL) {
    Qp <- if (is.null(perm)) Q else Q[order(perm), , drop = FALSE]
    sum(Qp * (G %*% Qp))
  }
  K <- length(terms); SS <- numeric(K); df <- numeric(K); prev <- 0; prev_rank <- 1
  for (k in seq_len(K)) { s <- ss_model(Qs[[k]]); SS[k] <- s - prev; df[k] <- ranks[k] - prev_rank; prev <- s; prev_rank <- ranks[k] }
  SS_res <- tot - prev; df_res <- n - ranks[K]
  Fv <- (SS / df) / (SS_res / df_res)
  pv <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    ids <- BASE$block_perm_ids(if (k == K && !is.null(blocks)) blocks else rep(1L, n), n_perm, seed + k)
    cnt <- 0
    for (b in seq_len(n_perm)) {
      pm <- ids[b, ]
      s_k  <- ss_model(Qs[[k]], pm); s_km <- if (k > 1) ss_model(Qs[[k - 1]], pm) else 0
      s_K  <- if (k == K) s_k else ss_model(Qs[[K]], pm)
      Fb <- ((s_k - s_km) / df[k]) / ((tot - s_K) / df_res)
      if (Fb >= Fv[k]) cnt <- cnt + 1
    }
    pv[k] <- (cnt + 1) / (n_perm + 1)
  }
  out <- data.frame(Df = c(df, df_res, n - 1), SumOfSqs = c(SS, SS_res, tot),
                    R2 = c(SS, SS_res, tot) / tot, F = c(Fv, NA, NA), `Pr(>F)` = c(pv, NA, NA),
                    check.names = FALSE)
  rownames(out) <- c(names(terms), "Residual", "Total")
  out
}
# PERMDISP (Anderson 2006): distance of each object to its group centroid,
# computed from the distance matrix, then a one-way test with label permutation.
BASE$dispersion <- function(D, g, n_perm, seed) {
  Dm <- as.matrix(D)^2; g <- factor(g); n <- nrow(Dm); z <- numeric(n)
  for (k in levels(g)) {
    i <- which(g == k); ng <- length(i)
    within <- sum(Dm[i, i]) / (2 * ng^2)
    z[i] <- sqrt(pmax(rowMeans(Dm[i, i, drop = FALSE]) - within, 0))
  }
  f_stat <- function(z, g) { m <- tapply(z, g, mean); nk <- table(g); K <- nlevels(g)
  ssb <- sum(nk * (m - mean(z))^2); ssw <- sum((z - m[as.character(g)])^2)
  (ssb / (K - 1)) / (ssw / (length(z) - K)) }
  F_obs <- f_stat(z, g); set.seed(seed); cnt <- 0
  for (b in seq_len(n_perm)) if (f_stat(z, sample(g)) >= F_obs) cnt <- cnt + 1
  out <- data.frame(Df = c(nlevels(g) - 1, n - nlevels(g)), F = c(F_obs, NA),
                    `Pr(>F)` = c((cnt + 1) / (n_perm + 1), NA), check.names = FALSE)
  rownames(out) <- c("Groups", "Residuals"); attr(out, "mean_dist") <- tapply(z, g, mean); out
}
BASE$nmds <- function(D, k = 2) {
  Dm <- as.matrix(D); Dm[Dm <= 0] <- 1e-6; diag(Dm) <- 0
  fit <- MASS::isoMDS(stats::as.dist(Dm), k = k, trace = FALSE)
  list(points = fit$points, stress = fit$stress / 100, engine = "base")
}

# ---- wrappers: packages when available, base otherwise --------------------
dist_bray <- function(m) if (ENGINE == "packages") vegan::vegdist(BASE$rel_profile(m), "bray") else BASE$bray_dist(BASE$rel_profile(m))

run_permanova <- function(D, terms, blocks, n_perm, label) {
  tryCatch({
    if (ENGINE == "packages") {
      df <- as.data.frame(lapply(terms, factor)); names(df) <- names(terms)
      f  <- stats::as.formula(paste("D ~", paste(names(terms), collapse = " + ")))
      vegan::adonis2(f, data = df, by = "terms",
                     permutations = if (is.null(blocks)) n_perm else
                       permute::how(blocks = factor(blocks), nperm = n_perm))
    } else BASE$permanova(D, terms, blocks, n_perm, CFG$seed)
  }, error = function(e) { log_failure(label, conditionMessage(e)); NULL })
}
run_dispersion <- function(D, g, n_perm, label) {
  tryCatch({
    if (ENGINE == "packages") stats::anova(vegan::betadisper(D, g)) else BASE$dispersion(D, g, n_perm, CFG$seed)
  }, error = function(e) { log_failure(label, conditionMessage(e)); NULL })
}
run_nmds <- function(m, label) {
  tryCatch({
    if (ENGINE == "packages") {
      nm <- vegan::metaMDS(BASE$rel_profile(m), distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = 0)
      list(points = vegan::scores(nm, display = "sites"), stress = nm$stress, engine = "vegan")
    } else BASE$nmds(dist_bray(m))
  }, error = function(e) { log_failure(label, conditionMessage(e)); NULL })
}
# Prey vectors on the ordination: correlation of each prey with the two axes.
fit_vectors <- function(nm_points, mat, dv) {
  cors <- t(sapply(dv, function(px) c(stats::cor(mat[, px], nm_points[, 1]), stats::cor(mat[, px], nm_points[, 2]))))
  cors[!is.finite(cors)] <- 0
  data.frame(prey = dv, NMDS1 = cors[, 1], NMDS2 = cors[, 2])
}

# =============================================================================
# 1. CHILD-PROCESS CORES (package fits) and the dispatchers
# =============================================================================
core_mglm_pkg <- function(comm, cell, period, lev, n_perm, seed, with_unadjusted) {
  suppressPackageStartupMessages({ library(vegan); library(mvabund) })
  cell <- factor(cell); period <- factor(as.character(period), levels = lev)
  mv <- mvabund::mvabund(as.matrix(vegan::decostand(comm, method = "hellinger")))
  set.seed(seed)
  idx <- seq_along(cell); spl <- split(idx, cell)
  ids <- t(vapply(seq_len(n_perm), function(i) { out <- idx
  for (s in spl) if (length(s) > 1) out[s] <- s[sample.int(length(s))]; out }, integer(length(cell))))
  fit <- mvabund::manylm(mv ~ cell + period); restricted <- TRUE
  an <- tryCatch(anova(fit, resamp = "perm.resid", test = "F", p.uni = "adjusted", bootID = ids),
                 error = function(e) { restricted <<- FALSE
                 anova(fit, resamp = "perm.resid", test = "F", p.uni = "adjusted", nBoot = n_perm) })
  uni_p <- an$uni.p["period", ]; uni_p <- uni_p[is.finite(uni_p)]
  mult_p <- tryCatch(an$table["period", "Pr(>F)"], error = function(e) NA_real_)
  cf <- tryCatch(coef(fit)[grep("^period", rownames(coef(fit)))[1], ], error = function(e) rep(NA_real_, ncol(comm)))
  out <- list(uni_p = uni_p, mult_p = mult_p, coef_period = cf, restricted_permutations = restricted, engine = "mvabund")
  if (isTRUE(with_unadjusted)) {
    set.seed(seed); fit0 <- mvabund::manylm(mv ~ period)
    an0 <- anova(fit0, resamp = "perm.resid", test = "F", p.uni = "adjusted", nBoot = n_perm)
    u0 <- an0$uni.p["period", ]; out$uni_p_unadj <- u0[is.finite(u0)]
    out$mult_p_unadj <- tryCatch(an0$table["period", "Pr(>F)"], error = function(e) NA_real_)
  }
  out
}
core_indval_pkg <- function(comm, cell, period, lev, n_perm, seed, min_rows) {
  suppressPackageStartupMessages({ library(vegan); library(labdsv) })
  set.seed(seed); cell <- as.character(cell); period <- factor(as.character(period), levels = lev)
  rel <- vegan::decostand(comm, "total"); out <- list()
  for (cl in unique(cell)) {
    i <- which(cell == cl); g <- droplevels(period[i])
    if (nlevels(g) != 2 || any(table(g) < min_rows)) next
    m <- rel[i, , drop = FALSE]; m <- m[, colSums(m > 0) >= 2, drop = FALSE]
    if (ncol(m) < 2) next
    iv <- tryCatch(labdsv::indval(as.data.frame(m), as.integer(g), numitr = n_perm), error = function(e) NULL)
    if (is.null(iv)) next
    out[[cl]] <- data.frame(cell = cl, prey = names(iv$indcls), indval = as.numeric(iv$indcls),
                            indval_p = as.numeric(iv$pval), indval_grp = levels(g)[iv$maxcls], stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

run_test_A <- function(comm, meta, n_perm, with_unadjusted, label) {
  t0 <- Sys.time()
  res <- tryCatch({
    if (ENGINE == "packages") {
      args <- list(comm = comm, cell = as.character(meta$cell), period = as.character(meta$period),
                   lev = levels(meta$period), n_perm = n_perm, seed = CFG$seed, with_unadjusted = with_unadjusted)
      if (isTRUE(CFG$isolate_tests) && HAVE$callr)
        callr::r(function(core, args) do.call(core, args), args = list(core = core_mglm_pkg, args = args),
                 timeout = CFG$test_timeout, show = FALSE)
      else do.call(core_mglm_pkg, args)
    } else {
      Y <- BASE$hellinger(comm)
      out <- BASE$mglm(Y, meta$cell, meta$period, levels(meta$period), n_perm, CFG$seed, adjusted = TRUE)
      if (isTRUE(with_unadjusted)) {
        u <- BASE$mglm(Y, meta$cell, meta$period, levels(meta$period), n_perm, CFG$seed, adjusted = FALSE)
        out$uni_p_unadj <- u$uni_p; out$mult_p_unadj <- u$mult_p
      }
      out
    }
  }, error = function(e) { log_failure(label, conditionMessage(e)); NULL })
  cat(sprintf("    %s: %.1f min\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}

run_cell_indval <- function(comm, meta, n_perm, label) {
  t0 <- Sys.time()
  res <- tryCatch({
    if (ENGINE == "packages") {
      args <- list(comm = comm, cell = as.character(meta$cell), period = as.character(meta$period),
                   lev = levels(meta$period), n_perm = n_perm, seed = CFG$seed, min_rows = CFG$min_rows_per_period)
      if (isTRUE(CFG$isolate_tests) && HAVE$callr)
        callr::r(function(core, args) do.call(core, args), args = list(core = core_indval_pkg, args = args),
                 timeout = CFG$test_timeout, show = FALSE)
      else do.call(core_indval_pkg, args)
    } else {
      rel <- BASE$rel_profile(comm); cell <- as.character(meta$cell); period <- meta$period; out <- list()
      for (cl in unique(cell)) {
        i <- which(cell == cl); g <- droplevels(period[i])
        if (nlevels(g) != 2 || any(table(g) < CFG$min_rows_per_period)) next
        m <- rel[i, , drop = FALSE]; m <- m[, colSums(m > 0) >= 2, drop = FALSE]
        if (ncol(m) < 2) next
        iv <- BASE$indval(m, as.integer(g), n_perm, CFG$seed)
        out[[cl]] <- data.frame(cell = cl, prey = iv$names, indval = iv$indcls, indval_p = iv$pval,
                                indval_grp = levels(g)[iv$maxcls], stringsAsFactors = FALSE)
      }
      do.call(rbind, out)
    }
  }, error = function(e) { log_failure(label, conditionMessage(e)); NULL })
  cat(sprintf("    %s: %.1f min\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}

# =============================================================================
# 2. LOAD AND PREPARE
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
n0 <- nrow(DT)
DT <- DT[!is.na(period) & !is.na(get(CFG$col_pred)) & !is.na(get(CFG$col_size)) & !is.na(get(CFG$col_area))]
if (nrow(DT) < n0) message(n0 - nrow(DT), " prey records dropped (outside the two periods or missing predator/size/Area).")
DT[, period := factor(period, levels = c(CFG$lab_early, CFG$lab_late))]
DT[, set_uid := paste(get(CFG$col_year), get(CFG$col_vessel), get(CFG$col_set), sep = "_")]
DT[, `:=`(predator = as.character(get(CFG$col_pred)), size_cl = as.character(get(CFG$col_size)),
          Area = as.character(get(CFG$col_area)))]
DT[, cell := paste(predator, size_cl, Area, sep = "|")]
DT[, row_id := paste(set_uid, predator, size_cl, sep = "#")]

# Resolutions available and the sweep
avail_x <- sort(as.integer(sub(paste0("^prey_category_(\\d+)", CFG$prey_family, "$"), "\\1",
                               grep(paste0("^prey_category_\\d+", CFG$prey_family, "$"), names(DT), value = TRUE))))
if (!CFG$res_col %in% names(DT)) stop("CFG$res_col '", CFG$res_col, "' not in the data.")
sweep_x <- intersect(CFG$sweep_x, avail_x)
if (length(sweep_x) < CFG$sweep_min && length(avail_x) >= CFG$sweep_min)
  sweep_x <- avail_x[unique(round(seq(1, length(avail_x), length.out = CFG$sweep_min)))]
sweep_x <- sort(union(sweep_x, 440L))

rank_cols <- character(0)
if (all(c("phylum", "class") %in% names(DT))) {
  for (r in CFG$rank_levels) {
    hier <- rev(intersect(CFG$rank_hierarchy[seq_len(match(r, CFG$rank_hierarchy))], names(DT)))
    nm <- paste0("prey_rank_", r); DT[, (nm) := NA_character_]
    for (h in hier) DT[is.na(get(nm)) & !is.na(get(h)) & get(h) != "", (nm) := as.character(get(h))]
    rank_cols <- c(rank_cols, nm)
  }
} else message("Rank columns (phylum, class, ...) not in dat_classed: rank-level anchors skipped.")

sweep_specs <- rbind(
  data.table(col = rank_cols, label = paste0("rank: ", CFG$rank_levels[seq_along(rank_cols)]),
             x_order = -rev(seq_along(rank_cols))),
  data.table(col = paste0("prey_category_", sweep_x, CFG$prey_family), label = paste0("x=", sweep_x), x_order = sweep_x))
cat("Resolution sweep: ", paste(sweep_specs$label, collapse = ", "), " (", nrow(sweep_specs), ")\n", sep = "")
cat("Prey records: ", nrow(DT), " | sets: ", uniqueN(DT$set_uid), " | rows (set x predator x size): ",
    uniqueN(DT$row_id), " | cells: ", uniqueN(DT$cell), "\n", sep = "")

# Identification depth / reach by period: the evidence behind the harmonisation.
depth_col <- intersect(c("tax_level_lowest", "tax_level"), names(DT))[1]
if (!is.na(depth_col)) {
  rank_order_full <- CFG$rank_hierarchy
  if (file.exists("R_helpers/taxonomic_rank_order.R")) {
    e2 <- new.env(); sys.source("R_helpers/taxonomic_rank_order.R", envir = e2)
    if (exists("ranknfile", envir = e2)) rank_order_full <- get("ranknfile", envir = e2)
  }
  ridx <- function(x) match(tolower(as.character(x)), rank_order_full)
  i_class <- match("class", rank_order_full); i_order <- match("order", rank_order_full)
  i_family <- match("family", rank_order_full); i_genus <- match("genus", rank_order_full)
  DT[, .depth := as.character(get(depth_col))]
  depth <- DT[!is.na(.depth), .N, by = .(period, tax_level = .depth)]
  depth[, pct := round(100 * N / sum(N), 1), by = period][, rank_order := ridx(tax_level)]
  setorder(depth, period, rank_order)
  write_csv(dcast(depth, tax_level + rank_order ~ period, value.var = "pct")[order(rank_order)],
            file.path(DIR_DRIVERS, "identification_depth_by_period.csv"))
  reach <- DT[!is.na(.depth), { r <- ridx(.depth)
  .(pct_class_or_finer = round(100 * mean(r >= i_class, na.rm = TRUE), 1),
    pct_order_or_finer = round(100 * mean(r >= i_order, na.rm = TRUE), 1),
    pct_family_or_finer = round(100 * mean(r >= i_family, na.rm = TRUE), 1),
    pct_genus_or_finer = round(100 * mean(r >= i_genus, na.rm = TRUE), 1), n_items = .N) }, by = period]
  reach_pred <- DT[!is.na(.depth), { r <- ridx(.depth)
  .(pct_order_or_finer = round(100 * mean(r >= i_order, na.rm = TRUE), 1),
    pct_family_or_finer = round(100 * mean(r >= i_family, na.rm = TRUE), 1), n_items = .N) }, by = .(predator, period)]
  write_csv(reach, file.path(DIR_DRIVERS, "identification_reach_by_period.csv"))
  write_csv(reach_pred, file.path(DIR_DRIVERS, "identification_reach_by_predator_period.csv"))
  DT[, .depth := NULL]
  cat("\nIdentification reach (% of prey items identified at least to the rank), by period:\n")
  print(as.data.frame(reach), row.names = FALSE)
  SUMMARY$identification_reach <- as.list(setNames(
    paste0("class ", reach$pct_class_or_finer, " | order ", reach$pct_order_or_finer,
           " | family ", reach$pct_family_or_finer, " | genus ", reach$pct_genus_or_finer), as.character(reach$period)))
}

# -----------------------------------------------------------------------------
# Row matrix builder: rows = set x predator x size class, columns = prey.
# -----------------------------------------------------------------------------
build_rows <- function(dat, currency, res_col) {
  d <- dat[, c("row_id", "cell", "set_uid", "predator", "size_cl", "Area", "period",
               CFG$col_stom, CFG$col_wt, res_col), with = FALSE]
  d[, prey := as.character(get(res_col))]
  d <- d[!is.na(prey) & !(prey %chin% CFG$drop_prey)]
  if (length(CFG$harmonise)) {
    map <- unlist(lapply(names(CFG$harmonise), function(to)
      setNames(rep(to, length(CFG$harmonise[[to]])), CFG$harmonise[[to]])))
    d[prey %chin% names(map), prey := map[prey]]
  }
  if (!nrow(d)) return(NULL)
  agg <- if (currency == "biomass") d[, .(val = sum(get(CFG$col_wt), na.rm = TRUE)), by = .(row_id, prey)] else
    d[, .(val = as.numeric(uniqueN(get(CFG$col_stom)))), by = .(row_id, prey)]
  m   <- dcast(agg, row_id ~ prey, value.var = "val", fill = 0)
  mat <- as.matrix(m[, -1]); rownames(mat) <- m$row_id
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  meta <- unique(d[, .(row_id, cell, set_uid, predator, size_cl, Area, period)])
  meta <- meta[match(rownames(mat), row_id)]
  meta[, n_stom := d[, uniqueN(get(CFG$col_stom)), by = row_id][match(meta$row_id, row_id), V1]]
  keep0 <- meta$n_stom >= CFG$min_stom_per_row
  mat <- mat[keep0, , drop = FALSE]; meta <- meta[keep0]
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

# =============================================================================
# 3. TEST B (cell-level consistency) and the driver table
# =============================================================================
cell_shares <- function(comm, meta) {
  rel <- BASE$rel_profile(comm)
  dt  <- as.data.table(rel); dt[, `:=`(cell = as.character(meta$cell), period = as.character(meta$period))]
  long <- melt(dt, id.vars = c("cell", "period"), variable.name = "prey", value.name = "share", variable.factor = FALSE)
  cp <- long[, .(share = mean(share), n_rows = .N), by = .(cell, period, prey)]
  w  <- dcast(cp, cell + prey ~ period, value.var = c("share", "n_rows"))
  setnames(w, paste0("share_", CFG$lab_early), "share_early"); setnames(w, paste0("share_", CFG$lab_late), "share_late")
  setnames(w, paste0("n_rows_", CFG$lab_early), "n_rows_early"); setnames(w, paste0("n_rows_", CFG$lab_late), "n_rows_late")
  w[, d_share := share_late - share_early]; w[]
}

cell_consistency <- function(cs, iv_cells, alpha = CFG$alpha) {
  x <- cs[, .(
    n_cells = .N, n_gain = sum(d_share > 0), n_loss = sum(d_share < 0), n_flat = sum(d_share == 0),
    share_early = mean(share_early), share_late = mean(share_late), d_share = mean(d_share),
    d_share_median = median(d_share),
    share_early_w = weighted.mean(share_early, pmin(n_rows_early, n_rows_late)),
    share_late_w  = weighted.mean(share_late,  pmin(n_rows_early, n_rows_late)),
    d_share_w     = weighted.mean(d_share,     pmin(n_rows_early, n_rows_late)),
    wilcox_p = if (.N >= CFG$min_cells_wilcox && any(d_share != 0))
      suppressWarnings(stats::wilcox.test(d_share, mu = 0, exact = FALSE)$p.value) else NA_real_
  ), by = prey]
  x[, direction := fifelse(d_share > 0, "gained", "lost")]
  x[, pct_modal := 100 * pmax(n_gain, n_loss) / n_cells]
  x[, wilcox_p_bh := p.adjust(wilcox_p, "BH")]
  if (!is.null(iv_cells) && nrow(iv_cells)) {
    ivs <- as.data.table(iv_cells)[indval_p <= alpha, .(iv_cells_early = sum(indval_grp == CFG$lab_early),
                                                        iv_cells_late = sum(indval_grp == CFG$lab_late)), by = prey]
    x <- merge(x, ivs, by = "prey", all.x = TRUE)
    x[is.na(iv_cells_early), iv_cells_early := 0L][is.na(iv_cells_late), iv_cells_late := 0L]
  } else x[, `:=`(iv_cells_early = NA_integer_, iv_cells_late = NA_integer_)]
  x[]
}

assemble_drivers <- function(A, B, label) {
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
  setcolorder(d, c("prey", "direction", "n_criteria", "A_sig", "B_sig", "mglm_p", "mglm_p_unadjusted",
                   "mglm_multivariate_p", "share_early", "share_late", "d_share", "d_share_median",
                   "share_early_w", "share_late_w", "d_share_w", "n_cells", "n_gain", "n_loss", "n_flat",
                   "pct_modal", "wilcox_p", "wilcox_p_bh", "iv_cells_early", "iv_cells_late"))
  d[]
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
      " | rows early/late: ", sum(meta$period == CFG$lab_early), "/", sum(meta$period == CFG$lab_late), "\n", sep = "")
  note(cur, "rows", nrow(comm)); note(cur, "prey", ncol(comm)); note(cur, "cells", R$n_cells)

  # ---- 4.1 Test A ------------------------------------------------------------
  cat("\nTest A: adjusted MGLM (", CFG$perms, " within-cell permutations)...\n", sep = "")
  A <- run_test_A(comm, meta, CFG$perms, TRUE, "Test A (overall)")
  if (!is.null(A)) {
    cat(sprintf("  multivariate period effect: p = %s (adjusted) | p = %s (unadjusted, pooled) | engine: %s\n",
                format(A$mult_p, digits = 3), format(A$mult_p_unadj, digits = 3), A$engine))
    note(cur, "A_engine", A$engine); note(cur, "A_within_cell_permutations", A$restricted_permutations)
    note(cur, "A_multivariate_p", A$mult_p); note(cur, "A_multivariate_p_unadjusted", A$mult_p_unadj)
  }

  # ---- 4.2 Test B ------------------------------------------------------------
  cat("\nTest B: cell-level consistency (IndVal in ", R$n_cells, " cells)...\n", sep = "")
  cs <- cell_shares(comm, meta)
  iv_cells <- run_cell_indval(comm, meta, CFG$perms_cell, "IndVal within cells")
  B <- cell_consistency(cs, iv_cells)
  cell_meta <- unique(meta[, .(cell = as.character(cell), predator, size_cl, Area)])
  cells_out <- merge(cs, if (!is.null(iv_cells)) as.data.table(iv_cells) else data.table(cell = character(), prey = character()),
                     by = c("cell", "prey"), all.x = TRUE)
  cells_out <- merge(cells_out, cell_meta, by = "cell", all.x = TRUE)
  write_csv(cells_out, file.path(DIR_DRIVERS, sprintf("drivers_cells_%s.csv", cur)))

  # ---- 4.3 the driver table -------------------------------------------------
  drv <- assemble_drivers(A, B, "overall"); DRV[[cur]] <- drv
  write_csv(drv, file.path(DIR_DRIVERS, sprintf("drivers_%s.csv", cur)))
  note(cur, "n_both_criteria", sum(drv$n_criteria == 2)); note(cur, "n_one_criterion", sum(drv$n_criteria == 1))
  cat("\nDrivers (>= 1 criterion). A = adjusted MGLM, B = cell consistency:\n")
  print(as.data.frame(head(drv[n_criteria >= 1], 25)[, .(
    prey, direction, A = A_sig, B = B_sig, mglm_p = signif(mglm_p, 2), unadj_p = signif(mglm_p_unadjusted, 2),
    share_early = round(100 * share_early, 2), share_late = round(100 * share_late, 2),
    d_pct = round(100 * d_share, 2), cells = n_cells, gain = n_gain, loss = n_loss,
    modal = round(pct_modal), wilcox_bh = signif(wilcox_p_bh, 2))]), row.names = FALSE)
  if (!is.null(A$uni_p_unadj)) {
    both <- drv[!is.na(mglm_p) & !is.na(mglm_p_unadjusted)]
    cat(sprintf("\nPrey significant unadjusted only (sampling-mix artefact candidates): %d | adjusted only: %d | both: %d\n",
                sum(both$mglm_p_unadjusted <= CFG$alpha & both$mglm_p > CFG$alpha),
                sum(both$mglm_p <= CFG$alpha & both$mglm_p_unadjusted > CFG$alpha),
                sum(both$mglm_p <= CFG$alpha & both$mglm_p_unadjusted <= CFG$alpha)))
  }

  # ---- 4.4 by axis --------------------------------------------------------------
  cs2 <- merge(cs, cell_meta, by = "cell")
  by_axis <- rbindlist(lapply(list(c("Area", "Ecoregion"), c("size_cl", "Size class"), c("predator", "Predator")), function(ax) {
    z <- cs2[, .(axis = ax[2], n_cells = .N, n_gain = sum(d_share > 0), n_loss = sum(d_share < 0), d_share = mean(d_share)),
             by = c("prey", ax[1])]
    setnames(z, ax[1], "axis_level"); z }))
  by_axis[, direction := fifelse(d_share > 0, "gained", "lost")]
  write_csv(by_axis, file.path(DIR_DRIVERS, sprintf("drivers_by_axis_%s.csv", cur)))

  # ---- 4.5 figures: turnover and cell consistency --------------------------------
  top <- head(drv[n_criteria >= 1], CFG$top_n_figs)
  if (nrow(top)) {
    top[, prey_f := factor(prey, levels = rev(prey))]
    tl <- melt(top[, .(prey_f, share_early, share_late)], id.vars = "prey_f", variable.name = "period", value.name = "share")
    tl[, period := factor(fifelse(period == "share_early", CFG$lab_early, CFG$lab_late), levels = c(CFG$lab_early, CFG$lab_late))]
    p_tv <- ggplot(top) +
      geom_segment(aes(x = 100 * share_early, xend = 100 * share_late, y = prey_f, yend = prey_f, colour = direction),
                   linewidth = 1.1, alpha = 0.6, show.legend = FALSE) +
      geom_point(data = tl, aes(x = 100 * share, y = prey_f, fill = period), shape = 21, size = 2.8, colour = "grey20") +
      geom_text(aes(x = 100 * pmax(share_early, share_late), y = prey_f,
                    label = ifelse(n_criteria == 2, "A+B", ifelse(A_sig, "A", "B"))), hjust = -0.35, size = 2.8, colour = "grey30") +
      scale_colour_manual(values = DIR_COL) + scale_fill_manual(values = PER_COL, name = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
      labs(title = paste0("Prey associated with the inter-decade turnover, ", CURRENCY_LAB[[cur]], " currency"),
           subtitle = paste0("Mean share of the diet per cell (predator x size x ecoregion),\naveraged over ", R$n_cells,
                             " cells sampled in both decades.\nA = adjusted MGLM significant; B = consistent across cells\n(Wilcoxon BH p <= ",
                             CFG$alpha, " and >= ", 100 * CFG$modal_share, " % of cells in the same direction)."),
           x = "Mean share of cell diet (%)", y = NULL) +
      theme_diag(base_size = 10) + theme(legend.position = "top")
    save_drv(p_tv, sprintf("Fig_prey_turnover_%s", cur), 7.8, 0.28 * nrow(top) + 2.8)

    cc <- top[, .(prey_f, gained = 100 * n_gain / n_cells, lost = -100 * n_loss / n_cells, n_cells)]
    ccl <- melt(cc, id.vars = c("prey_f", "n_cells"), variable.name = "dir", value.name = "pct")
    p_cc <- ggplot(ccl, aes(x = pct, y = prey_f, fill = dir)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = c(-100 * CFG$modal_share, 100 * CFG$modal_share), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
      geom_text(data = cc, aes(x = 104, y = prey_f, label = paste0("n=", n_cells)), inherit.aes = FALSE, size = 2.6, colour = "grey35", hjust = 0) +
      scale_fill_manual(values = DIR_COL, name = NULL) +
      scale_x_continuous(limits = c(-100, 125), breaks = seq(-100, 100, 50), labels = function(x) abs(x)) +
      labs(title = paste0("Consistency across cells, ", CURRENCY_LAB[[cur]], " currency"),
           subtitle = paste0("Share of the cells in which each prey lost (left) or gained (right)\nbetween decades. Dashed lines: the ",
                             100 * CFG$modal_share, " % consistency threshold of criterion B."),
           x = "Cells (%)", y = NULL) +
      theme_diag(base_size = 10) + theme(legend.position = "top")
    save_drv(p_cc, sprintf("Fig_prey_cell_consistency_%s", cur), 7.2, 0.28 * nrow(top) + 2.6)
  }

  # ---- 4.6 resolution sweep -----------------------------------------------------
  cat("\nResolution sweep (", nrow(sweep_specs), " resolutions; Test A with ",
      if (CFG$sweep_run_mglm) CFG$perms_sweep else "no", " permutations)...\n", sep = "")
  sweep <- rbindlist(lapply(seq_len(nrow(sweep_specs)), function(i) {
    rc <- sweep_specs$col[i]; lab <- sweep_specs$label[i]
    if (!rc %in% names(DT)) return(NULL)
    Rx <- build_rows(DT, cur, rc)
    if (is.null(Rx) || ncol(Rx$mat) < 2) return(NULL)
    Ax <- if (isTRUE(CFG$sweep_run_mglm)) run_test_A(Rx$mat, Rx$meta, CFG$perms_sweep, FALSE, paste0("Test A (", lab, ")")) else NULL
    Bx <- cell_consistency(cell_shares(Rx$mat, Rx$meta), NULL)
    dx <- assemble_drivers(Ax, Bx, lab)
    dx[, `:=`(x_threshold = lab, x_order = sweep_specs$x_order[i], n_categories = ncol(Rx$mat), n_cells = Rx$n_cells)]
    dx
  }), fill = TRUE)

  if (nrow(sweep)) {
    write_csv(sweep, file.path(DIR_DRIVERS, sprintf("drivers_resolution_sweep_%s.csv", cur)))
    n_res_run <- uniqueN(sweep$x_threshold)
    stability <- sweep[, .(
      n_res_present = uniqueN(x_threshold), n_res_sig = sum(n_criteria >= 1), n_res_both = sum(n_criteria == 2),
      direction_agree = uniqueN(direction[n_criteria >= 1]) <= 1,
      direction = if (any(n_criteria >= 1)) names(which.max(table(direction[n_criteria >= 1]))) else NA_character_,
      d_share_mean = mean(d_share), d_share_min = min(d_share), d_share_max = max(d_share)), by = prey]
    stability[, pct_res_sig := round(100 * n_res_sig / n_res_present)]
    stability[, robust := n_res_present >= 3 & pct_res_sig >= 60 & direction_agree]
    stability[, n_res_run := n_res_run]
    stability[, abs_d := abs(d_share_mean)]; setorder(stability, -robust, -n_res_sig, -abs_d); stability[, abs_d := NULL]
    write_csv(stability, file.path(DIR_DRIVERS, sprintf("drivers_resolution_stability_%s.csv", cur)))
    cat("Resolutions completed: ", n_res_run, "/", nrow(sweep_specs), " | prey robust across resolutions: ", sum(stability$robust), "\n", sep = "")
    print(as.data.frame(head(stability[robust == TRUE], 20)[, .(prey, direction, n_res_present, n_res_sig, n_res_both,
                                                                d_pct = round(100 * d_share_mean, 2))]), row.names = FALSE)
    note(cur, "resolutions_completed", paste0(n_res_run, "/", nrow(sweep_specs))); note(cur, "n_robust_across_resolutions", sum(stability$robust))

    show_prey <- union(head(drv[n_criteria >= 1], CFG$top_n_figs)$prey, head(stability[robust == TRUE], CFG$top_n_figs)$prey)
    show_prey <- show_prey[!is.na(show_prey)]
    tiles <- sweep[prey %in% show_prey]
    if (nrow(tiles)) {
      tiles[, prey := factor(prey, levels = rev(show_prey))]
      lv <- unique(sweep[order(x_order), x_threshold]); tiles[, xf := factor(x_threshold, levels = lv)]
      tiles[, mark := fifelse(n_criteria == 2, "A+B", fifelse(A_sig, "A", fifelse(B_sig, "B", "")))]
      lim <- max(abs(100 * tiles$d_share))
      p_st <- ggplot(tiles, aes(x = xf, y = prey, fill = 100 * d_share)) +
        geom_tile(colour = "white", linewidth = 0.4) + geom_text(aes(label = mark), size = 2.4, colour = "grey15") +
        scale_fill_gradient2(low = "#2c7bb6", mid = "grey96", high = "#d7191c", limits = c(-lim, lim), name = "Change in\ncell-mean share\n(% of diet)") +
        labs(title = paste0("Stability across taxonomic resolutions, ", CURRENCY_LAB[[cur]], " currency"),
             subtitle = paste0("Each column is the full cell-adjusted analysis at one grouping (", n_res_run,
                               " resolutions). A / B / A+B = criteria met.\nBlank cell = no category of that name at that resolution (merged)."),
             x = "Prey-grouping threshold", y = NULL) +
        theme_diag(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_drv(p_st, sprintf("Fig_resolution_stability_%s", cur), 9, 0.26 * nlevels(tiles$prey) + 2.8)
    }
  }

  # ---- 4.7 context: PERMANOVA on rows, PERMDISP -----------------------------------
  cat("\nPERMANOVA on rows (period after predator, size, ecoregion; permutations within cells)...\n")
  D <- dist_bray(comm)
  n_perm_rows <- if (ENGINE == "base") CFG$perms_rows_permanova else CFG$perms
  ado <- run_permanova(D, list(predator = meta$predator, size_class = meta$size_cl, Area = meta$Area, period = meta$period),
                       blocks = as.character(meta$cell), n_perm = n_perm_rows, label = paste("PERMANOVA rows", cur))
  if (!is.null(ado)) {
    print(ado); capture.output(print(ado), file = file.path(DIR_DRIVERS, sprintf("permanova_rows_%s.txt", cur)))
    note(cur, "permanova_period_R2", round(ado["period", "R2"], 4)); note(cur, "permanova_period_p", ado["period", "Pr(>F)"])
  }
  bd <- run_dispersion(D, meta$period, n_perm_rows, paste("PERMDISP", cur))
  if (!is.null(bd)) {
    print(bd); capture.output(print(bd), file = file.path(DIR_DRIVERS, sprintf("permdisp_period_%s.txt", cur)))
    note(cur, "permdisp_p", bd[["Pr(>F)"]][1])
  }

  # ---- 4.8 NMDS of cell centroids with decade arrows ------------------------------
  cat("\nNMDS (cell centroids)...\n")
  rel <- BASE$rel_profile(comm)
  cdt <- as.data.table(rel); cdt[, `:=`(cell = as.character(meta$cell), period = as.character(meta$period))]
  cm  <- cdt[, lapply(.SD, mean), by = .(cell, period), .SDcols = colnames(rel)]
  cmeta <- cm[, .(cell, period = factor(period, levels = c(CFG$lab_early, CFG$lab_late)))]
  cmat  <- as.matrix(cm[, colnames(rel), with = FALSE]); rownames(cmat) <- paste(cm$cell, cm$period, sep = "@@")
  cmat  <- cmat[, colSums(cmat) > 0, drop = FALSE]
  if (nrow(cmat) >= 8) {
    nmc <- run_nmds(cmat, paste("NMDS cells", cur))
    if (!is.null(nmc)) {
      s <- nmc$stress
      if (is.finite(s) && s > 0.2) warning("High NMDS stress for cell centroids, ", cur, ": ", round(s, 3), " (Clarke 1993).", call. = FALSE)
      note(cur, "nmds_cells_stress", round(s, 3)); note(cur, "nmds_engine", nmc$engine)
      sc <- as.data.table(nmc$points); setnames(sc, c("NMDS1", "NMDS2")); sc[, `:=`(cell = cmeta$cell, period = cmeta$period)]
      arr <- dcast(sc, cell ~ period, value.var = c("NMDS1", "NMDS2")); setnames(arr, c("cell", "x0", "x1", "y0", "y1"))
      p_cell <- ggplot(sc, aes(NMDS1, NMDS2, colour = period)) +
        geom_segment(data = arr, inherit.aes = FALSE, aes(x0, y0, xend = x1, yend = y1),
                     arrow = arrow(length = unit(0.12, "cm")), colour = "grey70", linewidth = 0.3, alpha = 0.6) +
        geom_point(size = 1.8, alpha = 0.85) +
        stat_ellipse(aes(group = period), type = "norm", level = 0.95, linewidth = 1)
      dv <- head(drv[n_criteria == 2], CFG$top_n_vectors)$prey; dv <- intersect(dv, colnames(cmat))
      if (length(dv) >= 2) {
        ev <- fit_vectors(as.matrix(sc[, .(NMDS1, NMDS2)]), cmat, dv)
        ev$NMDS1 <- ev$NMDS1 * 0.8 * max(abs(sc$NMDS1)); ev$NMDS2 <- ev$NMDS2 * 0.8 * max(abs(sc$NMDS2))
        p_cell <- p_cell +
          geom_segment(data = ev, inherit.aes = FALSE, aes(0, 0, xend = NMDS1, yend = NMDS2), arrow = arrow(length = unit(0.18, "cm")), colour = "grey25") +
          geom_text(data = ev, inherit.aes = FALSE, aes(NMDS1, NMDS2, label = prey), size = 2.5, colour = "grey15")
      }
      p_cell <- p_cell + scale_colour_manual(values = PER_COL, name = NULL) +
        labs(title = paste0("Cell centroids, ", CURRENCY_LAB[[cur]], " currency"),
             subtitle = sprintf(paste0("One point per predator x size class x ecoregion and decade;\nBray-Curtis, stress = %.2f. ",
                                       "Grey arrows link the same cell across decades;\ndark arrows: prey meeting both criteria ",
                                       "(direction = correlation with the axes)."), s)) +
        theme_diag(base_size = 11) + theme(panel.grid = element_blank(), legend.position = "top")
      save_drv(p_cell, sprintf("Fig_nmds_cells_%s", cur), 7.4, 6.4)

      ado_c <- run_permanova(dist_bray(cmat), list(period = cmeta$period), blocks = cmeta$cell, n_perm = CFG$perms,
                             label = paste("PERMANOVA cells", cur))
      if (!is.null(ado_c)) {
        print(ado_c); capture.output(print(ado_c), file = file.path(DIR_DRIVERS, sprintf("permanova_cellpaired_%s.txt", cur)))
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
  pick <- function(d, sfx) { z <- d[, .(prey, mglm_p, wilcox_p_bh, d_share, direction, n_criteria, pct_modal)]
  old <- setdiff(names(z), "prey"); setnames(z, old, paste0(old, "_", sfx)); z }
  xc <- merge(pick(DRV$biomass, "B"), pick(DRV$occurrence, "O"), by = "prey", all = TRUE)
  xc[, `:=`(sig_B = !is.na(n_criteria_B) & n_criteria_B >= 1, sig_O = !is.na(n_criteria_O) & n_criteria_O >= 1)]
  xc[, class := fcase(
    sig_B & sig_O & direction_B == direction_O & n_criteria_B == 2 & n_criteria_O == 2, "A+B in both currencies",
    sig_B & sig_O & direction_B == direction_O, "robust (both currencies)",
    sig_B & sig_O, "discordant", sig_B | sig_O, "one currency", default = "not significant")]
  xc[, d_share_mean := rowMeans(cbind(d_share_B, d_share_O), na.rm = TRUE)]
  xc[, class := factor(class, levels = c("A+B in both currencies", "robust (both currencies)", "discordant", "one currency", "not significant"))]
  xc[, abs_d := abs(d_share_mean)]; setorder(xc, class, -abs_d); xc[, abs_d := NULL]
  write_csv(xc, file.path(DIR_DRIVERS, "drivers_cross_currency.csv"))
  cat("\nPrey by robustness class:\n"); print(table(xc$class))
  cat("\nRobust in both currencies:\n")
  print(as.data.frame(xc[class %in% levels(class)[1:2], .(prey, class, direction = direction_B,
                                                          d_pct_B = round(100 * d_share_B, 2), d_pct_O = round(100 * d_share_O, 2))]), row.names = FALSE)
  SUMMARY$cross_currency <- as.list(table(xc$class))
} else message("Cross-currency join skipped: one currency has no driver table.")

# =============================================================================
# 6. RUN SUMMARY
# =============================================================================
con <- file(file.path(DIR_DRIVERS, "run_summary.txt"), "w")
writeLines(c("10_Prey_Drivers.R - run summary (cell-consistent design)", format(Sys.time()),
             paste0("Engine: ", ENGINE), paste0("Prey family: ", CFG$prey_family, " | manuscript resolution: ", CFG$res_col),
             paste0("Sweep: ", paste(sweep_specs$label, collapse = ", ")),
             paste0("Harmonised categories: ", if (length(CFG$harmonise)) paste(sapply(names(CFG$harmonise), function(to)
               paste0(paste(CFG$harmonise[[to]], collapse = "+"), " -> ", to)), collapse = "; ") else "none"),
             paste0("Permutations: ", CFG$perms, " (Test A) / ", CFG$perms_cell, " (IndVal in cells) / ", CFG$perms_sweep, " (sweep)"),
             paste0("Testable cell: >= ", CFG$min_rows_per_period, " rows per period | criterion B: BH p <= ", CFG$alpha,
                    " and >= ", 100 * CFG$modal_share, " % of cells in the same direction"), ""), con)
for (nm in setdiff(names(SUMMARY), "engine")) {
  writeLines(paste0("[", nm, "]"), con)
  for (k in names(SUMMARY[[nm]])) writeLines(paste0("  ", k, ": ", paste(unlist(SUMMARY[[nm]][[k]]), collapse = ", ")), con)
  writeLines("", con)
}
if (file.exists(FAIL_LOG)) writeLines(c("Failed tests:", paste0("  ", readLines(FAIL_LOG))), con) else writeLines("No test failed.", con)
close(con)

cat("\n", strrep("=", 70), "\nOutputs in ", normalizePath(DIR_DRIVERS), "\n\n", sep = "")
cat("Reporting guidance:\n")
cat("  * Association, not causation: 'prey associated with the transition'.\n")
cat("  * A prey earns a place in Results when it meets A AND B at x = 440 (harmonised),\n")
cat("    is consistent at the 'order' anchor, and is 'robust' or 'A+B' across currencies.\n")
cat("  * mglm_p vs mglm_p_unadjusted: prey significant only unadjusted are the\n")
cat("    sampling-mix artefacts a pooled design would have reported.\n")
cat("  * A label present in one decade only is an identification-level change:\n")
cat("    Methods caveat (with the identification-reach tables), not a Result.\n")
