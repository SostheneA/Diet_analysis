## =====================================================================
## 12_Network_Rewiring_v3.R
## Interaction beta-diversity of predator-prey networks from stomach data
##
##   Poisot axis   : prey turnover (beta_ST) vs rewiring (beta_OS)
##   Baselga axis  : balanced replacement vs intensity gradient inside beta_OS
##   Inference     : within-period split-half null at matched effort
##
## Self-contained. Supersedes v2 + the three patch files.
## Project: Diet_analysis | Input: data/dat_classed.rda (script 5)
##
## Replicate unit is the trawl set throughout, as in 6a. Link weight is the
## mean across sets of the per-set diet composition, so a set with 40
## stomachs and a set with 3 contribute equally.
## =====================================================================

## ---------------------------------------------------------------------
## 0. CONFIG
## ---------------------------------------------------------------------

if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"        # "_1" | "_2" | "_PP"
if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"  # all_gulf | ecoregion | stratum

DAT_PATH <- file.path("data", "dat_classed.rda")
P1_RUNS  <- sprintf("all_runs_%s%s.rda", SPATIAL_LEVEL, PREY_FAMILY)   # optional

## Run the matching diagnostic, then stop so the grain can be chosen on
## numbers. Set to FALSE for the real run.
STOP_AFTER_CHECK <- F

## --- resolved column names (from the real dat_classed inventory) ------
COL_PRED    <- "predator_species_common_name"
COL_SIZE    <- "size_class"
COL_PERIOD  <- "period"
COL_STOMACH <- "stomach_id"
COL_YEAR    <- "year"
COL_VESSEL  <- "vessel.code"
COL_SETNO   <- "set"
COL_AREA    <- "Area"
COL_STRATUM <- "stratum"        # 'strata' exists but is empty

## Biomass currency. Must match what 6a uses, or the two papers report
## "biomass" on different scales.
##   "pfi"          partial fullness index (prey weight / predator weight)
##   "somatic_wt_g" raw prey wet weight
BIOMASS_CURRENCY <- "somatic_wt_g"

## Spatial matching grain: "stratum" (strict), "Area" (gross shifts), "none"
MATCH_COL <- "Area"

## --- filters ----------------------------------------------------------
DROP_EMPTY       <- TRUE
NUTRITIONAL_ONLY <- TRUE
DROP_GEOMETRY    <- TRUE        # dat_classed is sf; geometry kills dplyr speed

## --- contrasts: period holds LABELS, not years ------------------------
CONTRASTS <- list(
  PT = list(p1 = "2004-2006", p2 = "2018-2019")
)
MODES <- c("biomass", "occurrence")

## --- design switches --------------------------------------------------
NODE_DEF      <- "species"              # "species" | "species_size"
RENORM_SHARED <- TRUE
X_THRESHOLDS  <- seq(10, 1000, by = 20)  # NULL = all 100 resolutions

MIN_SETS_PRED <- 3
B_RAREFY      <- 200
B_NULL        <- 200
N_NULL_NET    <- 199
SEED          <- 20260914
DO_MODULES    <- FALSE                  # slow; enable once verified
VERBOSE       <- TRUE

SPATIAL_COL <- switch(SPATIAL_LEVEL,
                      all_gulf = NULL, ecoregion = COL_AREA, stratum = COL_STRATUM,
                      stop("Unknown SPATIAL_LEVEL: ", SPATIAL_LEVEL))

OUT_DIR  <- sprintf("Network%s_%s", PREY_FAMILY, SPATIAL_LEVEL)
PLOT_DIR <- sprintf("Network_Plot%s_%s", PREY_FAMILY, SPATIAL_LEVEL)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})
set.seed(SEED)
msg <- function(...) if (VERBOSE) cat(sprintf(...), "\n")

## ---------------------------------------------------------------------
## 1. LOAD, FILTER, DERIVE
## ---------------------------------------------------------------------

stopifnot(file.exists(DAT_PATH))
.obj <- load(DAT_PATH); dat <- get(.obj[1])

if (DROP_GEOMETRY && inherits(dat, "sf")) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("sf needed to drop geometry.")
  dat <- sf::st_drop_geometry(dat); msg("Dropped sf geometry.")
}
dat <- as.data.frame(dat); nms <- names(dat)

need <- c(COL_PRED, COL_SIZE, COL_PERIOD, COL_STOMACH, COL_YEAR, COL_VESSEL, COL_SETNO)
if (length(setdiff(need, nms))) stop("Missing: ", paste(setdiff(need, nms), collapse = ", "))

n0 <- nrow(dat)
truthy <- function(v) v %in% c(TRUE, 1, "1", "Y", "y", "yes", "TRUE", "true", "O", "oui")

if (DROP_EMPTY && "is_empty" %in% nms) {
  drop <- truthy(dat$is_empty)
  msg("Empty-stomach rows dropped: %d", sum(drop, na.rm = TRUE))
  dat <- dat[!drop, , drop = FALSE]
}
if (NUTRITIONAL_ONLY && "is_nutritional_prey" %in% nms) {
  keep <- truthy(dat$is_nutritional_prey)
  msg("Non-nutritional prey rows dropped: %d", sum(!keep, na.rm = TRUE))
  if (mean(keep) < 0.01)
    warning("is_nutritional_prey matched almost nothing - check its coding before trusting this filter.")
  else dat <- dat[keep, , drop = FALSE]
}
msg("Rows: %d -> %d", n0, nrow(dat))

## biomass currency
col_wgt <- NULL
if (BIOMASS_CURRENCY %in% nms) {
  col_wgt <- ".wgt"
  dat$.wgt <- suppressWarnings(as.numeric(dat[[BIOMASS_CURRENCY]]))
  bad <- sum(is.na(dat$.wgt) | dat$.wgt < 0)
  dat$.wgt[is.na(dat$.wgt) | dat$.wgt < 0] <- 0
  msg("Biomass currency '%s': %d rows NA/negative set to 0.", BIOMASS_CURRENCY, bad)
  if (sum(dat$.wgt) <= 0) stop("Biomass currency is all zero.")
} else {
  warning(sprintf("Currency '%s' absent - occurrence mode only.", BIOMASS_CURRENCY))
  MODES <- setdiff(MODES, "biomass")
  if (!length(MODES)) stop("No usable mode left.")
}

## derived keys. The trawl set is the replicate unit - never stomach_id.
dat$.set     <- paste(dat[[COL_YEAR]], dat[[COL_VESSEL]], dat[[COL_SETNO]], sep = "_")
dat$.stom    <- as.character(dat[[COL_STOMACH]])
dat$.period  <- as.character(dat[[COL_PERIOD]])
dat$.unit    <- if (is.null(SPATIAL_COL)) "all_gulf" else as.character(dat[[SPATIAL_COL]])
dat$.node    <- if (NODE_DEF == "species_size")
  paste(dat[[COL_PRED]], dat[[COL_SIZE]], sep = " | ") else as.character(dat[[COL_PRED]])

if (MATCH_COL == "none" || !MATCH_COL %in% nms) {
  if (MATCH_COL != "none") warning("MATCH_COL '", MATCH_COL, "' absent - matching disabled.")
  use_strata <- FALSE; dat$.stratum <- "all"
  warning("No spatial matching: set locations differ between periods, ",
          "so part of the measured rewiring may be spatial rather than temporal.")
} else {
  use_strata <- TRUE
  dat$.stratum <- as.character(dat[[MATCH_COL]])
  dat$.stratum[is.na(dat$.stratum)] <- "unknown"
}
msg("Matching grain: %s (%d groups).", MATCH_COL, length(unique(dat$.stratum)))

## sanity
.chk <- tapply(dat$.stom, dat$.set, function(z) length(unique(z)))
msg("Sets: %d | stomachs/set: median %.1f, max %d", length(.chk), median(.chk), max(.chk))
if (median(.chk) <= 1) stop("Median 1 stomach per set - the set key is wrong.")

.pv <- sort(unique(dat$.period))
msg("Period labels: %s", paste(.pv, collapse = ", "))
for (ct in names(CONTRASTS)) {
  miss <- setdiff(unlist(CONTRASTS[[ct]]), .pv)
  if (length(miss)) stop(sprintf("Contrast %s references absent labels: %s",
                                 ct, paste(miss, collapse = ", ")))
}
msg("Sets per period: %s", paste(sprintf("%s=%d", names(table(dat$.period[!duplicated(dat$.set)])),
                                         table(dat$.period[!duplicated(dat$.set)])), collapse = ", "))
msg("Predators: %d | spatial units: %d", length(unique(dat$.node)), length(unique(dat$.unit)))

## prey resolution columns
prey_cols <- grep(sprintf("^prey_category_[0-9]+%s$", PREY_FAMILY), nms, value = TRUE)
if (!length(prey_cols)) stop("No prey_category_*", PREY_FAMILY, " columns.")
x_of <- function(cl) as.numeric(sub(sprintf("^prey_category_([0-9]+)%s$", PREY_FAMILY), "\\1", cl))
if (!is.null(X_THRESHOLDS)) {
  got <- intersect(sprintf("prey_category_%d%s", X_THRESHOLDS, PREY_FAMILY), prey_cols)
  if (!length(got)) stop("None of the requested thresholds exist.")
  prey_cols <- got
}
prey_cols <- prey_cols[order(x_of(prey_cols))]
msg("Resolutions: %d | modes: %s", length(prey_cols), paste(MODES, collapse = ", "))

## ---------------------------------------------------------------------
## 2. MATCHING-GRAIN DIAGNOSTIC
## ---------------------------------------------------------------------
## Matching drops any group sampled in only one period, per predator.
## This is what that control costs at each available grain.

.p1 <- CONTRASTS[[1]]$p1; .p2 <- CONTRASTS[[1]]$p2
grain_retention <- do.call(rbind, lapply(intersect(c("none", COL_AREA, COL_STRATUM), c("none", nms)),
                                         function(mc) {
                                           g <- if (mc == "none") rep("all", nrow(dat)) else as.character(dat[[mc]])
                                           g[is.na(g)] <- "unknown"
                                           u <- unique(data.frame(node = dat$.node, set = dat$.set, per = dat$.period,
                                                                  grp = g, stringsAsFactors = FALSE))
                                           u <- u[u$per %in% c(.p1, .p2), ]
                                           u$per <- ifelse(u$per %in% .p1, "p1", "p2")
                                           do.call(rbind, lapply(split(u, u$node), function(z) {
                                             t1 <- table(z$grp[z$per == "p1"]); t2 <- table(z$grp[z$per == "p2"])
                                             cm <- intersect(names(t1), names(t2))
                                             k  <- if (length(cm)) sum(pmin(t1[cm], t2[cm])) else 0
                                             n1 <- length(unique(z$set[z$per == "p1"])); n2 <- length(unique(z$set[z$per == "p2"]))
                                             data.frame(match_col = mc, node = z$node[1], n1 = n1, n2 = n2,
                                                        k_matched = as.numeric(k), k_unmatched = min(n1, n2),
                                                        loss = 1 - as.numeric(k) / max(min(n1, n2), 1),
                                                        stringsAsFactors = FALSE)
                                           }))
                                         }))

cat("\n=== effort retained per matching grain ===\n")
print(do.call(rbind, lapply(split(grain_retention, grain_retention$match_col), function(z)
  data.frame(match_col = z$match_col[1],
             predators_kept = sum(z$k_matched >= MIN_SETS_PRED), predators_total = nrow(z),
             sets_matched = sum(z$k_matched), median_k = median(z$k_matched),
             median_loss = round(median(z$loss), 3)))), row.names = FALSE)
write.csv(grain_retention, file.path(OUT_DIR, "matching_grain_retention.csv"), row.names = FALSE)

if (STOP_AFTER_CHECK)
  stop("STOP_AFTER_CHECK is TRUE. Pick MATCH_COL from the table above, set it to FALSE, rerun.")

## ---------------------------------------------------------------------
## 3. SET-LEVEL COMPOSITION AND NETWORK CONSTRUCTION
## ---------------------------------------------------------------------

set_prey_table <- function(d, prey_col, mode) {
  d <- d[!is.na(d[[prey_col]]) & d[[prey_col]] != "", , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$.prey <- as.character(d[[prey_col]])
  if (mode == "biomass") {
    tab <- d %>% group_by(.set, .node, .stratum, .prey) %>%
      summarise(v = sum(.wgt, na.rm = TRUE), .groups = "drop")
  } else {
    nst <- d %>% group_by(.set, .node) %>%
      summarise(n_stom = n_distinct(.stom), .groups = "drop")
    tab <- d %>% group_by(.set, .node, .stratum, .prey) %>%
      summarise(n_occ = n_distinct(.stom), .groups = "drop") %>%
      left_join(nst, by = c(".set", ".node")) %>%
      mutate(v = n_occ / n_stom) %>% select(-n_occ, -n_stom)
  }
  tab %>% group_by(.set, .node) %>% mutate(tot = sum(v, na.rm = TRUE)) %>%
    ungroup() %>% filter(tot > 0) %>% mutate(v = v / tot) %>% select(-tot)
}

net_from_sets <- function(tab, sets_by_node) {
  keep <- bind_rows(lapply(names(sets_by_node), function(p)
    data.frame(.node = p, .set = sets_by_node[[p]], stringsAsFactors = FALSE)))
  if (!nrow(keep)) return(NULL)
  sub <- inner_join(tab, keep, by = c(".node", ".set"))
  if (!nrow(sub)) return(NULL)
  ns  <- sub %>% distinct(.node, .set) %>% count(.node, name = "n_sets")
  agg <- sub %>% group_by(.node, .prey) %>% summarise(w = sum(v), .groups = "drop") %>%
    left_join(ns, by = ".node") %>% mutate(w = w / n_sets) %>% select(.node, .prey, w)
  m <- xtabs(w ~ .node + .prey, data = agg)
  matrix(as.numeric(m), nrow = nrow(m), dimnames = dimnames(m))
}

## ---------------------------------------------------------------------
## 4. EFFORT PLAN AND NULL DONORS
## ---------------------------------------------------------------------
## k_obs  full matched effort = min(n1, n2) sets per group per predator
## k_test effort at which the test runs; equals k_obs whenever one period
##        can supply 2*k_obs disjoint sets for the null
## donors which period(s) supply the within-period null

effort_plan2 <- function(idx1, idx2, nodes) {
  out <- list()
  for (nd in nodes) {
    d1 <- idx1[[nd]]; d2 <- idx2[[nd]]
    if (is.null(d1) || is.null(d2)) next
    if (use_strata) {
      t1 <- table(d1$.stratum); t2 <- table(d2$.stratum)
      cm <- intersect(names(t1), names(t2)); if (!length(cm)) next
      t1 <- t1[cm]; t2 <- t2[cm]
    } else { t1 <- c(all = nrow(d1)); t2 <- c(all = nrow(d2)) }

    k_obs <- pmin(t1, t2); k_obs <- k_obs[k_obs >= 1]
    if (!length(k_obs) || sum(k_obs) < 3) next
    nm <- names(k_obs)
    can1 <- all(floor(t1[nm] / 2) >= k_obs); can2 <- all(floor(t2[nm] / 2) >= k_obs)
    if (can1 || can2) {
      k_test <- k_obs; donors <- c(if (can1) 1L, if (can2) 2L)
    } else {
      k_test <- pmin(floor(t1[nm] / 2), floor(t2[nm] / 2)); k_test <- k_test[k_test >= 1]
      donors <- c(1L, 2L)
      if (!length(k_test) || sum(k_test) < 3) next
    }
    out[[nd]] <- list(k_obs = as.list(k_obs), k_test = as.list(k_test), donors = donors,
                      full_effort = identical(as.numeric(k_test),
                                              as.numeric(k_obs[names(k_test)])))
  }
  out
}

draw_two_periods <- function(idx1, idx2, plan, field) {
  s1 <- list(); s2 <- list()
  for (nd in names(plan)) {
    k <- plan[[nd]][[field]]; d1 <- idx1[[nd]]; d2 <- idx2[[nd]]
    for (s in names(k)) {
      kk <- k[[s]]; if (kk < 1) next
      a <- d1$.set[d1$.stratum == s]; b <- d2$.set[d2$.stratum == s]
      if (length(a) < kk || length(b) < kk) return(NULL)
      s1[[nd]] <- c(s1[[nd]], if (length(a) == kk) a else sample(a, kk))
      s2[[nd]] <- c(s2[[nd]], if (length(b) == kk) b else sample(b, kk))
    }
  }
  list(p1 = s1, p2 = s2)
}

draw_split_half <- function(idx, plan) {
  sa <- list(); sb <- list()
  for (nd in names(plan)) {
    k <- plan[[nd]]$k_test; d <- idx[[nd]]; if (is.null(d)) return(NULL)
    for (s in names(k)) {
      kk <- k[[s]]; if (kk < 1) next
      pool <- d$.set[d$.stratum == s]
      if (length(pool) < 2 * kk) return(NULL)
      pick <- sample(pool, 2 * kk)
      sa[[nd]] <- c(sa[[nd]], pick[seq_len(kk)])
      sb[[nd]] <- c(sb[[nd]], pick[(kk + 1):(2 * kk)])
    }
  }
  list(a = sa, b = sb)
}

## ---------------------------------------------------------------------
## 5. DISSIMILARITIES
## ---------------------------------------------------------------------

.pad <- function(m, rows, cols) {
  o <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
  rr <- intersect(rows, rownames(m)); cc <- intersect(cols, colnames(m))
  if (length(rr) && length(cc)) o[rr, cc] <- m[rr, cc, drop = FALSE]
  o
}

## Baselga (2013): BC = balanced variation + abundance gradient
bc_decomp <- function(v1, v2) {
  A <- sum(pmin(v1, v2)); B <- sum(v1 - pmin(v1, v2)); C <- sum(v2 - pmin(v1, v2))
  den <- 2 * A + B + C
  if (den <= 0) return(c(bc = NA_real_, bal = NA_real_, gra = NA_real_))
  bc <- (B + C) / den
  bal <- if ((A + min(B, C)) == 0) 0 else min(B, C) / (A + min(B, C))
  c(bc = bc, bal = bal, gra = bc - bal)
}

bin_diss <- function(v1, v2) {
  a <- sum(v1 > 0 & v2 > 0); b <- sum(v1 > 0 & v2 == 0); cc <- sum(v1 == 0 & v2 > 0)
  if ((2 * a + b + cc) == 0) return(NA_real_)
  (b + cc) / (2 * a + b + cc)
}

beta_full <- function(m1, m2) {
  nodes <- union(rownames(m1), rownames(m2)); prey <- union(colnames(m1), colnames(m2))
  M1 <- .pad(m1, nodes, prey); M2 <- .pad(m2, nodes, prey)
  wn  <- bc_decomp(as.vector(M1), as.vector(M2))
  wnb <- bin_diss(as.vector(M1), as.vector(M2))

  p1 <- colnames(m1); p2 <- colnames(m2)
  a <- length(intersect(p1, p2)); b <- length(setdiff(p1, p2)); cc <- length(setdiff(p2, p1))
  beta_S <- if ((2 * a + b + cc) == 0) NA_real_ else (b + cc) / (2 * a + b + cc)

  sp <- intersect(colnames(m1), colnames(m2)); sn <- intersect(rownames(m1), rownames(m2))
  os <- c(bc = NA_real_, bal = NA_real_, gra = NA_real_); osb <- NA_real_
  os_rn <- os; per_node <- NULL
  if (length(sp) && length(sn)) {
    S1 <- m1[sn, sp, drop = FALSE]; S2 <- m2[sn, sp, drop = FALSE]
    os  <- bc_decomp(as.vector(S1), as.vector(S2))
    osb <- bin_diss(as.vector(S1), as.vector(S2))
    if (RENORM_SHARED) {
      r1 <- S1 / pmax(rowSums(S1), .Machine$double.eps)
      r2 <- S2 / pmax(rowSums(S2), .Machine$double.eps)
      os_rn <- bc_decomp(as.vector(r1), as.vector(r2))
    }
    per_node <- t(vapply(sn, function(p) bc_decomp(S1[p, ], S2[p, ]), numeric(3)))
  }
  c(list(beta_WN = unname(wn["bc"]), beta_WN_bal = unname(wn["bal"]),
         beta_WN_gra = unname(wn["gra"]), beta_WN_bin = wnb,
         beta_OS = unname(os["bc"]), beta_OS_bal = unname(os["bal"]),
         beta_OS_gra = unname(os["gra"]), beta_OS_bin = osb,
         beta_OS_rn = unname(os_rn["bc"]), beta_OS_rn_bal = unname(os_rn["bal"]),
         beta_ST = unname(wn["bc"]) - unname(os["bc"]),
         beta_ST_bin = if (is.na(wnb) || is.na(osb)) NA_real_ else wnb - osb,
         beta_S = beta_S, n_prey_1 = ncol(m1), n_prey_2 = ncol(m2),
         n_shared_prey = length(sp), n_nodes = length(sn)),
    list(per_node = per_node))
}

BETA_FIELDS <- c("beta_WN", "beta_WN_bal", "beta_WN_gra", "beta_WN_bin",
                 "beta_OS", "beta_OS_bal", "beta_OS_gra", "beta_OS_bin",
                 "beta_OS_rn", "beta_OS_rn_bal", "beta_ST", "beta_ST_bin", "beta_S")

net_metrics <- function(m, do_modules = FALSE, n_null = 0) {
  bin <- (m > 0) * 1
  out <- list(connectance = sum(bin) / (nrow(bin) * ncol(bin)),
              links_per_sp = sum(bin) / (nrow(bin) + ncol(bin)),
              n_pred = nrow(m), n_prey = ncol(m),
              nodf = NA_real_, nodf_z = NA_real_,
              modularity = NA_real_, modularity_z = NA_real_)
  if (!requireNamespace("bipartite", quietly = TRUE)) return(out)
  nodf_of <- function(z) tryCatch(unname(bipartite::networklevel(z, index = "weighted NODF")),
                                  error = function(e) NA_real_)
  mod_of <- function(z) tryCatch(bipartite::computeModules(z)@likelihood,
                                 error = function(e) NA_real_)
  out$nodf <- nodf_of(m); if (do_modules) out$modularity <- mod_of(m)
  if (n_null > 0) {
    nl <- tryCatch(bipartite::vaznull(n_null, m), error = function(e) NULL)
    if (!is.null(nl)) {
      nd <- vapply(nl, nodf_of, numeric(1))
      if (sum(!is.na(nd)) > 2 && sd(nd, na.rm = TRUE) > 0)
        out$nodf_z <- (out$nodf - mean(nd, na.rm = TRUE)) / sd(nd, na.rm = TRUE)
      if (do_modules) {
        md <- vapply(nl, mod_of, numeric(1))
        if (sum(!is.na(md)) > 2 && sd(md, na.rm = TRUE) > 0)
          out$modularity_z <- (out$modularity - mean(md, na.rm = TRUE)) / sd(md, na.rm = TRUE)
      }
    }
  }
  out
}

## ---------------------------------------------------------------------
## 6. ONE CELL
## ---------------------------------------------------------------------

run_cell <- function(d_unit, prey_col, mode, p1_lab, p2_lab) {
  d1 <- d_unit[d_unit$.period %in% p1_lab, , drop = FALSE]
  d2 <- d_unit[d_unit$.period %in% p2_lab, , drop = FALSE]
  if (!nrow(d1) || !nrow(d2)) return(NULL)
  t1 <- set_prey_table(d1, prey_col, mode); t2 <- set_prey_table(d2, prey_col, mode)
  if (is.null(t1) || is.null(t2)) return(NULL)

  mk_idx <- function(tt) {
    u <- distinct(tt, .node, .set, .stratum)
    lapply(split(u, u$.node), function(z) as.data.frame(z[, c(".set", ".stratum")]))
  }
  idx1 <- mk_idx(t1); idx2 <- mk_idx(t2)
  core <- intersect(names(idx1)[vapply(idx1, nrow, 1L) >= MIN_SETS_PRED],
                    names(idx2)[vapply(idx2, nrow, 1L) >= MIN_SETS_PRED])
  if (length(core) < 2) return(NULL)
  plan <- effort_plan2(idx1, idx2, core)
  if (length(plan) < 2) return(NULL)
  core <- names(plan)

  obs <- vector("list", B_RAREFY)
  for (b in seq_len(B_RAREFY)) {
    dd <- draw_two_periods(idx1, idx2, plan, "k_test"); if (is.null(dd)) next
    ma <- net_from_sets(t1, dd$p1); mb <- net_from_sets(t2, dd$p2)
    if (is.null(ma) || is.null(mb)) next
    obs[[b]] <- beta_full(ma, mb)
  }
  obs <- Filter(Negate(is.null), obs); if (!length(obs)) return(NULL)

  full_same <- all(vapply(plan, function(z) isTRUE(z$full_effort), logical(1)))
  obs_full <- obs
  if (!full_same) {
    obs_full <- vector("list", max(30, B_RAREFY %/% 4))
    for (b in seq_along(obs_full)) {
      dd <- draw_two_periods(idx1, idx2, plan, "k_obs"); if (is.null(dd)) next
      ma <- net_from_sets(t1, dd$p1); mb <- net_from_sets(t2, dd$p2)
      if (is.null(ma) || is.null(mb)) next
      obs_full[[b]] <- beta_full(ma, mb)
    }
    obs_full <- Filter(Negate(is.null), obs_full)
  }

  donors <- sort(unique(unlist(lapply(plan, `[[`, "donors"))))
  nulls <- list(); by_donor <- list()
  for (p in donors) {
    idx <- if (p == 1L) idx1 else idx2; tt <- if (p == 1L) t1 else t2
    acc <- vector("list", B_NULL)
    for (b in seq_len(B_NULL)) {
      sp <- draw_split_half(idx, plan); if (is.null(sp)) next
      ma <- net_from_sets(tt, sp$a); mb <- net_from_sets(tt, sp$b)
      if (is.null(ma) || is.null(mb)) next
      acc[[b]] <- beta_full(ma, mb)
    }
    acc <- Filter(Negate(is.null), acc)
    by_donor[[as.character(p)]] <- acc; nulls <- c(nulls, acc)
  }

  grab <- function(lst, f) vapply(lst, function(z) as.numeric(z[[f]]), numeric(1))
  row <- list()
  for (f in BETA_FIELDS) {
    v <- grab(obs, f)
    row[[paste0(f, "_mean")]] <- mean(v, na.rm = TRUE)
    row[[paste0(f, "_lo")]]   <- unname(quantile(v, .025, na.rm = TRUE))
    row[[paste0(f, "_hi")]]   <- unname(quantile(v, .975, na.rm = TRUE))
    if (length(obs_full)) row[[paste0(f, "_full")]] <- mean(grab(obs_full, f), na.rm = TRUE)
    if (length(nulls)) {
      vn <- grab(nulls, f)
      row[[paste0(f, "_null")]]    <- mean(vn, na.rm = TRUE)
      row[[paste0(f, "_null_hi")]] <- unname(quantile(vn, .975, na.rm = TRUE))
      row[[paste0(f, "_p")]] <- (1 + sum(vn >= mean(v, na.rm = TRUE), na.rm = TRUE)) /
        (1 + sum(!is.na(vn)))
      row[[paste0(f, "_excess")]] <- mean(v, na.rm = TRUE) - mean(vn, na.rm = TRUE)
    }
  }
  for (p in names(by_donor)) if (length(by_donor[[p]]))
    row[[paste0("beta_OS_null_p", p)]] <- mean(grab(by_donor[[p]], "beta_OS"), na.rm = TRUE)

  dd <- draw_two_periods(idx1, idx2, plan, "k_obs")
  metrics <- NULL
  if (!is.null(dd)) {
    m1o <- net_from_sets(t1, dd$p1); m2o <- net_from_sets(t2, dd$p2)
    if (!is.null(m1o) && !is.null(m2o)) metrics <- bind_rows(
      data.frame(period = "p1", as.data.frame(net_metrics(m1o, DO_MODULES, N_NULL_NET))),
      data.frame(period = "p2", as.data.frame(net_metrics(m2o, DO_MODULES, N_NULL_NET))))
  }

  pn <- lapply(obs, function(z) z$per_node); pn <- pn[!vapply(pn, is.null, logical(1))]
  per_node <- NULL
  if (length(pn)) {
    allp <- unique(unlist(lapply(pn, rownames)))
    st <- lapply(c("bc", "bal", "gra"), function(f)
      do.call(rbind, lapply(pn, function(z) z[match(allp, rownames(z)), f])))
    per_node <- data.frame(node = allp, beta_OS = colMeans(st[[1]], na.rm = TRUE),
                           beta_OS_bal = colMeans(st[[2]], na.rm = TRUE),
                           beta_OS_gra = colMeans(st[[3]], na.rm = TRUE),
                           row.names = NULL, stringsAsFactors = FALSE)
    pnn <- lapply(nulls, function(z) z$per_node); pnn <- pnn[!vapply(pnn, is.null, logical(1))]
    if (length(pnn)) {
      nmz <- do.call(rbind, lapply(pnn, function(z) z[match(allp, rownames(z)), "bc"]))
      per_node$beta_OS_null <- colMeans(nmz, na.rm = TRUE)
      per_node$beta_OS_p <- vapply(seq_along(allp), function(i)
        (1 + sum(nmz[, i] >= per_node$beta_OS[i], na.rm = TRUE)) /
          (1 + sum(!is.na(nmz[, i]))), numeric(1))
    }
  }

  list(beta = as.data.frame(c(row, list(
    n_nodes = length(core), n_draws = length(obs), n_draws_full = length(obs_full),
    n_null = length(nulls), donors = paste(donors, collapse = "+"),
    full_effort = full_same,
    n_sets_test = sum(unlist(lapply(plan, function(z) unlist(z$k_test)))),
    n_sets_obs  = sum(unlist(lapply(plan, function(z) unlist(z$k_obs))))))),
    metrics = metrics, per_node = per_node)
}

## ---------------------------------------------------------------------
## 7. SWEEP
## ---------------------------------------------------------------------

res_beta <- list(); res_met <- list(); res_node <- list()
t_start <- Sys.time()

for (ctr_name in names(CONTRASTS)) {
  ctr <- CONTRASTS[[ctr_name]]
  for (mode in MODES) for (prey_col in prey_cols) {
    xth <- x_of(prey_col)
    for (u in sort(unique(dat$.unit))) {
      d_unit <- dat[dat$.unit == u, , drop = FALSE]
      r <- try(run_cell(d_unit, prey_col, mode, ctr$p1, ctr$p2), silent = TRUE)
      if (inherits(r, "try-error")) { msg("  ! %s x=%s unit=%s : %s", mode, xth, u,
                                          conditionMessage(attr(r, "condition"))); next }
      if (is.null(r)) next
      key <- data.frame(contrast = ctr_name, mode = mode, x_threshold = xth,
                        spatial_level = SPATIAL_LEVEL, unit = u, prey_family = PREY_FAMILY,
                        node_def = NODE_DEF, match_col = MATCH_COL, stringsAsFactors = FALSE)
      res_beta[[length(res_beta) + 1]] <- cbind(key, r$beta)
      if (!is.null(r$metrics))  res_met[[length(res_met) + 1]]   <- cbind(key, r$metrics)
      if (!is.null(r$per_node)) res_node[[length(res_node) + 1]] <- cbind(key, r$per_node)
    }
    msg("done | %s | %s | x=%s | %.1f min elapsed", ctr_name, mode, xth,
        as.numeric(difftime(Sys.time(), t_start, units = "mins")))
  }
}

beta_tab <- bind_rows(res_beta); met_tab <- bind_rows(res_met); node_tab <- bind_rows(res_node)
if (!nrow(beta_tab)) stop("No cell produced results - loosen MIN_SETS_PRED or the matching grain.")

## ---------------------------------------------------------------------
## 8. SIGNATURE AND OUTPUTS
## ---------------------------------------------------------------------

beta_tab <- beta_tab %>%
  mutate(os_share  = beta_OS_mean / pmax(beta_OS_mean + pmax(beta_ST_mean, 0), 1e-12),
         bal_share = beta_OS_bal_mean / pmax(beta_OS_mean, 1e-12),
         signal    = if ("beta_OS_p" %in% names(.)) beta_OS_p < 0.05 else NA,
         signature = case_when(
           !is.na(signal) & !signal            ~ "Indistinguishable from sampling",
           os_share <= 0.33                    ~ "Turnover-dominated (prey loss)",
           os_share >= 0.66 & bal_share >= .66 ~ "Rewiring, balanced (Ghost Shift)",
           os_share >= 0.66 & bal_share <  .66 ~ "Rewiring, intensity gradient",
           TRUE                                ~ "Mixed"))

saveRDS(list(beta = beta_tab, metrics = met_tab, per_node = node_tab,
             retention = grain_retention,
             config = list(PREY_FAMILY = PREY_FAMILY, SPATIAL_LEVEL = SPATIAL_LEVEL,
                           NODE_DEF = NODE_DEF, MATCH_COL = MATCH_COL,
                           BIOMASS_CURRENCY = BIOMASS_CURRENCY,
                           B_RAREFY = B_RAREFY, B_NULL = B_NULL, SEED = SEED)),
        file.path(OUT_DIR, sprintf("network_rewiring_v3%s_%s.rds", PREY_FAMILY, SPATIAL_LEVEL)))
write.csv(beta_tab, file.path(OUT_DIR, "beta_decomposition.csv"), row.names = FALSE)
write.csv(met_tab,  file.path(OUT_DIR, "network_metrics.csv"),    row.names = FALSE)
if (nrow(node_tab)) write.csv(node_tab, file.path(OUT_DIR, "per_node_rewiring.csv"), row.names = FALSE)

cat("\n=== signatures ===\n"); print(table(beta_tab$mode, beta_tab$signature))


if (nrow(node_tab)) {
  pred_summary_table <- node_tab %>%
    group_by(mode, unit, node) %>%
    summarise(
      mean_beta_OS = mean(beta_OS, na.rm = TRUE),
      mean_null_OS = if ("beta_OS_null" %in% names(node_tab)) mean(beta_OS_null, na.rm = TRUE) else NA_real_,
      prop_significant = if ("beta_OS_p" %in% names(node_tab)) mean(beta_OS_p < 0.05, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  write.csv(pred_summary_table, file.path(OUT_DIR, "predator_rewiring_significance_summary.csv"), row.names = FALSE)
  msg("Tableau récapitulatif de la significativité par prédateur exporté.")
}

## ---------------------------------------------------------------------
## 9. FIGURES
## ---------------------------------------------------------------------

theme_set(theme_bw(base_size = 11))

ggsave(file.path(PLOT_DIR, "Fig1_two_axes.png"),
       ggplot(beta_tab, aes(beta_ST_mean, beta_OS_mean, colour = bal_share, shape = mode)) +
         geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
         geom_point(size = 2, alpha = .85) +
         scale_colour_viridis_c(name = "balanced\nshare of OS", limits = c(0, 1)) +
         facet_wrap(~ unit) + coord_equal() +
         labs(x = expression(beta[ST]~"(prey turnover)"), y = expression(beta[OS]~"(rewiring)"),
              title = "Poisot axis (position), Baselga axis (colour)",
              subtitle = "Ghost Shift = upper-left and bright"),
       width = 9, height = 6.5, dpi = 300)

if ("beta_OS_null" %in% names(beta_tab))
  ggsave(file.path(PLOT_DIR, "Fig2_null_vs_observed.png"),
         ggplot(beta_tab, aes(x_threshold)) +
           geom_ribbon(aes(ymin = beta_OS_lo, ymax = beta_OS_hi), alpha = .2, fill = "steelblue") +
           geom_line(aes(y = beta_OS_mean), colour = "steelblue", linewidth = .8) +
           geom_line(aes(y = beta_OS_null), colour = "grey35", linetype = 2) +
           geom_line(aes(y = beta_OS_null_hi), colour = "grey65", linetype = 3) +
           facet_grid(mode ~ unit) +
           labs(x = "Prey-resolution threshold (X)", y = expression(beta[OS]),
                title = "Rewiring against the within-period null",
                subtitle = "Dashed = expected under stability at the same effort; dotted = null 97.5%"),
         width = 10, height = 6, dpi = 300)

ggsave(file.path(PLOT_DIR, "Fig3_stacked_partition.png"),
       beta_tab %>% select(mode, unit, x_threshold, Turnover = beta_ST_mean,
                           Balanced = beta_OS_bal_mean, Gradient = beta_OS_gra_mean) %>%
         pivot_longer(c(Turnover, Balanced, Gradient), names_to = "component", values_to = "value") %>%
         mutate(value = pmax(value, 0)) %>%
         ggplot(aes(factor(x_threshold), value, fill = component)) +
         geom_col(width = .9) + facet_grid(mode ~ unit, scales = "free_x") +
         labs(x = "Prey-resolution threshold (X)", y = expression(beta[WN]~"partitioned"),
              title = "Where the dissimilarity comes from, at every prey resolution") +
         theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 6)),
       width = 11, height = 6, dpi = 300)

if (nrow(node_tab)) {
  pn <- node_tab %>% group_by(mode, unit, node) %>%
    summarise(beta_OS = mean(beta_OS, na.rm = TRUE),
              null = if ("beta_OS_null" %in% names(node_tab))
                mean(beta_OS_null, na.rm = TRUE) else NA_real_,
              sig = if ("beta_OS_p" %in% names(node_tab)) mean(beta_OS_p < 0.05, na.rm = TRUE) > 0.5 else FALSE,
              .groups = "drop")

  p4 <- ggplot(pn, aes(reorder(node, beta_OS), beta_OS, fill = sig)) +
    geom_col(position = position_dodge(.8), width = .7) + coord_flip() +
    facet_wrap(~ unit) +
    scale_fill_manual(values = c("TRUE" = "#2b8cbe", "FALSE" = "grey70"),
                      labels = c("TRUE" = "Significant Ghost Shift", "FALSE" = "Not significant"),
                      name = "Null Model Test") +
    labs(x = NULL, y = expression(beta[OS]~"(Rewiring)"),
         title = "Per-predator rewiring among persisting prey",
         subtitle = "Crosses = within-period null at matched effort")

  if (all(!is.na(pn$null))) p4 <- p4 + geom_point(aes(y = null), shape = 4, size = 1.8, colour = "black", stroke = 1)

  ggsave(file.path(PLOT_DIR, "Fig4_per_predator_enhanced.png"), p4,
         width = 10, height = max(4, .25 * length(unique(pn$node))), dpi = 300)
}

msg("Done in %.1f min. Outputs in %s and %s",
    as.numeric(difftime(Sys.time(), t_start, units = "mins")), OUT_DIR, PLOT_DIR)

## =====================================================================
## READING THE OUTPUT
##   beta_OS_p     p against the within-period null. Not < 0.05 means the
##                 row is not interpretable, whatever the other columns say.
##   os_share      beta_OS / (beta_OS + beta_ST). High = rewiring.
##   bal_share     balanced share of beta_OS. High = prey A replaced by
##                 prey B. Low = the diet moved toward or away from the
##                 persisting prey as a block.
##   beta_OS_bin   link gain/loss among shared prey, weights ignored.
##   beta_OS_rn    rewiring after renormalising each predator's shared-prey
##                 row to 1: pure reshuffling.
##   *_full        same quantities at full matched effort, descriptive only.
##   beta_OS_null_p1 / _p2  per-donor nulls, to check the assumption that
##                 sampling noise is comparable between periods.
##
## CAVEATS FOR THE MANUSCRIPT
## 1. Effort is equalised in SETS, matched within MATCH_COL, not in
##    stomachs per set.
## 2. The predator node set is the both-period core, so beta_S is prey
##    turnover only; predator turnover here is a survey outcome.
## 3. Prey absence is absence in the stomachs examined. beta_ST is an upper
##    bound on true prey turnover.
## 4. The resolution sweep is the validity test, not an appendix: a
##    signature that vanishes at coarse grain was identification depth.
## 5. Nothing here separates behavioural from forced switching. Prey are
##    observed only through predators.
## =====================================================================
