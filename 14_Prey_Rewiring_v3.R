## =====================================================================
## 13_Prey_Rewiring_v3.R
## Interaction beta-diversity of predator-prey networks: Prey perspective
##
##   Poisot axis   : consumer turnover vs rewiring from the prey viewpoint
##   Baselga axis  : balanced replacement vs intensity gradient
##   Inference     : within-period split-half null at matched effort
##
## Project: Diet_analysis | Input: data/dat_classed.rda (script 5)
## =====================================================================

## ---------------------------------------------------------------------
## 0. CONFIG
## ---------------------------------------------------------------------

if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"        # "_1" | "_2" | "_PP"
if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"  # all_gulf | ecoregion | stratum

DAT_PATH <- file.path("data", "dat_classed.rda")

## --- resolved column names ---
COL_PRED    <- "predator_species_common_name"
COL_SIZE    <- "size_class"
COL_PERIOD  <- "period"
COL_STOMACH <- "stomach_id"
COL_YEAR    <- "year"
COL_VESSEL  <- "vessel.code"
COL_SETNO   <- "set"
COL_AREA    <- "Area"
COL_STRATUM <- "stratum"

BIOMASS_CURRENCY <- "pfi"
MATCH_COL        <- "Area"

DROP_EMPTY       <- TRUE
NUTRITIONAL_ONLY <- TRUE
DROP_GEOMETRY    <- TRUE

CONTRASTS <- list(
  PT = list(p1 = "2004-2006", p2 = "2018-2019")
)
MODES <- c("biomass", "occurrence")

RENORM_SHARED <- TRUE
X_THRESHOLDS  <- seq(10, 1000, by = 10)

MIN_SETS_PREY <- 3   # Minimum de sets par proie pour l'analyse
B_RAREFY      <- 200
B_NULL        <- 200
SEED          <- 20260914
VERBOSE       <- TRUE

SPATIAL_COL <- switch(SPATIAL_LEVEL,
                      all_gulf = NULL, ecoregion = COL_AREA, stratum = COL_STRATUM,
                      stop("Unknown SPATIAL_LEVEL: ", SPATIAL_LEVEL))

OUT_DIR  <- sprintf("PreyNetwork%s_%s", PREY_FAMILY, SPATIAL_LEVEL)
PLOT_DIR <- sprintf("PreyNetwork_Plot%s_%s", PREY_FAMILY, SPATIAL_LEVEL)
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
  dat <- sf::st_drop_geometry(dat); msg("Dropped sf geometry.")
}
dat <- as.data.frame(dat); nms <- names(dat)

truthy <- function(v) v %in% c(TRUE, 1, "1", "Y", "y", "yes", "TRUE", "true", "O", "oui")

if (DROP_EMPTY && "is_empty" %in% nms) {
  dat <- dat[!truthy(dat$is_empty), , drop = FALSE]
}
if (NUTRITIONAL_ONLY && "is_nutritional_prey" %in% nms) {
  dat <- dat[truthy(dat$is_nutritional_prey), , drop = FALSE]
}

if (BIOMASS_CURRENCY %in% nms) {
  dat$.wgt <- suppressWarnings(as.numeric(dat[[BIOMASS_CURRENCY]]))
  dat$.wgt[is.na(dat$.wgt) | dat$.wgt < 0] <- 0
} else {
  MODES <- setdiff(MODES, "biomass")
}


dat$.set      <- paste(dat[[COL_YEAR]], dat[[COL_VESSEL]], dat[[COL_SETNO]], sep = "_")
dat$.stom     <- as.character(dat[[COL_STOMACH]])
dat$.period   <- as.character(dat[[COL_PERIOD]])
dat$.unit     <- if (is.null(SPATIAL_COL)) "all_gulf" else as.character(dat[[SPATIAL_COL]])
dat$.consumer <- as.character(dat[[COL_PRED]])  # Le consommateur (prédateur)

if (MATCH_COL == "none" || !MATCH_COL %in% nms) {
  use_strata <- FALSE; dat$.stratum <- "all"
} else {
  use_strata <- TRUE
  dat$.stratum <- as.character(dat[[MATCH_COL]])
  dat$.stratum[is.na(dat$.stratum)] <- "unknown"
}

# --- Extraction indispensable des colonnes de résolutions de proies ---
prey_cols <- grep(sprintf("^prey_category_[0-9]+%s$", PREY_FAMILY), nms, value = TRUE)
if (!length(prey_cols)) stop("Aucune colonne prey_category_*", PREY_FAMILY, " trouvée.")
x_of <- function(cl) as.numeric(sub(sprintf("^prey_category_([0-9]+)%s$", PREY_FAMILY), "\\1", cl))
if (!is.null(X_THRESHOLDS)) {
  prey_cols <- intersect(sprintf("prey_category_%d%s", X_THRESHOLDS, PREY_FAMILY), prey_cols)
}
prey_cols <- prey_cols[order(x_of(prey_cols))]
msg("Résolutions de proies détectées : %d", length(prey_cols))

## ---------------------------------------------------------------------
## 2. SET-LEVEL COMPOSITION (Perspective Proie -> Consommateurs)
## ---------------------------------------------------------------------

set_prey_perspective_table <- function(d, prey_col, mode) {
  d <- d[!is.na(d[[prey_col]]) & d[[prey_col]] != "", , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$.prey <- as.character(d[[prey_col]])

  if (mode == "biomass") {
    tab <- d %>% group_by(.set, .prey, .stratum, .consumer) %>%
      summarise(v = sum(.wgt, na.rm = TRUE), .groups = "drop")
  } else {
    nst <- d %>% group_by(.set, .prey) %>%
      summarise(n_stom = n_distinct(.stom), .groups = "drop")
    tab <- d %>% group_by(.set, .prey, .stratum, .consumer) %>%
      summarise(n_occ = n_distinct(.stom), .groups = "drop") %>%
      left_join(nst, by = c(".set", ".prey")) %>%
      mutate(v = n_occ / n_stom) %>% select(-n_occ, -n_stom)
  }
  tab %>% group_by(.set, .prey) %>% mutate(tot = sum(v, na.rm = TRUE)) %>%
    ungroup() %>% filter(tot > 0) %>% mutate(v = v / tot) %>% select(-tot)
}

net_from_sets_prey <- function(tab, sets_by_prey) {
  keep <- bind_rows(lapply(names(sets_by_prey), function(p)
    data.frame(.prey = p, .set = sets_by_prey[[p]], stringsAsFactors = FALSE)))
  if (!nrow(keep)) return(NULL)
  sub <- inner_join(tab, keep, by = c(".prey", ".set"))
  if (!nrow(sub)) return(NULL)
  ns  <- sub %>% distinct(.prey, .set) %>% count(.prey, name = "n_sets")
  agg <- sub %>% group_by(.prey, .consumer) %>% summarise(w = sum(v), .groups = "drop") %>%
    left_join(ns, by = ".prey") %>% mutate(w = w / n_sets) %>% select(.prey, .consumer, w)
  m <- xtabs(w ~ .prey + .consumer, data = agg)
  matrix(as.numeric(m), nrow = nrow(m), dimnames = dimnames(m))
}

## ---------------------------------------------------------------------
## 3. EFFORT PLAN ET NULL DONORS (Côté Proies)
## ---------------------------------------------------------------------

effort_plan_prey <- function(idx1, idx2, prey_nodes) {
  out <- list()
  for (pr in prey_nodes) {
    d1 <- idx1[[pr]]; d2 <- idx2[[pr]]
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
    out[[pr]] <- list(k_obs = as.list(k_obs), k_test = as.list(k_test), donors = donors,
                      full_effort = identical(as.numeric(k_test), as.numeric(k_obs[names(k_test)])))
  }
  out
}

draw_two_periods_prey <- function(idx1, idx2, plan, field) {
  s1 <- list(); s2 <- list()
  for (pr in names(plan)) {
    k <- plan[[pr]][[field]]; d1 <- idx1[[pr]]; d2 <- idx2[[pr]]
    for (s in names(k)) {
      kk <- k[[s]]; if (kk < 1) next
      a <- d1$.set[d1$.stratum == s]; b <- d2$.set[d2$.stratum == s]
      if (length(a) < kk || length(b) < kk) return(NULL)
      s1[[pr]] <- c(s1[[pr]], if (length(a) == kk) a else sample(a, kk))
      s2[[pr]] <- c(s2[[pr]], if (length(b) == kk) b else sample(b, kk))
    }
  }
  list(p1 = s1, p2 = s2)
}

draw_split_half_prey <- function(idx, plan) {
  sa <- list(); sb <- list()
  for (pr in names(plan)) {
    k <- plan[[pr]]$k_test; d <- idx[[pr]]; if (is.null(d)) return(NULL)
    for (s in names(k)) {
      kk <- k[[s]]; if (kk < 1) next
      pool <- d$.set[d$.stratum == s]
      if (length(pool) < 2 * kk) return(NULL)
      pick <- sample(pool, 2 * kk)
      sa[[pr]] <- c(sa[[pr]], pick[seq_len(kk)])
      sb[[pr]] <- c(sb[[pr]], pick[(kk + 1):(2 * kk)])
    }
  }
  list(a = sa, b = sb)
}

## ---------------------------------------------------------------------
## 4. DISSIMILARITÉS
## ---------------------------------------------------------------------

.pad <- function(m, rows, cols) {
  o <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
  rr <- intersect(rows, rownames(m)); cc <- intersect(cols, colnames(m))
  if (length(rr) && length(cc)) o[rr, cc] <- m[rr, cc, drop = FALSE]
  o
}

bc_decomp <- function(v1, v2) {
  A <- sum(pmin(v1, v2)); B <- sum(v1 - pmin(v1, v2)); C <- sum(v2 - pmin(v1, v2))
  den <- 2 * A + B + C
  if (den <= 0) return(c(bc = NA_real_, bal = NA_real_, gra = NA_real_))
  bc <- (B + C) / den
  bal <- if ((A + min(B, C)) == 0) 0 else min(B, C) / (A + min(B, C))
  c(bc = bc, bal = bal, gra = bc - bal)
}

beta_full_prey_pers <- function(m1, m2) {
  preys <- union(rownames(m1), rownames(m2)); consumers <- union(colnames(m1), colnames(m2))
  M1 <- .pad(m1, preys, consumers); M2 <- .pad(m2, preys, consumers)
  wn  <- bc_decomp(as.vector(M1), as.vector(M2))

  p1 <- colnames(m1); p2 <- colnames(m2)
  a <- length(intersect(p1, p2)); b <- length(setdiff(p1, p2)); cc <- length(setdiff(p2, p1))
  beta_consumer_turnover <- if ((2 * a + b + cc) == 0) NA_real_ else (b + cc) / (2 * a + b + cc)

  sp <- intersect(colnames(m1), colnames(m2)); sn <- intersect(rownames(m1), rownames(m2))
  os <- c(bc = NA_real_, bal = NA_real_, gra = NA_real_)
  per_prey <- NULL
  if (length(sp) && length(sn)) {
    S1 <- m1[sn, sp, drop = FALSE]; S2 <- m2[sn, sp, drop = FALSE]
    os  <- bc_decomp(as.vector(S1), as.vector(S2))
    per_prey <- t(vapply(sn, function(p) bc_decomp(S1[p, ], S2[p, ]), numeric(3)))
  }
  c(list(beta_WN = unname(wn["bc"]), beta_WN_bal = unname(wn["bal"]),
         beta_WN_gra = unname(wn["gra"]),
         beta_OS = unname(os["bc"]), beta_OS_bal = unname(os["bal"]),
         beta_OS_gra = unname(os["gra"]),
         beta_ST = unname(wn["bc"]) - unname(os["bc"]),
         beta_Consumer_Turnover = beta_consumer_turnover,
         n_consumers_1 = ncol(m1), n_consumers_2 = ncol(m2),
         n_shared_consumers = length(sp), n_prey_nodes = length(sn)),
    list(per_prey = per_prey))
}

BETA_FIELDS_PREY <- c("beta_WN", "beta_WN_bal", "beta_WN_gra",
                      "beta_OS", "beta_OS_bal", "beta_OS_gra",
                      "beta_ST", "beta_Consumer_Turnover")

## ---------------------------------------------------------------------
## 5. EXÉCUTION CELLULAIRE
## ---------------------------------------------------------------------

run_cell_prey <- function(d_unit, prey_col, mode, p1_lab, p2_lab) {
  d1 <- d_unit[d_unit$.period %in% p1_lab, , drop = FALSE]
  d2 <- d_unit[d_unit$.period %in% p2_lab, , drop = FALSE]
  if (!nrow(d1) || !nrow(d2)) return(NULL)

  t1 <- set_prey_perspective_table(d1, prey_col, mode)
  t2 <- set_prey_perspective_table(d2, prey_col, mode)
  if (is.null(t1) || is.null(t2)) return(NULL)

  mk_idx <- function(tt) {
    u <- distinct(tt, .prey, .set, .stratum)
    lapply(split(u, u$.prey), function(z) as.data.frame(z[, c(".set", ".stratum")]))
  }
  idx1 <- mk_idx(t1); idx2 <- mk_idx(t2)
  core <- intersect(names(idx1)[vapply(idx1, nrow, 1L) >= MIN_SETS_PREY],
                    names(idx2)[vapply(idx2, nrow, 1L) >= MIN_SETS_PREY])
  if (length(core) < 2) return(NULL)

  plan <- effort_plan_prey(idx1, idx2, core)
  if (length(plan) < 2) return(NULL)
  core <- names(plan)

  obs <- vector("list", B_RAREFY)
  for (b in seq_len(B_RAREFY)) {
    dd <- draw_two_periods_prey(idx1, idx2, plan, "k_test"); if (is.null(dd)) next
    ma <- net_from_sets_prey(t1, dd$p1); mb <- net_from_sets_prey(t2, dd$p2)
    if (is.null(ma) || is.null(mb)) next
    obs[[b]] <- beta_full_prey_pers(ma, mb)
  }
  obs <- Filter(Negate(is.null), obs); if (!length(obs)) return(NULL)

  donors <- sort(unique(unlist(lapply(plan, `[[`, "donors"))))
  nulls <- list()
  for (p in donors) {
    idx <- if (p == 1L) idx1 else idx2; tt <- if (p == 1L) t1 else t2
    acc <- vector("list", B_NULL)
    for (b in seq_len(B_NULL)) {
      sp <- draw_split_half_prey(idx, plan); if (is.null(sp)) next
      ma <- net_from_sets_prey(tt, sp$a); mb <- net_from_sets_prey(tt, sp$b)
      if (is.null(ma) || is.null(mb)) next
      acc[[b]] <- beta_full_prey_pers(ma, mb)
    }
    acc <- Filter(Negate(is.null), acc)
    nulls <- c(nulls, acc)
  }

  grab <- function(lst, f) vapply(lst, function(z) as.numeric(z[[f]]), numeric(1))
  row <- list()
  for (f in BETA_FIELDS_PREY) {
    v <- grab(obs, f)
    row[[paste0(f, "_mean")]] <- mean(v, na.rm = TRUE)
    row[[paste0(f, "_lo")]]   <- unname(quantile(v, .025, na.rm = TRUE))
    row[[paste0(f, "_hi")]]   <- unname(quantile(v, .975, na.rm = TRUE))
    if (length(nulls)) {
      vn <- grab(nulls, f)
      row[[paste0(f, "_null")]]    <- mean(vn, na.rm = TRUE)
      row[[paste0(f, "_p")]] <- (1 + sum(vn >= mean(v, na.rm = TRUE), na.rm = TRUE)) /
        (1 + sum(!is.na(vn)))
    }
  }

  pp <- lapply(obs, function(z) z$per_prey); pp <- pp[!vapply(pp, is.null, logical(1))]
  per_prey_tab <- NULL
  if (length(pp)) {
    allpr <- unique(unlist(lapply(pp, rownames)))
    st <- lapply(c("bc", "bal", "gra"), function(f)
      do.call(rbind, lapply(pp, function(z) z[match(allpr, rownames(z)), f])))
    per_prey_tab <- data.frame(prey = allpr, beta_OS_consumer = colMeans(st[[1]], na.rm = TRUE),
                               beta_OS_consumer_bal = colMeans(st[[2]], na.rm = TRUE),
                               beta_OS_consumer_gra = colMeans(st[[3]], na.rm = TRUE),
                               row.names = NULL, stringsAsFactors = FALSE)
  }

  list(beta = as.data.frame(c(row, list(n_prey_nodes = length(core), n_draws = length(obs)))),
       per_prey = per_prey_tab)
}

## ---------------------------------------------------------------------
## 6. SWEEP ET SAUVEGARDE
## ---------------------------------------------------------------------

res_prey_beta <- list(); res_prey_node <- list()
msg("Lancement de l'analyse de la perspective des proies...")

for (ctr_name in names(CONTRASTS)) {
  ctr <- CONTRASTS[[ctr_name]]
  for (mode in MODES) for (prey_col in prey_cols) {
    xth <- x_of(prey_col)
    for (u in sort(unique(dat$.unit))) {
      d_unit <- dat[dat$.unit == u, , drop = FALSE]
      r <- try(run_cell_prey(d_unit, prey_col, mode, ctr$p1, ctr$p2), silent = TRUE)
      if (inherits(r, "try-error") || is.null(r)) next

      key <- data.frame(contrast = ctr_name, mode = mode, x_threshold = xth,
                        unit = u, prey_family = PREY_FAMILY, stringsAsFactors = FALSE)
      res_prey_beta[[length(res_prey_beta) + 1]] <- cbind(key, r$beta)
      if (!is.null(r$per_prey)) res_prey_node[[length(res_prey_node) + 1]] <- cbind(key, r$per_prey)
    }
    msg("  -> Fait pour mode=%s, x=%s", mode, xth)
  }
}

if (length(res_prey_beta)) {
  prey_beta_tab <- bind_rows(res_prey_beta)
  prey_node_tab <- bind_rows(res_prey_node)

  write.csv(prey_beta_tab, file.path(OUT_DIR, "prey_beta_decomposition.csv"), row.names = FALSE)
  if (nrow(prey_node_tab)) write.csv(prey_node_tab, file.path(OUT_DIR, "per_prey_vulnerability.csv"), row.names = FALSE)
  msg("Succès ! Résultats sauvegardés dans : %s", OUT_DIR)
} else {
  warning("Aucune cellule n'a généré de résultat. Essayez de baisser MIN_SETS_PREY à 3.")
}
