# =============================================================================
# 6a_Engine_Trophic.R - TROPHIC TRANSITION ENGINE
# -----------------------------------------------------------------------------
# Classifies the change in each predator x size-class cell between two periods
# from three signals computed with the trawl set as the replicate unit
# (Hurlbert 1984; Pennington & Volstad 1994): Shannon diversity H' (delete-one-
# set jackknife, Welch t-test; Zahl 1977), Levins niche breadth Bs (per set,
# Welch t-test) and prey composition (Bray-Curtis, PERMANOVA and PERMDISP;
# Anderson 2001, 2006). The three significance flags give nine diagnostic
# states, collapsed into four families by R_helpers/Config_Mappings.R.
#
# One engine, three spatial levels, chosen by SPATIAL_LEVEL before source():
#   all_gulf   every trawl set in one matrix per cell (Gulf-wide test)
#   ecoregion  cells = species x size class x ecoregion
#   stratum    cells = species x size class x survey stratum
# Each level runs in its own R session (6b, 6c, 6d): the engine fixes the
# spatial columns at source() time.
#
# This file defines objects only. run_all_scenarios() (section 8) runs the five
# contrasts in both currencies and saves the results.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(vegan); library(data.table); library(sf)
})
if ("igraph" %in% (.packages())) detach("package:igraph", unload = TRUE)

# -----------------------------------------------------------------------------
# SPATIAL LEVEL
# -----------------------------------------------------------------------------
#   source  column of diet_clean holding the spatial unit (NA for all_gulf)
#   out     spatial column written to the results ("Area" or "str")
#   pool    TRUE collapses every set into one unit before testing
#   label   value written in `out` when pool = TRUE

if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

SPATIAL_LEVELS <- list(
  all_gulf = list(
    source = NA_character_, out = "Area", pool = TRUE,
    label  = "All Gulf",
    title  = "L1 all_gulf - no spatial structure, one Gulf-wide test"
  ),
  ecoregion = list(
    source = "Area", out = "Area", pool = FALSE,
    label  = NA_character_,
    title  = "L2 ecoregion - disaggregated by ecoregion"
  ),
  stratum = list(
    source = "stratum", out = "str", pool = FALSE,
    label  = NA_character_,
    title  = "L3 stratum - disaggregated by survey stratum"
  )
)

if (!SPATIAL_LEVEL %in% names(SPATIAL_LEVELS)) {
  stop("SPATIAL_LEVEL must be one of: ",
       paste(names(SPATIAL_LEVELS), collapse = ", "),
       ". Got: ", SPATIAL_LEVEL)
}

.lvl <- SPATIAL_LEVELS[[SPATIAL_LEVEL]]
SPATIAL_SOURCE <- .lvl$source   # NA for all_gulf
SPATIAL_OUT    <- .lvl$out
SPATIAL_POOL   <- .lvl$pool
SPATIAL_LABEL  <- .lvl$label
SPATIAL_TITLE  <- .lvl$title

message("Spatial level: ", SPATIAL_TITLE)

# -----------------------------------------------------------------------------
# PARAMETERS (a run script may set them before sourcing this file)
# -----------------------------------------------------------------------------
if (!exists("PREY_FAMILY"))     PREY_FAMILY     <- "_1"     # prey-grouping family swept
if (!exists("N_MIN"))           N_MIN           <- 5        # min stomachs with prey per cell and period
if (!exists("N_STRICT"))        N_STRICT        <- 25       # min stomachs per period for the reliable flag
if (!exists("MIN_SETS"))        MIN_SETS        <- 3        # min sets per period for jackknife / Welch
if (!exists("ALPHA"))           ALPHA           <- 0.05
if (!exists("R_PERM"))          R_PERM          <- 999      # permutations, PERMANOVA and PERMDISP
if (!exists("BS_TEST"))         BS_TEST         <- "welch"  # "welch", "wilcox" or "perm"
if (!exists("COMP_RELATIVE"))   COMP_RELATIVE   <- TRUE     # relative profiles before Bray-Curtis
if (!exists("OCC_BINARY"))      OCC_BINARY      <- FALSE    # TRUE: Jaccard on presence/absence for occurrence
if (!exists("ADD_NICHE_PART"))  ADD_NICHE_PART  <- TRUE     # set-level TNW / WIC / BIC
if (!exists("ADD_IND_NICHE"))   ADD_IND_NICHE   <- FALSE    # stomach-level WIC / TNW (slow)

# =============================================================================
# 1. PERIOD AND SPATIAL PREPARATION
# =============================================================================
# make_dat_classed() applies the period map and the row filters and stores the
# spatial unit under `.spatial_src`; prepare_spatial_strata() writes the output
# spatial column from it, pooled or not according to the level.

make_dat_classed <- function(diet_clean, period_map,
                             spatial_var = SPATIAL_SOURCE) {
  years_keep <- unique(unlist(period_map)); period_labels <- names(period_map)

  dat <- diet_clean %>%
    st_drop_geometry() %>%
    filter(year %in% years_keep) %>%
    filter(!is.na(somatic_length_cm), somatic_length_cm > 0) %>%
    mutate(
      period = case_when(
        year %in% period_map[[1]] ~ period_labels[1],
        year %in% period_map[[2]] ~ period_labels[2],
        TRUE ~ NA_character_),
      period = factor(period, levels = period_labels)) %>%
    filter(!is.na(period))

  if (is.na(spatial_var)) {
    # all_gulf reads no spatial information, so diet_clean needs no spatial
    # column and no stomach is dropped for lacking one.
    dat$.spatial_src <- factor(SPATIAL_LABEL)
  } else {
    if (!spatial_var %in% names(diet_clean)) {
      stop("Spatial column '", spatial_var, "' not found in diet_clean. ",
           "Required by SPATIAL_LEVEL = '", SPATIAL_LEVEL, "'.")
    }
    n_before <- nrow(dat)
    dat$.spatial_src <- factor(dat[[spatial_var]])
    dat <- dat[!is.na(dat$.spatial_src), , drop = FALSE]
    dat$.spatial_src <- droplevels(dat$.spatial_src)
    n_dropped <- n_before - nrow(dat)
    if (n_dropped > 0) {
      message(n_dropped, " row(s) dropped for a missing '", spatial_var, "'.")
    }
  }

  as.data.frame(dat)
}

# Writes the output spatial column from `.spatial_src`, collapsing it to a
# single unit at the pooled levels. run_pipeline() calls this, so one dataset
# can be handed to any level without being rebuilt.
prepare_spatial_strata <- function(dat) {
  dat <- as.data.frame(dat)
  if (!".spatial_src" %in% names(dat)) {
    stop("Column '.spatial_src' is missing. Build the dataset with ",
         "make_dat_classed() under the same SPATIAL_LEVEL.")
  }
  dat[[SPATIAL_OUT]] <- if (isTRUE(SPATIAL_POOL)) {
    factor(SPATIAL_LABEL)
  } else {
    droplevels(dat$.spatial_src)
  }
  dat
}

# =============================================================================
# 2. SET x PREY MATRICES
# =============================================================================
# Stomach contents are summed to the trawl set: biomass = wet mass per prey
# category and set; occurrence = share of the set's stomachs containing the
# category.

.add_set_uid <- function(d) {
  d$set_uid <- paste(d$year, d[["vessel.code"]], d$set, sep = "_")
  d
}

# Fast subsetting using data.table semantics
.filter_cell <- function(d, sp = NULL, sz = NULL, ar = NULL, per = NULL) {
  if (!data.table::is.data.table(d)) d <- data.table::as.data.table(d)
  if (!is.null(sp))  d <- d[predator_species_common_name == sp]
  if (!is.null(sz))  d <- d[as.character(size_class) == sz]
  if (!is.null(ar))  d <- d[as.character(get(SPATIAL_OUT)) == ar]
  if (!is.null(per)) d <- d[as.character(period) == per]
  d
}

# Build Set x Prey matrix (Biomass)
build_biomass_set <- function(data, sp = NULL, sz = NULL, ar = NULL, per = NULL, col_prey) {
  d <- .filter_cell(data, sp, sz, ar, per)
  d <- d[!is.na(d[[col_prey]]) & !is.na(set)]
  if (nrow(d) == 0) return(NULL)

  d <- .add_set_uid(d)
  mat <- d %>%
    group_by(set_uid, prey = .data[[col_prey]]) %>%
    summarise(w = sum(somatic_wt_g, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = prey, values_from = w, values_fill = 0) %>%
    column_to_rownames("set_uid")

  mat <- as.matrix(mat)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)
  mat
}

# Build Set x Prey matrix (Frequency of Occurrence)
build_occ_set <- function(data, sp = NULL, sz = NULL, ar = NULL, per = NULL, col_prey) {
  d <- .filter_cell(data, sp, sz, ar, per)
  d <- d[!is.na(d[[col_prey]]) & !is.na(set) & !is.na(stomach_id)]
  if (nrow(d) == 0) return(NULL)

  d <- .add_set_uid(d)
  denom <- d %>% distinct(set_uid, stomach_id) %>% count(set_uid, name = "n_stomachs")

  mat <- d %>%
    distinct(set_uid, stomach_id, prey = .data[[col_prey]]) %>%
    count(set_uid, prey, name = "freq") %>%
    left_join(denom, by = "set_uid") %>%
    mutate(prop_occ = freq / n_stomachs) %>%
    select(set_uid, prey, prop_occ) %>%
    pivot_wider(names_from = prey, values_from = prop_occ, values_fill = 0) %>%
    column_to_rownames("set_uid")

  mat <- as.matrix(mat)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)
  mat
}

# Build Individual Stomach Matrix (for Optional Niche Partitioning)
build_biomass_stomach <- function(data, sp = NULL, sz = NULL, ar = NULL, per = NULL, col_prey) {
  d <- .filter_cell(data, sp, sz, ar, per)
  d <- d[!is.na(d[[col_prey]]) & !is.na(stomach_id) & !is.na(set)]
  if (nrow(d) == 0 || dplyr::n_distinct(d$stomach_id) < N_MIN) return(NULL)

  d <- .add_set_uid(d)
  key <- d %>% distinct(stomach_id, set_uid)

  mat <- d %>%
    group_by(stomach_id, prey = .data[[col_prey]]) %>%
    summarise(w = sum(somatic_wt_g, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = prey, values_from = w, values_fill = 0) %>%
    column_to_rownames("stomach_id")

  mat <- as.matrix(mat)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)

  attr(mat, "set_uid") <- key$set_uid[match(rownames(mat), key$stomach_id)]
  mat
}

# Matrix Alignment
align_mats <- function(m1, m2) {
  if (is.null(m1) || is.null(m2)) return(list(m1 = NULL, m2 = NULL))

  s1 <- attr(m1, "set_uid"); s2 <- attr(m2, "set_uid")
  m1 <- as.matrix(m1); m2 <- as.matrix(m2)
  all_cols <- union(colnames(m1), colnames(m2))

  add_missing <- function(m) {
    rn <- rownames(m)
    miss <- setdiff(all_cols, colnames(m))
    if (length(miss) > 0) {
      pad <- matrix(0, nrow(m), length(miss), dimnames = list(rn, miss))
      m <- cbind(m, pad)
    }
    m <- m[, all_cols, drop = FALSE]
    rownames(m) <- rn
    m
  }

  o1 <- add_missing(m1); o2 <- add_missing(m2)
  if (!is.null(s1)) attr(o1, "set_uid") <- s1
  if (!is.null(s2)) attr(o2, "set_uid") <- s2
  list(m1 = o1, m2 = o2)
}

# Guard: `somatic_wt_g` must hold the PREY weight (it varies among the prey
# rows of one stomach and does not scale with predator length), otherwise the
# biomass currency would sum predator mass. Warns, never stops.

check_prey_weight_column <- function(dat, wt_col = "somatic_wt_g",
                                     len_col = "somatic_length_cm",
                                     id_col = "stomach_id") {
  if (!all(c(wt_col, id_col) %in% names(dat))) {
    message("check_prey_weight_column: '", wt_col, "' or '", id_col,
            "' absent; check skipped.")
    return(invisible(NULL))
  }
  d <- as.data.frame(dat)[, intersect(c(wt_col, len_col, id_col), names(dat)), drop = FALSE]
  d <- d[stats::complete.cases(d[, c(wt_col, id_col)]), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))

  # (1) constant within a stomach?
  n_rows <- table(d[[id_col]])
  multi  <- names(n_rows)[n_rows > 1]
  const  <- NA_real_
  if (length(multi)) {
    sub <- d[d[[id_col]] %in% multi, , drop = FALSE]
    nd  <- tapply(sub[[wt_col]], sub[[id_col]], function(x) length(unique(x)))
    const <- mean(nd == 1)
  }

  # (2) allometric with predator length?
  slope <- r2 <- NA_real_
  if (len_col %in% names(d)) {
    ok <- is.finite(d[[wt_col]]) & d[[wt_col]] > 0 &
          is.finite(d[[len_col]]) & d[[len_col]] > 0
    if (sum(ok) > 30) {
      fit <- stats::lm(log(d[[wt_col]][ok]) ~ log(d[[len_col]][ok]))
      slope <- unname(stats::coef(fit)[2]); r2 <- summary(fit)$r.squared
    }
  }

  cat(sprintf(
    "Weight-column check on '%s': constant within stomach in %s of multi-prey stomachs; log-log slope vs length = %s (R2 = %s)\n",
    wt_col,
    if (is.na(const)) "n/a" else paste0(round(100 * const, 1), "%"),
    if (is.na(slope)) "n/a" else round(slope, 2),
    if (is.na(r2))    "n/a" else round(r2, 2)))

  looks_predator <- (!is.na(const) && const > 0.95) ||
                    (!is.na(slope) && !is.na(r2) && slope > 2.3 && r2 > 0.7)
  if (isTRUE(looks_predator)) {
    warning("'", wt_col, "' now behaves like a PREDATOR weight. It held the ",
            "PREY weight when this pipeline was validated, so the upstream ",
            "build has changed. The biomass currency would sum predator body ",
            "mass per prey category and would not be a prey biomass. Stop and ",
            "check the build before running the sweep.", call. = FALSE)
  }
  invisible(list(constant_within_stomach = const, allometric_slope = slope, r2 = r2))
}

# =============================================================================
# 3. DIVERSITY AND NICHE BREADTH
# =============================================================================

shannon_from_vec <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log(p))
}

levins_bs_from_vec <- function(x, n_cat) {
  x <- as.numeric(x); x <- x[is.finite(x) & x > 0]
  if (length(x) == 0 || n_cat < 2) return(NA_real_)
  if (length(x) < 2) return(0)
  p <- x / sum(x)
  ((1 / sum(p^2)) - 1) / (n_cat - 1)
}

row_indices_set <- function(mat, n_cat = ncol(mat)) {
  if (is.null(mat) || nrow(mat) == 0) {
    return(list(H = numeric(0), Bs = numeric(0)))
  }
  list(
    H  = apply(mat, 1, shannon_from_vec),
    Bs = apply(mat, 1, levins_bs_from_vec, n_cat = n_cat)
  )
}

shannon_jackknife_set <- function(mat_set1, mat_set2, min_sets = MIN_SETS) {
  H_from_mat <- function(m) {
    s <- sum(m)
    if (!is.finite(s) || s <= 0) return(NA_real_)
    p <- colSums(m) / s; p <- p[p > 0]
    if (length(p) == 0) return(NA_real_)
    -sum(p * log(p))
  }

  jackknife_H <- function(m) {
    n <- nrow(m)
    if (n < min_sets) return(list(H = NA_real_, var = NA_real_, n = n, ok = FALSE))
    H_full <- H_from_mat(m)
    if (!is.finite(H_full)) return(list(H = NA_real_, var = NA_real_, n = n, ok = FALSE))

    H_i <- vapply(seq_len(n), function(i) H_from_mat(m[-i, , drop = FALSE]), numeric(1))
    if (any(!is.finite(H_i))) return(list(H = NA_real_, var = NA_real_, n = n, ok = FALSE))

    pseudo <- n * H_full - (n - 1) * H_i
    list(H = mean(pseudo), var = stats::var(pseudo) / n, n = n, H_plugin = H_full, ok = TRUE)
  }

  j1 <- jackknife_H(mat_set1); j2 <- jackknife_H(mat_set2)
  base <- list(
    H1 = if (isTRUE(j1$ok)) round(j1$H_plugin, 3) else NA_real_,
    H2 = if (isTRUE(j2$ok)) round(j2$H_plugin, 3) else NA_real_,
    H1_jack = round(j1$H, 3), H2_jack = round(j2$H, 3),
    delta = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_,
    n1_sets = j1$n, n2_sets = j2$n
  )

  if (!isTRUE(j1$ok) || !isTRUE(j2$ok) || !is.finite(j1$var) || !is.finite(j2$var) || (j1$var + j2$var) <= 0) {
    return(base)
  }

  t_val <- (j2$H - j1$H) / sqrt(j1$var + j2$var)
  df <- (j1$var + j2$var)^2 / (j1$var^2 / (j1$n - 1) + j2$var^2 / (j2$n - 1))

  base$delta <- round(j2$H - j1$H, 3); base$t <- round(t_val, 3); base$df <- round(df, 1)
  base$p <- round(2 * stats::pt(-abs(t_val), df = max(df, 1)), 4)
  base
}

compare_index_sets <- function(v1, v2, test = BS_TEST, n_perm = R_PERM) {
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  out <- list(mean_P1 = NA_real_, mean_P2 = NA_real_, med_P1 = NA_real_, med_P2 = NA_real_, delta = NA_real_, p = NA_real_, test = test, n1 = length(v1), n2 = length(v2))
  if (length(v1) < 2 || length(v2) < 2) return(out)

  out$mean_P1 <- round(mean(v1), 3); out$mean_P2 <- round(mean(v2), 3)
  out$med_P1  <- round(stats::median(v1), 3); out$med_P2  <- round(stats::median(v2), 3)
  out$delta   <- round(mean(v2) - mean(v1), 3)

  p <- switch(
    test,
    welch = {
      if (stats::var(v1) == 0 && stats::var(v2) == 0) {
        if (isTRUE(all.equal(mean(v1), mean(v2)))) 1 else 0
      } else {
        suppressWarnings(stats::t.test(v2, v1, var.equal = FALSE)$p.value)
      }
    },
    wilcox = suppressWarnings(stats::wilcox.test(v1, v2, exact = FALSE)$p.value),
    perm = {
      allv <- c(v1, v2); g <- rep(c(1L, 2L), c(length(v1), length(v2)))
      obs <- mean(allv[g == 2L]) - mean(allv[g == 1L])
      null <- replicate(n_perm, {
        gg <- sample(g); mean(allv[gg == 2L]) - mean(allv[gg == 1L])
      })
      (sum(abs(null) >= abs(obs)) + 1) / (n_perm + 1)
    },
    stop("Unknown BS_TEST: ", test)
  )

  out$p <- round(as.numeric(p), 4)
  out
}

# =============================================================================
# 4. COMPOSITION
# =============================================================================
# PERMANOVA (adonis2) on Bray-Curtis dissimilarities between set profiles, with
# PERMDISP (betadisper) to flag cells where location and dispersion effects
# are confounded. Failures are counted so that a systematic error is visible.
.comp_fail_count <- 0L

run_composition_set <- function(s1, s2, mode = c("biomass", "occurrence"), perm = R_PERM, relative = COMP_RELATIVE) {
  mode <- match.arg(mode)
  empty <- list(BC = NA, R2 = NA, p_comp = NA, p_disp = NA, comp_dispersion = NA, n_set_P1 = 0L, n_set_P2 = 0L)

  al <- align_mats(s1, s2)
  if (is.null(al$m1) || is.null(al$m2)) return(empty)

  empty$n_set_P1 <- nrow(al$m1); empty$n_set_P2 <- nrow(al$m2)
  if (nrow(al$m1) < 2 || nrow(al$m2) < 2) return(empty)

  m1 <- al$m1; m2 <- al$m2
  combined <- rbind(m1, m2)
  groups <- factor(c(rep("P1", nrow(m1)), rep("P2", nrow(m2))))

  use_jaccard <- (mode == "occurrence") && isTRUE(OCC_BINARY)
  d <- if (use_jaccard) {
    vegdist(combined, method = "jaccard", binary = TRUE)
  } else {
    vegdist(if (relative) decostand(combined, "total") else combined, method = "bray")
  }

  tryCatch({
    pr <- adonis2(d ~ groups, permutations = perm)
    pmt <- permutest(suppressWarnings(betadisper(d, groups)), permutations = perm)

    if (use_jaccard) {
      c1 <- as.integer(colSums(m1) > 0); c2 <- as.integer(colSums(m2) > 0)
      BC <- as.numeric(vegdist(rbind(c1, c2), method = "jaccard", binary = TRUE))
    } else {
      c1 <- colSums(if (relative) decostand(m1, "total") else m1); c1 <- c1 / sum(c1)
      c2 <- colSums(if (relative) decostand(m2, "total") else m2); c2 <- c2 / sum(c2)
      BC <- as.numeric(vegdist(rbind(c1, c2), method = "bray"))
    }

    list(
      BC = round(BC, 3), R2 = round(pr$R2[1], 3),
      p_comp = round(pr$`Pr(>F)`[1], 4), p_disp = round(pmt$tab$`Pr(>F)`[1], 4),
      comp_dispersion = (!is.na(pr$`Pr(>F)`[1]) && pr$`Pr(>F)`[1] < ALPHA && !is.na(pmt$tab$`Pr(>F)`[1]) && pmt$tab$`Pr(>F)`[1] < ALPHA),
      n_set_P1 = nrow(m1), n_set_P2 = nrow(m2)
    )
  }, error = function(e) {
    .comp_fail_count <<- .comp_fail_count + 1L
    if (.comp_fail_count <= 5L) {
      message("run_composition_set failed (", .comp_fail_count, "): ",
              conditionMessage(e))
    }
    empty
  })
}

# =============================================================================
# 5. NICHE WIDTH PARTITION
# =============================================================================
# Shannon entropy decomposes additively over a nested grouping.
#   niche_partition_set(): unit = trawl set. TNW = niche width of the period,
#     WIC_set = mean width of one tow, BIC_set = among-tow component. This is
#     not the individual WIC of Roughgarden (1972) / Bolnick et al. (2002).
#   niche_partition_ind(): unit = stomach (Roughgarden 1972; Bolnick et al.
#     2002), with BIC split within and among sets. Biomass only.

niche_partition_set <- function(mat) {
  if (is.null(mat) || nrow(mat) == 0 || sum(mat) <= 0) return(NULL)
  w <- rowSums(mat)
  mat <- mat[w > 0, , drop = FALSE]; w <- w[w > 0]
  if (!nrow(mat)) return(NULL)

  p_i <- w / sum(w); q <- colSums(mat) / sum(mat)
  TNW <- shannon_from_vec(q)
  WIC <- sum(p_i * apply(mat, 1, shannon_from_vec), na.rm = TRUE)

  if (!is.finite(TNW) || TNW <= 0) return(NULL)
  list(TNW = round(TNW, 4), WIC_set = round(WIC, 4), BIC_set = round(TNW - WIC, 4), WIC_TNW_set = round(WIC / TNW, 4))
}

niche_partition_ind <- function(mat, set_uid) {
  if (is.null(mat) || nrow(mat) == 0 || sum(mat) <= 0) return(NULL)
  w <- rowSums(mat)
  keep <- w > 0
  mat <- mat[keep, , drop = FALSE]; w <- w[keep]
  set_uid <- as.character(set_uid)[keep]
  if (!nrow(mat)) return(NULL)

  p_i <- w / sum(w); q <- colSums(mat) / sum(mat)
  TNW <- shannon_from_vec(q)
  WIC <- sum(p_i * apply(mat, 1, shannon_from_vec), na.rm = TRUE)
  if (!is.finite(TNW) || TNW <= 0) return(NULL)

  sets <- unique(set_uid)
  w_s <- vapply(sets, function(s) sum(w[set_uid == s]), numeric(1))
  H_s <- vapply(sets, function(s) {
    m <- mat[set_uid == s, , drop = FALSE]
    shannon_from_vec(colSums(m))
  }, numeric(1))
  mean_H_set <- sum((w_s / sum(w_s)) * H_s, na.rm = TRUE)

  list(TNW_ind = round(TNW, 4), WIC_ind = round(WIC, 4), WIC_TNW_ind = round(WIC / TNW, 4),
       BIC_within_set = round((mean_H_set - WIC) / TNW, 4), BIC_among_set  = round((TNW - mean_H_set) / TNW, 4),
       n_sets_ind = length(sets))
}

design_effect <- function(y, set_uid) {
  d <- data.frame(y = as.numeric(y), g = factor(set_uid))
  d <- d[is.finite(d$y), ]
  d$g <- droplevels(d$g)
  if (nlevels(d$g) < 2 || nrow(d) <= nlevels(d$g)) return(NULL)

  tab <- summary(stats::aov(y ~ g, data = d))[[1]]
  ms_b <- tab[1, "Mean Sq"]; ms_w <- tab[2, "Mean Sq"]
  k <- nlevels(d$g); n <- nrow(d); ni <- as.numeric(table(d$g))
  m0 <- (n - sum(ni^2) / n) / (k - 1)
  var_b <- max((ms_b - ms_w) / m0, 0)
  icc <- if ((var_b + ms_w) > 0) var_b / (var_b + ms_w) else 0
  deff <- 1 + (mean(ni) - 1) * icc

  list(icc = round(icc, 3), deff = round(deff, 2), n_eff = round(n / deff, 1))
}

# =============================================================================
# 6. MAIN PIPELINE AND CLASSIFICATION
# =============================================================================

run_pipeline <- function(dat, mode = c("biomass", "occurrence"), period_1, period_2, seed = 1234) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  .comp_fail_count <<- 0L

  dat <- prepare_spatial_strata(dat)
  setDT(dat)

  prey_cols <- names(dat)[grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"), names(dat))]
  out <- list(); inventory_all <- list()

  for (target_col in prey_cols) {
    curr_x <- as.numeric(str_extract(target_col, "\\d+"))
    cat(sprintf("\n--- X = %d (%s) --- %s ", curr_x,
                format(Sys.time(), "%H:%M:%S"),
                paste0(which(prey_cols == target_col), "/", length(prey_cols))))

    inventory_x <- dat %>%
      filter(period %in% c(period_1, period_2), !is.na(.data[[target_col]]), !is.na(size_class)) %>%
      group_by(across(all_of(c("predator_species_common_name", SPATIAL_OUT, "size_class", "period")))) %>%
      summarise(n = n_distinct(stomach_id), .groups = "drop") %>%
      pivot_wider(names_from = period, values_from = n, values_fill = 0)

    if (!(period_1 %in% names(inventory_x)) || !(period_2 %in% names(inventory_x))) next

    inventory_x <- inventory_x %>%
      rename(n_P1 = all_of(period_1), n_P2 = all_of(period_2)) %>%
      filter(n_P1 >= N_MIN, n_P2 >= N_MIN)

    if (nrow(inventory_x) == 0) next
    inventory_x <- inventory_x %>% mutate(x_threshold = curr_x)
    inventory_all[[target_col]] <- inventory_x
    res_x <- vector("list", nrow(inventory_x))

    for (i in seq_len(nrow(inventory_x))) {

      # Seed per cell, so each cell draws the same permutations whatever ran before.
      if (!is.null(seed)) set.seed(seed + 1000L * curr_x + i)
      sp <- inventory_x$predator_species_common_name[i]
      sz <- as.character(inventory_x$size_class[i])
      ar <- as.character(inventory_x[[SPATIAL_OUT]][i])

      n_sto_P1 <- inventory_x$n_P1[i]; n_sto_P2 <- inventory_x$n_P2[i]

      if (mode == "biomass") {
        s1 <- build_biomass_set(dat, sp, sz, ar, period_1, target_col)
        s2 <- build_biomass_set(dat, sp, sz, ar, period_2, target_col)
      } else {
        s1 <- build_occ_set(dat, sp, sz, ar, period_1, target_col)
        s2 <- build_occ_set(dat, sp, sz, ar, period_2, target_col)
      }

      if (is.null(s1) || is.null(s2)) next
      al_set <- align_mats(s1, s2)
      s1a <- al_set$m1; s2a <- al_set$m2
      if (is.null(s1a) || is.null(s2a)) next

      n_cat <- ncol(s1a)
      htch <- shannon_jackknife_set(s1a, s2a)
      idx1 <- row_indices_set(s1a, n_cat); idx2 <- row_indices_set(s2a, n_cat)
      Bscmp <- compare_index_sets(idx1$Bs, idx2$Bs)
      comp <- run_composition_set(s1a, s2a, mode = mode, perm = R_PERM)

      np1 <- np2 <- NULL
      if (isTRUE(ADD_NICHE_PART)) {
        np1 <- niche_partition_set(s1a); np2 <- niche_partition_set(s2a)
      }

      ind1 <- ind2 <- de <- NULL
      if (isTRUE(ADD_IND_NICHE) && mode == "biomass") {
        b1 <- build_biomass_stomach(dat, sp, sz, ar, period_1, target_col)
        b2 <- build_biomass_stomach(dat, sp, sz, ar, period_2, target_col)
        if (!is.null(b1) && !is.null(b2)) {
          alb <- align_mats(b1, b2)
          ind1 <- niche_partition_ind(alb$m1, attr(alb$m1, "set_uid"))
          ind2 <- niche_partition_ind(alb$m2, attr(alb$m2, "set_uid"))
          de <- design_effect(
            c(apply(alb$m1, 1, shannon_from_vec), apply(alb$m2, 1, shannon_from_vec)),
            c(attr(alb$m1, "set_uid"), attr(alb$m2, "set_uid"))
          )
        }
      }

      H_sig      <- !is.na(htch$p) && htch$p < ALPHA
      Bs_sig     <- !is.na(Bscmp$p) && Bscmp$p < ALPHA
      comp_sig   <- !is.na(comp$p_comp) && comp$p_comp < ALPHA
      comp_trend <- !is.na(comp$p_comp) && comp$p_comp < 0.10

      diagnostic <- dplyr::case_when(
        !H_sig & !Bs_sig & !comp_sig & !comp_trend ~ "Stable Diet",
        !H_sig & !Bs_sig & !comp_sig &  comp_trend ~ "Emerging Shift",
        !H_sig & !Bs_sig &  comp_sig               ~ "Ghost Shift",
        H_sig & !Bs_sig & !comp_sig               ~ "Internal Rebalancing",
        !H_sig &  Bs_sig & !comp_sig               ~ "Niche Compression/Expansion",
        H_sig &  Bs_sig & !comp_sig               ~ "Niche Restructuring",
        H_sig & !Bs_sig &  comp_sig               ~ "Partial Diet Shift",
        !H_sig &  Bs_sig &  comp_sig               ~ "Structural Shift",
        H_sig &  Bs_sig &  comp_sig               ~ "Major Shift",
        .default = "Inconclusive"
      )

      testable <- !is.na(htch$p) && !is.na(Bscmp$p) && !is.na(comp$p_comp)
      if (!testable) diagnostic <- "Inconclusive"

      reliable <- (
        n_sto_P1 >= N_STRICT && n_sto_P2 >= N_STRICT &&
          comp$n_set_P1 >= MIN_SETS && comp$n_set_P2 >= MIN_SETS &&
          htch$n1_sets >= MIN_SETS && htch$n2_sets >= MIN_SETS &&
          Bscmp$n1 >= MIN_SETS && Bscmp$n2 >= MIN_SETS
      )

      gv <- function(o, f) if (is.null(o) || is.null(o[[f]])) NA_real_ else o[[f]]

      row <- tibble(
        x_threshold = curr_x,prey_family = PREY_FAMILY, species = sp, size_class = sz,
        period_1 = period_1, period_2 = period_2, spatial_level = SPATIAL_LEVEL,
        H_P1 = htch$H1, H_P2 = htch$H2, H_P1_jack = htch$H1_jack, H_P2_jack = htch$H2_jack,
        delta_H = htch$delta, t_H = htch$t, df_H = htch$df, p_H = htch$p,
        Bs_P1 = Bscmp$mean_P1, Bs_P2 = Bscmp$mean_P2, Bs_med_P1 = Bscmp$med_P1, Bs_med_P2 = Bscmp$med_P2,
        delta_Bs = Bscmp$delta, p_Bs = Bscmp$p, Bs_test = Bscmp$test, n_cat = n_cat,
        p_comp = comp$p_comp, R2_comp = comp$R2, BC = comp$BC,
        p_disp = comp$p_disp, comp_dispersion = comp$comp_dispersion,
        TNW_P1 = gv(np1, "TNW"), TNW_P2 = gv(np2, "TNW"),
        WIC_TNW_set_P1 = gv(np1, "WIC_TNW_set"), WIC_TNW_set_P2 = gv(np2, "WIC_TNW_set"),
        WIC_TNW_ind_P1 = gv(ind1, "WIC_TNW_ind"), WIC_TNW_ind_P2 = gv(ind2, "WIC_TNW_ind"),
        BIC_within_set_P1 = gv(ind1, "BIC_within_set"), BIC_within_set_P2 = gv(ind2, "BIC_within_set"),
        BIC_among_set_P1 = gv(ind1, "BIC_among_set"), BIC_among_set_P2 = gv(ind2, "BIC_among_set"),
        icc = gv(de, "icc"), deff = gv(de, "deff"), n_eff_sto = gv(de, "n_eff"),
        n_sto_P1 = n_sto_P1, n_sto_P2 = n_sto_P2,
        n_set_P1 = comp$n_set_P1, n_set_P2 = comp$n_set_P2,
        testable = testable, reliable = reliable,
        diagnostic = diagnostic, mode = mode
      )
      row[[SPATIAL_OUT]] <- ar
      res_x[[i]] <- row
      if (i %% 10 == 0) cat(".")
    }
    out[[target_col]] <- bind_rows(res_x)
  }
  res <- bind_rows(out)

  if (.comp_fail_count > 0L) {
    warning(.comp_fail_count, " composition test(s) failed and returned NA. ",
            "A count close to the number of cells points to a systematic cause ",
            "(package version, data structure) rather than sparse cells.",
            call. = FALSE)
  }
  cat(sprintf("\n[%s] %d cells | %.1f%% testable | %d composition failures\n",
              mode, nrow(res), 100 * mean(res$testable), .comp_fail_count))

  list(results = res, inventory = bind_rows(inventory_all),
       n_comp_failures = .comp_fail_count)
}

# =============================================================================
# 7. PER-RUN FIGURE AND SAVE
# =============================================================================

diet_diagnostics_table <- tribble(
  ~Diagnostic, ~Statistical_Logic, ~Ecological_Interpretation,
  "Stable Diet", "Shannon: ns | Niche breadth: ns | Composition: ns", "The diet is stable in prey consumed, their proportions, and overall resource breadth.",
  "Emerging Shift", "Shannon: ns | Niche breadth: ns | Composition: 0.05-0.10", "Early signals of taxonomic change; the prey list is beginning to deviate from baseline.",
  "Ghost Shift", "Shannon: ns | Niche breadth: ns | Composition: significant", "Species have been replaced by others occupying a similar functional space.",
  "Niche Compression/Expansion", "Shannon: ns | Niche breadth: significant | Composition: ns", "The range of prey consumed has changed, even though composition remains stable.",
  "Internal Rebalancing", "Shannon: significant | Niche breadth: ns | Composition: ns", "Relative abundance of prey has shifted within the same prey set.",
  "Niche Restructuring", "Shannon: significant | Niche breadth: significant | Composition: ns", "A deep reorganization of diet structure and breadth.",
  "Partial Diet Shift", "Shannon: significant | Niche breadth: ns | Composition: significant", "Species replacement accompanied by a change in diet diversity.",
  "Structural Shift", "Shannon: ns | Niche breadth: significant | Composition: significant", "Taxonomic turnover coupled with a change in foraging range.",
  "Major Shift", "Shannon: significant | Niche breadth: significant | Composition: significant", "A systemic overhaul: prey identity, diversity, and specialization all change.",
  "Inconclusive", "Sample size or data consistency insufficient", "No diagnostic can be reliably assigned."
)

make_diag_plot2 <- function(df, mode_label, period_1, period_2,
                            spatial_level = SPATIAL_LEVEL,
                            by_area = FALSE, file_out = NULL,
                            title_main = NULL, subtitle_main = NULL,
                            y_lab = NULL, inventory_df = NULL,
                            drop_untestable = TRUE) {

  if (drop_untestable && "testable" %in% names(df)) df <- dplyr::filter(df, testable)
  if (nrow(df) == 0) return(NULL)

  # by_area = facet by spatial unit; meaningless with a single unit.
  if (by_area && dplyr::n_distinct(df[[SPATIAL_OUT]]) < 2) by_area <- FALSE

  STABILITY_LABS <- c("Stable Diet", "Emerging Shift")
  TURNOVER_LABS  <- c("Major Shift", "Structural Shift", "Partial Diet Shift")

  if (!by_area) {
    n_eff <- dplyr::n_distinct(paste(df$species, df$size_class, df[[SPATIAL_OUT]]))

    df2 <- df %>%
      group_by(diagnostic) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        total = sum(count),
        percentage = 100 * count / total,
        std_error_pct = sqrt((percentage / 100) * (1 - percentage / 100) / n_eff) * 100
      ) %>%
      rename(Diagnostic = diagnostic) %>%
      left_join(diet_diagnostics_table, by = "Diagnostic")

    composites <- df2 %>%
      summarise(
        stability_pct = sum(percentage[Diagnostic %in% STABILITY_LABS]),
        turnover_func_pct = sum(percentage[Diagnostic %in% TURNOVER_LABS])
      )

    ymax <- max(df2$percentage, na.rm = TRUE)
    pct_x <- ymax + 4
    desc_x <- ymax + 22

    plot_data <- df2 %>%
      mutate(
        Stats_Label = paste0(round(percentage, 2), "% (±", round(std_error_pct, 2), ")"),
        Interpretation_Wrapped = stringr::str_wrap(Ecological_Interpretation, width = 40),
        Diagnostic = factor(Diagnostic, levels = rev(unique(Diagnostic)))
      )

    p <- ggplot(plot_data, aes(x = Diagnostic, y = percentage, fill = Diagnostic)) +
      geom_col(alpha = 0.85, width = 0.7) +
      geom_errorbar(
        aes(ymin = pmax(percentage - std_error_pct, 0),
            ymax = percentage + std_error_pct),
        width = 0.2
      ) +
      geom_text(aes(y = pct_x, label = Stats_Label),
                hjust = 0, size = 3.4, fontface = "bold", color = "black") +
      geom_text(aes(y = desc_x, label = Interpretation_Wrapped),
                hjust = 0, size = 3, lineheight = 0.9, color = "grey35") +
      coord_flip(clip = "off") +
      scale_y_continuous(limits = c(0, desc_x + 44),
                         breaks = seq(0, ceiling(ymax / 20) * 20, by = 20)) +
      scale_fill_brewer(palette = "Set3") +
      labs(
        title = title_main,
        subtitle = paste0(
          "Stability: ", round(composites$stability_pct, 2),
          "% | Turnover with functional change: ",
          round(composites$turnover_func_pct, 2), "%"
        ),
        x = NULL, y = y_lab,
        caption = paste0(
          "Source: trophic transition analysis (", mode_label, ", ",
          period_1, " vs ", period_2, ", ", spatial_level,
          "). Error bars: cluster-corrected SE, n = ", n_eff, " independent cells."
        )
      ) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.y = element_text(face = "bold"),
            plot.margin = margin(10, 25, 10, 10))

  } else {
    if (is.null(inventory_df)) stop("inventory_df is required when by_area = TRUE")

    unit_summary <- inventory_df %>%
      group_by(across(all_of(SPATIAL_OUT))) %>%
      summarise(n_total = sum(n_P1 + n_P2) / n_distinct(x_threshold), .groups = "drop")

    n_eff_unit <- df %>%
      group_by(across(all_of(SPATIAL_OUT))) %>%
      summarise(n_eff = dplyr::n_distinct(paste(species, size_class)), .groups = "drop")

    stats_by_unit <- df %>%
      group_by(across(all_of(c(SPATIAL_OUT, "diagnostic")))) %>%
      summarise(count = n(), .groups = "drop") %>%
      left_join(n_eff_unit, by = SPATIAL_OUT) %>%
      group_by(across(all_of(SPATIAL_OUT))) %>%
      mutate(
        total_area = sum(count),
        percentage = 100 * count / total_area,
        std_error_pct = sqrt((percentage / 100) * (1 - percentage / 100) / n_eff) * 100,
        is_stable_group = diagnostic %in% STABILITY_LABS,
        is_turnover_func = diagnostic %in% TURNOVER_LABS,
        stability_score = sum(percentage[is_stable_group]),
        turnover_func_score = sum(percentage[is_turnover_func]),
        Unit_Label = paste0(
          .data[[SPATIAL_OUT]],
          "\n(Stability: ", round(stability_score, 2),
          "% | Turnover with functional change: ",
          round(turnover_func_score, 2), "%)"
        )
      ) %>%
      ungroup() %>%
      rename(Diagnostic = diagnostic) %>%
      left_join(unit_summary, by = SPATIAL_OUT)

    ann_unit <- stats_by_unit %>%
      distinct(across(all_of(SPATIAL_OUT)), Unit_Label, n_total) %>%
      mutate(x = 0.5, y = 50, label = paste0("n = ", round(n_total)))

    plot_data <- stats_by_unit %>%
      left_join(diet_diagnostics_table, by = "Diagnostic") %>%
      mutate(
        Full_Label = paste0(round(percentage, 2), "% (±", round(std_error_pct, 2), ")"),
        Diagnostic = factor(Diagnostic, levels = rev(unique(Diagnostic)))
      )

    ymax_u <- max(plot_data$percentage, na.rm = TRUE)

    p <- ggplot(plot_data, aes(x = Diagnostic, y = percentage, fill = Diagnostic)) +
      geom_col(alpha = 0.85, width = 0.7) +
      geom_errorbar(
        aes(ymin = pmax(percentage - std_error_pct, 0),
            ymax = percentage + std_error_pct),
        width = 0.2, color = "black"
      ) +
      geom_text(aes(y = ymax_u + 6, label = Full_Label),
                hjust = 0, size = 2.8, fontface = "bold", color = "black") +
      geom_text(data = ann_unit, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, color = "black", size = 3.5, fontface = "bold") +
      coord_flip(clip = "off") +
      facet_wrap(~Unit_Label, scales = "free_y") +
      scale_y_continuous(limits = c(0, ymax_u + 22)) +
      scale_fill_brewer(palette = "Set3") +
      labs(
        title = title_main, subtitle = subtitle_main, x = NULL, y = y_lab,
        caption = paste0(
          "Source: trophic transition analysis (", mode_label, ", ",
          period_1, " vs ", period_2, ", ", spatial_level,
          "). Error bars: cluster-corrected SE on independent cells."
        )
      ) +
      theme_minimal() +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold", size = 9, color = "#1a3a5c"),
            strip.background = element_rect(fill = "#f0f0f0", color = NA),
            panel.spacing.x = unit(3, "lines"),
            plot.margin = margin(10, 60, 10, 10))
  }

  if (!is.null(file_out)) ggsave(file_out, plot = p, width = 12, height = 8, dpi = 300)
  p
}

run_save_plot <- function(dat, sc, mode, period_1, period_2,
                          drop_untestable = TRUE,
                          out_dir  = paste0("Sensitivity_", SPATIAL_LEVEL),
                          plot_dir = paste0("Sensitivity_Plot_", SPATIAL_LEVEL)) {

  dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  label <- paste0(
    sc,PREY_FAMILY, "_", SPATIAL_LEVEL, "_",
    gsub("[^A-Za-z0-9]+", "_", period_1), "_vs_",
    gsub("[^A-Za-z0-9]+", "_", period_2)
  )

  res <- run_pipeline(dat = dat, mode = mode, period_1 = period_1, period_2 = period_2)

  save(res, file = file.path(out_dir, paste0(mode, "_", label, ".rda")))

  make_diag_plot2(
    df = res$results, inventory_df = res$inventory, mode_label = mode,
    period_1 = period_1, period_2 = period_2,
    by_area = FALSE,
    file_out = file.path(plot_dir, paste0(mode, "_summary_", label, ".png")),
    title_main = paste0(
      ifelse(mode == "biomass", "Biomass", "Occurrence"),
      "-based dietary transitions: ", period_1, " vs ", period_2
    ),
    y_lab = "Frequency (%)"
  )

  # Per-unit figure, only where more than one spatial unit exists.
  if (!isTRUE(SPATIAL_POOL)) {
    make_diag_plot2(
      df = res$results, inventory_df = res$inventory, mode_label = mode,
      period_1 = period_1, period_2 = period_2,
      by_area = TRUE,
      file_out = file.path(plot_dir, paste0(mode, "_by_unit_", label, ".png")),
      title_main = paste0(
        ifelse(mode == "biomass", "Biomass", "Occurrence"),
        "-based dietary transitions by spatial unit: ", period_1, " vs ", period_2
      ),
      subtitle_main = "Stable + Emerging vs Major + Partial + Structural shifts",
      y_lab = "Within-unit frequency (%)"
    )
  }

  invisible(res)
}

# =============================================================================
# 8. ALL SCENARIOS FOR THE CURRENT LEVEL
# =============================================================================
# Five contrasts x two currencies, results in Sensitivity_<level><family>/ and
# figures in Sensitivity_Plot_<level><family>/, then one .rda with all runs.
SCENARIOS <- list(
  P1 = list("2004"      = 2004,             "2006"      = 2006),
  P2 = list("2018"      = 2018,             "2019"      = 2019),
  P3 = list("2004-2005" = c(2004, 2005),    "2006"      = 2006),
  P4 = list("2006"      = 2006,             "2018"      = 2018),
  PT = list("2004-2006" = c(2004, 2005, 2006), "2018-2019" = c(2018, 2019))
)

run_all_scenarios <- function(data_path = "data/dat_classed.rda") {
  load(data_path)
  diet_clean <- as.data.frame(dat_classed)
  if (!is.na(SPATIAL_SOURCE)) {
    if (!SPATIAL_SOURCE %in% names(diet_clean))
      stop("Column '", SPATIAL_SOURCE, "' not found in dat_classed.")
    diet_clean[[SPATIAL_SOURCE]] <- as.factor(diet_clean[[SPATIAL_SOURCE]])
  }

  cat("\n--- Stomachs per year ---\n")
  print(diet_clean %>% filter(year %in% c(2004:2006, 2018:2019)) %>%
          distinct(stomach_id, year) %>% count(year))
  check_prey_weight_column(diet_clean)

  out_dir_data  <- paste0("Sensitivity_", SPATIAL_LEVEL, PREY_FAMILY)
  out_dir_plots <- paste0("Sensitivity_Plot_", SPATIAL_LEVEL, PREY_FAMILY)

  runs <- list()
  for (sc in names(SCENARIOS)) {
    pm  <- SCENARIOS[[sc]]
    dat <- make_dat_classed(diet_clean, period_map = pm)
    for (mode in c("biomass", "occurrence")) {
      cat("\n=== ", SPATIAL_LEVEL, PREY_FAMILY, " | ", sc, " | ", mode, " ===\n", sep = "")
      runs[[paste0("res_", if (mode == "biomass") "biomass" else "occ", "_", sc)]] <-
        run_save_plot(dat, sc, mode, names(pm)[1], names(pm)[2],
                      out_dir = out_dir_data, plot_dir = out_dir_plots)
    }
  }

  dir.create("data/Sensitivity", recursive = TRUE, showWarnings = FALSE)
  list2env(runs, envir = environment())
  save(list = names(runs),
       file = paste0("data/Sensitivity/all_runs_", SPATIAL_LEVEL, PREY_FAMILY, ".rda"))
  cat("\nDone: ", SPATIAL_LEVEL, PREY_FAMILY, " (", length(runs), " runs saved)\n", sep = "")
  invisible(runs)
}

# =============================================================================
# LITERATURE CITED
# =============================================================================
# Anderson, M. J. (2001). A new method for non-parametric multivariate analysis
#   of variance. Austral Ecology, 26, 32-46.
# Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate
#   dispersions. Biometrics, 62, 245-253.
# Bolnick, D. I., et al. (2002). Measuring individual-level resource
#   specialization. Ecology, 83, 2936-2941.
# Hurlbert, S. H. (1978). The measurement of niche overlap and some relatives.
#   Ecology, 59, 67-77.
# Hurlbert, S. H. (1984). Pseudoreplication and the design of ecological field
#   experiments. Ecological Monographs, 54, 187-211.
# Kish, L. (1965). Survey Sampling. Wiley.
# Levins, R. (1968). Evolution in Changing Environments. Princeton Univ. Press.
# Pennington, M., & Volstad, J. H. (1994). Assessing the effect of intra-haul
#   correlation and variable density on estimates of population characteristics
#   from marine surveys. Biometrics, 50, 725-732.
# Roughgarden, J. (1972). Evolution of niche width. Am. Nat., 106, 683-718.
# Welch, B. L. (1947). The generalization of "Student's" problem when several
#   different population variances are involved. Biometrika, 34, 28-35.
# Zahl, S. (1977). Jackknifing an index of diversity. Ecology, 58, 907-913.
# =============================================================================
