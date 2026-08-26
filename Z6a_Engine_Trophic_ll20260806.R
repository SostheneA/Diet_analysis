# =============================================================================
# TROPHIC TRANSITION ANALYSIS ENGINE
# =============================================================================
#
# This framework classifies dietary change in marine populations between two
# temporal periods. The analytical unit throughout is the trawl set rather than
# the individual stomach, because stomachs from one tow are not independent
# (Hurlbert 1984; Pennington & Volstad 1994). What the engine produces is a
# classification of observed change, not a causal account of it.
#
# The pipeline evaluates three axes of dietary change—alpha diversity, niche
# breadth, and multivariate composition—to classify trophic transitions into
# discrete ecological states.
#
# -----------------------------------------------------------------------------
# THREE LEVELS OF SPATIAL DISAGGREGATION
# -----------------------------------------------------------------------------
# One engine, three spatial levels. The run script sets SPATIAL_LEVEL before
# sourcing this file; everything downstream is derived from it. Each level is a
# separate pipeline run and the three are independent, so they can be launched
# in three parallel R sessions.
#
#   L1  all_gulf    No spatial structure. Every trawl set enters a single
#                   PERMANOVA. Answers: did the Gulf-wide diet change?
#
#   L2  ecoregion   Cells are species x size class x ecoregion, tested within
#                   each ecoregion. Answers: where did it change?
#
#   L3  stratum     Same, at the survey-stratum grain — the finest available.
#
# All three run on the same stomachs: the source data carries no missing Area
# and no missing stratum, so no level drops rows the others keep. Differences
# between levels are therefore attributable to the spatial treatment alone.
# make_dat_classed() still reports any row it has to drop, so a future version
# of the data that does contain gaps will announce itself rather than shift the
# comparison silently.
#
# What each comparison isolates:
#   L1 vs L2   aggregation: pooling the diet matrices before testing, versus
#              testing within units
#   L2 vs L3   spatial grain, at constant sample
#
# TWO KINDS OF POOLING — DO NOT CONFLATE THEM IN THE METHODS
#   Pooling the DATA (level L1) merges every trawl set into one matrix and runs
#   a single test. Spatial heterogeneity becomes within-group variance, which
#   inflates dispersion; betadisper detects it, and it can either mask a real
#   shift or manufacture one (Warton et al. 2012).
#
#   Pooling the RESULTS (the `_summary_` figure emitted by L2 and L3) counts how
#   many local units fall in each diagnostic class. It is immune to the
#   dispersion problem, but its n counts cells rather than independent
#   replicates, which is what the n_eff correction in make_diag_plot2() handles.
#
#   The two answer different questions and both are produced. L1 is the
#   Gulf-wide test; the `_summary_` figures are the proportion of local units
#   showing each pattern.
# =============================================================================

# =============================================================================
# 0. ENVIRONMENT SETUP AND GLOBAL PARAMETERS
# =============================================================================
# The engine relies on `data.table` for high-performance data manipulation,
# `vegan` for multivariate ecology metrics, and `mvabund` for resampling-based
# linear models.
#
# Analytical parameters are declared with `if (!exists(...))`, so a run script
# that sets them before sourcing this file keeps its own values; the assignments
# here only supply defaults. SPATIAL_SOURCE and SPATIAL_OUT are deliberate
# exceptions: they identify which spatial variant of the engine this file is,
# and are fixed by plain assignment.

packages <- c(
  "tidyverse", "vegan", "DT", "plotly", "kableExtra",
  "RColorBrewer", "htmltools", "sf", "data.table", "stringr",
  "mvabund", "labdsv"
)

new_pkg <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(new_pkg)) install.packages(new_pkg, dependencies = TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(vegan); library(data.table)
  library(sf); library(mvabund); library(labdsv)
})

if ("igraph" %in% (.packages())) detach("package:igraph", unload = TRUE)

# -----------------------------------------------------------------------------
# SPATIAL LEVEL
# -----------------------------------------------------------------------------
# SPATIAL_LEVEL selects one of the four levels. A run script sets it before
# sourcing this file; the default below makes an unattended source() run the
# broadest level rather than fail.
#
# Each entry declares:
#   source  column in diet_clean holding the spatial unit; NA for L1, which
#           needs no spatial information at all
#   out     name of the spatial column written to the results, kept as "Area"
#           and "str" so existing mapping scripts keep working
#   pool    TRUE collapses every stomach into one unit before testing
#   label   value written in the output column when pool = TRUE; distinct per
#           level so L1 and L2 results stay distinguishable if row-bound
#   title   human-readable name used in plot captions and messages

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

# Analytical Thresholds
if (!exists("OCC_BINARY"))      OCC_BINARY      <- FALSE    # If TRUE, Jaccard for occurrence; else Bray-Curtis
if (!exists("MIN_SETS"))        MIN_SETS        <- 3        # Minimum sets per period to run jackknife
if (!exists("N_MIN"))           N_MIN           <- 5        # Minimum stomachs per cell
if (!exists("N_STRICT"))        N_STRICT        <- 25       # Minimum stomachs per period for reliable inference
if (!exists("ALPHA"))           ALPHA           <- 0.05     # Significance threshold
if (!exists("R_PERM"))          R_PERM          <- 999      # Permutations for distance-based tests

# Methodological Switches
if (!exists("BS_TEST"))         BS_TEST         <- "welch"    # Welch avoids the combinatorial floor of rank tests
if (!exists("ADD_NICHE_PART"))  ADD_NICHE_PART  <- TRUE       # Set-level niche width decomposition
if (!exists("ADD_IND_NICHE"))   ADD_IND_NICHE   <- FALSE      # Individual-level WIC/TNW (computationally heavy)
if (!exists("DRIVER_PUNI"))     DRIVER_PUNI     <- "adjusted" # Step-down adjusted p-values across taxa
if (!exists("DRIVER_RELATIVE")) DRIVER_RELATIVE <- TRUE       # Standardize abundance for Indicator Value
if (!exists("COMP_RELATIVE"))   COMP_RELATIVE   <- TRUE       # Standardize row sums before Bray-Curtis

# NEW line 156****** -----------------------------------------------------------
# Column cleaning (over-resolution). Complements N_MIN: N_MIN filters CELLS
# (rows / analysis units); clean_matrix() filters PREY CATEGORIES (columns) that
# are too sparse within a cell. A category is kept only if it clears BOTH floors:
#   S_{g,c}  >= K_COL   (stomachs in the cell containing prey g)
#   FO_{g,c} >= FO_MIN  (frequency of occurrence = S_{g,c} / stomachs in cell)
# Master switch OFF by default, so existing results are unchanged.
if (!exists("CLEAN_MATRIX"))    CLEAN_MATRIX    <- FALSE
if (!exists("K_COL"))           K_COL           <- 5L         # min stomachs per prey per cell
if (!exists("FO_MIN"))          FO_MIN          <- 0.05       # min frequency of occurrence per cell

# =========================================================
# 1. PERIOD & SPATIAL PREPARATION
# =========================================================
# make_dat_classed() applies the period map and the row filters, and stores the
# spatial unit under the reserved name `.spatial_src`. Keeping it under its own
# name means pooling never overwrites the source column, so a pooled dataset can
# still be audited against the disaggregated one.
#
# Rows with a missing spatial value are dropped at every level except all_gulf,
# which is what makes L2, L3 and L4 run on the same stomachs.

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
# 2. MATRIX CONSTRUCTION AND THE REPLICATE UNIT
# =============================================================================
# The core mechanism to prevent pseudoreplication is aggregating stomach
# contents up to the **trawl set** level. The set becomes the definitive
# replicate for assessing ecological shifts across periods.
#
# For occurrence matrices, the data is normalized by the number of stomachs
# sampled in that specific set. This generates a genuine proportional gradient
# rather than a coarse binary presence/absence grid.

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

# NEW line 251*****-------------------------------------------------------------
# --- Column cleaning helpers (over-resolution) --------------------------------
# Prey categories to KEEP within a cell. Support is always computed from stomach
# counts (S = distinct stomachs containing the prey; FO = S / stomachs in cell),
# then the surviving names are used to drop columns from whichever matrix (set-
# or stomach-level) was built. Applied only when CLEAN_MATRIX = TRUE.
.keep_prey_cols <- function(d, col_prey, k = K_COL, fo = FO_MIN) {
  if (!data.table::is.data.table(d)) d <- data.table::as.data.table(d)
  d2 <- d[!is.na(get(col_prey)) & !is.na(stomach_id)]
  n_sto <- data.table::uniqueN(d2$stomach_id)
  if (n_sto == 0) return(character(0))
  supp <- unique(d2[, .(stomach_id, prey = get(col_prey))])
  supp <- supp[, .(S = .N), by = prey][, FO := S / n_sto]
  supp[S >= k & FO >= fo, prey]
}

# Restrict a matrix to the kept columns, then drop rows emptied by the cut.
clean_matrix <- function(mat, keep) {
  if (is.null(mat)) return(NULL)
  cols <- intersect(colnames(mat), keep)
  if (length(cols) == 0) return(NULL)
  mat <- mat[, cols, drop = FALSE]
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)
  mat
}


# Build Set x Prey matrix (Biomass)
build_biomass_set <- function(data,
                              sp = NULL,
                              sz = NULL,
                              ar = NULL,
                              per = NULL,
                              col_prey) {
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
  # NEW line 267 ******---------------------------------------------------------
  if (isTRUE(CLEAN_MATRIX)) { mat <- clean_matrix(mat, .keep_prey_cols(d, col_prey)); if (is.null(mat)) return(NULL) }
  mat
}

# Build Set x Prey matrix (Frequency of Occurrence) - it's a proportion of occurrence
build_occ_set <- function(data,
                          sp = NULL,
                          sz = NULL,
                          ar = NULL,
                          per = NULL,
                          col_prey) {
  d <- .filter_cell(data, sp, sz, ar, per)
  d <- d[!is.na(d[[col_prey]]) & !is.na(set) & !is.na(stomach_id)]
  if (nrow(d) == 0) return(NULL)
  if (dplyr::n_distinct(d$stomach_id) < N_MIN) return(NULL)

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
  # NEW line 292 ---------------------------------------------------------------
  if (isTRUE(CLEAN_MATRIX)) { mat <- clean_matrix(mat, .keep_prey_cols(d, col_prey)); if (is.null(mat)) return(NULL) }
  mat
}

# Build Individual Stomach Matrix (for Optional Niche Partitioning)
build_biomass_stomach <- function(data,
                                  sp = NULL,
                                  sz = NULL,
                                  ar = NULL,
                                  per = NULL,
                                  col_prey) {
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
  # NEW line 313 ********-------------------------------------------------------
  if (isTRUE(CLEAN_MATRIX)) { mat <- clean_matrix(mat, .keep_prey_cols(d, col_prey)); if (is.null(mat)) return(NULL) }

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

# =============================================================================#
# =============================================================================#
# =============================================================================#
# =============================================================================#


# =============================================================================#
# SCENARIOS ----------
# =============================================================================#
period_map_P1 <- list(
  "2004" = c(2004),
  "2006" = c(2006)
)

period_map_P2 <- list(
  "2018" = c(2018),
  "2019" = c(2019)
)

period_map_P3 <- list(
  "2004-2005" = c(2004, 2005),
  "2006" = c(2006)
)

period_map_P4 <- list(
  "2006" = c(2006),
  "2018" = c(2018)
)

period_map_PT <- list(
  "2004-2006" = c(2004, 2005, 2006),
  "2018-2019" = c(2018, 2019)
)


# =============================================================================#
# TOY DATA SETS -----
# =============================================================================#
path <- here::here()
load(paste0(path, "/data/dat_classed.rda"))

# build dataset
dat_classed_PT_area <- make_dat_classed(diet_clean = dat_classed,
                                        period_map = period_map_PT)

verif <- dat_classed_PT_area %>%
  distinct(stomach_id, year) %>%
  count(stomach_id) %>%
  filter(n > 1)

print(verif)




# build a %W  matrix on trawl set as unit  for the period 1 ----------------
m1_bio_set <- build_biomass_set(
  data     = dat_classed_PT_area,
  sp       = "Atlantic cod",
  sz       = "adult",
  ar       = "Magdalen Shallows",
  per      = "2004-2006",
  col_prey = "prey_category_100_2"
)

dim(m1_bio_set)
#head(m1_bio_set)
print(m1_bio_set[1:6, 1:6])
#View(m1_bio_set)

# build a %W  matrix on trawl set as unit  for the period 2 ----------------
m2_bio_set <- build_biomass_set(
  data     = dat_classed_PT_area,
  sp       = "Atlantic cod",
  sz       = "adult",
  ar       = "Magdalen Shallows",
  per      = "2018-2019",
  col_prey = "prey_category_100_2"
)

dim(m2_bio_set)
#head(m2_bio_set)
m2_bio_set[1:6, 1:6]
#View(m2_bio_set)

# build a %W  matrix on stomach as unit  for the period 1 ----------------
m1_bio_sto <- build_biomass_stomach(
  data     = dat_classed_PT_area,
  sp       = "Atlantic cod",
  sz       = "adult",
  ar       = "Magdalen Shallows",
  per      = "2004-2006",
  col_prey = "prey_category_100_2"
)

dim(m1_bio_sto)
#head(m1_bio_sto)
m1_bio_sto[1:6, 1:6]
#View(m1_bio_sto)

# build a %W  matrix on stomach as unit  for the period 2 ----------------
m2_bio_sto <- build_biomass_stomach(
  data     = dat_classed_PT_area,
  sp       = "Atlantic cod",
  sz       = "adult",
  ar       = "Magdalen Shallows",
  per      = "2018-2019",
  col_prey = "prey_category_100_2"
)

dim(m2_bio_sto)
#head(m2_bio_sto)
m2_bio_sto[1:6, 1:6]
#View(m2_bio_sto)

# build a %FO  matrix on set as unit  for the period 1 ----------------
m1_occ_set <- build_occ_set(
  data     = dat_classed_PT_area,
  sp       = "Atlantic cod",
  sz       = "adult",
  ar       = "Magdalen Shallows",
  per      = "2004-2006",
  col_prey = "prey_category_100_2"
)

dim(m1_occ_set)
#head(m1_occ_set)
m1_occ_set[1:6, 1:6]
#View(m1_occ_set)

# build a %FO  matrix on set as unit  for the period 2 ----------------
m2_occ_set <- build_occ_set(
  data     = dat_classed_PT_area,
  sp       = "Atlantic cod",
  sz       = "adult",
  ar       = "Magdalen Shallows",
  per      = "2018-2019",
  col_prey = "prey_category_100_2"
)

dim(m2_occ_set)
#head(m2_occ_set)
m2_occ_set[1:6, 1:6]
#View(m2_occ_set)

# align mat ----------------------------------------------------------------

m_bio_set <- align_mats(m1 = m1_bio_set, m2 = m2_bio_set)
m_occ_set <- align_mats(m1 = m1_occ_set, m2 = m2_occ_set)

m_bio_sto <- align_mats(m1 = m1_bio_sto, m2 = m2_bio_sto)


# =============================================================================#
# =============================================================================#
# =============================================================================#
# =============================================================================#

# -----------------------------------------------------------------------------
# Guard: the weight column must hold the PREY weight
# -----------------------------------------------------------------------------
# build_biomass_set() sums `somatic_wt_g` per set and prey category, which is
# only meaningful if that column carries the weight of the prey item on each
# row. It does: confirmed against the upstream build.
#
# The check below is kept as a regression guard, not an open question. The
# column name invites the opposite reading — "somatic weight" conventionally
# means the gutted weight of a fish, and the neighbouring `somatic_length_cm`
# IS a predator measurement — so a future change to the build could substitute
# a predator weight without anything downstream complaining. It would not
# error; it would quietly turn the biomass currency into a count of prey
# records weighted by predator size.
#
# NOTE FOR THE DATA DICTIONARY AND THE METHODS
#   `somatic_length_cm` describes the PREDATOR (it sets the size class relative
#   to length at maturity); `somatic_wt_g` describes the PREY item. Two columns
#   sharing the `somatic_` prefix describe different entities. State this
#   explicitly wherever the data are archived.
#
# Two signatures separate the two cases:
#   * a predator attribute is constant across the several prey rows of one
#     stomach; a prey weight varies between them;
#   * a predator weight scales allometrically with predator length, with a
#     log-log slope near 3.
#
# Called by the run scripts before the sweep. It warns rather than stops, and
# under the current build it should print a low constant-within-stomach share
# and a flat slope, and raise nothing.
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
# 3. ALPHA DIVERSITY AND NICHE BREADTH
# =============================================================================
# To test differences in Shannon diversity (H'), the pipeline utilizes a
# delete-one-set jackknife resampling procedure combined with a Welch's t-test
# (Zahl, 1977). For Levins' niche breadth (Bs), the metric is standardized
# using a fixed n_cat representing the total prey categories available in the
# pooled temporal cell.

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


# =============================================================================#
# =============================================================================#
# =============================================================================#
# =============================================================================#
# Shannon (H) : jackknife entre périodes ---
H_bio_p <- shannon_jackknife_set(m_bio_set$m1, m_bio_set$m2)
H_occ_p <- shannon_jackknife_set(m_occ_set$m1, m_occ_set$m2)

# Levins standardisé (Bs) entre les périodes ----
levins_jackknife_set <- function(mat_set1, mat_set2, n_cat = NULL, min_sets = MIN_SETS) {
  if (is.null(n_cat)) n_cat <- ncol(mat_set1)
  Bs_from_mat <- function(m) {
    s <- sum(m)
    if (!is.finite(s) || s <= 0 || n_cat < 2) return(NA_real_)
    p <- colSums(m) / s; p <- p[p > 0]
    if (length(p) == 0) return(NA_real_)
    ((1 / sum(p^2)) - 1) / (n_cat - 1)
  }
  jackknife_B <- function(m) {
    n <- nrow(m)
    if (n < min_sets) return(list(B = NA_real_, var = NA_real_, n = n, ok = FALSE))
    B_full <- Bs_from_mat(m)
    if (!is.finite(B_full)) return(list(B = NA_real_, var = NA_real_, n = n, ok = FALSE))
    B_i <- vapply(seq_len(n), function(i) Bs_from_mat(m[-i, , drop = FALSE]), numeric(1))
    if (any(!is.finite(B_i))) return(list(B = NA_real_, var = NA_real_, n = n, ok = FALSE))
    pseudo <- n * B_full - (n - 1) * B_i
    list(B = mean(pseudo), var = stats::var(pseudo) / n, n = n, B_plugin = B_full, ok = TRUE)
  }
  j1 <- jackknife_B(mat_set1); j2 <- jackknife_B(mat_set2)
  base <- list(
    Bs1 = if (isTRUE(j1$ok)) round(j1$B_plugin, 3) else NA_real_,
    Bs2 = if (isTRUE(j2$ok)) round(j2$B_plugin, 3) else NA_real_,
    Bs1_jack = round(j1$B, 3), Bs2_jack = round(j2$B, 3),
    delta = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_,
    n1_sets = j1$n, n2_sets = j2$n
  )
  if (!isTRUE(j1$ok) || !isTRUE(j2$ok) || !is.finite(j1$var) ||
      !is.finite(j2$var) || (j1$var + j2$var) <= 0) return(base)
  t_val <- (j2$B - j1$B) / sqrt(j1$var + j2$var)
  df <- (j1$var + j2$var)^2 / (j1$var^2 / (j1$n - 1) + j2$var^2 / (j2$n - 1))
  base$delta <- round(j2$B - j1$B, 3); base$t <- round(t_val, 3); base$df <- round(df, 1)
  base$p <- round(2 * stats::pt(-abs(t_val), df = max(df, 1)), 4)
  base
}

Bs_bio_p <- levins_jackknife_set(m_bio_set$m1, m_bio_set$m2)
Bs_occ_p <- levins_jackknife_set(m_occ_set$m1, m_occ_set$m2)



# Levins standardisé (Bs) : indices par set, puis comparaison ---
ncat_bio <- ncol(m_bio_set$m1)
ncat_occ <- ncol(m_occ_set$m1)

Bs_bio1 <- row_indices_set(m_bio_set$m1, n_cat = ncat_bio)$Bs
Bs_bio2 <- row_indices_set(m_bio_set$m2, n_cat = ncat_bio)$Bs
Bs_occ1 <- row_indices_set(m_occ_set$m1, n_cat = ncat_occ)$Bs
Bs_occ2 <- row_indices_set(m_occ_set$m2, n_cat = ncat_occ)$Bs

Bs_bio_set <- compare_index_sets(Bs_bio1, Bs_bio2)
Bs_occ_set <- compare_index_sets(Bs_occ1, Bs_occ2)

# H par trait, comparé par Welch (même cadre que Bs)
H_bio_set <- compare_index_sets(
  row_indices_set(m_bio_set$m1)$H,
  row_indices_set(m_bio_set$m2)$H
)
H_occ_set <- compare_index_sets(
  row_indices_set(m_occ_set$m1)$H,
  row_indices_set(m_occ_set$m2)$H
)


# LL20260806 choisir un type ou justifier dans le texte en expliquant pourquoi on choisit un plutôt que l'autre selon le test.
list(H_W_p = H_bio_p,
     H_FO_p = H_occ_p,
     Bs_W_p = Bs_bio_p,
     Bs_FO_p = Bs_occ_p)

list(H_W_set = H_bio_set,
     H_FO_set = H_occ_set,
     Bs_W_set = Bs_bio_set,
     Bs_FO_set = Bs_occ_set)


# =============================================================================#
# Que mesurent H (Shannon) et Bs (Levins standardisé) ?
# -----------------------------------------------------------------------------

# ------------------------------------------------------------------
# Indices de diversité du régime alimentaire (par trawl set)
# ------------------------------------------------------------------
#
# Indice de Shannon (H') et indice de Levins (B) : deux mesures
# distinctes calculees sur les proportions d'abondance de proies (p_i).
#
# Shannon : H' = -sum(p_i * ln(p_i))
#   -> pondere la richesse et l'equitabilite des proportions.
#   -> sensible aux proies rares.
#
# Levins : B = 1 / sum(p_i^2)   (inverse de l'indice de Simpson)
#   -> domine par les proies abondantes, peu sensible aux rares.
#   -> souvent standardise : B_A = (B - 1) / (n - 1), borne entre 0 et 1.
#
#
# Shannon = diversite du regime (variete + equilibre des proies
#           consommees dans un trawl set).

# Levins  = largeur de niche trophique :
#             B_A proche de 1 = generaliste (proies reparties uniformement)
#             B_A proche de 0 = specialiste (concentre sur peu de proies)
#
# ------------------------------------------------------------------


# Bloc explicatif
# Il réutilise EXACTEMENT les deux fonctions du moteur (shannon_from_vec et
# levins_bs_from_vec), donc les figures décrivent bien ce que la section 3
# calcule, pas une réimplémentation.
#
# Rappel des trois quantités :
#   H  = -Σ p_j ln p_j            diversité de Shannon (en nats),
#   B  = 1 / Σ p_j^2              Levins = inverse de Simpson,
#   Bs = (B - 1) / (n_cat - 1)    Levins STANDARDISÉ, borné 0..1  (ce que renvoie
#                                 levins_bs_from_vec)

# "Largeur de niche de la population" = ces indices calculés sur la diète POOLÉE
# de la période (colSums de la matrice sets x proies). H et Bs décrivent la MÊME
# chose (la largeur de niche), mais avec un poids différent : Shannon (H) compte
# davantage les proies RARES, Levins/Simpson (Bs) les proies ABONDANTES.
# Sur les figures, H et Bs sont mis sur une échelle comparable 0..1 : Shannon est
# normalisé par ln(S) (H / ln S) et Bs est déjà borné à [0, 1] par construction.




EXPLAIN_NICHE <- TRUE
if (exists("EXPLAIN_NICHE") && isTRUE(EXPLAIN_NICHE)) {
  library(ggplot2)
  set.seed(4127)
  S <- 6

  lab_H  <- "H (Shannon, / ln S) "
  lab_Bs <- "Bs (Levins standardisé)"

  # FIGURE 1 : que mesurent H et Bs ? =====================================

  # Panneau A : une même diète, deux lectures ----------------------------
  # 6 catégories de proies. On fait varier la DOMINANCE de la proie principale
  # (p1), le reste étant réparti également sur les 5 autres.
  dom <- seq(1 / S, 0.95, length.out = 40)          # uniforme -> spécialiste
  panelA <- rbindlist(lapply(dom, function(p1) {
    p <- c(p1, rep((1 - p1) / (S - 1), S - 1))
    data.table(dominance = p1,
               H  = shannon_from_vec(p) / log(S),
               Bs = levins_bs_from_vec(p, n_cat = S))
  }))

  panelA <- melt(panelA, id.vars = "dominance",
                 variable.name = "indice", value.name = "valeur")
  panelA[, indice := factor(indice, levels = c("H", "Bs"),
                            labels = c(lab_H, lab_Bs))]

  pA <- ggplot(panelA, aes(dominance, valeur, colour = indice)) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Une même diète, deux lectures",
         subtitle = "Vers la droite : la proie principale domine (spécialisation)",
         x = "Dominance de la proie principale (p1)",
         y = "Largeur de niche (0 = spécialiste, 1 = uniforme)",
         colour = NULL) +
    theme(legend.position = "bottom")

  # Panneau B : sensibilité aux proies rares -----------------------------
  # Même nombre de catégories (6) partout ; seule l'ÉVENNESS change. H (normalisé)
  # reste plus haut que Bs tant qu'une proie domine, parce que Shannon compte plus
  # les proies rares.
  diets <- list(
    "Spécialiste"          = c(0.90, 0.06, 0.02, 0.01, 0.005, 0.005),
    "Généraliste inégal"   = c(0.50, 0.25, 0.13, 0.07, 0.03, 0.02),
    "Généraliste uniforme" = rep(1 / 6, 6)
  )
  panelB <- rbindlist(lapply(names(diets), function(nm) {
    p <- diets[[nm]] / sum(diets[[nm]])
    data.table(diete = nm,
               H  = shannon_from_vec(p) / log(S),
               Bs = levins_bs_from_vec(p, n_cat = S))
  }))
  panelB[, diete := factor(diete, levels = names(diets))]
  panelB <- melt(panelB, id.vars = "diete",
                 variable.name = "indice", value.name = "valeur")
  panelB[, indice := factor(indice, levels = c("H", "Bs"),
                            labels = c(lab_H, lab_Bs))]

  pB <- ggplot(panelB, aes(diete, valeur, fill = indice)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = round(valeur, 2)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_y_continuous(limits = c(0, 1.05)) +
    labs(title = "Sensibilité aux proies rares",
         subtitle = "H (normalisé) dépasse Bs quand une proie domine",
         x = NULL, y = "Largeur de niche (0..1)", fill = NULL) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 20, hjust = 1))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pA, pB, ncol = 2) +
            patchwork::plot_annotation(
              title = "Figure 1 — H (Shannon) et Bs (Levins standardisé)"))
  } else {
    print(pA); print(pB)
  }

  # FIGURE 2 : par set vs poolé, pour H puis pour Bs ======================
  # Une période = plusieurs trawl sets (lignes de la matrice). Chaque indice se
  # calcule de DEUX façons : (1) sur la diète POOLÉE de la période (colSums) =
  # valeur "plugin" du jackknife ; (2) par set, puis moyenne des sets
  # (compare_index_sets). Les deux diffèrent quand les sets sont hétérogènes :
  # le poolé est en général PLUS large que la moyenne des sets.
  n_sets <- 12
  mat <- t(sapply(seq_len(n_sets), function(i) {
    a <- rgamma(S, shape = runif(1, 0.2, 1.5)); a / sum(a) * rpois(1, 60)
  }))

  H_set   <- apply(mat, 1, shannon_from_vec)
  H_pool  <- shannon_from_vec(colSums(mat))
  Bs_set  <- apply(mat, 1, levins_bs_from_vec, n_cat = S)
  Bs_pool <- levins_bs_from_vec(colSums(mat), n_cat = S)

  pool_set_plot <- function(vals, pooled, ylab, ttl, ann) {
    ggplot(data.table(v = vals), aes(x = "P1", y = v)) +
      geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.7) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.3,
                   colour = "steelblue", linewidth = 0.4) +
      geom_hline(yintercept = pooled, linetype = "dashed", colour = "firebrick") +
      annotate("text", x = 1.4, y = pooled,
               label = paste0(ann, " = ", round(pooled, 2)),
               colour = "firebrick", vjust = -0.4, size = 3) +
      labs(title = ttl, x = NULL, y = ylab)
  }

  pH  <- pool_set_plot(H_set,  H_pool,  "H (Shannon, nats)",
                       "H : par set vs poolé", "H poolé")
  pBs <- pool_set_plot(Bs_set, Bs_pool, "Bs (Levins standardisé, 0..1)",
                       "Bs : par set vs poolé", "Bs poolé")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pH, pBs, ncol = 2) +
            patchwork::plot_annotation(
              title = "Figure 2 — par set vs poolé",
              subtitle = "points = valeur par set | bleu = moyenne des sets | rouge = valeur poolée"))
  } else {
    print(pH); print(pBs)
  }
}


# =============================================================================#
# FIGURE 3 (VRAIES DONNÉES) — par set vs poolé sur m_bio_set, P1 vs P2
# -----------------------------------------------------------------------------
# Même idée que la Figure 2, mais sur les matrices de la
# période (m_bio_set$m1 / $m2, en biomasse), et pour les deux périodes côte à
# côte. Ce bloc réutilise les objets DÉJÀ calculés par la section 3 :
#   H_bio_p  <- shannon_jackknife_set(m_bio_set$m1, m_bio_set$m2)
#   Bs_bio_p <- levins_jackknife_set(m_bio_set$m1, m_bio_set$m2)
#
# donc les repères tracés SONT les valeurs que le moteur rapporte :
#   - plugin poolé      : H1/H2  et Bs1/Bs2      (diète poolée, colSums)
#   - jackknife corrigé : H1_jack/H2_jack et Bs1_jack/Bs2_jack (delete-one-set)
#   - moyenne des sets  : moyenne des indices par set (voie compare_index_sets)
# Le t de Welch et la p-value affichés sous chaque panneau proviennent du
# jackknife (H_bio_p$t/$df/$p et Bs_bio_p$t/$df/$p).

EXPLAIN_NICHE_REAL <- TRUE
if (exists("EXPLAIN_NICHE_REAL") && isTRUE(EXPLAIN_NICHE_REAL) && exists("m_bio_set")) {
  library(ggplot2)
  m1 <- m_bio_set$m1
  m2 <- m_bio_set$m2
  ncat_bio <- ncol(m1)

  H_p1  <- apply(m1, 1, shannon_from_vec)
  H_p2  <- apply(m2, 1, shannon_from_vec)
  Bs_p1 <- apply(m1, 1, levins_bs_from_vec, n_cat = ncat_bio)
  Bs_p2 <- apply(m2, 1, levins_bs_from_vec, n_cat = ncat_bio)

  panel_real <- function(v1, v2, jr, plug1, plug2, jack1, jack2, ylab, ttl) {
    raw <- rbind(data.table(period = "P1", value = v1[is.finite(v1)]),
                 data.table(period = "P2", value = v2[is.finite(v2)]))

    summ <- rbind(
      data.table(period = "P1", type = "Sets (moyenne)",   value = mean(v1[is.finite(v1)])),
      data.table(period = "P2", type = "Sets (moyenne)",   value = mean(v2[is.finite(v2)])),

      data.table(period = "P1", type = "Periodes",     value = plug1),
      data.table(period = "P2", type = "Periodes",     value = plug2),

      data.table(period = "P1", type = "Periodes (jack corr)",  value = jack1),
      data.table(period = "P2", type = "Periodes (jack corr)",  value = jack2))

    summ[, type := factor(type, levels = c("Periodes", "Periodes (jack corr)",
                                           "Sets (moyenne)"))]
    ggplot() +
      geom_jitter(data = raw, aes(period, value),
                  width = 0.12, height = 0, colour = "grey65", alpha = 0.5, size = 1.4) +
      geom_point(data = summ, aes(period, value, colour = type, shape = type),
                 size = 2, stroke = 1.1, position = position_dodge(width = 0.6)) +
      scale_colour_manual(values = c("Periodes" = "firebrick",
                                     "Periodes (jack corr)" = "darkorange2",
                                     "Sets (moyenne)" = "steelblue")) +
      scale_shape_manual(values = c("Periodes" = 18, "Periodes (jack corr)" = 17,
                                    "Sets (moyenne)" = 15)) +
      labs(#title = ttl,
           subtitle = paste0("Welch (jackknife) : t = ", jr$t, "\n, df = ", jr$df,
                             ", p = ", jr$p, "  |  \u0394 = ", jr$delta),
           x = NULL, y = ylab, colour = NULL, shape = NULL) +
      theme(legend.position = "bottom")
  }

  pH_real  <- panel_real(H_p1,
                         H_p2,
                         H_bio_p,
                         H_bio_p$H1,
                         H_bio_p$H2,
                         H_bio_p$H1_jack,
                         H_bio_p$H2_jack,
                         "H (Shannon, nats)", "H")
  pBs_real <- panel_real(Bs_p1,
                         Bs_p2,
                         Bs_bio_p,
                         Bs_bio_p$Bs1,
                         Bs_bio_p$Bs2,
                         Bs_bio_p$Bs1_jack,
                         Bs_bio_p$Bs2_jack,
                         "Bs (Levins standardisé, 0..1)", "Bs")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(pH_real, pBs_real, ncol = 2, guides = "collect") +
            patchwork::plot_annotation(
              title = "Figure 3 (vraies données, biomasse)") &
            theme(legend.position = "bottom"))
  } else {
    print(pH_real); print(pBs_real)
  }
}

library(data.table); library(ggplot2)
m1 <- m_bio_set$m1; m2 <- m_bio_set$m2
ncat_bio <- ncol(m1)

# reproductions EXACTES des fonctions internes du jackknife (section 3)
H_from_mat <- function(m) {
  s <- sum(m); if (!is.finite(s) || s <= 0) return(NA_real_)
  p <- colSums(m) / s; p <- p[p > 0]; if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}
Bs_from_mat <- function(m, n_cat) {
  s <- sum(m); if (!is.finite(s) || s <= 0 || n_cat < 2) return(NA_real_)
  p <- colSums(m) / s; p <- p[p > 0]; if (length(p) == 0) return(NA_real_)
  ((1 / sum(p^2)) - 1) / (n_cat - 1)
}
pseudo <- function(m, fun) {
  n <- nrow(m); full <- fun(m)
  n * full - (n - 1) * vapply(seq_len(n), function(i) fun(m[-i, , drop = FALSE]), numeric(1))
}

ps <- rbind(
  data.table(indice = "H",  periode = "P1", val = pseudo(m1, H_from_mat)),
  data.table(indice = "H",  periode = "P2", val = pseudo(m2, H_from_mat)),
  data.table(indice = "Bs", periode = "P1", val = pseudo(m1, function(m) Bs_from_mat(m, ncat_bio))),
  data.table(indice = "Bs", periode = "P2", val = pseudo(m2, function(m) Bs_from_mat(m, ncat_bio)))
)

# Shapiro-Wilk par groupe (+ n, asymétrie, aplatissement)
sw <- ps[, {
  x <- val[is.finite(val)]
  m_ <- mean(x); s_ <- sd(x)
  .(n = length(x),
    W = round(unname(shapiro.test(x)$statistic), 3),
    p_shapiro = signif(shapiro.test(x)$p.value, 3),
    skew = round(mean(((x - m_) / s_)^3), 2),
    kurt = round(mean(((x - m_) / s_)^4) - 3, 2))
}, by = .(indice, periode)]
print(sw)


ggplot(ps, aes(sample = val)) +
  stat_qq(size = 1, alpha = 0.6) + stat_qq_line(colour = "firebrick") +
  facet_wrap(~ indice + periode, scales = "free", ncol = 2) +
  labs(title = "Normalité des pseudo-valeurs jackknife (biomasse)",
       subtitle = "QQ-plots ; ligne rouge = normale théorique",
       x = "Quantiles théoriques", y = "Pseudo-valeurs")

# =============================================================================#
# =============================================================================#
# =============================================================================#
# =============================================================================#



# =============================================================================
# 4. MULTIVARIATE COMPOSITION AND PREY ASSOCIATED WITH THE CHANGE
# =============================================================================
# Change in overall diet structure is tested by PERMANOVA (`adonis2`), with
# PERMDISP (`betadisper`) to check whether a significant result reflects a shift
# in location or a difference in dispersion.
#
# Prey associated with that change are identified by multivariate generalized
# linear models (`manylm`) on Hellinger-transformed data, together with
# Indicator Value analysis (`indval`). Both are associative: they identify prey
# whose relative abundance differs between periods, or which characterise one
# period. Neither establishes that a prey caused the compositional change. The
# `driver_prey` column name is retained for schema continuity with earlier
# outputs; read it as "prey associated with the transition".

# Failure counters. Both tests run inside tryCatch() and return NA on error, so
# without these a systematic failure (package API change, unexpected data shape)
# would leave every cell blank while the sweep finishes normally. run_pipeline()
# resets them at the start of a run and reports them at the end.
.comp_fail_count   <- 0L
.driver_fail_count <- 0L

run_composition_set <- function(s1, s2, mode = c("biomass", "occurrence"), perm = R_PERM, relative = COMP_RELATIVE) {
  mode <- match.arg(mode)
  empty <- list(BC = NA, R2 = NA, p_comp = NA, p_disp = NA, comp_dispersion = NA, n_set_P1 = 0L, n_set_P2 = 0L)

  al <- align_mats(s1, s2)
  if (is.null(al$m1) || is.null(al$m2)) return(empty)

  # Set counts are recorded before the early return, so a cell with too few tows
  # still reports its true depth. audit_set_depth() divides by
  # (n_set_P1 + n_set_P2) and would return Inf on a zero.
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
    # Failures are counted and the first few are printed, so a systematic cause
    # (vegan API change, adonis2 output layout change) is visible instead of
    # surfacing as a sweep where every cell is Inconclusive.
    .comp_fail_count <<- .comp_fail_count + 1L
    if (.comp_fail_count <= 5L) {
      message("run_composition_set failed (", .comp_fail_count, "): ",
              conditionMessage(e))
    }
    empty
  })
}

run_mglm_indval <- function(m1, m2, alpha = ALPHA, p_uni = DRIVER_PUNI, relative = DRIVER_RELATIVE, n_perm = R_PERM) {
  if (is.null(m1) || is.null(m2)) return(NA_character_)
  al <- align_mats(m1, m2)
  combined <- rbind(al$m1, al$m2)
  groups <- factor(c(rep("P1", nrow(al$m1)), rep("P2", nrow(al$m2))))

  keep_col <- colSums(combined > 0) >= 2
  combined <- combined[, keep_col, drop = FALSE]
  if (ncol(combined) < 2 || nrow(combined) < 4 || any(table(groups) < 2)) return(NA_character_)

  tryCatch({
    mv <- mvabund(as.matrix(decostand(combined, method = "hellinger")))
    fit <- manylm(mv ~ groups)
    an <- anova(fit, test = "F", p.uni = p_uni, nBoot = n_perm)

    # Non-finite entries are dropped from the p-value vector before its names
    # are read, so taxon labels stay aligned with the p-values they belong to.
    uni_p <- an$uni.p["groups", ]
    uni_p <- uni_p[is.finite(uni_p)]
    mglm_sig <- names(uni_p)[uni_p <= alpha]

    # numitr is passed explicitly: indval() otherwise runs its own default of
    # 1000 permutations, independent of R_PERM. as.integer() on the factor gives
    # the level indices, which is what maxcls indexes back into below.
    #
    # API check, worth running once after any labdsv update:
    #   args(labdsv:::indval.data.frame)
    # expects (x, clustering, numitr = 1000, ...). If that signature has
    # changed, this call will either error or silently ignore numitr, and the
    # per-cell permutation count would stop tracking R_PERM.
    iv_input <- if (relative) decostand(combined, "total") else combined
    iv <- indval(as.data.frame(iv_input), as.integer(groups), numitr = n_perm)

    indval_sig <- tibble(
      prey = names(iv$indcls), indval = as.numeric(iv$indcls),
      group = levels(groups)[iv$maxcls], p_indval = as.numeric(iv$pval)
    ) %>% filter(p_indval <= alpha)

    c1 <- colSums(al$m1[, colnames(combined), drop = FALSE])
    c2 <- colSums(al$m2[, colnames(combined), drop = FALSE])
    c1 <- c1 / sum(c1); c2 <- c2 / sum(c2)

    sig <- intersect(colnames(combined), union(mglm_sig, indval_sig$prey))
    if (length(sig) == 0) return(NA_character_)

    keep <- tibble(prey = sig) %>% mutate(dProp = unname(c2[prey] - c1[prey])) %>% arrange(desc(abs(dProp))) %>% pull(prey)

    out <- vapply(keep, function(px) {
      iv_row <- indval_sig[indval_sig$prey == px, , drop = FALSE]
      has_iv <- nrow(iv_row) > 0
      dprop <- unname(c2[px] - c1[px])
      sprintf("%s(MGLM=%s, IndVal=%s, grp=%s, p=%s, dProp=%+0.3f)",
              px, px %in% mglm_sig, ifelse(has_iv, round(iv_row$indval[1], 2), NA),
              ifelse(has_iv, iv_row$group[1], NA), ifelse(has_iv, round(iv_row$p_indval[1], 3), NA), dprop)
    }, character(1))

    paste(out, collapse = " | ")
  }, error = function(e) {
    # Same reporting logic as run_composition_set(): a driver_prey column that
    # is empty because manylm or indval failed on every cell is otherwise
    # indistinguishable from one where no prey reached significance.
    .driver_fail_count <<- .driver_fail_count + 1L
    if (.driver_fail_count <= 5L) {
      message("run_mglm_indval failed (", .driver_fail_count, "): ",
              conditionMessage(e))
    }
    NA_character_
  })
}

# =============================================================================#
# =============================================================================#
# =============================================================================#
# =============================================================================#

# -----------------------------------------------------------------------------
# remove_extreme_sets() — retire les TRAWL SETS à composition atypique
# -----------------------------------------------------------------------------
# Un trait peut avoir une composition très éloignée du reste de son groupe (p.
# ex. 2018_T_103) et étirer une ordination ou gonfler la dispersion. Cette
# fonction repère ces lignes en MULTIVARIÉ, de façon cohérente avec la section 4
# (même distance que run_composition_set : Bray-Curtis sur diète relative, ou
# Jaccard binaire si binary = TRUE).
#
# Critère : distance de Bray-Curtis au centroïde du GROUPE (via betadisper),
# puis seuil robuste PAR période = médiane + k * MAD (k = 3 par défaut). Le seuil
# est par groupe pour ne pas juger P1 et P2 à la même aune.
#
# Retourne toujours un `report` (une ligne par set, avec distance, seuil et
# drapeau `extreme`). Avec drop = TRUE (défaut) il renvoie aussi m1/m2 filtrées,
# colonnes d'origine conservées ; avec drop = FALSE il se contente de SIGNALER.
remove_extreme_sets <- function(m1, m2, relative = COMP_RELATIVE,
                                binary = OCC_BINARY, k = 3, drop = FALSE) {
  al <- align_mats(m1, m2)
  combined <- rbind(al$m1, al$m2)
  groups <- factor(c(rep("P1", nrow(al$m1)), rep("P2", nrow(al$m2))))

  d <- if (isTRUE(binary)) {
    vegan::vegdist(combined, method = "jaccard", binary = TRUE)
  } else {
    vegan::vegdist(if (relative) vegan::decostand(combined, "total") else combined,
                   method = "bray")
  }
  dist_c <- vegan::betadisper(d, groups)$distances

  # Seuil robuste par groupe ; MAD nul (groupe homogène) -> aucun retrait
  thr_by_grp <- tapply(dist_c, groups, function(x) {
    md <- mad(x)
    if (md == 0) Inf else median(x) + k * md
  })
  thr <- thr_by_grp[as.character(groups)]
  is_out <- dist_c > thr

  report <- tibble::tibble(
    set_uid       = rownames(combined),
    period        = groups,
    dist_centroid = round(dist_c, 3),
    threshold     = round(unname(thr), 3),
    extreme       = is_out
  )

  out <- list(report = report,
              removed = report[report$extreme, , drop = FALSE],
              n_removed = sum(is_out))
  if (drop) {
    keep1 <- rownames(al$m1)[!is_out[groups == "P1"]]
    keep2 <- rownames(al$m2)[!is_out[groups == "P2"]]
    out$m1 <- m1[rownames(m1) %in% keep1, , drop = FALSE]
    out$m2 <- m2[rownames(m2) %in% keep2, , drop = FALSE]
  }
  out
}

res_bio_set <- run_composition_set(m_bio_set$m1, m_bio_set$m2, mode = "biomass")
res_occ_set <- run_composition_set(m_occ_set$m1, m_occ_set$m2, mode = "occurrence")

# =============================================================================#
# APERÇU / INTERPRÉTATION — PERMANOVA + PERMDISP (composition par trawl set)
# -----------------------------------------------------------------------------
# La section 4 lance DEUX tests de permutation complémentaires sur la même
# matrice de distances (Bray-Curtis sur diète relative ; Jaccard binaire si
# OCC_BINARY = TRUE en occurrence) :
#
#   1. PERMANOVA  (adonis2)               -> res$p_comp
#      Teste si la POSITION (le centroïde) du nuage diffère entre P1 et P2,
#      c.-à-d. si la composition MOYENNE du régime a changé.
#
#   2. PERMDISP   (permutest(betadisper)) -> res$p_disp
#      Teste si la DISPERSION (variabilité inter-traits autour du centroïde)
#      diffère entre P1 et P2, c.-à-d. si un régime est plus hétérogène d'un
#      trait à l'autre que l'autre.
#
# POURQUOI LES DEUX. adonis2 est sensible à l'hétérogénéité de dispersion : un
# p_comp significatif SEUL ne prouve pas un déplacement de composition, il peut
# refléter une simple différence de dispersion. Le drapeau `comp_dispersion`
# (= TRUE quand p_comp ET p_disp sont significatifs) signale ce cas. Le signal
# combine alors un déplacement de centroïde ET une différence d'hétérogénéité :
# il faut rapporter LES DEUX plutôt que de conclure "la composition moyenne a
# changé". La dispersion accrue est elle-même interprétable (régime plus
# opportuniste/variable, hétérogénéité spatio-temporelle des proies, régime en
# transition), et n'est pas qu'un artéfact d'échantillonnage : l'effort par
# trait (nombre d'estomacs) corrèle avec la richesse détectée, mais il ne
# diffère pas systématiquement entre périodes.
#
# NOTE SUR LE DÉSÉQUILIBRE (Anderson & Walsh 2013). Quand le plus PETIT groupe
# est le plus dispersé, la PERMANOVA devient CONSERVATRICE ; un p_comp
# significatif dans ce cas est donc peu susceptible d'être un faux positif dû à
# la dispersion.
# -----------------------------------------------------------------------------

# Il réutilise les matrices m_occ_set (occurrence) et dat_classed (effort), et
# reproduit la distance EXACTEMENT comme run_composition_set() (mode occurrence).

EXPLAIN_COMPOSITION <- TRUE
if (exists("EXPLAIN_COMPOSITION") && isTRUE(EXPLAIN_COMPOSITION) &&
    exists("m_occ_set") && exists("dat_classed")) {
  library(vegan); library(ggplot2); library(dplyr); library(tibble)

  combined <- rbind(m_occ_set$m1, m_occ_set$m2)
  groups <- factor(c(rep("P1", nrow(m_occ_set$m1)), rep("P2", nrow(m_occ_set$m2))))

  use_jaccard <- isTRUE(OCC_BINARY)
  d <- if (use_jaccard) {
    vegdist(combined, method = "jaccard", binary = TRUE)
  } else {
    vegdist(decostand(combined, "total"), method = "bray")
  }

  set.seed(6820)
  ado <- adonis2(d ~ groups, permutations = R_PERM)
  bd  <- betadisper(d, groups)
  pmt <- permutest(bd, permutations = R_PERM)
  print(ado); print(pmt)

  # Effort par trait : nombre d'estomacs distincts par set_uid (même clé que
  # .add_set_uid : year_vessel.code_set)
  eff <- dat_classed |>
    mutate(set_uid = paste(year, vessel.code, set, sep = "_")) |>
    filter(!is.na(stomach_id), !is.na(set)) |>
    distinct(set_uid, stomach_id) |>
    count(set_uid, name = "n_stomachs")

  info <- tibble(
    set_uid       = rownames(combined),
    period        = groups,
    dist_centroid = bd$distances,
    richesse      = rowSums(combined > 0)
  ) |>
    left_join(eff, by = "set_uid")

  # Figure 1 : dispersion (PERMDISP) + effet de l'effort ------------------
  p_disp_fig <- ggplot(info, aes(period, dist_centroid)) +
    geom_boxplot(outlier.alpha = 0.3) +
    labs(title = "Dispersion (PERMDISP)",
         subtitle = sprintf("p = %.3f", pmt$tab$`Pr(>F)`[1]),
         x = NULL, y = "Distance au centroïde")

  p_eff_fig <- ggplot(info, aes(n_stomachs, richesse)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(#title = "Effort vs richesse détectée",
         subtitle = sprintf("Spearman rho = %.2f",
                            cor(info$n_stomachs, info$richesse,
                                method = "spearman", use = "complete.obs")),
         x = "Estomacs par trait", y = "Nb de proies détectées")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(p_disp_fig, p_eff_fig, ncol = 2))
  } else {
    print(p_disp_fig); print(p_eff_fig)
  }

  # Figure 2 : NMDS avec ellipses par période -----------------------------
  # Ellipses de type "t" (robustes) ; un trait très atypique peut étirer les
  # axes. Retirer les traits extrêmes en amont donne une ordination plus lisible.
  set.seed(6820)
  nmds <- metaMDS(d, k = 2, trymax = 100, trace = FALSE)
  print(
    vegan::scores(nmds, display = "sites") |>
      as_tibble(rownames = "set_uid") |>
      mutate(period = groups) |>
      ggplot(aes(NMDS1, NMDS2, colour = period)) +
      geom_point(alpha = 0.6) +
      stat_ellipse(aes(fill = period), geom = "polygon", alpha = 0.12, type = "t") +
      labs(title = "NMDS composition par trawl set",
           subtitle = sprintf("Stress = %.3f | PERMANOVA p = %.3f | PERMDISP p = %.3f",
                              nmds$stress, ado$`Pr(>F)`[1], pmt$tab$`Pr(>F)`[1]),
           colour = "Période", fill = "Période") +
      theme(legend.position = "top")
  )

  # Repèrer les outliers -------------------------------------------------------
  library(plotly)

  scores_df <- vegan::scores(nmds, display = "sites") |>
    as_tibble(rownames = "set_uid") |>
    mutate(period = groups)

  g <- ggplot(scores_df, aes(NMDS1, NMDS2, colour = period)) +
    geom_point(aes(text = paste0("set_uid: ", set_uid,
                                 "<br>NMDS1: ", round(NMDS1, 3),
                                 "<br>NMDS2: ", round(NMDS2, 3))),
               alpha = 0.6) +
    stat_ellipse(aes(fill = period), geom = "polygon", alpha = 0.12, type = "t") +
    labs(title = "NMDS composition par trawl set",
         colour = "Période", fill = "Période")

  ggplotly(g, tooltip = "text")



  # Retirer les outliers --------------------------------------------------------
  m2_bio_set <- m2_bio_set[rownames(m2_bio_set) != "2018_T_103", , drop = FALSE]
  m2_occ_set <- m2_occ_set[rownames(m2_occ_set) != "2018_T_103", , drop = FALSE]

  m_bio_set <- align_mats(m1 = m1_bio_set, m2 = m2_bio_set)
  m_occ_set <- align_mats(m1 = m1_occ_set, m2 = m2_occ_set)
  res_bio_set <- run_composition_set(m_bio_set$m1, m_bio_set$m2, mode = "biomass")
  res_occ_set <- run_composition_set(m_occ_set$m1, m_occ_set$m2, mode = "occurrence")


  # occurence ---
  combined <- rbind(m_occ_set$m1, m_occ_set$m2)
  groups <- factor(c(rep("P1", nrow(m_occ_set$m1)), rep("P2", nrow(m_occ_set$m2))))

  use_jaccard <- isTRUE(OCC_BINARY)
  d <- if (use_jaccard) {
    vegdist(combined, method = "jaccard", binary = TRUE)
  } else {
    vegdist(decostand(combined, "total"), method = "bray")
  }

  set.seed(6820)
  ado <- adonis2(d ~ groups, permutations = R_PERM)
  bd  <- betadisper(d, groups)
  pmt <- permutest(bd, permutations = R_PERM)

  nmds <- metaMDS(d, k = 2, trymax = 100, trace = FALSE)
  vegan::scores(nmds, display = "sites") |>
    as_tibble(rownames = "set_uid") |>
    mutate(period = groups) |>
    ggplot(aes(NMDS1, NMDS2, colour = period)) +
    geom_point(alpha = 0.6) +
    stat_ellipse(aes(fill = period), geom = "polygon", alpha = 0.12, type = "t") +
    labs(title = "NMDS composition par trawl set",
         subtitle = sprintf("Stress = %.3f | PERMANOVA p = %.3f | PERMDISP p = %.3f",
                            nmds$stress, ado$`Pr(>F)`[1], pmt$tab$`Pr(>F)`[1]),
         colour = "Période", fill = "Période") +
    theme(legend.position = "top")

  # biomass ------------
  combined_bio <- rbind(m_bio_set$m1, m_bio_set$m2)
  groups_bio <- factor(c(rep("P1", nrow(m_bio_set$m1)),
                         rep("P2", nrow(m_bio_set$m2))))

  d_bio <- vegdist(decostand(combined_bio, "total"), method = "bray")

  set.seed(6820)
  ado_bio <- adonis2(d_bio ~ groups_bio, permutations = R_PERM)
  bd_bio  <- betadisper(d_bio, groups_bio)
  pmt_bio <- permutest(bd_bio, permutations = R_PERM)

  nmds_bio <- metaMDS(d_bio, k = 2, trymax = 100, trace = FALSE)
  vegan::scores(nmds_bio, display = "sites") |>
    as_tibble(rownames = "set_uid") |>
    mutate(period = groups_bio) |>
    ggplot(aes(NMDS1, NMDS2, colour = period)) +
    geom_point(alpha = 0.6) +
    stat_ellipse(aes(fill = period), geom = "polygon", alpha = 0.12, type = "t") +
    labs(title = "NMDS composition par trawl set (biomasse)",
         subtitle = sprintf("Stress = %.3f | PERMANOVA p = %.3f | PERMDISP p = %.3f",
                            nmds_bio$stress, ado_bio$`Pr(>F)`[1], pmt_bio$tab$`Pr(>F)`[1]),
         colour = "Période", fill = "Période") +
    theme(legend.position = "top")
}

# =============================================================================#
# =============================================================================#
# =============================================================================#
# =============================================================================#


# =============================================================================
# LITERATURE CITED
# =============================================================================
# Anderson, M. J. (2001). A new method for non-parametric multivariate analysis
#   of variance. Austral Ecology, 26(1), 32-46.                    [PERMANOVA]
# Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate
#   dispersions. Biometrics, 62(1), 245-253.                      [betadisper]
# Bolnick, D. I., Yang, L. H., Fordyce, J. A., Davis, J. M., & Svanback, R.
#   (2002). Measuring individual-level resource specialization.
#   Ecology, 83(10), 2936-2941.                       [WIC/TNW, discrete form]
# Bolnick, D. I., Svanback, R., Fordyce, J. A., Yang, L. H., Davis, J. M.,
#   Hulsey, C. D., & Forister, M. L. (2003). The ecology of individuals:
#   incidence and implications of individual specialization.
#   The American Naturalist, 161(1), 1-28.
# Clarke, K. R. (1993). Non-parametric multivariate analyses of changes in
#   community structure. Australian Journal of Ecology, 18(1), 117-143.
#                                                        [SIMPER, superseded]
# Dufrene, M., & Legendre, P. (1997). Species assemblages and indicator
#   species: the need for a flexible asymmetrical approach.
#   Ecological Monographs, 67(3), 345-366.                            [IndVal]
# Hurlbert, S. H. (1978). The measurement of niche overlap and some relatives.
#   Ecology, 59(1), 67-77.                        [standardisation of Levins B]
# Hurlbert, S. H. (1984). Pseudoreplication and the design of ecological field
#   experiments. Ecological Monographs, 54(2), 187-211.
# Kish, L. (1965). Survey Sampling. Wiley.                    [design effect]
# Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful
#   transformations for ordination of species data.
#   Oecologia, 129(2), 271-280.                          [Hellinger transform]
# Levins, R. (1968). Evolution in Changing Environments. Princeton University
#   Press.                                                  [niche breadth B]
# Roughgarden, J. (1972). Evolution of niche width. The American Naturalist,
#   106(952), 683-718.                                             [WIC / TNW]
# Wang, Y., Naumann, U., Wright, S. T., & Warton, D. I. (2012). mvabund - an R
#   package for model-based analysis of multivariate abundance data.
#   Methods in Ecology and Evolution, 3(3), 471-474.
# Warton, D. I., Wright, S. T., & Wang, Y. (2012). Distance-based multivariate
#   analyses confound location and dispersion effects.
#   Methods in Ecology and Evolution, 3(1), 89-101.
#                                     [rationale for replacing SIMPER by MGLM]
# Welch, B. L. (1947). The generalization of "Student's" problem when several
#   different population variances are involved. Biometrika, 34(1/2), 28-35.
# Zahl, S. (1977). Jackknifing an index of diversity. Ecology, 58(4), 907-913.
#                                              [delete-one-set jackknife of H']
# =============================================================================
