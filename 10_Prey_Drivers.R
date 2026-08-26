# =============================================================================
# 10_Prey_Drivers.R - PREY ASSOCIATED WITH THE INTER-DECADE TURNOVER
# -----------------------------------------------------------------------------
# Identifies the prey taxa whose relative abundance differs between the two
# decades. These are associations, not mechanisms: nothing here holds predator
# identity, prey availability or environmental change constant, so a taxon
# appearing below is a candidate for interpretation against independent
# evidence on the Gulf prey field, not a demonstrated cause of the diet shift.
# The word "driver" is kept only in file and column names, for continuity with
# the engine's output schema.
#
# Replaces the exploratory SIMPER of the earlier draft. SIMPER ranks taxa by
# their contribution to the average Bray-Curtis dissimilarity between groups,
# which confounds a difference in mean abundance with a difference in variance:
# a highly variable prey scores high whether or not its mean changed (Warton et
# al. 2012). The primary analysis here is therefore model-based.
#
#   PRIMARY   manylm on Hellinger-transformed set x prey matrices, with
#             resampling-based univariate p-values adjusted step-down across
#             taxa. Answers: for which prey does the period term matter, once
#             the mean-variance relationship is modelled?
#
#   PRIMARY   Indicator Value (Dufrene & Legendre 1997). Answers: which prey
#             are both faithful to and abundant in one period. Complements
#             the model by ranking on ecological indicator strength rather
#             than on statistical significance alone.
#
#   LEGACY    SIMPER, retained and clearly labelled as descriptive, so the
#             convergence with the earlier draft can be checked and reported.
#             It is not the basis for any claim.
#
# These are the same two tests the engine runs per cell in run_mglm_indval();
# this script applies them at the assemblage level, which is what the
# Discussion needs.
#
# WHAT THIS SCRIPT PRODUCES
# -----------------------------------------------------------------------------
#   drivers_mglm_indval_<currency>.csv   The driver table. One row per prey:
#                                        MGLM p-value, IndVal and its p, the
#                                        period each prey indicates, mean share
#                                        per period, and the change in share.
#                                        Sorted by absolute change.
#
#   drivers_by_axis_<currency>.csv       The same, computed separately within
#                                        each ecoregion, size class and
#                                        predator, so the consistency of the
#                                        drivers across the assemblage is
#                                        visible rather than assumed.
#
#   drivers_resolution_sweep_<cur>.csv   The driver ranking at several
#                                        taxonomic resolutions. A prey that
#                                        only appears at one resolution is a
#                                        coding artefact, not a signal.
#
#   simper_legacy_<currency>.csv         The descriptive cross-check.
#
#   Fig8_nmds_sets_<currency>            Set-level NMDS, one 95 % ellipse per
#                                        period, driver prey as fitted vectors.
#
#   Fig8b_nmds_cells_<currency>          Cell centroids (predator x size x
#                                        ecoregion) with an arrow linking each
#                                        cell across decades. Controls for
#                                        predator identity; this is the version
#                                        for the manuscript.
#
#   permanova_aggregate_<currency>.txt   The whole-dataset test described in
#                                        Methods §2.4: D ~ period * Area +
#                                        size_class, with PERMDISP on period.
#                                        Indicative only — it pools predators,
#                                        sizes and areas.
#
# CAUTION TO CARRY INTO THE DISCUSSION
#   None of this holds predator identity, body size or area constant in the
#   pooled version. It is a convergence check against independently documented
#   changes in the Gulf prey field, not a mechanistic demonstration. The
#   per-axis table exists so that a driver claimed for the assemblage can be
#   shown to hold across its parts.
# =============================================================================

rm(list = ls())

# labdsv is attached last and masks vegan's scores(), calibrate(), density(),
# loadings() and predict() generics. labdsv's scores() has no method for
# metaMDS objects, so every ordination call below is qualified as
# vegan::scores(). Reordering the library() calls would fix it for this script
# but silently break any script that attaches them the other way round;
# qualifying is order-independent.
suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ggplot2)
  library(dplyr); library(tidyr); library(readr)
  library(mvabund); library(labdsv)
})

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
  # The PREY item weight, one row per prey item. Confirmed against the upstream
  # build; note that the neighbouring somatic_length_cm is a PREDATOR
  # measurement despite the shared prefix.
  col_wt     = "somatic_wt_g",

  # Resolution retained for the manuscript, and the sweep used to check that a
  # reported taxon is not an artefact of one aggregation threshold.
  PREY_FAMILY <- "_1",
  res_col     = paste0("prey_category_440", PREY_FAMILY),
  sweep_cols  = paste0("prey_category_", c(100, 440, 750), PREY_FAMILY),
  yr_early  = 2004:2006,
  yr_late   = 2018:2019,
  lab_early = "2004-2006",
  lab_late  = "2018-2019",

  currencies = c("biomass", "occurrence"),
  drop_prey  = c("Empty", "empty", "Unidentified", "unidentified",
                 "Digested", "digested", "NA"),
  min_sets_per_prey = 2,    # a prey must occur in at least this many sets
  seed  = 42,
  perms = 999,
  alpha = 0.05,
  top_n_vectors = 12,       # driver arrows on the NMDS
  run_legacy_simper = TRUE
)

set.seed(CFG$seed)
PER_COL <- c("2004-2006" = "#2C6E91", "2018-2019" = "#D2691E")

save_drv <- function(p, stem, w, h, dpi = 320) save_fig(p, stem, w, h, dpi, dir = DIR_DRIVERS)

# =============================================================================
# 1. LOAD AND PREPARE
# =============================================================================
cat("\nLoading ", CFG$rda_path, "\n", sep = "")
e <- new.env(); loaded <- load(CFG$rda_path, envir = e)
obj <- if (CFG$obj_name %in% loaded) CFG$obj_name else loaded[1]
DT  <- as.data.table(get(obj, envir = e))

DT[, period := fifelse(get(CFG$col_year) %in% CFG$yr_early, CFG$lab_early,
                fifelse(get(CFG$col_year) %in% CFG$yr_late,  CFG$lab_late,
                        NA_character_))]
DT <- DT[!is.na(period)]
DT[, period := factor(period, levels = c(CFG$lab_early, CFG$lab_late))]

# The trawl set is the independent compositional unit, as everywhere else.
DT[, set_uid := paste(get(CFG$col_year), get(CFG$col_vessel), get(CFG$col_set),
                      sep = "_")]

# Fail early and informatively if the configured resolution column is absent:
# otherwise the first get() deep inside a matrix builder raises an opaque
# "object not found" with no indication of what is actually available.
avail_res <- grep("^prey_category_", names(DT), value = TRUE)
if (!CFG$res_col %in% names(DT)) {
  stop("CFG$res_col = '", CFG$res_col, "' is not a column of the data.\n",
       "Resolution columns available: ",
       if (length(avail_res)) paste(avail_res, collapse = ", ") else "(none)",
       "\nSet CFG$res_col to the resolution retained for the manuscript.")
}
missing_sweep <- setdiff(CFG$sweep_cols, names(DT))
if (length(missing_sweep)) {
  message("Sweep columns absent and skipped: ",
          paste(missing_sweep, collapse = ", "))
  CFG$sweep_cols <- intersect(CFG$sweep_cols, names(DT))
}
for (nm in c(CFG$col_pred, CFG$col_size, CFG$col_area, CFG$col_stom,
             CFG$col_year, CFG$col_vessel, CFG$col_set, CFG$col_wt)) {
  if (!nm %in% names(DT)) stop("Column '", nm, "' declared in CFG is not in the data.")
}

cat("Records: ", nrow(DT), " | sets: ", uniqueN(DT$set_uid),
    " | ", sum(DT$period == CFG$lab_early), " early / ",
    sum(DT$period == CFG$lab_late), " late\n", sep = "")

# -----------------------------------------------------------------------------
# Matrix builders. Rows are trawl sets, columns prey categories.
#   biomass    summed prey weight per set
#   occurrence number of stomachs in the set that contained the prey
# No square-root transform here: the multivariate model uses a Hellinger
# transform, and Bray-Curtis for the ordination is computed on row profiles.
# Applying both would down-weight abundance twice.
# -----------------------------------------------------------------------------
build_set_matrix <- function(dat, currency, res_col) {
  d <- copy(dat)
  d[, prey := as.character(get(res_col))]
  d <- d[!is.na(prey) & !(prey %chin% CFG$drop_prey)]
  if (!nrow(d)) return(NULL)

  agg <- if (currency == "biomass") {
    d[, .(val = sum(get(CFG$col_wt), na.rm = TRUE)), by = .(set_uid, prey)]
  } else {
    d[, .(val = uniqueN(get(CFG$col_stom))), by = .(set_uid, prey)]
  }

  m   <- dcast(agg, set_uid ~ prey, value.var = "val", fill = 0)
  mat <- as.matrix(m[, -1]); rownames(mat) <- m$set_uid
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  # A prey seen in a single set cannot support a period contrast.
  mat[, colSums(mat > 0) >= CFG$min_sets_per_prey, drop = FALSE]
}

# Set-level metadata. A set can hold several predators; the modal level is
# tagged for grouping and the predator count is kept so the mix is visible.
set_meta <- function(dat) {
  d <- dat[, .(set_uid,
               period = period,
               area = as.character(get(CFG$col_area)),
               size = as.character(get(CFG$col_size)),
               pred = as.character(get(CFG$col_pred)),
               stom = get(CFG$col_stom))]

  modal_by <- function(d, col) {
    z <- d[!is.na(get(col)), .(n = .N), by = c("set_uid", col)]
    setorderv(z, c("set_uid", "n"), c(1L, -1L))
    z <- unique(z, by = "set_uid")
    setnames(z, col, "value")
    z[, .(set_uid, value)]
  }
  m_area <- modal_by(d, "area"); setnames(m_area, "value", "Area")
  m_size <- modal_by(d, "size"); setnames(m_size, "value", "size_cl")
  m_pred <- modal_by(d, "pred"); setnames(m_pred, "value", "predator")

  base <- d[, .(period = period[1L], n_pred = uniqueN(pred),
                n_stom = uniqueN(stom)), by = set_uid]

  Reduce(function(a, b) merge(a, b, by = "set_uid", all.x = TRUE),
         list(base, m_area, m_size, m_pred))[]
}

# =============================================================================
# 2. THE ASSOCIATION TEST
# =============================================================================
# Returns one row per prey with both tests plus the descriptive change in share.
# Neither test alone is sufficient: MGLM says the period term matters, IndVal
# says the prey characterises a period, and dProp says by how much its share of
# the diet moved. A taxon worth reporting satisfies at least two of the three.
# Satisfying all three still establishes association only.
drivers_mglm_indval <- function(comm, groups, n_perm = CFG$perms,
                                alpha = CFG$alpha, label = "overall") {
  groups <- droplevels(factor(groups))
  if (nlevels(groups) != 2 || any(table(groups) < 3) ||
      ncol(comm) < 2 || nrow(comm) < 8) return(NULL)

  out <- tryCatch({
    mv  <- mvabund(as.matrix(decostand(comm, method = "hellinger")))
    fit <- manylm(mv ~ groups)
    an  <- anova(fit, test = "F", p.uni = "adjusted", nBoot = n_perm)

    uni_p <- an$uni.p["groups", ]
    uni_p <- uni_p[is.finite(uni_p)]

    iv <- indval(as.data.frame(decostand(comm, "total")),
                 as.integer(groups), numitr = n_perm)

    # Mean share of each prey in each period, the descriptive quantity.
    rel <- decostand(comm, "total")
    m1  <- colMeans(rel[groups == levels(groups)[1], , drop = FALSE])
    m2  <- colMeans(rel[groups == levels(groups)[2], , drop = FALSE])

    tibble(
      prey        = colnames(comm),
      mglm_p      = unname(uni_p[colnames(comm)]),
      indval      = as.numeric(iv$indcls)[match(colnames(comm), names(iv$indcls))],
      indval_p    = as.numeric(iv$pval)[match(colnames(comm), names(iv$indcls))],
      indval_grp  = levels(groups)[iv$maxcls][match(colnames(comm), names(iv$indcls))],
      share_early = unname(m1[colnames(comm)]),
      share_late  = unname(m2[colnames(comm)])
    ) %>%
      mutate(
        d_share   = share_late - share_early,
        direction = ifelse(d_share > 0, "gained", "lost"),
        mglm_sig  = !is.na(mglm_p)   & mglm_p   <= alpha,
        iv_sig    = !is.na(indval_p) & indval_p <= alpha,
        n_criteria = as.integer(mglm_sig) + as.integer(iv_sig),
        level     = label
      ) %>%
      arrange(desc(abs(d_share)))
  }, error = function(err) {
    message("  driver test failed for '", label, "': ", conditionMessage(err))
    NULL
  })
  out
}

# Legacy SIMPER, descriptive only.
simper_legacy <- function(comm, groups, label = "overall") {
  groups <- droplevels(factor(groups))
  if (nlevels(groups) != 2 || any(table(groups) < 2)) return(NULL)
  sm <- try(simper(comm, group = groups, permutations = CFG$perms), silent = TRUE)
  if (inherits(sm, "try-error")) return(NULL)
  ss <- summary(sm)[[1]]
  if (is.null(ss) || !nrow(ss)) return(NULL)
  out <- data.frame(prey = rownames(ss), ss, row.names = NULL,
                    stringsAsFactors = FALSE)
  # vegan has renamed these columns across versions; map only what is present.
  ren <- c(average = "contrib", sd = "contrib_sd", ava = "mean_early",
           avb = "mean_late", cumsum = "cum")
  for (old_nm in names(ren)) {
    if (old_nm %in% names(out)) names(out)[names(out) == old_nm] <- ren[[old_nm]]
  }
  out$level <- label
  tibble::as_tibble(out) %>% arrange(desc(.data[[intersect(c("contrib", "average"),
                                                           names(out))[1]]]))
}


# Stress diagnostics. A stress near zero on more than a handful of points is a
# degenerate solution, not a perfect fit: the ordination has collapsed onto a
# line or a pair of points and carries no within-group information. Above 0.2
# the configuration is not a reliable summary of the distances (Clarke 1993).
check_stress <- function(nm, label, n_points) {
  s <- nm$stress
  if (is.finite(s) && s < 0.001 && n_points > 10) {
    warning("Degenerate NMDS for ", label, ": stress = ", signif(s, 3),
            " on ", n_points, " points. The configuration has collapsed; do ",
            "not read within-period structure from it.", call. = FALSE)
  } else if (is.finite(s) && s > 0.2) {
    warning("High NMDS stress for ", label, ": ", round(s, 3),
            ". Treat the ordination as indicative only.", call. = FALSE)
  }
  invisible(s)
}

# =============================================================================
# 3. RUN, PER CURRENCY
# =============================================================================
for (cur in CFG$currencies) {

  cat("\n", strrep("=", 70), "\n", cur, "\n", strrep("=", 70), "\n", sep = "")

  comm <- build_set_matrix(DT, cur, CFG$res_col)
  if (is.null(comm) || !ncol(comm)) { message("No matrix for ", cur); next }

  meta <- set_meta(DT)
  meta <- meta[match(rownames(comm), set_uid)]
  stopifnot(identical(meta$set_uid, rownames(comm)))
  meta[, period := factor(period, levels = c(CFG$lab_early, CFG$lab_late))]

  cat("Sets: ", nrow(comm), " | prey categories: ", ncol(comm), "\n", sep = "")

  # ---- overall drivers ------------------------------------------------------
  cat("\nOverall driver test...\n")
  drv <- drivers_mglm_indval(comm, meta$period, label = "overall")

  if (!is.null(drv)) {
    write_csv(drv, file.path(DIR_DRIVERS,
                             sprintf("drivers_mglm_indval_%s.csv", cur)))

    top <- drv %>% filter(n_criteria >= 1) %>% head(20)
    cat("\nTop drivers (at least one test significant):\n")
    print(as.data.frame(
      top %>% transmute(prey,
                        mglm_p = signif(mglm_p, 3),
                        indval = round(indval, 3),
                        indval_p = signif(indval_p, 3),
                        indicates = indval_grp,
                        share_early = round(100 * share_early, 2),
                        share_late  = round(100 * share_late, 2),
                        d_share_pct = round(100 * d_share, 2),
                        direction)),
      row.names = FALSE)

    cat("\nAgreement between the two tests:\n")
    print(table(MGLM = drv$mglm_sig, IndVal = drv$iv_sig))
  }

  # ---- drivers within each axis --------------------------------------------
  # Consistency across the assemblage: a driver that only appears in one
  # ecoregion or one predator is not an assemblage-level result.
  cat("\nPer-axis driver tests...\n")
  run_axis <- function(axis_col, axis_name) {
    lv <- sort(unique(meta[[axis_col]]))
    bind_rows(lapply(lv, function(l) {
      idx <- which(meta[[axis_col]] == l)
      if (length(unique(meta$period[idx])) < 2) return(NULL)
      cm <- comm[idx, , drop = FALSE]
      cm <- cm[, colSums(cm > 0) >= CFG$min_sets_per_prey, drop = FALSE]
      d  <- drivers_mglm_indval(cm, meta$period[idx],
                                label = paste(axis_name, l, sep = ": "))
      if (!is.null(d)) d$axis <- axis_name
      d
    }))
  }

  drv_axis <- bind_rows(
    run_axis("Area",     "Ecoregion"),
    run_axis("size_cl",  "Size class"),
    run_axis("predator", "Predator")
  )

  if (nrow(drv_axis)) {
    write_csv(drv_axis, file.path(DIR_DRIVERS,
                                  sprintf("drivers_by_axis_%s.csv", cur)))

    # How many sub-levels each prey is a driver in: the consistency score.
    consistency <- drv_axis %>%
      filter(n_criteria >= 1) %>%
      count(axis, prey, direction, name = "n_levels") %>%
      pivot_wider(names_from = axis, values_from = n_levels, values_fill = 0) %>%
      arrange(desc(rowSums(across(where(is.numeric)))))
    write_csv(consistency, file.path(DIR_DRIVERS,
                                     sprintf("drivers_consistency_%s.csv", cur)))
    cat("\nMost consistent drivers across axes:\n")
    print(as.data.frame(head(consistency, 15)), row.names = FALSE)
  }

  # ---- resolution sweep -----------------------------------------------------
  # A prey appearing at one resolution only is a coding artefact. The
  # manuscript already flags part of the Amphipoda signal on these grounds.
  cat("\nResolution sweep...\n")
  sweep <- bind_rows(lapply(intersect(CFG$sweep_cols, names(DT)), function(rc) {
    cm <- build_set_matrix(DT, cur, rc)
    if (is.null(cm) || !ncol(cm)) return(NULL)
    mm <- set_meta(DT); mm <- mm[match(rownames(cm), set_uid)]
    d  <- drivers_mglm_indval(cm, factor(mm$period,
                                         levels = c(CFG$lab_early, CFG$lab_late)),
                              label = rc)
    if (!is.null(d)) d$resolution <- rc
    d
  }))
  if (nrow(sweep)) {
    write_csv(sweep, file.path(DIR_DRIVERS,
                               sprintf("drivers_resolution_sweep_%s.csv", cur)))
    stable <- sweep %>%
      filter(n_criteria >= 1) %>%
      count(prey, direction, name = "n_resolutions") %>%
      arrange(desc(n_resolutions))
    write_csv(stable, file.path(DIR_DRIVERS,
                                sprintf("drivers_resolution_stability_%s.csv", cur)))
    cat("Drivers recovered at all ", dplyr::n_distinct(sweep$resolution),
        " resolutions: ", sum(stable$n_resolutions ==
                                dplyr::n_distinct(sweep$resolution)), "\n", sep = "")
  }

  # ---- legacy SIMPER --------------------------------------------------------
  if (isTRUE(CFG$run_legacy_simper)) {
    cat("\nLegacy SIMPER (descriptive cross-check)...\n")
    sp <- simper_legacy(comm, meta$period)
    if (!is.null(sp)) {
      write_csv(sp, file.path(DIR_DRIVERS, sprintf("simper_legacy_%s.csv", cur)))
      if (!is.null(drv)) {
        agree <- inner_join(
          drv %>% filter(n_criteria >= 1) %>% select(prey, direction),
          sp  %>% head(25) %>% select(prey),
          by = "prey")
        cat("Prey in both the model-based set and the SIMPER top 25: ",
            nrow(agree), "\n", sep = "")
      }
    }
  }

  # ---- aggregate PERMANOVA (Methods 2.4) ------------------------------------
  cat("\nAggregate PERMANOVA...\n")
  D <- vegdist(decostand(comm, "total"), method = "bray")
  perm_df <- data.frame(period = meta$period,
                        Area = factor(meta$Area),
                        size_class = factor(meta$size_cl))
  ado <- try(adonis2(D ~ period * Area + size_class, data = perm_df,
                     permutations = CFG$perms, by = "terms"), silent = TRUE)
  if (!inherits(ado, "try-error")) {
    print(ado)
    capture.output(ado, file = file.path(DIR_DRIVERS,
                                         sprintf("permanova_aggregate_%s.txt", cur)))
  }

  bd <- betadisper(D, meta$period)
  bd_a <- anova(bd)
  print(bd_a)
  capture.output(bd_a, file = file.path(DIR_DRIVERS,
                                        sprintf("permdisp_period_%s.txt", cur)))

  # ---- NMDS, set level ------------------------------------------------------
  cat("\nNMDS (sets)...\n")
  nm <- try(metaMDS(decostand(comm, "total"), distance = "bray", k = 2,
                    trymax = 100, autotransform = FALSE, trace = 0),
            silent = TRUE)

  if (inherits(nm, "try-error")) {
    message("  Set-level NMDS failed for ", cur, ": ",
            conditionMessage(attr(nm, "condition")),
            "  -> Figure 8 skipped.")
  } else {
    cat(sprintf("  stress = %.3f\n", nm$stress))
    check_stress(nm, paste0("set-level ordination, ", cur), nrow(comm))
    scr <- as.data.frame(vegan::scores(nm, display = "sites"))
    scr$period <- meta$period

    p_sets <- ggplot(scr, aes(NMDS1, NMDS2, colour = period, fill = period)) +
      geom_point(alpha = 0.4, size = 1.5) +
      stat_ellipse(aes(group = period), type = "norm", level = 0.95, linewidth = 0.9)

    # Driver prey as fitted vectors: which taxa pull the ellipses apart.
    if (!is.null(drv)) {
      dv <- drv %>% filter(n_criteria >= 1) %>% head(CFG$top_n_vectors) %>% pull(prey)
      dv <- intersect(dv, colnames(comm))
      if (length(dv) >= 2) {
        ef <- envfit(nm, decostand(comm, "total")[, dv, drop = FALSE],
                     permutations = CFG$perms)
        ev <- as.data.frame(vegan::scores(ef, "vectors")); ev$prey <- dv
        ev$NMDS1 <- ev$NMDS1 * 0.8 * max(abs(scr$NMDS1))
        ev$NMDS2 <- ev$NMDS2 * 0.8 * max(abs(scr$NMDS2))
        p_sets <- p_sets +
          geom_segment(data = ev, inherit.aes = FALSE,
                       aes(0, 0, xend = NMDS1, yend = NMDS2),
                       arrow = arrow(length = unit(0.18, "cm")), colour = "grey30") +
          geom_text(data = ev, inherit.aes = FALSE,
                    aes(NMDS1, NMDS2, label = prey), size = 2.5, colour = "grey20")
      }
    }

    p_sets <- p_sets +
      scale_colour_manual(values = PER_COL, name = NULL) +
      scale_fill_manual(values = PER_COL, name = NULL) +
      labs(title = paste0("Set-level ordination, ", CURRENCY_LAB[[cur]], " currency"),
           subtitle = sprintf(paste0("Bray-Curtis on row profiles, stress = %.2f; ",
                                     "ellipses = 95%% per period.\nArrows are ",
                                     "prey associated with the period contrast."),
                              nm$stress)) +
      theme_classic(base_size = 11) + theme(legend.position = "top")

    save_drv(p_sets, sprintf("Fig8_nmds_sets_%s", cur), 7.4, 6.2)
  }

  # ---- NMDS, cell centroids with shift arrows -------------------------------
  # One point per predator x size x ecoregion x period, on relative diet. This
  # removes both the within-set noise and the predator-mix confound of the
  # set-level ordination, and is the version to put in the manuscript.
  cat("\nNMDS (cell centroids)...\n")
  dd <- copy(DT)
  dd[, prey := as.character(get(CFG$res_col))]
  dd <- dd[!is.na(prey) & !(prey %chin% CFG$drop_prey)]
  dd[, cell := paste(get(CFG$col_pred), get(CFG$col_size), get(CFG$col_area),
                     sep = "|")]

  agg <- if (cur == "biomass") {
    dd[, .(val = sum(get(CFG$col_wt), na.rm = TRUE)), by = .(cell, period, prey)]
  } else {
    # as.numeric matters: uniqueN() returns an integer, and data.table's `:=`
    # coerces the right-hand side to the existing column type, so dividing an
    # integer column in place truncates every proportion below 1 to zero and
    # the whole occurrence matrix silently becomes empty.
    dd[, .(val = as.numeric(uniqueN(get(CFG$col_stom)))), by = .(cell, period, prey)]
  }
  agg[, val := val / sum(val), by = .(cell, period)]
  agg[, key := paste(cell, period, sep = "@@")]

  mw  <- dcast(agg, key ~ prey, value.var = "val", fill = 0)
  cm  <- as.matrix(mw[, -1]); rownames(cm) <- mw$key
  ks  <- tstrsplit(rownames(cm), "@@", fixed = TRUE)
  cmeta <- data.table(cell = ks[[1]],
                      period = factor(ks[[2]],
                                      levels = c(CFG$lab_early, CFG$lab_late)))

  both <- cmeta[, .N, by = cell][N == 2, cell]     # arrows need both periods
  keep <- cmeta$cell %chin% both
  cm <- cm[keep, , drop = FALSE]; cmeta <- cmeta[keep]

  # A cell x period whose entire diet fell into drop_prey divides by a zero
  # total above and comes back as NaN; Bray-Curtis then returns missing values
  # and metaMDS fails with an uninformative error. Clean both cases here and
  # say how much was removed, because dropping a cell silently would change
  # which cells the arrows are drawn for.
  cm[!is.finite(cm)] <- 0
  cm <- cm[, colSums(cm) > 0, drop = FALSE]
  empty_rows <- rowSums(cm) <= 0
  if (any(empty_rows)) {
    message("  ", sum(empty_rows), " cell-period centroid(s) empty after ",
            "filtering and dropped.")
    cm <- cm[!empty_rows, , drop = FALSE]; cmeta <- cmeta[!empty_rows]
    # An arrow needs both periods; a cell that lost one of them must go too.
    both2 <- cmeta[, .N, by = cell][N == 2, cell]
    k2 <- cmeta$cell %chin% both2
    cm <- cm[k2, , drop = FALSE]; cmeta <- cmeta[k2]
  }

  if (nrow(cm) < 8) {
    message("  Only ", nrow(cm), " cell-period centroids: too few for an ",
            "ordination. Figure 8b skipped for ", cur, ".")
  } else {
    nmc <- try(metaMDS(cm, distance = "bray", k = 2, trymax = 100,
                       autotransform = FALSE, trace = 0), silent = TRUE)
    if (inherits(nmc, "try-error")) {
      message("  Cell-centroid NMDS failed for ", cur, ": ",
              conditionMessage(attr(nmc, "condition")),
              "  -> Figure 8b and the paired PERMANOVA are skipped.")
    } else {
      check_stress(nmc, paste0("cell centroids, ", cur), nrow(cm))
      sc <- as.data.table(vegan::scores(nmc, display = "sites"))
      sc[, `:=`(cell = cmeta$cell, period = cmeta$period)]
      arr <- dcast(sc, cell ~ period, value.var = c("NMDS1", "NMDS2"))
      setnames(arr, c("cell", "x0", "x1", "y0", "y1"))

      p_cell <- ggplot(sc, aes(NMDS1, NMDS2, colour = period)) +
        geom_segment(data = arr, inherit.aes = FALSE,
                     aes(x0, y0, xend = x1, yend = y1),
                     arrow = arrow(length = unit(0.12, "cm")),
                     colour = "grey70", linewidth = 0.3, alpha = 0.6) +
        geom_point(size = 1.8, alpha = 0.85) +
        stat_ellipse(aes(group = period), type = "norm", level = 0.95, linewidth = 1) +
        scale_colour_manual(values = PER_COL, name = NULL) +
        labs(title = paste0("Cell centroids, ", CURRENCY_LAB[[cur]], " currency"),
             subtitle = sprintf(paste0("One point per predator x size class x ",
                                       "ecoregion and period; stress = %.2f.\n",
                                       "Arrows link the same cell across decades."),
                                nmc$stress)) +
        theme_classic(base_size = 11) + theme(legend.position = "top")

      save_drv(p_cell, sprintf("Fig8b_nmds_cells_%s", cur), 7.4, 6.2)

      # Paired test: the period effect with each cell as its own block.
      Dc <- vegdist(cm, method = "bray")
      ado_c <- try(adonis2(Dc ~ period,
                           data = data.frame(period = cmeta$period),
                           permutations = how(blocks = factor(cmeta$cell),
                                              nperm = CFG$perms)),
                   silent = TRUE)
      if (!inherits(ado_c, "try-error")) {
        print(ado_c)
        capture.output(ado_c,
                       file = file.path(DIR_DRIVERS,
                                        sprintf("permanova_cellpaired_%s.txt", cur)))
      }
    }
  }
}

# =============================================================================
# 4. A NOTE FOR THE DISCUSSION
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("Drivers written to ", normalizePath(DIR_DRIVERS), "\n\n", sep = "")
cat("Reporting guidance:\n")
cat("  * Phrase every claim as association, not causation: 'prey associated\n")
cat("    with the transition', not 'prey driving the transition'.\n")
cat("  * Lead with the model-based result (MGLM + IndVal). SIMPER is written\n")
cat("    to disk as a cross-check only and should not carry a claim.\n")
cat("  * A prey is worth naming when it is significant on both tests AND is\n")
cat("    recovered at more than one taxonomic resolution. Check\n")
cat("    drivers_resolution_stability_*.csv before writing a taxon into the text.\n")
cat("  * Report the per-axis table alongside the overall one. A driver present\n")
cat("    in one ecoregion only is not an assemblage-level result.\n")
cat("  * The aggregate PERMANOVA pools predators, sizes and areas. Keep it\n")
cat("    secondary to the cell-based framework, as Methods 2.4 already states.\n")
