# =============================================================================
# 6g_Compare_families.R - APPENDIX B: SENSITIVITY TO THE PREY-GROUPING RULE
# -----------------------------------------------------------------------------
# Compares the three prey-grouping families at each spatial level:
#   pooled_q1     "_1"   pooled across predators, q = 1 (main analysis)
#   pooled_q2     "_2"   pooled across predators, q = 2
#   per_predator  "_PP"  grouping rebuilt inside each predator (3PP)
# read from data/Sensitivity/all_runs_<level>_1 / _2 / _PP.rda (or the
# Sensitivity_<level><family>/ folders). A family or level without results
# is skipped.
#
# Outputs (Output_Appendices/):
#   FigB1_grouping_families_<level>          4 families x 2 currencies, the
#                                            three families in colour
#   FigB2_grouping_pct_significant_<level>   % of testable cells significant
#                                            per test
#   TableB1_grouping_families_<level>.csv (+ _wide), TableB2_..._<level>.csv
# Percentages are computed as in Tables 2-3: families from STATE_TO_FAMILY,
# within each resolution on the full family grid, then averaged across the
# resolutions common to the families present.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

PREY_FAMILY <- "_1"
source("R_helpers/Config_Mappings.R")

MAIN_PROJECT <- normalizePath(getwd())

ALPHA  <- 0.05
LEVELS <- SPATIAL_LEVELS_ORD
P_COLS <- c("p_comp", "p_H", "p_Bs", "p_disp")
TEST_LAB <- c(p_comp = "Composition (PERMANOVA)", p_H = "Diversity H'",
              p_Bs = "Niche breadth Bs", p_disp = "Dispersion (PERMDISP)")

APPROACHES <- list(
  pooled_q1    = list(root = MAIN_PROJECT, suffix = "_1"),
  pooled_q2    = list(root = MAIN_PROJECT, suffix = "_2"),
  per_predator = list(root = MAIN_PROJECT, suffix = "_PP")
)
APPROACH_LAB <- c(pooled_q1    = "Pooled grouping, q = 1 (main analysis)",
                  pooled_q2    = "Pooled grouping, q = 2",
                  per_predator = "Per-predator grouping")
APPROACH_PAL <- c(pooled_q1 = "#2F4A5A", pooled_q2 = "#7F9CAB", per_predator = "#D96C4A")

# Scenario -> contrast code, from the period labels stored in the results.
SCEN_MAP <- data.table(
  period_1 = c("2004", "2004-2005", "2018", "2006", "2004-2006"),
  period_2 = c("2006", "2006",      "2019", "2018", "2018-2019"),
  contrast = c("P1a",  "P1b",       "P2",   "PTb",  "PTa")
)
FAM_MAP <- data.table(diagnostic = names(STATE_TO_FAMILY),
                      family     = unname(STATE_TO_FAMILY))

save_app <- function(p, stem, w, h) save_fig(p, stem, w, h, dir = DIR_APPEND)

# --- Loading -----------------------------------------------------------------

.is_results <- function(x)
  is.data.frame(x) && all(c("x_threshold", "p_comp", "mode") %in% names(x))

.collect_results <- function(x) {
  if (.is_results(x)) return(list(as.data.table(x)))
  if (is.list(x) && !is.data.frame(x))
    return(unlist(lapply(x, .collect_results), recursive = FALSE))
  list()
}

.results_from_rda <- function(f) {
  e <- new.env(); load(f, envir = e)
  unlist(lapply(ls(e), function(o) .collect_results(get(o, envir = e))),
         recursive = FALSE)
}

load_approach <- function(root, suffix, spatial) {
  if (!dir.exists(root)) return(NULL)
  hit <- list.files(root, pattern = paste0("^all_runs_", spatial, suffix, "\\.rda$"),
                    recursive = TRUE, full.names = TRUE)
  if (length(hit)) {
    tabs <- .results_from_rda(hit[1])
    if (length(tabs)) return(list(dt = rbindlist(tabs, fill = TRUE), src = hit[1]))
  }
  d <- file.path(root, paste0("Sensitivity_", spatial, suffix))
  if (dir.exists(d)) {
    fs <- list.files(d, pattern = "\\.rda$", full.names = TRUE)
    if (length(fs)) {
      tabs <- unlist(lapply(fs, .results_from_rda), recursive = FALSE)
      if (length(tabs)) return(list(dt = rbindlist(tabs, fill = TRUE),
                                    src = paste0(d, " (", length(fs), " files)")))
    }
  }
  NULL
}

# Shaded bands: within-period contrasts on the left, between-period on the right.
period_bands <- function(n_within, n_total) {
  list(
    annotate("rect", xmin = 0.5, xmax = n_within + 0.5, ymin = -Inf, ymax = Inf,
             fill = "#D7E6EC", alpha = 0.55),
    annotate("rect", xmin = n_within + 0.5, xmax = n_total + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#F6E0D6", alpha = 0.65),
    geom_vline(xintercept = n_within + 0.5, linetype = "dashed",
               colour = "grey60", linewidth = 0.3)
  )
}

# =============================================================================
# LOOP OVER SPATIAL LEVELS
# =============================================================================
for (SPATIAL in LEVELS) {

  cat("\n=============", SPATIAL, "=============\n")

  res_list <- list()
  for (ap in names(APPROACHES)) {
    r <- load_approach(APPROACHES[[ap]]$root, APPROACHES[[ap]]$suffix, SPATIAL)
    if (!is.null(r)) {
      dt <- merge(r$dt[, approach := ap], SCEN_MAP,
                  by = c("period_1", "period_2"), all.x = TRUE)
      res_list[[ap]] <- dt
      xs <- sort(unique(dt$x_threshold))
      message(sprintf("[ok]     %-13s %s | %d rows | X = %d-%d",
                      ap, r$src, nrow(dt), min(xs), max(xs)))
    } else message("[absent] ", ap)
  }
  if (length(res_list) < 2) {
    message("Fewer than two approaches available for ", SPATIAL, " -> skipped")
    next
  }

  both <- rbindlist(res_list, fill = TRUE)
  both[, approach := factor(approach, levels = names(APPROACHES))]
  both[, contrast := factor(contrast, levels = CONTRAST_LEVELS)]

  # Resolutions common to every approach present (the PP sweep is coarser).
  common_x <- Reduce(intersect, lapply(res_list, function(d) unique(d$x_threshold)))
  sig_x    <- if (length(common_x)) common_x else unique(both$x_threshold)
  message(length(sig_x), " common resolutions used (",
          min(sig_x), "-", max(sig_x), ")")

  # ===========================================================================
  # B2. % SIGNIFICANT per approach x contrast x currency x test
  # ===========================================================================
  long <- melt(both[testable == TRUE & x_threshold %in% sig_x & !is.na(contrast),
                    c("approach", "contrast", "mode", "x_threshold",
                      intersect(P_COLS, names(both))), with = FALSE],
               id.vars = c("approach", "contrast", "mode", "x_threshold"),
               variable.name = "test", value.name = "p")

  pct_by_x <- long[!is.na(p), .(pct_sig = 100 * mean(p < ALPHA), n_cells = .N),
                   by = .(approach, contrast, mode, test, x_threshold)]

  sig_summary <- pct_by_x[, .(pct_mean = round(mean(pct_sig), 1),
                              pct_sd = round(sd(pct_sig), 1),
                              n_thresholds = uniqueN(x_threshold),
                              n_cells_mean = round(mean(n_cells), 1)),
                          by = .(approach, contrast, mode, test)]
  setorder(sig_summary, test, contrast, mode, approach)

  fwrite(sig_summary[, .(approach = APPROACH_LAB[as.character(approach)],
                         contrast = contrast_label(contrast),
                         currency = CURRENCY_LAB[as.character(mode)],
                         test = TEST_LAB[as.character(test)],
                         pct_mean, pct_sd, n_thresholds, n_cells_mean)],
         file.path(DIR_APPEND, paste0("TableB2_grouping_pct_significant_", SPATIAL, ".csv")))

  pB2 <- ggplot(sig_summary, aes(x = contrast, y = pct_mean, fill = approach)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(pct_mean - pct_sd, 0),
                      ymax = pmin(pct_mean + pct_sd, 100)),
                  position = position_dodge(width = 0.8), width = 0.25) +
    facet_grid(test ~ mode, labeller = labeller(test = TEST_LAB, mode = CURRENCY_LAB)) +
    scale_x_discrete(labels = CONTRAST_LAB) +
    scale_fill_manual(values = APPROACH_PAL, labels = APPROACH_LAB, name = NULL) +
    labs(title = paste0("Share of testable cells significant (p < ", ALPHA,
                        ") by prey-grouping approach, ", LEVEL_SHORT[[SPATIAL]], " scale"),
         subtitle = paste0("Mean +/- SD across ", length(sig_x),
                           " common taxonomic resolutions"),
         x = NULL, y = "Cells significant (%)") +
    theme_diag(base_size = 10)
  save_app(pB2, paste0("FigB2_grouping_pct_significant_", SPATIAL), 10, 11)

  # ===========================================================================
  # B1. FAMILY COMPOSITION (Fig5 grammar) : 4 families x 3 approaches
  # ===========================================================================
  db <- both[testable == TRUE & x_threshold %in% sig_x &
               !is.na(diagnostic) & diagnostic != INCONCLUSIVE_LAB &
               !is.na(contrast)]
  db <- merge(db, FAM_MAP, by = "diagnostic", all.x = TRUE)
  if (anyNA(db$family))
    warning("Diagnostics outside STATE_TO_FAMILY: ",
            paste(unique(db[is.na(family), diagnostic]), collapse = ", "))
  db <- db[!is.na(family)]

  # Family percentages within each resolution on the full family grid (absent
  # family = 0), then mean and SD across resolutions.
  keys <- c("approach", "mode", "contrast", "x_threshold")
  tot  <- db[, .(n_tot = .N), by = keys]
  cnt  <- db[, .(n = .N), by = c(keys, "family")]
  grid <- tot[, CJ(family = FAMILY_LEVELS), by = keys]
  fam_by_x <- merge(grid, cnt, by = c(keys, "family"), all.x = TRUE)
  fam_by_x <- merge(fam_by_x, tot, by = keys)
  fam_by_x[is.na(n), n := 0L]
  fam_by_x[, pct := 100 * n / n_tot]

  fam_summary <- fam_by_x[, .(pct_mean = round(mean(pct), 1),
                              pct_sd   = round(sd(pct), 1),
                              n_thresholds = uniqueN(x_threshold)),
                          by = .(approach, mode, contrast, family)]
  fam_summary[, family := factor(family, levels = FAMILY_LEVELS)]
  setorder(fam_summary, approach, mode, contrast, family)

  out_long <- fam_summary[, .(approach = APPROACH_LAB[as.character(approach)],
                              currency = CURRENCY_LAB[as.character(mode)],
                              contrast = contrast_label(contrast),
                              family, pct_mean, pct_sd, n_thresholds)]
  fwrite(out_long, file.path(DIR_APPEND, paste0("TableB1_grouping_families_", SPATIAL, ".csv")))

  fam_wide <- dcast(out_long, approach + currency + contrast ~ family, value.var = "pct_mean")
  setcolorder(fam_wide, c("approach", "currency", "contrast", FAMILY_LEVELS))
  fwrite(fam_wide, file.path(DIR_APPEND, paste0("TableB1_grouping_families_wide_", SPATIAL, ".csv")))

  # Figure: one row of panels per family, the approaches in colour.
  dplot <- copy(fam_summary)[, x := as.integer(contrast)]
  present  <- levels(droplevels(dplot$contrast))
  n_within <- sum(INTRA_CONTRASTS %in% present)
  n_total  <- length(CONTRAST_LEVELS)

  pB1 <- ggplot(dplot, aes(x = x, y = pct_mean, colour = approach, group = approach)) +
    period_bands(n_within, n_total) +
    geom_errorbar(aes(ymin = pmax(pct_mean - pct_sd, 0), ymax = pct_mean + pct_sd),
                  width = 0.12, linewidth = 0.35, alpha = 0.8) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = seq_len(n_total), labels = unname(CONTRAST_LAB),
                       limits = c(0.5, n_total + 0.5), expand = c(0, 0)) +
    scale_colour_manual(values = APPROACH_PAL, labels = APPROACH_LAB, name = NULL) +
    facet_grid(family ~ mode, scales = "free_y",
               labeller = labeller(mode = CURRENCY_LAB)) +
    labs(title = paste0("Diagnostic families by prey-grouping approach, ",
                        LEVEL_SHORT[[SPATIAL]], " scale"),
         subtitle = paste0("Mean +/- SD across ", length(sig_x),
                           " common taxonomic resolutions.\nShaded bands: ",
                           "within-period (left) and between-period (right) contrasts."),
         x = NULL, y = "Frequency (%)") +
    theme_diag(base_size = 10)
  save_app(pB1, paste0("FigB1_grouping_families_", SPATIAL), 9, 10)

  message("Written for ", SPATIAL, ": FigB1, FigB2, TableB1 (long + wide), TableB2")
}

cat("\nDone. Appendix outputs in ", normalizePath(DIR_APPEND), "\n", sep = "")
