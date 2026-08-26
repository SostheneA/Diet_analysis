# =============================================================================
# 6f_Compare_PP_vs_pooled.R — Comparaison des 3 approches de regroupement
#   (pooled_2 = sans suffixe, pooled_1 = _1, per_predator = _PP, projet frere)
# pour CHAQUE niveau spatial (all_gulf, ecoregion, stratum), selon les
# resultats disponibles. Deux sorties par niveau :
#
#   A. pct_significant_<niveau>.png / .csv
#      % de cellules testables significatives par approche x scenario x mode
#      x test (p_comp, p_H, p_Bs, p_disp), moyenne +/- ecart-type entre seuils.
#
#   B. regime_break_compare_<niveau>.png / .csv   (style Fig5 du manuscrit)
#      Frequence des 4 familles diagnostiques (Stability, Substitution,
#      Functional change, Reorganisation) par contraste (P1a, P1b, P2, PTb,
#      PTa) x currency (biomass/occurrence), moyennee entre les resolutions —
#      UNE facette par approche, pour comparer les 3 en un coup d'oeil.
#
#   C. regime_break_by_family_<niveau>.png
#      Meme contenu, transpose : une LIGNE de facettes par famille diagnostique
#      et, dans chaque panneau, les 3 approches en couleur — le changement du
#      a l'approche se lit directement dans chaque figure.
#
# Tolerant aux campagnes en cours : assemblage depuis les dossiers
# Sensitivity_* si le all_runs consolide manque ; approche/niveau absent saute.
# A lancer depuis le PROJET PRINCIPAL.
# =============================================================================

library(data.table)
library(ggplot2)

MAIN_PROJECT <- normalizePath(here::here())
PP_PROJECT   <- file.path(dirname(MAIN_PROJECT), paste0(basename(MAIN_PROJECT), "_PP"))

ALPHA  <- 0.05
LEVELS <- c("all_gulf", "ecoregion", "stratum")
P_COLS <- c("p_comp", "p_H", "p_Bs", "p_disp")

APPROACHES <- list(
  pooled_2     = list(root = MAIN_PROJECT, suffix = ""),
  pooled_1     = list(root = MAIN_PROJECT, suffix = "_1"),
  per_predator = list(root = PP_PROJECT,   suffix = "_PP")
)

# Scenario (run) -> contraste (label manuscrit, ordre de la Fig5)
SCEN_MAP <- data.table(
  period_1 = c("2004", "2004-2005", "2018", "2006", "2004-2006"),
  period_2 = c("2006", "2006",      "2019", "2018", "2018-2019"),
  scenario = c("P1",   "P3",        "P2",   "P4",   "PT"),
  contrast = c("P1a",  "P1b",       "P2",   "PTb",  "PTa")
)
CONTRAST_ORD <- c("P1a", "P1b", "P2", "PTb", "PTa")

# Diagnostics (9 classes du moteur) -> 4 familles du manuscrit.
# !! A VERIFIER contre add_families() de R_helpers/Config_Mappings.R :
#    regle utilisee ici — Stability = pas de signal ; Functional change =
#    indices (H/Bs) sans composition ; Substitution = composition SANS Bs ;
#    Reorganisation = composition AVEC Bs.
DIAG_FAMILY_MAP <- data.table(
  diagnostic = c("Stable Diet", "Emerging Shift",
                 "Ghost Shift", "Partial Diet Shift",
                 "Internal Rebalancing", "Niche Compression/Expansion",
                 "Niche Restructuring",
                 "Structural Shift", "Major Shift"),
  dfam = c("Stability", "Stability",
           "Substitution", "Substitution",
           "Functional change", "Functional change", "Functional change",
           "Reorganisation", "Reorganisation")
)
FAMILY_ORD <- c("Stability", "Substitution", "Functional change", "Reorganisation")
FAMILY_COL <- c("Stability"         = "#1f3b57",
                "Substitution"      = "#e8804f",
                "Functional change" = "#8aa8b8",
                "Reorganisation"    = "#a01f1f")

dir.create(file.path(MAIN_PROJECT, "data/Sensitivity"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(MAIN_PROJECT, "Sensitivity_Plot_PPcompare"), showWarnings = FALSE)

# --- Chargement ---------------------------------------------------------------

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
                                    src = paste0(d, " (", length(fs), " fichiers)")))
    }
  }
  NULL
}

# =============================================================================
# BOUCLE SUR LES NIVEAUX SPATIAUX
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
      message(sprintf("[ok]     %-13s %s | %d lignes | X = %d-%d",
                      ap, r$src, nrow(dt), min(xs), max(xs)))
    } else message("[absent] ", ap)
  }
  if (!length(res_list)) { message("Aucune approche pour ", SPATIAL, " -> saute"); next }

  both <- rbindlist(res_list, fill = TRUE)

  common_x <- Reduce(intersect, lapply(res_list, function(d) unique(d$x_threshold)))
  sig_x    <- if (length(res_list) >= 2 && length(common_x)) common_x else
    unique(both$x_threshold)
  message(length(sig_x), " seuils (resolutions) utilises")

  # ===========================================================================
  # A. % SIGNIFICATIF par approche x scenario x mode x test
  # ===========================================================================
  long <- melt(both[testable == TRUE & x_threshold %in% sig_x,
                    c("approach", "scenario", "mode", "x_threshold",
                      intersect(P_COLS, names(both))), with = FALSE],
               id.vars = c("approach", "scenario", "mode", "x_threshold"),
               variable.name = "test", value.name = "p")

  pct_by_x <- long[!is.na(p), .(pct_sig = 100 * mean(p < ALPHA), n_cells = .N),
                   by = .(approach, scenario, mode, test, x_threshold)]

  sig_summary <- pct_by_x[, .(pct_mean = round(mean(pct_sig), 1),
                              pct_sd = round(sd(pct_sig), 1),
                              n_thresholds = uniqueN(x_threshold),
                              n_cells_mean = round(mean(n_cells), 1)),
                          by = .(approach, scenario, mode, test)]
  setorder(sig_summary, test, scenario, mode, approach)

  fwrite(sig_summary, file.path(MAIN_PROJECT, "data/Sensitivity",
                                paste0("pct_significant_", SPATIAL, ".csv")))
  fwrite(pct_by_x, file.path(MAIN_PROJECT, "data/Sensitivity",
                             paste0("pct_significant_by_threshold_", SPATIAL, ".csv")))

  pA <- ggplot(sig_summary, aes(x = scenario, y = pct_mean, fill = approach)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(pct_mean - pct_sd, 0),
                      ymax = pmin(pct_mean + pct_sd, 100)),
                  position = position_dodge(width = 0.8), width = 0.25) +
    facet_grid(test ~ mode) +
    labs(title = paste0("% de cellules testables significatives (p < ", ALPHA,
                        ") — moyenne +/- ecart-type entre seuils (", SPATIAL, ")"),
         x = "Scenario (comparaison de periodes)", y = "% significatif", fill = NULL) +
    theme_minimal(base_size = 10)
  ggsave(file.path(MAIN_PROJECT, "Sensitivity_Plot_PPcompare",
                   paste0("pct_significant_", SPATIAL, ".png")),
         pA, width = 11, height = 12, dpi = 200)

  # ===========================================================================
  # B. REGIME BREAK (style Fig5) : 4 familles diagnostiques x 3 approches
  # ===========================================================================
  db <- both[testable == TRUE & x_threshold %in% sig_x &
               !is.na(diagnostic) & diagnostic != "Inconclusive" &
               !is.na(contrast)]
  db <- merge(db, DIAG_FAMILY_MAP, by = "diagnostic", all.x = TRUE)
  if (anyNA(db$dfam))
    warning("Diagnostics hors mapping : ",
            paste(unique(db[is.na(dfam), diagnostic]), collapse = ", "))

  # % de chaque famille AU SEIN de chaque resolution (somme = 100)...
  fam_by_x <- db[!is.na(dfam),
                 .(n = .N), by = .(approach, mode, contrast, x_threshold, dfam)]
  fam_by_x[, pct := 100 * n / sum(n), by = .(approach, mode, contrast, x_threshold)]

  # ...moyenne (+ ecart-type) entre les resolutions
  fam_summary <- fam_by_x[, .(pct_mean = round(mean(pct), 1),
                              pct_sd   = round(sd(pct), 1),
                              n_thresholds = uniqueN(x_threshold)),
                          by = .(approach, mode, contrast, dfam)]
  fam_summary[, `:=`(contrast = factor(contrast, levels = CONTRAST_ORD),
                     dfam     = factor(dfam, levels = FAMILY_ORD))]
  setorder(fam_summary, approach, mode, contrast, dfam)

  fwrite(fam_summary, file.path(MAIN_PROJECT, "data/Sensitivity",
                                paste0("regime_break_compare_", SPATIAL, ".csv")))

  # Table large "manuscrit" : une colonne par famille
  fam_wide <- dcast(fam_summary, approach + mode + contrast ~ dfam,
                    value.var = "pct_mean")
  fwrite(fam_wide, file.path(MAIN_PROJECT, "data/Sensitivity",
                             paste0("regime_break_compare_wide_", SPATIAL, ".csv")))

  # Figure : meme grammaire que la Fig5 (lignes par famille, linetype par
  # currency, bandes intra/inter-decennies), UNE FACETTE PAR APPROCHE.
  dplot <- copy(fam_summary)[, x := as.integer(contrast)]
  n_within <- sum(c("P1a", "P1b", "P2") %in% levels(droplevels(dplot$contrast)))

  pB <- ggplot(dplot, aes(x = x, y = pct_mean,
                          colour = dfam, linetype = mode,
                          group = interaction(dfam, mode))) +
    annotate("rect", xmin = 0.5, xmax = n_within + 0.5, ymin = -Inf, ymax = Inf,
             fill = "#dfe9ee", alpha = 0.5) +
    annotate("rect", xmin = n_within + 0.5, xmax = length(CONTRAST_ORD) + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#f7e2d6", alpha = 0.5) +
    geom_vline(xintercept = n_within + 0.5, linetype = "dashed", colour = "grey55") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = seq_along(CONTRAST_ORD), labels = CONTRAST_ORD,
                       limits = c(0.5, length(CONTRAST_ORD) + 0.5),
                       expand = c(0, 0)) +
    scale_colour_manual(values = FAMILY_COL) +
    facet_wrap(~ approach, ncol = 1) +
    labs(title = paste0("Diet regime break — comparaison des 3 approches (",
                        SPATIAL, ")"),
         subtitle = paste0("Frequence des 4 familles diagnostiques, moyenne sur ",
                           length(sig_x), " resolutions communes ; ",
                           "bande bleue = intra-decennie, orange = inter-decennies"),
         x = NULL, y = "Frequency (%)", colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  ggsave(file.path(MAIN_PROJECT, "Sensitivity_Plot_PPcompare",
                   paste0("regime_break_compare_", SPATIAL, ".png")),
         pB, width = 10, height = 12, dpi = 200)

  # ===========================================================================
  # C. REGIME BREAK "par famille" : une ligne de facettes PAR FAMILLE
  #    diagnostique, et DANS chaque panneau les 3 approches en couleur —
  #    l'effet du choix d'approche se lit directement, famille par famille.
  #    Lignes = approches ; colonnes = currency ; barres fines = +/- ecart-type
  #    entre resolutions.
  # ===========================================================================
  pC <- ggplot(dplot, aes(x = x, y = pct_mean,
                          colour = approach, group = approach)) +
    annotate("rect", xmin = 0.5, xmax = n_within + 0.5, ymin = -Inf, ymax = Inf,
             fill = "#dfe9ee", alpha = 0.5) +
    annotate("rect", xmin = n_within + 0.5, xmax = length(CONTRAST_ORD) + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#f7e2d6", alpha = 0.5) +
    geom_vline(xintercept = n_within + 0.5, linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(ymin = pmax(pct_mean - pct_sd, 0),
                      ymax = pct_mean + pct_sd),
                  width = 0.12, linewidth = 0.35, alpha = 0.8) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = seq_along(CONTRAST_ORD), labels = CONTRAST_ORD,
                       limits = c(0.5, length(CONTRAST_ORD) + 0.5),
                       expand = c(0, 0)) +
    facet_grid(dfam ~ mode, scales = "free_y") +
    labs(title = paste0("Diet regime break par famille diagnostique — ",
                        "effet de l'approche de regroupement (", SPATIAL, ")"),
         subtitle = paste0("Chaque panneau : les 3 approches pour une famille ; ",
                           "moyenne +/- ecart-type sur ", length(sig_x),
                           " resolutions communes ; bande bleue = intra-decennie, ",
                           "orange = inter-decennies"),
         x = NULL, y = "Frequency (%)", colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  ggsave(file.path(MAIN_PROJECT, "Sensitivity_Plot_PPcompare",
                   paste0("regime_break_by_family_", SPATIAL, ".png")),
         pC, width = 10, height = 12, dpi = 200)

  message("Ecrit pour ", SPATIAL,
          " : pct_significant + regime_break_compare + regime_break_by_family")
}

cat("\nTermine. Sorties dans Sensitivity_Plot_PPcompare/ et data/Sensitivity/\n")
