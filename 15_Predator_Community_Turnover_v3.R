## =====================================================================
## 15_Predator_Community_Turnover_v3.R
## Spatial community beta-diversity of predators (Volet 3)
## Partitioning total dissimilarity into turnover (SIM) vs nestedness (SNE)
##
## Project: Diet_analysis | Input: data/dat_classed.rda
## =====================================================================

## ---------------------------------------------------------------------
## 0. CONFIG
## ---------------------------------------------------------------------

if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"
DAT_PATH <- file.path("data", "dat_classed.rda")

COL_PRED    <- "predator_species_common_name"
COL_PERIOD  <- "period"
COL_YEAR    <- "year"
COL_VESSEL  <- "vessel.code"
COL_SETNO   <- "set"
COL_AREA    <- "Area"
COL_STRATUM <- "stratum"

CONTRASTS <- list(
  PT = list(p1 = "2004-2006", p2 = "2018-2019")
)

OUT_DIR  <- sprintf("CommunityTurnover_%s", SPATIAL_LEVEL)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(betapart)
})

msg <- function(...) cat(sprintf(...), "\n")

## ---------------------------------------------------------------------
## 1. CHARGEMENT ET PRÉPARATION DES DONNÉES DE COMMUNAUTÉ
## ---------------------------------------------------------------------

stopifnot(file.exists(DAT_PATH))
.obj <- load(DAT_PATH); dat <- get(.obj[1])

if (inherits(dat, "sf")) {
  dat <- sf::st_drop_geometry(dat)
  msg("Géométrie sf retirée pour accélérer le traitement.")
}
dat <- as.data.frame(dat)

# Clé unique pour chaque station de chalutage (trawl set)
dat$.set    <- paste(dat[[COL_YEAR]], dat[[COL_VESSEL]], dat[[COL_SETNO]], sep = "_")
dat$.period <- as.character(dat[[COL_PERIOD]])
dat$.pred   <- as.character(dat[[COL_PRED]])
dat$.unit   <- if (SPATIAL_LEVEL == "all_gulf") "all_gulf" else as.character(dat[[SPATIAL_LEVEL]])

msg("Données chargées : %d lignes, %d stations uniques.", nrow(dat), length(unique(dat$.set)))

## ---------------------------------------------------------------------
## 2. FONCTION DE DÉCOMPOSITION SPATIALE (Baselga / Sørensen)
## ---------------------------------------------------------------------

run_community_decomposition <- function(d_unit, p1_lab, p2_lab) {
  # Filtrer pour les deux périodes du contraste
  d_sub <- d_unit[d_unit$.period %in% c(p1_lab, p2_lab), , drop = FALSE]
  if (!nrow(d_sub)) return(NULL)

  # Créer la matrice station x espèce de prédateur (Présence / Absence)
  comm_mat <- d_sub %>%
    distinct(.set, .period, .pred) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = .pred, values_from = present, values_fill = 0)

  sets_info <- comm_mat[, c(".set", ".period")]
  species_mat <- as.matrix(comm_mat[, setdiff(names(comm_mat), c(".set", ".period"))])

  # Séparer les sites par période pour calculer la bêta-diversité inter-périodes
  # (On compare le vecteur moyen ou l'incidence globale de la communauté entre p1 et p2)
  idx_p1 <- sets_info$.period == p1_lab
  idx_p2 <- sets_info$.period == p2_lab

  if (sum(idx_p1) < 2 || sum(idx_p2) < 2) {
    warning("Pas assez de stations pour comparer les périodes.")
    return(NULL)
  }

  # Incidence globale (présence de l'espèce dans au moins un set de la période)
  inc_p1 <- as.numeric(colSums(species_mat[idx_p1, , drop = FALSE]) > 0)
  inc_p2 <- as.numeric(colSums(species_mat[idx_p2, , drop = FALSE]) > 0)

  mat_comparison <- rbind(p1 = inc_p1, p2 = inc_p2)
  colnames(mat_comparison) <- colnames(species_mat)

  # Nettoyer les colonnes ne contenant que des zéros sur les deux périodes
  keep_cols <- colSums(mat_comparison) > 0
  mat_comparison <- mat_comparison[, keep_cols, drop = FALSE]

  if (ncol(mat_comparison) < 2) return(NULL)

  # Calcul de la décomposition de Baselga (Indice de Sørensen)
  # beta.sor = Dissimilarité totale
  # beta.sim = Turnover pur (remplacement d'espèces)
  # beta.sne = Imbrication (nestedness / perte-gain net)
  bp <- betapart::beta.pair(mat_comparison, index.family = "sorensen")

  data.frame(
    beta_sor = as.numeric(bp$beta.sor),
    beta_sim = as.numeric(bp$beta.sim),
    beta_sne = as.numeric(bp$beta.sne),
    n_species_p1 = sum(inc_p1),
    n_species_p2 = sum(inc_p2),
    n_shared_species = sum(inc_p1 == 1 & inc_p2 == 1),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------
## 3. EXÉCUTION DU SWEEP SPATIAL ET SAUVEGARDE
## ---------------------------------------------------------------------

results_list <- list()

for (ctr_name in names(CONTRASTS)) {
  ctr <- CONTRASTS[[ctr_name]]
  for (u in sort(unique(dat$.unit))) {
    d_unit <- dat[dat$.unit == u, , drop = FALSE]
    res <- try(run_community_decomposition(d_unit, ctr$p1, ctr$p2), silent = TRUE)

    if (inherits(res, "try-error") || is.null(res)) {
      msg("  ! Erreur ou données insuffisantes pour l'unité spatiale : %s", u)
      next
    }

    row_out <- data.frame(
      contrast = ctr_name,
      spatial_level = SPATIAL_LEVEL,
      unit = u,
      res,
      stringsAsFactors = FALSE
    )
    results_list[[length(results_list) + 1]] <- row_out
  }
}

if (length(results_list)) {
  community_tab <- bind_rows(results_list)
  out_file <- file.path(OUT_DIR, "predator_community_turnover_baselga.csv")
  write.csv(community_tab, out_file, row.names = FALSE)

  msg("\n=== Résultats du Volet 3 (Turnover des Prédateurs) ===")
  print(community_tab)
  msg("\nFichier sauvegardé avec succès dans : %s", out_file)
} else {
  warning("Aucun résultat généré pour le Volet 3.")
}
