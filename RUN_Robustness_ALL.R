# =============================================================================
# 00_RUN_Robustness_ALL.R
# -----------------------------------------------------------------------------
# A lancer APRES que 6d (stratum, RUN_TAG = "_T") soit termine.
#
# Ce fichier n'est PAS a sourcer d'un bloc. Chaque BLOC se lance separement,
# dans une session R neuve, parce que 6i et 6k chargent le moteur, qui fige le
# niveau spatial au moment du source().
#
#   BLOC 0  verification : tout est-il bien sur le disque
#   BLOC 1  6h aux trois niveaux, sur les runs patches (_T)      ~30 min
#   BLOC 2  comparaison avant / apres patch (non-regression)     ~1 min
#   BLOC 3  6k calibration sous permutation                      heures
#   BLOC 4  6i + 6j effort equilibre en traits                   heures
#   BLOC 5  recapitulatif : les sept chiffres a lire
#
# Dans RStudio : selectionner le bloc, Ctrl+Entree. Ctrl+Shift+F10 entre les
# blocs 3, 4 et tout ce qui precede.
# =============================================================================


# =============================================================================
# BLOC 0 - VERIFICATION      (session neuve)
# =============================================================================

setwd("C:/Users/SOSTHENEA/Desktop/Diet_analysis")
source("R_helpers/Load_all_runs.R")

# Les dix runs de chaque niveau sont-ils la ?
for (lvl in c("all_gulf", "ecoregion", "stratum")) {
  d <- sprintf("Sensitivity_%s_1_T", lvl)
  n <- if (dir.exists(d)) length(list.files(d, pattern = "\\.rda$")) else 0L
  cat(sprintf("%-28s %2d / 10 runs%s\n", d, n,
              if (n == 10L) "  OK" else "  INCOMPLET -> relancer ce niveau"))
}

# Les .rda groupes (absents = pas grave, le chargeur retombe sur le dossier)
print(inventory_runs()$pooled)

# Lecture reelle des trois niveaux + presence de se_Bs (le patch)
for (lvl in c("all_gulf", "ecoregion", "stratum")) {
  r <- try(load_all_runs(lvl, "_1", tag = "_T", verbose = FALSE), silent = TRUE)
  if (inherits(r, "try-error")) { cat(lvl, ": ECHEC DE LECTURE\n"); next }
  cat(sprintf("%-10s %6d lignes | %3d cellules | se_Bs rempli : %s\n",
              lvl, nrow(r), dplyr::n_distinct(r$cell),
              if ("se_Bs" %in% names(r) && any(is.finite(r$se_Bs))) "OUI" else "NON"))
}

# Ne continuer que si les trois lignes affichent se_Bs = OUI.


# =============================================================================
# BLOC 1 - 6h AUX TROIS NIVEAUX, SUR LES RUNS PATCHES
# =============================================================================
# Ne relance pas le moteur. Produit M0 a M6 dans Robustness_<niveau>_1_T/.
# Environ 30 min au total, le stratum etant le plus long.

rm(list = ls())
setwd("C:/Users/SOSTHENEA/Desktop/Diet_analysis")

for (lvl in c("all_gulf", "ecoregion", "stratum")) {
  PREY_FAMILY <- "_1"
  SPATIAL_LEVEL <- lvl
  RUN_TAG <- "_T"
  B_BOOT <- 2000
  cat("\n\n########## 6h :", lvl, "##########\n")
  source("6h_Robustness_PostHoc.R")
}


# =============================================================================
# BLOC 2 - NON-REGRESSION : AVANT vs APRES PATCH
# =============================================================================
# Le patch n'ajoute que trois colonnes ; aucun pourcentage ne doit bouger.
# Si un ecart depasse 0.01 point, il faut comprendre pourquoi AVANT de publier.

library(dplyr); library(readr)

verif <- lapply(c("all_gulf", "ecoregion", "stratum"), function(lvl) {
  f_old <- sprintf("Robustness_%s_1/M0_family_frequencies_reference.csv", lvl)
  f_new <- sprintf("Robustness_%s_1_T/M0_family_frequencies_reference.csv", lvl)
  if (!file.exists(f_old) || !file.exists(f_new)) {
    cat(lvl, ": M0 manquant d'un cote, comparaison impossible\n"); return(NULL)
  }
  full_join(read_csv(f_old, show_col_types = FALSE) %>% rename(pct_old = pct_mean),
            read_csv(f_new, show_col_types = FALSE) %>% rename(pct_new = pct_mean),
            by = c("contrast", "contrast_type", "mode", "family")) %>%
    mutate(level = lvl, ecart = abs(pct_new - pct_old))
}) %>% bind_rows()

cat("\nEcart maximal entre les deux versions :", round(max(verif$ecart, na.rm = TRUE), 6), "\n")
print(verif %>% filter(ecart > 0.01) %>% select(level, contrast, mode, family, pct_old, pct_new))
write_csv(verif, "VERIF_non_regression_patch.csv")

# Et le gain : pct_equiv_Bs passe de NA a une valeur
print(read_csv("Robustness_all_gulf_1_T/M1_TOST_substitution.csv", show_col_types = FALSE) %>%
        select(contrast, mode, n_substitution, pct_equiv_H, pct_equiv_Bs, pct_equiv_both))


# =============================================================================
# BLOC 3 - 6k : CALIBRATION SOUS PERMUTATION      (session neuve obligatoire)
# =============================================================================
# Permute les etiquettes de periode a l'interieur de chaque cellule : sous cette
# permutation il n'y a rien a trouver, donc le taux de rejet EST le taux de faux
# positifs de chaque test, a l'effectif reel des cellules.
#
# ETAPE 3a : chronometrer. Editer 6k_Null_Calibration.R, mettre B_PERM <- 2,
#            lancer, lire le temps affiche, multiplier par 100.
# ETAPE 3b : remettre B_PERM <- 200 et lancer pour de bon.

rm(list = ls())
setwd("C:/Users/SOSTHENEA/Desktop/Diet_analysis")
PREY_FAMILY <- "_1"
SPATIAL_LEVEL <- "all_gulf"
source("6k_Null_Calibration.R")

# Le niveau Gulf-wide suffit pour l'argument. Pour les autres niveaux, session
# neuve a chaque fois :
#   rm(list = ls()); SPATIAL_LEVEL <- "ecoregion"; source("6k_Null_Calibration.R")


# =============================================================================
# BLOC 4 - 6i + 6j : EFFORT EQUILIBRE EN TRAITS   (session neuve obligatoire)
# =============================================================================
# Plafonne chaque contraste au meme nombre de traits par periode, cellule par
# cellule, et repete B fois. Repond a : la separation intra/inter vient-elle de
# l'ecologie ou du fait que le contraste decennal a plus de traits ?
#
# ETAPE 4a : chronometrer avec B_REP <- 3 dans le fichier.
#            Si le total projete depasse une nuit, reduire X_GRID a
#            seq(100, 1000, by = 200) AVANT de toucher a B_REP.
# ETAPE 4b : B_REP <- 300, puis lancer.

rm(list = ls())
setwd("C:/Users/SOSTHENEA/Desktop/Diet_analysis")
PREY_FAMILY <- "_1"
SPATIAL_LEVEL <- "all_gulf"
source("6i_Balanced_Resampling.R")

# LIRE resample_design_all_gulf_1.csv AVANT d'aller plus loin : si l'equilibrage
# ecarte plus d'un quart des cellules, il faudra le declarer dans le texte.

source("6j_Summarise_Balanced.R")


# =============================================================================
# BLOC 5 - RECAPITULATIF : LES SEPT CHIFFRES A LIRE
# =============================================================================

rm(list = ls())
setwd("C:/Users/SOSTHENEA/Desktop/Diet_analysis")
library(dplyr); library(readr)

show <- function(path, label, cols = NULL) {
  if (!file.exists(path)) { cat("\n[", label, "] absent\n"); return(invisible()) }
  d <- read_csv(path, show_col_types = FALSE)
  if (!is.null(cols)) d <- d %>% select(any_of(cols))
  cat("\n=====", label, "=====\n"); print(as.data.frame(d))
}

for (lvl in c("all_gulf", "ecoregion", "stratum")) {
  R <- sprintf("Robustness_%s_1_T", lvl)
  cat("\n\n##################", lvl, "##################\n")

  # 1. les pourcentages reproduisent-ils le manuscrit
  show(file.path(R, "M0_family_frequencies_reference.csv"),
       "M0 - frequences de familles")

  # 2. equivalence par cellule (attendu : faible)
  show(file.path(R, "M1_TOST_substitution.csv"), "M1 - TOST par cellule")

  # 3. equivalence agregee (LE chiffre qui rattrape M1)
  show(file.path(R, "M6_aggregate_equivalence.csv"),
       "M6 - equivalence agregee",
       c("contrast", "mode", "family", "n_cells", "mean_dH", "dH_lo", "dH_hi",
         "equiv_H_at_margin", "mean_dBs", "dBs_lo", "dBs_hi", "equiv_Bs_at_margin"))

  # 4. l'ordre des familles est-il du signal
  show(file.path(R, "M6b_family_ordering.csv"), "M6b - ordre des familles")

  # 5. la multiplicite detruit-elle le resultat
  show(file.path(R, "M2_FDR_family_frequencies.csv"), "M2 - apres BH")

  # 6. les bandes intra et inter se recouvrent-elles
  show(file.path(R, "M3_CI_overlap_verdict.csv"), "M3 - verdict bootstrap")

  # 7. la classification tient-elle sur les resolutions
  show(file.path(R, "M5_stability_summary.csv"), "M5 - stabilite inter-seuils")
}

# 8. calibration des trois tests
show("NullCalibration_all_gulf_1/null_false_positive_rates.csv",
     "6k - taux de faux positifs (nominal 5 %)")
show("NullCalibration_all_gulf_1/null_family_distribution.csv",
     "6k - distribution NULLE des familles")

# 9. effort equilibre
show("Robustness_all_gulf_1_BAL/balanced_separation_verdict.csv",
     "6j - verdict effort equilibre")

# =============================================================================
# COMMENT LIRE
# =============================================================================
# M0        doit reproduire le manuscrit. Sinon, tout s'arrete la.
# M6        dH_lo et dH_hi dans +/-0.20 pour Substitution = tu peux ecrire
#           l'equivalence de la diversite. C'est le chiffre qui sauve M1.
# M6b       ordered = TRUE : Substitution a bien un |delta H'| plus faible que
#           Reorganisation, et ce n'est pas du bruit.
# M2        shift_pp positif sur Substitution = la correction te renforce.
# M3        disjoint = TRUE ET gap_pp large = separation solide. Un gap_pp
#           sous 2 points signifie que les bandes se touchent : ne revendique
#           pas ce contraste-la.
# 6k        les trois taux proches de 5 % = typologie calibree. Un ecart entre
#           tests = Substitution partiellement structurelle, a declarer.
# 6j        pct_replicates_separated proche de 100 = la separation ne vient pas
#           de l'asymetrie d'effort.
# =============================================================================
