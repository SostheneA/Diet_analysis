# =============================================================================
# PATCH_6a_export_Bs_stats.R
# -----------------------------------------------------------------------------
# UNE seule modification demandee dans 6a_Engine_Trophic.R, en deux endroits.
#
# POURQUOI
#   run_pipeline() ecrit t_H et df_H pour la diversite, donc l'erreur type de
#   delta_H se reconstruit par |delta_H| / |t_H|. Pour l'amplitude de niche,
#   compare_index_sets() ne renvoie que delta et p : ni t, ni erreur type. Le
#   test d'equivalence (TOST) sur Bs est donc impossible, alors que le test de
#   Welch qui produit p calcule deja l'erreur type et la jette.
#
#   Sans ce patch, le manuscrit ne peut affirmer que "H' et Bs n'ont pas ete
#   detectes comme differents" dans les cellules Substitution. Avec, il peut
#   affirmer qu'ils sont equivalents a une marge declaree pres -- ce qui est
#   l'argument central du papier.
#
# EFFET SUR L'EXISTANT : aucun. Deux colonnes s'ajoutent aux resultats, rien
# n'est renomme ni supprime, les scripts 7 a 9 ne bougent pas.
# =============================================================================


# -----------------------------------------------------------------------------
# MODIFICATION 1 -- remplacer integralement compare_index_sets() (ligne ~390)
# -----------------------------------------------------------------------------
# Trois ajouts, signales par ##### : t et se dans la liste de sortie, capture de
# l'objet t.test dans la branche welch, et ecriture des deux valeurs.

compare_index_sets <- function(v1, v2, test = BS_TEST, n_perm = R_PERM) {
  v1 <- v1[is.finite(v1)]; v2 <- v2[is.finite(v2)]
  out <- list(mean_P1 = NA_real_, mean_P2 = NA_real_,
              med_P1 = NA_real_, med_P2 = NA_real_,
              delta = NA_real_, p = NA_real_,
              t = NA_real_, df = NA_real_, se = NA_real_,          ##### AJOUT
              test = test, n1 = length(v1), n2 = length(v2))
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
        tt <- suppressWarnings(stats::t.test(v2, v1, var.equal = FALSE))
        out$t  <- round(unname(tt$statistic), 3)                   ##### AJOUT
        out$df <- round(unname(tt$parameter), 1)                   ##### AJOUT
        # erreur type de la difference des moyennes, telle qu'utilisee par Welch
        out$se <- round(unname(tt$stderr), 5)                      ##### AJOUT
        tt$p.value
      }
    },
    wilcox = suppressWarnings(stats::wilcox.test(v1, v2, exact = FALSE)$p.value),
    perm = {
      allv <- c(v1, v2); g <- rep(c(1L, 2L), c(length(v1), length(v2)))
      obs <- mean(allv[g == 2L]) - mean(allv[g == 1L])
      null <- replicate(n_perm, {
        gg <- sample(g); mean(allv[gg == 2L]) - mean(allv[gg == 1L])
      })
      # erreur type empirique, pour que le TOST reste possible hors welch
      out$se <- round(stats::sd(null), 5)                          ##### AJOUT
      (sum(abs(null) >= abs(obs)) + 1) / (n_perm + 1)
    },
    stop("Unknown BS_TEST: ", test)
  )

  out$p <- round(as.numeric(p), 4)
  out
}


# -----------------------------------------------------------------------------
# MODIFICATION 2 -- dans run_pipeline(), tibble `row` (ligne ~669)
# -----------------------------------------------------------------------------
# Remplacer cette ligne :
#
#   delta_Bs = Bscmp$delta, p_Bs = Bscmp$p, Bs_test = Bscmp$test, n_cat = n_cat,
#
# par celle-ci :
#
#   delta_Bs = Bscmp$delta, p_Bs = Bscmp$p,
#   t_Bs = Bscmp$t, df_Bs = Bscmp$df, se_Bs = Bscmp$se,
#   Bs_test = Bscmp$test, n_cat = n_cat,
#
# -----------------------------------------------------------------------------
# APRES LE PATCH
#   Relancer 6b/6c/6d (et 6e pour _2 et _PP). C'est le seul coup de calcul
#   complet a repayer : environ le temps d'un run standard par niveau.
#   6h_Robustness_PostHoc.R detecte alors t_Bs et remplit pct_equiv_Bs et
#   pct_equiv_both, qui sortent actuellement en NA.
#
#   Si tu ne veux pas relancer tout de suite : le volet H' du TOST fonctionne
#   deja sur les sorties existantes, et c'est le plus important des deux,
#   puisque c'est H' que "Substitution" affirme inchange en premier lieu.
# =============================================================================
