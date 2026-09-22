# =============================================================================
# 6m_Effort_Descriptive.R
# -----------------------------------------------------------------------------
# Traitement de l'objection sur l'effort d'echantillonnage SANS relancer le
# moteur. Quelques minutes au lieu de plusieurs heures.
#
# L'OBJECTION
#   Le contraste decennal pool trois annees contre deux. Son plus petit groupe
#   compte plus de traits que celui de tout contraste intra-periode, donc les
#   trois tests y ont plus de puissance. La separation observee pourrait venir
#   du design plutot que de l'ecologie.
#
# TROIS REPONSES, de la plus faible a la plus forte
#   A. L'effort effectivement realise, par contraste. Quantifie l'asymetrie.
#   B. Le gradient d'effort A L'INTERIEUR de chaque contraste. Si Stability
#      etait un artefact de puissance, elle devrait chuter quand le nombre de
#      traits augmente. C'est testable sans rien relancer.
#   C. La restriction au support commun : on ne garde, dans chaque contraste,
#      que les cellules dont l'effort tombe dans la plage couverte par TOUS les
#      contrastes, puis on recalcule les familles. Comparaison a effort
#      comparable, sans reechantillonnage.
#
# C est l'equivalent observationnel de 6i. Plus faible, parce que la restriction
# selectionne des cellules au lieu de les egaliser par tirage, mais elle repond
# a la meme question et coute quelques minutes.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(tibble)
})

if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"
if (!exists("RUN_TAG"))       RUN_TAG       <- "_T"

OUT <- sprintf("Effort_%s%s", SPATIAL_LEVEL, PREY_FAMILY)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

source("R_helpers/Load_all_runs.R")

runs <- load_all_runs(SPATIAL_LEVEL, PREY_FAMILY, tag = RUN_TAG)
UNIT <- attr(runs, "unit_col")

runs <- runs %>%
  filter(testable, !is.na(family)) %>%
  mutate(n_min = pmin(n_set_P1, n_set_P2),
         n_sto_min = pmin(n_sto_P1, n_sto_P2))

# Agregation identique a 6h : unite spatiale, zeros inclus, puis seuils.
pct_fam <- function(d) {
  d %>%
    count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT),
          family, name = "n") %>%
    group_by(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT)) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    tidyr::complete(
      tidyr::nesting(contrast, contrast_type, mode, x_threshold,
                     !!rlang::sym(UNIT)),
      family = FAMILY_LEVELS, fill = list(n = 0L, pct = 0)) %>%
    group_by(contrast, contrast_type, mode, x_threshold, family) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(contrast, contrast_type, mode, family) %>%
    summarise(pct = mean(pct), .groups = "drop")
}

# =============================================================================
# A. L'EFFORT REALISE, PAR CONTRASTE
# =============================================================================
# Une cellule apparait a 100 resolutions ; on la compte une fois, a la mediane
# de ses valeurs, sinon les cellules bien echantillonnees pesent autant de fois
# qu'il y a de resolutions.

per_cell <- runs %>%
  group_by(contrast, contrast_type, mode, cell) %>%
  summarise(n_set_P1 = median(n_set_P1), n_set_P2 = median(n_set_P2),
            n_min = median(n_min), n_sto_min = median(n_sto_min),
            .groups = "drop")

effort <- per_cell %>%
  group_by(contrast, contrast_type, mode) %>%
  summarise(
    n_cells    = n(),
    sets_P1_med = median(n_set_P1),
    sets_P2_med = median(n_set_P2),
    sets_min_med = median(n_min),
    sets_min_q25 = quantile(n_min, .25),
    sets_min_q75 = quantile(n_min, .75),
    sto_min_med  = median(n_sto_min),
    .groups = "drop")

write_csv(effort, file.path(OUT, "A_effort_by_contrast.csv"))
cat("\n=== A. Effort realise par contraste (traits par periode) ===\n")
print(as.data.frame(effort %>%
  select(contrast, contrast_type, mode, n_cells,
         sets_P1_med, sets_P2_med, sets_min_med, sets_min_q25, sets_min_q75)))

# =============================================================================
# B. GRADIENT D'EFFORT A L'INTERIEUR DE CHAQUE CONTRASTE
# =============================================================================
# Si Stability n'etait qu'un deficit de puissance, elle devrait decroitre quand
# le nombre de traits augmente, A L'INTERIEUR d'un meme contraste. On le mesure.

gradient <- runs %>%
  mutate(n_bin = cut(n_min, breaks = c(2, 4, 6, 9, 14, 25, Inf), right = FALSE)) %>%
  filter(!is.na(n_bin)) %>%
  count(contrast, contrast_type, mode, n_bin, family, name = "n") %>%
  group_by(contrast, contrast_type, mode, n_bin) %>%
  mutate(pct = 100 * n / sum(n), n_tot = sum(n)) %>%
  ungroup()

write_csv(gradient, file.path(OUT, "B_family_by_effort_class.csv"))

# Pente : regression logistique de "cellule = Stability" sur log(n traits),
# ajustee separement dans chaque contraste et devise.
slopes <- runs %>%
  mutate(is_stab = as.integer(family == "Stability")) %>%
  group_by(contrast, contrast_type, mode) %>%
  group_modify(function(d, k) {
    if (n_distinct(d$is_stab) < 2 || n_distinct(d$n_min) < 3) return(tibble())
    fit <- try(glm(is_stab ~ log(n_min), family = binomial, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(tibble())
    s <- summary(fit)$coefficients
    tibble(slope_log_nsets = s[2, 1], se = s[2, 2],
           p = s[2, 4], n_obs = nrow(d))
  }) %>% ungroup()

write_csv(slopes, file.path(OUT, "B_stability_vs_effort_slope.csv"))
cat("\n=== B. Pente de Stability en fonction de log(traits) ===\n")
cat("Une pente NEGATIVE et significative indiquerait que Stability decroit\n")
cat("quand l'effort augmente, donc un effet de puissance.\n\n")
print(as.data.frame(slopes))

p_grad <- gradient %>%
  filter(family %in% c("Stability", "Substitution")) %>%
  ggplot(aes(n_bin, pct, colour = contrast, group = contrast,
             linetype = contrast_type)) +
  geom_line() + geom_point(size = 1.2) +
  facet_grid(mode ~ family) +
  labs(x = "Traits par periode (minimum des deux)",
       y = "% des cellules classifiables", colour = NULL, linetype = NULL,
       title = "Frequence des familles en fonction de l'effort, dans chaque contraste",
       subtitle = "Des courbes plates indiquent que l'effort ne pilote pas la classification") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "B_effort_gradient.png"), p_grad,
       width = 10, height = 6, dpi = 300)

# =============================================================================
# C. RESTRICTION AU SUPPORT COMMUN
# =============================================================================
# On ne garde que les cellules dont l'effort tombe dans la plage couverte par
# TOUS les contrastes. Les familles sont alors comparees a effort comparable.

# Le support commun est l'intersection des etendues d'effort des cinq
# contrastes. On elargit progressivement les centiles jusqu'a obtenir un
# intervalle non vide : d'abord [10e, 90e], puis [5e, 95e], puis l'etendue
# complete. Un intervalle toujours vide signifie que les contrastes ne
# partagent aucune plage d'effort, ce qui est en soi le resultat a rapporter.

common_support <- function(q_lo, q_hi) {
  per_cell %>%
    group_by(contrast, mode) %>%
    summarise(lo = quantile(n_min, q_lo), hi = quantile(n_min, q_hi),
              .groups = "drop") %>%
    group_by(mode) %>%
    summarise(lo_common = max(lo), hi_common = min(hi), .groups = "drop")
}

support <- NULL
for (qq in list(c(.10, .90), c(.05, .95), c(0, 1))) {
  cand <- common_support(qq[1], qq[2])
  if (all(cand$lo_common <= cand$hi_common)) {
    support <- cand %>% mutate(q_lo = qq[1], q_hi = qq[2]); break
  }
}

cat("\n=== C. Support commun en traits par periode ===\n")
if (is.null(support)) {
  cat("AUCUN support commun : les contrastes ne partagent aucune plage d'effort.\n")
  cat("A rapporter tel quel ; la restriction a effort comparable est impossible\n")
  cat("et seul un reechantillonnage (6i) peut trancher.\n")
  runs_common <- runs[0, ]
} else {
  print(as.data.frame(support))
  cat(sprintf("Centiles utilises : %.0f%% - %.0f%%\n",
              100 * support$q_lo[1], 100 * support$q_hi[1]))
  runs_common <- runs %>%
    left_join(support %>% select(mode, lo_common, hi_common), by = "mode") %>%
    filter(n_min >= lo_common, n_min <= hi_common)
}

if (!nrow(runs_common)) {
  message("Module C non concluant : aucune cellule dans le support commun.")
}

n_kept <- runs_common %>%
  group_by(contrast, mode) %>%
  summarise(cells_kept = n_distinct(cell), .groups = "drop") %>%
  left_join(per_cell %>% count(contrast, mode, name = "cells_total"),
            by = c("contrast", "mode")) %>%
  mutate(pct_kept = round(100 * cells_kept / cells_total, 1))
write_csv(n_kept, file.path(OUT, "C_cells_retained.csv"))
cat("\nCellules retenues dans le support commun :\n")
print(as.data.frame(n_kept))

comparison <- pct_fam(runs) %>% rename(pct_all = pct) %>%
  left_join(pct_fam(runs_common) %>% rename(pct_common = pct),
            by = c("contrast", "contrast_type", "mode", "family")) %>%
  mutate(pct_common = coalesce(pct_common, 0),
         shift_pp = pct_common - pct_all,
         family = factor(family, levels = FAMILY_LEVELS)) %>%
  arrange(mode, contrast, family)
write_csv(comparison, file.path(OUT, "C_common_support_vs_all.csv"))

verdict <- comparison %>%
  filter(family %in% c("Stability", "Substitution")) %>%
  group_by(mode, family, contrast_type) %>%
  summarise(lo_all = min(pct_all), hi_all = max(pct_all),
            lo_com = min(pct_common), hi_com = max(pct_common), .groups = "drop") %>%
  pivot_wider(names_from = contrast_type,
              values_from = c(lo_all, hi_all, lo_com, hi_com)) %>%
  mutate(
    gap_all = if_else(family == "Stability",
                      lo_all_within - hi_all_between, lo_all_between - hi_all_within),
    gap_common = if_else(family == "Stability",
                         lo_com_within - hi_com_between, lo_com_between - hi_com_within))
write_csv(verdict, file.path(OUT, "C_separation_verdict.csv"))

cat("\n=== VERDICT : ecart entre bandes intra et inter (points) ===\n")
print(as.data.frame(verdict %>% select(mode, family, gap_all, gap_common)))
cat("\nSi gap_common reste large et positif, la separation ne vient pas de\n")
cat("l'asymetrie d'effort entre contrastes.\n")

message("6m termine -> ", OUT)
