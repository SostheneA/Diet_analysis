# =============================================================================
# 6i_Robustness_PostHoc.R   (v2)
# -----------------------------------------------------------------------------
# Robustesse SANS relance du moteur, assemblee depuis les dossiers Sensitivity_*
# (les sauvegardes all_runs groupees n'ayant jamais ete executees).
#
#   M1  Tailles d'effet par famille + TOST sur les cellules Substitution  [CRITIQUE]
#   M2  Multiplicite (BH-FDR) et recalcul des familles
#   M3  IC bootstrap sur les pourcentages de familles (cellule = unite)
#   M4  Sensibilite a alpha et a la bande de tendance
#   M5  Stabilite de la classification a travers les seuils
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(rlang); library(tibble)
})

PROJECT_DIR <- "C:/Users/SOSTHENEA/Desktop/Diet_analysis"
if (dir.exists(PROJECT_DIR)) setwd(PROJECT_DIR)

if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

ALPHA <- 0.05; TREND_UPPER <- 0.10
B_BOOT <- 2000; SEED <- 20260911
MARGIN_H <- 0.20; MARGIN_BS <- 0.05      # marges d'equivalence, a justifier

OUT <- sprintf("Robustness_%s%s", SPATIAL_LEVEL, PREY_FAMILY)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

source("R_helpers/Load_all_runs.R")
if (file.exists("R_helpers/Config_Mappings.R")) source("R_helpers/Config_Mappings.R")

# =============================================================================
# 0. CHARGEMENT
# =============================================================================

FAM_ORDER <- c("Stability", "Substitution", "Functional change", "Reorganisation")
S2F <- c("Stable Diet" = "Stability", "Emerging Shift" = "Stability",
         "Ghost Shift" = "Substitution",
         "Niche Compression / Expansion" = "Functional change",
         "Internal Rebalancing" = "Functional change",
         "Niche Restructuring" = "Functional change",
         "Partial Diet Shift" = "Reorganisation",
         "Structural Shift" = "Reorganisation", "Major Shift" = "Reorganisation",
         "Inconclusive" = "Inconclusive")

classify_from_p <- function(p_comp, p_H, p_Bs, alpha = ALPHA, trend = TREND_UPPER) {
  cs <- !is.na(p_comp) & p_comp <  alpha
  ct <- !is.na(p_comp) & p_comp >= alpha & p_comp < trend
  hs <- !is.na(p_H)  & p_H  < alpha
  bs <- !is.na(p_Bs) & p_Bs < alpha
  o <- rep(NA_character_, length(p_comp))
  o[ cs & !hs & !bs] <- "Ghost Shift"
  o[ cs &  hs & !bs] <- "Partial Diet Shift"
  o[ cs & !hs &  bs] <- "Structural Shift"
  o[ cs &  hs &  bs] <- "Major Shift"
  o[!cs & !hs &  bs] <- "Niche Compression / Expansion"
  o[!cs &  hs & !bs] <- "Internal Rebalancing"
  o[!cs &  hs &  bs] <- "Niche Restructuring"
  o[!cs & !hs & !bs &  ct] <- "Emerging Shift"
  o[!cs & !hs & !bs & !ct] <- "Stable Diet"
  o[is.na(p_comp) | is.na(p_H) | is.na(p_Bs)] <- "Inconclusive"
  o
}

add_fam <- function(d, alpha = ALPHA, trend = TREND_UPPER,
                    pc = "p_comp", ph = "p_H", pb = "p_Bs") {
  d$diagnostic <- classify_from_p(d[[pc]], d[[ph]], d[[pb]], alpha, trend)
  # mapping canonique si Config_Mappings l'a fourni
  d$family <- if (exists("STATE_TO_FAMILY")) unname(STATE_TO_FAMILY[d$diagnostic]) else
    unname(S2F[d$diagnostic])
  d
}

runs <- load_all_runs(SPATIAL_LEVEL, PREY_FAMILY, materialise = TRUE)
UNIT <- attr(runs, "unit_col")

runs <- runs %>%
  filter(testable %in% c(TRUE, "TRUE")) %>%
  add_fam() %>%
  mutate(cell = paste(species, size_class, .data[[UNIT]], sep = "|"))

cls <- runs %>% filter(family %in% FAM_ORDER)

message(sprintf("%d lignes testables | %d cellules | %d seuils | scenarios: %s",
                nrow(runs), n_distinct(cls$cell), n_distinct(cls$x_threshold),
                paste(sort(unique(cls$scenario)), collapse = ", ")))

# Helper unique : % par famille, moyenne sur unites spatiales puis sur seuils.
pct_fam <- function(d, fam = "family") {
  d %>%
    filter(.data[[fam]] %in% FAM_ORDER) %>%
    count(scenario, contrast_type, mode, x_threshold, .data[[UNIT]],
          family = .data[[fam]], name = "n") %>%
    group_by(scenario, contrast_type, mode, x_threshold, .data[[UNIT]]) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    group_by(scenario, contrast_type, mode, x_threshold, family) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(scenario, contrast_type, mode, family) %>%
    summarise(sd_across_x = sd(pct), pct = mean(pct), .groups = "drop") %>%
    select(scenario, contrast_type, mode, family, pct, sd_across_x)
}

write_csv(pct_fam(cls), file.path(OUT, "M0_family_frequencies_reference.csv"))

# =============================================================================
# M1. TAILLES D'EFFET + EQUIVALENCE
# =============================================================================

if (!"se_delta_H"  %in% names(cls) && "t_H"  %in% names(cls))
  cls$se_delta_H  <- abs(cls$delta_H)  / pmax(abs(cls$t_H),  1e-9)
if (!"se_delta_Bs" %in% names(cls) && "t_Bs" %in% names(cls))
  cls$se_delta_Bs <- abs(cls$delta_Bs) / pmax(abs(cls$t_Bs), 1e-9)

eff_summary <- cls %>%
  group_by(scenario, contrast_type, mode, family) %>%
  summarise(n = n(),
            dH_med  = median(abs(delta_H),  na.rm = TRUE),
            dH_q90  = quantile(abs(delta_H),  .90, na.rm = TRUE),
            dBs_med = median(abs(delta_Bs), na.rm = TRUE),
            dBs_q90 = quantile(abs(delta_Bs), .90, na.rm = TRUE),
            BC_med  = median(BC, na.rm = TRUE),
            R2_med  = median(R2_comp, na.rm = TRUE), .groups = "drop")
write_csv(eff_summary, file.path(OUT, "M1_effect_sizes_by_family.csv"))

if (all(c("se_delta_H", "se_delta_Bs") %in% names(cls))) {
  tost <- cls %>%
    filter(family == "Substitution") %>%
    mutate(lo_H  = delta_H  - 1.645 * se_delta_H,  hi_H  = delta_H  + 1.645 * se_delta_H,
           lo_Bs = delta_Bs - 1.645 * se_delta_Bs, hi_Bs = delta_Bs + 1.645 * se_delta_Bs,
           equiv_H  = lo_H  > -MARGIN_H  & hi_H  < MARGIN_H,
           equiv_Bs = lo_Bs > -MARGIN_BS & hi_Bs < MARGIN_BS,
           equiv_both = equiv_H & equiv_Bs)

  tost_summary <- tost %>%
    group_by(scenario, contrast_type, mode) %>%
    summarise(n_substitution = n(),
              pct_equiv_H    = 100 * mean(equiv_H,    na.rm = TRUE),
              pct_equiv_Bs   = 100 * mean(equiv_Bs,   na.rm = TRUE),
              pct_equiv_both = 100 * mean(equiv_both, na.rm = TRUE), .groups = "drop")
  write_csv(tost_summary, file.path(OUT, "M1_TOST_substitution.csv"))
  print(as.data.frame(tost_summary))

  tost_curve <- expand_grid(margin_H = seq(0.05, 0.50, by = 0.05),
                            margin_Bs = seq(0.02, 0.15, by = 0.02)) %>%
    mutate(pct_equiv = map2_dbl(margin_H, margin_Bs, function(mH, mBs)
      100 * mean(tost$lo_H > -mH & tost$hi_H < mH &
                 tost$lo_Bs > -mBs & tost$hi_Bs < mBs, na.rm = TRUE)))
  write_csv(tost_curve, file.path(OUT, "M1b_TOST_margin_sensitivity.csv"))
} else {
  warning("Ni se_delta_* ni t_* dans les sorties du moteur : ",
          "exporter les erreurs types de Welch depuis 6a. Sans elles, pas de TOST.")
}

p_eff <- cls %>%
  filter(contrast_type == "between") %>%
  select(mode, family, delta_H, delta_Bs) %>%
  mutate(family = factor(family, levels = FAM_ORDER),
         `|delta H'| (nats)` = abs(delta_H), `|delta Bs|` = abs(delta_Bs)) %>%
  pivot_longer(c(`|delta H'| (nats)`, `|delta Bs|`),
               names_to = "index", values_to = "value") %>%
  ggplot(aes(family, value, fill = family)) +
  geom_violin(scale = "width", alpha = .6, colour = NA) +
  geom_boxplot(width = .13, outlier.size = .3) +
  facet_grid(index ~ mode, scales = "free_y") +
  labs(x = NULL, y = "Taille d'effet absolue",
       title = "Tailles d'effet univariees par famille (contrastes inter-periodes)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "M1_effect_sizes.png"), p_eff, width = 9, height = 6, dpi = 300)

# =============================================================================
# M2. MULTIPLICITE (BH-FDR)
# =============================================================================

runs_fdr <- runs %>%
  group_by(scenario, mode, x_threshold) %>%
  mutate(q_comp = p.adjust(p_comp, "BH"),
         q_H    = p.adjust(p_H,    "BH"),
         q_Bs   = p.adjust(p_Bs,   "BH")) %>%
  ungroup() %>%
  add_fam(pc = "q_comp", ph = "q_H", pb = "q_Bs")

fdr_comparison <- pct_fam(cls) %>% select(-sd_across_x) %>% rename(pct_raw = pct) %>%
  left_join(pct_fam(runs_fdr) %>% select(-sd_across_x) %>% rename(pct_bh = pct),
            by = c("scenario", "contrast_type", "mode", "family")) %>%
  mutate(shift_pp = pct_bh - pct_raw)
write_csv(fdr_comparison, file.path(OUT, "M2_FDR_family_frequencies.csv"))

# =============================================================================
# M3. IC BOOTSTRAP  (la cellule est l'unite reechantillonnee)
# =============================================================================
# On tire les cellules avec remise UNE FOIS par replicat, puis on applique le
# meme tirage a tous les seuils : la cellule reste l'unite d'echantillonnage a
# travers le balayage, ce qui est la structure reelle des donnees.

set.seed(SEED)
cell_list <- sort(unique(cls$cell))
K <- length(cell_list)

base_counts <- cls %>%
  count(scenario, contrast_type, mode, x_threshold, cell, family, name = "n_cell")

boot_one <- function(b) {
  w <- tibble(cell = cell_list,
              w = as.integer(rmultinom(1, K, rep(1 / K, K))))
  base_counts %>%
    inner_join(w, by = "cell") %>%
    filter(w > 0) %>%
    count(scenario, contrast_type, mode, x_threshold, family, wt = n_cell * w, name = "n") %>%
    group_by(scenario, contrast_type, mode, x_threshold) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    group_by(scenario, contrast_type, mode, family) %>%
    summarise(pct = mean(pct), .groups = "drop")
}

boot <- map_dfr(seq_len(B_BOOT), boot_one)

boot_ci <- boot %>%
  group_by(scenario, contrast_type, mode, family) %>%
  summarise(pct_mean = mean(pct),
            ci_lo = quantile(pct, .025), ci_hi = quantile(pct, .975), .groups = "drop")
write_csv(boot_ci, file.path(OUT, "M3_bootstrap_CI_family_pct.csv"))

overlap <- boot_ci %>%
  filter(family %in% c("Stability", "Substitution")) %>%
  group_by(mode, family, contrast_type) %>%
  summarise(lo = min(ci_lo), hi = max(ci_hi), .groups = "drop") %>%
  pivot_wider(names_from = contrast_type, values_from = c(lo, hi)) %>%
  mutate(disjoint = if_else(family == "Stability",
                            hi_between < lo_within, lo_between > hi_within),
         gap_pp = if_else(family == "Stability",
                          lo_within - hi_between, lo_between - hi_within))
write_csv(overlap, file.path(OUT, "M3_CI_overlap_verdict.csv"))
print(as.data.frame(overlap))

# =============================================================================
# M4. SENSIBILITE A ALPHA ET A LA BANDE DE TENDANCE
# =============================================================================

alpha_sens <- expand_grid(a = c(.01, .025, .05, .075, .10),
                          tr = c(.05, .10, .15, .20)) %>%
  filter(tr > a) %>%
  pmap_dfr(function(a, tr)
    pct_fam(add_fam(runs, alpha = a, trend = tr)) %>%
      select(-sd_across_x) %>% mutate(alpha = a, trend = tr))
write_csv(alpha_sens, file.path(OUT, "M4_alpha_trend_sensitivity.csv"))

p_alpha <- alpha_sens %>% filter(trend == .10) %>%
  mutate(family = factor(family, levels = FAM_ORDER)) %>%
  ggplot(aes(alpha, pct, colour = contrast_type, group = scenario)) +
  geom_line(alpha = .8) + geom_point(size = 1) +
  facet_grid(mode ~ family) +
  labs(x = "alpha", y = "% des cellules classifiables", colour = NULL,
       title = "Sensibilite des frequences de familles au seuil de significativite") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "M4_alpha_sensitivity.png"), p_alpha, width = 10, height = 5.5, dpi = 300)

# =============================================================================
# M5. STABILITE DE LA CLASSIFICATION SUR LES SEUILS
# =============================================================================

stability_cell <- cls %>%
  group_by(scenario, mode, cell) %>%
  summarise(n_thresholds = n(),
            modal_family = names(which.max(table(family))),
            support      = max(table(family)) / n(),
            n_families   = n_distinct(family), .groups = "drop")
write_csv(stability_cell, file.path(OUT, "M5_cell_classification_stability.csv"))

stability_summary <- stability_cell %>%
  group_by(scenario, mode) %>%
  summarise(pct_support_ge_80 = 100 * mean(support >= .80),
            pct_support_ge_60 = 100 * mean(support >= .60),
            median_support = median(support), .groups = "drop")
write_csv(stability_summary, file.path(OUT, "M5_stability_summary.csv"))
print(as.data.frame(stability_summary))

message("6i termine -> ", OUT)
