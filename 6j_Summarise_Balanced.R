# =============================================================================
# 6j_Summarise_Balanced.R
# -----------------------------------------------------------------------------
# Agrege la sortie de 6i et produit :
#   (A) frequences de familles par contraste, en distribution sur les replicats,
#       avec les valeurs non equilibrees en regard
#   (B) le VERDICT : dans quelle proportion des replicats la bande inter-periode
#       reste-t-elle entierement disjointe de la bande intra-periode
#   (C) la restriction aux cellules communes (denominateur constant)
#   (D) un objet consensus au schema standard, pour que 7/8/9 tournent dessus
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(tibble)
})

if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

OUT <- sprintf("Robustness_%s%s_BAL", SPATIAL_LEVEL, PREY_FAMILY)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

source("R_helpers/Load_all_runs.R")

tidy_runs <- function(r) {
  r$contrast <- contrast_from_periods(r$period_1, r$period_2)
  r$contrast_type <- ifelse(as.character(r$contrast) %in% INTRA_CONTRASTS,
                            "within", "between")
  r$family <- as.character(family_of(r$diagnostic))
  unit <- intersect(c("Area", "str"), names(r))[1]
  r$cell <- paste(r$species, r$size_class, r[[unit]], sep = "|")
  attr(r, "unit_col") <- unit
  r
}

load(sprintf("data/Sensitivity/all_runs_%s%s_BAL.rda", SPATIAL_LEVEL, PREY_FAMILY))
bal <- tidy_runs(all_runs_balanced)
UNIT <- attr(bal, "unit_col")
CFG  <- attr(all_runs_balanced, "config")
bal_cls <- bal %>% filter(testable, !is.na(family))

n_rep <- n_distinct(bal$rep_id)
message(sprintf("%d replicats | %d lignes | %d cellules",
                n_rep, nrow(bal), n_distinct(bal_cls$cell)))

# =============================================================================
# (A) FREQUENCES PAR REPLICAT
# =============================================================================
# Meme ordre d'agregation que le manuscrit : par unite spatiale, puis moyenne
# sur les unites, puis moyenne sur les resolutions.

pct_by_rep <- bal_cls %>%
  count(rep_id, contrast, contrast_type, mode, x_threshold,
        !!rlang::sym(UNIT), family, name = "n") %>%
  group_by(rep_id, contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT)) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  group_by(rep_id, contrast, contrast_type, mode, x_threshold, family) %>%
  summarise(pct = mean(pct), .groups = "drop") %>%
  group_by(rep_id, contrast, contrast_type, mode, family) %>%
  summarise(pct = mean(pct), .groups = "drop")

fam_summary <- pct_by_rep %>%
  group_by(contrast, contrast_type, mode, family) %>%
  summarise(mean_pct = mean(pct), sd_pct = sd(pct),
            q025 = quantile(pct, .025), q500 = median(pct),
            q975 = quantile(pct, .975), .groups = "drop") %>%
  mutate(family = factor(family, levels = FAMILY_LEVELS)) %>%
  arrange(mode, contrast, family)

write_csv(fam_summary, file.path(OUT, "balanced_family_frequencies.csv"))
write_csv(pct_by_rep,  file.path(OUT, "balanced_family_per_replicate.csv"))

inconclusive <- bal %>%
  filter(testable) %>%
  group_by(rep_id, contrast, mode) %>%
  summarise(pct_inconclusive = 100 * mean(is.na(family)), .groups = "drop")
write_csv(inconclusive, file.path(OUT, "balanced_inconclusive_rate.csv"))

# =============================================================================
# (B) VERDICT
# =============================================================================
# Pour chaque replicat : la bande inter-periode (min-max sur les contrastes
# inter) est-elle entierement disjointe de la bande intra-periode ?

verdict <- pct_by_rep %>%
  filter(family %in% c("Stability", "Substitution")) %>%
  group_by(rep_id, mode, family, contrast_type) %>%
  summarise(lo = min(pct), hi = max(pct), .groups = "drop") %>%
  pivot_wider(names_from = contrast_type, values_from = c(lo, hi)) %>%
  mutate(separated = if_else(family == "Stability",
                             hi_between < lo_within, lo_between > hi_within),
         gap_pp    = if_else(family == "Stability",
                             lo_within - hi_between, lo_between - hi_within)) %>%
  group_by(mode, family) %>%
  summarise(n_rep = n(),
            pct_replicates_separated = 100 * mean(separated, na.rm = TRUE),
            median_gap_pp = median(gap_pp, na.rm = TRUE),
            q025_gap_pp   = quantile(gap_pp, .025, na.rm = TRUE), .groups = "drop")

write_csv(verdict, file.path(OUT, "balanced_separation_verdict.csv"))
cat("\n=================== VERDICT ===================\n")
print(as.data.frame(verdict))
cat("===============================================\n\n")

# =============================================================================
# Comparaison avec la serie non equilibree
# =============================================================================

ref_ok <- tryCatch({
  ref <- load_all_runs(SPATIAL_LEVEL, PREY_FAMILY, verbose = FALSE) %>%
    filter(testable, !is.na(family), x_threshold %in% CFG$X_GRID)
  fam_ref <- ref %>%
    count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT), family, name = "n") %>%
    group_by(contrast, contrast_type, mode, x_threshold, !!rlang::sym(UNIT)) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    group_by(contrast, contrast_type, mode, x_threshold, family) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    group_by(contrast, contrast_type, mode, family) %>%
    summarise(pct_unbalanced = mean(pct), .groups = "drop")

  comparison <- fam_summary %>%
    left_join(fam_ref, by = c("contrast", "contrast_type", "mode", "family")) %>%
    mutate(shift_pp = mean_pct - pct_unbalanced)
  write_csv(comparison, file.path(OUT, "balanced_vs_unbalanced.csv"))
  fam_ref
}, error = function(e) { message("Serie non equilibree indisponible : ",
                                 conditionMessage(e)); NULL })

p <- ggplot(pct_by_rep %>% mutate(family = factor(family, levels = FAMILY_LEVELS)),
            aes(contrast, pct, fill = contrast_type)) +
  geom_violin(scale = "width", alpha = .55, colour = NA) +
  geom_boxplot(width = .14, outlier.size = .3, alpha = .9) +
  facet_grid(mode ~ family, scales = "free_y") +
  labs(x = NULL, y = "% des cellules classifiables", fill = "Contraste",
       title = sprintf("Effort equilibre en traits - %s%s", SPATIAL_LEVEL, PREY_FAMILY),
       subtitle = sprintf("%d replicats%s", n_rep,
                          if (!is.null(ref_ok)) " ; croix = valeur non equilibree" else "")) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (!is.null(ref_ok)) {
  p <- p + geom_point(data = ref_ok %>% mutate(family = factor(family, levels = FAMILY_LEVELS)),
                      aes(contrast, pct_unbalanced), inherit.aes = FALSE,
                      shape = 4, size = 2.6, stroke = 1)
}
ggsave(file.path(OUT, "balanced_family_distributions.png"), p,
       width = 11, height = 6, dpi = 300)

# =============================================================================
# (C) CELLULES COMMUNES
# =============================================================================
# Cellules classifiables dans TOUS les contrastes et TOUS les replicats :
# denominateur constant, ce qui retire la seconde objection sur les pourcentages.

n_sc <- n_distinct(bal_cls$contrast)
common <- bal_cls %>%
  distinct(cell, x_threshold, mode, contrast, rep_id) %>%
  count(cell, x_threshold, mode, name = "n_obs") %>%
  filter(n_obs == n_rep * n_sc)

fam_common <- bal_cls %>%
  semi_join(common, by = c("cell", "x_threshold", "mode")) %>%
  count(rep_id, contrast, contrast_type, mode, x_threshold, family, name = "n") %>%
  group_by(rep_id, contrast, contrast_type, mode, x_threshold) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  group_by(contrast, contrast_type, mode, family) %>%
  summarise(mean_pct = mean(pct), q025 = quantile(pct, .025),
            q975 = quantile(pct, .975), .groups = "drop")
write_csv(fam_common, file.path(OUT, "balanced_common_cells.csv"))
message(sprintf("Cellules communes : %d", n_distinct(common$cell)))

# =============================================================================
# (D) OBJET CONSENSUS POUR 7 / 8 / 9
# =============================================================================
# Une ligne par cellule x resolution x devise x contraste, au schema attendu.
# Conventions a declarer en legende : diagnostic et famille = valeur modale sur
# les replicats ; p-values et tailles d'effet = mediane ; rep_support = part des
# replicats soutenant la modale ; testable = classifiable dans >= 50 % d'entre eux.

modal <- function(x) { t <- table(x); names(t)[which.max(t)] }
UNITQ <- rlang::sym(UNIT)

consensus <- bal %>%
  filter(testable) %>%
  group_by(species, size_class, !!UNITQ, x_threshold, mode, contrast) %>%
  summarise(
    prey_family = dplyr::first(prey_family),
    period_1 = dplyr::first(period_1), period_2 = dplyr::first(period_2),
    spatial_level = dplyr::first(spatial_level),
    n_rep_cell = dplyr::n(),
    testable = dplyr::n() >= 0.5 * n_rep,
    diagnostic = modal(diagnostic),
    rep_support = mean(diagnostic == modal(diagnostic)),
    across(any_of(c("p_comp", "p_H", "p_Bs", "p_disp", "BC", "R2_comp",
                    "delta_H", "delta_Bs", "t_H", "t_Bs", "se_Bs",
                    "n_sto_P1", "n_sto_P2", "n_set_P1", "n_set_P2")),
           ~ median(.x, na.rm = TRUE)),
    reliable = mean(reliable) >= 0.5,
    .groups = "drop")

res <- list(results = as.data.frame(consensus),
            inventory = attr(all_runs_balanced, "design"),
            n_comp_failures = 0L)
assign("res_balanced_consensus", res)
save(res_balanced_consensus,
     file = sprintf("data/Sensitivity/all_runs_%s%s_BALCONS.rda",
                    SPATIAL_LEVEL, PREY_FAMILY))

message("6j termine -> ", OUT,
        " | consensus : data/Sensitivity/all_runs_", SPATIAL_LEVEL,
        PREY_FAMILY, "_BALCONS.rda")
