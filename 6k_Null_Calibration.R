# =============================================================================
# 6k_Null_Calibration.R
# -----------------------------------------------------------------------------
# Calibration empirique des trois tests sous une hypothese nulle construite a
# partir des donnees reelles.
#
# PRINCIPE
#   Dans chaque cellule, les etiquettes de periode sont permutees ENTRE LES
#   TRAITS, en conservant le nombre de traits par periode et toute la structure
#   interne (composition des traits, nombre d'estomacs par trait, richesse).
#   Sous cette permutation il n'y a, par construction, aucune difference entre
#   periodes. Le taux de rejet a alpha = 0.05 est donc le taux de faux positifs
#   reel de chaque test, a l'effectif de CETTE cellule.
#
# CE QUE CA REGLE
#   L'objection la plus lourde du papier : "Substitution est definie par un test
#   significatif et deux non significatifs, c'est donc exactement ce qu'on
#   attend si PERMANOVA est plus sensible que les deux tests univaries."
#   Si les trois tests rejettent autour de 5 %, la typologie est calibree et
#   l'objection tombe. Si PERMANOVA rejette a 12 % et H' a 2 %, il faut le dire
#   et calibrer les seuils par test.
#
#   La distribution NULLE des familles est aussi une meilleure reference interne
#   que les contrastes intra-periode, qui contiennent du vrai signal ecologique.
#
# AUCUNE modification de 6a : la permutation agit sur la colonne `period` du
# jeu construit par make_dat_classed(), avant l'appel a run_pipeline().
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(ggplot2); library(tibble)
})

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

N_MIN <- 5; N_STRICT <- 25; MIN_SETS <- 3; R_PERM <- 999; ALPHA <- 0.05

B_PERM <- 200            # permutations des etiquettes de periode
SEED   <- 20260911
X_GRID <- seq(100, 1000, by = 100)   # grille de resolutions (10 valeurs)
SC     <- "PT"           # contraste dont on calibre les tests (decennal)

OUT <- sprintf("NullCalibration_%s%s", SPATIAL_LEVEL, PREY_FAMILY)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

source("6a_Engine_Trophic.R")
for (.p in c("R_helpers/Config_Mappings.R", "Config_Mappings.R"))
  if (file.exists(.p)) { source(.p); break }

# =============================================================================
# 1. DONNEES
# =============================================================================

load("data/dat_classed.rda")
diet_clean <- as.data.frame(dat_classed)
if (!is.na(SPATIAL_SOURCE)) {
  if (!SPATIAL_SOURCE %in% names(diet_clean))
    stop("Colonne '", SPATIAL_SOURCE, "' absente de dat_classed.")
  diet_clean[[SPATIAL_SOURCE]] <- as.factor(diet_clean[[SPATIAL_SOURCE]])
}

# Ne garder que les resolutions de la grille : run_pipeline() balaie toutes les
# colonnes prey_category_<x><family>, il n'a pas d'argument de seuil.
prey_all  <- grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"),
                  names(diet_clean), value = TRUE)
prey_keep <- intersect(paste0("prey_category_", X_GRID, PREY_FAMILY), prey_all)
if (!length(prey_keep))
  stop("Aucune colonne de X_GRID presente. Disponibles : ",
       paste(head(prey_all, 15), collapse = ", "))
diet_grid <- diet_clean[, setdiff(names(diet_clean), setdiff(prey_all, prey_keep)),
                        drop = FALSE]
cat(sprintf("Resolutions : %d colonnes retenues sur %d\n",
            length(prey_keep), length(prey_all)))

pm  <- SCENARIOS[[SC]]
dat <- make_dat_classed(diet_grid, period_map = pm)
dat <- prepare_spatial_strata(dat)
dat$set_uid <- paste(dat$year, dat[["vessel.code"]], dat$set, sep = "_")

PERIOD_LEVELS <- levels(dat$period)

# =============================================================================
# 2. PERMUTATION DES ETIQUETTES DE PERIODE
# =============================================================================
# La permutation est faite A L'INTERIEUR DE CHAQUE CELLULE, en conservant n1 et
# n2. Une permutation globale melangerait les deficits d'effort entre cellules
# et fabriquerait un nul plus permissif que la realite.
#
# Un meme trait peut recevoir des etiquettes differentes selon la cellule : le
# moteur traite chaque cellule independamment, donc c'est legitime et c'est ce
# qui preserve exactement l'effectif de chaque cellule.

key <- dat %>%
  distinct(sp = predator_species_common_name,
           sz = as.character(size_class),
           un = as.character(.data[[SPATIAL_OUT]]),
           set_uid, period)

permute_labels <- function() {
  key %>%
    group_by(sp, sz, un) %>%
    mutate(period_perm = sample(period)) %>%
    ungroup()
}

apply_labels <- function(kp) {
  dat %>%
    select(-period) %>%
    inner_join(
      kp %>% select(sp, sz, un, set_uid, period = period_perm),
      by = c("predator_species_common_name" = "sp",
             "size_class" = "sz", "set_uid" = "set_uid")) %>%
    filter(as.character(.data[[SPATIAL_OUT]]) == un) %>%
    select(-un) %>%
    mutate(period = factor(as.character(period), levels = PERIOD_LEVELS))
}

run_one_perm <- function(b) {
  set.seed(SEED + 1e5 + b)
  d_p <- apply_labels(permute_labels())
  out <- list()
  for (md in c("biomass", "occurrence")) {
    r <- tryCatch(
      run_pipeline(d_p, mode = md,
                   period_1 = names(pm)[1], period_2 = names(pm)[2])$results,
      error = function(e) { warning("perm ", b, " / ", md, " : ",
                                    conditionMessage(e), call. = FALSE); NULL })
    if (is.null(r) || !nrow(r)) next
    r$perm_id <- b
    out[[length(out) + 1L]] <- r
  }
  bind_rows(out)
}

# =============================================================================
# 3. EXECUTION
# =============================================================================

cat(sprintf("\nNul : %d permutations x 2 devises x %d resolutions\n",
            B_PERM, length(prey_keep)))
cat("Chronometre d'abord avec B_PERM <- 2.\n\n")

t0 <- Sys.time()

N_CORES <- max(1L, parallel::detectCores() - 1L)

run_parallel <- function(idx, fun, n_cores) {
  if (n_cores <= 1L) return(lapply(idx, fun))
  wd <- getwd()
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("wd", "SPATIAL_LEVEL", "PREY_FAMILY",
                                "N_MIN", "N_STRICT", "MIN_SETS", "R_PERM", "ALPHA"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    setwd(wd)
    suppressPackageStartupMessages({
      library(dplyr); library(tidyr); library(purrr); library(tibble)
    })
    source("6a_Engine_Trophic.R")
    TRUE
  })
  parallel::clusterExport(cl, c("dat", "key", "pm", "SEED", "PERIOD_LEVELS",
                                "permute_labels", "apply_labels"),
                          envir = environment())
  parallel::parLapplyLB(cl, idx, fun)
}

null_runs <- bind_rows(run_parallel(seq_len(B_PERM), run_one_perm, N_CORES))
cat(sprintf("\nTermine en %.1f min | %d lignes\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(null_runs)))

save(null_runs, file = sprintf("data/null_runs_%s%s.rda", SPATIAL_LEVEL, PREY_FAMILY))

# =============================================================================
# 4. RESULTATS
# =============================================================================

nr <- null_runs %>% filter(testable)

# ---- (1) taux de faux positifs par test -------------------------------------
fpr <- nr %>%
  group_by(mode) %>%
  summarise(n = n(),
            fpr_comp = 100 * mean(p_comp < ALPHA, na.rm = TRUE),
            fpr_H    = 100 * mean(p_H    < ALPHA, na.rm = TRUE),
            fpr_Bs   = 100 * mean(p_Bs   < ALPHA, na.rm = TRUE),
            fpr_disp = 100 * mean(p_disp < ALPHA, na.rm = TRUE), .groups = "drop")
write_csv(fpr, file.path(OUT, "null_false_positive_rates.csv"))
cat("\n--- Taux de rejet sous le nul (nominal 5 %) ---\n")
print(as.data.frame(fpr))

# ---- (2) calibration en fonction de l'effectif de la cellule ----------------
fpr_by_n <- nr %>%
  mutate(n_min = pmin(n_set_P1, n_set_P2),
         n_bin = cut(n_min, breaks = c(2, 3, 5, 8, 12, 20, Inf), right = FALSE)) %>%
  group_by(mode, n_bin) %>%
  summarise(n = n(),
            `Composition (PERMANOVA)` = 100 * mean(p_comp < ALPHA, na.rm = TRUE),
            `Diversite (H')`          = 100 * mean(p_H    < ALPHA, na.rm = TRUE),
            `Amplitude (Bs)`          = 100 * mean(p_Bs   < ALPHA, na.rm = TRUE),
            .groups = "drop")
write_csv(fpr_by_n, file.path(OUT, "null_fpr_by_sample_size.csv"))

p_cal <- fpr_by_n %>%
  pivot_longer(c(`Composition (PERMANOVA)`, `Diversite (H')`, `Amplitude (Bs)`),
               names_to = "test", values_to = "fpr") %>%
  ggplot(aes(n_bin, fpr, colour = test, group = test)) +
  geom_hline(yintercept = 100 * ALPHA, linetype = 2, colour = "grey40") +
  geom_line() + geom_point() + facet_wrap(~ mode) +
  labs(x = "Traits par periode (minimum des deux)",
       y = "Taux de rejet sous le nul (%)", colour = NULL,
       title = "Calibration des trois tests sous permutation des etiquettes de periode",
       subtitle = "Pointilles : nominal 5 %. Un ecart entre courbes = sensibilite differentielle") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "null_calibration.png"), p_cal, width = 9, height = 4.5, dpi = 300)

# ---- (3) distribution NULLE des familles ------------------------------------
nf <- nr %>%
  mutate(family = as.character(family_of(diagnostic))) %>%
  filter(!is.na(family)) %>%
  count(perm_id, mode, x_threshold, family) %>%
  group_by(perm_id, mode, x_threshold) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  group_by(perm_id, mode, family) %>%
  summarise(pct = mean(pct), .groups = "drop")

null_fam <- nf %>%
  group_by(mode, family) %>%
  summarise(null_mean = mean(pct), null_q025 = quantile(pct, .025),
            null_q975 = quantile(pct, .975), .groups = "drop")
write_csv(null_fam, file.path(OUT, "null_family_distribution.csv"))
cat("\n--- Distribution nulle des familles ---\n")
print(as.data.frame(null_fam))

cat("\nA rapporter : la frequence OBSERVEE de Substitution (M0) contre sa borne\n",
    "superieure nulle a 97.5 %, et l'ecart entre les trois courbes de calibration.\n", sep = "")
