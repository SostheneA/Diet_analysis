# =============================================================================
# 6i_Balanced_Resampling.R
# -----------------------------------------------------------------------------
# Analyse de sensibilite a effort d'echantillonnage equilibre.
#
# LE PROBLEME
#   Le contraste decennal (2004-2006 vs 2018-2019) met en regard trois annees
#   de campagne contre deux. Son plus petit groupe est plus grand que celui de
#   tout contraste intra-periode, donc les trois tests y ont plus de puissance.
#   La separation observee entre familles intra- et inter-periodes est en partie
#   confondue avec la replication.
#
# LA CORRECTION
#   Relancer le MEME moteur sur des donnees ou chaque contraste est plafonne,
#   cellule par cellule, au meme nombre de traits par periode. Rien dans 6a
#   n'est modifie : on filtre les traits avant make_dat_classed().
#
# SORTIE
#   Sensitivity_<level><family>_BAL/rep_XXXX.rds      un fichier par replicat
#   data/Sensitivity/all_runs_<level><family>_BAL.rda un data.frame + rep_id
#   resample_design_<level><family>.csv               la cible par cellule
#
# Numerote 6i pour ne pas heurter 6g_Compare_families.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble)
})

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
if (!exists("PREY_FAMILY"))   PREY_FAMILY   <- "_1"
if (!exists("SPATIAL_LEVEL")) SPATIAL_LEVEL <- "all_gulf"

N_MIN <- 5; N_STRICT <- 25; MIN_SETS <- 3; R_PERM <- 999; ALPHA <- 0.05

B_REP <- 10
SEED  <- 20260911

# Le balayage complet a 100 seuils x B replicats est hors budget. Grille grossiere
# couvrant la meme plage, a declarer telle quelle dans le manuscrit.
X_GRID <- seq(50, 1000, by = 100)

# Equilibrer les CINQ contrastes, pas seulement le decennal : comparer une valeur
# equilibree a une bande intra-periode non equilibree reintroduit l'asymetrie
# que cette analyse existe pour supprimer.
BALANCE_ALL <- TRUE

# Regle de la cible par cellule :
#   "min_within" le minimum, sur les trois contrastes intra-periode et les deux
#                periodes, du nombre de traits observe dans cette cellule
#   "min_all"    idem sur les cinq contrastes
TARGET_RULE <- "min_within"

# Coeurs. detectCores() renvoie les coeurs LOGIQUES sous Windows ; le master
# reste inactif pendant l'attente, donc on peut les prendre tous. Baisser si la
# memoire sature : chaque worker detient une copie des donnees.
if (!exists("N_CORES")) N_CORES <- max(1L, parallel::detectCores())

# -----------------------------------------------------------------------------
# MOTEUR
# -----------------------------------------------------------------------------
source("6a_Engine_Trophic.R")   # definit SCENARIOS, make_dat_classed, run_pipeline,
# prepare_spatial_strata, SPATIAL_OUT, SPATIAL_SOURCE

OUT_DIR <- sprintf("Sensitivity_%s%s_BAL", SPATIAL_LEVEL, PREY_FAMILY)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create("data/Sensitivity", showWarnings = FALSE, recursive = TRUE)

WITHIN_SC <- c("P1", "P2", "P3")     # codes moteur ; P1a, P2, P1b au manuscrit
RUN_SC    <- if (BALANCE_ALL) names(SCENARIOS) else c("P4", "PT")

# =============================================================================
# 1. DONNEES ET INVENTAIRE DES TRAITS
# =============================================================================

load("data/dat_classed.rda")
diet_clean <- as.data.frame(dat_classed)
if (!is.na(SPATIAL_SOURCE)) {
  if (!SPATIAL_SOURCE %in% names(diet_clean))
    stop("Colonne '", SPATIAL_SOURCE, "' absente de dat_classed.")
  diet_clean[[SPATIAL_SOURCE]] <- as.factor(diet_clean[[SPATIAL_SOURCE]])
}

# set_uid reproduit .add_set_uid() du moteur
diet_clean$set_uid <- paste(diet_clean$year, diet_clean[["vessel.code"]],
                            diet_clean$set, sep = "_")

# Unite spatiale de la cellule, telle que le moteur l'ecrira
unit_vec <- if (isTRUE(SPATIAL_POOL)) rep(SPATIAL_LABEL, nrow(diet_clean)) else
  as.character(diet_clean[[SPATIAL_SOURCE]])
diet_clean$.unit <- unit_vec
diet_clean <- diet_clean[!is.na(diet_clean$.unit), , drop = FALSE]

# Le moteur filtre sur une longueur valide : l'inventaire doit filtrer pareil,
# sinon on plafonne sur des traits que le moteur n'utilisera pas.
inv_src <- diet_clean %>%
  filter(!is.na(somatic_length_cm), somatic_length_cm > 0, !is.na(size_class)) %>%
  distinct(species = predator_species_common_name,
           size = as.character(size_class), unit = .unit,
           set_uid, year)

n_sets_for <- function(sc) {
  pm <- SCENARIOS[[sc]]
  bind_rows(
    inv_src %>% filter(year %in% pm[[1]]) %>%
      count(species, size, unit, name = "n_sets") %>% mutate(period = 1L),
    inv_src %>% filter(year %in% pm[[2]]) %>%
      count(species, size, unit, name = "n_sets") %>% mutate(period = 2L)
  ) %>% mutate(sc = sc)
}

availability <- map_dfr(names(SCENARIOS), n_sets_for)

targets <- availability %>%
  filter(sc %in% if (TARGET_RULE == "min_within") WITHIN_SC else names(SCENARIOS)) %>%
  group_by(species, size, unit) %>%
  summarise(n_target = min(n_sets), .groups = "drop")

cap <- availability %>%
  filter(sc %in% RUN_SC) %>%
  group_by(species, size, unit) %>%
  summarise(n_min_available = min(n_sets), n_sc = n_distinct(sc), .groups = "drop")

design <- targets %>%
  left_join(cap, by = c("species", "size", "unit")) %>%
  mutate(n_min_available = coalesce(n_min_available, 0L),
         complete_sc     = coalesce(n_sc, 0L) == length(RUN_SC),
         n_target_final  = pmin(n_target, n_min_available),
         keep            = complete_sc & n_target_final >= MIN_SETS,
         drop_reason     = case_when(
           !complete_sc                   ~ "cellule absente d'au moins un contraste",
           n_target_final < MIN_SETS      ~ sprintf("cible %d < %d traits", n_target_final, MIN_SETS),
           TRUE                           ~ NA_character_))

write_csv(design, sprintf("resample_design_%s%s.csv", SPATIAL_LEVEL, PREY_FAMILY))

cat(sprintf("\nCellules : %d au total | %d conservees | %d ecartees\n",
            nrow(design), sum(design$keep), sum(!design$keep)))
print(design %>% filter(!keep) %>% count(drop_reason))
n_keep <- sum(design$keep)
if (n_keep == 0L)
  stop("Aucune cellule ne survit a l'equilibrage a ce niveau spatial.\n",
       "Les cibles tombent toutes sous ", MIN_SETS, " traits : les contrastes ",
       "intra-periode n'ont pas assez de traits par cellule ici.\n",
       "Voir resample_design_", SPATIAL_LEVEL, PREY_FAMILY, ".csv.")

cat(sprintf("Cible en traits : mediane %d (min %d, max %d)\n",
            median(design$n_target_final[design$keep]),
            min(design$n_target_final[design$keep]),
            max(design$n_target_final[design$keep])))

pct_drop <- 100 * (1 - n_keep / nrow(design))
cat(sprintf("Part des cellules ecartees : %.1f %%\n", pct_drop))
if (pct_drop > 40)
  warning(sprintf("%.0f %% des cellules sont ecartees. Le resultat porte sur ",
                  pct_drop),
          "un sous-ensemble mieux echantillonne et doit etre declare comme tel.",
          call. = FALSE)

design_keep <- design %>% filter(keep) %>% select(species, size, unit, n_target_final)

# =============================================================================
# 2. UN REPLICAT
# =============================================================================
# Le tirage se fait par paire (cellule, trait) : le moteur traite chaque cellule
# independamment, donc rien n'impose de garder le meme jeu de traits pour tous
# les predateurs -- et un tirage global detruirait les cellules rares.

# Ne garder que les resolutions de la grille : run_pipeline() balaie TOUTES les
# colonnes prey_category_<x><family>, il n'a pas d'argument de seuil.
prey_cols_all <- grep(paste0("^prey_category_\\d+", PREY_FAMILY, "$"),
                      names(diet_clean), value = TRUE)
keep_prey <- paste0("prey_category_", X_GRID, PREY_FAMILY)
keep_prey <- intersect(keep_prey, prey_cols_all)
if (!length(keep_prey))
  stop("Aucune colonne de la grille X_GRID dans dat_classed. Presentes : ",
       paste(head(prey_cols_all, 20), collapse = ", "))
drop_prey <- setdiff(prey_cols_all, keep_prey)
cat(sprintf("Resolutions retenues : %s\n", paste(X_GRID[paste0("prey_category_", X_GRID, PREY_FAMILY) %in% keep_prey], collapse = ", ")))

diet_grid <- diet_clean[, setdiff(names(diet_clean), drop_prey), drop = FALSE]

draw_sets <- function(sc) {
  pm <- SCENARIOS[[sc]]
  inv_src %>%
    mutate(period = case_when(year %in% pm[[1]] ~ 1L,
                              year %in% pm[[2]] ~ 2L,
                              TRUE ~ NA_integer_)) %>%
    filter(!is.na(period)) %>%
    inner_join(design_keep, by = c("species", "size", "unit")) %>%
    group_by(species, size, unit, period) %>%
    # rang aleatoire par groupe puis coupe : equivaut a un tirage sans remise,
    # sans passer par slice_sample(n = <variable>), refuse par dplyr, ni par
    # first(), masque par data.table.
    mutate(.rank = sample.int(dplyr::n())) %>%
    filter(.rank <= n_target_final) %>%
    ungroup() %>%
    select(species, size, unit, set_uid)
}

# --- Grille de taches ---------------------------------------------------------
# On parallelise au niveau (replicat x contraste x devise) et non au niveau du
# replicat. Avec B_REP = 100 cela fait 1000 taches courtes au lieu de 100 taches
# longues, donc un equilibrage de charge bien meilleur : aucun coeur ne reste
# inactif en fin de run a attendre le dernier replicat.

TASKS <- expand.grid(rep_id = seq_len(B_REP), sc = RUN_SC,
                     mode = c("biomass", "occurrence"),
                     stringsAsFactors = FALSE)
TASKS <- TASKS[order(TASKS$rep_id, TASKS$sc, TASKS$mode), ]
rownames(TASKS) <- NULL

run_one_task <- function(k) {
  b  <- TASKS$rep_id[k]
  sc <- TASKS$sc[k]
  md <- TASKS$mode[k]
  pm <- SCENARIOS[[sc]]

  # Re-seeder AVANT chaque tirage. run_pipeline() refixe le grain a chaque appel
  # (set.seed(1234), puis seed + 1000*x + i par cellule), donc l'etat du
  # generateur en sortie est identique d'un replicat a l'autre. Le germe depend
  # du replicat ET du contraste, jamais de la devise, pour que les deux devises
  # d'un meme replicat voient exactement les memes traits.
  set.seed(SEED + 1000L * b + match(sc, RUN_SC))
  keep <- draw_sets(sc)

  d_b <- diet_grid %>%
    inner_join(keep,
               by = c("predator_species_common_name" = "species",
                      "size_class" = "size", ".unit" = "unit",
                      "set_uid" = "set_uid"))
  if (!nrow(d_b)) return(NULL)

  dat_b <- make_dat_classed(d_b, period_map = pm)

  r <- tryCatch(
    run_pipeline(dat_b, mode = md,
                 period_1 = names(pm)[1], period_2 = names(pm)[2])$results,
    error = function(e) {
      warning(sprintf("rep %d | %s | %s : %s", b, sc, md, conditionMessage(e)),
              call. = FALSE); NULL })
  if (is.null(r) || !nrow(r)) return(NULL)

  r$rep_id <- b; r$sc_code <- sc
  saveRDS(r, file.path(OUT_DIR, sprintf("task_%04d_%s_%s.rds", b, sc, md)))
  r
}

# =============================================================================
# 3. EXECUTION
# =============================================================================

cat(sprintf("\n%d replicats x %d contrastes x 2 devises x %d resolutions\n",
            B_REP, length(RUN_SC), length(keep_prey)))
cat(sprintf("%d taches reparties sur %d coeurs\n", nrow(TASKS), N_CORES))
cat("Chronometre d'abord avec B_REP <- 10.\n\n")

t0 <- Sys.time()

run_parallel <- function(idx, fun, n_cores) {
  if (n_cores <= 1L) return(lapply(idx, fun))

  wd <- getwd()
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Chaque worker charge le moteur une seule fois, avec le bon niveau spatial.
  parallel::clusterExport(cl, c("wd", "SPATIAL_LEVEL", "PREY_FAMILY",
                                "N_MIN", "N_STRICT", "MIN_SETS",
                                "R_PERM", "ALPHA"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    setwd(wd)
    suppressPackageStartupMessages({
      library(dplyr); library(tidyr); library(purrr); library(tibble)
    })
    sink(nullfile())            # le moteur ecrit sa progression ; on la coupe
    source("6a_Engine_Trophic.R")
    sink()
    TRUE
  })

  parallel::clusterExport(cl, c("diet_grid", "inv_src", "design_keep",
                                "RUN_SC", "SEED", "OUT_DIR", "draw_sets",
                                "TASKS", "SCENARIOS"),
                          envir = environment())

  # Taches longues en premier : meilleur remplissage en fin de run.
  parallel::parLapplyLB(cl, idx, fun)
}

res_list <- run_parallel(seq_len(nrow(TASKS)), run_one_task, N_CORES)
reps <- Filter(Negate(is.null), res_list)

all_runs_balanced <- bind_rows(reps)
cat(sprintf("\nTermine en %.1f min | %d lignes\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            nrow(all_runs_balanced)))

attr(all_runs_balanced, "design") <- design
attr(all_runs_balanced, "config") <- list(
  B_REP = B_REP, SEED = SEED, X_GRID = X_GRID, TARGET_RULE = TARGET_RULE,
  BALANCE_ALL = BALANCE_ALL, MIN_SETS = MIN_SETS, R_PERM = R_PERM,
  SPATIAL_LEVEL = SPATIAL_LEVEL, PREY_FAMILY = PREY_FAMILY)

# Controle : chaque contraste doit varier d'un replicat a l'autre. Une variance
# nulle signale que le tirage n'a pas ete re-seede.
if (B_REP > 1L) {
  chk <- all_runs_balanced %>%
    filter(testable) %>%
    count(sc_code, mode, rep_id, diagnostic) %>%
    group_by(sc_code, mode, diagnostic) %>%
    summarise(sd_n = sd(n), .groups = "drop") %>%
    group_by(sc_code, mode) %>%
    summarise(varies = any(sd_n > 0, na.rm = TRUE), .groups = "drop")
  cat("\nVariation entre replicats, par contraste :\n")
  print(as.data.frame(chk))
  if (any(!chk$varies))
    warning("Certains contrastes ne varient pas entre replicats : le tirage ",
            "n'est pas aleatoire pour eux.", call. = FALSE)
}

save(all_runs_balanced,
     file = sprintf("data/Sensitivity/all_runs_%s%s_BAL.rda",
                    SPATIAL_LEVEL, PREY_FAMILY))
cat("Sortie : data/Sensitivity/all_runs_", SPATIAL_LEVEL, PREY_FAMILY, "_BAL.rda\n", sep = "")
