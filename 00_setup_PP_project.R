# =============================================================================
# 00_setup_PP_project.R — Cree le projet PP separe, a lancer DEPUIS le projet
# principal (working directory = projet principal).
# =============================================================================
# Philosophie : DONNEES ET RESULTATS SEPARES, MOTEUR PARTAGE.
#   - Le projet PP a son propre data/ et ses propres Sensitivity_* : aucun
#     risque d'ecraser ou de melanger avec le pooled.
#   - 6a n'est PAS copie : les run scripts du projet PP le sourcent depuis le
#     projet principal -> une seule version du moteur, toute correction
#     profite aux deux projets.
#   - Seuls 3PP, 4PP et 5 vivent dans le projet PP (5 est une copie conforme).
# =============================================================================

MAIN_PROJECT <- normalizePath(here::here())
PP_PROJECT   <- file.path(dirname(MAIN_PROJECT), paste0(basename(MAIN_PROJECT), "_PP"))

dir.create(PP_PROJECT, showWarnings = FALSE)
for (d in c("data", "output/Figures")) {
  dir.create(file.path(PP_PROJECT, d), recursive = TRUE, showWarnings = FALSE)
}

# 1) Donnees d'entree : copie de tout data/ SAUF les sorties de runs ----------
#    (les .rda d'entree + shapefiles Spatial_data ; on exclut Sensitivity et
#     les gros produits regenerables par 4PP/5)
src_data <- list.files("data", recursive = TRUE, full.names = TRUE)
# unique_samples.rda est copie volontairement : il ne depend que du cote
# predateur et son recalcul (appels weight par ligne) est long.
excl     <- grepl("Sensitivity|all_runs|diet_clean|dat_classed|Database\\.rda",
                  src_data)
for (f in src_data[!excl]) {
  dest <- file.path(PP_PROJECT, f)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  file.copy(f, dest, overwrite = FALSE)
}

# 2) Helpers : copie de R_helpers/ (Config_Mappings.R inclus) -----------------
dir.create(file.path(PP_PROJECT, "R_helpers"), showWarnings = FALSE)
file.copy(list.files("R_helpers", full.names = TRUE),
          file.path(PP_PROJECT, "R_helpers"), overwrite = FALSE)

# 3) Scripts propres au projet PP ---------------------------------------------
#    3PP + 4PP doivent etre places dans le projet principal avant de lancer ce
#    setup ; 5 est copie tel quel (il selectionne par pattern, rien a changer).
for (f in c("3PP_taxonomic_groups.R", "4PP_Database_with_others_factors.R",
            "5_From_cleandiet_to_dat_classed.R")) {
  if (file.exists(f)) file.copy(f, file.path(PP_PROJECT, f), overwrite = TRUE)
  else warning("Introuvable dans le projet principal : ", f)
}

# 4) Run scripts PP : copies de 6b/6c/6d qui (a) fixent PREY_FAMILY = "_PP",
#    (b) sourcent le moteur 6a DU PROJET PRINCIPAL (pas de copie du moteur).
for (f in c("6b_Run_L1_all_gulf.R", "6c_Run_L2_ecoregion.R", "6d_Run_L3_stratum.R")) {
  txt <- readLines(f, warn = FALSE)

  # remplace toute affectation PREY_FAMILY existante par "_PP"
  txt <- sub('^(\\s*PREY_FAMILY\\s*(<-|=)\\s*)"_[12]"', '\\1"_PP"', txt)

  # source du moteur : chemin absolu vers le projet principal
  txt <- sub('source\\("6a_Engine_Trophic.R"\\)',
             paste0('MAIN_PROJECT <- "', MAIN_PROJECT, '"\n',
                    'PREY_FAMILY <- "_PP"   # projet PP\n',
                    'source(file.path(MAIN_PROJECT, "6a_Engine_Trophic.R"))'),
             txt)

  writeLines(txt, file.path(PP_PROJECT, f))
}

# 5) Fichier projet RStudio ----------------------------------------------------
writeLines(c("Version: 1.0", "", "RestoreWorkspace: No", "SaveWorkspace: No"),
           file.path(PP_PROJECT, paste0(basename(PP_PROJECT), ".Rproj")))

cat("\nProjet PP cree :", PP_PROJECT,
    "\nOrdre d'execution dans le projet PP : 3PP -> 4PP -> 5 -> 6b/6c/6d",
    "\n(6a reste dans le projet principal et y est source ; PREY_FAMILY = \"_PP\")\n")
