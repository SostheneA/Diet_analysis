#------------------------------------------------------------------------------#
# Objective: Find a logic to do groups and
#           reduce redundancy in taxonomic hierarchy for the analysis
#
#          - First group are the lowest possible with the available taxonomic
#               hierarchy
#          - Second are a function that may be based on the following rules
#             - keep the lowest if n_stomach_per_prey >= 100 &
#                                   n_predator_per_prey >= 2
#             - otherwise keep if it exists or meaningful/not much redundancy
#               (in this order) : "phylum", "subphylum", "class", "subclass", "order",
#                                 "infraorder", "family", "genus", "species"
#
#           **Either way, prey life stages are kept, but prey_size_cat
#             are discarded.
#             Intestinal parasites are not considered prey.
#             Plantae-Chromista mixed identification is considerate as
#               a dual Kingdom (may be to removed.)
#             General prey groups are "identified_prey_species", "plantae_chromista"
#                     "parasite", "egg_larvae", "empty", "digested"
#
# *****This must be run carefully since it may require manual execution****
#
# Input : /data/prey_taxon_data.rda
#         /R_helpers/taxonomic_rank_order.R"
#         /R_helpers/RemoveAllNasCol.R"
#
# Output : /data/prey_groups.RData
#
# Documentation :
#  With the help of this way of classifications (may be to update)
#  https://publications.gc.ca/collections/collection_2018/mpo-dfo/fs70-5/Fs70-5-2018-003-fra.pdf
#  WoRMS, ITIS
# see https://publications.gc.ca/collections/collection_2020/mpo-dfo/Fs97-6-3383-eng.pdf
#     https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2022.963039/full
#------------------------------------------------------------------------------#

#rm(list = ls())
gc()

# 1) Environment & Data Loading ------------------------------------------------
#------------------------------------------------------------------------------#
project_path <- here::here()

## a) load libraries and helpers -----------------------------------------------
library(data.table)

source(paste0(project_path, "/R_helpers/RemoveAllNasCol.R"))
source(paste0(project_path, "/R_helpers/taxonomic_rank_order.R"))
load(paste0(project_path, "/data/prey_taxon_data.rda"))
prey_groups <-  rbind(species_prey_id, species_prey_noid, fill = T)

setDT(prey_groups)

unique(prey_groups$prey_category)


# 2) prey_category_lowest (the lowest possible classification) -----------------
#------------------------------------------------------------------------------#
prey_groups[, ':='(
  prey_category_lowest = fcase(
    !is.na(verified_name) & prey_category == "identified_prey_species", as.character(verified_name),
    !is.na(prey_category) & prey_category != "identified_prey_species", as.character(prey_category),
    default = as.character(prey_category)
  ),
  tax_level_lowest = fcase(
    !is.na(verified_name) & prey_category == "identified_prey_species", as.character(tax_level),
    !is.na(prey_category)& prey_category != "identified_prey_species", "prey_category",
    default = as.character("prey_category")
  )
)]

setDT(prey_groups)

unique(prey_groups$prey_category)
table(prey_groups$prey_category)
length(unique(prey_groups$prey_category_lowest))
nrow(prey_groups)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#



## b) Identified Prey handling function ----------------------------------------
#------------------------------------------------------------------------------#
# Improved prey categorization function for a data.table
PreyCategory <- function(dt,
                                   target_ranks = c("phylum", "subphylum", "class", "subclass", "order",
                                                    "infraorder", "family", "genus", "species"),
                                   full_taxonomy,
                                   target_rows = NULL,
                                   n_stomach_threshold = 100,
                                   n_predator_threshold = 2,
                                   agg_level_var = "agg_level",
                                   agg_name_var = "agg_name") {

  dt <- copy(as.data.table(dt))

  # A. all_target_idx is the FULL set of rows in scope and never shrinks; it is
  #    what the thresholds are evaluated on.
  #    rem_idx shrinks as rows get assigned to a taxa level (rank) and drives the assignment only.
  all_target_idx <- if (is.null(target_rows)) seq_len(nrow(dt)) else which(target_rows)
  rem_idx <- all_target_idx

  if (length(rem_idx) == 0) return(dt)

  ## Initialize columns to avoid 'set' masking conflicts
  if (!agg_level_var %in% names(dt)) dt[, (agg_level_var) := NA_character_]
  if (!agg_name_var %in% names(dt))  dt[, (agg_name_var) := NA_character_]

  ## Reset targeted rows to NA
  dt[rem_idx, c(agg_level_var, agg_name_var) := NA_character_]

  ## Filter ranks existing in data and sort Bottom-Up (Specific to General)
  ordered_target_ranks <- full_taxonomy[full_taxonomy %in% target_ranks]
  valid_cols <- intersect(ordered_target_ranks, names(dt))
  spec_to_gen <- rev(valid_cols)

  # Bi. Threshold evaluation, on ALL target rows, BEFORE any assignment.
  valid_names_per_rank <- list()
  for (rank in spec_to_gen) {
    rank_vec <- dt[[rank]]
    rows <- all_target_idx[!is.na(rank_vec[all_target_idx]) & rank_vec[all_target_idx] != ""]
    if (length(rows) == 0) next

    agg <- dt[rows, .(
      n_stomach  = uniqueN(unlist(list_stomach_id)),
      n_predator = uniqueN(unlist(list_predator_id))
    ), by = c(rank)]

    valid_names_per_rank[[rank]] <- agg[
      n_stomach >= n_stomach_threshold & n_predator >= n_predator_threshold, get(rank)]
  }

  # Bii. Assignment, still bottom-up, but using the counts computed above.
  for (rank in spec_to_gen) {
    if (length(rem_idx) == 0) break

    valid_names <- valid_names_per_rank[[rank]]
    if (is.null(valid_names) || length(valid_names) == 0) next

    rank_vec  <- dt[[rank]]
    valid_idx <- rem_idx[!is.na(rank_vec[rem_idx]) & rank_vec[rem_idx] != ""]
    update_idx <- valid_idx[rank_vec[valid_idx] %chin% valid_names]

    if (length(update_idx) > 0) {
      dt[update_idx, c(agg_level_var, agg_name_var) := .(rank, rank_vec[update_idx])]
      rem_idx <- setdiff(rem_idx, update_idx)
    }
  }

  # C. Fallback logic for remaining residual rows
  if (length(rem_idx) > 0) {
    char_list <- lapply(dt[rem_idx, .SD, .SDcols = valid_cols], function(x) {
      x <- as.character(x)
      x[x == ""] <- NA_character_
      x
    })

    dt[rem_idx, c(agg_level_var, agg_name_var) := .(
      valid_cols[1],
      do.call(fcoalesce, char_list)
    )]
  }

  return(dt)
}


## c) Identified Prey handling loop n_predator = 2, species reduce to genus ----
#------------------------------------------------------------------------------#

# Defining target ranks (taxa levels) to run
chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")

# Création de la liste des phylums (Animalia non NA), on which the loop work
phyla_list <- unique(prey_groups[kingdom == "Animalia" & !is.na(phylum), phylum])

# Séquence de 10 à 1000 (incrément de 10), on what the loop run
thresholds <- seq(10, 1000, by = 10)

for (x in thresholds) {

  # A. Noms dynamiques des colonnes
  col_agg_level <- paste0("tax_level_", x, "_2")
  col_agg_name  <- paste0("prey_category_", x, "_2")

  # B. Appel de la fonction PreyCategory
  prey_groups <- PreyCategory(
    dt = prey_groups,
    target_ranks = chosen_ranks,
    full_taxonomy = ranknfile,
    target_rows = prey_groups$phylum %in% phyla_list, # Utilise la phyla_list définie plus haut
    n_stomach_threshold = x,
    n_predator_threshold = 2,
    agg_level_var = col_agg_level,
    agg_name_var = col_agg_name
  )

  # C. Harmonisation selon prey_category (fcase)
  # On utilise (col) := pour l'assignation dynamique
  prey_groups[, (col_agg_name) := fcase(
    prey_category == "identified_prey_species", get(col_agg_name),
    default = as.character(prey_category)
  )]

  prey_groups[, (col_agg_level) := fcase(
    prey_category == "identified_prey_species", get(col_agg_level),
    default = "prey_category"
  )]

  # D. Nettoyage Genus -> Species
  # Récupère les genres qui sont devenus des catégories finales
  current_genuses <- unique(prey_groups[get(col_agg_level) == "genus"][[col_agg_name]])

  if (length(current_genuses) > 0) {
    # Correctif
    prey_groups[get(col_agg_level) == "species" & genus %in% current_genuses,
                c(col_agg_name, col_agg_level) := .(genus, "genus")]
  }

  # E. Gestion des Phyla sous le seuil (Invalid Phyla)
  # On recalcule les invalides spécifiquement pour le seuil x actuel
  invalid_phyla <- prey_groups[prey_category == "identified_prey_species",
                               .(total_stomach = uniqueN(unlist(list_stomach_id)),
                                 total_pred    = uniqueN(unlist(list_predator_id))),
                               by = phylum
  ][total_stomach < x | total_pred < 2, phylum]

  if (length(invalid_phyla) > 0) {
    prey_groups[get(col_agg_level) == "phylum" & phylum %in% invalid_phyla,
                (col_agg_name) := "other_phyla"]
  }

  # Message de suivi pour voir l'avancement dans la console
  message(paste("[_2] species reduced to genus - Calcul terminé pour le seuil x =", x))
  gc()
}


#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#

## D) Identified Prey handling loop n_predator = 1, no species reduction to genus ----
#------------------------------------------------------------------------------#
# Notes : Lorsqu'on a écrit la fonction avec un n_predateur = 2, on n'avait en tête
# la représentaion de toutes les proies de tous les prédateurs ensembles et pas l'analyse.
# Pour faire l'analyse, la version par prédateur est préférable (voir 3b_question2_resolution_sensitivity.R)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Suggestion 1 - correctifs de sur-résolution dû à la boucle avec plusieurs prédateur
#
#   1) n_predator_threshold = 1  -> une proie propre à un seul prédateur
#      n'est plus remontée à un rang grossier.
#   2) Suppression du collapse espèce -> genre (l'étape 4 de la boucle _2).
#      La règle des frères (prey_category_new) le produit naturellement
#
# Les colonnes _2 restent INCHANGÉES dans le fichier : les deux familles
#  (n_predator = _1, predator = _2) coexistent dans
# prey_groups pour permettre la comparaison _1 vs _2 en aval.
#------------------------------------------------------------------------------#

for (x in thresholds) {

  col_agg_level <- paste0("tax_level_", x, "_1")
  col_agg_name  <- paste0("prey_category_", x, "_1")

  prey_groups <- PreyCategory(
    dt = prey_groups,
    target_ranks = chosen_ranks,
    full_taxonomy = ranknfile,
    target_rows = prey_groups$phylum %in% phyla_list,
    n_stomach_threshold = x,
    n_predator_threshold = 1,          # correctif
    agg_level_var = col_agg_level,
    agg_name_var = col_agg_name
  )

  prey_groups[, (col_agg_name) := fcase(
    prey_category == "identified_prey_species", get(col_agg_name),
    default = as.character(prey_category)
  )]

  prey_groups[, (col_agg_level) := fcase(
    prey_category == "identified_prey_species", get(col_agg_level),
    default = "prey_category"
  )]

  # correctif : PAS de collapse Genus -> Species pour les colonnes _1.

  # Phyla sous le seuil : n_predator étant retiré, seul le seuil d'estomacs
  # (total_stomach < x) déclasse un phylum vers "other_phyla".
  invalid_phyla <- prey_groups[prey_category == "identified_prey_species",
                               .(total_stomach = uniqueN(unlist(list_stomach_id))),
                               by = phylum
  ][total_stomach < x, phylum]

  if (length(invalid_phyla) > 0) {
    prey_groups[get(col_agg_level) == "phylum" & phylum %in% invalid_phyla,
                (col_agg_name) := "other_phyla"]
  }

  message(paste("[_1] Calcul terminé pour le seuil x =", x))
  gc()
}


#------------------------------------------------------------------------------#
# Suggestion 2 - correctifs de sur-résolution dû à la boucle avec plusieurs prédateur ----
# ** pooled vs per-predator grouping -----------------------
#------------------------------------------------------------------------------#
# The loops above build the taxonomic grouping with ALL PREDATORS POOLED: the
# thresholds are evaluated on uniqueN(list_stomach_id) summed over every
# predator. But the downstream analyses (script 6a onward) are done PER
# PREDATOR. A fine category that clears the threshold when pooled can be nearly
# empty for a single predator -> per-predator over-resolution (see the diagnostic
# in 3c_question2_resolution_sensitivity.R).
#
# We rerun the SAME grouping (the _1 rule: n_predator = 1, no species->genus
# collapse) but SEPARATELY within each predator's own stomachs, then compare the
# number of prey categories obtained per predator: pooled vs separate.in 3c_question2_resolution_sensitivity.R).


# ============================================================================ #
# prey_groups_by_predator : script 3's _1 grouping sweep, rebuilt PER PREDATOR
# ---------------------------------------------------------------------------- #
# Script 3 runs the _1 loop (for x in thresholds) with ALL PREDATORS POOLED.
# Because within ONE predator n_predator is always 1, the corrected _1 rule
# (n_predator_threshold = 1) is exactly what a per-predator grouping needs.
# So we run the SAME threshold loop separately inside each predator's stomachs
# and stack the results. Output has one row per (predator, prey taxon) and, for
# every threshold x, the pair prey_category_<x>_1 / tax_level_<x>_1.
# ============================================================================ #
if (!exists("diet_corr")) {
  load(paste0(project_path, "/data/diet_corr.RData"))
  setDT(diet_corr)
}
# Taxonomy + pooled _1 labels, one row per prey taxon (keyed by common name).
rank_cols <- intersect(
  c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","parvphylum",
    "gigaclass","superclass","class","subclass","infraclass","subterclass","superorder",
    "order","suborder","infraorder","parvorder","section","subsection","superfamily",
    "family","subfamily","tribe","genus","subgenus","species"),
  names(prey_groups))
pcols1      <- paste0("prey_category_", thresholds, "_1")   # created by the _1 loop above
lookup_cols <- c("prey_species_common_name", "prey_category", "verified_name",
                 "tax_level", rank_cols, pcols1)
tax_lookup  <- unique(prey_groups[, ..lookup_cols], by = "prey_species_common_name")

# Main predators = the study set: those with >= 100 stomachs sampled (whole diet,
# regardless of prey), as in script 2. This is 40 predators.
MAIN_PRED_MIN_STOM <- 100L
main_predators <- diet_corr[!is.na(stomach_id),
                            .(n_sto = uniqueN(stomach_id)),
                            by = predator_species_common_name
][n_sto >= MAIN_PRED_MIN_STOM, predator_species_common_name]

# One row per (predator, stomach, prey taxon) for the 40 main predators, with
# taxonomy attached. rec_all already carries EVERY prey category (identified +
# digested/parasite/egg_larvae/etc.), just like the pooled prey_groups.
rec_all <- diet_corr[!is.na(stomach_id) & predator_species_common_name %in% main_predators,
                     .(predator = predator_species_common_name, stomach_id, prey_species_common_name)]
rec_all <- tax_lookup[rec_all, on = "prey_species_common_name", nomatch = 0L]

# Columns that define one unique prey taxon (identity + full taxonomy).
by_cols <- c("prey_species_common_name", "prey_category", "verified_name", "tax_level", rank_cols)

# BUILD: predator loop wrapping script 3's _1 threshold loop ------------- #
prey_groups_by_predator <- rbindlist(lapply(main_predators, function(p) {

  # Collapse this predator's stomach-level records to ONE row per prey taxon.
  # list_stomach_id feeds PreyCategory's stomach-count test; list_predator_id is
  # a single element (p), so n_predator = 1 by construction.
  pg <- rec_all[predator == p,
                .(list_stomach_id  = list(unique(stomach_id)),
                  list_predator_id = list(p)),
                by = by_cols]
  if (!nrow(pg)) return(NULL)

  # identical to script 3's _1 loop, but on this predator's table ----
  for (x in thresholds) {

    col_agg_level <- paste0("tax_level_", x, "_1")
    col_agg_name  <- paste0("prey_category_", x, "_1")

    # Keep the finest rank reaching x stomachs; n_predator_threshold = 1 makes
    # the pooled ">= 2 predators" test vacuous (correctif).
    pg <- PreyCategory(
      dt = pg,
      target_ranks = chosen_ranks,
      full_taxonomy = ranknfile,
      target_rows = pg$phylum %in% phyla_list,
      n_stomach_threshold = x,
      n_predator_threshold = 1,
      agg_level_var = col_agg_level,
      agg_name_var  = col_agg_name
    )

    # Non-identified prey keep their own prey_category (coalesce step).
    pg[, (col_agg_name) := fcase(
      prey_category == "identified_prey_species", get(col_agg_name),
      default = as.character(prey_category)
    )]
    pg[, (col_agg_level) := fcase(
      prey_category == "identified_prey_species", get(col_agg_level),
      default = "prey_category"
    )]

    # No species -> genus collapse for _1 columns (correctif).

    # A whole phylum still under x stomachs (within THIS predator) -> other_phyla.
    invalid_phyla <- pg[prey_category == "identified_prey_species",
                        .(total_stomach = uniqueN(unlist(list_stomach_id))),
                        by = phylum
    ][total_stomach < x, phylum]

    if (length(invalid_phyla) > 0) {
      pg[get(col_agg_level) == "phylum" & phylum %in% invalid_phyla,
         (col_agg_name) := "other_phyla"]
    }
  }

  pg[, predator := p]
  pg
}), fill = TRUE)

setcolorder(prey_groups_by_predator,
            c("predator", "prey_species_common_name", "prey_category", "tax_level"))

list(
  n_rows       = nrow(prey_groups_by_predator),
  n_predators  = uniqueN(prey_groups_by_predator$predator),
  n_thresh_col = sum(grepl("^prey_category_\\d+_1$", names(prey_groups_by_predator))),
  n_cols_total = ncol(prey_groups_by_predator)
)
#==============================================================================#
#==============================================================================#
#==============================================================================#
#==============================================================================#
