#------------------------------------------------------------------------------#
# 3PP_taxonomic_groups.R
#
# Objective : PER-PREDATOR prey grouping database ("PP" family).
#   Same logic as script 3's `_1` family (n_predator_threshold = 1, NO species
#   -> genus collapse, invalid_phyla on the stomach threshold only), but the
#   whole threshold sweep is rebuilt INSIDE each main predator's own stomachs :
#   thresholds are evaluated on that predator's support only, so each predator
#   gets its own grouping vectors (over-resolution fix, cf. 3b/3c diagnostics).
#
#   - Predators : the study set, >= 100 stomachs sampled (expected : 40).
#   - Thresholds: X = seq(5, 250, by = 5)  -> 50 vectors per predator.
#   - Columns   : prey_category_<x>_PP / tax_level_<x>_PP  (suffix "_PP" so the
#                 pooled _1/_2 columns are never overwritten downstream).
#
# Input  : data/prey_groups.RData   (script 3 output ; taxonomy source)
#          data/diet_corr.RData
#          R_helpers/taxonomic_rank_order.R
# Output : data/prey_groups_PP.RData  (object : prey_groups_PP ; long table,
#          one row per predator x prey taxon)
#
# Downstream : script 4 merges it on (predator, prey_species_common_name) ;
#              script 5 unchanged ; 6+ run with PREY_FAMILY <- "_PP".
#------------------------------------------------------------------------------#

gc()

# 1) Environment & data --------------------------------------------------------
project_path <- here::here()
library(data.table)

source(paste0(project_path, "/R_helpers/taxonomic_rank_order.R"))
load(paste0(project_path, "/data/prey_groups.RData"))   # taxonomie (script 3)
load(paste0(project_path, "/data/diet_corr.RData"))
setDT(prey_groups); setDT(diet_corr)

# 2) PreyCategory (copie conforme du script 3) ---------------------------------
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

  # 1. Use integer indices for speed.
  #    all_target_idx is the FULL set of rows in scope and never shrinks; it is
  #    what the thresholds are evaluated on. rem_idx shrinks as rows get
  #    assigned and drives the assignment only.
  all_target_idx <- if (is.null(target_rows)) seq_len(nrow(dt)) else which(target_rows)
  rem_idx <- all_target_idx

  if (length(rem_idx) == 0) return(dt)

  # Initialize columns safely to avoid 'set' masking conflicts
  if (!agg_level_var %in% names(dt)) dt[, (agg_level_var) := NA_character_]
  if (!agg_name_var %in% names(dt))  dt[, (agg_name_var) := NA_character_]

  # Reset targeted rows to NA
  dt[rem_idx, c(agg_level_var, agg_name_var) := NA_character_]

  # Filter ranks existing in data and sort Bottom-Up (Specific to General)
  ordered_target_ranks <- full_taxonomy[full_taxonomy %in% target_ranks]
  valid_cols <- intersect(ordered_target_ranks, names(dt))
  spec_to_gen <- rev(valid_cols)

  # 2. Threshold evaluation, on ALL target rows, BEFORE any assignment.
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

  # 2b. Assignment, still bottom-up, but using the counts computed above.
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

  # 3. Fallback logic for remaining residual rows
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


# 3) Parametres ----------------------------------------------------------------
PP_SUFFIX          <- "_PP"
thresholds_PP      <- seq(5, 250, by = 5)      # 50 seuils
MAIN_PRED_MIN_STOM <- 100L                     # jeu d'etude (attendu : 40 pred.)

chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")
phyla_list   <- unique(prey_groups[kingdom == "Animalia" & !is.na(phylum), phylum])

# 4) Lookup taxonomique : une ligne par taxon-proie ----------------------------
rank_cols <- intersect(
  c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","parvphylum",
    "gigaclass","superclass","class","subclass","infraclass","subterclass","superorder",
    "order","suborder","infraorder","parvorder","section","subsection","superfamily",
    "family","subfamily","tribe","genus","subgenus","species"),
  names(prey_groups))

lookup_cols <- c("prey_species_common_name", "prey_category", "verified_name",
                 "tax_level", rank_cols)
tax_lookup  <- unique(prey_groups[, ..lookup_cols], by = "prey_species_common_name")

# 5) Predateurs principaux -----------------------------------------------------
main_predators <- diet_corr[!is.na(stomach_id),
                            .(n_sto = uniqueN(stomach_id)),
                            by = predator_species_common_name
][n_sto >= MAIN_PRED_MIN_STOM, predator_species_common_name]

message(length(main_predators), " predateurs principaux (attendu : 40)")

# 6) Enregistrements estomac x proie pour ces predateurs, taxonomie attachee ---
rec_all <- diet_corr[!is.na(stomach_id) &
                       predator_species_common_name %in% main_predators,
                     .(predator = predator_species_common_name,
                       stomach_id, prey_species_common_name)]
rec_all <- tax_lookup[rec_all, on = "prey_species_common_name", nomatch = 0L]

by_cols <- c("prey_species_common_name", "prey_category", "verified_name",
             "tax_level", rank_cols)

# 7) BOUCLE PREDATEUR (enveloppe le sweep de seuils, regles famille _1) --------
prey_groups_PP <- rbindlist(lapply(main_predators, function(p) {

  # Une ligne par taxon-proie DE CE predateur ; list_stomach_id alimente le
  # test de seuil de PreyCategory ; list_predator_id = {p}, donc n_predator = 1
  # par construction (le critere >= 2 predateurs est retire).
  pg <- rec_all[predator == p,
                .(list_stomach_id  = list(unique(stomach_id)),
                  list_predator_id = list(p)),
                by = by_cols]
  if (!nrow(pg)) return(NULL)

  for (x in thresholds_PP) {

    col_agg_level <- paste0("tax_level_", x, PP_SUFFIX)
    col_agg_name  <- paste0("prey_category_", x, PP_SUFFIX)

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

    # Harmonisation selon prey_category (identique aux boucles du script 3)
    pg[, (col_agg_name) := fcase(
      prey_category == "identified_prey_species", get(col_agg_name),
      default = as.character(prey_category)
    )]
    pg[, (col_agg_level) := fcase(
      prey_category == "identified_prey_species", get(col_agg_level),
      default = "prey_category"
    )]

    # PAS de collapse espece -> genre (regle famille _1).

    # Phyla sous le seuil AU SEIN de ce predateur : seuil d'estomacs seul.
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
  message("[PP] termine : ", p, "  (", nrow(pg), " taxons)")
  pg
}), fill = TRUE)

setcolorder(prey_groups_PP,
            c("predator", "prey_species_common_name", "prey_category", "tax_level"))

# 8) Controles & sauvegarde ----------------------------------------------------
n_pp_cols <- sum(grepl(paste0("^prey_category_\\d+", PP_SUFFIX, "$"),
                       names(prey_groups_PP)))
stopifnot(n_pp_cols == length(thresholds_PP))          # attendu : 50

cat(sprintf("\nprey_groups_PP : %d lignes | %d predateurs | %d paires de colonnes %s\n",
            nrow(prey_groups_PP), uniqueN(prey_groups_PP$predator),
            n_pp_cols, PP_SUFFIX))

save(prey_groups_PP, file = paste0(project_path, "/data/prey_groups_PP.RData"))
message("Sauve : data/prey_groups_PP.RData")
