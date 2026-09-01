#------------------------------------------------------------------------------#
# 3PP_taxonomic_groups.R - PER-PREDATOR PREY GROUPING (FAMILY "_PP")
#
# Same rule as the _1 family of script 3 (stomach threshold only, no species
# -> genus collapse), but the threshold sweep is rebuilt inside each main
# predator's own stomachs, so each predator gets its own grouping vectors.
#   - Predators : >= 100 stomachs sampled (the study set).
#   - Thresholds: x = 5, 10, ..., 250 (50 vectors per predator).
#   - Columns   : prey_category_<x>_PP / tax_level_<x>_PP.
#
# Input  : data/prey_groups.RData (taxonomy, from script 3), data/diet_corr.RData
#          R_helpers/taxonomic_rank_order.R, R_helpers/PreyCategory.R
# Output : data/prey_groups_PP.RData (prey_groups_PP: one row per predator x
#          prey taxon). Script 4 merges it on (predator, prey_species_common_name);
#          6b-6d run with PREY_FAMILY <- "_PP" (see 6f_Run_families.R).
#------------------------------------------------------------------------------#

project_path <- here::here()
library(data.table)

source(paste0(project_path, "/R_helpers/taxonomic_rank_order.R"))
source(paste0(project_path, "/R_helpers/PreyCategory.R"))
load(paste0(project_path, "/data/prey_groups.RData"))
load(paste0(project_path, "/data/diet_corr.RData"))
setDT(prey_groups); setDT(diet_corr)

# 1) Parameters ----------------------------------------------------------------
PP_SUFFIX          <- "_PP"
thresholds_PP      <- seq(5, 250, by = 5)
MAIN_PRED_MIN_STOM <- 100L

chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")
phyla_list   <- unique(prey_groups[kingdom == "Animalia" & !is.na(phylum), phylum])

# 2) Taxonomic lookup: one row per prey name ----------------------------------
rank_cols <- intersect(
  c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","parvphylum",
    "gigaclass","superclass","class","subclass","infraclass","subterclass","superorder",
    "order","suborder","infraorder","parvorder","section","subsection","superfamily",
    "family","subfamily","tribe","genus","subgenus","species"),
  names(prey_groups))

lookup_cols <- c("prey_species_common_name", "prey_category", "verified_name",
                 "tax_level", rank_cols)
tax_lookup  <- unique(prey_groups[, ..lookup_cols], by = "prey_species_common_name")

# 3) Main predators and their stomach x prey records ---------------------------
main_predators <- diet_corr[!is.na(stomach_id),
                            .(n_sto = uniqueN(stomach_id)),
                            by = predator_species_common_name
][n_sto >= MAIN_PRED_MIN_STOM, predator_species_common_name]
message(length(main_predators), " main predators")

rec_all <- diet_corr[!is.na(stomach_id) & predator_species_common_name %in% main_predators,
                     .(predator = predator_species_common_name, stomach_id, prey_species_common_name)]
rec_all <- tax_lookup[rec_all, on = "prey_species_common_name", nomatch = 0L]

by_cols <- c("prey_species_common_name", "prey_category", "verified_name",
             "tax_level", rank_cols)

# 4) Threshold sweep inside each predator --------------------------------------
prey_groups_PP <- rbindlist(lapply(main_predators, function(p) {

  # One row per prey taxon of this predator; list_predator_id = {p}, so the
  # predator criterion of PreyCategory is 1 by construction.
  pg <- rec_all[predator == p,
                .(list_stomach_id  = list(unique(stomach_id)),
                  list_predator_id = list(p)),
                by = by_cols]
  if (!nrow(pg)) return(NULL)

  is_prey     <- !is.na(pg$prey_category) & pg$prey_category == "identified_prey_species"
  pc_chr      <- as.character(pg$prey_category)
  phylum_vec  <- as.character(pg$phylum)
  target_rows <- pg$phylum %in% phyla_list
  valid_cols  <- intersect(ranknfile[ranknfile %in% chosen_ranks], names(pg))
  counts      <- rank_counts(pg, rev(valid_cols), target_rows)
  phylum_tot  <- pg[is_prey, .(total_stomach = uniqueN(unlist(list_stomach_id))), by = phylum]

  for (x in thresholds_PP) {
    r <- PreyCategory(pg, chosen_ranks, ranknfile, target_rows = target_rows,
                      n_stomach_threshold = x, n_predator_threshold = 1, counts = counts)
    name <- r$name; level <- r$level
    name[!is_prey]  <- pc_chr[!is_prey]
    level[!is_prey] <- "prey_category"
    bad <- phylum_tot[total_stomach < x, phylum]
    idx <- which(level == "phylum" & phylum_vec %in% bad)
    name[idx] <- "other_phyla"
    set(pg, j = paste0("prey_category_", x, PP_SUFFIX), value = name)
    set(pg, j = paste0("tax_level_", x, PP_SUFFIX),     value = level)
  }

  pg[, predator := p]
  message("[PP] ", p, " (", nrow(pg), " taxa)")
  pg
}), fill = TRUE)

setcolorder(prey_groups_PP,
            c("predator", "prey_species_common_name", "prey_category", "tax_level"))

# 5) Checks and export ----------------------------------------------------------
n_pp_cols <- sum(grepl(paste0("^prey_category_\\d+", PP_SUFFIX, "$"), names(prey_groups_PP)))
stopifnot(n_pp_cols == length(thresholds_PP))

cat(sprintf("\nprey_groups_PP: %d rows | %d predators | %d resolution columns\n",
            nrow(prey_groups_PP), uniqueN(prey_groups_PP$predator), n_pp_cols))

save(prey_groups_PP, file = paste0(project_path, "/data/prey_groups_PP.RData"))
