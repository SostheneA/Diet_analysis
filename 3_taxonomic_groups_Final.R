#------------------------------------------------------------------------------#
# 3_taxonomic_groups_Final.R - PREY AGGREGATION SWEEP (POOLED FAMILIES _2 AND _1)
#
# Input  : data/prey_taxon_data.rda (from 2_prey_to_taxon.R)
#          R_helpers/taxonomic_rank_order.R, R_helpers/PreyCategory.R
# Output : data/prey_groups.RData (prey_groups: one row per prey name with
#          prey_category_<x>_2 / tax_level_<x>_2 and prey_category_<x>_1 /
#          tax_level_<x>_1 for x = 10, 20, ..., 1000, plus prey_category_keep
#          (x = 440) and prey_category_new (sibling rule applied))
#
# Rule: each prey record keeps its lowest rank, or walks up the hierarchy
# (species -> genus -> family -> infraorder -> order -> subclass -> class ->
# subphylum -> phylum) until it reaches a rank whose name occurs in at least x
# stomachs, thresholds being evaluated on the complete data before assignment.
# Phyla below the threshold are pooled as "other_phyla".
#   _2 : name must also occur in >= 2 predator species; species collapse to
#        their genus when the genus is itself a category.
#   _1 : no predator criterion, no species -> genus collapse (manuscript).
# Life stages (egg/larvae) are kept as categories; parasites, digested and
# empty records keep their prey_category label at every resolution.
#
# The name counts behind the thresholds are computed once (they do not depend
# on x); each resolution then only assigns and writes two columns with set().
#------------------------------------------------------------------------------#

project_path <- here::here()

library(data.table)

source(paste0(project_path, "/R_helpers/taxonomic_rank_order.R"))
source(paste0(project_path, "/R_helpers/PreyCategory.R"))
load(paste0(project_path, "/data/prey_taxon_data.rda"))

prey_groups <- rbind(species_prey_id, species_prey_noid, fill = TRUE)
data.table::setDT(prey_groups)

# 1) Lowest available classification -------------------------------------------
is_prey <- !is.na(prey_groups$prey_category) & prey_groups$prey_category == "identified_prey_species"
pc_chr  <- as.character(prey_groups$prey_category)

data.table::set(prey_groups, j = "prey_category_lowest",
                value = ifelse(is_prey & !is.na(prey_groups$verified_name), as.character(prey_groups$verified_name), pc_chr))
data.table::set(prey_groups, j = "tax_level_lowest",
                value = ifelse(is_prey & !is.na(prey_groups$verified_name), as.character(prey_groups$tax_level), "prey_category"))

# 2) Shared inputs of the sweep ---------------------------------------------------
chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")
valid_cols   <- intersect(ranknfile[ranknfile %in% chosen_ranks], names(prey_groups))
phyla_list   <- unique(prey_groups[kingdom == "Animalia" & !is.na(phylum), phylum])
target_rows  <- prey_groups$phylum %in% phyla_list
thresholds   <- seq(10, 1000, by = 10)

counts <- rank_counts(prey_groups, rev(valid_cols), target_rows)

phylum_tot <- prey_groups[is_prey,
                          .(total_stomach = uniqueN(unlist(list_stomach_id)),
                            total_pred    = uniqueN(unlist(list_predator_id))),
                          by = phylum]
genus_vec  <- as.character(prey_groups$genus)
phylum_vec <- as.character(prey_groups$phylum)

# One resolution: assignment, non-prey labels, optional genus collapse,
# phyla below the threshold, then two columns written with data.table::set().
build_resolution <- function(x, n_pred, collapse_genus, col_name, col_level) {
  r <- PreyCategory(prey_groups, chosen_ranks, ranknfile, target_rows = target_rows,
                    n_stomach_threshold = x, n_predator_threshold = n_pred, counts = counts)
  name <- r$name; level <- r$level

  name[!is_prey]  <- pc_chr[!is_prey]
  level[!is_prey] <- "prey_category"

  if (collapse_genus) {
    g   <- unique(name[!is.na(level) & level == "genus"])
    idx <- which(level == "species" & genus_vec %in% g)
    name[idx]  <- genus_vec[idx]
    level[idx] <- "genus"
  }

  bad <- if (n_pred >= 2) phylum_tot[total_stomach < x | total_pred < n_pred, phylum] else
    phylum_tot[total_stomach < x, phylum]
  idx <- which(level == "phylum" & phylum_vec %in% bad)
  name[idx] <- "other_phyla"

  data.table::set(prey_groups, j = col_name,  value = name)
  data.table::set(prey_groups, j = col_level, value = level)
  invisible(NULL)
}

# 3) Family _2: >= x stomachs and >= 2 predators, species collapse to genus ------
for (x in thresholds) {
  build_resolution(x, n_pred = 2, collapse_genus = TRUE,
                   col_name = paste0("prey_category_", x, "_2"), col_level = paste0("tax_level_", x, "_2"))
  message("[_2] x = ", x, " done")
}

# 4) Family _1: >= x stomachs only, no species -> genus collapse ----------------
for (x in thresholds) {
  build_resolution(x, n_pred = 1, collapse_genus = FALSE,
                   col_name = paste0("prey_category_", x, "_1"), col_level = paste0("tax_level_", x, "_1"))
  message("[_1] x = ", x, " done")
}

n_cat <- rbindlist(lapply(thresholds, function(x) data.table(
  x = x,
  n_2 = uniqueN(prey_groups[[paste0("prey_category_", x, "_2")]]),
  n_1 = uniqueN(prey_groups[[paste0("prey_category_", x, "_1")]]))))
cat("\nNumber of prey categories by threshold:\n")
print(n_cat[x %in% c(10, 50, 100, 200, 440, 750, 1000)])

# 5) prey_category_keep (x = 440, family _2 rule) -------------------------------
build_resolution(440, n_pred = 2, collapse_genus = TRUE,
                 col_name = "prey_category_keep", col_level = "tax_level_keep")

# 6) Sibling rule: a parent label that also has finer labels under it becomes
#    "<parent>_others" (prey_category_new).
keep_name  <- prey_groups$prey_category_keep
keep_level <- prey_groups$tax_level_keep
is_taxo    <- !is.na(keep_level) & keep_level != "prey_category"
has_child  <- rep(FALSE, nrow(prey_groups))

rank_mat <- lapply(c("phylum", "class", "order", "family", "genus"), function(cc) as.character(prey_groups[[cc]]))
for (cat_name in unique(keep_name[is_taxo])) {
  under <- Reduce(`|`, lapply(rank_mat, function(v) !is.na(v) & v == cat_name)) & is_taxo
  if (length(unique(keep_name[under])) > 1) has_child[!is.na(keep_name) & keep_name == cat_name] <- TRUE
}
data.table::set(prey_groups, j = "has_child", value = has_child)
data.table::set(prey_groups, j = "prey_category_new",
                value = ifelse(has_child, paste0(keep_name, "_others"), keep_name))

cat("\nCategories at x = 440: ", uniqueN(prey_groups$prey_category_keep),
    " (", uniqueN(prey_groups$prey_category_new), " with the sibling rule)\n", sep = "")

# 7) Export --------------------------------------------------------------------
save(prey_groups, file = paste0(project_path, "/data/prey_groups.RData"))
