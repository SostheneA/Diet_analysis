#------------------------------------------------------------------------------#
# Objective: Find a logic to do groups and
#           reduce redundancy in taxonomic hierarchy for the analysis
#
#          - First group are the lowest possible with the available taxonomic
#               hierarchy
#          - Second are a funtion that may be based on the folowing rules
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



## b) Identified Prey handling -------------------------------------------------
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


############  loop


# --- PRÉPARATION INITIALE (Une seule fois) ---

# Définition des rangs cibles
chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")

# Création de la liste des phylums (Animalia non NA)
# On la définit ici pour qu'elle soit disponible pour la boucle
phyla_list <- unique(prey_groups[kingdom == "Animalia" & !is.na(phylum), phylum])

# --- BOUCLE D'AUTOMATISATION (100 variables) ---

# Séquence de 10 à 1000 (incrément de 10)
thresholds <- seq(10, 1000, by = 10)

for (x in thresholds) {

  # 1. Noms dynamiques des colonnes
  col_agg_level <- paste0("tax_level_", x, "_2")
  col_agg_name  <- paste0("prey_category_", x, "_2")

  # 2. Appel de la fonction PreyCategory
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

  # 3. Harmonisation selon prey_category (fcase)
  # On utilise (col) := pour l'assignation dynamique
  prey_groups[, (col_agg_name) := fcase(
    prey_category == "identified_prey_species", get(col_agg_name),
    default = as.character(prey_category)
  )]

  prey_groups[, (col_agg_level) := fcase(
    prey_category == "identified_prey_species", get(col_agg_level),
    default = "prey_category"
  )]

  # 4. Nettoyage Genus -> Species
  # Récupère les genres qui sont devenus des catégories finales
  current_genuses <- unique(prey_groups[get(col_agg_level) == "genus"][[col_agg_name]])

  if (length(current_genuses) > 0) {
    pattern <- paste0("\\b(", paste(current_genuses, collapse = "|"), ")\\b")

    # Mise à jour des espèces appartenant à ces genres
    prey_groups[get(col_agg_level) == "species" & grepl(pattern, get(col_agg_name)),
                c(col_agg_name, col_agg_level) := .(genus, "genus")]
  }

  # 5. Gestion des Phyla sous le seuil (Invalid Phyla)
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
  message(paste("Calcul terminé pour le seuil x =", x))
  gc()
}



for (x in thresholds) {
  col_name <- paste0("prey_category_", x, "_2")

  # Calcul du nombre de catégories uniques
  n_categories <- uniqueN(prey_groups[[col_name]])

  cat("\n--- Seuil x =", x, "(Colonne:", col_name, ") ---")
  cat("\nNombre total de catégories :", n_categories)

  # Affichage du top 5 des catégories les plus fréquentes
  print(prey_groups[, .N, by = col_name][order(-N)][1:5])
}


library(ggplot2)

# Création d'un tableau résumé
summary_list <- lapply(thresholds, function(x) {
  col_name <- paste0("prey_category_", x, "_2")
  data.table(
    threshold = x,
    n_unique_categories = uniqueN(prey_groups[[col_name]]),
    n_items_total = nrow(prey_groups)
  )
})

summary_dt <- rbindlist(summary_list)

# Affichage du tableau de synthèse
print(summary_dt)

# Petit bonus : Visualiser la perte de résolution taxonomique
ggplot(summary_dt, aes(x = threshold, y = n_unique_categories)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point() +
  labs(title = "Évolution du nombre de catégories de proies",
       subtitle = "Plus le seuil est haut, moins on a de catégories différentes",
       x = "Seuil (n_stomach_threshold)",
       y = "Nombre de catégories uniques") +
  theme_minimal()+
  scale_x_continuous(breaks = seq(0, 1000, by = 40))+scale_y_continuous(breaks = seq(0, 1000, by = 10))



ggplot(summary_dt, aes(x = threshold, y = n_unique_categories)) +
  # Utilisation de geom_col pour un bar plot (identique à geom_bar(stat="identity"))
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(
    title = "Évolution du nombre de catégories de proies",
    subtitle = "Réduction de la diversité taxonomique en fonction du seuil d'agrégation",
    x = "Seuil (n_stomach_threshold)",
    y = "Nombre de catégories uniques"
  ) +
  theme_minimal() +
  # Garder tes échelles personnalisées
  scale_x_continuous(breaks = seq(0, 1000, by = 40)) +
  scale_y_continuous(breaks = seq(0, 1000, by = 10)) +
  # Optionnel : incliner les étiquettes x si elles se chevauchent
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# 4) prey_category_keep-----------------------------------------------------------
#------------------------------------------------------------------------------#
chosen_ranks <- c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                  "family", "genus", "species")
phyla_list <- unique(prey_groups[prey_groups$kingdom == "Animalia" &
                                   !is.na(prey_groups$phylum),]$phylum)
dput(phyla_list)

prey_groups <- PreyCategory(dt = prey_groups,
                            target_ranks = c("phylum", "subphylum", "class", "subclass", "order", "infraorder",
                                             "family", "genus", "species"),
                            full_taxonomy = ranknfile,
                            target_rows = prey_groups$phylum %in% phyla_list,
                            n_stomach_threshold = 440,
                            n_predator_threshold = 2,
                            agg_level_var = "tax_level_keep",
                            agg_name_var = "prey_category_keep")


prey_groups[, ':='(
  prey_category_keep = fcase(
    prey_category == "identified_prey_species", as.character(prey_category_keep),
    prey_category != "identified_prey_species", as.character(prey_category),
    default = as.character(prey_category)
  ),
  tax_level_keep = fcase(
    prey_category == "identified_prey_species", as.character(tax_level_keep),
    prey_category != "identified_prey_species", "prey_category",
    default = as.character("prey_category")
  )
)]


genuses <- unique(prey_groups[prey_groups$tax_level_keep == "genus",]$prey_category_keep)
pattern <- paste0("\\b(", paste(genuses, collapse = "|"), ")\\b")
prey_groups[tax_level_keep == "species" & grepl(pattern, prey_category_keep),
            ":="(
            prey_category_keep = genus,
            tax_level_keep = "genus"
            )]

# if the phylum did not respect the threshold
invalid_phyla <- prey_groups[prey_category == "identified_prey_species",
                             .(total_stomach = uniqueN(unlist(list_stomach_id)),
                               total_pred    = uniqueN(unlist(list_predator_id))),
                             by = phylum
][total_stomach < 440 | total_pred < 2, phylum]

prey_groups[tax_level_keep == "phylum" & phylum %in% invalid_phyla,
            prey_category_keep := "other_phyla"]



# 5) --- Dynamic Detection: Unique Label Rule for "Others" Groups ---

# 1. Initialization
# 'has_child' flags groups that need a suffix because they have specific sub-groups.
# 'prey_category_new' will hold the final publication-ready names.
prey_groups[, has_child := FALSE]
prey_groups[, prey_category_new := prey_category_keep]

# 2. Identify category names based on taxonomy
# We exclude generic categories (e.g., "digested", "parasite") to focus on biological ranks.
taxo_labels <- unique(prey_groups[tax_level_keep != "prey_category", prey_category_keep])

# 3. Verification Loop
for (cat_name in taxo_labels) {

  # CORE LOGIC: Find all final labels that belong to the taxonomic group 'cat_name'.
  # We check across all taxonomic levels from Phylum to Genus.
  # We only count labels that are taxonomic (tax_level_keep != "prey_category").
  labels_under_parent <- prey_groups[
    (phylum == cat_name | class == cat_name | order == cat_name |
       family == cat_name | genus == cat_name) &
      tax_level_keep != "prey_category",
    unique(prey_category_keep)
  ]

  # THE SIBLING RULE:
  # If only 1 unique label exists (e.g., only "Gastropoda"), then 'has_child' remains FALSE.
  # If 2 or more exist (e.g., "Mollusca" AND "Gastropoda"), then "Mollusca" acts
  # as a "catch-all" for the remainder and must become "Mollusca_others".
  if (length(labels_under_parent) > 1) {
    prey_groups[prey_category_keep == cat_name, has_child := TRUE]
  }
}

# 4. Suffix Application
# Apply the "_others" suffix only to the "parent" rows that were flagged.
prey_groups[has_child == TRUE, prey_category_new := paste0(prey_category_keep, "_others")]

# --- Verification ---
unique(prey_groups[has_child == TRUE, .(prey_category_keep, prey_category_new)])

table(prey_groups$has_child, prey_groups$tax_level_keep)
table(prey_groups$has_child, prey_groups$tax_level_keep, prey_groups$phylum)


unique(prey_groups$prey_category)
table(prey_groups$prey_category)
length(unique(prey_groups$prey_category_lowest))
length(unique(prey_groups$prey_category_keep))
nrow(prey_groups)

sort(unique(prey_groups$prey_category_new))




# 5) Final export---------------------------------------------------------------
save(prey_groups, file = paste0(project_path, "/data/prey_groups.RData"))
