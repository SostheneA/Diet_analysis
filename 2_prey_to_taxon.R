#------------------------------------------------------------------------------#
# 2_prey_to_taxon.R - PREY CATEGORIES AND TAXONOMIC HIERARCHY (WoRMS)
#
# Input  : data/diet_corr.RData (from 1_diet_corrections.R)
#          R_helpers/taxonomic_rank_order.R
# Output : data/prey_taxon_data.rda (species_prey_id: prey with a WoRMS
#          classification; species_prey_noid: prey left unclassified)
#
# Prey records are first sorted into broad categories (identified prey,
# egg/larvae, digested, parasite, plantae/chromista, empty), then the
# identified ones are matched to WoRMS by AphiaID (gulf::species.names) or,
# failing that, by latin name, and their full hierarchy is retrieved.
#------------------------------------------------------------------------------#

rm(list = ls())

project_path <- here::here()

library(data.table)
library(taxize)
library(worrms)
library(gulf)
library(here)

source(paste0(project_path, "/R_helpers/taxonomic_rank_order.R"))
load(paste0(project_path, "/data/diet_corr.RData"))
setDT(diet_corr)

# 1. One row per prey name, with its sampling footprint -------------------------
species_prey <- diet_corr[,
                          .(
                            n_prey_count_per_prey = sum(number_of_prey, na.rm = TRUE),
                            n_stomach_per_prey = uniqueN(stomach_id, na.rm = TRUE),
                            n_predator_per_prey = uniqueN(predator_species_latin_name, na.rm = TRUE),
                            list_stomach_id       = list(unique(na.omit(stomach_id))),
                            list_predator_id      = list(unique(na.omit(predator_species_code)))
                          ),
                          by = .(
                            code = prey_species_code,
                            prey_species_common_name,
                            prey_species_latin_name
                          )
]


# 2. Broad Categorization & Manual Aphia Overrides -----------------------------
#------------------------------------------------------------------------------#

## a) define prey_category lists --------------------------------------------------
Not_list <- c("")


Inorganic_list <- c("Sand",
                  "Stones and rocks",
                  "Foreign articles,
                  garbage",
                  "Mud",
                  "Foreign articles,garbage"
)

Organic_list <- c("Mucus",
               "Organic debris",
               "Unid remains,digested",
               "Unid fish and invertebrates",
               "Unid. species",
               "Unidentified",
               "no prey information",
               "Odobenidae (f.)"   # walrus code in a sea raven stomach: not a prey
)

Unid_fish_list <- c("Unid fish and remains",
               "Fish remains",
               "Unid. Finfish",
               "Unid. fish (larvae,juvenile and adults)",
               "Unid fish and eggs")

Unid_invert_list <- c("Unid. marine invertebrata",
              "Ctenophora,coelenterata,porifera",
              "Protochordata sp.",
              "Zooplankton",
              "Invertebrate"
              )

Eggs_larvae_list <- list("Eggs - invertebrate" = c("Invertebrate", "egg"),
                         "Eggs - Invertebrate" = c("Invertebrate", "egg"),
                         "Eggs - unid." = c(NA, "egg"),
                         "Eggs - Unid." = c(NA, "egg"),
                         "Unid. Fish eggs" = c("Teleostei", "egg"),
                         "Eggs - decapoda" = c("Decapoda", "egg"),
                         "Eggs - Decapoda" = c("Decapoda", "egg"),
                         "Hermit crab eggs" = c("Paguridae", "egg"),
                         "Buccinidae sp. eggs" = c("Buccinidae", "egg"),
                         "Eggs - Atlantic herring"= c("Clupea harengus", "egg"),
                         "Larvae - unid." = c(NA, "larvae"),
                         "Larvae - Unid." = c(NA, "larvae"),
                         "Larvae - decapoda (o.)" = c("Decapoda", "larvae"),
                         "Larvae - Decapoda (o.)" = c("Decapoda", "larvae"),
                         "Larvae - crustacea (s.p.)" = c("Crustacea", "larvae"),
                         "Larvae - Crustacea (s.p.)" = c("Crustacea", "larvae"),
                         "Larvae - mollusca (p.)" = c("Mollusca", "larvae"),
                         "Larvae - Mollusca (p.)" = c("Mollusca", "larvae"),
                         "Larvae - polychaeta (c.)" = c("Polychaeta", "larvae"),
                         "Larvae - Polychaeta (c.)" = c("Polychaeta", "larvae"),
                         "Larvae - crab" = c("Pleocyemata", "larvae"),
                         "Larvae - Crab" = c("Pleocyemata", "larvae"),
                         "Larvae - homarus americanus" = c("Homarus americanus",
                                                           "larvae"),
                         "Larvae - Homarus americanus" = c("Homarus americanus",
                                                           "larvae"),
                         "Calappa megalops" = c("Calappa", "megalopae")
)

Parasite_list <- c("Parasites,round worms",
                   "Cestoda (c.)",
                   "Trematoda")

Rest_list <- c("Scales - unid.",
               "Scales - unid. fish",
               "Shells - unid. fish",
               "Mollusca sp. empty",
               "Casing - unidentified invertebrate",
               "Casing - Unidentified Invertebrate",
               "Scales - Unid.",
               "Scales - Unid. fish",
               "Shells - Unid. fish")

Algae_list <- c("Seaweed/Algae/Kelp unidentified",
                "Brown seaweed unidentified",
                "foraminiferans / forams / hole bearers",
                "Brown rockweed unidentified",
                "Red seaweed unidentified",
                "Eel grass"
                )


## b) apply the logic to the dataframe --------------------------------------------

### egg and larvae association to have the taxonomic hierarchy -----
map_latin <- setNames(vapply(Eggs_larvae_list, `[`, character(1L), 1L),
                      names(Eggs_larvae_list))
map_stage <- setNames(vapply(Eggs_larvae_list, `[`, character(1L), 2L),
                      names(Eggs_larvae_list))

setDT(species_prey)
species_prey[prey_species_latin_name %in% names(Eggs_larvae_list),
             `:=`(
               prey_species_latin_name = map_latin[prey_species_latin_name],
               life_stage              = map_stage[prey_species_latin_name],
               code = NA_integer_
             )]

### all fish unidentified in Teleostei -----
species_prey[prey_species_latin_name %in% Unid_fish_list,
             `:=`(
               prey_species_latin_name = "Teleostei",
               code = NA_integer_
             )]

### prey category -----
species_prey[, prey_category := fcase(
  prey_species_latin_name %in% Not_list, "empty",
  prey_species_latin_name %in% Inorganic_list, "digested",
  prey_species_latin_name %in% Organic_list, "digested",
  prey_species_latin_name %in% Unid_invert_list, "digested",
  (!is.na(life_stage)), "egg_larvae",
  (!is.na(life_stage) & is.na(prey_species_latin_name)), "digested",
  prey_species_latin_name %in% Rest_list, "digested",
  prey_species_latin_name %in% Parasite_list, "parasite",
  default = "identified_prey_species"
)]


species_prey[is.na(prey_species_latin_name) & is.na(prey_species_common_name),
          `:=`(
            prey_category  = "empty")
]

### Plantae_chromista particularity --------------
species_prey$prey_category <- ifelse(species_prey$prey_species_common_name %in%  Algae_list,
                                     "plantae_chromista",
                                     species_prey$prey_category)


## c) apply manual aphia_id overrides for some groups -----------------------------
all.sp <- as.data.table(gulf::species.names)
all.sp[, aphia_id_origin := aphia_id]

### i) aphia ID for the survey's composite codes ---------------------------------
all.sp[code == 290, aphia_id_origin := 126417] # Clupeidae/Osmeridae -> Clupea harengus
all.sp[code == 2699, aphia_id_origin := 1086] # Euphausiacea/Mysida -> Eumalacostraca
all.sp[code == 2902, aphia_id_origin := 1080] # copepoda small -> Copepoda
all.sp[code == 2901, aphia_id_origin := 1080] # copepoda large -> Copepoda
all.sp[code == 4998, aphia_id_origin := 325342] # squid beak -> Decapodiformes
all.sp[code == 5200, aphia_id_origin := 382158] # Limpet unidentified -> Patellogastropoda
all.sp[code == 2507, aphia_id_origin := 106670] # Crabs (could be Anomura or Brachyura) -> Pleocyemata
all.sp[code == 3199, aphia_id_origin := 883] # Polychaete sp. remains -> Polychaeta
all.sp[code == 3101, aphia_id_origin := 883] # Polychaeta (c.) large
all.sp[code == 3102, aphia_id_origin := 883] # Polychaeta (c.) small
all.sp[code == 4310, aphia_id_origin := 105] # Protobranchia/heteroconchia <- Bivalvia
all.sp[code == 6700, aphia_id_origin := 123111] # Psoluse/thyone <- Dendrochirotida

### ii) merge aphia_id in species_prey -----------------------------------------
species_prey <- merge(all.sp[, .(code, aphia_id_origin)],
                      species_prey,
                      by = "code",
                      all.y = T
)

#### manage prey small and large
species_prey[, prey_size_cat := NA_character_]
species_prey[, prey_size_cat := fifelse(code %in% c(3102, 2902), "small", prey_size_cat)]
species_prey[, prey_size_cat := fifelse(code %in% c(3101, 2901), "large", prey_size_cat)]


# 3. Taxonomic Lookup (WoRMS via AphiaID) -------------------------------------
#------------------------------------------------------------------------------#
## find those that have aphia_id -------------------------
aphia_id_origins  <- unique(species_prey[aphia_id_origin > 0, aphia_id_origin])

## look if all aphia_ids are valid.----------------------------

# worrms accepts short vectors only: query in batches of 50.
aphia_id_origins_batches <- split(aphia_id_origins, ceiling(seq_along(aphia_id_origins) / 50))


aphia_ids_verified_list <- lapply(aphia_id_origins_batches, function(batche) {
  tryCatch({
    Sys.sleep(0.2)
    return(wm_record(id = batche))
  }, error = function(e) {
    warning("Error on this batch of ids : ", paste(head(batche), collapse=", "))
    return(NULL)
  })
})

aphia_ids_verified_df <- rbindlist(aphia_ids_verified_list, fill = TRUE)
setDT(aphia_ids_verified_df)
aphia_ids_verified_df[, valid_AphiaID := fifelse(is.na(valid_AphiaID), parentNameUsageID, valid_AphiaID)]
aphia_id_origins_verified <- unique(aphia_ids_verified_df$valid_AphiaID)


## a) Batch retrieve classification aphia_id, WoRMS----------------------------------
if(!exists("class_list_prey_w")){
class_list_prey_w <- taxize::classification(aphia_id_origins_verified, db = "worms", return_id = TRUE)
}
class_tbl_prey  <- rbindlist(class_list_prey_w[!is.na(class_list_prey_w)],
                             fill = TRUE, idcol = "aphia_id")

### i) Hierarchy processing ----------------------------------------------------
class_tbl_prey[, rank := tolower(rank)]
class_tbl_prey[rank == "phylum (division)", rank := "phylum"]
class_tbl_prey[rank == "subphylum (subdivision)", rank := "subphylum"]
class_tbl_prey[, rank := factor(rank, levels = tolower(ranknfile), ordered = TRUE)]

#### Extract lowest rank per AphiaID
class_tbl_prey[, ":="(verified_name = .SD[rank == max(rank), name],
                      tax_level = .SD[, max(rank)]),
               by = aphia_id]

#### Pivot to wide format
class_tbl_prey <- dcast(class_tbl_prey,
                        aphia_id + verified_name + tax_level ~ rank,
                        value.var = "name")

setDT(class_tbl_prey)

### ii) merge valid code and classification with maybe not valid former id.
class_tbl_prey <- merge(
  x = class_tbl_prey[, aphia_id := as.integer(aphia_id)],
  y = aphia_ids_verified_df[, .(aphia_id_origin = as.integer(AphiaID), aphia_id = as.integer(valid_AphiaID))],
  by = "aphia_id",
  all.y = TRUE
)


### iii) merge WoRMS aphia_id and hierarchy -------------------------------------
species_prey_aphia <- merge(
  x = species_prey,
  y = class_tbl_prey,
  by = "aphia_id_origin",
  all.y = TRUE
)

aphia_id_origins_rest <- setdiff(aphia_id_origins, unique(class_tbl_prey$aphia_id_origin))
if (length(aphia_id_origins_rest)) message("AphiaIDs without classification: ",
                                           paste(aphia_id_origins_rest, collapse = ", "))


## b) Batch retrieve classification latin names, WoRMS --------------------------
## find those that have'nt aphia_id but latin_names -------------------------
latin_names  <- unique(species_prey[is.na(aphia_id_origin) &
                                                !is.na(prey_species_latin_name) &
                                                !(prey_species_latin_name %in% c("Obsolete", "Parasites,round worms")) &
                                                prey_category %in% c("identified_prey_species", "plantae_chromista", "parasite", "egg_larvae"),
                                              prey_species_latin_name])


latin_names_batches <- split(latin_names, ceiling(seq_along(latin_names) / 50))

latin_names_verified_list <- lapply(latin_names_batches, function(batche) {
  tryCatch({
    Sys.sleep(0.2)
    return(wm_records_names(name = batche))
  }, error = function(e) {
    warning("Error on this batch of names : ", paste(head(batche), collapse=", "))
    return(NULL)
  })
})

latin_names_verified_df <- dplyr::bind_rows(latin_names_verified_list)
setDT(latin_names_verified_df)
latin_names_verified_df[, valid_AphiaID := fifelse(is.na(valid_AphiaID), parentNameUsageID, valid_AphiaID)]
latin_names_id_verified <- unique(latin_names_verified_df$valid_AphiaID)


if(!exists("class_list_prey_w_latin")){
  class_list_prey_w_latin <- taxize::classification(latin_names_id_verified,
                                                    db = "worms",
                                                    return_id = TRUE)
}


class_tbl_prey2  <- rbindlist(class_list_prey_w_latin[!is.na(class_list_prey_w_latin)],
                             fill = TRUE, idcol = "aphia_id")


### i) Hierarchy processing ----------------------------------------------------
class_tbl_prey2[, rank := tolower(rank)]
class_tbl_prey2[rank == "phylum (division)", rank := "phylum"]
class_tbl_prey2[rank == "subphylum (subdivision)", rank := "subphylum"]
class_tbl_prey2[, rank := factor(rank, levels = tolower(ranknfile), ordered = TRUE)]

#### Extract lowest rank per AphiaID
class_tbl_prey2[, ":="(verified_name = .SD[rank == max(rank), name],
                      tax_level = .SD[, max(rank)]),
               by = aphia_id]

#### Pivot to wide format
class_tbl_prey2 <- dcast(class_tbl_prey2,
                        aphia_id + verified_name + tax_level ~ rank,
                        value.var = "name")

setDT(class_tbl_prey2)

### ii) merge valid code and classification with maybe not valid former latin_names.
class_tbl_prey2 <- merge(
  x = class_tbl_prey2[, aphia_id := as.integer(aphia_id)],
  y = latin_names_verified_df[, .(prey_species_latin_name = as.character(scientificname),
                                  aphia_id = as.integer(valid_AphiaID))],
  by = "aphia_id",
  all.y = TRUE
)


### iii) merge WoRMS aphia_id and hierarchy ------------------------------------
species_prey_latin <- merge(
  x = subset(species_prey, is.na(aphia_id_origin)),
  y = class_tbl_prey2,
  by = "prey_species_latin_name",
  all.y = T
)


latin_names_rest <- setdiff(latin_names, unique(class_tbl_prey2$prey_species_latin_name))
if (length(latin_names_rest)) message("Latin names without classification: ",
                                      paste(latin_names_rest, collapse = ", "))


# 4) Merge identified preys together "Id", and find the remain ones "NoId" -----
species_prey_id <- rbind(species_prey_aphia, species_prey_latin, fill = TRUE)

species_prey_noid_list <- setdiff(unique(species_prey$prey_species_latin_name), unique(species_prey_id$prey_species_latin_name))
species_prey_noid <- species_prey[species_prey$prey_species_latin_name %in% species_prey_noid_list,]

cat("Classified prey names: ", uniqueN(species_prey_id$verified_name),
    " | unclassified: ", nrow(species_prey_noid), "\n", sep = "")

# 5) Export --------------------------------------------------------------------
save(species_prey_id, species_prey_noid, file = paste0(project_path, "/data/prey_taxon_data.rda"))

