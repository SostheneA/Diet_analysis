#------------------------------------------------------------------------------#
# objectives : assigning prey_catagory and/or taxonomic hierarchy to preys
#
# author : Sosthene Akia / LL
#
# Input :
#   - data/diet_corr from 1_diet_corrections.R
#   - R/taxonomic_rank_order.R", similar to source("R/helpers/taxonomy_levels.R")
#
# Output :
#  - species_prey_id, and species_prey_noid in /data/prey_taxon_data.rda
#
#
#
# TODO (assignment):
# Better grouping of preys or significant grouping for preys
# 40 mains predators (based on more than 100 stomachs threshold)
#
# Main prey or stomach elements groups
#   A - garbage (parasites, plants or old algae groups, land insects)
#   B - well identified preys
#   C - unidentified preys
#
#  Target date : February 20th 2026.
#  Adds to TODO
#  - functional groups (not possible in the allocated time),
#  - define empty,
#  - what to remove before analysis
#  - corrections to diet_raw
#  - larvae, plantae/chromista problem, big vs small item
#
#------------------------------------------------------------------------------#

rm(list = ls())
gc()


# 1. Environment & Data Loading ------------------------------------------------
#------------------------------------------------------------------------------#
project_path <- here::here()

## a) load libraries and helpers -----------------------------------------------
library(data.table)
library(taxize)
library(worrms)
library(gulf)
library(here)

# Probably similar to source("R/helpers/taxonomy_levels.R")
source(paste0(project_path, "/R_helpers/taxonomic_rank_order.R"))

# data - run 1_diet_corrections.R to obtain this
load(paste0(project_path, "/data/diet_corr.RData"))
setDT(diet_corr)
## c) Create a species_prey file  ----------------------------------------------

### ii) ---

# A CORRIGER CHEZ LYSANDRE
species_prey <- diet_corr[,
                          .(
                            #  n = .N, # no need of this see below for exploration
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


### i) exploration to know if my n n_... ok ----
#species_prey2 <- diet_corr[
#  ,
#  .(
#    n = .N,
#    n_prey_count_per_prey = sum(number_of_prey, na.rm = TRUE),
#    n_stomach_per_prey = uniqueN(stomach_id, na.rm = TRUE),
#    n_predator_per_prey = uniqueN(predator_species_latin_name, na.rm = TRUE)
#  ),
#  by = .(
#    cruise_number,
#    cruise_year,
#    set,
#    stomach_id,
#    predator_species_common_name,
#    code = prey_species_code,
#    prey_species_common_name,
#    prey_species_latin_name
#  )
#]
#
#species_prey3 <- diet_corr[
#  ,
#  .(
#    n = .N,
#    n_prey_count_per_prey = sum(number_of_prey, na.rm = TRUE),
#    n_stomach_per_prey = uniqueN(stomach_id, na.rm = TRUE),
#    n_predator_per_prey = uniqueN(predator_species_latin_name, na.rm = TRUE)
#  ),
#  by = .(
#    #cruise_number,
#    #cruise_year,
#    #set,
#    stomach_id,
#    predator_species_common_name,
#    code = prey_species_code,
#    prey_species_common_name,
#    prey_species_latin_name
#  )
#]

#rm(species_prey2, species_prey3)





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
               "Odobenidae (f.)" # error ?? walrus in "Hemitripterus americanus" stomach unless it's genetic
)

Unid_fish_list <- c("Unid fish and remains",
               "Fish remains",
               "Unid. Finfish",
               "Unid. fish (larvae,juvenile and adults)",
               "Unid fish and eggs")

Unid_invert_list <- c("Unid. marine invertebrata",
              "Ctenophora,coelenterata,porifera",
              #"Hyalellidae (f.)", # ?? see comment in diet file, one observation (Amphipoda (Order), but it has a problem in WoRMS with it.)
              "Protochordata sp.", # could be tunicata or something else, one observation
              "Zooplankton",
              "Invertebrate"
              )

Mixed_prey_list <- c(#"Clupeidae/osmeridae (f.)",
                    #"Euphausiacea/mysidacea o.",
                    #"Protobranchia/heteroconchia", # Bivalvia
                     #"Ctenophora,coelenterata,porifera"#, # old grouping Ctenophora/Cnidaria/Porifera Phylum
                    # "Psoluse/thyone" # Dendrochirotida
                    ) # just one stay here for now, the other are now classified, see below

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
                ) #c("Obsolete") - for latin name, Kingdom: Plantae, Chromista



## b) apply the logic to the dataframe --------------------------------------------

### egg and larvae association to have the taxonomic hierarchy -----
map_latin <- setNames(vapply(Eggs_larvae_list, `[`, character(1L), 1L),
                      names(Eggs_larvae_list))
map_stage <- setNames(vapply(Eggs_larvae_list, `[`, character(1L), 2L),
                      names(Eggs_larvae_list))

# Update only matching rows
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
  prey_species_latin_name %in% Inorganic_list, "digested",# "scale_shell_debris",
  prey_species_latin_name %in% Organic_list, "digested",
  # prey_species_latin_name %in% Unid_fish_list, "Teleostei",
  prey_species_latin_name %in% Unid_invert_list, "digested",
  (!is.na(life_stage)), "egg_larvae",
  (!is.na(life_stage) & is.na(prey_species_latin_name)), "digested",
  prey_species_latin_name %in% Rest_list, "digested",# "scale_shell_debris",
  prey_species_latin_name %in% Parasite_list, "parasite",
  # prey_species_latin_name %in% Mixed_prey_list, "Mixed_Prey",
  default = "identified_prey_species"
)]


# just to be sure but it's an error if it's match something.
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

### i) where necessary, add aphia ID from Andrew code -----------------------------
all.sp[code == 290, aphia_id_origin := 126417] # Clupeidae/Osmeridae -> Teleostei 293496 or as François 126417 (Clupea harengus)
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

# The worrm:: function don't work with long vector
# we have to create batches
aphia_id_origins_batches <- split(aphia_id_origins, ceiling(seq_along(aphia_id_origins) / 50))


aphia_ids_verified_list <- lapply(aphia_id_origins_batches, function(batche) {
  tryCatch({
    Sys.sleep(0.2) # small pause for the server to think if we had too much data to send
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

#rank_levels <- tolower(ranknfile)

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

#species_prey_aphia <- species_prey_aphia[aphia_id_origin %in% aphia_id_origins,]
aphia_id_origins_rest <- setdiff(aphia_id_origins,  unique(class_tbl_prey$aphia_id_origin))
aphia_id_origins_rest
#species_prey_aphiax <- species_prey_aphia[, .(.N), by = "aphia_id_origin"]





## b) Batch retrieve classification latin names, WoRMS --------------------------
## find those that have'nt aphia_id but latin_names -------------------------
latin_names  <- unique(species_prey[is.na(aphia_id_origin) &
                                                !is.na(prey_species_latin_name) &
                                                !(prey_species_latin_name %in% c("Obsolete", "Parasites,round worms")) &
                                                prey_category %in% c("identified_prey_species", "plantae_chromista", "parasite", "egg_larvae"),
                                              prey_species_latin_name])


# The worrm:: function don't work with long vector
# we have to create batches
latin_names_batches <- split(latin_names, ceiling(seq_along(latin_names) / 50))

latin_names_verified_list <- lapply(latin_names_batches, function(batche) {
  tryCatch({
    Sys.sleep(0.2) # small pause for the server to think if we had too much data to send
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

#rank_levels <- tolower(ranknfile)

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


### ii) merge WoRMS aphia_id and hierarchy -------------------------------------------
# A CORRIGER CHEZ LYSANDRE
species_prey_latin <- merge(
  x = subset(species_prey, is.na(aphia_id_origin)),
  y = class_tbl_prey2,
  by = "prey_species_latin_name",
  all.y = T
)


latin_names_rest <- setdiff(latin_names,  unique(class_tbl_prey2$prey_species_latin_name))
latin_names_rest


# 4) Merge identified preys together "Id", and find the remain ones "NoId" -----
species_prey_id <- rbind(species_prey_aphia, species_prey_latin, fill = T)
#species_prey_id <- unique(species_prey_id)

species_prey_noid_list <- setdiff(unique(species_prey$prey_species_latin_name), unique(species_prey_id$prey_species_latin_name))
species_prey_noid <- species_prey[species_prey$prey_species_latin_name %in% species_prey_noid_list,]
#species_prey_noid <- unique(species_prey_noid)

length(unique(species_prey_id$verified_name))
length(unique(species_prey_noid$prey_category))


# 5) Final export---------------------------------------------------------------
#NoIdentify_tobe_classify <- copy(species_prey_noid)
#Identify_tobe_group <- copy(species_prey_id)

save(species_prey_id, species_prey_noid, file = paste0(project_path, "/data/prey_taxon_data.rda"))

