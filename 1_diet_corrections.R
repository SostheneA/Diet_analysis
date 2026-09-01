#------------------------------------------------------------------------------#
# 1_diet_corrections.R - READ THE SURVEY DIET FILES AND CORRECT KNOWN ISSUES
#
# Input  : data/Survey_data/4T_*.csv (one file per vessel x year)
# Output : data/diet_corr.RData (diet_corr, one row per prey record, with
#          vessel.code, year and the is_empty flag per stomach)
#------------------------------------------------------------------------------#
rm(list = ls())

project_path <- here::here()

library(data.table)
library(dplyr)

# 1) Read ----------------------------------------------------------------------

survey_files <- list(
  list(path = "data/Survey_data/4T_Needler_2004.csv", vessel = "N"),
  list(path = "data/Survey_data/4T_Needler_2005.csv", vessel = "N"),
  list(path = "data/Survey_data/4T_TEL_2004.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2005.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2006.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2018.csv", vessel = "T"),
  list(path = "data/Survey_data/4T_TEL_2019.csv", vessel = "T")
)

diet_list <- lapply(survey_files, function(f) {
  dt <- fread(f$path, encoding = "UTF-8")
  dt[, vessel.code := f$vessel]
  return(dt)
})

diet_raw <- rbindlist(diet_list, fill = TRUE)
diet_raw <- diet_raw %>% mutate(year = cruise_year)

setDT(diet_raw)

diet_corr <- copy(diet_raw)


# 2) corrections ---------------------------------------------------------------
#------------------------------------------------------------------------------#

## empty but not empty ---------------------------------------------------------
diet_corr$prey_species_latin_name <- ifelse((is.na(diet_corr$prey_species_latin_name) | diet_corr$prey_species_latin_name == "") &
                                             diet_corr$prey_species_common_name  %in%
                                             c("Fish remains", "Eggs - Atlantic herring", "Mud"),
                                           diet_corr$prey_species_common_name,
                                           diet_corr$prey_species_latin_name)

## Amphipod but just state in comments or or in WoRMs but not retrieved---------
diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_id == 60840, "Amphipoda", diet_corr$prey_species_latin_name)
diet_corr$prey_species_common_name <- ifelse(diet_corr$prey_id == 60840, "Amphipod unidentified", diet_corr$prey_species_common_name)
diet_corr$prey_species_code <- ifelse(diet_corr$prey_id == 60840, 2800, diet_corr$prey_species_code)


## 68342 10044 Hyalellidae Hyalellidae (f.) ----
diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_id == 68342, "Hyaloidea", diet_corr$prey_species_latin_name)
diet_corr$prey_species_common_name <- ifelse(diet_corr$prey_id == 68342, "Hyaloidea", diet_corr$prey_species_common_name)
diet_corr$prey_species_code <- ifelse(diet_corr$prey_id == 68342, NA_integer_, diet_corr$prey_species_code)

## white baracudina (one species but two different names) ----------------------
diet_corr$predator_species_latin_name <- ifelse(diet_corr$predator_species_latin_name == "Notolepis rissoi",
                                               "Arctozenus risso",
                                               diet_corr$predator_species_latin_name)

diet_corr$predator_species_code <- ifelse(diet_corr$predator_species_code == 727,
                                         712,
                                         diet_corr$predator_species_code)


## code 910 (Odobenidae, walruses) is not a prey code ---------------------------
diet_corr$prey_species_code <- ifelse(diet_corr$prey_species_code == 910,
                                     NA_integer_,
                                     diet_corr$prey_species_code)


## code corrections ------------------------------------------------------------
diet_corr$prey_species_code <- ifelse(diet_corr$prey_species_code == 2500,
                                     2507,
                                     diet_corr$prey_species_code) # Crabs (could be Anomura or Brachyura) -> Pleocyemata


## latin name corrections not in the Gulf gulf::species.names table) -----------
taxa_not_in_all_sp <- c("Nemertea (p.)", "Malacostraca (c.)", "Acteonidae (f.)",
                        "Eusiridae (f.)",
                        "Eunicida (o.)", "Liljeborgiidae (f.)", "Ampithoidae (f.)")

diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_species_latin_name %in% taxa_not_in_all_sp,
                                           gsub(" .*", "", trimws(diet_corr$prey_species_latin_name)),
                                           diet_corr$prey_species_latin_name)




#3) redefine emptyness ---------------------------------------------------------
#------------------------------------------------------------------------------#
diet_corr$prey_species_latin_name <- ifelse(diet_corr$prey_species_latin_name %in% c(""),
                                            NA_character_,
                                            diet_corr$prey_species_latin_name)

diet_corr$prey_species_common_name <- ifelse(diet_corr$prey_species_common_name %in% c(""),
                                            NA_character_,
                                            diet_corr$prey_species_common_name)

# could also be parasites but not intestinal: Caligus sp., Gnathia cerina

potential_empty_elements <- c(
  "Parasites,round worms",
  "Cestoda (c.)",
  "Trematoda",
  "Scales - unid.",
  "Scales - unid. fish",
  "Shells - unid. fish",
  "Mollusca sp. empty",
  "Casing - unidentified invertebrate",
  "Mucus",
  "Sand",
  "Stones and rocks",
  "Foreign articles,
                  garbage",
  "Mud",
  "Foreign articles,garbage"
)

emptyness <- diet_corr[,
                       .(
                         prey_list = list({
                           u <- unique(trimws(prey_species_latin_name))
                           u <- u[!is.na(u) & nzchar(u)]  # remove NA and ""
                           u
                         }),
                         codes = paste0(unique(na.omit(prey_species_code)), collapse = ","),
                         n_prey_names = uniqueN(prey_species_latin_name, na.rm = TRUE)
                       ),
                       by = .(
                         stomach_id
                       )
][
  ,
  keep := vapply(prey_list, function(x) length(x) == 0L || all(x %chin% potential_empty_elements), logical(1L))
][keep == T,

][
  ,
  `:=`(
    prey_names = vapply(prey_list, paste, character(1L), collapse = ";"),
    prey_list = NULL,
    keep = NULL
  )
]

empty_stomach_ids <- unique(emptyness$stomach_id)

diet_corr[, is_empty := fifelse(stomach_id %in% empty_stomach_ids, TRUE, FALSE)]
setDT(diet_corr)


# Stomachs with an empty-looking row that are not empty.
problematic_stomach <- diet_corr[!is_empty & is.na(prey_species_common_name),]$stomach_id
problematic_stomach_table1 <- diet_corr[stomach_id %in% c(problematic_stomach),]
problematic_stomach_table <- diet_corr[stomach_id %in% c(problematic_stomach),
                       .(
                         preys = paste0(unique(prey_species_code), collapse = ",")
                       ),
                       by = .(
                         stomach_id
                       )
]

diet_corr[stomach_id %in% problematic_stomach &
            is.na(prey_species_code) &
            is.na(prey_species_common_name) &
            is.na(prey_species_latin_name),

          `:=`(
            prey_species_latin_name  = "no prey information",
            prey_species_common_name = "no prey information",
            prey_species_code        = NA_integer_
          )
]
setDT(diet_corr)

# 4) Save ----------------------------------------------------------------------
#------------------------------------------------------------------------------#
save(diet_corr, file = paste0(project_path, "/data/diet_corr.RData"))

rm(emptyness, problematic_stomach_table, problematic_stomach_table1, taxa_not_in_all_sp)

