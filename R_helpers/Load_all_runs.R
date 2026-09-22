# =============================================================================
# R_helpers/Load_all_runs.R
# -----------------------------------------------------------------------------
# Assemble un data.frame all_runs a partir des sorties reelles de 6a.
#
# OU SONT LES RESULTATS (verifie contre run_all_scenarios() et run_save_plot())
#   1. data/Sensitivity/all_runs_<level><family>.rda
#      -> objets res_biomass_<sc> et res_occ_<sc>, sc dans P1 P2 P3 P4 PT,
#         chacun = list(results, inventory, n_comp_failures)
#   2. Sensitivity_<level><family>/<mode>_<sc><family>_<level>_<p1>_vs_<p2>.rda
#      -> un objet `res`, meme structure (repli si le .rda groupe manque)
#
# Le contraste est derive des colonnes period_1/period_2 via
# contrast_from_periods() de Config_Mappings.R, jamais du nom de fichier.
# La famille vient de family_of(diagnostic) : le moteur a deja classe la cellule.
# =============================================================================

suppressPackageStartupMessages({library(dplyr); library(purrr); library(tibble)})

`%||%` <- function(a, b) if (is.null(a)) b else a

for (.p in c("R_helpers/Config_Mappings.R", "Config_Mappings.R")) {
  if (file.exists(.p)) { source(.p); break }
}
if (!exists("contrast_from_periods"))
  stop("Config_Mappings.R introuvable : place-le dans R_helpers/.")

.pooled_path <- function(level, fam, tag = "", base = "Sensitivity") {
  file.path("data", base, sprintf("all_runs_%s%s%s.rda", level, fam, tag))
}

.extract <- function(o) {
  if (is.data.frame(o)) return(o)
  if (is.list(o) && !is.null(o$results)) return(o$results)
  NULL
}

#' @param level "all_gulf" | "ecoregion" | "stratum"
#' @param prey_family "_1" | "_2" | "_PP"
#' @param tag suffixe de run : "" pour le run canonique, "_T" pour un run d'essai
#' @param base prefixe des dossiers : "Sensitivity" (canonique) ou "SensitivityT" (essai)
load_all_runs <- function(level = "all_gulf", prey_family = "_1", tag = "",
                          base = "Sensitivity",
                          dir_override = NULL, verbose = TRUE) {

  pooled <- .pooled_path(level, prey_family, tag, base)
  res <- NULL

  if (file.exists(pooled)) {
    e <- new.env(); nm <- load(pooled, envir = e)
    res <- bind_rows(map(nm, ~ .extract(get(.x, envir = e))))
    if (verbose) message("Charge depuis ", pooled, " (", length(nm), " objets)")
  } else {
    sens_dir <- dir_override %||% sprintf("%s_%s%s%s", base, level, prey_family, tag)
    if (!dir.exists(sens_dir)) {
      dirs <- grep("^Sensitivity", list.dirs(".", recursive = FALSE,
                                             full.names = FALSE), value = TRUE)
      stop("Ni ", pooled, " ni le dossier ", sens_dir, ".",
           "\nDossiers Sensitivity presents : ",
           if (length(dirs)) paste(dirs, collapse = ", ") else "aucun")
    }
    f <- list.files(sens_dir, pattern = "\\.rda$", full.names = TRUE)
    if (!length(f)) stop("Aucun .rda dans ", sens_dir)
    if (verbose) message("Assemblage depuis ", sens_dir, " (", length(f), " fichiers)")
    res <- bind_rows(map(f, function(ff) {
      e <- new.env(); nm <- load(ff, envir = e)
      bind_rows(map(nm, ~ .extract(get(.x, envir = e))))
    }))
  }

  if (is.null(res) || !nrow(res)) stop("Aucune ligne de resultats recuperee.")

  # --- contraste, depuis les periodes ----------------------------------------
  res$contrast <- contrast_from_periods(res$period_1, res$period_2)
  if (anyNA(res$contrast)) {
    bad <- res %>% filter(is.na(contrast)) %>% distinct(period_1, period_2)
    warning("Paires de periodes non reconnues : ",
            paste(sprintf("%s|%s", bad$period_1, bad$period_2), collapse = ", "))
  }
  res$contrast_type <- ifelse(as.character(res$contrast) %in% INTRA_CONTRASTS,
                              "within", "between")
  res$contrast_label <- contrast_label(res$contrast)

  # --- famille, depuis le diagnostic ecrit par le moteur ---------------------
  res$family <- as.character(family_of(res$diagnostic))

  unit <- intersect(c("Area", "str"), names(res))[1]
  res$cell <- paste(res$species, res$size_class, res[[unit]], sep = "|")
  attr(res, "unit_col") <- unit

  if (verbose) {
    message(sprintf("%d lignes | %d cellules | %d seuils | contrastes: %s",
                    nrow(res), dplyr::n_distinct(res$cell),
                    dplyr::n_distinct(res$x_threshold),
                    paste(levels(droplevels(res$contrast)), collapse = ", ")))
  }
  res
}

# Que contient reellement le projet ?
inventory_runs <- function() {
  pooled <- unlist(lapply(c("Sensitivity", "SensitivityT"), function(b) {
    f <- list.files(file.path("data", b), pattern = "^all_runs_.*\\.rda$")
    if (length(f)) file.path("data", b, f) else character(0)
  }))
  dirs <- grep("^Sensitivity_", list.dirs(".", recursive = FALSE,
                                          full.names = FALSE), value = TRUE)
  dirs <- dirs[!grepl("^Sensitivity_Plot", dirs)]
  list(
    pooled   = tibble(fichier = pooled),
    dossiers = map_dfr(dirs, ~ tibble(dossier = .x,
                                      n_rda = length(list.files(.x, pattern = "\\.rda$"))))
  )
}
