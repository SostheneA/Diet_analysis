# =============================================================================
# R_helpers/paper_helpers.R
# -----------------------------------------------------------------------------
# Shared by Paper_Tables.qmd and Paper_Appendix.qmd: readers, table display,
# labels, and a finder for the tables written by the main pipeline (6f, 6g, 7,
# 8, 9), whose exact file names may vary.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
})

USE_FT  <- requireNamespace("flextable", quietly = TRUE)
OUTDIRS <- c("Output_Appendix_rebuilt", "Output_Tables", "Output_Appendices",
             "Output_Figures", "Output_Figures_extra")

FAMILIES     <- c("Stability", "Substitution", "Functional change", "Reorganization")
CONTRAST_LAB <- c(P1a = "2004 vs 2006", P1b = "2004–2005 vs 2006",
                  P2 = "2018 vs 2019", PTb = "2006 vs 2018",
                  PTa = "2004–2006 vs 2018–2019")
MODE_LAB     <- c(biomass = "Biomass", occurrence = "Occurrence")
LEVEL_LAB    <- c(all_gulf = "Gulf-wide", ecoregion = "Ecoregion", stratum = "Stratum")

fmt <- function(x, d = 1) formatC(x, format = "f", digits = d)

# Writes a table straight into the document (works for several per chunk).
# The table spans the page width and wraps text inside cells.
show_tbl <- function(d, digits = 1) {
  d <- d %>% mutate(across(where(is.double), ~ round(.x, digits)))
  if (USE_FT) {
    ft <- flextable::flextable(d)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
    ft <- flextable::set_table_properties(ft, layout = "autofit", width = 1)
    flextable::flextable_to_rmd(ft)
  } else {
    cat(knitr::kable(d), sep = "\n")
  }
  cat("\n\n")
  invisible(NULL)
}

rd <- function(path) if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL

absent <- function(what) {
  cat("\n*[Missing: ", what, " — run the corresponding analysis.]*\n\n", sep = "")
}

relabel <- function(d) {
  if ("family" %in% names(d))
    d <- d %>% mutate(family = recode(family, Reorganisation = "Reorganization"),
                      family = factor(family, levels = FAMILIES))
  if ("contrast" %in% names(d))
    d <- d %>% mutate(contrast = factor(contrast, levels = names(CONTRAST_LAB),
                                        labels = unname(CONTRAST_LAB)))
  if ("mode" %in% names(d))
    d <- d %>% mutate(mode = factor(mode, levels = names(MODE_LAB),
                                    labels = unname(MODE_LAB)))
  d
}

# ---- Tables written by the main pipeline -----------------------------------
# Looks for a file named after the table identifier, e.g. "TableA2_...csv",
# "Table_A2.xlsx" or "A2_effect_sizes.csv", in the output folders.
find_pipeline <- function(id) {
  dirs <- OUTDIRS[dir.exists(OUTDIRS)]
  if (!length(dirs)) return(character(0))
  # folders are searched in order; the first folder holding a match wins, so
  # the rebuilt appendix tables take precedence over older pipeline copies
  pat <- sprintf("(?i)(^|[^A-Za-z0-9])(table|tab)?[_ .-]?%s([^0-9A-Za-z]|_|$).*\\.(csv|xlsx|rds)$",
                 id)
  for (d in dirs) {
    f <- list.files(d, full.names = TRUE, recursive = TRUE)
    f <- f[grepl(pat, basename(f), perl = TRUE)]
    if (length(f)) return(f)
  }
  character(0)
}

read_any <- function(f) {
  ext <- tolower(tools::file_ext(f))
  if (ext == "csv")  return(read_csv(f, show_col_types = FALSE))
  if (ext == "rds")  { o <- readRDS(f); return(if (is.data.frame(o)) o else NULL) }
  if (ext == "xlsx" && requireNamespace("readxl", quietly = TRUE))
    return(readxl::read_excel(f))
  NULL
}

show_pipeline <- function(id, digits = 1) {
  hits <- find_pipeline(id)
  if (!length(hits)) {
    absent(sprintf("pipeline table %s (no file matching '%s' in %s)",
                   id, id, paste(OUTDIRS, collapse = ", ")))
    return(invisible(NULL))
  }
  if (length(hits) > 1)
    warning("Several files for ", id, ": ", paste(basename(hits), collapse = ", "),
            " -- first used", call. = FALSE)
  d <- read_any(hits[1])
  if (is.null(d)) absent(hits[1]) else show_tbl(relabel(d), digits)
}

# ---- Figures ---------------------------------------------------------------
fig <- function(pattern) {
  dirs <- OUTDIRS[dir.exists(OUTDIRS)]
  hits <- list.files(dirs, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) stop("introuvable : ", pattern, call. = FALSE)
  if (length(hits) > 1)
    warning("plusieurs fichiers pour '", pattern, "' : ",
            paste(basename(hits), collapse = ", "), " -- premier retenu", call. = FALSE)
  knitr::include_graphics(normalizePath(hits[1]))
}

# ---- Family frequencies, the primary summary -------------------------------
# Mean (SD across the 100 resolutions) per contrast and currency, one column
# per family, read from the M0 output of 6h.
freq_wide <- function(level, tag = "_T") {
  d <- rd(sprintf("Robustness_%s_1%s/M0_family_frequencies_reference.csv", level, tag))
  if (is.null(d)) return(NULL)
  d %>% relabel() %>%
    mutate(cell = sprintf("%s (%s)", fmt(pct_mean), fmt(sd_across_x))) %>%
    select(Contrast = contrast, Currency = mode, family, cell) %>%
    pivot_wider(names_from = family, values_from = cell) %>%
    select(Contrast, Currency, any_of(FAMILIES)) %>%
    arrange(Currency, Contrast)
}

# In the appendix a missing figure shows a visible placeholder rather than
# stopping the whole render, since pipeline outputs may be produced separately.
fig_safe <- function(pattern) {
  dirs <- OUTDIRS[dir.exists(OUTDIRS)]
  hits <- character(0)
  for (d in dirs) {
    hits <- list.files(d, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) break
  }
  if (!length(hits)) { absent(sprintf("figure matching '%s'", pattern)); return(invisible(NULL)) }
  knitr::include_graphics(normalizePath(hits[1]))
}
