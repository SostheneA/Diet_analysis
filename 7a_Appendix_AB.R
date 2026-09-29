# =============================================================================
# 7a_Appendix_AB.R
# -----------------------------------------------------------------------------
# Rebuilds the Appendix A and B tables and figures from the engine runs, so the
# appendices no longer depend on file names produced elsewhere, and checks every
# number the manuscript cites against the recomputed values.
#
#   TableA2_effect_sizes.csv        compositional effect sizes (Gulf-wide)
#   TableA5_resolution_range.csv    range of family frequencies across resolutions
#   TableA6_sparse_units.csv        ecoregion / stratum after excluding sparse units
#   TableB1/B2/B3_grouping_rules.csv family frequencies under the three rules
#   FigB1/B2/B3_grouping_rules.png
#   VERIF_appendix_vs_manuscript.csv
#
# Output folder: Output_Appendix_rebuilt/ (read first by Paper_Appendix.qmd).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2)
})
source("R_helpers/Load_all_runs.R")

OUT <- "Output_Appendix_rebuilt"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FAM  <- c("Stability", "Substitution", "Functional change", "Reorganisation")
FAM_EN <- c(FAM[1:3], "Reorganization")
CLAB <- c(P1a = "2004 vs 2006", P1b = "2004\u20132005 vs 2006", P2 = "2018 vs 2019",
          PTb = "2006 vs 2018", PTa = "2004\u20132006 vs 2018\u20132019")
MLAB <- c(biomass = "Biomass", occurrence = "Occurrence")
LLAB <- c(all_gulf = "Gulf-wide", ecoregion = "Ecoregion", stratum = "Stratum")
lab <- function(d) {
  d %>% mutate(
    across(any_of("family"), ~ factor(recode(.x, Reorganisation = "Reorganization"), levels = FAM_EN)),
    across(any_of("contrast"), ~ factor(.x, levels = names(CLAB), labels = unname(CLAB))),
    across(any_of("mode"), ~ factor(.x, levels = names(MLAB), labels = unname(MLAB))))
}
f1 <- function(x, d = 1) formatC(x, format = "f", digits = d)
MIN_CELLS_UNIT <- 5        # sparse-unit rule of Section 2.5
SHARED_MAX     <- 250      # thresholds shared by the three grouping rules

# Primary rule: the patched run (_T) if present, otherwise the original.
load_rule <- function(level, fam) {
  for (tag in if (fam == "_1") c("_T", "") else "") {
    r <- try(load_all_runs(level, fam, tag = tag, verbose = FALSE), silent = TRUE)
    if (!inherits(r, "try-error")) {
      message(sprintf("%-9s %-3s <- tag '%s' (%d rows)", level, fam, tag, nrow(r)))
      return(r)
    }
  }
  NULL
}

# Family frequencies with the paper's aggregation: per spatial unit (absent
# families = 0), mean over units, then per resolution or mean over resolutions.
pct_unit <- function(d, unit) {
  d %>% filter(testable, !is.na(family)) %>%
    count(contrast, contrast_type, mode, x_threshold, !!rlang::sym(unit),
          family, name = "n") %>%
    group_by(contrast, contrast_type, mode, x_threshold, !!rlang::sym(unit)) %>%
    mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
    complete(nesting(contrast, contrast_type, mode, x_threshold, !!rlang::sym(unit)),
             family = FAM, fill = list(n = 0L, pct = 0)) %>%
    group_by(contrast, contrast_type, mode, x_threshold, family) %>%
    summarise(pct = mean(pct), .groups = "drop")
}
mean_over_x <- function(p) {
  p %>% group_by(contrast, contrast_type, mode, family) %>%
    summarise(pct = mean(pct), .groups = "drop")
}

verif <- list()
chk <- function(what, got, expected) {
  verif[[length(verif) + 1L]] <<- tibble(
    item = what, manuscript = expected, recomputed = round(got, 3),
    match = abs(round(got, 1) - expected) < 0.15)
}

# =============================================================================
# APPENDIX A
# =============================================================================
g1 <- load_rule("all_gulf", "_1")
if (is.null(g1)) stop("Runs Gulf-wide '_1' introuvables.")
UNIT_G <- attr(g1, "unit_col")

# ---- A2 : compositional effect sizes ---------------------------------------
a2 <- g1 %>% filter(testable) %>%
  group_by(contrast, contrast_type, mode) %>%
  summarise(BCmed = median(BC, na.rm = TRUE),
            R2med = median(R2_comp, na.rm = TRUE),
            nrows = n(), .groups = "drop") %>%
  arrange(mode, contrast)
write_csv(a2 %>% lab() %>% transmute(Contrast = contrast, Currency = mode,
             BCmed = round(BCmed, 3), R2med = round(R2med, 3), nrows = nrows) %>%
            arrange(Currency, Contrast) %>%
            rename(!!"Median Bray\u2013Curtis" := BCmed,
                   !!"Median PERMANOVA R\u00b2" := R2med,
                   !!"Cell \u00d7 resolution" := nrows),
          file.path(OUT, "TableA2_effect_sizes.csv"))

bc <- function(ct, md) a2$BCmed[a2$contrast == ct & a2$mode == md]
chk("A2 BC PTa biomass",    bc("PTa", "biomass"),    0.576)
chk("A2 BC PTa occurrence", bc("PTa", "occurrence"), 0.564)

# ---- A5 : range across the 100 resolutions ---------------------------------
pg <- pct_unit(g1, UNIT_G)
a5 <- pg %>% group_by(contrast, contrast_type, mode, family) %>%
  summarise(Mean = mean(pct), Minimum = min(pct), Maximum = max(pct),
            SD = sd(pct), .groups = "drop") %>%
  arrange(mode, family, contrast)
write_csv(a5 %>% lab() %>%
            mutate(cell = sprintf("%s [%s\u2013%s]", f1(Mean), f1(Minimum), f1(Maximum))) %>%
            select(Contrast = contrast, Currency = mode, family, cell) %>%
            pivot_wider(names_from = family, values_from = cell) %>%
            select(Contrast, Currency, any_of(FAM_EN)) %>% arrange(Currency, Contrast),
          file.path(OUT, "TableA5_resolution_range.csv"))

st <- function(ct, md, f) a5 %>% filter(contrast == ct, mode == md, family == "Stability") %>% pull(f)
chk("A5 Stability PTa biomass min",    st("PTa", "biomass", "Minimum"),    5.1)
chk("A5 Stability PTa biomass max",    st("PTa", "biomass", "Maximum"),    13.6)
chk("A5 Stability PTa occurrence min", st("PTa", "occurrence", "Minimum"), 3.4)
chk("A5 Stability PTa occurrence max", st("PTa", "occurrence", "Maximum"), 10.2)
lw <- a5 %>% filter(contrast_type == "within", family == "Stability") %>%
  group_by(mode) %>% summarise(m = min(Minimum))
chk("A5 lowest within Stability biomass",    lw$m[lw$mode == "biomass"],    43.6)
chk("A5 lowest within Stability occurrence", lw$m[lw$mode == "occurrence"], 44.4)

# ---- A6 : excluding sparse spatial units -----------------------------------
a6 <- map_dfr(c("ecoregion", "stratum"), function(l) {
  r <- load_rule(l, "_1"); if (is.null(r)) return(NULL)
  u <- attr(r, "unit_col")
  full <- mean_over_x(pct_unit(r, u)) %>% rename(allunits = pct)
  # a unit is kept at a given resolution when it holds >= 5 classifiable cells
  keep <- r %>% filter(testable, !is.na(family)) %>%
    count(contrast, mode, x_threshold, !!rlang::sym(u), name = "n_cells") %>%
    filter(n_cells >= MIN_CELLS_UNIT)
  rr <- r %>% semi_join(keep, by = c("contrast", "mode", "x_threshold", u))
  restr <- mean_over_x(pct_unit(rr, u)) %>% rename(restricted = pct)
  nu <- keep %>% count(contrast, mode, x_threshold, name = "n_units") %>%
    group_by(contrast, mode) %>%
    summarise(unitsret = round(mean(n_units), 1), .groups = "drop")
  full %>% left_join(restr, by = c("contrast", "contrast_type", "mode", "family")) %>%
    left_join(nu, by = c("contrast", "mode")) %>% mutate(Scale = l, .before = 1)
})
if (nrow(a6)) write_csv(a6 %>% lab() %>%
            mutate(Scale = LLAB[Scale],
                   cell = sprintf("%s \u2192 %s", f1(allunits), ifelse(is.na(restricted), "\u2014", f1(restricted))),
                   unitsret = coalesce(unitsret, 0)) %>%
            select(Scale, Contrast = contrast, Currency = mode,
                   unitsret, family, cell) %>%
            pivot_wider(names_from = family, values_from = cell) %>%
            select(Scale, Contrast, Currency, unitsret, any_of(FAM_EN)) %>%
            arrange(Scale, Currency, Contrast) %>%
            rename(!!"Units retained" := unitsret),
          file.path(OUT, "TableA6_sparse_units.csv"))

if (nrow(a6)) {
  s6 <- function(l, ct, md, col) a6 %>% filter(Scale == l, contrast == ct, mode == md,
                                               family == "Stability") %>% pull(col)
  chk("A6 ecoregion 2018v2019 Stability biomass, all",   s6("ecoregion", "P2", "biomass", "allunits"), 69.5)
  chk("A6 ecoregion 2018v2019 Stability biomass, >=5",   s6("ecoregion", "P2", "biomass", "restricted"), 59.4)
  chk("A6 ecoregion 2018v2019 Stability occurrence, >=5", s6("ecoregion", "P2", "occurrence", "restricted"), 71.1)
  chk("A6 stratum PTa Stability biomass, >=5",           s6("stratum", "PTa", "biomass", "restricted"), 35.8)
  chk("A6 stratum PTa Stability occurrence, >=5",        s6("stratum", "PTa", "occurrence", "restricted"), 24.0)
}

# =============================================================================
# APPENDIX B : the three grouping rules on the shared thresholds
# =============================================================================
RULES <- c(`_1` = "Pooled (primary)", `_2` = "Pooled cross-predator",
           `_PP` = "Per-predator")

b_level <- function(l) {
  runs <- imap(RULES, function(lab, fam) load_rule(l, fam)) %>% compact()
  if (length(runs) < 2) return(NULL)
  xs <- reduce(map(runs, ~ unique(.x$x_threshold)), intersect)
  xs <- xs[xs <= SHARED_MAX & xs %% 10 == 0]
  map_dfr(names(runs), function(fam) {
    r <- runs[[fam]] %>% filter(x_threshold %in% xs)
    mean_over_x(pct_unit(r, attr(runs[[fam]], "unit_col"))) %>%
      mutate(Rule = RULES[[fam]], n_thresholds = length(xs))
  }) %>% mutate(Scale = l, .before = 1)
}

B <- list()
for (i in seq_along(c("all_gulf", "ecoregion", "stratum"))) {
  l <- c("all_gulf", "ecoregion", "stratum")[i]
  b <- b_level(l)
  if (is.null(b)) { message("B", i, " : regles alternatives absentes pour ", l); next }
  B[[l]] <- b
  write_csv(b %>% lab() %>%
              mutate(Rule = factor(Rule, levels = RULES), pct = round(pct, 1)) %>%
              select(Rule, Contrast = contrast, Currency = mode, family, pct) %>%
              pivot_wider(names_from = family, values_from = pct) %>%
              select(Rule, Contrast, Currency, any_of(FAM_EN)) %>%
              arrange(Currency, Rule, Contrast),
            file.path(OUT, sprintf("TableB%d_grouping_rules.csv", i)))

  p <- b %>%
    mutate(family = factor(recode(family, Reorganisation = "Reorganization"),
                           levels = c(FAM[1:3], "Reorganization")),
           Rule = factor(Rule, levels = RULES),
           contrast = factor(contrast, levels = c("P1a", "P1b", "P2", "PTb", "PTa"),
                             labels = c("2004 vs 2006", "2004\u201305 vs 2006",
                                        "2018 vs 2019", "2006 vs 2018",
                                        "2004\u201306 vs\n2018\u201319")),
           mode = factor(mode, levels = c("biomass", "occurrence"),
                         labels = c("Biomass", "Occurrence"))) %>%
    ggplot(aes(contrast, pct, fill = Rule)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_grid(mode ~ family) +
    scale_fill_manual(values = c("#1F5FA8", "#7FA7D4", "#C8553D")) +
    labs(x = NULL, y = "Classifiable cells (%)", fill = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(OUT, sprintf("FigB%d_grouping_rules.png", i)), p,
         width = 7.5, height = 5, dpi = 300)
}

if (!is.null(B$all_gulf)) {
  gb <- function(rule, md, f) B$all_gulf %>% filter(Rule == rule, contrast == "PTa",
                                                    mode == md, family == f) %>% pull(pct)
  chk("B1 per-predator Stability biomass",       gb("Per-predator", "biomass", "Stability"), 31.7)
  chk("B1 per-predator Stability occurrence",    gb("Per-predator", "occurrence", "Stability"), 21.7)
  chk("B1 primary Stability biomass (shared x)", gb("Pooled (primary)", "biomass", "Stability"), 8.3)
  chk("B1 primary Stability occurrence",         gb("Pooled (primary)", "occurrence", "Stability"), 4.5)
  chk("B1 per-predator Substitution biomass",    gb("Per-predator", "biomass", "Substitution"), 39.0)
  chk("B1 per-predator Substitution occurrence", gb("Per-predator", "occurrence", "Substitution"), 37.9)
  chk("B1 primary Substitution biomass",         gb("Pooled (primary)", "biomass", "Substitution"), 48.5)
  chk("B1 primary Substitution occurrence",      gb("Pooled (primary)", "occurrence", "Substitution"), 44.3)
}

# =============================================================================
# VERIFICATION AGAINST THE MANUSCRIPT
# =============================================================================
v <- bind_rows(verif)
write_csv(v, file.path(OUT, "VERIF_appendix_vs_manuscript.csv"))
cat("\n================ CONTROLE CONTRE LE MANUSCRIT ================\n")
print(as.data.frame(v), row.names = FALSE)
cat(sprintf("\n%d / %d valeurs concordent a 0.1 pres.\n", sum(v$match), nrow(v)))
if (any(!v$match))
  cat("Les ecarts signalent une convention de calcul differente de celle qui a produit\n",
      "le chiffre du texte : a examiner avant de retenir l'une ou l'autre.\n")
