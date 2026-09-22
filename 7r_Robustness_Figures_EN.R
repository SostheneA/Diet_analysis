# =============================================================================
# 7r_Robustness_Figures_EN.R
# -----------------------------------------------------------------------------
# Rebuilds the robustness figures in English, publication-ready, from the CSV
# outputs of 6h, 6j, 6k and 6l. Writes them to Output_Appendices/ under the
# FigS<n>_ names that Paper_Figures.qmd picks up.
#
#   FigS8_null_calibration.png           6k  (NullCalibration_all_gulf_1/)
#   FigS9_calibrated_thresholds.png      6l  (Calibrated_all_gulf_1/)
#   FigS10_balanced_effort_gulf.png      6j  (Robustness_all_gulf_1_BAL/)
#   FigS11_balanced_effort_ecoregion.png 6j  (Robustness_ecoregion_1_BAL/)
#
# Nothing is recomputed: the script only reads and redraws.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

OUT <- "Output_Appendices"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Shared conventions -----------------------------------------------------

FAMILIES <- c("Stability", "Substitution", "Functional change", "Reorganization")

CONTRAST_LAB <- c(
  P1a = "2004 vs 2006",
  P1b = "2004\u201305 vs 2006",
  P2  = "2018 vs 2019",
  PTb = "2006 vs 2018",
  PTa = "2004\u201306 vs\n2018\u201319"
)

MODE_LAB <- c(biomass = "Biomass", occurrence = "Occurrence")
TYPE_COL <- c(within = "#4C8CB5", between = "#C8553D")

theme_paper <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1))
}

tidy_labels <- function(d) {
  d %>%
    mutate(
      family   = recode(family, Reorganisation = "Reorganization"),
      family   = factor(family, levels = FAMILIES),
      contrast = factor(contrast, levels = names(CONTRAST_LAB),
                        labels = unname(CONTRAST_LAB)),
      mode     = factor(mode, levels = names(MODE_LAB), labels = unname(MODE_LAB)))
}

read_req <- function(path) {
  if (!file.exists(path)) stop("Fichier introuvable : ", path, call. = FALSE)
  read_csv(path, show_col_types = FALSE)
}

save_fig <- function(p, name, w = 6.5, h = 4.5) {
  f <- file.path(OUT, name)
  ggsave(f, p, width = w, height = h, dpi = 300)
  message("Ecrit : ", f)
}

# =============================================================================
# FIGURE S8 — Null calibration of the three tests
# =============================================================================

fpr <- read_req("NullCalibration_all_gulf_1/null_fpr_by_sample_size.csv")

# The CSV carries French column names from 6k. Rename to ASCII keys: non-ASCII
# column names break across encodings, notably under Windows. Typographic
# labels are applied only at display time, through factor labels.
names(fpr)[grepl("^Composition", names(fpr))] <- "comp"
names(fpr)[grepl("^Diversit",    names(fpr))] <- "div"
names(fpr)[grepl("^Amplitude",   names(fpr))] <- "bs"

TEST_LAB <- c(comp = "Composition (PERMANOVA)",
              div  = "Diversity (H\u2032)",
              bs   = "Niche breadth (Bs)")

s8 <- fpr %>%
  pivot_longer(c(comp, div, bs), names_to = "test", values_to = "rate") %>%
  mutate(test  = factor(test, levels = names(TEST_LAB), labels = unname(TEST_LAB)),
         mode  = factor(mode, levels = names(MODE_LAB), labels = unname(MODE_LAB)),
         n_bin = factor(n_bin, levels = unique(n_bin)))

p8 <- ggplot(s8, aes(n_bin, rate, colour = test, group = test)) +
  geom_hline(yintercept = 5, linetype = 2, colour = "grey45") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  facet_wrap(~ mode) +
  scale_colour_manual(values = c("#2E7D32", "#1F5FA8", "#C8553D")) +
  labs(x = "Trawl sets per period (smaller of the two)",
       y = "Rejection rate under the null (%)") +
  theme_paper() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
save_fig(p8, "FigS8_null_calibration.png", h = 3.6)

# =============================================================================
# FIGURE S9 — Family frequencies under nominal and calibrated thresholds
# =============================================================================

cal <- read_req("Calibrated_all_gulf_1/calibrated_vs_nominal.csv") %>%
  tidy_labels() %>%
  pivot_longer(c(pct_nominal, pct_calibrated), names_to = "threshold",
               values_to = "pct") %>%
  mutate(threshold = factor(threshold, levels = c("pct_nominal", "pct_calibrated"),
                            labels = c("\u03b1 = 0.05", "Calibrated \u03b1")))

p9 <- ggplot(cal, aes(contrast, pct, fill = threshold)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  facet_grid(mode ~ family) +
  scale_fill_manual(values = c("grey70", "#1F5FA8")) +
  labs(x = NULL, y = "Classifiable cells (%)") +
  theme_paper()
save_fig(p9, "FigS9_calibrated_thresholds.png", w = 7.5, h = 5)

# =============================================================================
# FIGURES S10 / S11 — Balanced sampling effort
# =============================================================================

balanced_plot <- function(dir) {
  rep_d <- read_req(file.path(dir, "balanced_family_per_replicate.csv")) %>% tidy_labels()
  unb   <- read_req(file.path(dir, "balanced_vs_unbalanced.csv")) %>%
    select(contrast, contrast_type, mode, family, pct_unbalanced) %>%
    tidy_labels()

  ggplot(rep_d, aes(contrast, pct, fill = contrast_type)) +
    geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.3) +
    geom_point(data = unb, aes(contrast, pct_unbalanced), inherit.aes = FALSE,
               shape = 4, size = 2.2, stroke = 0.9) +
    facet_grid(mode ~ family) +
    scale_fill_manual(values = TYPE_COL, breaks = c("within", "between"),
                      labels = c(within = "Within-period contrast",
                                 between = "Between-period contrast")) +
    labs(x = NULL, y = "Classifiable cells (%)") +
    theme_paper()
}

if (dir.exists("Robustness_all_gulf_1_BAL")) {
  save_fig(balanced_plot("Robustness_all_gulf_1_BAL"),
           "FigS10_balanced_effort_gulf.png", w = 7.5, h = 5)
} else message("Robustness_all_gulf_1_BAL absent : FigS10 non produite.")

if (dir.exists("Robustness_ecoregion_1_BAL")) {
  save_fig(balanced_plot("Robustness_ecoregion_1_BAL"),
           "FigS11_balanced_effort_ecoregion.png", w = 7.5, h = 5)
} else message("Robustness_ecoregion_1_BAL absent : FigS11 non produite.")

message("7r termine.")
