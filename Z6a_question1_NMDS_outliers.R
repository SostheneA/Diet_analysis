# =============================================================================
# 6f — NMDS INTERACTIVES POUR REPÉRAGE VISUEL DES OUTLIERS (par prédateur)
# =============================================================================
# Objectif : pour CHAQUE prédateur x classe de taille, produire deux NMDS
# (occurrence %FO et biomasse %W) en plotly, afin d'identifier visuellement les
# trawl sets atypiques (survol -> set_uid). Une fois les outliers listés, la
# fonction rebuild_aligned_cell() reconstruit les matrices alignées filtrées.
#
# Choix de cadrage (validés avec l'analyste) :
#   * périodes : period_map_PT   -> P1 = 2004-2006, P2 = 2018-2019
#   * aires    : POOLÉES          -> ar = NULL
#   * tailles  : SÉPARÉES         -> une figure par classe de taille
#   * proies   : prey_category_100_2
#
# Réutilise telles quelles les fonctions du moteur (6a_Engine_Trophic_ll…) :
#   make_dat_classed(), build_biomass_set(), build_occ_set(), align_mats().
# Ces fonctions doivent être en mémoire (source du moteur) OU le sont déjà.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(vegan); library(plotly); library(htmlwidgets)
})

stopifnot(exists("dat_classed"),
          all(sapply(c("make_dat_classed", "build_biomass_set",
                       "build_occ_set", "align_mats"), exists)))

# --- Paramètres ---------------------------------------------------------------
if (!exists("OCC_BINARY")) OCC_BINARY <- FALSE
COL_PREY  <- "prey_category_500_2"
PREY_LAB  <- sub("^prey_category_", "", COL_PREY)   # ex. "100_2" -> affiché dans les titres
DATE_TAG  <- format(Sys.Date(), "%Y%m%d")           # ex. "20260810" -> nom de fichier
OUT_DIR   <- file.path(here::here(), "nmds_out")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

period_map_PT <- list(
  "2004-2006" = c(2004, 2005, 2006),
  "2018-2019" = c(2018, 2019)
)
PER_LABELS <- names(period_map_PT)   # c("2004-2006", "2018-2019")

# Jeu de données au niveau période (aires conservées dans la colonne Area,
# mais on ne filtre PAS dessus -> ar = NULL dans les build_*).
dat_PT <- make_dat_classed(diet_clean = dat_classed, period_map = period_map_PT)

# =============================================================================
# Cœur : construire une NMDS plotly pour une matrice alignée
# =============================================================================
# Distance identique à run_composition_set() du moteur : Bray-Curtis sur diète
# relative (decostand "total"), ou Jaccard binaire si OCC_BINARY = TRUE en occ.
nmds_plotly_from_aligned <- function(m_al, mode = c("biomass", "occurrence"),
                                     title = "", seed = 6820,
                                     n_sto = c(NA_integer_, NA_integer_)) {
  mode <- match.arg(mode)
  m1 <- m_al$m1; m2 <- m_al$m2
  if (is.null(m1) || is.null(m2)) return(NULL)
  if (nrow(m1) < 2 || nrow(m2) < 2) return(NULL)       # NMDS ininterprétable
  if (nrow(m1) + nrow(m2) < 4)      return(NULL)

  combined <- rbind(m1, m2)
  groups <- factor(c(rep(PER_LABELS[1], nrow(m1)),
                     rep(PER_LABELS[2], nrow(m2))), levels = PER_LABELS)

  use_jaccard <- (mode == "occurrence") && isTRUE(OCC_BINARY)
  d <- if (use_jaccard) {
    vegdist(combined, method = "jaccard", binary = TRUE)
  } else {
    vegdist(decostand(combined, "total"), method = "bray")
  }

  set.seed(seed)
  nmds <- tryCatch(metaMDS(d, k = 2, trymax = 100, trace = FALSE),
                   error = function(e) NULL)
  if (is.null(nmds)) return(NULL)

  scores_df <- vegan::scores(nmds, display = "sites") |>
    as_tibble(rownames = "set_uid") |>
    mutate(period = groups)

  # ggplotly() IGNORE le subtitle d'un ggplot : on met donc toutes les infos
  # (stress, n sets, n estomacs) dans le TITRE plotly via un <sub> multi-lignes.
  subtitle_txt <- sprintf(
    "Stress = %.3f  |  n sets: %s = %d, %s = %d  |  n estomacs: %s = %s, %s = %s",
    nmds$stress, PER_LABELS[1], nrow(m1), PER_LABELS[2], nrow(m2),
    PER_LABELS[1], n_sto[1], PER_LABELS[2], n_sto[2])

  g <- ggplot(scores_df, aes(NMDS1, NMDS2, colour = period)) +
    geom_point(aes(text = paste0("set_uid: ", set_uid,
                                 "<br>NMDS1: ", round(NMDS1, 3),
                                 "<br>NMDS2: ", round(NMDS2, 3))),
               alpha = 0.75, size = 2) +
    stat_ellipse(aes(fill = period), geom = "polygon", alpha = 0.12, type = "t") +
    labs(colour = "Période", fill = "Période")

  wid <- ggplotly(g, tooltip = "text") |>
    layout(title = list(
             text = paste0(title, "<br><sub>", subtitle_txt, "</sub>"),
             x = 0, xanchor = "left"),
           margin = list(t = 90))

  list(widget = wid,
       stress = nmds$stress, n1 = nrow(m1), n2 = nrow(m2),
       scores = scores_df)
}

# Construit les matrices occ/bio alignées d'une cellule prédateur x taille
build_cell_aligned <- function(sp, sz, data = dat_PT) {
  b1 <- build_biomass_set(data, sp = sp, sz = sz, ar = NULL,
                          per = PER_LABELS[1], col_prey = COL_PREY)
  b2 <- build_biomass_set(data, sp = sp, sz = sz, ar = NULL,
                          per = PER_LABELS[2], col_prey = COL_PREY)
  o1 <- build_occ_set(data, sp = sp, sz = sz, ar = NULL,
                      per = PER_LABELS[1], col_prey = COL_PREY)
  o2 <- build_occ_set(data, sp = sp, sz = sz, ar = NULL,
                      per = PER_LABELS[2], col_prey = COL_PREY)
  list(bio = align_mats(b1, b2), occ = align_mats(o1, o2),
       raw = list(b1 = b1, b2 = b2, o1 = o1, o2 = o2))
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# Nombre d'estomacs distincts (effort) par période pour la cellule sp x taille,
# aires poolées. Indépendant du mode occ/bio (propriété du prédateur).
count_stomachs <- function(sp, sz, data = dat_PT) {
  tab <- data |>
    filter(predator_species_common_name == sp,
           as.character(size_class) == sz, !is.na(stomach_id)) |>
    distinct(period, stomach_id) |>
    count(period)
  n <- tab$n[match(PER_LABELS, tab$period)]
  ifelse(is.na(n), 0L, n)
}

# =============================================================================
# Boucle : 40 prédateurs x classes de taille -> HTML dans nmds_out/
# =============================================================================
predators <- sort(unique(dat_PT$predator_species_common_name))
sizes      <- c("adult", "juvenile", "one_class")

log_rows <- list()
for (sp in predators) {
  for (sz in sizes) {
    cell <- build_cell_aligned(sp, sz)
    n_sto <- count_stomachs(sp, sz)
    for (mode in c("occurrence", "biomass")) {
      m_al <- if (mode == "occurrence") cell$occ else cell$bio
      tag  <- if (mode == "occurrence") "occ" else "bio"
      ttl  <- sprintf("%s — %s — %s — %s (aires poolées)", sp, sz,
                      if (mode == "occurrence") "%FO occurrence" else "%W biomasse",
                      PREY_LAB)
      res <- nmds_plotly_from_aligned(m_al, mode = mode, title = ttl, n_sto = n_sto)
      if (is.null(res)) next
      fname <- sprintf("%s__%s__%s__%s__%s.html",
                       safe_name(sp), sz, tag, safe_name(PREY_LAB), DATE_TAG)
      # Bibliothèque JS partagée (nmds_lib_assets/) -> fichiers HTML légers.
      saveWidget(res$widget, file.path(OUT_DIR, fname),
                 selfcontained = FALSE, libdir = "nmds_lib_assets", title = ttl)
      log_rows[[length(log_rows) + 1]] <- tibble(
        predator = sp, size = sz, mode = tag,
        n_P1 = res$n1, n_P2 = res$n2,
        stress = round(res$stress, 3), file = fname)
    }
  }
}

nmds_index <- bind_rows(log_rows) |> arrange(predator, size, mode)
write_csv(nmds_index, file.path(OUT_DIR, sprintf("nmds_index__%s__%s.csv", safe_name(PREY_LAB), DATE_TAG)))
message("Figures écrites : ", nrow(nmds_index), " -> ", OUT_DIR)




# =============================================================================
# Reconstruction après repérage visuel des outliers
# =============================================================================
# Donner une liste nommée : outliers$<safe_name(sp)>$<size> = c("2018_T_103", …)
# rebuild_aligned_cell() retire ces set_uid des matrices brutes puis ré-aligne.
rebuild_aligned_cell <- function(sp, sz, drop_ids = character(0), data = dat_PT) {
  cell <- build_cell_aligned(sp, sz, data)
  drop_rows <- function(m) if (is.null(m)) m else
    m[!rownames(m) %in% drop_ids, , drop = FALSE]
  list(
    bio = align_mats(drop_rows(cell$raw$b1), drop_rows(cell$raw$b2)),
    occ = align_mats(drop_rows(cell$raw$o1), drop_rows(cell$raw$o2))
  )
}

# 1. Après avoir survolé la NMDS de la morue adulte et repéré des sets atypiques,
#    tu listes leurs set_uid (ceux lus dans l'infobulle plotly) :
cell_cod <- rebuild_aligned_cell(
  sp       = "Atlantic cod",
  sz       = "adult",
  drop_ids = c("2018_T_103", "2019_T_11")   # <- outliers repérés visuellement
)

# 2. Tu récupères les matrices alignées SANS ces sets, prêtes pour la suite
#    (PERMANOVA, jackknife, etc.) — une version biomasse et une occurrence :
cell_cod$bio$m1   # %W, période 2004-2006, outliers retirés
cell_cod$bio$m2   # %W, période 2018-2019
cell_cod$occ$m1   # %FO, période 2004-2006
cell_cod$occ$m2   # %FO, période 2018-2019

# 3. Vérifier que les outliers ont bien disparu :
c("2018_T_103", "2019_N_58") %in% rownames(cell_cod$bio$m2)   # attendu : FALSE FALSE
