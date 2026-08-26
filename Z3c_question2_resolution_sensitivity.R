#------------------------------------------------------------------------------#
# 3c_question2_resolution_sensitivity.R
# POOLED vs PER-PREDATOR taxonomic grouping -- comparison in the style of 3b.
#
# Two grouping schemes, both built by 3_taxonomic_groups_ll20260812.R and both
# using the SAME _1 rule (n_predator = 1, no species->genus collapse) :
#   * prey_groups             -> POOLED  : thresholds evaluated on stomach counts
#                                summed over ALL predators (prey_category_X_1,
#                                keyed by prey_species_common_name).
#   * prey_groups_by_predator -> PER PRED: the _1 sweep rerun SEPARATELY inside
#                                each predator's own stomachs (prey_category_X_1,
#                                keyed by predator x prey taxon).
#
# We attach BOTH labelings to the SAME dat_classed diet records (pooled labels
# joined from prey_groups by common name, per-predator labels joined from
# prey_groups_by_predator by predator x common name). The only difference between
# schemes is therefore "pooled vs per-predator"; outside identified prey in
# phyla_list the two schemes coincide.
#
# Because a category needs X stomachs WITHIN a single predator to survive the
# per-predator sweep (vs X across all predators when pooled), the per-predator
# scheme is generally COARSER for an individual predator -- the analogue of
# _1 (fine) -> _2 (coarse) in 3b.
#
# Outputs (as in 3b) :
#   A) resolution_cmp        - # prey categories per predator: pooled vs per-pred
#   B) idx / idx_summ        - n_cat, H, Bs (pooled diet) + BC, R2 (P1 vs P2),
#                              per predator x threshold x scheme x mode
#   C) loss / loss_summ      - dominant-prey reassignment pooled -> per-predator
#   + figures for each
#------------------------------------------------------------------------------#

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(vegan)
})

# 0) Prerequisites : objects built by script 3 --------------------------------
#    prey_groups, prey_groups_by_predator, thresholds, chosen_ranks, phyla_list
#    (script 3) and dat_classed (script 5). Source only if something is missing.
needed <- c("prey_groups", "prey_groups_by_predator", "thresholds",
            "phyla_list", "dat_classed")
if (!all(vapply(needed, exists, logical(1)))) {
  source("3_taxonomic_groups_ll20260812.R")
  if (!exists("dat_classed")) load("data/dat_classed.rda")
}
setDT(prey_groups); setDT(prey_groups_by_predator); setDT(dat_classed)

pcols1 <- paste0("prey_category_", thresholds, "_1")   # _1 columns in BOTH objects

# Predator set + parameters (as in 3b) ----------------------------------------
top10 <- c("Atlantic cod", "American plaice", "Yellowtail flounder",
           "Thorny skate", "Capelin", "Winter flounder", "White hake",
           "Longhorn sculpin", "Atlantic poacher", "Witch flounder")
top5          <- top10[1:5]
floor_stom    <- 5L
period_map_PT <- list("2004-2006" = c(2004, 2005, 2006), "2018-2019" = c(2018, 2019))
P1 <- names(period_map_PT)[1]
P2 <- names(period_map_PT)[2]

# Pooled label lookup (one row per prey taxon), reused throughout.
pooled_lu <- unique(prey_groups[, .SD, .SDcols = c("prey_species_common_name", pcols1)],
                    by = "prey_species_common_name")

# ============================================================================ #
# A) RESOLUTION : number of prey categories per predator, pooled vs per-predator
# ---------------------------------------------------------------------------- #
# Scope where the two schemes can differ: identified prey within phyla_list.
# For each (predator, threshold) we count DISTINCT pooled labels vs DISTINCT
# per-predator labels, using ONLY the two script-3 objects.
sep_lab_long <- melt(
  prey_groups_by_predator[prey_category == "identified_prey_species" & phylum %in% phyla_list,
                          .SD, .SDcols = c("predator", "prey_species_common_name", pcols1)],
  id.vars      = c("predator", "prey_species_common_name"),
  measure.vars = pcols1, variable.name = "col", value.name = "sep_lab")
sep_lab_long[, threshold := as.integer(sub("prey_category_(\\d+)_1", "\\1", as.character(col)))]
sep_lab_long[, col := NULL]

pooled_lab_long <- melt(pooled_lu, id.vars = "prey_species_common_name",
                        measure.vars = pcols1, variable.name = "col", value.name = "pooled_lab")
pooled_lab_long[, threshold := as.integer(sub("prey_category_(\\d+)_1", "\\1", as.character(col)))]
pooled_lab_long[, col := NULL]

# Attach the pooled label to each (predator, taxon, threshold) the predator eats.
lab_cmp <- pooled_lab_long[sep_lab_long, on = c("prey_species_common_name", "threshold")]

cmp <- lab_cmp[, .(n_cat_pooled = uniqueN(pooled_lab),
                   n_cat_sep    = uniqueN(sep_lab)),
               by = .(predator, threshold)]

resolution_cmp <- cmp[, .(n_predators  = .N,
                          pooled_med   = median(n_cat_pooled),
                          sep_med      = median(n_cat_sep),
                          pooled_total = sum(n_cat_pooled),
                          sep_total    = sum(n_cat_sep),
                          gap_med      = median(n_cat_pooled - n_cat_sep)),
                      by = threshold][order(threshold)]
print(resolution_cmp[threshold %in% c(10, 50, 100, 200, 300, 440, 600, 1000)])

# Figure A : per-predator resolution, pooled vs per-predator -------------------
plot_dt <- melt(cmp, id.vars = c("predator", "threshold"),
                measure.vars = c("n_cat_pooled", "n_cat_sep"),
                variable.name = "scheme", value.name = "n_cat")
plot_dt[, scheme := factor(scheme, levels = c("n_cat_pooled", "n_cat_sep"),
                           labels = c("Pooled (all predators)", "Per predator"))]
plot_sum <- plot_dt[, .(med = median(n_cat),
                        q25 = quantile(n_cat, 0.25),
                        q75 = quantile(n_cat, 0.75)),
                    by = .(scheme, threshold)]
print(
  ggplot(plot_sum, aes(threshold, med, color = scheme, fill = scheme)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 1000, 100)) +
    labs(x = "Threshold X (n_stomach)", y = "Prey categories per predator",
         color = "Grouping scheme", fill = "Grouping scheme",
         title = "Fig A - Per-predator taxonomic resolution: pooled vs per predator",
         subtitle = "Line = median predator; ribbon = interquartile range across predators") +
    theme_minimal()
)
print(
  ggplot(plot_dt[predator %in% top5], aes(threshold, n_cat, color = scheme)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ predator) +
    scale_x_continuous(breaks = seq(0, 1000, 200)) +
    labs(x = "Threshold X (n_stomach)", y = "Prey categories", color = "Grouping scheme",
         title = "Fig A2 - Resolution for the top 5 predators: pooled vs per predator",
         subtitle = "One panel per predator") +
    theme_minimal()
)

# ============================================================================ #
# B/C) TROPHIC INDICES + DOMINANT-PREY LOSS
# ---------------------------------------------------------------------------- #
# Same downstream metrics as 3b, computed per predator x threshold x scheme x
# mode. Tools are COPIED from the 6a engine (autonomous, exactly as in 3b).
shannon_from_vec <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  p <- x / sum(x); -sum(p * log(p))
}
levins_bs_from_vec <- function(x, n_cat) {
  x <- as.numeric(x); x <- x[is.finite(x) & x > 0]
  if (length(x) == 0 || n_cat < 2) return(NA_real_)
  if (length(x) < 2) return(0)
  p <- x / sum(x); ((1 / sum(p^2)) - 1) / (n_cat - 1)
}
build_bio <- function(d, per, col) {
  x <- if (is.null(per)) d[!is.na(get(col)) & !is.na(set)]
  else d[period == per & !is.na(get(col)) & !is.na(set)]
  if (nrow(x) == 0) return(NULL)
  x[, set_uid := paste(year, `vessel.code`, set, sep = "_")]
  m <- x[, .(w = sum(somatic_wt_g, na.rm = TRUE)), by = .(set_uid, prey = get(col))]
  m <- dcast(m, set_uid ~ prey, value.var = "w", fill = 0)
  rn <- m$set_uid; mm <- as.matrix(m[, -1, with = FALSE]); rownames(mm) <- rn
  mm <- mm[rowSums(mm) > 0, , drop = FALSE]
  if (nrow(mm) == 0) NULL else mm
}
build_occ <- function(d, per, col, n_min = 5L) {
  x <- if (is.null(per)) d[!is.na(get(col)) & !is.na(set) & !is.na(stomach_id)]
  else d[period == per & !is.na(get(col)) & !is.na(set) & !is.na(stomach_id)]
  if (nrow(x) == 0 || uniqueN(x$stomach_id) < n_min) return(NULL)
  x[, set_uid := paste(year, `vessel.code`, set, sep = "_")]
  denom <- unique(x[, .(set_uid, stomach_id)])[, .(n_sto = .N), by = set_uid]
  occ <- unique(x[, .(set_uid, stomach_id, prey = get(col))])[, .(freq = .N), by = .(set_uid, prey)]
  occ <- merge(occ, denom, by = "set_uid"); occ[, prop := freq / n_sto]
  m <- dcast(occ, set_uid ~ prey, value.var = "prop", fill = 0)
  rn <- m$set_uid; mm <- as.matrix(m[, -1, with = FALSE]); rownames(mm) <- rn
  mm <- mm[rowSums(mm) > 0, , drop = FALSE]
  if (nrow(mm) == 0) NULL else mm
}
align_mats <- function(m1, m2) {
  if (is.null(m1) || is.null(m2)) return(list(m1 = NULL, m2 = NULL))
  m1 <- as.matrix(m1); m2 <- as.matrix(m2); allc <- union(colnames(m1), colnames(m2))
  pad <- function(m) {
    miss <- setdiff(allc, colnames(m))
    if (length(miss)) m <- cbind(m, matrix(0, nrow(m), length(miss),
                                           dimnames = list(rownames(m), miss)))
    m[, allc, drop = FALSE]
  }
  list(m1 = pad(m1), m2 = pad(m2))
}
# run_composition_set allege : BC + R2 seulement (perm = 1, effets seuls).
run_composition_set <- function(s1, s2, mode = c("biomass", "occurrence"),
                                 perm = 1L, relative = TRUE) {
  mode <- match.arg(mode)
  empty <- list(BC = NA_real_, R2 = NA_real_)
  al <- align_mats(s1, s2)
  if (is.null(al$m1) || is.null(al$m2) || nrow(al$m1) < 2 || nrow(al$m2) < 2) return(empty)
  m1 <- al$m1; m2 <- al$m2
  combined <- rbind(m1, m2)
  groups <- factor(c(rep("P1", nrow(m1)), rep("P2", nrow(m2))))
  d <- vegdist(if (relative) decostand(combined, "total") else combined, "bray")
  tryCatch({
    pr <- adonis2(d ~ groups, permutations = perm)
    c1 <- colSums(if (relative) decostand(m1, "total") else m1); c1 <- c1 / sum(c1)
    c2 <- colSums(if (relative) decostand(m2, "total") else m2); c2 <- c2 / sum(c2)
    BC <- as.numeric(vegdist(rbind(c1, c2), "bray"))
    list(BC = round(BC, 3), R2 = round(pr$R2[1], 3))
  }, error = function(e) empty)
}
# n_cat / H / Bs a partir d'une matrice set x proie (diete poolee = colSums).
div_from_mat <- function(m) {
  if (is.null(m)) return(list(n_cat = NA_integer_, H = NA_real_, Bs = NA_real_))
  v <- colSums(m); v <- v[v > 0]; nc <- length(v)
  list(n_cat = nc, H = shannon_from_vec(v), Bs = levins_bs_from_vec(v, n_cat = nc))
}

# Analysis data : dat_classed already carries period, size_class, somatic_wt_g,
# stomach_id, set, year, vessel.code and prey_species_common_name.
dt0 <- copy(dat_classed)
if (inherits(dt0, "sf")) dt0 <- sf::st_drop_geometry(dt0)
dt0 <- as.data.table(dt0)
dt0 <- dt0[!is.na(period) & !is.na(somatic_length_cm) & somatic_length_cm > 0]

# Boucle : predateurs x seuils x schema x mode --------------------------------
res_idx <- list(); res_loss <- list()
for (p in top10) {
  dp <- dt0[predator_species_common_name == p & size_class == "adult"]
  if (nrow(dp) == 0) next

  # Per-predator label lookup for THIS predator (one row per prey taxon).
  bp <- prey_groups_by_predator[predator == p]

  for (x in thresholds) {
    pcol <- paste0("prey_category_", x, "_1")
    if (!pcol %in% names(pooled_lu) || !pcol %in% names(bp)) next

    # Attach BOTH labelings to the SAME diet records: pooled from prey_groups,
    # per-predator from prey_groups_by_predator (fall back to pooled if a taxon
    # is missing from this predator's per-predator table).
    dp[, pooled_lab := pooled_lu[[pcol]][match(prey_species_common_name,
                                               pooled_lu$prey_species_common_name)]]
    dp[, sep_lab := bp[[pcol]][match(prey_species_common_name, bp$prey_species_common_name)]]
    dp[is.na(sep_lab), sep_lab := pooled_lab]

    # (B) indices H / Bs / n_cat / composition, par schema et par mode
    for (scheme in c("pooled", "separate")) {
      col <- if (scheme == "pooled") "pooled_lab" else "sep_lab"
      for (md in c("biomass", "occurrence")) {
        bf <- if (md == "biomass") build_bio else build_occ
        dg <- div_from_mat(bf(dp, per = NULL, col = col))
        comp <- run_composition_set(bf(dp, per = P1, col = col),
                                    bf(dp, per = P2, col = col), mode = md, perm = 1L)
        res_idx[[length(res_idx) + 1L]] <- data.table(
          predator = p, threshold = x, scheme = scheme, mode = md,
          n_cat = dg$n_cat, H = dg$H, Bs = dg$Bs, BC = comp$BC, R2 = comp$R2)
      }
    }

    # (C) perte de la proie dominante : dominante sous POOLED, survie du label
    #     sous PER-PREDATEUR (loss = 1 - part conservee).
    dd <- dp[!is.na(pooled_lab) & !is.na(sep_lab)]
    for (md in c("biomass", "occurrence")) {
      w <- if (md == "biomass")
        dd[!is.na(somatic_wt_g), .(w = sum(somatic_wt_g)), by = .(cat1 = pooled_lab, cat2 = sep_lab)]
      else
        unique(dd[!is.na(stomach_id), .(stomach_id, cat1 = pooled_lab, cat2 = sep_lab)])[, .(w = .N), by = .(cat1, cat2)]
      if (nrow(w) == 0) next
      tot1 <- w[, .(tot = sum(w)), by = cat1]
      L1 <- tot1[which.max(tot), cat1]
      wl <- w[cat1 == L1]
      preserved <- sum(wl[cat2 == L1, w]) / sum(wl$w)
      res_loss[[length(res_loss) + 1L]] <- data.table(
        predator = p, threshold = x, mode = md, dom_prey = L1,
        dom_share = tot1[cat1 == L1, tot] / sum(tot1$tot),
        loss = 1 - preserved, is_lost = preserved < 0.5)
    }
  }
  message("[3c] pooled vs per-predateur termine pour : ", p)
}
idx <- rbindlist(res_idx); loss <- rbindlist(res_loss)

# Resumes inter-predateurs ----------------------------------------------------
idx_long <- melt(idx, id.vars = c("predator", "threshold", "scheme", "mode"),
                 measure.vars = c("n_cat", "H", "Bs", "BC", "R2"),
                 variable.name = "metric", value.name = "value")
idx_summ <- idx_long[!is.na(value),
                     .(med = median(value),
                       q25 = as.numeric(quantile(value, .25)),
                       q75 = as.numeric(quantile(value, .75)),
                       n_pred = .N),
                     by = .(mode, metric, scheme, threshold)]
idx_summ[, scheme := factor(scheme, levels = c("pooled", "separate"),
                            labels = c("Pooled (all predators)", "Per predator"))]

loss_summ <- loss[!is.na(loss),
                  .(loss_med = median(loss),
                    loss_q25 = as.numeric(quantile(loss, .25)),
                    loss_q75 = as.numeric(quantile(loss, .75)),
                    lost_rate = mean(is_lost), n_pred = .N),
                  by = .(mode, threshold)]

cat("\n-- Perte de la proie dominante (pooled -> per-predator) : mediane inter-predateurs --\n")
print(dcast(loss_summ[threshold %in% c(10, 100, 200, 500, 1000)],
            threshold ~ mode, value.var = c("loss_med", "lost_rate")))

cat("\n-- Indices : mediane inter-predateurs (pooled vs per-predator, seuils choisis) --\n")
print(dcast(idx_summ[metric %in% c("n_cat", "H", "Bs", "BC", "R2") &
                       threshold %in% c(10, 100, 500, 1000)],
            mode + metric + threshold ~ scheme, value.var = "med"))

# Figure C : perte de la proie dominante (pooled -> per-predator) --------------
print(
  ggplot(loss_summ, aes(threshold, loss_med, color = mode, fill = mode)) +
    geom_ribbon(aes(ymin = loss_q25, ymax = loss_q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    labs(x = "Seuil X (n_stomach)",
         y = "Part de la proie dominante (pooled) reassignee sous per-predator",
         color = "Monnaie", fill = "Monnaie",
         title = "Fig C - Perte de la proie dominante : pooled (fin) -> per-predator (grossier)",
         subtitle = "Mediane inter-predateurs (ruban = interquartile) ; dominante = 1re proie de la diete poolee") +
    theme_minimal()
)

# Figures B : effet du schema (pooled vs per-predator) sur n_cat, H, Bs --------
idx_fig <- function(mt, ylab, ttl) {
  ggplot(idx_summ[metric == mt], aes(threshold, med, color = scheme, fill = scheme)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ mode, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    labs(x = "Seuil X (n_stomach)", y = ylab, color = "Regroupement", fill = "Regroupement",
         title = ttl,
         subtitle = "Mediane inter-predateurs (ruban = interquartile) ; 10 predateurs ; adult ; toutes regions ; PT") +
    theme_minimal()
}
print(idx_fig("n_cat", "n_cat - diete poolee",
              "Fig B1 - Effet du schema sur n_cat (pooled vs per-predator)"))
print(idx_fig("H",  "H (Shannon, nats) - diete poolee",
              "Fig B2 - Effet du schema sur H (pooled vs per-predator)"))
print(idx_fig("Bs", "Bs (Levins standardise, 0..1)",
              "Fig B3 - Effet du schema sur Bs (pooled vs per-predator)"))
print(
  ggplot(idx_summ[metric %in% c("BC", "R2")],
         aes(threshold, med, color = scheme, fill = scheme)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_grid(metric ~ mode, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    labs(x = "Seuil X (n_stomach)", y = "Mediane inter-predateurs (ruban = interquartile)",
         color = "Regroupement", fill = "Regroupement",
         title = "Fig B4 - Effet du schema sur la composition P1 vs P2 (BC / R2)",
         subtitle = "run_composition_set entre 2004-2006 et 2018-2019 ; adult ; toutes regions") +
    theme_minimal()
)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# OBSERVATION — POOLÉ vs PAR PRÉDATEUR POUR LES ANALYSES EN AVAL (6a_Engine_Trophic.R)
#------------------------------------------------------------------------------#
#
# Le moteur (6a_Engine_Trophic.R) est déjà construit PAR PRÉDATEUR : run_pipeline()
# définit des cellules predator x SPATIAL_OUT x size_class x period (~ ligne 807)
# puis bâtit la matrice pour UN SEUL prédateur (sp) avant de tester P1 vs P2 par
# PERMANOVA (~ lignes 829, 836-840). Chaque adonis2(d ~ période) compare donc les
# deux périodes À L'INTÉRIEUR d'une espèce.
#
# RECOMMANDATION : garder les matrices SÉPARÉES PAR PRÉDATEUR pour comparer les
# deux périodes. Raisons :
#   1. Écologique — les espèces ont des diètes fondamentalement différentes ;
#      une matrice tous prédateurs confondus serait dominée par les différences
#      ENTRE espèces, pas par le signal temporel recherché.
#   2. Confusion / paradoxe de Simpson — si la proportion de prédateurs
#      échantillonnés diffère entre P1 et P2, un "effet période" sur la matrice
#      poolée serait confondu avec un changement de composition d'échantillonnage
#      des prédateurs (faux changement de diète).
#   3. Dispersion multivariée — pooler fusionne plusieurs nuages de diète
#      distincts, ce qui gonfle la dispersion et fait confondre à la PERMANOVA
#      décalage de localisation et différence de dispersion (Anderson 2006 ;
#      Warton et al. 2012).
#
# POOLER n'est légitime que si la question est explicitement au niveau du
# GUILDE/communauté de prédateurs. Dans ce cas, ne PAS tester la période seule :
# mettre le prédateur comme facteur/strate, p. ex.
#   adonis2(d ~ predator + period, by = "terms")   # ou strata = predator
# afin de tester l'effet période À L'INTÉRIEUR de chaque prédateur. Une matrice
# poolée testée sur ~ period seul n'est pas valide (pseudoréplication,
# Hurlbert 1984).
#
# NB : axe DIFFÉRENT de la section A/B/C ci-dessus. Ces sections portent sur la
# DÉFINITION des catégories de proies (regroupement taxonomique poolé vs par
# prédateur) ; ici il s'agit de la MATRICE D'ANALYSE. Les deux devraient rester
# par prédateur pour rester cohérents.
#
# Sources :
#   Anderson (2001) Austral Ecology 26:32-46            [PERMANOVA]
#   Anderson (2006) Biometrics 62:245-253               [PERMDISP]
#   Warton, Wright & Wang (2012) Methods Ecol Evol 3:89-101
#   Hurlbert (1984) Ecological Monographs 54:187-211    [pseudoréplication]
#   (bloc LITERATURE CITED de 6a_Engine_Trophic_ll20260806.R, ~ ligne 1577)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# OBSERVATION de l'IA non entièrement révisé par moi LL20260813 —
# H ET Bs : CALCULER PAR PRÉDATEUR, PAS SUR UNE MATRICE POOLÉE
#------------------------------------------------------------------------------#
#
# Contrairement à BC/R2 (qui COMPARENT deux périodes), H (Shannon) et surtout Bs
# (largeur de niche de Levins) sont des descripteurs AU NIVEAU DU CONSOMMATEUR :
# ils décrivent comment UN prédateur répartit son alimentation entre les
# catégories de proies.
#
#   * Bs (Levins) — la largeur de niche est DÉFINIE pour un consommateur donné
#     (Levins 1968 ; standardisation de Hurlbert 1978). Une valeur poolée décrit
#     la niche d'un "super-prédateur" fictif et mélange la largeur de niche
#     INTRA-prédateur avec les différences de diète INTER-prédateurs — la logique
#     même de la décomposition WIC/TNW (Roughgarden 1972 ; Bolnick et al. 2002).
#   * H (Shannon) — même confusion intra/inter, et la matrice serait dominée par
#     l'espèce la plus échantillonnée (la morue ~ 40 % des estomacs), donc le H
#     "global" refléterait surtout la diète de la morue, pas celle du guilde.
#
# CONCLUSION : garder TOUT par prédateur.
#     Bs, H, n_cat  -> par consommateur, par période -> PAR PRÉDATEUR
#     BC, R2        -> comparaison entre périodes     -> PAR PRÉDATEUR
# C'est déjà la structure du moteur (cellules predator x ... x period), donc Bs
# et H sont calculés par prédateur ET par période -> comparables entre P1 et P2.
# Une valeur poolée n'a de sens que comme statistique explicite DE GUILDE, jamais
# en substitut des Bs/H par espèce.
#
# Sources :
#   Levins (1968) Evolution in Changing Environments   [largeur de niche B]
#   Hurlbert (1978) Ecology 59:67-77                    [standardisation de B]
#   Roughgarden (1972) Am. Nat. 106:683-718             [WIC / TNW]
#   Bolnick et al. (2002) Ecology 83:2936-2941          [WIC/TNW, forme discrète]
#   (bloc LITERATURE CITED de 6a_Engine_Trophic_ll20260806.R, ~ ligne 1577)
#------------------------------------------------------------------------------#
