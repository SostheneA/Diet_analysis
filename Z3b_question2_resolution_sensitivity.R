#------------------------------------------------------------------------------#
# 3b_question2_resolution_sensitivity.R - Diagnostic de resolution taxonomique : pool vs analyse par predateur
#
# Objectif :
#   Documenter et quantifier le compromis de resolution introduit par le
#   regroupement taxonomique (script 3), qui est calcule TOUS PREDATEURS
#   CONFONDUS (poole), alors que les analyses (6a et suivantes) sont faites
#   PAR PREDATEUR.
#
# Trois constats :
#
#   1) n_predator_threshold = 2.
#      La regle "garder le rang fin seulement si >= 2 predateurs consomment
#      la proie". Une proie propre a un seul predateur
#      (meme dominante) est remontee a un rang grossier uniquement parce
#      qu'aucun autre predateur ne la mange. Cela efface l'information la
#      plus discriminante entre predateurs. Mettre le seuil a 1 rend le
#      critere vide et supprime ce biais.
#
#   2) Collapse espece -> genre.
#      Il rabaisse au genre une espece qui avait DEJA passe le
#      seuil au rang espece, ce qui contredit la logique "garder le rang le
#      plus fin qui atteint le seuil". Sans cette etape, la regle
#      produit naturellement "Genus species" + "Genus_others".
#
#   3) Le regroupement de prédateurs est une probablement une SUR-resolution
#      par predateur : des categories fines mais quasi vides
#      (peu d'estomacs) pour un predateur donné.
#
# Input  : data/dat_classed.rda   (produit par le script 5 ; contient les
#          colonnes prey_category_X_1 (corrigee) ET _X_2 (origine) pour
#          X = 10..1000, ainsi que stomach_id, predator_species_common_name,
#          size_class et year)
#
# Output : Fig A  p_over_res         - sur-resolution PAR PREDATEUR (_1 vs _2)
#          Fig B  p_cell_removal     - retrait PAR CELLULE (famille _1)
#          Fig C  p_cell_removal_sp  - retrait par cellule ET par predateur
#          objets : res_resolution, cell_removal, cell_removal_sp, per_cell
#------------------------------------------------------------------------------#

library(data.table)
library(ggplot2)
library(here)
library(tidyverse)

# 1) Donnees -------------------------------------------------------------------
lines <- readLines(paste0(here::here(), "/3_taxonomic_groups_ll20260812.R"),
                   warn = FALSE)[39:310]
eval(parse(text = lines), envir = .GlobalEnv)

lines <- readLines(paste0(here::here(), "/4_Database_with_others_factors.R"),
                   warn = FALSE)[17:93]
eval(parse(text = lines), envir = .GlobalEnv)

diet_clean <- copy(Database)

lines <- readLines(paste0(here::here(), "/5_From_cleandiet_to_dat_classed.R"),
                   warn = FALSE)[85:129]
eval(parse(text = lines), envir = .GlobalEnv)

dat_classed_2 <- dat_classed

# 2) Parametres du diagnostic --------------------------------------------------
if(!exists("thresholds")) thresholds <- seq(10, 1000, by = 10)   # seuils X presents dans les colonnes _2
floor_stom <- 5L                       # "quasi vide" : < floor_stom estomacs
top10 <- c("Atlantic cod", "American plaice", "Yellowtail flounder",
                   "Thorny skate", "Capelin",
                   "Winter flounder", "White hake", "Longhorn sculpin",
                   "Atlantic poacher", "Witch flounder")
top5 <- c("Atlantic cod", "American plaice", "Yellowtail flounder",
           "Thorny skate", "Capelin")

top15 <- c("Atlantic cod", "American plaice", "Yellowtail flounder",
           "Thorny skate", "Capelin",
           "Winter flounder", "White hake", "Longhorn sculpin",
           "Atlantic poacher", "Witch flounder", "Greenland halibut", "Atlantic herring",
           "Alewife", "Rainbow smelt", "Atlantic mackerel")

x <- setdiff(unique(dat_classed_2$predator_species_common_name), top10)
# 3) Taux de colonnes quasi vides PAR PREDATEUR, selon le seuil et la famille ---------
# Pour chaque seuil, chaque famille (_1 corrigee, _2 origine) et chaque predateur :
#   n_stomachs     = nb d'estomacs uniques du predateur contenant la categorie de proie
#   low_count_rate = part des categories PRESENTES (n_stomachs > 0) avec n_stomachs < floor_stom
#   singleton_rate = part des categories PRESENTES vues dans un seul estomac (n_stomachs == 1)
sparse_by_pred <- function(diet_family) {
  results <- lapply(thresholds, function(threshold) {
    category_col <- paste0("prey_category_", threshold, "_", diet_family)

    stomach_categories <- unique(dat_classed_2[, .(
      predator_species_common_name,
      stomach_id,
      prey_category = get(category_col)
    )])
    stomach_categories <- stomach_categories[!is.na(prey_category)]

    category_counts <- stomach_categories[, .(n_stomachs = .N),
                                          by = .(predator_species_common_name, prey_category)]

    predator_rates <- category_counts[, .(
      low_count_rate = mean(n_stomachs < floor_stom),
      singleton_rate = mean(n_stomachs == 1)
    ), by = predator_species_common_name]

    # On identifie la famille/le seuil pour pouvoir les distinguer une fois empilees
    predator_rates[, `:=`(family = paste0("_", diet_family), threshold = threshold)]

    summary_row <- data.table(
      family      = paste0("_", diet_family),
      threshold   = threshold,
      n_cat_total = uniqueN(category_counts$prey_category),
      rate_median = median(predator_rates$low_count_rate),
      rate_q90    = as.numeric(quantile(predator_rates$low_count_rate, 0.90)),
      sing_median = median(predator_rates$singleton_rate)
    )

    list(summary = summary_row, predator_rates = predator_rates)
  })

  list(
    summary        = rbindlist(lapply(results, `[[`, "summary")),
    predator_rates = rbindlist(lapply(results, `[[`, "predator_rates"))
  )
}
resolution_1 <- sparse_by_pred(1)
resolution_2 <- sparse_by_pred(2)

res_resolution_s <- rbind(resolution_1$summary, resolution_2$summary)
res_resolution_p <- rbind(resolution_1$predator_rates, resolution_2$predator_rates)

print(res_resolution_s[threshold %in% c(10, 50, 100, 200, 300, 440, 600, 1000)][
        order(family , threshold)])

# 4) Figure Ai : sur-resolution par predateur, famille _1 vs _2 -----------------
res_resolution_s <- res_resolution_s |>
  mutate(family = forcats::fct_relevel(family, "_2", after = Inf))


p_over_res <- ggplot(res_resolution_s, aes(x = threshold, color = family, fill = family )) +
  geom_ribbon(aes(ymin = rate_median, ymax = rate_q90), alpha = 0.12, color = NA) +
  geom_line(aes(y = rate_median), linewidth = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(0, 1000, 100)) +
  labs(
    x = "Seuil X (n_stomach)",
    y = paste0("Part des categories presentes avec < ", floor_stom, " estomacs"),
    color = "Famille", fill = "Famille",
    title = "Fig A - Sur-resolution par predateur selon le seuil et la famille (n_predator)",
    subtitle = "Trait = mediane entre predateurs ; ruban = mediane -> 90e centile"
  ) +
  theme_minimal()

print(p_over_res)


# 4) Figure Aii : sur-resolution par predateur, famille _1 vs _2 -----------------
p_over_res_p <- ggplot(res_resolution_p[predator_species_common_name %in% top15,],
                       aes(x = threshold, color = family, fill = family, linetype = family)) +

  geom_ribbon(aes(ymin = singleton_rate, ymax = low_count_rate), alpha = 0.08, color = NA) +
    geom_line(aes(y = low_count_rate), linewidth = 0.75, alpha = 0.7) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(0, 1000, 100)) +
  #facet_grid(~ predator_species_common_name, scales = "free_y") +
  facet_wrap(~ predator_species_common_name, nrow = 3, scales = "free_y") +
  labs(
    x = "Seuil X (n_stomach)",
    y = paste0("Part des categories presentes avec < ", floor_stom, " estomacs"),
    title = "Fig Aii - Sur-resolution par predateur selon le seuil et la famille",
    subtitle = paste0("Trait = part des categories <", floor_stom, " ; ruban = part des categories singletons (n=1) -> < ", floor_stom)
  ) +

  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_over_res_p)


#if (!dir.exists("figures")) dir.create("figures")
#ggsave("figures/6g_over_resolution_vs_threshold.png", p_over_res,
#       width = 8, height = 5, dpi = 150)

#------------------------------------------------------------------------------#
# Lecture :
#   - Augmenter X reduit la sur-resolution lentement:
#     le plancher du 90e centile reste eleve car X ne borne pas S_{g,p}.
#   - Seul un seuil applique PAR PREDATEUR (S_{g,p} >= plancher), ou un
#     nettoyage/regroupement des colonnes trop creuses au moment de construire
#     chaque matrice dans 6a, corrige reellement le probleme -- tout en gardant
#     le schema poole comme jeu de colonnes commun (comparabilite pour 6a, mais voir la suite ci-bas).
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# 8) PERTE DE LA PROIE DOMINANTE (_1 -> _2) + effet sur Bs / H / composition
#
#   Question : en passant de la famille _1 (corrigee : n_predator = 1, PAS de
#   collapse espece->genre -> plus FINE) a la famille _2 (origine : n_predator
#   = 2, avec collapse -> plus GROSSIERE), (a) perd-on la proie DOMINANTE de
#   chaque predateur, et (b) quel est l'effet sur les indices trophiques ?
#
#   (a) Perte de la proie dominante -- Option B "survie de la dominante fine".
#       Pour chaque cellule (predateur x seuil X), diete poolee (2 periodes
#       confondues) ; on identifie la proie DOMINANTE sous _1 (1re proie par la
#       biomasse ou occurence du mode) puis la part de cette proie reassignee a une AUTRE
#       etiquette sous _2 (loss = 1 - part conservee). Monnaie : biomasse =
#       somme somatic_wt_g ; occurrence = nb d'estomacs distincts.
#
#   (b) Effet sur Bs / H / composition, par famille (_1 vs _2) :
#       H = shannon_from_vec, Bs = levins_bs_from_vec sur la diete poolee ;
#       BC / R2 = run_composition_set entre P1 (2004-2006) et P2 (2018-2019).
#       p-values non rapportees (effets seuls, perm = 1).
#
#   Cadre 3b  : top10 (predators) ; adult ; TOUTES regions
#   (ar = NULL) ; period_map_PT ; 100 seuils. Outils COPIES du moteur 6a mais
#   AUTONOMES. Les
#   fonctions d'indices sont reprises EXACTEMENT du moteur.
#------------------------------------------------------------------------------#
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(vegan) })

# Outils du moteur 6a (COPIES ; autonomes) -----------
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

# Donnees : dat_classed (porte _1, _2 ET les colonnes d'analyse) ----------
#if (!exists("dat_classed")) load(file.path(here::here(), "data/dat_classed.rda"))
period_map_PT <- list("2004-2006" = c(2004, 2005, 2006), "2018-2019" = c(2018, 2019))
P1 <- names(period_map_PT)[1]
P2 <- names(period_map_PT)[2]

dt0 <- copy(dat_classed_2)
if (inherits(dt0, "sf")) dt0 <- sf::st_drop_geometry(dt0)

dt0 <- as.data.table(dt0)
dt0 <- dt0[year %in% unlist(period_map_PT) & !is.na(somatic_length_cm) & somatic_length_cm > 0]
dt0[, period := fifelse(year %in% period_map_PT[[1]], P1,
                        fifelse(year %in% period_map_PT[[2]], P2, NA_character_))]
dt0 <- dt0[!is.na(period)]
if (!exists("thresholds")) thresholds <- seq(10, 1000, by = 10)

# Boucle : predateurs x seuils --------------------------------------------
res_idx <- list(); res_loss <- list()
for (p in top10) {
  dp <- dt0[predator_species_common_name == p & size_class == "adult"]
  if (nrow(dp) == 0) next
  for (x in thresholds) {
    c1c <- paste0("prey_category_", x, "_1")
    c2c <- paste0("prey_category_", x, "_2")
    if (!all(c(c1c, c2c) %in% names(dp))) next

    # indices H / Bs / composition, par famille et par mode
    for (fam in c("1", "2")) {
      col <- paste0("prey_category_", x, "_", fam)

      for (md in c("biomass", "occurrence")) {

        bf <- if (md == "biomass") build_bio else build_occ
        dg <- div_from_mat(bf(dp, per = NULL, col = col))
        comp <- run_composition_set(bf(dp, per = P1, col = col),
                                    bf(dp, per = P2, col = col), mode = md, perm = 1L)

        res_idx[[length(res_idx) + 1L]] <- data.table(
          predator = p,
          threshold = x,
          family = paste0("_", fam),
          mode = md,
          n_cat = dg$n_cat,
          H = dg$H,
          Bs = dg$Bs,
          BC = comp$BC,
          R2 = comp$R2)
      }
    }

    # perte de la proie dominante : dominante sous _1, survie du label sous _2
    dd <- dp[!is.na(get(c1c)) & !is.na(get(c2c))]
    dd[, `:=`(cat1 = get(c1c), cat2 = get(c2c))]

    for (md in c("biomass", "occurrence")) {
      w <- if (md == "biomass")
        dd[!is.na(somatic_wt_g),
           .(w = sum(somatic_wt_g)),
           by = .(cat1, cat2)]
      else
        unique(dd[!is.na(stomach_id),
                  .(stomach_id, cat1, cat2)])[, .(w = .N), by = .(cat1, cat2)]

      if (nrow(w) == 0) next

      tot1 <- w[, .(tot = sum(w)),
                by = cat1]
      L1 <- tot1[which.max(tot),
                 cat1]
      wl <- w[cat1 == L1]

      preserved <- sum(wl[cat2 == L1, w]) / sum(wl$w)

      res_loss[[length(res_loss) + 1L]] <- data.table(
        predator = p,
        threshold = x,
        mode = md,
        dom_prey = L1,
        dom_share = tot1[cat1 == L1, tot] / sum(tot1$tot),
        loss = 1 - preserved,
        is_lost = preserved < 0.5)
    }
  }
  message("[8] perte dominante + indices termine pour : ", p)
}
idx <- rbindlist(res_idx); loss <- rbindlist(res_loss)

# Resumes inter-predateurs ------------------------------------------------
idx_long <- melt(idx, id.vars = c("predator","threshold","family","mode"),
                 measure.vars = c("n_cat","H","Bs","BC","R2"),
                 variable.name = "metric", value.name = "value")
idx_summ <- idx_long[!is.na(value),
                     .(med = median(value),
                       q25 = as.numeric(quantile(value,.25)),
                       q75 = as.numeric(quantile(value,.75)),
                       n_pred = .N),
                     by = .(mode, metric, family, threshold)]

loss_summ <- loss[!is.na(loss),
                  .(loss_med = median(loss),
                    loss_q25 = as.numeric(quantile(loss,.25)),
                    loss_q75 = as.numeric(quantile(loss,.75)),
                    lost_rate = mean(is_lost), n_pred = .N),
                  by = .(mode, threshold)]

cat("\n-- Perte de la proie dominante : mediane inter-predateurs (seuils choisis) --\n")
print(dcast(loss_summ[threshold %in% c(10,100,200,500,1000)],
            threshold ~ mode, value.var = c("loss_med","lost_rate")))

cat("\n-- Indices : mediane inter-predateurs (_1 vs _2, seuils choisis) --\n")
print(dcast(idx_summ[metric %in% c("n_cat","H","Bs","BC","R2") &
                       threshold %in% c(10,100,500,1000)],
            mode + metric + threshold ~ family, value.var = "med"))

# Figure 8a : perte de la proie dominante (_1 -> _2) -----------------------
print(
  ggplot(loss_summ, aes(threshold, loss_med, color = mode, fill = mode)) +
    geom_ribbon(aes(ymin = loss_q25, ymax = loss_q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    labs(x = "Seuil X (n_stomach)",
         y = "Part de la proie dominante _1 reassignee sous _2",
         color = "Monnaie", fill = "Monnaie",
         title = "Fig 8a - Perte de la proie dominante en passant de _1 (fin) a _2 (grossier)",
         subtitle = "Mediane inter-predateurs (ruban = interquartile) ; dominante = 1re proie de la diete _1") +
    theme_minimal()
)

# Figures 8b / 8c / 8d : effet de la famille sur H, Bs, composition -------
idx_fig <- function(mt, ylab, ttl) {
  ggplot(idx_summ[metric == mt], aes(threshold, med, color = family, fill = family)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ mode, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    labs(x = "Seuil X (n_stomach)", y = ylab, color = "Famille", fill = "Famille",
         title = ttl,
         subtitle = "Mediane inter-predateurs (ruban = interquartile) ; 10 predateurs ; adult ; toutes regions ; PT") +
    theme_minimal()
}
print(idx_fig("H",  "H (Shannon, nats) - diete poolee",
              "Fig 8b - Effet de la famille sur H (_1 vs _2)"))
print(idx_fig("Bs", "Bs (Levins standardise, 0..1)",
              "Fig 8c - Effet de la famille sur Bs (_1 vs _2)"))
print(
  ggplot(idx_summ[metric %in% c("BC","R2")],
         aes(threshold, med, color = family, fill = family)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_grid(metric ~ mode, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    labs(x = "Seuil X (n_stomach)", y = "Mediane inter-predateurs (ruban = interquartile)",
         color = "Famille", fill = "Famille",
         title = "Fig 8d - Effet de la famille sur la composition P1 vs P2 (BC / R2)",
         subtitle = "run_composition_set entre 2004-2006 et 2018-2019 ; adult ; toutes regions") +
    theme_minimal()
)
#------------------------------------------------------------------------------#


# Perte de la proie dominante : 0 % partout (0 cellule sur 1800, biomasse et occurrence, tous seuils). La dominante de chaque prédateur est toujours une catégorie commune (p. ex. Clupea harengus en biomasse, Mysida en occurrence) qui franchit les seuils et n'est jamais absorbée par le collapse _2. Le regroupement _1→_2 ne touche que la queue rare, pas la dominante. La Fig 8a est donc une ligne plate à 0 — c'est le résultat, pas un bug.
# H : _1 légèrement supérieur à _2 (sensible aux proies rares que _1 conserve), écart faible et qui se referme quand X monte.
# Bs et composition (BC, R²) : _1 et _2 quasi superposés — l'effet de la famille est négligeable sur ces indices.
