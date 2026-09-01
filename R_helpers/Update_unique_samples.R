# =============================================================================
# Update_unique_samples.R - COMPLETE THE PREDATOR-WEIGHT CACHE
# -----------------------------------------------------------------------------
# data/unique_samples.rda holds one estimated weight per (species, length,
# year). After a change in the length correction, new combinations appear in
# diet_clean; this script estimates only those and appends them to the cache.
# Run once from the project root, then rerun 4_Database_with_others_factors.R.
# =============================================================================
library(data.table)
library(gulf)

load("data/unique_samples.rda"); setDT(unique_samples)
load("data/diet_clean.rda");     setDT(diet_clean)

combos <- diet_clean[!is.na(predator_species_code) & !is.na(somatic_length_cm) & !is.na(year),
                     .(n_fish = .N), by = .(predator_species_code, somatic_length_cm, year)]
missing <- combos[!unique_samples, on = .(predator_species_code, somatic_length_cm, year)]
cat("Cached:", nrow(unique_samples), "| in diet_clean:", nrow(combos),
    "| to estimate:", nrow(missing), "\n")

if (nrow(missing) > 0) {
  missing[, pred_weight_est := tryCatch({
    w <- weight(somatic_length_cm, species = predator_species_code, year = year) * 1000
    if (length(w) > 0) w[1] else NA_real_
  }, error = function(e) NA_real_), by = seq_len(nrow(missing))]
  missing[year %in% c(2018, 2019) & predator_species_code == 60,
          pred_weight_est := weight(somatic_length_cm, species = predator_species_code,
                                    year = 2017) * 1000]
  unique_samples <- rbind(unique_samples, missing, fill = TRUE)
  save(unique_samples, file = "data/unique_samples.rda")
  cat("Cache updated:", nrow(unique_samples), "combinations\n")
}

na_left <- unique_samples[is.na(pred_weight_est), .N, by = predator_species_code][order(-N)]
if (nrow(na_left)) { cat("Combinations without a weight estimate, by species code:\n"); print(na_left) }
