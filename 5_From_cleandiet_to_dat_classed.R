# =============================================================================
# 5_From_cleandiet_to_dat_classed.R - SIZE CLASSES AND THE ANALYSIS TABLE
# -----------------------------------------------------------------------------
# Input  : data/diet_clean.rda (from script 4)
# Output : data/dat_classed.rda (dat_classed: diet_clean with size_class)
#
# Size class relative to length at maturity (lmat): juvenile below lmat, adult
# at or above; predators without a documented lmat (lmat = lmax in
# lmat_predators.csv) form a single class.
# =============================================================================
rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)
library(sf)

load("data/diet_clean.rda")
setDT(diet_clean)

diet_clean[, `:=`(lmat = as.numeric(lmat), lmax = as.numeric(lmax))]
diet_clean[, size_class := fcase(
  lmat == lmax,               "one_class",
  somatic_length_cm <  lmat,  "juvenile",
  somatic_length_cm >= lmat,  "adult"
)]

size_share <- diet_clean %>%
  st_drop_geometry() %>%
  count(predator_species_common_name, period, size_class) %>%
  group_by(predator_species_common_name, period) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = size_class, values_from = pct, values_fill = 0)
print(size_share, n = Inf)

dat_classed <- diet_clean
save(dat_classed, file = "data/dat_classed.rda")
