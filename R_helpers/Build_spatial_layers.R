# =============================================================================
# Build_spatial_layers.R - GENERATES THE THREE "new_*" SPATIAL LAYERS
# -----------------------------------------------------------------------------
# Self-contained: starts from the source shapefiles and writes
#   data/Spatial_data/new_ecoregions_final.shp   (col. EcoZone, 6 zones)
#   data/Spatial_data/new_regions_v2_final.shp   (col. Zone,    2nd proposal)
#   data/Spatial_data/new_gulf_final.shp         (col. Area,    THE 4 manuscript
#                                                 ecoregions - CORRECTED build)
#
# Sources (relative to the project root):
#   data/Spatial_data/EAR_map/EAR_map.shp        ecological regions (EAR)
#     PID 1 northwest_estuary . PID 3 centre . PID 5 magdalen_shallows
#     PID 6 northumberland_strait . PID 50 baie_des_chaleurs
#   data/Spatial_data/NAFO/nafo_2014_02.shp      NAFO unit areas (level_2 = 4TX)
#
# Corrections vs the original construction of new_gulf_final:
#   1. everything is clipped to the UNION OF THE 4T UNIT AREAS
#      (before: level_0 == 4 = all of division 4 -> 70% of Central outside 4T)
#   2. the EAR strip 'northwest_estuary' (Gaspe north shore + estuary mouth,
#      strata 415-416 area) is included EXPLICITLY (before: 17% of 4T had no
#      ecoregion and its sets were silently attached to the nearest polygon)
#   3. no S2/planar mixing: one geometry engine for the whole script
#   4. loud checks: the script stops if an overlap or a hole appears in the
#      surveyed area (only the unsampled upper estuary, west of -66.5, may
#      remain uncovered)
#
# After regenerating new_gulf_final.shp, rerun 4 -> 5 -> 6c (the ecoregion
# runs), then 7/8/9. 6b (Gulf-wide) and 6d (stratum) do not depend on Area.
# Run from the project root:  source("R_helpers/Build_spatial_layers.R")
# =============================================================================

library(sf)
library(dplyr)

# ---- CHOICE: destination of the 4T part of 'northwest_estuary' (PID 1).
# "Central" keeps the deep-channel logic; "Magdalen Shallows" is the
# alternative; a new name (e.g. "Gaspe") creates a 5th ecoregion.
GASPE_TO <- "Central"

SPATIAL_DIR <- "data/Spatial_data"

sf_use_s2(FALSE)

EAR_sf     <- st_make_valid(st_read(file.path(SPATIAL_DIR, "EAR_map/EAR_map.shp"),   quiet = TRUE))
NAFO_4T_sf <- st_make_valid(st_read(file.path(SPATIAL_DIR, "NAFO/nafo_2014_02.shp"), quiet = TRUE))

# Union of the 4T unit areas: the study domain used to clip every layer.
z_4t <- NAFO_4T_sf %>% filter(grepl("^4T", level_2)) %>% st_union()

# Generic helper: intersection -> merged single-row sf with the given id column.
zone_sf <- function(geom_a, geom_b, name, col) {
  g <- st_union(st_geometry(suppressWarnings(st_intersection(geom_a, geom_b))))
  out <- st_sf(x = name, geometry = g)
  names(out)[1] <- col
  out
}

# =============================================================================
# 1. new_ecoregions_final.shp  (EcoZone - first proposal, 6 zones)
# =============================================================================
z_4tmn  <- NAFO_4T_sf %>% filter(level_2 %in% c("4TM", "4TN")) %>% st_union()
z_4thlg <- NAFO_4T_sf %>% filter(level_2 %in% c("4TH", "4TL", "4TG")) %>% st_union()

eco_chaleurs <- zone_sf(EAR_sf %>% filter(PID == 50), z_4tmn,  "Chaleur Bay", "EcoZone")
eco_northumb <- zone_sf(EAR_sf %>% filter(PID == 6),  z_4thlg, "Northumberland",    "EcoZone")

nafo_4tgj <- NAFO_4T_sf %>% filter(level_2 %in% c("4TG", "4TJ")) %>% st_union()
eco_4tgj  <- st_sf(EcoZone = "4TGJ clean",
                   geometry = st_union(suppressWarnings(
                     st_difference(nafo_4tgj, st_geometry(eco_northumb)))))

z_nafo_otl <- NAFO_4T_sf %>% filter(level_2 %in% c("4TO", "4TN", "4TL")) %>% st_union()
z_ear_mag  <- EAR_sf %>% filter(regn_nm %in% c("centre", "magdalen_shallows")) %>% st_union()
eco_shediac <- zone_sf(z_ear_mag, z_nafo_otl, "Shediac-Magdalen", "EcoZone") %>%
  st_difference(st_geometry(eco_chaleurs)) %>% st_as_sf()

eco_pure <- NAFO_4T_sf %>%
  filter(level_2 %in% c("4TF", "4TK")) %>%
  mutate(EcoZone = paste(level_2, "Only")) %>%
  select(EcoZone)

new_ecoregions_clean <- bind_rows(eco_chaleurs, eco_northumb, eco_4tgj,
                                  eco_shediac, eco_pure) %>%
  filter(!st_is_empty(.)) %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

st_write(new_ecoregions_clean, file.path(SPATIAL_DIR, "new_ecoregions_final.shp"),
         delete_layer = TRUE, quiet = TRUE)
cat("Written:", file.path(SPATIAL_DIR, "new_ecoregions_final.shp"), "\n")

# =============================================================================
# 2. new_regions_v2_final.shp  (Zone - second proposal, with the custom box)
# =============================================================================
box <- st_sfc(st_polygon(list(matrix(c(-64, 45.5,  -61, 45.5,  -61, 46.47,
                                       -64, 46.47, -64, 45.5),
                                     ncol = 2, byrow = TRUE))), crs = 4326)

v2_chaleurs <- zone_sf(EAR_sf %>% filter(PID == 50), z_4tmn, "Chaleur Bay", "Zone")

z_4thg <- NAFO_4T_sf %>% filter(level_2 %in% c("4TH", "4TG")) %>% st_union()
v2_northcoast <- zone_sf(st_as_sf(z_4thg), box, "North Coast", "Zone")

z_4tjg <- NAFO_4T_sf %>% filter(level_2 %in% c("4TJ", "4TG")) %>% st_union()
v2_4tgj <- st_sf(Zone = "4TGJ coast",
                 geometry = st_union(suppressWarnings(st_difference(z_4tjg, box))))

z_ear_mag2 <- EAR_sf %>%
  filter(regn_nm %in% c("centre", "magdalen_shallows", "northumberland_strait")) %>%
  st_union()
v2_shediac <- zone_sf(z_ear_mag2, z_nafo_otl, "Shediac-Magdalen", "Zone") %>%
  st_difference(st_geometry(v2_chaleurs)) %>% st_as_sf()

v2_pure <- NAFO_4T_sf %>%
  filter(level_2 %in% c("4TF", "4TK")) %>%
  mutate(Zone = paste(level_2, "Only")) %>%
  select(Zone)

new_regions_v2_clean <- bind_rows(v2_chaleurs, v2_northcoast, v2_4tgj,
                                  v2_shediac, v2_pure) %>%
  filter(!st_is_empty(.)) %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

st_write(new_regions_v2_clean, file.path(SPATIAL_DIR, "new_regions_v2_final.shp"),
         delete_layer = TRUE, quiet = TRUE)
cat("Written:", file.path(SPATIAL_DIR, "new_regions_v2_final.shp"), "\n")

# =============================================================================
# 3. new_gulf_final.shp  (Area - THE 4 manuscript ecoregions, CORRECTED)
# =============================================================================
clip_ear <- function(pid, name) zone_sf(EAR_sf %>% filter(PID == pid), z_4t, name, "Area")

area_chaleurs  <- clip_ear(50, "Chaleur Bay")
area_northumb  <- clip_ear(6,  "Northumberland")
area_central   <- clip_ear(3,  "Central")
area_gaspe     <- clip_ear(1,  GASPE_TO)            # ex-'northwest_estuary' in 4T

area_magdalen  <- clip_ear(5, "Magdalen Shallows") %>%
  st_difference(st_geometry(area_chaleurs)) %>%
  st_difference(st_geometry(area_northumb)) %>%
  st_as_sf()

new_gulf_clean <- bind_rows(area_chaleurs, area_northumb, area_central,
                            area_magdalen, area_gaspe) %>%
  filter(!st_is_empty(.)) %>%
  group_by(Area) %>% summarise(.groups = "drop") %>%   # merges gaspe into GASPE_TO
  st_make_valid() %>%
  st_collection_extract("POLYGON")

# ---- Checks (the script STOPS if the build is wrong) ------------------------
cat("\nAreas (km2):\n")
print(data.frame(Area = new_gulf_clean$Area,
                 km2 = round(as.numeric(st_area(new_gulf_clean)) / 1e6)))

for (i in seq_len(nrow(new_gulf_clean) - 1)) {
  for (j in (i + 1):nrow(new_gulf_clean)) {
    ov <- suppressWarnings(st_intersection(st_geometry(new_gulf_clean)[i],
                                           st_geometry(new_gulf_clean)[j]))
    a <- if (length(ov)) sum(as.numeric(st_area(ov))) / 1e6 else 0
    if (a > 1) stop(sprintf("Overlap %s x %s: %.0f km2",
                            new_gulf_clean$Area[i], new_gulf_clean$Area[j], a))
  }
}

# Nothing may spill outside 4T.
out_each <- vapply(seq_len(nrow(new_gulf_clean)), function(i) {
  d <- suppressWarnings(st_difference(st_geometry(new_gulf_clean)[i], z_4t))
  if (length(d)) sum(as.numeric(st_area(d))) / 1e6 else 0
}, numeric(1))
stopifnot(all(out_each < 200))

# Only the unsampled upper estuary (west of -66.5) may remain uncovered.
gap <- suppressWarnings(st_difference(z_4t, st_union(new_gulf_clean)))
gap_km2 <- if (length(gap)) sum(as.numeric(st_area(gap))) / 1e6 else 0
cat(sprintf("4T not covered: %.0f km2 (upper estuary only is allowed)\n", gap_km2))
if (gap_km2 > 0) {
  gap_east <- suppressWarnings(st_crop(gap, xmin = -66.5, xmax = -55,
                                       ymin = 40, ymax = 55))
  gap_east_km2 <- if (length(gap_east)) sum(as.numeric(st_area(gap_east))) / 1e6 else 0
  cat(sprintf("  of which east of -66.5 (surveyed sGSL): %.0f km2\n", gap_east_km2))
  stopifnot(gap_east_km2 < 200)
}

st_write(new_gulf_clean, file.path(SPATIAL_DIR, "new_gulf_final.shp"),
         delete_layer = TRUE, quiet = TRUE)
cat("Written:", file.path(SPATIAL_DIR, "new_gulf_final.shp"), "\n")

sf_use_s2(TRUE)
