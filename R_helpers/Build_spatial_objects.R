### Spatial Pre-processing for EcoRegion Classification and Diet Data Spatial Join

#This report documents the spatial processing of the Southern Gulf of St. Lawrence ecoregions.
#The objective is to merge ecological boundaries (EAR) with administrative fishing zones (NAFO) to create a standardized spatial framework for diet analysis.


# Set global spatial engine preference
library(sf)
sf_use_s2(TRUE)

# Ensure geometries are valid before processing
EAR_sf <- st_read("data/Spatial_data/EAR_map/EAR_map.shp")
NAFO_4T_sf <- st_read("data/Spatial_data/NAFO/nafo_2014_02.shp")
NAFO_4A_sf <- st_read("data/Spatial_data/nafo_4T/nafo_2014_02.shp")

EAR_sf      <- st_make_valid(EAR_sf)
NAFO_4T_sf  <- st_make_valid(NAFO_4T_sf)

# Helper function to handle intersections and return a clean sf object
quick_intersect <- function(sf_a, sf_b, name) {
  inter <- st_intersection(sf_a, sf_b)

  # Merge fragments into a single geometry to avoid row-replacement errors
  geom_combined <- st_union(st_geometry(inter))

  # Reconstruct as a clean sf object with a single row
  result <- st_sf(EcoZone = name, geometry = geom_combined)
  return(result)
}


quick_intersect0 <- function(sf_a, sf_b, name) {
  inter <- st_intersection(sf_a, sf_b)

  # Merge fragments into a single geometry to avoid row-replacement errors
  geom_combined <- st_union(st_geometry(inter))

  # Reconstruct as a clean sf object with a single row
  result <- st_sf(Area = name, geometry = geom_combined)
  return(result)
}



### Defining Ecoregions: We define specific zones where ecological boundaries intersect with NAFO units.


# 1. Baie des Chaleurs (PID 50 + 4TM/4TN)
z_4tmn <- NAFO_4T_sf %>% filter(level_2 %in% c("4TM", "4TN")) %>% st_union()
baie_chaleurs <- quick_intersect(EAR_sf %>% filter(PID == 50), z_4tmn, "Baie des Chaleurs")

# 2. Northumberland (PID 6 + 4TH/4TL/4TG)
# Planar processing (S2 off) is safer for complex coastal intersections
sf_use_s2(FALSE)
z_4thlg <- NAFO_4T_sf %>% filter(level_2 %in% c("4TH", "4TL", "4TG")) %>% st_union()
northumberland <- quick_intersect(EAR_sf %>% filter(PID == 6), z_4thlg, "Northumberland")
sf_use_s2(TRUE)


#### Cleaning Overlaps: To ensure no overlapping geometries, we subtract high-priority zones from broader administrative units.

sf_use_s2(FALSE)

# 1. Clean 4TGJ (Subtracting Northumberland overlap)
nafo_4tgj <- NAFO_4T_sf %>% filter(level_2 %in% c("4TG", "4TJ")) %>% st_union()
nafo_4tgj_clean <- st_difference(nafo_4tgj, northumberland) %>%
  st_as_sf() %>%
  mutate(EcoZone = "4TGJ clean") %>%
  rename(geometry = x)

# 2. Clean Shediac-Magdalen (Subtracting Baie des Chaleurs overlap)
# Note: z_ear_mag and z_nafo_otl assumed pre-defined in environment

z_nafo_otl <- NAFO_4T_sf %>%
  filter(level_2 %in% c("4TO", "4TN", "4TL")) %>%
  st_union()

z_ear_mag <- EAR_sf %>%
  filter(regn_nm %in% c("centre", "magdalen_shallows")) %>%
  st_union()

shediac_mag_brut <- quick_intersect(z_ear_mag, z_nafo_otl, "Shediac-Magdalen")
shediac_clean <- st_difference(shediac_mag_brut, st_geometry(baie_chaleurs)) %>%
  st_as_sf() %>%
  mutate(EcoZone = "Shediac-Magdalen")

sf_use_s2(TRUE)

# 3. Pure NAFO Zones (4TF, 4TK)
zones_pure <- NAFO_4T_sf %>%
  filter(level_2 %in% c("4TF", "4TK")) %>%
  mutate(EcoZone = paste(level_2, "Only")) %>%
  select(EcoZone)



### Final Assembly and Export: We bind all processed layers into a single master ecoregion layer.


new_ecoregions <- bind_rows(
  baie_chaleurs %>% select(EcoZone),
  northumberland %>% select(EcoZone),
  nafo_4tgj_clean %>% select(EcoZone),
  shediac_clean %>% select(EcoZone),
  zones_pure %>% select(EcoZone)
) %>%
  filter(!st_is_empty(.))

# Final geometry cleanup to ensure polygon-only output (required for Shapefiles)
new_ecoregions_clean <- new_ecoregions %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

# Export to Shapefile
st_write(new_ecoregions_clean, "data/Spatial_data/new_ecoregions_final.shp", delete_layer = TRUE)



#### Intersection Zones for a 2nd Ecozone proposal: We define specific zones where ecological boundaries intersect with NAFO units.

# --- 1. SETTINGS & UTILITIES ---
sf_use_s2(TRUE)

# Updated utility function to ensure consistent column naming ('Zone')
quick_intersect_v2 <- function(sf_a, sf_b, name) {
  inter <- st_intersection(sf_a, sf_b)
  geom_combined <- st_union(st_geometry(inter))
  result <- st_sf(Zone = name, geometry = geom_combined)
  return(result)
}

# --- 2. CUSTOM POLYGON DEFINITION ---
# Defining a box from Longitude -65 to -61 and Latitude 45.5 to 46.47
coords <- matrix(c(
  -64, 45.5,
  -61, 45.5,
  -61, 46.47,
  -64, 46.47,
  -64, 45.5
), ncol = 2, byrow = TRUE)

new_poly_zone <- st_polygon(list(coords)) %>%
  st_sfc(crs = 4326) %>%
  st_sf(Zone = "New Custom Zone", geometry = .)

# --- 3. COMPONENT GENERATION ---

# A. Baie des Chaleurs (Intersection of EAR PID 50 and NAFO 4TM/N)
z_4tmn_union <- NAFO_4T_sf %>% filter(level_2 %in% c("4TM", "4TN")) %>% st_union()
baie_chaleurs <- quick_intersect_v2(EAR_sf %>% filter(PID == 50), z_4tmn_union, "Baie des Chaleurs")

# B. North Coast (Custom Box intersected with NAFO 4TH/G)
sf_use_s2(FALSE) # Switch off S2 for planar intersection with box
z_4thg_union <- NAFO_4T_sf %>% filter(level_2 %in% c("4TH", "4TG")) %>% st_union()
nort_coast <- quick_intersect_v2(z_4thg_union, new_poly_zone, "North Coast")
sf_use_s2(TRUE)

# C. 4TGJ Coast (NAFO 4TJ/G minus the custom box area)
z_4tjg_union <- NAFO_4T_sf %>% filter(level_2 %in% c("4TJ", "4TG")) %>% st_union()
JG4T_coast <- st_difference(z_4tjg_union, st_geometry(new_poly_zone)) %>%
  st_as_sf() %>%
  mutate(Zone = "4TGJ coast") %>%
  # Ensure geometry column name is standardized
  st_set_geometry("geometry")

# D. Shediac-Magdalen (Cleaned of Baie des Chaleurs overlap)
sf_use_s2(FALSE)
z_nafo_otl_union <- NAFO_4T_sf %>% filter(level_2 %in% c("4TO", "4TN", "4TL")) %>% st_union()
z_ear_mag_union  <- EAR_sf %>% filter(regn_nm %in% c("centre", "magdalen_shallows","northumberland_strait")) %>% st_union()
shediac_brut <- quick_intersect_v2(z_ear_mag_union, z_nafo_otl_union, "Shediac-Magdalen")

shediac_cl <- st_difference(shediac_brut, st_geometry(baie_chaleurs)) %>%
  st_as_sf() %>%
  mutate(Zone = "Shediac-Magdalen")
sf_use_s2(TRUE)

# E. Pure NAFO Zones
zones_pure2 <- NAFO_4T_sf %>%
  filter(level_2 %in% c("4TF", "4TK")) %>%
  mutate(Zone = paste(level_2, "Only")) %>%
  select(Zone)

# --- 4. FINAL ASSEMBLY ---

# Standardizing all components to have only the 'Zone' column
new_regions_v2 <- bind_rows(
  baie_chaleurs  %>% select(Zone),
  nort_coast     %>% select(Zone),
  JG4T_coast     %>% select(Zone),
  shediac_cl     %>% select(Zone),
  zones_pure2    %>% select(Zone)
) %>%
  filter(!st_is_empty(.))

# Clean geometries for Shapefile export
new_regions_v2_clean <- new_regions_v2 %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

# --- 5. EXPORT ---
st_write(new_regions_v2_clean, "data/Spatial_data/new_regions_v2_final.shp", delete_layer = TRUE)



# ==============================================================================
# SCRIPT: Spatial Boundary Definition and Multi-Period Mapping
# ==============================================================================

library(sf)
library(dplyr)
library(ggplot2)

# --- 1. Base Boundary Definition ---
# Combine NAFO Level 4 divisions to create the broad study area boundary
z_nafo4_otl <- NAFO_4A_sf %>%
  filter(level_0 == 4) %>%
  st_union()

# --- 2. Geometric Intersection & Regional Cleaning ---
# Note: We toggle s2 geometry engine off for certain planar intersections
# to avoid topology errors.

sf_use_s2(FALSE)

# Define Magdalen Shallows ecoregion
z_ear_mag        <- EAR_sf %>% filter(regn_nm == "magdalen_shallows") %>% st_union()
shediac_mag_brut <- quick_intersect0(z_ear_mag, z_nafo4_otl, "Magdalen Shallows")
shediac_clean <- st_difference(shediac_mag_brut, st_geometry(baie_chaleurs)) %>%
  st_as_sf()

# Define Central region
z_ear_cent        <- EAR_sf %>% filter(regn_nm == "centre") %>% st_union()
shediac_cent_brut <- quick_intersect0(z_ear_cent, z_nafo4_otl, "Central")

# Define Baie des Chaleurs (Intersection of EAR PID 50 and NAFO 4TM/N)
sf_use_s2(TRUE)
z_4tmn         <- NAFO_4T_sf %>% filter(level_2 %in% c("4TM", "4TN")) %>% st_union()
baie_chaleurs  <- quick_intersect0(EAR_sf %>% filter(PID == 50), z_4tmn, "Baie des Chaleurs")

# Define Northumberland Strait (Intersection of EAR PID 6 and NAFO 4TH/L/G)
sf_use_s2(FALSE)
z_4thlg        <- NAFO_4T_sf %>% filter(level_2 %in% c("4TH", "4TL", "4TG")) %>% st_union()
northumberland <- quick_intersect0(EAR_sf %>% filter(PID == 6), z_4thlg, "Northumberland")

# --- 3. Data Integration & Topology Cleanup ---
# Merge all custom regions into a single sf object
new_gulf <- bind_rows(
  baie_chaleurs  %>% select(Area),
  northumberland %>% select(Area),
  shediac_cent_brut %>% select(Area),
  shediac_clean  %>% select(Area)
) %>%
  filter(!st_is_empty(.))

# Ensure valid geometry and extract only polygons (removes lines/points from edges)
new_gulf_clean <- new_gulf %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

# Export refined spatial layers
st_write(new_gulf_clean, "data/Spatial_data/new_gulf_final.shp", delete_layer = TRUE)



