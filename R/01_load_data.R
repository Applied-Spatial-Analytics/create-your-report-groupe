source("R/00_config.R")
library(sf)
library(terra)
library(dplyr)

sf_use_s2(FALSE)

message("Loading Yuexiu data...")

# ── Yuexiu vector layers ──────────────────────────────────────────────────────
yx_bnd       <- read_sf(YX$boundary)
yx_sub       <- read_sf(YX$subdistricts)
yx_grn       <- read_sf(YX$green)
yx_grn <- filter_green_space(yx_grn)
yx_rds       <- read_sf(YX$roads)
yx_water_poly <- read_sf(YX$water_poly)
yx_water_line <- read_sf(YX$water_line)
yx_gbif      <- read_sf(YX$gbif)
yx_viirs_pts <- read_sf(YX$viirs_pts)

# ── Yuexiu rasters ────────────────────────────────────────────────────────────
yx_pop   <- rast(YX$worldpop)
yx_ndvi  <- rast(YX$ndvi)
yx_cover <- rast(YX$worldcover)
yx_viirs <- rast(YX$viirs_rast)

message("Loading Delft data...")


# ── Delft vector layers ───────────────────────────────────────────────────────
dl_bnd   <- read_sf(DL$boundary)
dl_wijk  <- read_sf(DL$wijken)

# ── FIX: drop sliver polygons from neighbouring municipalities ──────────────
# dl_wijk as loaded contains sliver polygons from neighbouring municipalities
# (e.g. Midden-Delfland: "Schipluiden", "Delfgauw") that leak in at the
# boundary clip. Confirmed via area check: these slivers are ~10-40 m²
# (vs. hundreds of km² for genuine Delft wijken) and carry wijkcode prefixes
# outside the real Delft range (WK0503xx). Previously this filter only
# existed in script 02, AFTER dl_wijk_access.rds (used by script 04) had
# already been computed on the unfiltered data. Filtering here, immediately
# after read_sf(), ensures every downstream script inherits clean geometry.
message(sprintf("Delft wijken before gemeentecode filter: %d", nrow(dl_wijk)))
dl_wijk  <- dl_wijk[dl_wijk$gemeentecode == "GM0503", ]
message(sprintf("Delft wijken after gemeentecode filter:  %d", nrow(dl_wijk)))

dl_brt   <- read_sf(DL$buurten)
dl_inc   <- read_sf(DL$income)
dl_grn   <- read_sf(DL$green)
dl_grn <- filter_green_space(dl_grn)
dl_rds   <- read_sf(DL$roads)
dl_water <- read_sf(DL$water)
dl_gbif  <- read_sf(DL$gbif)

# Load Delft rasters
dl_pop   <- rast(DL$worldpop)
dl_ndvi  <- rast(DL$ndvi)
dl_cover <- rast(DL$worldcover)

message("Validating geometries and CRS...")

# ── Fix geometries ────────────────────────────────────────────────────────────
yx_grn        <- st_make_valid(yx_grn)
yx_water_poly <- st_make_valid(yx_water_poly)
yx_water_line <- st_make_valid(yx_water_line)
yx_sub        <- st_make_valid(yx_sub)
dl_grn        <- st_make_valid(dl_grn)
dl_water      <- st_make_valid(dl_water)
dl_wijk       <- st_make_valid(dl_wijk)

# ── Drop green patches below minimum size ─────────────────────────────────────
yx_grn <- yx_grn |>
  mutate(area_ha = as.numeric(st_area(
    st_transform(st_geometry(yx_grn), CRS_YX))) / 10000) |>
  filter(area_ha >= MIN_PATCH_HA)

dl_grn <- dl_grn |>
  mutate(area_ha = as.numeric(st_area(
    st_transform(st_geometry(dl_grn), CRS_DELFT))) / 10000) |>
  filter(area_ha >= MIN_PATCH_HA)

message(sprintf("Yuexiu green patches: %d (>= %.1f ha)", nrow(yx_grn), MIN_PATCH_HA))
message(sprintf("Delft green patches:  %d (>= %.1f ha)", nrow(dl_grn), MIN_PATCH_HA))
message(sprintf("Yuexiu subdistricts:  %d", nrow(yx_sub)))
message(sprintf("Delft wijken:         %d", nrow(dl_wijk)))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(
  list(
    yx_bnd        = yx_bnd,
    yx_sub        = yx_sub,
    yx_grn        = yx_grn,
    yx_rds        = yx_rds,
    yx_water_poly = yx_water_poly,
    yx_water_line = yx_water_line,
    yx_gbif       = yx_gbif,
    yx_viirs_pts  = yx_viirs_pts,
    yx_pop        = wrap(yx_pop),
    yx_ndvi       = wrap(yx_ndvi),
    yx_cover      = wrap(yx_cover),
    yx_viirs      = wrap(yx_viirs)
  ),
  file.path(OUT_ROOT, "yuexiu_data.rds")
)

saveRDS(
  list(
    dl_bnd       = dl_bnd,
    dl_wijk      = dl_wijk,
    dl_brt       = dl_brt,
    dl_inc       = dl_inc,
    dl_grn       = dl_grn,
    dl_rds       = dl_rds,
    dl_water     = dl_water,
    dl_gbif      = dl_gbif,
    dl_pop       = wrap(dl_pop),
    dl_ndvi      = wrap(dl_ndvi),
    dl_cover     = wrap(dl_cover)
  ),
  file.path(OUT_ROOT, "delft_data.rds")
)

message("Script 01 complete — data loaded and saved.")

