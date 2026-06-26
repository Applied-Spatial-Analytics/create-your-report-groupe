source("R/00_config.R")
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

sf_use_s2(FALSE)

d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

# Unwrap rasters (safe: only unwrap if still PackedSpatRaster)
safe_rast <- function(x) if (inherits(x, "PackedSpatRaster")) rast(x) else x
d$yx_ndvi   <- safe_rast(d$yx_ndvi)
d$yx_cover  <- safe_rast(d$yx_cover)
dl$dl_ndvi  <- safe_rast(dl$dl_ndvi)
dl$dl_cover <- safe_rast(dl$dl_cover)

# Pre-clean admin polygons that are known to have topology issues
d$yx_sub  <- st_make_valid(d$yx_sub)
d$yx_sub  <- st_collection_extract(d$yx_sub, "POLYGON")

# ── Helper: scrub an sf object of every geometry problem ─────────────────────
scrub_sf <- function(x, label = "") {
  n0   <- nrow(x)
  keep <- seq_len(n0)

  ok <- !st_is_empty(x)
  x  <- x[ok, ]; keep <- keep[ok]

  x <- st_make_valid(x)
  # Extra precision repair — fixes TopologyException from BGT/OSM data
  if (st_is_longlat(x)) {
    crs_orig <- st_crs(x)
    x <- st_transform(x, 3857)
    x <- st_buffer(x, 0)
    x <- st_transform(x, crs_orig)
  } else {
    x <- st_buffer(x, 0)
  }

  ok <- !st_is_empty(x)
  x  <- x[ok, ]; keep <- keep[ok]

  ok <- st_is_valid(x)
  ok[is.na(ok)] <- FALSE
  x  <- x[ok, ]; keep <- keep[ok]

  geom_types <- as.character(st_geometry_type(x))
  is_coll    <- geom_types == "GEOMETRYCOLLECTION"
  if (any(is_coll)) {
    st_geometry(x)[is_coll] <- lapply(st_geometry(x)[is_coll], function(g) {
      parts <- st_collection_extract(st_sfc(g, crs = st_crs(x)), "POLYGON")
      if (length(parts) == 0) st_geometrycollection() else st_union(parts)[[1]]
    })
    x <- st_set_crs(x, st_crs(x))
  }

  has_na_coords <- function(g) {
    coords <- tryCatch(st_coordinates(g), error = function(e) matrix(NA_real_, 1, 2))
    anyNA(coords)
  }
  ok <- !vapply(st_geometry(x), has_na_coords, logical(1))
  x  <- x[ok, ]; keep <- keep[ok]

  if (any(st_geometry_type(x) %in% c("POLYGON","MULTIPOLYGON","GEOMETRYCOLLECTION"))) {
    ok <- suppressWarnings(as.numeric(st_area(x))) > 0
    ok[is.na(ok)] <- FALSE
    x  <- x[ok, ]; keep <- keep[ok]
  }

  dropped <- n0 - nrow(x)
  if (dropped > 0)
    message(sprintf("  scrub_sf [%s]: dropped %d / %d geometries", label, dropped, n0))

  list(sf = x, kept = keep)
}

# ── 2A. Green space typology ──────────────────────────────────────────────────
classify_green_type <- function(green_sf) {
  geom_col <- attr(green_sf, "sf_column")
  cols     <- tolower(names(green_sf))

  is_bgt <- "fysiek_voorkomen" %in% cols

  #put a filter on which green spaces to count as publicly accessible urban green space

  if (is_bgt) {
    green_sf |>
      mutate(
        tag_raw = dplyr::coalesce(
          as.character(plus_fysiek_voorkomen),
          as.character(fysiek_voorkomen),
          "unknown"
        ),
        green_type = case_when(
          grepl("loofbos|gemengd bos|houtwal|bosplantsoen", tag_raw, ignore.case = TRUE) ~ "Forest / Woodland",
          grepl("grasland|gras- en kruid|braakliggend",      tag_raw, ignore.case = TRUE) ~ "Grass / Meadow",
          grepl("groenvoorziening|planten|heesters|struiken", tag_raw, ignore.case = TRUE) ~ "Park / Recreation",
          grepl("rietland",                                  tag_raw, ignore.case = TRUE) ~ "Nature Reserve / Scrub",
          grepl("bouwland|fruitteelt|hoogstam",              tag_raw, ignore.case = TRUE) ~ "Allotment / Agriculture",
          TRUE                                                                            ~ "Other / Unclassified"
        )
      )
  } else {
    green_sf |>
      mutate(
        across(-all_of(geom_col), ~ as.character(.x)),
        tag_raw = dplyr::coalesce(
          if ("leisure"    %in% cols) .data[["leisure"]]    else NA_character_,
          if ("landuse"    %in% cols) .data[["landuse"]]    else NA_character_,
          if ("natural"    %in% cols) .data[["natural"]]    else NA_character_,
          if ("other_tags" %in% cols) .data[["other_tags"]] else NA_character_,
          "unknown"
        ),
        green_type = case_when(
          grepl("park|recreation_ground|garden|pleasure_ground", tag_raw, ignore.case = TRUE) ~ "Park / Recreation",
          grepl("forest|wood|tree_row",                          tag_raw, ignore.case = TRUE) ~ "Forest / Woodland",
          grepl("grass|meadow|village_green|common",             tag_raw, ignore.case = TRUE) ~ "Grass / Meadow",
          grepl("cemetery|grave",                                tag_raw, ignore.case = TRUE) ~ "Cemetery",
          grepl("allotment|farm|farmland",                       tag_raw, ignore.case = TRUE) ~ "Allotment / Agriculture",
          grepl("nature_reserve|wetland|scrub|heath",            tag_raw, ignore.case = TRUE) ~ "Nature Reserve / Scrub",
          grepl("pitch|track|sports_centre",                     tag_raw, ignore.case = TRUE) ~ "Sports Facility",
          TRUE                                                                                ~ "Other / Unclassified"
        )
      )
    recreational_green_types <- c(
      "Park / Recreation",
      "Forest / Woodland",
      "Grass / Meadow",
      "Nature Reserve / Scrub"
    )

    green_sf <- green_sf |>
      filter(green_type %in% recreational_green_types)
  }
}

message("Classifying Yuexiu green typology...")
yx_grn_typed <- classify_green_type(d$yx_grn)
yx_grn_typed <- yx_grn_typed |>
  mutate(area_ha = as.numeric(st_area(st_transform(yx_grn_typed, CRS_YX))) / 10000)

message("Classifying Delft green typology...")
dl_grn_typed <- classify_green_type(dl$dl_grn)
dl_grn_typed <- dl_grn_typed |>
  mutate(area_ha = as.numeric(st_area(st_transform(dl_grn_typed, CRS_DELFT))) / 10000)

# Summary tables
yx_type_summary <- yx_grn_typed |>
  st_drop_geometry() |>
  group_by(green_type) |>
  summarise(n_patches = n(),
            total_ha  = sum(area_ha, na.rm = TRUE),
            mean_ha   = mean(area_ha, na.rm = TRUE),
            .groups   = "drop") |>
  mutate(pct_area = total_ha / sum(total_ha) * 100, city = "Yuexiu")

dl_type_summary <- dl_grn_typed |>
  st_drop_geometry() |>
  group_by(green_type) |>
  summarise(n_patches = n(),
            total_ha  = sum(area_ha, na.rm = TRUE),
            mean_ha   = mean(area_ha, na.rm = TRUE),
            .groups   = "drop") |>
  mutate(pct_area = total_ha / sum(total_ha) * 100, city = "Delft")

type_summary <- bind_rows(yx_type_summary, dl_type_summary)
message("Green typology summary:")
print(type_summary)

# ── 2B. NDVI stats per green patch ───────────────────────────────────────────
safe_ndvi_patches <- function(grn_typed, ndvi_rast, label) {
  message(sprintf("Extracting NDVI per patch — %s...", label))

  proj    <- st_transform(grn_typed, st_crs(ndvi_rast))
  sc      <- scrub_sf(proj, label)
  proj_ok <- sc$sf

  grn_typed$ndvi_mean <- NA_real_
  grn_typed$ndvi_sd   <- NA_real_

  if (nrow(proj_ok) == 0) {
    warning(sprintf("safe_ndvi_patches [%s]: no valid geometries.", label))
    return(grn_typed)
  }

  sv  <- tryCatch(vect(proj_ok), error = function(e) NULL)
  if (is.null(sv)) return(grn_typed)
  raw <- tryCatch(terra::extract(ndvi_rast, sv), error = function(e) NULL)
  if (is.null(raw)) return(grn_typed)

  val_col <- names(raw)[2]
  agg <- raw |>
    group_by(ID) |>
    summarise(ndvi_mean = mean(.data[[val_col]], na.rm = TRUE),
              ndvi_sd   = sd(.data[[val_col]],   na.rm = TRUE),
              .groups   = "drop")

  ndvi_mean <- rep(NA_real_, nrow(proj_ok))
  ndvi_sd   <- rep(NA_real_, nrow(proj_ok))
  ndvi_mean[agg$ID] <- agg$ndvi_mean
  ndvi_sd[agg$ID]   <- agg$ndvi_sd

  grn_typed$ndvi_mean[sc$kept] <- ndvi_mean
  grn_typed$ndvi_sd[sc$kept]   <- ndvi_sd
  grn_typed
}

yx_grn_typed <- safe_ndvi_patches(yx_grn_typed, d$yx_ndvi,  "Yuexiu")
dl_grn_typed <- safe_ndvi_patches(dl_grn_typed, dl$dl_ndvi, "Delft")

# ── 2C. NDVI zonal stats per admin unit ──────────────────────────────────────
calc_ndvi_zonal <- function(admin_sf, ndvi_rast, local_crs, label = "") {
  result <- admin_sf |>
    mutate(ndvi_mean = NA_real_, ndvi_p25 = NA_real_, ndvi_p75 = NA_real_)

  admin_proj <- st_transform(admin_sf, st_crs(ndvi_rast))
  sc         <- scrub_sf(admin_proj, label)
  admin_ok   <- sc$sf

  if (nrow(admin_ok) == 0) return(result)

  not_empty <- !vapply(st_geometry(admin_ok), st_is_empty, logical(1))
  sc$kept   <- sc$kept[not_empty]
  admin_ok  <- admin_ok[not_empty, ]
  if (nrow(admin_ok) == 0) return(result)

  admin_sv <- tryCatch(vect(admin_ok), error = function(e) NULL)
  if (is.null(admin_sv)) return(result)

  raw <- tryCatch(terra::extract(ndvi_rast, admin_sv), error = function(e) NULL)
  if (is.null(raw)) return(result)

  val_col <- names(raw)[2]
  agg <- raw |>
    group_by(ID) |>
    summarise(ndvi_mean = mean(.data[[val_col]],           na.rm = TRUE),
              ndvi_p25  = quantile(.data[[val_col]], 0.25, na.rm = TRUE),
              ndvi_p75  = quantile(.data[[val_col]], 0.75, na.rm = TRUE),
              .groups   = "drop")

  ndvi_mean <- rep(NA_real_, nrow(admin_ok))
  ndvi_p25  <- rep(NA_real_, nrow(admin_ok))
  ndvi_p75  <- rep(NA_real_, nrow(admin_ok))
  ndvi_mean[agg$ID] <- agg$ndvi_mean
  ndvi_p25[agg$ID]  <- agg$ndvi_p25
  ndvi_p75[agg$ID]  <- agg$ndvi_p75

  result$ndvi_mean[sc$kept] <- ndvi_mean
  result$ndvi_p25[sc$kept]  <- ndvi_p25
  result$ndvi_p75[sc$kept]  <- ndvi_p75
  result
}

message("NDVI zonal stats — Yuexiu subdistricts...")
yx_sub_ndvi <- calc_ndvi_zonal(d$yx_sub,    d$yx_ndvi,  CRS_YX,    "Yuexiu-sub")

message("NDVI zonal stats — Delft wijken...")
dl_wijk_ndvi <- calc_ndvi_zonal(dl$dl_wijk, dl$dl_ndvi, CRS_DELFT, "Delft-wijk")

# ── 2D. Species observation density (GBIF) ───────────────────────────────────
count_pts_in_polys <- function(green_sf, gbif_sf, local_crs) {
  grn_m  <- st_transform(green_sf, local_crs)
  gbif_m <- st_transform(gbif_sf,  local_crs)
  gbif_m <- gbif_m[!st_is_empty(gbif_m) & !is.na(st_is_valid(gbif_m)) & st_is_valid(gbif_m), ]

  idx <- st_intersects(grn_m, gbif_m)
  grn_m |>
    mutate(n_gbif_obs      = lengths(idx),
           gbif_obs_per_ha = n_gbif_obs / pmax(area_ha, 0.01))
}

message("Species density — Yuexiu...")
yx_grn_bio <- count_pts_in_polys(yx_grn_typed, d$yx_gbif,  CRS_YX)

message("Species density — Delft...")
dl_grn_bio <- count_pts_in_polys(dl_grn_typed, dl$dl_gbif, CRS_DELFT)

# ── 2E. Blue-green ratio / Environmental Justice Metrics ──────────────────────
calc_blue_green_ratio <- function(admin_sf, green_sf, water_sf,
                                  water_line_sf  = NULL,
                                  local_crs,
                                  water_line_buf = WATER_LINE_BUFFER,
                                  label          = "") {

  admin_m <- scrub_sf(st_transform(admin_sf, local_crs), paste0(label, "-bg-admin"))$sf
  green_m <- scrub_sf(st_transform(green_sf, local_crs), paste0(label, "-bg-green"))$sf
  water_m <- scrub_sf(st_transform(water_sf, local_crs), paste0(label, "-bg-water"))$sf

  if (!is.null(water_line_sf)) {
    wl_poly <- st_buffer(st_transform(water_line_sf, local_crs), water_line_buf)
    wl_poly <- scrub_sf(wl_poly, paste0(label, "-bg-waterline"))$sf
    water_m <- rbind(
      water_m[, intersect(names(water_m), names(wl_poly))],
      wl_poly[, intersect(names(water_m), names(wl_poly))]
    )
  }

  # Add admin row id
  admin_m$.adm_id <- seq_len(nrow(admin_m))

  area_in_admin <- function(feat_m, admin_m, feat_label) {
    areas <- rep(0, nrow(admin_m))
    feat_m <- st_make_valid(feat_m)
    if (st_is_longlat(feat_m)) {
      crs_orig <- st_crs(feat_m)
      feat_m <- st_transform(feat_m, 3857)
      feat_m <- st_buffer(feat_m, 0)
      feat_m <- st_transform(feat_m, crs_orig)
    } else {
      feat_m <- st_buffer(feat_m, 0)
    }
    feat_geom <- st_as_sf(st_geometry(feat_m))

    inter <- tryCatch(
      suppressWarnings(st_intersection(
        feat_geom,
        admin_m[, c(".adm_id", attr(admin_m, "sf_column"))]
      )),
      error = function(e) { warning(feat_label, ": ", conditionMessage(e)); NULL }
    )
    if (is.null(inter) || nrow(inter) == 0) return(areas)

    inter <- inter[!st_is_empty(inter), ]
    inter <- inter[st_geometry_type(inter) %in%
                     c("POLYGON","MULTIPOLYGON","GEOMETRYCOLLECTION"), ]
    if (nrow(inter) == 0) return(areas)

    inter$area_m2 <- as.numeric(st_area(inter))
    agg <- tapply(inter$area_m2, inter$.adm_id, sum, na.rm = TRUE)
    areas[as.integer(names(agg))] <- as.numeric(agg)
    areas
  }

  message("   computing green areas...")
  green_areas <- area_in_admin(green_m, admin_m, "green")
  message("   computing water areas...")
  water_areas <- area_in_admin(water_m, admin_m, "water")

  out <- admin_sf |>
    mutate(green_m2 = NA_real_, water_m2 = NA_real_,
           green_pc_m2 = NA_real_, water_pc_m2 = NA_real_,
           pop_density = NA_real_, green_pressure_idx = NA_real_,
           blue_green_balance = NA_real_)

  sc_idx <- scrub_sf(st_transform(admin_sf, local_crs), paste0(label, "-bg-admin-idx"))$kept
  idx    <- if (length(sc_idx) == length(green_areas)) sc_idx else seq_along(green_areas)

  out$green_m2[idx] <- green_areas
  out$water_m2[idx] <- water_areas

  # SAFELY HANDLE MISSING POPULATION COLUMN BEFORE SLICING
  if (!"pop_count" %in% names(admin_sf)) {
    warning(sprintf("pop_count column missing in study context [%s]. Falling back to area metrics (Pop = 1).", label))
    pop_vec_orig <- rep(1, nrow(admin_sf))
  } else {
    pop_vec_orig <- admin_sf$pop_count
  }

  # Clean alignment logic with execution geometry subsets
  pop_vec <- pop_vec_orig[idx]
  pop_vec[is.na(pop_vec)] <- 1

  green_pc <- green_areas / pmax(pop_vec, 1)
  water_pc <- water_areas / pmax(pop_vec, 1)

  pop_density <- pop_vec / (as.numeric(st_area(st_transform(admin_sf[idx, ], local_crs))) / 10000)

  out$green_pc_m2[idx]        <- green_pc
  out$water_pc_m2[idx]        <- water_pc
  out$pop_density[idx]        <- pop_density
  out$green_pressure_idx[idx] <- pop_density / (green_pc + 0.01)
  out$blue_green_balance[idx] <- log1p(water_pc) - log1p(green_pc)

  return(out)
}

# Deduplicate wijken by wijkcode to avoid phantom rows from CBS join
dl$dl_wijk <- dl$dl_wijk[!duplicated(dl$dl_wijk$wijkcode), ]
dl$dl_wijk <- dl$dl_wijk[!is.na(dl$dl_wijk$wijknaam) & dl$dl_wijk$wijknaam != "", ]
dl$dl_wijk <- dl$dl_wijk[dl$dl_wijk$gemeentecode == "GM0503", ]

message("Blue-green ratio — Yuexiu subdistricts...")
yx_sub_bg <- calc_blue_green_ratio(
  d$yx_sub, d$yx_grn, d$yx_water_poly,
  water_line_sf = d$yx_water_line,
  local_crs     = CRS_YX,
  label         = "Yuexiu"
)

message("Blue-green ratio — Delft wijken...")
dl_wijk_bg <- calc_blue_green_ratio(
  dl$dl_wijk, dl$dl_grn, dl$dl_water,
  water_line_sf = NULL,
  local_crs     = CRS_DELFT,
  label         = "Delft"
)

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(yx_grn_bio,   file.path(OUT_ROOT, "yx_grn_bio.rds"))
saveRDS(dl_grn_bio,   file.path(OUT_ROOT, "dl_grn_bio.rds"))
saveRDS(yx_sub_ndvi,  file.path(OUT_ROOT, "yx_sub_ndvi.rds"))
saveRDS(dl_wijk_ndvi, file.path(OUT_ROOT, "dl_wijk_ndvi.rds"))
saveRDS(yx_sub_bg,    file.path(OUT_ROOT, "yx_sub_bg.rds"))
saveRDS(dl_wijk_bg,   file.path(OUT_ROOT, "dl_wijk_bg.rds"))
saveRDS(type_summary, file.path(OUT_ROOT, "type_summary.rds"))

# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 2A: Fixed assignment pipeline
# Figure 2AA - typology quantity
type_summary <- type_summary |>
  dplyr::filter(green_type != "Other / Unclassified")
type_summary$green_type <- factor(type_summary$green_type, levels = c(
    "Forest / Woodland",
    "Park / Recreation",
    "Nature Reserve / Scrub",
    "Grass / Meadow"))

p_type <- ggplot(type_summary, aes(x = city, y = pct_area, fill = green_type)) +
geom_col() +
  scale_fill_manual(values = c("Forest / Woodland" = COLORS$green_dark, "Nature Reserve / Scrub" = COLORS$green_mid, "Park / Recreation" = COLORS$green_light, "Grass / Meadow" = COLORS$beige ), name = "Typology") +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_y_continuous(labels = scales::label_number(suffix = "%")) +
  theme_minimal(base_size = 12) +
  labs(title    = "Green Space Composition by Typology", subtitle = "OSM/BGT tag-based classification — % of total green area", x = NULL, y = "% of total green area") +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(OUT_ROOT, "fig_green_typology.png"), p_type, width = 9, height = 6, dpi = 300)

# Figure 2AB - NDVI

# transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_ndvi, CRS_YX)
dl_plot <- st_transform(dl_wijk_ndvi, 28992)

# calculate real maps widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

p_ndvi_yx <- ggplot(yx_plot) +
  geom_sf(aes(fill = ndvi_mean), linewidth = 0.25, color = "white") +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  scale_fill_gradientn( colours = pal_green, limits = c(0.1, 0.8), na.value = COLORS$grey85, name = "Mean NDVI") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs( title = "Yuexiu", subtitle = "Sentinel-2 NDVI imagery")

p_ndvi_dl <- ggplot(dl_plot) +
  geom_sf(aes(fill = ndvi_mean), linewidth = 0.25, color = "white") +
  scale_fill_gradientn( colours = pal_green, limits = c(0.1, 0.8), na.value = COLORS$grey85, name = "Sentinel-2 NDVI imagery", guide = "none") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "Sentinel-2 NDVI imagery")

fig_ndvi <- p_ndvi_yx + p_ndvi_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation(title = "Mean NDVI", theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_ndvi

ggsave(file.path(OUT_ROOT, "fig_ndvi.png"), fig_ndvi, width = 14, height = 6, dpi = 300)

# Figure 2D: NDVI violin
ndvi_df <- bind_rows(
  yx_grn_bio |> st_drop_geometry() |> transmute(city = "Yuexiu", ndvi_mean, green_type),
  dl_grn_bio |> st_drop_geometry() |> transmute(city = "Delft",  ndvi_mean, green_type)
) |> filter(!is.na(ndvi_mean))

p_ndvi_violin <- ggplot(ndvi_df, aes(x = green_type, y = ndvi_mean, fill = city)) +
  geom_violin(alpha = 0.7, position = position_dodge(0.8), draw_quantiles = 0.5) +
  scale_fill_manual(values = c("Yuexiu" = COLORS$orange, "Delft" = COLORS$blue)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title    = "NDVI Distribution by Green Typology",
       subtitle = "Patch-level mean NDVI (vegetation health proxy)",
       x = "Green Type", y = "Mean NDVI", fill = "City") +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(OUT_ROOT, "fig_ndvi_violin.png"),
  p_ndvi_violin,        # was p_type
  width = 12, height = 6, dpi = 300)

# Figure 2E: GBIF species density
bio_df <- bind_rows(
  yx_grn_bio |> st_drop_geometry() |>
    transmute(city = "Yuexiu", green_type, gbif_obs_per_ha, area_ha),
  dl_grn_bio |> st_drop_geometry() |>
    transmute(city = "Delft",  green_type, gbif_obs_per_ha, area_ha)
) |> filter(!is.na(gbif_obs_per_ha), gbif_obs_per_ha > 0, area_ha > 0)

p_gbif <- ggplot(bio_df, aes(x = area_ha, y = gbif_obs_per_ha,
                             colour = green_type, shape = city)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_colour_manual(
    values = c(
      "Forest / Woodland" = COLORS$green_dark,
      "Grass / Meadow" = COLORS$beige,
      "Nature Reserve / Scrub" = COLORS$green_mid,
      "Park / Recreation" = COLORS$green_light)) +
  theme_minimal(base_size = 11) +
  labs(title    = "Biodiversity Observation Density vs. Patch Size",
       subtitle = "GBIF occurrences per ha — log-log scale",
       x = "Patch area (ha, log)", y = "Obs. per ha (log)", shape = "City") +
    theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_ROOT, "fig_gbif_density.png"),
  p_gbif,               # was p_type
  width = 10, height = 6, dpi = 300)


# Figure 2F: Blue-Green Balance

# Transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_bg, CRS_YX)
dl_plot <- st_transform(dl_wijk_bg, 28992)

# Calculate real map widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

# Use shared limits so both maps are comparable
bg_limits <- range(c(yx_plot$blue_green_balance, dl_plot$blue_green_balance), na.rm = TRUE)

bg_values <- scales::rescale(
  c(bg_limits[1], 0, bg_limits[2]), from = bg_limits)

# recreate blue_green_balance before plotting
yx_sub_bg <- yx_sub_bg %>%
  mutate(blue_green_balance = log1p(water_pc_m2) - log1p(green_pc_m2))

dl_wijk_bg <- dl_wijk_bg %>%
  mutate(blue_green_balance = log1p(water_pc_m2) - log1p(green_pc_m2))

p_bg_yx <- ggplot(yx_plot) +
  geom_sf(aes(fill = blue_green_balance), linewidth = 0.25, color = "white") +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  scale_fill_gradientn(colours = pal_green_blue, values = bg_values, limits = bg_limits, na.value = COLORS$grey85, name = "BG Balance\n(log scale)") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs( title = "Yuexiu", subtitle = "Source: green, water and population datasets")

p_bg_dl <- ggplot(dl_plot) +
  geom_sf(aes(fill = blue_green_balance), linewidth = 0.25, color = "white") +
  scale_fill_gradientn( colours = pal_green_blue, values = bg_values, limits = bg_limits, na.value = COLORS$grey85, name = "BG Balance\n(log scale)", guide = "none"
  ) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "Source: green, water and population datasets")

fig_blue_green_ratio <- p_bg_yx + p_bg_dl +
  plot_layout( widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation( title = "Blue-Green Balance", theme = theme( plot.title = element_text(face = "bold", size = 14))) &
  theme( legend.position = "right")

fig_blue_green_ratio

ggsave(file.path(OUT_ROOT, "fig_blue_green_ratio.png"), fig_blue_green_ratio, width = 14, height = 6, dpi = 300)

message("Script 03 complete — typology & biodiversity figures saved.")
