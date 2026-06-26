source("R/00_config.R")
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

sf_use_s2(FALSE)

d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

# Unwrapping raster
d$yx_pop   <- rast(d$yx_pop)
d$yx_ndvi  <- rast(d$yx_ndvi)
d$yx_cover <- rast(d$yx_cover)
d$yx_viirs <- rast(d$yx_viirs)

dl$dl_pop   <- rast(dl$dl_pop)
dl$dl_ndvi  <- rast(dl$dl_ndvi)
dl$dl_cover <- rast(dl$dl_cover)

# ── 1A. Green space per capita per subdistrict ────────────────────────────────
calc_green_per_capita <- function(admin_sf, green_sf, pop_rast, local_crs) {
  # 1. Ensure geometries are valid and non-empty
  admin_sf <- admin_sf[!st_is_empty(admin_sf), ]
  admin_sf <- st_make_valid(admin_sf)

  admin_m <- st_transform(admin_sf, local_crs)
  green_m <- st_transform(green_sf, local_crs)
  green_union <- st_union(green_m)

  results <- lapply(seq_len(nrow(admin_m)), function(i) {
    unit    <- admin_m[i, ]
    clipped <- tryCatch(
      st_intersection(green_union, st_geometry(unit)),
      error = function(e) st_sfc(crs = local_crs)
    )
    area_m2 <- if (length(clipped) > 0 && !st_is_empty(clipped)) {
      as.numeric(st_area(clipped))
    } else { 0 }
    area_m2
  })

  # 2. FIX: Project the extraction boundary to the RASTER's native CRS (EPSG:4326)
  admin_for_extract <- st_transform(admin_sf, st_crs(pop_rast))

  admin_sf |>
    mutate(
      green_area_m2 = unlist(results),
      pop_count     = exact_extract(pop_rast, admin_for_extract, "sum"),
      pop_count = ifelse(is.na(pop_count), NA_real_, pop_count),
      green_pc_m2   = green_area_m2 / pop_count
    )
}

message("Calculating Yuexiu green per capita...")
yx_sub_access <- calc_green_per_capita(d$yx_sub, d$yx_grn, d$yx_pop, CRS_YX)

message("Calculating Delft green per capita...")
dl_wijk_access <- calc_green_per_capita(dl$dl_wijk, dl$dl_grn, dl$dl_pop, CRS_DELFT)


# ── 1B. % population within 300m and 500m of any green space ─────────────────
pct_pop_within_buffer <- function(green_sf, pop_rast, boundary_sf, dist_m, local_crs) {
  # Perform distance buffering in the metric projected CRS
  green_m <- st_transform(green_sf, local_crs)
  bnd_m   <- st_transform(boundary_sf, local_crs)
  buf     <- st_union(st_buffer(green_m, dist_m))

  # Clip inside the local metric projection framework
  clipped_m <- tryCatch(st_intersection(bnd_m, buf), error = function(e) NULL)

  # FIX: Transform vectors to match the raster's native CRS right before exact_extract
  bnd_for_extract <- st_transform(boundary_sf, st_crs(pop_rast))

  pop_total  <- exact_extract(pop_rast, bnd_for_extract, "sum")

  pop_within <- if (!is.null(clipped_m) && nrow(clipped_m) > 0) {
    # Transform clipped boundary to raster CRS
    clipped_for_extract <- st_transform(clipped_m, st_crs(pop_rast))
    exact_extract(pop_rast, clipped_for_extract, "sum")
  } else { 0 }

  round(sum(pop_within, na.rm = TRUE) / sum(pop_total, na.rm = TRUE) * 100, 1)
}

message("Calculating buffer coverage...")
yx_300 <- pct_pop_within_buffer(d$yx_grn,  d$yx_pop, d$yx_bnd, 300, CRS_YX)
yx_500 <- pct_pop_within_buffer(d$yx_grn,  d$yx_pop, d$yx_bnd, 500, CRS_YX)
dl_300 <- pct_pop_within_buffer(dl$dl_grn, dl$dl_pop, dl$dl_bnd, 300, CRS_DELFT)
dl_500 <- pct_pop_within_buffer(dl$dl_grn, dl$dl_pop, dl$dl_bnd, 500, CRS_DELFT)

message(sprintf("Yuexiu: %.1f%% within 300m, %.1f%% within 500m", yx_300, yx_500))
message(sprintf("Delft:  %.1f%% within 300m, %.1f%% within 500m", dl_300, dl_500))


# ── 1C. Mean nearest green space distance per subdistrict ────────────────────
nearest_park_dist <- function(admin_sf, green_sf, local_crs) {
  admin_m  <- st_transform(admin_sf, local_crs)
  green_m  <- st_transform(green_sf, local_crs)
  cents    <- st_centroid(admin_m)
  idx      <- st_nearest_feature(cents, green_m)
  dists    <- st_distance(cents, green_m[idx, ], by_element = TRUE)
  admin_sf |> mutate(nearest_green_m = as.numeric(dists))
}

yx_sub_access  <- nearest_park_dist(yx_sub_access,  d$yx_grn,  CRS_YX)
dl_wijk_access <- nearest_park_dist(dl_wijk_access, dl$dl_grn, CRS_DELFT)


# ── Save Data Matrices ────────────────────────────────────────────────────────
buffer_summary <- data.frame(
  city     = c("Yuexiu", "Yuexiu", "Delft", "Delft"),
  buffer_m = c(300, 500, 300, 500),
  pct_pop  = c(yx_300, yx_500, dl_300, dl_500)
)

saveRDS(yx_sub_access,   file.path(OUT_ROOT, "yx_sub_access.rds"))
saveRDS(dl_wijk_access,  file.path(OUT_ROOT, "dl_wijk_access.rds"))
saveRDS(buffer_summary,  file.path(OUT_ROOT, "buffer_summary.rds"))


# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 1A: Population density maps
density_breaks <- c(0, 1000, 2500, 5000, 10000, 20000, 40000, 80000, Inf)

density_labels <- c("0–1k", "1k–2.5k", "2.5k–5k", "5k–10k","10k–20k", "20k–40k", "40k–80k", ">80k")

density_cols <- c(
  COLORS$beige,
  "#e7d6bc",
  "#d8b48a",
  "#c98f5f",
  COLORS$red_light,
  "#7f3f3e",
  COLORS$red,
  "#4a1f1f")

yx_sub_access <- yx_sub_access %>%
  mutate(area_km2 = as.numeric(st_area(st_transform(., CRS_YX))) / 1e6, pop_density_km2 = pop_count / area_km2, pop_density_class = cut( pop_density_km2, breaks = density_breaks, labels = density_labels, include.lowest = TRUE))

dl_wijk_access <- dl_wijk_access %>%
  mutate(pop_density_km2 = bevolkingsdichtheidInwonersPerKm2, pop_density_class = cut( pop_density_km2, breaks = density_breaks, labels = density_labels,include.lowest = TRUE))

# transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_access, CRS_YX)
dl_plot <- st_transform(dl_wijk_access, 28992)  # RD New

# calculate real maps widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin


p_den_yx <- ggplot(yx_plot) +
  geom_sf(aes(fill = pop_density_class), linewidth = 0.25, color = "white")  +
  geom_sf(aes(fill = pop_density_class), color = NA, alpha = 0.001, show.legend = TRUE) +
  annotation_scale( location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4 ) +
  scale_fill_manual( name = "residents/km²", values = setNames(density_cols, density_labels), limits = density_labels, breaks = density_labels, drop = FALSE, na.value = COLORS$grey85, guide = guide_legend( reverse = TRUE, override.aes = list( fill = rev(density_cols), alpha = 1,color = NA))) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs( title = "Yuexiu", subtitle = "Source: WorldPop + administrative polygons")

p_den_dl <- ggplot(dl_plot) +
  geom_sf(aes(fill = pop_density_class), linewidth = 0.25, color = "white") +
  scale_fill_manual( name = "residents/km²", values = setNames(density_cols, density_labels), breaks = density_labels, drop = FALSE, na.value = COLORS$grey85, guide = "none") +
  coord_sf(expand = FALSE, datum = NA) +
theme_map_clean() +
  labs( title = "Delft", subtitle = "Source: CBS wijk data")

fig_population_density <- p_den_yx + p_den_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect")+
  plot_annotation( title = "Population Density", theme = theme( plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_population_density

ggsave(file.path(OUT_ROOT, "fig_population_density.png"), fig_population_density, width = 14, height = 6, dpi = 300)

# Figure 1B: Green space per capita maps
green_colours <- c("#f4ecd0", "#e6dfbb", "#d5d0a5", "#c3c193", "#b1b181", "#92975f", "#5c612f", "#252f18", "#111809")
green_breaks  <- c(0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 100, 500)

# transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_access, CRS_YX)
dl_plot <- st_transform(dl_wijk_access, 28992)  # RD New

# calculate real maps widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin


p1 <- ggplot(yx_sub_access) + geom_sf(aes(fill = green_pc_m2), linewidth = 0.25, color = "white") +
  annotation_scale( location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  scale_fill_stepsn(name = "m² per person", colours = green_colours, values = scales::rescale(green_breaks), breaks = green_breaks, limits = c(0, 500), na.value = "grey85") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs( title = "Yuexiu", subtitle = "source")

p2 <- ggplot(dl_wijk_access) + geom_sf(aes(fill = green_pc_m2)) +
  scale_fill_stepsn(name = "m² per person", colours = green_colours, values = scales::rescale(green_breaks), breaks = green_breaks, limits = c(0, 500), na.value = "grey85") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs( title = "Delft", subtitle = "source" )

fig_green_per_capita <- p1 + p2 +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation (title = "Green space per capita", theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_green_per_capita

ggsave(file.path(OUT_ROOT, "fig_access_per_capita.png"), fig_green_per_capita, width = 14, height = 6, dpi = 300)


# Figure 1C: Mean nearest green distance maps
# transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_access, CRS_YX)
dl_plot <- st_transform(dl_wijk_access, 28992)  # RD New

# calculate real maps widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

p3 <- ggplot(yx_plot) +
  geom_sf(aes(fill = nearest_green_m), linewidth = 0.25, color = "white") +
  annotation_scale( location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4 ) +
  scale_fill_gradient(name = "Distance (m)", low = COLORS$pink_light, high = COLORS$red, na.value = COLORS$grey85) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Yuexiu", subtitle = "Source")

p4 <- ggplot(dl_plot) +
  geom_sf(aes(fill = nearest_green_m), linewidth = 0.25, color = "white") +
  scale_fill_gradient(name = "Distance (m)", low = COLORS$pink_light, high = COLORS$red, na.value = COLORS$grey85, guide = "none") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "Source")
fig_nearest_green_distance <- p3 + p4 +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") + plot_annotation(title = "Mean Nearest Distance to Urban Green Space", theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_nearest_green_distance

ggsave(file.path(OUT_ROOT, "fig_nearest_green_distance.png"), fig_nearest_green_distance, width = 14, height = 6, dpi = 300)

# Figure 1D: Buffer coverage line chart
#changed to linechart to show the differences in a more readable way
p5 <- ggplot(buffer_summary, aes( x = buffer_m, y = pct_pop, colour = city, group = city)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("Yuexiu" = COLORS$orange, "Delft" = COLORS$blue), name = NULL ) +
  scale_x_continuous(breaks = sort(unique(buffer_summary$buffer_m))) +
  scale_y_continuous( limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs( title = "Population within Walking Distance of Green Space", x = "Buffer distance (m)", y = "Population covered (%)") +
  theme_map_clean() +
  theme( plot.title = element_text( face = "bold", size = 14), legend.position = "right", axis.text = element_text(size = 10), axis.title = element_text(size = 11))

p5

ggsave(file.path(OUT_ROOT, "fig_buffer_coverage.png"), p5, width = 8, height = 5, dpi = 300)

# ── Figure 1E: Network Entry Gates & Multi-Ring Walking Buffers Panel ───────────
generate_network_access <- function(green_sf, road_sf, local_crs) {
  green_m <- st_transform(green_sf, local_crs)
  road_m  <- st_transform(road_sf, local_crs)

  buf_500  <- st_union(st_buffer(green_m, 500))
  buf_300  <- st_union(st_buffer(green_m, 300))
  buf_100  <- st_union(st_buffer(green_m, 100))

  entry_buffer <- st_buffer(st_boundary(green_m), 3)
  road_m_filtered <- road_m[road_m$highway %in% c("primary", "secondary", "tertiary", "residential", "living_street", "unclassified", "service", "pedestrian", "cycleway", "footway"), ]

  green_core <- st_buffer(st_union(green_m), -10)
  roads_outside_green <- st_difference(road_m_filtered, green_core)
  entry_segments <- st_intersection(roads_outside_green, entry_buffer)

  access_pts_raw <- st_centroid(entry_segments)
  access_pts <- access_pts_raw |> st_buffer(35) |> st_union() |> st_cast("POLYGON") |> st_centroid() |> st_as_sf()

  list(b100 = buf_100, b300 = buf_300, b500 = buf_500, entry_buffer = entry_buffer, entry_segments = entry_segments, roads = road_m, green = green_m, pts = access_pts)
}

# transform both maps to metre-based CRS
yx_plot <- st_transform(yx_sub_access, CRS_YX)
dl_plot <- st_transform(dl_wijk_access, 28992)  # RD New


yx_net <- generate_network_access(d$yx_grn,  d$yx_rds, CRS_YX)
dl_net <- generate_network_access(dl$dl_grn, dl$dl_rds, CRS_DELFT)


# Calculate real map widths from network layers
yx_bb <- st_bbox(yx_net$green)
dl_bb <- st_bbox(dl_net$green)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin


p_yx_network <- ggplot() +
  geom_sf(data = yx_net$roads, aes(color = "Roads"), linewidth = 0.2) +
  geom_sf(data = yx_net$b300, aes(fill = "300 m buffer"), color = NA, alpha = 0.18) +
  geom_sf(data = yx_net$green, aes(fill = "Green patches"), color = NA) +
  geom_sf(data = yx_net$pts, aes(color = "Access points"), size = 0.5, alpha = 0.8) +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  scale_color_manual( name = NULL, values = c( "Roads" = "grey75", "Access points" = COLORS$red)) +
  scale_fill_manual( name = NULL, values = c( "300 m buffer" = COLORS$pink_light, "Green patches" = COLORS$green_dark)) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Yuexiu", subtitle = "Network entry thresholds")

p_dl_network <- ggplot() +
  geom_sf(data = dl_net$roads, aes(color = "Roads"), linewidth = 0.2) +
  geom_sf(data = dl_net$b300, aes(fill = "300 m buffer"), color = NA, alpha = 0.18) +
  geom_sf(data = dl_net$green, aes(fill = "Green patches"), color = NA) +
  geom_sf(data = dl_net$pts, aes(color = "Access points"), size = 0.5, alpha = 0.8) +
  scale_color_manual(name = NULL, values = c("Roads" = "grey75", "Access points" = COLORS$red), guide = "none") +
  scale_fill_manual(name = NULL, values = c("300 m buffer" = COLORS$pink_light, "Green patches" = COLORS$green_dark), guide = "none") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "Network entry thresholds")

fig_network_walk_access <- p_yx_network + p_dl_network +
  plot_layout( widths = c(yx_w, dl_w),
    guides = "collect") +
  plot_annotation(title = "Network Entry Gates & Walking Buffer", theme = theme( plot.title = element_text(face = "bold", size = 14))) &
  theme( legend.position = "right")

ggsave(
  file.path(OUT_ROOT, "fig_network_walk_access.png"), fig_network_walk_access, width = 14, height = 6, dpi = 300)

message("All figures corrected and saved successfully!")
