# Context maps for presentation and report introduction
# Setup ------------------------------------------------------
source("R/00_config.R")

library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggspatial)

# ── Load boundaries ───────────────────────────────────────────────────────────
gz_bnd  <- st_read("data/Yuexiu/vector/guangzhou_boundary_proj.gpkg")
yx_bnd  <- st_read(YX$boundary)
dl_bnd  <- st_read(DL$boundary)

# ── Study area units (from config globals) ────────────────────────────────────
dl_wijk_context <- dl_wijk
yx_sub_context  <- yx_sub

# ── 0A: Guangzhou + Delft context ────────────────────────────────────────────
p_context_gz <- ggplot() +
  geom_sf(data = gz_bnd,  fill = COLORS$orange_light, color = NA) +
  geom_sf(data = yx_bnd,  fill = COLORS$orange, color = "white", linewidth = 0.12) +
  geom_sf(data = gz_bnd,  fill = NA, color = COLORS$red_light, linewidth = 0.55) +
  annotation_scale(location = "bl", style = "ticks",
                   line_width = 0.6, text_cex = 0.7, width_hint = 0.2) +
  theme_map_clean() +
  labs(title = "Guangzhou", subtitle = "Source: GADM / OSM")

p_context_dl <- ggplot() +
  geom_sf(data = dl_bnd,          fill = COLORS$blue_light, color = NA) +
  geom_sf(data = dl_wijk_context, fill = NA, color = "white", linewidth = 0.12) +
  geom_sf(data = dl_bnd,          fill = NA, color = COLORS$red_light, linewidth = 0.55) +
  annotation_scale(location = "bl", style = "ticks",
                   unit_category = "metric", width_hint = 0.3,
                   text_cex = 0.7, line_width = 0.6) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "Source: CBS wijk data")

fig_context_gz_dl <- p_context_gz + p_context_dl +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Context",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

fig_context_gz_dl
ggsave(file.path(OUT_ROOT, "fig_context_gz_dl.png"),
       fig_context_gz_dl, width = 14, height = 6, dpi = 300)

# ── 0B: Spatial units — Yuexiu jiēdào + Delft wijken ─────────────────────────
yx_plot <- st_transform(yx_sub_context, CRS_YX)
dl_plot <- st_transform(dl_wijk_context, CRS_DELFT)

yx_w <- diff(st_bbox(yx_plot)[c("xmin","xmax")])
dl_w <- diff(st_bbox(dl_plot)[c("xmin","xmax")])

p_su_yx <- ggplot() +
  geom_sf(data = yx_bnd,          fill = COLORS$orange_light, color = NA) +
  geom_sf(data = yx_sub_context,  fill = NA, color = "white", linewidth = 0.12) +
  geom_sf(data = yx_bnd,          fill = NA, color = COLORS$red_light, linewidth = 0.55) +
  annotation_scale(location = "bl", style = "ticks",
                   width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Yuexiu", subtitle = "Source: administrative polygons")

p_su_dl <- ggplot() +
  geom_sf(data = dl_bnd,          fill = COLORS$blue_light, color = NA) +
  geom_sf(data = dl_wijk_context, fill = NA, color = "white", linewidth = 0.12) +
  geom_sf(data = dl_bnd,          fill = NA, color = COLORS$red_light, linewidth = 0.55) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "Source: CBS wijk data")

fig_context_yx_dl <- p_su_yx + p_su_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation(title = "Spatial units",
                  theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_context_yx_dl
ggsave(file.path(OUT_ROOT, "fig_context_yx_dl.png"),
       fig_context_yx_dl, width = 14, height = 6, dpi = 300)

# ── 0C: Urban green space maps ────────────────────────────────────────────────
dl_roads <- st_intersection(dl_net$roads, dl_bnd)
dl_green <- st_intersection(dl_net$green, dl_bnd)

yx_roads <- st_transform(yx_net$roads, st_crs(yx_bnd)) |> st_intersection(yx_bnd)
yx_green <- st_transform(yx_net$green, st_crs(yx_bnd)) |> st_intersection(yx_bnd)

dl_water_context <- st_intersection(dl_water, dl_bnd)
yx_water_context <- st_transform(yx_water_poly, st_crs(yx_bnd)) |> st_intersection(yx_bnd)

p_ugs_yx <- ggplot() +
  geom_sf(data = yx_bnd,           aes(fill = "Urban fabric"), color = NA) +
  geom_sf(data = yx_water_context, aes(fill = "Water"), color = NA, alpha = 0.6) +
  geom_sf(data = yx_roads,         aes(color = "Roads"), linewidth = 0.15) +
  geom_sf(data = yx_green,         aes(fill = "Urban green space"), color = NA, alpha = 0.9) +
  geom_sf(data = yx_sub_context,   aes(color = "Administrative boundaries"), fill = NA, linewidth = 0.16) +
  geom_sf(data = yx_bnd,           aes(color = "City boundary"), fill = NA, linewidth = 0.55) +
  annotation_scale(location = "bl", style = "ticks",
                   width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_fill_manual(name = NULL,
                    values = c("Urban fabric"      = COLORS$beige,
                               "Water"             = COLORS$blue_light,
                               "Urban green space" = COLORS$green_dark)) +
  scale_color_manual(name = NULL,
                     values = c("Roads"                    = COLORS$grey85,
                                "Administrative boundaries" = COLORS$red,
                                "City boundary"            = COLORS$red_light)) +
  theme_map_clean() +
  labs(title = "Yuexiu — Jiēdào", subtitle = "Source: OSM + administrative polygons")

p_ugs_dl <- ggplot() +
  geom_sf(data = dl_bnd,           aes(fill = "Urban fabric"), color = NA) +
  geom_sf(data = dl_water_context, aes(fill = "Water"), color = NA, alpha = 0.6) +
  geom_sf(data = dl_roads,         aes(color = "Roads"), linewidth = 0.10) +
  geom_sf(data = dl_green,         aes(fill = "Urban green space"), color = NA, alpha = 0.95) +
  geom_sf(data = dl_wijk_context,  aes(color = "Administrative boundaries"), fill = NA, linewidth = 0.16) +
  geom_sf(data = dl_bnd,           aes(color = "City boundary"), fill = NA, linewidth = 0.55) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_fill_manual(name = NULL,
                    values = c("Urban fabric"      = COLORS$beige,
                               "Water"             = COLORS$blue_light,
                               "Urban green space" = COLORS$green_dark),
                    guide = "none") +
  scale_color_manual(name = NULL,
                     values = c("Roads"                    = COLORS$grey85,
                                "Administrative boundaries" = COLORS$red,
                                "City boundary"            = COLORS$red_light),
                     guide = "none") +
  theme_map_clean() +
  labs(title = "Delft — Wijken", subtitle = "Source: CBS wijk data")

fig_ugs_yx_dl <- p_ugs_yx + p_ugs_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation(title = "Urban green space",
                  theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_ugs_yx_dl
ggsave(file.path(OUT_ROOT, "fig_ugs_yx_dl.png"),
       fig_ugs_yx_dl, width = 14, height = 6, dpi = 300)

message("Script 08 complete — context figures saved.")