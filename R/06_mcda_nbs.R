source("R/00_config.R")
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

sf_use_s2(FALSE)

# ── Load all prior outputs ────────────────────────────────────────────────────
d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

yx_sub_access  <- readRDS(file.path(OUT_ROOT, "yx_sub_access.rds"))
dl_wijk_access <- readRDS(file.path(OUT_ROOT, "dl_wijk_access.rds"))

yx_sub_ndvi    <- readRDS(file.path(OUT_ROOT, "yx_sub_ndvi.rds"))
dl_wijk_ndvi   <- readRDS(file.path(OUT_ROOT, "dl_wijk_ndvi.rds"))

yx_bivar       <- readRDS(file.path(OUT_ROOT, "yx_bivar.rds"))
dl_bivar       <- readRDS(file.path(OUT_ROOT, "dl_bivar.rds"))

yx_graph       <- readRDS(file.path(OUT_ROOT, "yx_graph.rds"))
dl_graph       <- readRDS(file.path(OUT_ROOT, "dl_graph.rds"))

# ── FIX: pop-floor-artifact filter (see comment from earlier) ──
POP_MIN <- 100

yx_sub_access  <- yx_sub_access  |> filter(!is.na(pop_count) & pop_count >= POP_MIN)
dl_wijk_access <- dl_wijk_access |> filter(!is.na(pop_count) & pop_count >= POP_MIN)
message(sprintf("MCDA input after pop filter — Yuexiu: %d | Delft: %d",
                nrow(yx_sub_access), nrow(dl_wijk_access)))

yx_sub_ndvi <- yx_sub_ndvi |> filter(full_id %in% yx_sub_access$full_id)

# MCDA weights from config
w <- MCDA_WEIGHTS


# ── Helper: min-max normalise to [0, 1] ───────────────────────────────────────
norm_minmax <- function(x, invert = FALSE) {
  mn <- min(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  if (mx == mn) return(rep(0.5, length(x)))
  s <- (x - mn) / (mx - mn)
  if (invert) 1 - s else s
}

# ── 5A. Build MCDA score table per admin unit ─────────────────────────────────
# For each city separately (different admin units), assemble four sub-scores
# then compute the weighted MCDA composite.

build_mcda <- function(access_sf, ndvi_sf, bivar_sf, graph_obj,
                        local_crs, city_label) {

  # --- Accessibility sub-score ---
  # Higher green_pc_m2 → better; lower nearest_green_m → better
  acc <- access_sf |>
    st_drop_geometry() |>
    mutate(
      s_green_pc   = norm_minmax(green_pc_m2, invert = TRUE),
      s_nearest    = norm_minmax(nearest_green_m, invert = FALSE),
      score_access = (s_green_pc + s_nearest) / 2
    ) |>
    select(matches("^(name|naam|BU_|WK_|GEO_|jiedao|subdistrict)", ignore.case = TRUE),
           score_access)

  # Use geometry from access_sf as the spatial backbone
  backbone <- access_sf |>
    mutate(row_id = row_number())

  # --- Biodiversity sub-score ---
  # Use mean NDVI per admin unit as proxy (already computed in 03)
  bio <- ndvi_sf |>
    st_drop_geometry() |>
    # invert: high NDVI = already vegetated = low intervention urgency
    mutate(score_bio = norm_minmax(ndvi_mean, invert = TRUE)) |>
    select(score_bio)

  # --- Equity sub-score ---
  # Pop density in a LOW-green tertile → HIGH equity need → invert so high score = needs intervention
  # For MCDA we treat equity score as "intervention urgency from equity lens":
  #   High pop + Low green → high score
  equity <- bivar_sf |>
    st_drop_geometry() |>
    mutate(
      equity_score_raw = case_when(
        bivar_class == "Low-High"  ~ 1.00,
        bivar_class == "Mid-High"  ~ 0.75,
        bivar_class == "Low-Mid"   ~ 0.65,
        bivar_class == "Mid-Mid"   ~ 0.40,
        bivar_class == "Low-Low"   ~ 0.30,
        bivar_class == "High-Mid"  ~ 0.20,
        bivar_class == "High-High" ~ 0.15,
        bivar_class == "Mid-Low"   ~ 0.10,
        bivar_class == "High-Low"  ~ 0.10,
        TRUE ~ 0.10
      ),
      score_equity = equity_score_raw
    ) |>
    select(score_equity)

  # --- Connectivity sub-score ---
  # For each admin unit: mean betweenness of patches whose centroid falls within it
  admin_m  <- st_transform(access_sf, local_crs)
  nodes_m  <- st_transform(graph_obj$nodes, local_crs)
  node_cents <- st_centroid(nodes_m)

  # Betweenness is ~0 for >75% of patches in both cities — too sparse to
  # produce meaningful per-admin-unit variation after averaging. Degree
  # (number of patch connections within dispersal threshold) varies more
  # continuously and is a better proxy for local connectivity density.
  joined_conn <- st_join(admin_m, node_cents["degree"],
                        join = st_contains) |>
    st_drop_geometry() |>
    group_by(row_number()) |>
    summarise(mean_degree = mean(degree, na.rm = TRUE), .groups = "drop")

  # Invert: low mean degree = poorly connected patches = high intervention need.
  conn_scores <- joined_conn |>
    mutate(
      mean_degree = replace_na(mean_degree, 0),
      score_conn = 1 - norm_minmax(mean_degree)
    ) |>
    select(score_conn)
  message(city_label, " connectivity mean degree:")
  print(summary(joined_conn$mean_degree))

  # --- Assemble & compute weighted composite ---
  n <- nrow(access_sf)
  # Ensure all sub-score vectors are length n
  pad <- function(x, len) { if (length(x) < len) c(x, rep(NA_real_, len - length(x))) else x[seq_len(len)] }

  score_access  <- pad(acc$score_access,    n)
  score_bio     <- pad(bio$score_bio,       n)
  score_equity  <- pad(replace_na(equity$score_equity, 0.10), n)
  score_equity  <- replace_na(score_equity, 0.10)
  score_conn    <- pad(conn_scores$score_conn, n)

  access_sf |>
    mutate(
      score_access   = score_access,
      score_bio      = score_bio,
      score_equity   = score_equity,
      score_conn     = score_conn,
      mcda_composite = w["accessibility"] * score_access  +
                       w["biodiversity"]  * score_bio     +
                       w["connectivity"]  * score_conn    +
                       w["equity"]        * score_equity,
      city           = city_label,
      # Priority tier
      priority_tier  = case_when(
        mcda_composite >= quantile(mcda_composite, 0.67, na.rm = TRUE) ~ "High priority",
        mcda_composite >= quantile(mcda_composite, 0.33, na.rm = TRUE) ~ "Medium priority",
        TRUE                                                            ~ "Low priority"
      )
    )
}

message("Building MCDA scores — Yuexiu...")
yx_mcda <- build_mcda(yx_sub_access, yx_sub_ndvi, yx_bivar, yx_graph,
                       CRS_YX, "Yuexiu")

message("Building MCDA scores — Delft...")
dl_mcda <- build_mcda(dl_wijk_access, dl_wijk_ndvi, dl_bivar, dl_graph,
                       CRS_DELFT, "Delft")

# ── 5B. NbS corridor prioritisation ───────────────────────────────────────────
# Identify potential NbS intervention zones:
#   Criteria: HIGH priority MCDA tier + admin units with low connectivity patches
#             → candidate corridors between isolated patches

identify_nbs_corridors <- function(mcda_sf, graph_obj, local_crs, thresh_m = DISPERSAL_THRESH_M) {
  high_priority <- mcda_sf |>
    filter(priority_tier == "High priority") |>
    st_transform(local_crs)

  # Isolated patches: enn_m > dispersal threshold (not currently connected)
  isolated_patches <- graph_obj$nodes |>
    st_transform(local_crs) |>
    filter(enn_m > thresh_m)

  if (nrow(isolated_patches) == 0) {
    message("  No isolated patches at threshold ", thresh_m, " m")
    return(list(high_priority = high_priority, isolated_patches = isolated_patches,
                corridor_lines = NULL))
  }

  # Potential corridors: lines from each isolated patch to its nearest neighbour
  # regardless of threshold — showing the gap that needs bridging
  cents <- st_centroid(isolated_patches)

  # Exclude isolated patches from targets so each patch finds a *different* neighbour
  connected_patches <- graph_obj$nodes |>
    st_transform(local_crs) |>
    filter(enn_m <= thresh_m)

  all_cents <- st_centroid(connected_patches)
  nearest_idx <- st_nearest_feature(cents, all_cents)
  corridor_lines <- lapply(seq_len(nrow(cents)), function(i) {
    a <- st_coordinates(cents[i, ])
    b <- st_coordinates(all_cents[nearest_idx[i], ])
    st_linestring(rbind(a, b))
  })
  corridors_sf <- st_sf(
    from_patch   = isolated_patches$area_ha,
    gap_m        = isolated_patches$enn_m,
    geometry     = st_sfc(corridor_lines, crs = local_crs)
  )

  list(high_priority   = high_priority,
       isolated_patches = isolated_patches,
       corridor_lines  = corridors_sf)
}

message("Identifying NbS corridors — Yuexiu...")
yx_nbs <- identify_nbs_corridors(yx_mcda, yx_graph, CRS_YX)

message("Identifying NbS corridors — Delft...")
dl_nbs <- identify_nbs_corridors(dl_mcda, dl_graph, CRS_DELFT)

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(yx_mcda, file.path(OUT_ROOT, "yx_mcda.rds"))
saveRDS(dl_mcda, file.path(OUT_ROOT, "dl_mcda.rds"))
saveRDS(yx_nbs,  file.path(OUT_ROOT, "yx_nbs.rds"))
saveRDS(dl_nbs,  file.path(OUT_ROOT, "dl_nbs.rds"))

# MCDA weight table for reporting
weight_df <- data.frame(
  Criterion    = c("Accessibility", "Biodiversity / NDVI", "Connectivity", "Equity"),
  Weight       = unname(w),
  Sub_scores   = c("Green pc + nearest distance",
                   "Mean NDVI per unit",
                   "Inverse mean betweenness",
                   "Pop × green density matrix")
)
write.csv(weight_df, file.path(OUT_ROOT, "mcda_weights.csv"), row.names = FALSE)
message("MCDA weights saved to mcda_weights.csv")

# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 5A: MCDA composite score maps

# Transform both maps to metre-based CRS
yx_plot <- st_transform(yx_mcda, CRS_YX)
dl_plot <- st_transform(dl_mcda, CRS_DELFT)

# Calculate real map widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

# Shared limits so colours are comparable
mcda_limits <- c(0, 1)

p_mcda_yx <- ggplot(yx_plot) +
  geom_sf(aes(fill = mcda_composite), colour = "white", linewidth = 0.3) +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  scale_fill_gradientn(colours = pal_full,limits = mcda_limits, name = "MCDA score\n(0–1)", na.value = COLORS$grey85) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Yuexiu")

p_mcda_dl <- ggplot(dl_plot) +
  geom_sf(aes(fill = mcda_composite), colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = pal_full, limits = mcda_limits, name = "MCDA score\n(0–1)", na.value = COLORS$grey85, guide = "none") +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft")

fig_mcda_maps <- p_mcda_yx + p_mcda_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation(
    title = "MCDA Composite Score",
    subtitle = sprintf("Weights: Acc %.0f%% | Bio %.0f%% | Conn %.0f%% | Eq %.0f%%",
      w["accessibility"] * 100,
      w["biodiversity"] * 100,
      w["connectivity"] * 100,
      w["equity"] * 100),
    theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_mcda_maps

ggsave(file.path(OUT_ROOT, "fig_mcda_maps.png"), fig_mcda_maps, width = 14, height = 7, dpi = 300)

# Figure 5B: Priority tier maps
tier_colours <- c("High priority"   = COLORS$red_light, "Medium priority" = COLORS$beige, "Low priority"    = COLORS$green_mid)

# Transform both maps to metre-based CRS
yx_plot <- st_transform(yx_mcda, CRS_YX)
dl_plot <- st_transform(dl_mcda, CRS_DELFT)

# Calculate real map widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

p_tier_yx <- ggplot(yx_plot) +
  geom_sf(aes(fill = priority_tier), colour = "white", linewidth = 0.3) +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  scale_fill_manual(values = tier_colours, name = "NbS Priority",  breaks = c("High priority", "Medium priority", "Low priority"), drop = FALSE) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Yuexiu")

p_tier_dl <- ggplot(dl_plot) +
  geom_sf(aes(fill = priority_tier), colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = tier_colours,
    name = "NBS Priority",
    breaks = c("High priority", "Medium priority", "Low priority"),
    drop = FALSE,
    guide = "none"
  ) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(
    title = "Delft"
  )

fig_priority_tiers <- p_tier_yx + p_tier_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation( title = "NBS Intervention Priority", subtitle = "Tertile classification of MCDA composite score", theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_priority_tiers

ggsave(file.path(OUT_ROOT, "fig_priority_tiers.png"), fig_priority_tiers, width = 14, height = 7, dpi = 300)

# Figure 5C: NbS corridor maps (isolated patches + proposed links)
# Transform both maps to metre-based CRS
yx_plot <- st_transform(d$yx_bnd, CRS_YX)
dl_plot <- st_transform(dl$dl_bnd, CRS_DELFT)

# Calculate real map widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

plot_nbs_map <- function(nbs_obj, bnd_sf, local_crs, city_label, show_scale = FALSE) {
  # Pre-transform all layers to local_crs before passing to ggplot.
  # Without this, ggplot inherits the CRS of the first geom_sf layer (bnd in
  # WGS84) and renders subsequent layers — especially corridor lines already
  # in UTM/RD New — outside the map extent, making them invisible.
    bnd           <- st_transform(bnd_sf, local_crs)
    high_priority <- st_transform(nbs_obj$high_priority, local_crs)
    isolated      <- st_transform(nbs_obj$isolated_patches, local_crs)
    corridors     <- if (!is.null(nbs_obj$corridor_lines))
      st_transform(nbs_obj$corridor_lines, local_crs)
    else NULL

    p <- ggplot() +
      geom_sf(data = bnd, aes(fill = "Urban fabric"), colour = COLORS$grey50, linewidth = 0.5) +
      geom_sf(data = high_priority, aes(fill = "High-priority areas"), colour = NA, alpha = 0.55)

    if (!is.null(corridors) && nrow(corridors) > 0) {
      p <- p +
        geom_sf(data = corridors, aes(colour = "Proposed corridors"), linewidth = 1.4, linetype = "11")
    }

    p <- p +
      geom_sf(data = isolated, aes(fill = "Isolated patches"), colour = COLORS$green_mid, linewidth = 0.1, alpha = 1) +
      scale_fill_manual(name = NULL, values = c( "Urban fabric" = COLORS$beige, "High-priority areas" = COLORS$red_light, "Isolated patches" = COLORS$green_light)) +
      scale_colour_manual(name = NULL, values = c("Proposed corridors" = COLORS$green_dark)) +
      coord_sf(expand = FALSE, datum = NA) +
      theme_map_clean() +
      labs( title = city_label)

    if (show_scale) {
      p <- p +
        annotation_scale( location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4)
    }
    p}


p_nbs_yx <- plot_nbs_map( yx_nbs, d$yx_bnd, CRS_YX, "Yuexiu", show_scale = TRUE)

p_nbs_dl <- plot_nbs_map(dl_nbs, dl$dl_bnd, CRS_DELFT, "Delft", show_scale = FALSE) +
  guides(fill = "none", colour = "none")

fig_nbs_corridors <- p_nbs_yx + p_nbs_dl +
  plot_layout(widths = c(yx_w, dl_w), guides = "collect") +
  plot_annotation(title = "Nature-based Solution Corridor Prioritisation", subtitle = "High-priority areas, isolated patches and proposed corridors", theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_nbs_corridors

ggsave(file.path(OUT_ROOT, "fig_nbs_corridors.png"), fig_nbs_corridors, width = 14, height = 7, dpi = 300)

# Figure 5D: Sub-score radar / spider chart (city-level means)
score_means <- bind_rows(
  yx_mcda |> st_drop_geometry() |>
    summarise(across(c(score_access, score_bio, score_equity, score_conn), mean, na.rm = TRUE)) |>
    mutate(city = "Yuexiu"),
  dl_mcda |> st_drop_geometry() |>
    summarise(across(c(score_access, score_bio, score_equity, score_conn), mean, na.rm = TRUE)) |>
    mutate(city = "Delft")
) |>
  pivot_longer(-city, names_to = "criterion", values_to = "score") |>
  mutate(criterion = recode(criterion,
    "score_access"  = "Accessibility",
    "score_bio"     = "Biodiversity",
    "score_equity"  = "Equity",
    "score_conn"    = "Connectivity"
  ))

# coord_polar() removed: polar layout compresses small differences between
# cities and makes directional comparison (which city scores higher on each
# criterion) very difficult to read. Standard dodged bar chart shows the
# same information clearly with readable absolute values on the y-axis.
p_radar_bar <- ggplot(score_means, aes(x = criterion, y = score, fill = city)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dotted", colour = "grey70") +
  geom_text(aes(label = sprintf("%.2f", score)),
            position = position_dodge(0.7), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Yuexiu" = COLORS$orange_light, "Delft" = COLORS$blue_light)) +
  scale_y_continuous(limits = c(0, 1.05), labels = label_percent()) +
  theme_minimal(base_size = 12) +
  labs(title = "MCDA Sub-score Profile — City Comparison",
       subtitle = "Mean normalised intervention urgency score (higher = greater need)",
       x = NULL, y = "Mean normalised score", fill = "City") +
  theme(plot.title = element_text(face = "bold"))

# Wider aspect ratio suits a horizontal bar comparison better than the
# square format used for the polar chart.
ggsave(file.path(OUT_ROOT, "fig_mcda_radar.png"),
       p_radar_bar, width = 10, height = 6, dpi = 300)

# Figure 5E: MCDA score distribution comparison
score_dist <- bind_rows(
  yx_mcda |> st_drop_geometry() |> transmute(city = "Yuexiu", mcda_composite),
  dl_mcda |> st_drop_geometry() |> transmute(city = "Delft",  mcda_composite)
)

p_dist <- ggplot(score_dist, aes(x = mcda_composite, fill = city)) +
  geom_density(alpha = 0.6) +
  geom_rug(aes(colour = city), alpha = 0.6) +
  scale_fill_manual(values  = c("Yuexiu" = COLORS$orange_light, "Delft" = COLORS$blue_light)) +
  scale_colour_manual(values = c("Yuexiu" = COLORS$orange_light, "Delft" = COLORS$blue_light),
                      guide = "none") +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of MCDA Composite Scores",
       subtitle = "Density plot — all admin units",
       x = "MCDA composite score (0–1)", y = "Density", fill = "City") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_ROOT, "fig_mcda_distribution.png"),
       p_dist, width = 8, height = 5, dpi = 300)

message("Script 06 complete — MCDA & NbS figures saved.")
message("")
message("=== All scripts complete. Output files in: ", OUT_ROOT, " ===")
