source("R/00_config.R")
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

sf_use_s2(FALSE)

d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))

yx_grn_bio <- readRDS(file.path(OUT_ROOT, "yx_grn_bio.rds"))
dl_grn_bio <- readRDS(file.path(OUT_ROOT, "dl_grn_bio.rds"))

# ── 4A. Fragmentation indices (landscape metrics) ─────────────────────────────
# Computed directly from the green patch vectors; mirrors FRAGSTATS metrics:
#   NP  = Number of Patches
#   MPS = Mean Patch Size (ha)
#   ENN = Euclidean Nearest-Neighbour distance (mean inter-patch gap, m)

calc_fragmentation <- function(green_sf, local_crs, city_label) {
  grn_m <- st_transform(green_sf, local_crs)

  # NP and MPS
  np  <- nrow(grn_m)
  mps <- mean(grn_m$area_ha, na.rm = TRUE)

  # ENN: centroid-to-centroid distance to nearest other patch
  cents <- st_centroid(grn_m)
  dist_mat <- st_distance(cents)
  diag(dist_mat) <- NA
  enn_per_patch <- apply(dist_mat, 1, function(r) min(r, na.rm = TRUE))
  mean_enn <- mean(as.numeric(enn_per_patch), na.rm = TRUE)
  sd_enn   <- sd(as.numeric(enn_per_patch), na.rm = TRUE)

  # Patch size distribution stats
  total_area_ha <- sum(grn_m$area_ha, na.rm = TRUE)
  pct_large     <- sum(grn_m$area_ha >= 1, na.rm = TRUE) / np * 100  # patches >= 1 ha

  data.frame(
    city          = city_label,
    n_patches     = np,
    mean_patch_ha = mps,
    total_area_ha = total_area_ha,
    mean_enn_m    = mean_enn,
    sd_enn_m      = sd_enn,
    pct_large_patches = pct_large
  )
}

message("Fragmentation indices — Yuexiu...")
yx_frag <- calc_fragmentation(yx_grn_bio, CRS_YX, "Yuexiu")

message("Fragmentation indices — Delft...")
dl_frag <- calc_fragmentation(dl_grn_bio, CRS_DELFT, "Delft")

frag_summary <- bind_rows(yx_frag, dl_frag)
message("Fragmentation summary:")
print(frag_summary)

# ── 4B. Per-patch ENN for visualisation ───────────────────────────────────────
add_enn <- function(green_sf, local_crs) {
  grn_m  <- st_transform(green_sf, local_crs)
  cents  <- st_centroid(grn_m)
  dist_mat <- st_distance(cents)
  diag(dist_mat) <- NA
  enn <- apply(dist_mat, 1, function(r) as.numeric(min(r, na.rm = TRUE)))
  grn_m |> mutate(enn_m = enn)
}

message("Per-patch ENN — Yuexiu...")
yx_grn_conn <- add_enn(yx_grn_bio, CRS_YX)

message("Per-patch ENN — Delft...")
dl_grn_conn <- add_enn(dl_grn_bio, CRS_DELFT)

# ── 4C. Functional connectivity graph ─────────────────────────────────────────
# Nodes  = green patch centroids
# Edges  = pairs of patches whose centroid distance <= DISPERSAL_THRESH_M
# Metric = betweenness centrality (identifies stepping-stone patches)

build_connectivity_graph <- function(green_sf, local_crs,
                                      thresh_m = DISPERSAL_THRESH_M) {
  grn_m  <- st_transform(green_sf, local_crs)
  cents  <- st_centroid(grn_m)
  coords <- st_coordinates(cents)
  n      <- nrow(grn_m)

  # Build edge list: all pairs within dispersal threshold
  dist_mat <- units::drop_units(st_distance(cents))
  edges <- which(dist_mat <= thresh_m & dist_mat > 0, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]  # upper triangle only

  if (nrow(edges) == 0) {
    message("  No edges at threshold ", thresh_m, " m — consider increasing DISPERSAL_THRESH_M")
    return(list(nodes = grn_m |> mutate(betweenness = 0, degree = 0),
                edges_sf = NULL, n_components = n))
  }

  # Degree centrality (simple: number of connections per node)
  deg <- tabulate(c(edges[, 1], edges[, 2]), nbins = n)

  # Betweenness centrality — implemented via BFS shortest paths (pure R, no igraph dep)
  # If igraph is available, use it; else fall back to degree as proxy
  betw <- tryCatch({
    if (!requireNamespace("igraph", quietly = TRUE)) stop("no igraph")
    # Build graph with explicit vertex count so isolated patches are retained.
    # graph_from_edgelist only creates vertices it sees in edges — any patch
    # with no neighbours within DISPERSAL_THRESH_M gets dropped, causing a
    # length mismatch with betweenness output (n_vertices != n_patches).
    g <- igraph::make_empty_graph(n = n, directed = FALSE)
    g <- igraph::add_edges(g, t(edges))
    # graph_from_edgelist creates vertices 1..n in input order — no $name attribute.
    # Assign directly; no index mapping needed.
    igraph::betweenness(g, normalized = TRUE)

  }, error = function(e) {
    message("  igraph not available – using degree as betweenness proxy")
    as.numeric(deg) / max(deg, 1)
  })

  # Number of connected components (union-find / BFS approach)
  adj <- vector("list", n)
  for (k in seq_len(nrow(edges))) {
    i <- edges[k, 1]; j <- edges[k, 2]
    adj[[i]] <- c(adj[[i]], j)
    adj[[j]] <- c(adj[[j]], i)
  }
  visited  <- logical(n)
  n_comp   <- 0
  comp_id  <- integer(n)
  for (start in seq_len(n)) {
    if (!visited[start]) {
      n_comp <- n_comp + 1
      queue  <- start
      while (length(queue) > 0) {
        v <- queue[1]; queue <- queue[-1]
        if (visited[v]) next
        visited[v]  <- TRUE
        comp_id[v]  <- n_comp
        # unique() prevents duplicate queue entries if parallel edges exist
        queue <- c(queue, unique(adj[[v]][!visited[adj[[v]]]]))
      }}}

  # Build edge lines for visualisation
  edge_lines <- lapply(seq_len(nrow(edges)), function(k) {
    i <- edges[k, 1]; j <- edges[k, 2]
    st_linestring(rbind(coords[i, ], coords[j, ]))
  })
  edges_sf <- st_sf( from   = edges[, 1], to     = edges[, 2], length_m = as.numeric(dist_mat[edges]), geometry = st_sfc(edge_lines, crs = local_crs))

  nodes_out <- grn_m |>
    mutate( degree       = deg, betweenness  = betw, component_id = comp_id)

  list(nodes = nodes_out, edges_sf = edges_sf, n_components = n_comp, n_edges = nrow(edges))
}

message(sprintf("Building connectivity graph — Yuexiu (threshold: %d m)...", DISPERSAL_THRESH_M))
yx_graph <- build_connectivity_graph(yx_grn_conn, CRS_YX, DISPERSAL_THRESH_M)
message(sprintf("  Yuexiu: %d edges, %d components", yx_graph$n_edges, yx_graph$n_components))

message(sprintf("Building connectivity graph — Delft (threshold: %d m)...", DISPERSAL_THRESH_M))
dl_graph <- build_connectivity_graph(dl_grn_conn, CRS_DELFT, DISPERSAL_THRESH_M)
message(sprintf("  Delft:  %d edges, %d components", dl_graph$n_edges, dl_graph$n_components))

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(frag_summary,    file.path(OUT_ROOT, "frag_summary.rds"))
saveRDS(yx_grn_conn,     file.path(OUT_ROOT, "yx_grn_conn.rds"))
saveRDS(dl_grn_conn,     file.path(OUT_ROOT, "dl_grn_conn.rds"))
saveRDS(yx_graph,        file.path(OUT_ROOT, "yx_graph.rds"))
saveRDS(dl_graph,        file.path(OUT_ROOT, "dl_graph.rds"))

# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 4A: ENN distribution — how isolated are patches?
enn_df <- bind_rows(
  yx_grn_conn |> st_drop_geometry() |> transmute(city = "Yuexiu", enn_m, area_ha),
  dl_grn_conn |> st_drop_geometry() |> transmute(city = "Delft",  enn_m, area_ha))

p_enn <- ggplot(enn_df, aes(x = enn_m, fill = city)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = DISPERSAL_THRESH_M, linetype = "dashed", colour = COLORS$red, linewidth = 0.8) +
  annotate("text", x = DISPERSAL_THRESH_M + 5, y = Inf,
           label = paste0("Dispersal\nthreshold\n", DISPERSAL_THRESH_M, " m"),
           hjust = 0, vjust = 1.2, size = 3.5, colour = COLORS$red) +
  scale_fill_manual(values = c("Yuexiu" = COLORS$orange, "Delft" = COLORS$blue)) +
  facet_wrap(~city, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "Euclidean Nearest-Neighbour Distance Between Green Patches",
       subtitle = "Dashed line = dispersal threshold for connectivity graph",
       x = "Distance to nearest patch (m)", y = "Count", fill = "City") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_ROOT, "fig_enn_distribution.png"),
       p_enn, width = 12, height = 5, dpi = 300)

# Figure 4B: Fragmentation metrics comparison table-chart
frag_long <- frag_summary |>
  select(city, n_patches, mean_patch_ha, mean_enn_m, pct_large_patches) |>
  pivot_longer(-city, names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
    "n_patches"          = "No. of patches",
    "mean_patch_ha"      = "Mean patch size (ha)",
    "mean_enn_m"         = "Mean ENN (m)",
    "pct_large_patches"  = "% patches ≥ 1 ha"))

p_frag <- ggplot(frag_long, aes(x = city, y = value, fill = city)) +
  geom_col(width = 0.55) +
  scale_fill_manual(values = c("Yuexiu" = COLORS$orange, "Delft" = COLORS$blue),
                    guide = "none") +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 11) +
  labs(title = "Landscape Fragmentation Metrics",
       subtitle = "Computed from OSM green patch vectors", x = NULL, y = NULL) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_ROOT, "fig_fragmentation_metrics.png"), p_frag, width = 10, height = 7, dpi = 300)

# Figure 4C: Connectivity maps — betweenness centrality of patches

# plotting function
plot_connectivity_map <- function(graph_obj, bnd_sf, local_crs, city_label, show_scale = FALSE) {
  nodes <- st_transform(graph_obj$nodes, local_crs)
  bnd   <- st_transform(bnd_sf, local_crs)

  has_edges <- !is.null(graph_obj$edges_sf) && nrow(graph_obj$edges_sf) > 0

  p <- ggplot() +
    geom_sf(data = bnd, fill = COLORS$beige, colour = COLORS$grey85, linewidth = 0.5)

  if (has_edges) {
    edges <- st_transform(graph_obj$edges_sf, local_crs)

    p <- p +
      geom_sf( data = edges, aes(colour = "Dispersal edge"), linewidth = 0.6, alpha = 0.8)}

  p <- p +
    geom_sf( data = nodes, aes(fill = area_ha), shape = 21, size = 3, colour = COLORS$green_dark, linewidth = 0.25, alpha = 0.9) +
    scale_colour_manual( name = NULL, values = c("Dispersal edge" = COLORS$orange)) +
    scale_fill_gradientn(colours = pal_green, name = "Patch size (ha)", na.value = COLORS$grey85) +
    guides(fill = guide_legend(title = "Patch size (ha)"), colour = guide_legend(override.aes = list(linewidth = 1.5))) +
    coord_sf(expand = FALSE, datum = NA) +
    theme_map_clean() +
    labs( title = city_label, subtitle = sprintf("Dispersal threshold: %d m | Nodes coloured by patch area", DISPERSAL_THRESH_M))

  if (show_scale) {
    p <- p +
      annotation_scale(location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4)}
p}

# Create plots
p_conn_yx <- plot_connectivity_map(
  yx_graph, d$yx_bnd, CRS_YX, "Yuexiu", show_scale = TRUE) +
  guides(fill = guide_legend(title = "Patch size (ha)\nYuexiu", order = 1), colour = guide_legend(title = "Network connection", order = 2))

p_conn_dl <- plot_connectivity_map(dl_graph, dl$dl_bnd, CRS_DELFT, "Delft", show_scale = FALSE) +
  guides(fill = guide_legend(title = "Patch size (ha)\nDelft", order = 1), colour = guide_legend(title = "Network connection", order = 2))

# Calculate real map widths
yx_conn_bb <- st_bbox(st_transform(d$yx_bnd, CRS_YX))
dl_conn_bb <- st_bbox(st_transform(dl$dl_bnd, CRS_DELFT))

yx_conn_w <- yx_conn_bb$xmax - yx_conn_bb$xmin
dl_conn_w <- dl_conn_bb$xmax - dl_conn_bb$xmin

# Combine
fig_connectivity_maps <- p_conn_yx + p_conn_dl +
  plot_layout(widths = c(yx_conn_w, dl_conn_w)) +
  plot_annotation(title = "Green Patch Network Connectivity", theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "right")

fig_connectivity_maps

ggsave(file.path(OUT_ROOT, "fig_connectivity_maps.png"), fig_connectivity_maps, width = 14, height = 7, dpi = 300)

# Figure 4D: Patch size vs. betweenness scatter
conn_df <- bind_rows(
  yx_graph$nodes |> st_drop_geometry() |>
    transmute(city = "Yuexiu", area_ha, betweenness, degree),
  dl_graph$nodes |> st_drop_geometry() |>
    transmute(city = "Delft",  area_ha, betweenness, degree))

p_bw <- ggplot(conn_df, aes(x = area_ha, y = betweenness, colour = city)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_x_log10(labels = label_comma()) +
  scale_colour_manual(values = c("Yuexiu" = COLORS$orange, "Delft" = COLORS$blue)) +
  theme_minimal(base_size = 12) +
  labs(title = "Patch Area vs. Betweenness Centrality", subtitle = "High betweenness = stepping-stone patch critical for connectivity", x = "Patch area (ha, log)", y = "Betweenness centrality (normalised)", colour = "City") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_ROOT, "fig_betweenness_vs_area.png"), p_bw, width = 9, height = 6, dpi = 300)

message("Script 05 complete — connectivity figures saved.")
