# ── 04_spatial_justice.R (v2 — cleaned) ──────────────────────────────────────
# Spatial justice analysis: Gini, bivariate choropleth,
# green space vs. income / VIIRS correlation.
#
# v2 changes on top of v1 fixes:
#   8. Filter pop_count < 100 from both cities before ALL analyses.
#      These are WorldPop placeholder values (no real population signal),
#      not real residential units. Leaving them in inflates green_pc_m2
#      to thousands of m²/person and corrupts every per-capita metric.
#   9. Delft Wijk 16 (Delftse Hout) and Wijk 26 (Abtswoude) are park/polder
#      wijken with real pop but no CBS income disclosure — noted in output.
#  10. Lorenz / Gini now computed on the cleaned set for consistency.
# ─────────────────────────────────────────────────────────────────────────────

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

# ── Load data ─────────────────────────────────────────────────────────────────
d  <- readRDS(file.path(OUT_ROOT, "yuexiu_data.rds"))
dl <- readRDS(file.path(OUT_ROOT, "delft_data.rds"))


d$yx_pop   <- unwrap(d$yx_pop)
d$yx_viirs <- unwrap(d$yx_viirs)
dl$dl_pop  <- unwrap(dl$dl_pop)

# ── FIX :
# Threshold catches Yuexiu's pop_count == 1 jiēdào (genuine WorldPop floor
# artifacts on parks/water/non-residential zones — confirmed via sorted
# pop_count inspection: next value above 1 jumps to 28,253, a hard cliff).
#
# NOTE: Delft's equivalent pop==1 rows (sliver polygons from neighbouring
# municipalities, e.g. Schipluiden/Delfgauw) are now fixed upstream in
# script 01 via the gemeentecode == "GM0503" filter, applied immediately
# after dl_wijk is read. This POP_MIN filter is no longer acting as a
# backstop for Delft — there should be zero rows left for it to catch
# on that side. It remains necessary for Yuexiu.

POP_MIN <- 100

yx_sub_access_all  <- yx_sub_access
dl_wijk_access_all <- dl_wijk_access

yx_sub_access  <- yx_sub_access  |> filter(!is.na(pop_count) & pop_count >= POP_MIN)
dl_wijk_access <- dl_wijk_access |> filter(!is.na(pop_count) & pop_count >= POP_MIN)
message(sprintf("After pop filter — Yuexiu: %d subdistricts | Delft: %d wijken",
                nrow(yx_sub_access), nrow(dl_wijk_access)))

# ── 3A. Gini coefficient ──────────────────────────────────────────────────────
gini_index <- function(x) {
  x <- x[!is.na(x) & x >= 0]
  if (length(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

lorenz_df <- function(x, city_label) {
  x <- sort(x[!is.na(x) & x >= 0])
  cumx <- cumsum(x) / sum(x)
  data.frame(
    city    = city_label,
    cum_pop = seq_along(x) / length(x),
    cum_grn = cumx
  )
}

yx_gini <- gini_index(yx_sub_access$green_pc_m2)
dl_gini <- gini_index(dl_wijk_access$green_pc_m2)

message(sprintf("Gini — Yuexiu: %.3f  |  Delft: %.3f", yx_gini, dl_gini))

lorenz_data <- bind_rows(
  lorenz_df(yx_sub_access$green_pc_m2, "Yuexiu"),
  lorenz_df(dl_wijk_access$green_pc_m2, "Delft")
)
gini_labels <- data.frame(
  city  = c("Yuexiu", "Delft"),
  gini  = c(yx_gini, dl_gini),
  label = sprintf("Gini = %.3f", c(yx_gini, dl_gini))
)

# ── 3B. Bivariate choropleth ──────────────────────────────────────────────────

bivar_palette <- c(
  "Low-Low"   = COLORS$pink_light,
  "Low-Mid"   = COLORS$pink,
  "Low-High"  = COLORS$red_light,
  "Mid-Low"   = COLORS$beige,
  "Mid-Mid"   = COLORS$blue_light,
  "Mid-High"  = COLORS$blue,
  "High-Low"  = COLORS$green_light,
  "High-Mid"  = COLORS$green_mid,
  "High-High" = COLORS$green_dark
)

bivar_order <- c(
  "Low-Low", "Low-Mid", "Low-High",
  "Mid-Low", "Mid-Mid", "Mid-High",
  "High-Low", "High-Mid", "High-High"
)

tertile_label <- function(x, labels = c("Low", "Mid", "High")) {
  x_num <- as.numeric(x)
  if (all(is.na(x_num))) {
    return(factor(rep("Mid", length(x_num)), levels = labels))
  }
  brks <- quantile(x_num, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if (length(unique(brks)) < 4) {
    return(factor(ifelse(is.na(x_num), "Mid", "Mid"), levels = labels))
  }
  cut(x_num, breaks = brks, labels = labels, include.lowest = TRUE)
}

yx_bivar <- yx_sub_access |>
  mutate(
    admin_area_m2 = as.numeric(st_area(st_transform(yx_sub_access, CRS_YX))),
    green_density = green_area_m2 / admin_area_m2,
    pop_density   = pop_count / (admin_area_m2 / 10000),
    grn_class     = tertile_label(green_density),
    pop_class     = tertile_label(pop_density),
    bivar_class   = factor(paste(grn_class, pop_class, sep = "-"),
                           levels = names(bivar_palette))
  )

dl_bivar <- dl_wijk_access |>
  mutate(
    admin_area_m2 = as.numeric(st_area(st_transform(dl_wijk_access, CRS_DELFT))),
    green_density = green_area_m2 / admin_area_m2,
    pop_density   = pop_count / (admin_area_m2 / 10000),
    grn_class     = tertile_label(green_density),
    pop_class     = tertile_label(pop_density),
    bivar_class   = factor(paste(grn_class, pop_class, sep = "-"),
                           levels = names(bivar_palette))
  )

# ── 3C. Income / VIIRS correlations ──────────────────────────────────────────

# --- Delft income (buurt → wijk aggregation) ---
dl_inc_sf <- dl$dl_inc
INC_COL   <- "mean_income_per_resident_x1000eur"

if (!INC_COL %in% names(dl_inc_sf)) {
  cbs <- c("gemiddeldInkomenPerInwoner",
           "gemiddeldInkomenPerInkomensontvanger",
           "gemiddeldGestandaardiseerdInkomenVanHuishoudens")
  INC_COL <- cbs[cbs %in% names(dl_inc_sf)][1]
  if (is.na(INC_COL)) stop("No income column found in dl_inc.")
  dl_inc_sf[[INC_COL]] <- ifelse(dl_inc_sf[[INC_COL]] < 0, NA_real_, dl_inc_sf[[INC_COL]])
  message("Fell back to CBS column: ", INC_COL)
} else {
  message("Using income column: ", INC_COL)
}

dl_inc_wijk <- st_drop_geometry(dl_inc_sf) |>
  filter(!is.na(wijkcode)) |>
  group_by(wijkcode) |>
  summarise(income_wijk = mean(.data[[INC_COL]], na.rm = TRUE), .groups = "drop") |>
  mutate(income_wijk = ifelse(is.nan(income_wijk), NA_real_, income_wijk))

dl_inc_joined <- dl_wijk_access |>
  left_join(dl_inc_wijk, by = "wijkcode")

n_income <- sum(!is.na(dl_inc_joined$income_wijk))
message(sprintf("Delft: %d / %d wijken have income data after pop filter",
                n_income, nrow(dl_inc_joined)))

# Note park/polder wijken with real pop but suppressed income
park_wijken <- dl_inc_joined |>
  st_drop_geometry() |>
  filter(is.na(income_wijk)) |>
  select(wijknaam, pop_count, green_pc_m2)
if (nrow(park_wijken) > 0) {
  message("Wijken with population but no income disclosure (likely park/polder):")
  print(park_wijken)
}

if (n_income < 3) stop("Too few Delft wijken with income data (n = ", n_income, ") — check CBS join.")
dl_corr <- cor(dl_inc_joined$green_pc_m2, dl_inc_joined$income_wijk,
               use = "complete.obs", method = "pearson")
message(sprintf("Delft Pearson r (green pc ~ income): %.3f  [n = %d]",
                dl_corr, n_income))

# --- Yuexiu VIIRS ---
yx_viirs_rast <- d$yx_viirs
yx_sub_reproj <- st_transform(yx_sub_access, terra::crs(yx_viirs_rast))

yx_viirs_zonal <- yx_sub_access |>
  mutate(viirs_mean = exact_extract(yx_viirs_rast, yx_sub_reproj, "mean"))

yx_corr <- cor(yx_viirs_zonal$green_pc_m2, yx_viirs_zonal$viirs_mean,
               use = "complete.obs", method = "pearson")
message(sprintf("Yuexiu Pearson r (green pc ~ VIIRS): %.3f  [n = %d]",
                yx_corr, sum(!is.na(yx_viirs_zonal$viirs_mean))))

# # --- Delft VIIRS ---

# # Fetch from the list and immediately unwrap it
# dl_viirs_raster_final <- unwrap(dl$dl_viirs)

# # Validate that it is a SpatRaster
# if (!inherits(dl_viirs_raster_final, "SpatRaster")) {
#   stop("The object is still not a SpatRaster. Please check dl$dl_viirs content.")
# }

# # Set CRS if missing
# if (is.na(terra::crs(dl_viirs_raster_final))) {
#   terra::crs(dl_viirs_raster_final) <- "EPSG:4326"
# }

# # Project the polygons
# dl_sub_reproj <- st_transform(dl_wijk_access, terra::crs(dl_viirs_raster_final))

# # Perform the extraction
# dl_wijk_access$viirs_mean <- exact_extract(dl_viirs_raster_final, dl_sub_reproj, "mean")
# # 1. Access the raster object from the 'dl' list, unwrap it, and check it
# dl_viirs_rast <- unwrap(dl$dl_viirs)

# if(is.na(terra::crs(dl_viirs_rast))) {
#   terra::crs(dl_viirs_rast) <- "EPSG:4326"
# }

# # 2. Project your vector data to match the raster's CRS
# dl_sub_reproj <- st_transform(dl_wijk_access, terra::crs(dl_viirs_rast))

# # 3. Now run the extraction on the confirmed SpatRaster
# dl_wijk_access$viirs_mean <- exact_extract(dl_viirs_rast, dl_sub_reproj, "mean")

# # 4. Calculate correlation
# dl_viirs_zonal <- dl_wijk_access
# dl_corr <- cor(dl_viirs_zonal$green_pc_m2, dl_viirs_zonal$viirs_mean,
#                use = "complete.obs", method = "pearson")

# message(sprintf("Delft Pearson r (green pc ~ VIIRS): %.3f  [n = %d]",
#                 dl_corr, sum(!is.na(dl_viirs_zonal$viirs_mean))))

# ── Save outputs ──────────────────────────────────────────────────────────────
saveRDS(lorenz_data,    file.path(OUT_ROOT, "lorenz_data.rds"))
saveRDS(gini_labels,    file.path(OUT_ROOT, "gini_labels.rds"))
saveRDS(yx_bivar,       file.path(OUT_ROOT, "yx_bivar.rds"))
saveRDS(dl_bivar,       file.path(OUT_ROOT, "dl_bivar.rds"))
saveRDS(dl_inc_joined,  file.path(OUT_ROOT, "dl_inc_joined.rds"))
saveRDS(yx_viirs_zonal, file.path(OUT_ROOT, "yx_viirs_zonal.rds"))

corr_summary <- data.frame(
  city      = c("Delft", "Yuexiu"),
  proxy     = c("Household income (CBS)", "Nighttime light (VIIRS)"),
  pearson_r = c(dl_corr, yx_corr),
  n         = c(n_income, sum(!is.na(yx_viirs_zonal$viirs_mean)))
)
saveRDS(corr_summary, file.path(OUT_ROOT, "corr_summary.rds"))

# ── Figures ───────────────────────────────────────────────────────────────────

# Figure 3A: Lorenz curves
p_lorenz <- ggplot(lorenz_data, aes(x = cum_pop, y = cum_grn, colour = city)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text(data = gini_labels,
            aes(x = 0.25, y = c(0.82, 0.72), label = label, colour = city),
            size = 4, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = c("Yuexiu" = COLORS$red, "Delft" = COLORS$blue)) +
  scale_x_continuous(labels = label_percent()) +
  scale_y_continuous(labels = label_percent()) +
  theme_minimal(base_size = 12) +
  labs(title = "Lorenz Curve — Green Space per Capita Inequality",
       subtitle = "Dashed diagonal = perfect equality; further below = less equal",
       x = "Cumulative share of neighbourhoods (ranked by green pc)",
       y = "Cumulative share of total green space",
       colour = "City") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_ROOT, "fig_lorenz_gini.png"),
       p_lorenz, width = 8, height = 6, dpi = 300)

# Figure 3B: Bivariate choropleth — single shared legend
# Key fix: convert bivar_class to a factor with IDENTICAL levels in both
# sf objects. patchwork's guides = "collect" only merges legends when the
# underlying scale (values + breaks + name) is byte-for-byte identical.
# Sorting + factoring guarantees that.
all_classes <- sort(names(bivar_palette))   # all 9 possible cells, always same order

yx_bivar <- yx_bivar |>
  mutate(bivar_class = factor(paste(grn_class, pop_class, sep = "-"),
                     levels = bivar_order))
dl_bivar <- dl_bivar |>
  mutate(bivar_class = factor(paste(grn_class, pop_class, sep = "-"),
                     levels = bivar_order))

# 2. Transform both maps to metre-based CRS
yx_plot <- st_transform(yx_bivar, CRS_YX)
dl_plot <- st_transform(dl_bivar, 28992)

# 3. Calculate real map widths
yx_bb <- st_bbox(yx_plot)
dl_bb <- st_bbox(dl_plot)

yx_w <- yx_bb$xmax - yx_bb$xmin
dl_w <- dl_bb$xmax - dl_bb$xmin

# 4. Shared legend scale
shared_fill_no_legend <- scale_fill_manual(
  values = bivar_palette,
  breaks = bivar_order,
  limits = bivar_order,
  drop = FALSE,
  guide = "none"
)

# 3×3 legend tile — legend_df must be defined first
legend_df <- expand.grid(grn_class = c("Low","Mid","High"),
                         pop_class = c("Low","Mid","High")) |>
  mutate(bivar_class = paste(grn_class, pop_class, sep = "-"),
         grn_class = factor(grn_class, c("Low","Mid","High")),
         pop_class = factor(pop_class, c("Low","Mid","High")))

p_legend <- ggplot(legend_df, aes(x = grn_class, y = pop_class, fill = bivar_class)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = bivar_palette, guide = "none") +
  coord_fixed(clip = "off") +
  theme_minimal(base_size = 9) +
  labs(
    title = "Bivariate legend",
    x = "← Green density →",
    y = "← Population density →"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )


p_bv_yx <- ggplot(yx_bivar) +
  geom_sf(aes(fill = bivar_class), colour = "white", linewidth = 0.3) +
  annotation_scale(location = "bl", style = "ticks", width_hint = 0.25, text_cex = 0.7, line_width = 0.4) +
  shared_fill_no_legend +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Yuexiu", subtitle = "3×3 tertile classification")

p_bv_dl <- ggplot(dl_bivar) +
  geom_sf(aes(fill = bivar_class), colour = "white", linewidth = 0.3) +
  shared_fill_no_legend +
  coord_sf(expand = FALSE, datum = NA) +
  theme_map_clean() +
  labs(title = "Delft", subtitle = "3×3 tertile classification")

# give legend more space
fig_bivariate_choropleth <- p_bv_yx + p_bv_dl + p_legend +
  plot_layout(
    widths = c(yx_w, dl_w, max(yx_w, dl_w) * 0.35)
  ) +
  plot_annotation(
    title = "Bivariate: Green Density × Population Density",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  )

fig_bivariate_choropleth

ggsave(file.path(OUT_ROOT, "fig_bivariate_choropleth.png"), fig_bivariate_choropleth, width = 14, height = 6, dpi = 300)

ggsave(file.path(OUT_ROOT, "fig_bivariate_legend.png"),
       p_legend, width = 3.5, height = 3.5, dpi = 300)

# Figure 3C: Scatter plots
p_corr_dl <- ggplot(dl_inc_joined |> st_drop_geometry(),
                    aes(x = income_wijk, y = green_pc_m2)) +
  geom_point(colour = COLORS$blue, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = COLORS$blue,
              fill = COLORS$blue, alpha = 0.15) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.3f\n(n = %d)", dl_corr, n_income),
           hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold") +
  scale_x_continuous(labels = label_comma(suffix = "k€")) +
  scale_y_continuous(labels = label_comma()) +
  theme_minimal(base_size = 12) +
  labs(title = "Green Space per Capita vs. Household Income — Delft",
       subtitle = "Pearson correlation; wijk level (residential wijken only)",
       x = "Mean disposable income (€1 000 / resident)",
       y = "Green space (m² per person)") +
  theme(plot.title = element_text(face = "bold"))

n_yx <- sum(!is.na(yx_viirs_zonal$viirs_mean))
p_corr_yx <- ggplot(yx_viirs_zonal |> st_drop_geometry(),
                    aes(x = viirs_mean, y = green_pc_m2)) +
  geom_point(colour = COLORS$orange, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = COLORS$orange,
              fill = COLORS$orange, alpha = 0.15) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("r = %.3f\n(n = %d)", yx_corr, n_yx),
           hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold") +
  scale_y_continuous(labels = label_comma()) +
  theme_minimal(base_size = 12) +
  labs(title = "Green Space per Capita vs. VIIRS Nighttime Light — Yuexiu",
       subtitle = "Pearson correlation; subdistrict level (residential only)",
       x = "Mean VIIRS radiance (nW/cm²/sr)",
       y = "Green space (m² per person)") +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_ROOT, "fig_equity_correlations.png"),
       p_corr_dl + p_corr_yx, width = 14, height = 6, dpi = 300)

message("Script 04 v2 complete — spatial justice figures saved.")

