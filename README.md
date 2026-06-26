# Spatial Justice through Urban Green Space Connectivity

Project report of Group E — Applied Spatial Analytics 2026, TU Delft.
1. Roxanne Vuijk
2. Belina Aileen
3. Daman Dogra

---

## Overview

This project analyses urban green space accessibility, biodiversity, and spatial justice across two cities — **Delft, Netherlands** and **Yuexiu, Guangzhou, China** — using open spatial datasets. The analysis addresses five sub-questions (SQ1–SQ5) spanning accessibility, typology, equity, connectivity, and multi-criteria prioritisation.

---

## Prerequisites

- [R](https://cran.r-project.org/) ≥ 4.3
- [RStudio](https://posit.co/download/rstudio-desktop/) or [VS Code](https://code.visualstudio.com/) with the Quarto extension
- [Quarto CLI](https://quarto.org/docs/get-started/) ≥ 1.4

Required R packages (install once):

```r
install.packages(c(
  "sf", "terra", "exactextractr", "dplyr", "tidyr",
  "ggplot2", "patchwork", "scales", "ggrepel",
  "ggspatial", "ineq", "here"
))
```

---

## Data

Data files are **not tracked in this repo**. Place them at the paths expected by `R/00_config.R`.

Download the data folder from: [Google Drive](https://drive.google.com/file/d/1Pk3aUc24ta8mHQLxxRYyFL_q4feXcNhA/view?usp=drive_link)

Expected structure after download:

```
data/
├── delft/
│   ├── raster/    # delft_worldpop_proj.tif, ndvi_delft_proj.tif, …
│   └── vector/    # delft_boundary_proj.gpkg, delft_osm_green_proj.gpkg, …
│   └── README.md
└── Yuexiu/
    ├── raster/    # Yuexiu_viirs_proj.tif, …
    └── vector/    # yuexiu_boundary.gpkg, guangzhou_osm_green_proj.gpkg, …
    └── README.md
```

> All paths are centralised in `R/00_config.R` — edit that file if your data lives elsewhere.

---

## How to run

### 1. Clone and open

```bash
git clone https://github.com/Applied-Spatial-Analytics/create-your-report-groupe.git
```

Open `asa2025-report.Rproj` in RStudio, or open the folder in VS Code.

### 2. Run the full pipeline

```r
source("run_pipeline.R")
```

This sources all scripts in order (`R/01_load_data.R` → `R/07_label_reference_mcda_map.R`). Intermediate `.rds` objects and all `fig_*.png` figures are saved automatically to `outputs/` (path controlled by `OUT_ROOT` in `R/00_config.R`).

### 3. Render the report

In RStudio: open `report.qmd` and click **Render**.
In VS Code / terminal:

```bash
quarto render
```

Output goes to `docs/index.html` (set via `output-dir: docs` in `_quarto.yml`). This is the folder GitHub Pages serves the live report from — no manual copying needed.

---

## File structure

```
create-your-report-groupe/
│
├── data/                               # not tracked — see Data section above
│   ├── delft/
│   │   ├── raster/                     # worldpop, ndvi, worldcover, viirs
│   │   └── vector/                     # boundary, wijken, buurten, income, OSM green/roads/water, GBIF
│   └── Yuexiu/
│       ├── raster/                     # worldpop, ndvi, worldcover, viirs
│       └── vector/                     # boundary, subdistricts, OSM green/roads/water, GBIF, VIIRS points
│
├── R/
│   ├── 00_config.R                     # All paths, CRS constants, buffer distances, MCDA weights, colour palettes
│   ├── 01_load_data.R                  # Load & validate all vector + raster layers; saves yuexiu_data.rds / delft_data.rds
│   ├── 02_accessibility.R              # SQ1 — green space per capita, 300/500 m buffers, nearest-green distance
│   ├── 03_typology_biodiversity.R      # SQ2 — OSM green typology, NDVI zonal stats, GBIF species density
│   ├── 04_spatial_justice.R            # SQ3 — Lorenz curve, Gini coefficient, bivariate choropleth, income/VIIRS correlation
│   ├── 05_connectivity.R               # SQ4 — fragmentation metrics (NP, MPS, ENN), graph connectivity, betweenness
│   ├── 06_mcda_nbs.R                   # SQ5 — MCDA scoring, NbS corridor prioritisation
│   ├── 07_label_reference_mcda_map.R   # Adds transliterated district labels to MCDA maps (Yuexiu)
│   └── 08_context.R                    # Context/overview maps for report introduction
│
├── outputs/                            # Auto-generated: .rds intermediates + fig_*.png figures
│
├── docs/                               # Rendered report served by GitHub Pages
│   ├── index.html                      # Rendered Quarto report
│   ├── outputs/                        # Figure PNGs referenced by the HTML report
│   └── report_files/libs/              # JS/CSS assets (Bootstrap, Quarto HTML)
│
├── pages/                              # Split pages for future multi-page website version
│
├── run_pipeline.R                      # Sources all R/ scripts in order
├── _quarto.yml                         # Quarto project config (output-dir: docs)
├── report.qmd                          # Main report narrative — pulls figures from outputs/
├── report.qmd.bak                      # Backup of original single-page report
├── asa2025-report.Rproj                # RStudio project file
├── references.bib                      # BibTeX bibliography
├── methodology.md                      # Extended methodology notes
└── README.md
```

---

## Outputs

All figures are written to `outputs/` by the pipeline scripts and copied to `docs/outputs/` on render:

| Figure | Script | Content |
|---|---|---|
| `fig_context_gz_dl.png` | 08 | Study site context map — Guangzhou and Delft |
| `fig_context_yx_dl.png` | 08 | Study site context map — Yuexiu and Delft |
| `fig_ugs_yx_dl.png` | 08 | Urban green space overview map |
| `fig_population_density.png` | 02 | Population density by subdistrict/wijk |
| `fig_access_per_capita.png` | 02 | Green space m² per capita by subdistrict/wijk |
| `fig_buffer_coverage.png` | 02 | Population within 300/500 m of green space |
| `fig_nearest_green_distance.png` | 02 | Distance to nearest green space |
| `fig_network_walk_access.png` | 02 | Walk-network accessibility map |
| `fig_green_typology.png` | 03 | OSM green space typology breakdown |
| `fig_ndvi_zonal.png` | 03 | NDVI zonal statistics by admin unit |
| `fig_ndvi_violin.png` | 03 | NDVI distribution violin plots by typology |
| `fig_gbif_density.png` | 03 | Biodiversity (GBIF) observation density per hectare |
| `fig_blue_green_ratio.png` | 03 | Blue vs. green space balance index |
| `fig_lorenz_gini.png` | 04 | Lorenz curve and Gini coefficient |
| `fig_bivariate_choropleth.png` | 04 | Bivariate map: green density × population density |
| `fig_equity_correlations.png` | 04 | Scatter plots: green access vs. socioeconomic proxy |
| `fig_enn_distribution.png` | 05 | Euclidean nearest-neighbour distance distribution |
| `fig_fragmentation_metrics.png` | 05 | NP, MPS, ENN fragmentation metrics |
| `fig_connectivity_maps.png` | 05 | Green space connectivity graph maps |
| `fig_betweenness_vs_area.png` | 05 | Patch betweenness centrality vs. area |
| `fig_mcda_maps.png` | 06 | MCDA composite urgency scores by subdistrict/wijk |
| `fig_mcda_radar.png` | 06 | City-mean MCDA sub-scores by criterion |
| `fig_priority_tiers.png` | 06 | NbS intervention priority tiers (High/Medium/Low) |
| `fig_nbs_corridors.png` | 06 | Proposed green corridors overlaid on priority zones |
| `fig_mcda_distribution.png` | 06 | Distribution of MCDA composite scores |

---

## Reproducibility notes

- All file paths live in `R/00_config.R` — no hardcoded paths elsewhere.
- Uses only open, globally available datasets (OSM, WorldPop, GBIF, VIIRS, NDVI / WorldCover).
- The workflow is adaptable to any city/district by updating `R/00_config.R`.
- The report is published via GitHub Pages from the `docs/` folder on `main`. Every render-and-push automatically updates the live version.
- Population floor artifact filter (`pop_count < 100`) is applied before all per-capita calculations to exclude WorldPop placeholder values on non-residential zones.
