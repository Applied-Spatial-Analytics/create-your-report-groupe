## Data Requirements #1

| Analysis | Method | Types of data | Data source | SQ |
| :--- | :--- | :--- | :--- | :--- |
| **Accessibility** | | | | |
| Green space per capita | sum polygon area ÷ Zonal statistics population per neighbourhood | Green space polygons, population grid | [cite_start]OSM, WorldPop / GHS-POP, CBS (Delft) [cite: 2] | [cite_start]SQ1 SQ3 [cite: 2] |
| Buffer/catchment analysis | 300 m & 500 m Euclidean buffer; % population within buffer | Green polygons, population grid, admin boundaries | [cite_start]OSM, WorldPop, GADM [cite: 2] | [cite_start]SQ1 SQ4 [cite: 2] |
| Network walking isochrones | Road/path network analysis - 5, 10, 15 min isochrones from park entrances | Street network, green polygon entrances | [cite_start]OSM (QGIS network analysis) [cite: 2] | [cite_start]SQ1 [cite: 2] |
| Mean nearest park distance | Near tool-centroid of each population cell to nearest green polygon | Population grid, green polygons | [cite_start]WorldPop, OSM [cite: 2] | [cite_start]SQ1 SQ3 [cite: 2] |
| **Typology and Biodiversity** | | | | |
| Green space typology | Tag-based classification of OSM polygons into types | OSM landuse / leisure attributes | [cite_start]OSM (Overpass API) [cite: 2] | [cite_start]SQ2 [cite: 2] |
| NDVI vegetation density | Raster band calculation (NIR - Red) / (NIR + Red); zonal stats per patch | Multispectral satellite imagery (10 m) | [cite_start]Sentinel-2/gscloud.cn (China) [cite: 2] | [cite_start]SQ2 SQ4 [cite: 2] |
| Species observation density | Spatial join occurrence points to green polygons; count per ha | Biodiversity occurrence points | [cite_start]GBIF [cite: 2] | [cite_start]SQ2 [cite: 2] |
| Blue-green ratio | Separate water vs. vegetation area; ratio per neighbourhood | Land cover + OSM water polygons | [cite_start]ESA World Cover, OSM [cite: 2] | [cite_start]SQ2 SQ4 [cite: 2] |

---

## Data Requirements #2

| Analysis | Method | Types of data / Data source | SQ |
| :--- | :--- | :--- | :--- |
| **Spatial Justice** | | | |
| Green space Gini coefficient | Lorenz curve + Gini index of per-capita green space across neighbourhoods | [cite_start]Green space per capita per neighbourhood (Derived from OSM + WorldPop) [cite: 5] | [cite_start]SQ3 SQ4 [cite: 5] |
| Bivariate choropleth | Cross-tabulate green density vs. population density; $3\times3$ colour matrix | [cite_start]Admin boundaries, green area, population (GADM/CBS, OSM, World Pop) [cite: 5] | [cite_start]SQ3 [cite: 5] |
| Green space-income correlation | Pearson r - green space per capita vs. income/nighttime light proxy | [cite_start]Income data or VIIRS night lights (CBS (Delft); VIIRS (Guangzhou)) [cite: 5] | [cite_start]SQ3 [cite: 5] |
| **Connectivity and NbS** | | | |
| Fragmentation indices | Landscape metrics: NP, MPS, ENN via LecoS / FRAGSTATS (landscape ecology research tool) | [cite_start]Green space polygons + land cover raster (OSM, ESA World Cover) [cite: 5] | [cite_start]SQ2 SQ4 [cite: 5] |
| Functional connectivity graph | Graph model - nodes = patches; edges at dispersal threshold (50-200 m); betweenness centrality | [cite_start]Green polygons centroids, dispersal thresholds (literature) (OSM + NetworkX (Python)) [cite: 5] | [cite_start]SQ2 SQ5 [cite: 5] |
| NBS corridor prioritisation | MULTI CRITERIA DECISION ANALYSIS (MCDA) - : low-access zones + disconnected patches = priority intervention areas | [cite_start]All prior outputs (Derived layer) [cite: 5] | [cite_start]SQ5 [cite: 5] |