# Delft, Netherlands

This folder contains the processed geospatial data for the **Municipality of Delft** (~24.06 km²). All layers have been clipped to the municipal boundary and reprojected to the local Dutch coordinate system to ensure geometric accuracy for comparison with the Yuexiu District (Guangzhou) dataset.

---

## 1. Vector Data (GeoPackage - .gpkg)
*Discrete features representing administrative divisions and socioeconomic infrastructure.*

| File Name | Description | Source |
| :--- | :--- | :--- |
| **delft_boundary_proj** | The official municipal boundary of Delft. Master mask used for all data clipping. | PDOK / Basisregistratie Gemeenten |
| **delft_buurten_proj** | The 91 Neighborhoods (*Buurten*) of Delft. Provides high-resolution statistical granularity. | CBS (Statistics Netherlands) |
| **delft_income_proj** | Socioeconomic data representing the average disposable income per household/district. | CBS StatLine |
| **delft_osm_roads_proj** | Complete road network, including the historic center's pedestrian paths and cycling infrastructure. | OpenStreetMap |
| **delft_gbif_proj** | Biodiversity occurrence points indicating recorded species within city limits. | GBIF.org |

---

## 2. Raster Data (GeoTIFF - .tif)
*Continuous grid data representing environmental and demographic variables.*

| File Name | Description | Source |
| :--- | :--- | :--- |
| **delft_worldpop_proj** | Gridded population density estimates (people per pixel). | WorldPop |
| **delft_worldcover_proj** | 10m resolution Land Use / Land Cover classification (ESA WorldCover). | European Space Agency (ESA) |
| **ndvi_delft_proj** | Normalized Difference Vegetation Index indicating urban biomass and green health. | Sentinel-2 / Landsat |

---

## 3. Data Processing Workflow
1. **Projection:** All layers reprojected to **EPSG:28992** (Amersfoort / RD New), the official projected coordinate system for the Netherlands, to ensure sub-meter accuracy.
2. **Clipping:** All layers clipped to the `delft_boundary_proj` extent to remove peripheral data from neighboring municipalities (e.g., Rijswijk, Schiedam).
3. **Alignment:** Demographic and socioeconomic vectors are aligned with the 2025 CBS neighborhood boundaries for consistency.

---

## 4. Usage & Comparison
This dataset is designed for a direct "scale-to-scale" comparison with **Yuexiu, Guangzhou**.
- **Wijken/Buurten vs. Jiedao:** Note that Delft’s 91 buurten offer a higher spatial resolution than the 18 Jiedao in Yuexiu, making them more comparable to Yuexiu's *Shèqū* level.
- **Income vs. Nightlight:** Used to evaluate the correlation between household wealth in Delft and luminous intensity (economic proxy) in Yuexiu.
