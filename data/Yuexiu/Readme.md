# Project: Urban Comparative Analysis – Yuexiu District, Guangzhou

This folder contains the processed geospatial data for the **Yuexiu District** (~33.64 km²). All data has been clipped to the district boundary and reprojected to a metric Coordinate Reference System (CRS) to ensure accurate spatial analysis and comparison with the Delft (Netherlands) dataset.

---

## 1. Vector Data (GeoPackage - .gpkg)
*Discrete features representing administrative boundaries and urban infrastructure.*

| File Name | Description | Source |
| :--- | :--- | :--- |
| **Yuexiu_boundary_proj** | The master administrative boundary of Yuexiu. Used as the clipping mask for all layers. | GADM / Local Admin |
| **Yuexiu_subdistricts_proj** | The 18 Sub-districts (*Jiēdào*) of Yuexiu. Primary unit for administrative comparison. | Civil Affairs Bureau |
| **Yuexiu_VIIRS_Points_Joined_proj** | Point cloud of Nighttime Light centroids joined with Sub-district names for zonal statistics. | NASA / Processed |
| **Yuexiu_osm_roads_proj** | Full road network including highways, primary roads, and residential streets. | OpenStreetMap |
| **Yuexiu_osm_green_proj** | Urban parks, gardens, and significant green spaces. | OpenStreetMap |
| **Yuexiu_osm_water_polygon_proj** | Standing water bodies (Lakes, ponds, and wide river basins). | OpenStreetMap |
| **Yuexiu_osm_water_proj** | Linear water features (Canals, streams, and river centerlines). | OpenStreetMap |
| **Yuexiu_gbif_proj** | Biodiversity occurrence points indicating species observations in the urban area. | GBIF.org |

---

## 2. Raster Data (GeoTIFF - .tif)
*Continuous grid data representing environmental and demographic variables.*

| File Name | Description | Source |
| :--- | :--- | :--- |
| **Yuexiu_viirs_proj** | Nighttime Light intensity (Radiance) used as a proxy for economic activity. | NOAA/NASA |
| **Yuexiu_worldpop_proj** | Gridded population estimates (people per pixel). | WorldPop |
| **Yuexiu_worldcover_proj** | 10m resolution Land Use / Land Cover (LULC) classification. | ESA WorldCover |
| **Yuexiu_ndvi_proj** | Normalized Difference Vegetation Index indicating plant health and density. | Sentinel-2 / Landsat |

---

## 3. Data Processing Workflow
1. **Projection:** All layers reprojected to a local metric CRS (e.g., UTM Zone 49N or EPSG:4490) to support `$area` calculations.
2. **Clipping:** All layers clipped to the `Yuexiu_boundary_proj` extent.
3. **Geometry Fix:** Vector layers passed through "Fix Geometries" to ensure topological integrity.
4. **Integration:** Nighttime Light rasters converted to points and spatially joined to administrative polygons for cross-regional density analysis.

---

## 4. Usage
This data is intended for a side-by-side comparison with **Delft, Netherlands**. 
- Use the `subdistricts` layer to compare with Delft `wijken`.
- Use the `VIIRS_Points` to compare urban luminous intensity.