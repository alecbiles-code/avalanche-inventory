import geopandas as gpd
import pandas as pd
import os

# CONFIG 
INPUT_DIR  = r"C:\Users\User\Desktop\avalanche_project"
OUTPUT_DIR = r"C:\Users\User\Desktop\avalanche_project\data\processed"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Seasonal shapefiles
SEASONAL = [
    ("combined_z_anomaly_21_22_edit.shp", "2021-2022", "avalanche_candidates_21_22"),
    ("combined_z_anomaly_22_23_edit.shp", "2022-2023", "avalanche_candidates_22_23"),
    ("combined_z_anomaly_23_24_edit.shp", "2023-2024", "avalanche_candidates_23_24"),
]

HOTSPOT_SHP  = "mean_combined_z_cross_season_hotspot_edit.shp"
HOTSPOT_OUT  = "mean_combined_z_cross_season_hotspot"

#  clean and standardize a GeoDataFrame 
def clean_gdf(gdf, season_label=None, is_hotspot=False):
    # Reproject to WGS84  required for GeoParquet interoperability
    gdf = gdf.to_crs("EPSG:4326")

    # Keep only what matters  drop ArcGIS noise columns
    keep_cols = ["geometry", "confidence", "Area_m2"]
    gdf = gdf[[c for c in keep_cols if c in gdf.columns]]

    # Add season label
    if season_label:
        gdf["season"] = season_label

    if is_hotspot:
        gdf["type"] = "cross_season_hotspot"
    else:
        gdf["type"] = "seasonal_candidate"

    # Add centroid coordinate for DuckDB spatial queries later
    centroids = gdf.geometry.centroid
    gdf["centroid_lon"] = centroids.x
    gdf["centroid_lat"] = centroids.y

    # Add unique ID
    prefix = season_label.replace("-", "_") if season_label else "hotspot"
    gdf = gdf.reset_index(drop=True)
    gdf["event_id"] = [f"AVA_{prefix}_{i:05d}" for i in gdf.index]

    # Drop null geometries
    gdf = gdf[gdf.geometry.notna()]

    return gdf

#  PROCESS SEASONAL FILES 
seasonal_gdfs = []

for shp_name, season_label, out_name in SEASONAL:
    shp_path = os.path.join(INPUT_DIR, shp_name)
    out_path  = os.path.join(OUTPUT_DIR, f"{out_name}.parquet")

    print(f"\nProcessing: {shp_name}")
    gdf = gpd.read_file(shp_path)
    print(f"  Features loaded: {len(gdf)}")
    print(f"  CRS: {gdf.crs}")

    gdf = clean_gdf(gdf, season_label=season_label)

    gdf.to_parquet(out_path, index=False)
    print(f"  Saved → {out_path}")
    print(f"  Confidence breakdown:\n{gdf['confidence'].value_counts().to_string()}")

    seasonal_gdfs.append(gdf)

#  PROCESS HOTSPOT FILE 
hotspot_path = os.path.join(INPUT_DIR, HOTSPOT_SHP)
hotspot_out  = os.path.join(OUTPUT_DIR, f"{HOTSPOT_OUT}.parquet")

print(f"\nProcessing: {HOTSPOT_SHP}")
gdf_hotspot = gpd.read_file(hotspot_path)
print(f"  Features loaded: {len(gdf_hotspot)}")

gdf_hotspot = clean_gdf(gdf_hotspot, is_hotspot=True)
gdf_hotspot.to_parquet(hotspot_out, index=False)
print(f"  Saved → {hotspot_out}")

#  COMBINE ALL SEASONS INTO ONE FILE
print(f"\nBuilding all_seasons combined file...")
all_seasons = gpd.GeoDataFrame(
    pd.concat(seasonal_gdfs, ignore_index=True),
    crs="EPSG:4326"
)
all_out = os.path.join(OUTPUT_DIR, "all_seasons.parquet")
all_seasons.to_parquet(all_out, index=False)
print(f"  Total features: {len(all_seasons)}")
print(f"  Saved → {all_out}")

# SUMMARY
print("\n" + "=" * 55)
print("GeoParquet export complete.")
print(f"Output directory: {OUTPUT_DIR}")
print()
print("Files created:")
for _, _, out_name in SEASONAL:
    print(f"  {out_name}.parquet")
print(f"  {HOTSPOT_OUT}.parquet")
print(f"  all_seasons.parquet")
print()
print("Season summary:")
for gdf, (_, label, _) in zip(seasonal_gdfs, SEASONAL):
    hc = (gdf['confidence'] == 'higher confidence').sum()
    lc = (gdf['confidence'] == 'lower confidence').sum()
    print(f"  {label}: {len(gdf)} total | {hc} higher confidence | {lc} lower confidence")
print("=" * 55)