#!/usr/bin/env python3
from pathlib import Path
import geopandas as gpd
import pandas as pd
import sys

BASE_DIR = Path(r"C:\Users\User\Desktop\avalanche_project\data\processed")
ASSUME_CRS = "EPSG:4326"

PARQUET_FILES = [
    "all_seasons.parquet",
    "avalanche_candidates_21_22.parquet",
    "avalanche_candidates_22_23.parquet",
    "avalanche_candidates_23_24.parquet",
    "mean_combined_z_cross_season_hotspot.parquet",
]


def confirm_geometry(gdf, src_name):
    if "geometry" in gdf.columns:
        try:
            gdf = gpd.GeoDataFrame(gdf, geometry="geometry")
        except Exception:
            pass
    else:
        # Try creating geometry from centroid columns
        if {"centroid_lon", "centroid_lat"}.issubset(gdf.columns):
            print(f" -> {src_name}: creating Point geometry from centroid_lon/centroid_lat")
            gdf["geometry"] = gpd.points_from_xy(
                gdf["centroid_lon"], gdf["centroid_lat"]
            )
            gdf = gpd.GeoDataFrame(gdf, geometry="geometry")
        else:
            print(f"WARNING: {src_name} has no geometry.", file=sys.stderr)
            return None

    # If CRS missing, assume EPSG:4326
    if gdf.crs is None:
        print(f" -> {src_name}: CRS missing, assuming {ASSUME_CRS}")
        gdf.set_crs(ASSUME_CRS, inplace=True)

    return gdf


def process_file(p):
    src_name = p.name

    try:
        gdf = gpd.read_parquet(p)
    except Exception as e:
        print(f"ERROR reading {src_name}: {e}", file=sys.stderr)
        return

    gdf = confirm_geometry(gdf, src_name)
    if gdf is None:
        return

    # Drop null geometries
    gdf = gdf[gdf.geometry.notna()].copy()
    if len(gdf) == 0:
        print(f"{src_name}: no valid geometries — skipping.")
        return

    # Reproject to WGS84 for GeOSON
    try:
        gdf = gdf.to_crs("EPSG:4326")
    except Exception:
        pass

    out_path = p.with_suffix(".geojson")

    try:
        gdf.to_file(out_path, driver="GeoJSON")
        print(f"WROTE: {out_path} ({len(gdf)} features)")
    except Exception as e:
        print(f"ERROR writing {out_path}: {e}", file=sys.stderr)


def main():
    for fname in PARQUET_FILES:
        p = BASE_DIR / fname
        if not p.exists():
            print(f"WARNING: file not found: {p}", file=sys.stderr)
            continue
        process_file(p)


if __name__ == "__main__":
    main()




'''
####sIMPLIFY

from pathlib import Path
import geopandas as gpd
import pandas as pd
import 


gdf['geometry'] = gdf.geometry.simplify(tolerance=0.0001)  # ~10m at this scale
gdf.to_file(r'C:\Users\User\Desktop\avalanche_project\data\processed\avalanche_inventory.geojson', driver='GeoJSON')
'''