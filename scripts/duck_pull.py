import duckdb
import os

# CONFIG
DATA_DIR = r"C:\Users\User\Desktop\avalanche_project\data\processed"

PARQUET_FILES = {
    "all_seasons"       : os.path.join(DATA_DIR, "all_seasons.parquet"),
    "season_21_22"      : os.path.join(DATA_DIR, "avalanche_candidates_21_22.parquet"),
    "season_22_23"      : os.path.join(DATA_DIR, "avalanche_candidates_22_23.parquet"),
    "season_23_24"      : os.path.join(DATA_DIR, "avalanche_candidates_23_24.parquet"),
    "cross_season"      : os.path.join(DATA_DIR, "mean_combined_z_cross_season_hotspot.parquet"),
}

DILLON_RESERVOIR = (39.61, -106.06)   # lat, lon
RADIUS_M         = 50000              # 50km

#  INITIALIZEEE
con = duckdb.connect()
con.execute("INSTALL spatial;")
con.execute("LOAD spatial;")

print("=" * 60)
print("DuckDB Avalanche Inventory — Spatial Query Validation")
print(f"Spatial filter: {RADIUS_M/1000:.0f}km radius from Dillon Reservoir")
print(f"  ({DILLON_RESERVOIR[0]}°N, {abs(DILLON_RESERVOIR[1])}°W)")
print("=" * 60)

#QUERY 1: Season summary across all files
print("\n── Query 1: Feature count and area summary per file ──\n")

for name, path in PARQUET_FILES.items():
    result = con.execute(f"""
        SELECT
            COUNT(*)                        AS total_features,
            ROUND(AVG(Area_m2), 1)          AS avg_area_m2,
            ROUND(SUM(Area_m2) / 10000, 2)  AS total_area_ha,
            COUNT(*) FILTER (WHERE confidence = 'higher confidence') AS higher_confidence,
            COUNT(*) FILTER (WHERE confidence = 'lower confidence')  AS lower_confidence
        FROM read_parquet('{path}')
    """).df()
    print(f"  {name}")
    print(result.to_string(index=False))
    print()

# QUERY 2: Spatial filter — within 50km of Dillon Reservoir
print("\n── Query 2: Features within 50km of Dillon Reservoir ──\n")

for name, path in PARQUET_FILES.items():
    season_col = "season," if name != "cross_season" else ""
    result = con.execute(f"""
        SELECT
            event_id,
            {season_col}
            confidence,
            ROUND(Area_m2, 1)       AS area_m2,
            ROUND(centroid_lat, 5)  AS centroid_lat,
            ROUND(centroid_lon, 5)  AS centroid_lon,
            ROUND(ST_Distance_Sphere(
                ST_Point(centroid_lon, centroid_lat),
                ST_Point({DILLON_RESERVOIR[1]}, {DILLON_RESERVOIR[0]})
            ) / 1000, 2)            AS distance_km
        FROM read_parquet('{path}')
        WHERE ST_Distance_Sphere(
            ST_Point(centroid_lon, centroid_lat),
            ST_Point({DILLON_RESERVOIR[1]}, {DILLON_RESERVOIR[0]})
        ) <= {RADIUS_M}
        ORDER BY distance_km ASC
    """).df()
    print(f"  {name}: {len(result)} features within {RADIUS_M/1000:.0f}km")
    if len(result) > 0:
        print(result.to_string(index=False))
    print()

#─── QUERY 3: Confidence breakdown near Dillon resy
print("\n── Query 3: Confidence breakdown within 50km of Dillon Reservoir ──\n")

result = con.execute(f"""
    SELECT
        season,
        confidence,
        COUNT(*)                       AS feature_count,
        ROUND(AVG(Area_m2), 1)         AS avg_area_m2,
        ROUND(SUM(Area_m2) / 10000, 2) AS total_area_ha
    FROM read_parquet('{PARQUET_FILES["all_seasons"]}')
    WHERE ST_Distance_Sphere(
        ST_Point(centroid_lon, centroid_lat),
        ST_Point({DILLON_RESERVOIR[1]}, {DILLON_RESERVOIR[0]})
    ) <= {RADIUS_M}
    GROUP BY season, confidence
    ORDER BY season, confidence
""").df()

print(result.to_string(index=False))

# QUERY 4: Largest features near Dillon — all seasons
print("\n── Query 4: Top 10 largest features within 50km of Dillon Reservoir ──\n")

result = con.execute(f"""
    SELECT
        event_id,
        season,
        confidence,
        ROUND(Area_m2, 1)      AS area_m2,
        ROUND(ST_Distance_Sphere(
            ST_Point(centroid_lon, centroid_lat),
            ST_Point({DILLON_RESERVOIR[1]}, {DILLON_RESERVOIR[0]})
        ) / 1000, 2)           AS distance_km
    FROM read_parquet('{PARQUET_FILES["all_seasons"]}')
    WHERE ST_Distance_Sphere(
        ST_Point(centroid_lon, centroid_lat),
        ST_Point({DILLON_RESERVOIR[1]}, {DILLON_RESERVOIR[0]})
    ) <= {RADIUS_M}
    ORDER BY Area_m2 DESC
    LIMIT 10
""").df()

print(result.to_string(index=False))

#QUERY 5: Cross-season hotspot check near lake Dillon
print("\n── Query 5: Cross-season hotspot features within 50km of Dillon Reservoir ──\n")

result = con.execute(f"""
    SELECT
        event_id,
        confidence,
        ROUND(Area_m2, 1)      AS area_m2,
        ROUND(centroid_lat, 5) AS centroid_lat,
        ROUND(centroid_lon, 5) AS centroid_lon,
        ROUND(ST_Distance_Sphere(
            ST_Point(centroid_lon, centroid_lat),
            ST_Point({DILLON_RESERVOIR[1]}, {DILLON_RESERVOIR[0]})
        ) / 1000, 2)           AS distance_km
    FROM read_parquet('{PARQUET_FILES["cross_season"]}')
    WHERE ST_Distance_Sphere(
        ST_Point(centroid_lon, centroid_lat),
        ST_Point({DILLON_RESERVOIR[1]}, {DILLON_RESERVOIR[0]})
    ) <= {RADIUS_M}
    ORDER BY Area_m2 DESC
""").df()

print(f"  Cross-season hotspots near Dillon: {len(result)}")
if len(result) > 0:
    print(result.to_string(index=False))

#CLEANUP 
con.close()

print("\n" + "=" * 60)
print("All queries complete.")
print("=" * 60)