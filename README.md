# Colorado Front Range SAR Avalanche Inventory

This is an end-to-end satellite-derived avalanche detection pipeline for Colorado's Front Range, built entirely on free, cloud-native tools. Three winters of Sentinel-1 C-band SAR data processed through Google Earth Engine, served as an interactive 3D map via PMTiles and MapLibre GL JS.

**[CLICK HERE TO VIEW WEBMAP](https://alecbiles-code.github.io/avalanche-inventory/)**

---

## Overview

Avalanche debris leaves a measurable signature in SAR backscatter. Disturbed snow and exposed terrain scatter C-band energy differently than undisturbed snowpack. This project detects those anomalies by comparing winter backscatter against a per-pixel summer baseline across three seasons (2021–22, 2022–23, 2023–24) and vectorizing statistically significant changes into an avalanche inventory.

Cross-season hotspots highlight terrain that triggered detections in two or more winters. The most operationally useful layer for backcountry travel planning.

---

## General Pipeline

```
Sentinel-1 GRD (GEE)
        │
        v
Speckle filter + Gamma-naught terrain normalization
        │
        v
Per-pixel summer baseline (median + MAD)
        │
        v
Multi-channel z-score anomaly (VV × 0.50 + VH × 0.35 + ratio × 0.15)
        │
        v
Terrain + water masking --> candidate polygons
        │
        v
Confidence classification (z > 6.0 = higher, 3.0–6.0 = lower)
        │
        v
GeoJSON --> GeoParquet --> DuckDB --> PMTiles --> GitHub Pages
```

---

## Detection Methodology

| Parameter | Value |
|---|---|
| Sensor | Sentinel-1 IW GRD, descending orbit |
| Polarizations | VV + VH |
| Baseline period | July–September (per season) |
| Anomaly metric | Robust z-score (MAD-normalized) |
| Detection threshold | Combined z-score ≥ 3.0 |
| Slope range | 20–70° |
| Aspect filter | None (all aspects retained) |
| Water mask | JRC Global Surface Water (occurrence < 50%) |
| Forest mask | Not applied (known limitation) |
| Resolution | 20 m |
| Higher confidence | Peak z-score > 6.0 |
| Lower confidence | Peak z-score 3.0–6.0 |

---

## Tech Stack

| Stage | Tool |
|---|---|
| SAR processing | Google Earth Engine (Python API) |
| Vectorization | GEE `connectedComponents` + `reduceToVectors` |
| Spatial analysis | GeoPandas, DuckDB spatial |
| Tile generation | Tippecanoe --> PMTiles |
| Visualization | MapLibre GL JS + PMTiles protocol |
| Hosting | GitHub Pages |

---

## Scripts

| File | Description |
|---|---|
| `scripts/gee_processing.py` | Full GEE pipeline  speckle filter, terrain normalization, z-score anomaly scoring, multi-season export |
| `scripts/vectorize.py` | Server-side vectorization and confidence classification |
| `scripts/parq_to_geojson.py` | Shapefile --> GeoParquet export with standardized schema |
| `scripts/duck_pull.py` | DuckDB spatial queries  50km radius filter, confidence breakdowns, area stats |

---

## Known Limitations

- **No forest mask applied**  detections in forested areas may include false positives from canopy snow loading/unloading
- **C-band amplitude only**  coherence-based detection (requires SLC data) would improve accuracy for small or shallow events
- **12-day revisit**  avalanches that settle and compact between Sentinel-1 passes problay not be detected
- **South-facing aspects retained**  increased false positive risk on sun-exposed slopes due to melt/refreeze cycles

---

## Validation

Detection results are intended for comparison against the [Colorado Avalanche Information Center (CAIC)](https://avalanche.state.co.us) observation database. Published SAR-based avalanche detection rates in comparable terrain range from 40–70% (Eckerstorfer et al., 2016; Vickers et al., 2016).

---

## References

1. Eckerstorfer, M. et al. (2016). Remote sensing of snow avalanches: Recent advances, potential, and limitations. *Cold Regions Science and Technology*, 121, 126–140.
2. Vickers, H. et al. (2016). Changes in SAR backscatter as an indicator of avalanche activity. *Cold Regions Science and Technology*, 135, 1–15.
3. Leinss, S. et al. (2020). Coherence analysis of SAR data for snow avalanche detection. *Remote Sensing of Environment*, 240, 111695.
4. Gorelick, N. et al. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18–27.
5. Colorado Avalanche Information Center (CAIC). Avalanche accident and observation database. colorado.gov/caic. Accessed 2024.

---

## Author

**Alec Biles**  
[LinkedIn](linkedin.com/in/alec-biles) · [GitHub](https://github.com/alecbiles-code)
