#!/usr/bin/env python3
"""
SAR Avalanche Inventory v2 Colorado Front Range

"""

#
# 1. IMPORTS AND INITIALIZZATION
# 
import ee
import sys

ee.Initialize(project='alec-avalanche')
print(" Initialized (project: alec-avalanche)")

# 
# 2. PARAM section
# 
# (coodrs for central spot. righ now dillion resovoir)
PARAMS = {
    'bbox': [-106.43, 39.41, -105.34, 39.95],

    'seasons': [
        ('21_22', '2021-07-01', '2021-09-30', '2022-01-01', '2022-04-30'),
        ('22_23', '2022-07-01', '2022-09-30', '2023-01-01', '2023-04-30'),
        ('23_24', '2023-07-01', '2023-09-30', '2024-01-01', '2024-04-30'),
    ],

    'orbit_direction': 'DESCENDING',
    'instrument_mode': 'IW',
    'enl': 4.4,                     # Equivalent Number of looks for S1 GRD IW
    'speckle_kernel_radius': 1,     # pixels (1 is 3×3 window)

    # Radar geometry for LIA computation
    'radar_look_azimuth_deg': 283.0,
    'ellipsoid_incidence_deg': 38.0,

    # Channel weights for combined z-score
    'w_vv': 0.45,
    'w_vh': 0.35,
    'w_ratio': 0.20,

    # Export settings
    'export_folder': 'avalanche_project',
    'export_crs': 'EPSG:32613',
    'export_scale_m': 20,
}

aoi = ee.Geometry.BBox(*PARAMS['bbox'])
export_tasks = []

print(f" Parameters loaded   {len(PARAMS['seasons'])} seasons")


# 
# 3. STATIC LAYERS
# 
dem = ee.Image('USGS/3DEP/10m').clip(aoi)
slope_deg = ee.Terrain.slope(dem).rename('slope')
aspect_deg = ee.Terrain.aspect(dem).rename('aspect')

###### Local incidence angle from DEM 
# All math uses single-band images renamed to a common name to avoid
# 
_s = slope_deg.multiply(3.14159265 / 180).rename('x')     # slope in rad
_a = aspect_deg.multiply(3.14159265 / 180).rename('x')    # aspect in rad
_te = ee.Image.constant(PARAMS['ellipsoid_incidence_deg']).multiply(3.14159265 / 180).rename('x')
_pa = ee.Image.constant(PARAMS['radar_look_azimuth_deg']).multiply(3.14159265 / 180).rename('x')
##
#equation with specal characters
# cos(LIA) = cos(slope)*cos(θ_el) + sin(slope)*sin(θ_el)*cos(φ_radar - aspect)
# With all bands named 'x', every arithmetic op matches
cos_lia = (
    _s.cos().multiply(_te.cos())
    .add(
        _s.sin().multiply(_te.sin()).multiply(
            _pa.subtract(_a).cos()
        )
    )
).clamp(-1, 1).rename('x')  # clamp to valid acos domain

lia_deg = cos_lia.acos().multiply(180 / 3.14159265).rename('lia')

#Terrain correction factor (dB) 
# γ0 = σ0 − 10·log10(cos(LIA))
# apply this as: corrected = sigma0_dB - correction_dB
# Floor cos(LIA) at 0.05 to avoid log(0) in shadow zones.
correction_db = cos_lia.max(0.05).log10().multiply(10).rename('correction')

#  quick diagnostic on LIA 
lia_sample = lia_deg.reduceRegion(
    reducer=ee.Reducer.minMax(),
    geometry=aoi, scale=500, maxPixels=1e6, bestEffort=True
).getInfo()
print(f"  LIA range: {lia_sample}")

print(" Static layers built")


# 
# 4. PREPROCESSING FUNCTIONS
# 

def load_s1(date_start, date_end):
    """Load S1 GRD IW collection for the study area and date range."""
    col = (ee.ImageCollection('COPERNICUS/S1_GRD')
           .filterBounds(aoi)
           .filterDate(date_start, date_end)
           .filter(ee.Filter.eq('instrumentMode', PARAMS['instrument_mode']))
           .filter(ee.Filter.eq('orbitProperties_pass', PARAMS['orbit_direction']))
           .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
           .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
           .select(['VV', 'VH', 'angle']))
    return col


def lee_filter(image):
    """Lee adaptive speckle filter in the linear power domain.

    Operates ONLY on VV and VH.  Angle band is never touched.
    All intermediate images are renamed to matching band names
    to guarantee GEE arithmetic operates correctly.
    """
    enl = PARAMS['enl']
    r = PARAMS['speckle_kernel_radius']
    kernel = ee.Kernel.square(r, 'pixels')

    # Process each polarization indepentent to avoid any multi-band
    # name-matching issues
    filtered_bands = []
    for pol in ['VV', 'VH']:
        # dB --> linear power
        band_db = image.select(pol).rename('b')
        band_lin = ee.Image(10).rename('b').pow(band_db.divide(10))

        # Local stas
        local_mean = band_lin.reduceNeighborhood(
            reducer=ee.Reducer.mean(), kernel=kernel
        ).rename('b')
        local_var = band_lin.reduceNeighborhood(
            reducer=ee.Reducer.variance(), kernel=kernel
        ).rename('b')

        # Lee weight: k = max(0, 1 - Cu²/Ci²)
        cu_sq = 1.0 / enl   # scalarnoise variance coefficient
        ci_sq = local_var.divide(local_mean.pow(2)).max(1e-10)
        k = ee.Image(1).rename('b').subtract(
            ee.Image(cu_sq).rename('b').divide(ci_sq)
        ).clamp(0, 1)

        # Apply: filtered = mean + k * (pixel - mean)
        filt_lin = local_mean.add(k.multiply(band_lin.subtract(local_mean)))

        # Linear --> dB
        filt_db = filt_lin.max(1e-10).log10().multiply(10).rename(pol)
        filtered_bands.append(filt_db)

    # Reassemble: filtered VV + filtered VH + original angle
    return (filtered_bands[0]
            .addBands(filtered_bands[1])
            .addBands(image.select('angle')))


def terrain_flatten(image):
    """Apply gamma-nought terrain correction using DEM-derived LIA.

    γ0_dB = σ0_dB − 10·log10(cos(LIA))

    The correction_db image (single band 'correction') is subtracted
    from each polarization band independently.
    """
    vv = image.select('VV').subtract(correction_db.rename('VV')).rename('VV')
    vh = image.select('VH').subtract(correction_db.rename('VH')).rename('VH')
    return vv.addBands(vh).addBands(image.select('angle'))


def add_ratio(image):
    """Add VV/VH ratio band (dB difference)."""
    ratio = image.select('VV').subtract(image.select('VH')).rename('ratio')
    return image.addBands(ratio)


def preprocess(image):
    """Full preprocessing: speckle filter --> terrain flatten --> add ratio."""
    return add_ratio(terrain_flatten(lee_filter(image)))


print(" Preprocessing functions defined")


# 
# 5. BASELINE computTION
# 

def build_baseline(summer_start, summer_end, label):
    """Build per-pixel temporal median and MAD from summer imagery.

    The temporal median of ~7-9 summer images acts as a multi-temporal
    speckle filter AND a robust central tendency estimator.
    """
    col = load_s1(summer_start, summer_end).map(preprocess)

    # how many images
    count = col.size().getInfo()
    print(f"  {label} summer images: {count}")
    if count == 0:
        print(f"  ⚠ ZERO summer images for {label}! Check date range / filters.")
        sys.exit(1)

    # Process each band independently for clarity
    bands = ['VV', 'VH', 'ratio']
    median = col.select(bands).median()   # per-pixel temporal median
    # median bands: ['VV', 'VH', 'ratio']   median() preserves names

    # MAD = median(|x_i - median|)
    def abs_dev(img):
        return img.select(bands).subtract(median).abs()

    mad = col.map(abs_dev).median()
    # Floor MAD to avoid /0   0.1 dB is well below measurement noise
    mad = mad.max(0.1)

    #  sample baseline values
    sample = median.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi, scale=500, maxPixels=1e6, bestEffort=True
    ).getInfo()
    print(f"  {label} baseline mean values: {sample}")

    mad_sample = mad.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi, scale=500, maxPixels=1e6, bestEffort=True
    ).getInfo()
    print(f"  {label} MAD mean values: {mad_sample}")

    return {'median': median, 'mad': mad}


print(" Baseline builder defined")


# 
# 6. ANOMALY SCORING / CHANGE DETECTION
# 

def score_winter_image(image, baseline):
    """Score one winter image against the summer baseline.

    For each band: z = (winter - baseline_median) / MAD
    Positive z --> backscatter INCREASED relative to summer.
    """
    bands = ['VV', 'VH', 'ratio']
    diff = image.select(bands).subtract(baseline['median'])
    z = diff.divide(baseline['mad'])
    # Also keep raw VV dB change for diagnostics
    vv_db = diff.select('VV').rename('VV_db_change')
    return z.addBands(vv_db)


def compute_anomaly(winter_start, winter_end, baseline, label):
    """Compute per-pixel max anomaly across all winter images.

    Returns a multi-band image with:
      VV_z, VH_z, ratio_z, combined_z, VV_db_max
    """
    col = load_s1(winter_start, winter_end).map(preprocess)

    #  Diagnostic 
    count = col.size().getInfo()
    print(f"  {label} winter images: {count}")
    if count == 0:
        print(f"  ⚠ ZERO winter images for {label}!")
        sys.exit(1)

    z_col = col.map(lambda img: score_winter_image(img, baseline))

    # Per-pixel MAX z across winter stack (captures peak debris signal)
    max_vv_z = z_col.select('VV').max().rename('VV_z')
    max_vh_z = z_col.select('VH').max().rename('VH_z')
    max_ratio_z = z_col.select('ratio').max().rename('ratio_z')
    max_vv_db = z_col.select('VV_db_change').max().rename('VV_db_max')

    #  Combined z-score (weighted average of channels) 
    #     ######  rename all channels to 'z' before arithmetic so
    # ee.Image.add() matches bands correctly.
    w_vv = PARAMS['w_vv']
    w_vh = PARAMS['w_vh']
    w_r = PARAMS['w_ratio']
    w_total = w_vv + w_vh + w_r

    combined_z = (
        max_vv_z.rename('z').multiply(w_vv)
        .add(max_vh_z.rename('z').multiply(w_vh))
        .add(max_ratio_z.rename('z').multiply(w_r))
    ).divide(w_total).rename('combined_z')

    # test sample anomaly values
    sample = combined_z.reduceRegion(
        reducer=ee.Reducer.percentile([5, 50, 95, 99]),
        geometry=aoi, scale=500, maxPixels=5e7, bestEffort=True
    ).getInfo()
    print(f"  {label} combined_z percentiles: {sample}")

    vv_db_sample = max_vv_db.reduceRegion(
        reducer=ee.Reducer.percentile([5, 50, 95, 99]),
        geometry=aoi, scale=500, maxPixels=5e7, bestEffort=True
    ).getInfo()
    print(f"  {label} VV dB change percentiles: {vv_db_sample}")

    # Stack all channels into one multi-band image for export
    result = (combined_z
              .addBands(max_vv_z)
              .addBands(max_vh_z)
              .addBands(max_ratio_z)
              .addBands(max_vv_db))

    return result


print(" Anomaly scoring defined")


# 
# 7–8. NOT USING ANYMORE!!!!
# 
print(" Thresholding/vectorization SKIPPED   exporting raw anomalies")


# 
# 9. PER-SEASON PROCESSING + CROSS-SEASON
# 

all_combined_z = []  # accumulate for cross-season max

for (label, s_start, s_end, w_start, w_end) in PARAMS['seasons']:
    print(f"\n{'='*60}")
    print(f"  SEASON {label}")
    print(f"{'='*60}")

    #  Baseline 
    baseline = build_baseline(s_start, s_end, label)

    #  Anomaly 
    anomaly = compute_anomaly(w_start, w_end, baseline, label)

    # Extract combined_z for cross-season stacking
    season_cz = anomaly.select('combined_z')
    all_combined_z.append(season_cz)

    #  EXPORTS 

    # (a) Full anomaly stack (combined_z + VV_z + VH_z + ratio_z + VV_db_max)
    task_anom = ee.batch.Export.image.toDrive(
        image=anomaly.toFloat().clip(aoi),
        description=f'anomaly_{label}_v2',
        folder=PARAMS['export_folder'],
        region=aoi,
        scale=PARAMS['export_scale_m'],
        crs=PARAMS['export_crs'],
        maxPixels=1e10,
        fileFormat='GeoTIFF',
        formatOptions={'cloudOptimized': True},
    )
    task_anom.start()
    export_tasks.append(f'anomaly_{label}_v2.tif')

    # (b) Baseline median (VV + VH + ratio)
    task_base = ee.batch.Export.image.toDrive(
        image=baseline['median'].toFloat().clip(aoi),
        description=f'baseline_{label}_v2',
        folder=PARAMS['export_folder'],
        region=aoi,
        scale=PARAMS['export_scale_m'],
        crs=PARAMS['export_crs'],
        maxPixels=1e10,
        fileFormat='GeoTIFF',
        formatOptions={'cloudOptimized': True},
    )
    task_base.start()
    export_tasks.append(f'baseline_{label}_v2.tif')

    # (c) Baseline MADnoise floor)
    task_mad = ee.batch.Export.image.toDrive(
        image=baseline['mad'].toFloat().clip(aoi),
        description=f'mad_{label}_v2',
        folder=PARAMS['export_folder'],
        region=aoi,
        scale=PARAMS['export_scale_m'],
        crs=PARAMS['export_crs'],
        maxPixels=1e10,
        fileFormat='GeoTIFF',
        formatOptions={'cloudOptimized': True},
    )
    task_mad.start()
    export_tasks.append(f'mad_{label}_v2.tif')

    print(f"  --> {label}: 3 exports submitted")


#  CROSS-SEASON: max combined_z across all seasons 
print(f"\n{'='*60}")
print("  CROSS-SEASON")
print(f"{'='*60}")

# Stack and take the max this is a continuous "hotspot intensity" raster.
# Pixels with high values in multiple seasons are best candidates.
cross_max = ee.ImageCollection(all_combined_z).max().rename('max_combined_z')

#  compute mean across seasonshotspot "persistence" indicator
cross_mean = ee.ImageCollection(all_combined_z).mean().rename('mean_combined_z')

# Stack both into one export
cross_stack = cross_max.addBands(cross_mean)

task_cross = ee.batch.Export.image.toDrive(
    image=cross_stack.toFloat().clip(aoi),
    description='cross_season_hotspot_v2',
    folder=PARAMS['export_folder'],
    region=aoi,
    scale=PARAMS['export_scale_m'],
    crs=PARAMS['export_crs'],
    maxPixels=1e10,
    fileFormat='GeoTIFF',
    formatOptions={'cloudOptimized': True},
)
task_cross.start()
export_tasks.append('cross_season_hotspot_v2.tif')

#  Export analysis mask and terrain for reference 
terrain_stack = (
    dem.rename('elevation')
    .addBands(slope_deg)
    .addBands(aspect_deg)
    .addBands(lia_deg)
)

task_terrain = ee.batch.Export.image.toDrive(
    image=terrain_stack.toFloat().clip(aoi),
    description='terrain_reference_v2',
    folder=PARAMS['export_folder'],
    region=aoi,
    scale=PARAMS['export_scale_m'],
    crs=PARAMS['export_crs'],
    maxPixels=1e10,
    fileFormat='GeoTIFF',
    formatOptions={'cloudOptimized': True},
)
task_terrain.start()
export_tasks.append('terrain_reference_v2.tif')

print("  --> Cross-season + terrain exports submitted")


#
# 10. EXPORT sECTUION)
# 


# 
# 11. SUMMARY
# 
print(f"\n{'='*60}")
print("  EXPORT SUMMARY")
print(f"{'='*60}")
print(f"  Folder : Google Drive / {PARAMS['export_folder']}")
print(f"  CRS    : {PARAMS['export_crs']}")
print(f"  Scale  : {PARAMS['export_scale_m']} m")
print(f"  Tasks  : {len(export_tasks)}")
print()
for i, name in enumerate(export_tasks, 1):
    print(f"    {i:2d}. {name}")
print()
print("  Per-season files (×3):")
print("    anomaly_{season}_v2.tif   5-band: combined_z, VV_z, VH_z, ratio_z, VV_db_max")
print("    baseline_{season}_v2.tif   3-band: VV, VH, ratio (summer median)")
print("    mad_{season}_v2.tif   3-band: VV, VH, ratio (summer MAD)")
print()
print("  Cross-season files:")
print("    cross_season_hotspot_v2.tif   2-band: max_combined_z, mean_combined_z")
print("    terrain_reference_v2.tif   4-band: elevation, slope, aspect, lia")
print()
print("  In ArcGIS Pro:")
print("    - combined_z > ~2.5 are candidate detections")
print("    - VV_db_max > ~1.0 dB confirms physical backscatter increase")
print("    - Use terrain_reference to mask slope < 25° or > 60°")
print("    - Cross-season max_combined_z shows persistent hotspots")
print()
print("  Monitor: ee.batch.Task.list()")
print(f"{'='*60}")
print(" Done.\n")


# 
# 12. COHERENCE-BASED DETECTION (FUTURE CHANGES)
# 
########
########  T
#########
# InSAR coherence measures phase correlation between two SAR acquisitions.
# Avalanche debris destroys scattering geometry  sharp coherence drop.
# This is especially powerful for dry continental snowpack where
# backscatter change is subtle but coherence loss is dramatic.
#
# GEE does not expose Sentinel-1 SLC data, so coherence requires
# offline processing:
#
# Workflow:
#   1. Download SLC pairs from ASF DAAC (same relative orbit, 6-12 day baseline)
#   2. Process with SNAP: Read --> Apply-Orbit-File  Back-Geocoding
#      --> Coherence (5×20 window) --> TOPSAR-Deburst  Terrain-Correction (3DEP)
#   3. Compute: coherence_change = summer_mean_coh − winter_coh
#   4. Fuse with backscatter z: combined = 0.5 * backscatter_z + 0.5 * coherence_z
#
# References:
#   Eckerstorfer et al. (2019) Remote Sensing 11(23):2863
#   Leinss et al. (2020) NHESS 20:1783-1803