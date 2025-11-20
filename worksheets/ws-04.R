############################################################
# 02_download_and_prepare_data.R
# 
# Homework Script – Module 2: Geodata Preprocessing
# 
# Tasks:
# 1) Download DGM1 (1 m) for Hesse for the Burgwald region
# 2) Download CORINE Land Cover (CLC5 2018) for Germany
# 3) Download OSM data for Hesse (Geofabrik)
# 4) Clip all datasets to the Burgwald area of interest (AOI)
#
# IMPORTANT:
# - This script assumes you are inside an RStudio Project.
# - Paths are handled via here::here(), no setwd() is used.
# - Some URLs must be adapted by each group (see TODO markers).
############################################################

# -----------------------------
# 0) Packages & project root
# -----------------------------
library(here)
library(fs)
library(terra)
library(sf)

# Ensure that the directories exist (idempotent)
fs::dir_create(c(
  here::here("data/raw"),
  here::here("data/processed"),
  here::here("outputs/figures"),
  here::here("metadata",
  here::here("docs")           )
))

message("Project root: ", here::here())

# -----------------------------
# 1) Define Burgwald AOI
# -----------------------------
# Approximate bounding box for Burgwald (WGS84)
# Adapt if you want to be more precise.
burgwald_bbox <- c(
  xmin = 8.7,
  xmax = 9.0,
  ymin = 50.85,
  ymax = 51.05
)

aoi_burgwald_wgs <- st_as_sfc(
  st_bbox(burgwald_bbox, crs = 4326)
)

# Save AOI as GeoPackage for later use / inspection
aoi_file <- here::here("data", "processed", "aoi_burgwald.gpkg")
st_write(st_sf(geometry = aoi_burgwald_wgs), aoi_file, delete_dsn = TRUE)

source("worksheets/ws-04-helper_func.R")

# -------------------------------------------------------
# 3) DGM1 (1 m) for Hesse (Burgwald subset)
# -------------------------------------------------------
# NOTE: Each group must obtain a *direct* download URL for a DGM1 tile /
#       subset that fully covers the Burgwald AOI.
# 
# Sources (examples to search manually, DO NOT hard-code portals):
# - hoehendaten.de (DGM1, Open Data, 1 m)
# - HVBG Geodaten online (DGM1 Hesse, Open Data)
#
# TODO: Replace the placeholder URL below with your real DGM1 file URL.
#       Prefer GeoTIFF (.tif) or ZIP containing a GeoTIFF.
# ------------------------------------
# DGM1 Burgwald: mehrere Zips mosaiken
# ------------------------------------

options(timeout = 600)

# ---- 1) tagesaktuellen Ordner erzeugen ----
today <- format(Sys.Date(), "%Y%m%d")

base_url_wf <- paste0(
  "https://gds.hessen.de/downloadcenter/",
  today,
  "/3D-Daten/Digitales%20Gel%C3%A4ndemodell%20(DGM1)/Landkreis%20Waldeck-Frankenberg/"
)

base_url_mr <- paste0(
  "https://gds.hessen.de/downloadcenter/",
  today,
  "/3D-Daten/Digitales%20Gel%C3%A4ndemodell%20(DGM1)/Landkreis%20Marburg-Biedenkopf/"
)

# ---- 2) URLs für alle Teilgebiete (nur hier definieren!) ----
dgm1_urls <- c(
  burgwald    = paste0(base_url_wf, "Burgwald%20-%20DGM1.zip"),
  gemuenden   = paste0(base_url_wf, "Gem%C3%BCnden%20(Wohra)%20-%20DGM1.zip"),
  rosenthal   = paste0(base_url_wf, "Rosenthal%20-%20DGM1.zip"),
  muenchhausen = paste0(base_url_mr, "M%C3%BCnchhausen%20-%20DGM1.zip"),
  rauschenberg = paste0(base_url_mr, "Rauschenberg%20-%20DGM1.zip"),
  coelbe       = paste0(base_url_mr, "C%C3%B6lbe%20-%20DGM1.zip"),
  lahntal      = paste0(base_url_mr,"Lahntal%20-%20DGM1.zip"),
  wohra        = paste0(base_url_mr,"Wohratal%20-%20DGM1.zip"),
  wetter      = paste0(base_url_mr,"Wetter%20(Hessen)%20-%20DGM1.zip"),
  haina      = paste0(base_url_wf,"Haina%20(Kloster)%20-%20DGM1.zip")
)
# ---- 3) Zips herunterladen und entpacken ----
all_tif_files <- c()

for (nm in names(dgm1_urls)) {
  url_i <- dgm1_urls[[nm]]
  
  zip_file_i     <- here::here("data", "raw", paste0("dgm1_", nm, ".zip"))
  unzipped_dir_i <- here::here("data", "raw", paste0("dgm1_", nm))
  
  download_if_missing(url_i, zip_file_i, mode = "wb")
  
  if (!fs::dir_exists(unzipped_dir_i)) {
    fs::dir_create(unzipped_dir_i)
    unzip(zip_file_i, exdir = unzipped_dir_i)
  }
  
  tif_i <- dir(
    unzipped_dir_i,
    pattern = "\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  all_tif_files <- c(all_tif_files, tif_i)
}

if (length(all_tif_files) == 0) {
  stop("Keine TIF-Dateien in den DGM1-Zips gefunden.")
}

message("Gefundene DGM1-TIFs:\n", paste(all_tif_files, collapse = "\n"))

# ---- 4) Mosaic (terra) ----
dem_list <- lapply(all_tif_files, function(f) {
  r <- terra::rast(f)
  terra::crs(r) <- "EPSG:25832"   # HVBG Standard: ETRS89 / UTM 32N
  r
})

dem_merged <- do.call(terra::mosaic, c(dem_list, fun = "min"))

# ---- 5) AOI clip ----
# AOI (WGS84) -> CRS des DEM, dann als SpatVector
aoi_dem_sf <- sf::st_transform(aoi_burgwald_wgs, terra::crs(dem_merged))
aoi_dem_v  <- terra::vect(aoi_dem_sf)

dem_burgwald <- dem_merged |>
  terra::crop(aoi_dem_v) |>
  terra::mask(aoi_dem_v)

plot(dem_burgwald)
mb =aggregate(dem_burgwald,100)
mapview::mapview(mb)
dem_out_file <- here::here("data", "processed", "dem_dgm1_burgwald.tif")
terra::writeRaster(dem_burgwald, dem_out_file, overwrite = TRUE)

message("Gespeichertes mosaikiertes & geclippes DEM: ", dem_out_file)


# -------------------------------------------------------
# 1) Haupt-ZIP entpacken
clc_zip  <- here("data/raw/31916.zip")
clc_root <- here("data/raw/clc5_2018_copernicus")
unzip(clc_zip, exdir = clc_root)

# 2) Ordner unter Results/ finden (egal wie er heißt)
results_dir <- dir(clc_root, pattern = "Results", full.names = TRUE)[1]

# DATA-Ordner (wie bei dir im Filesystem vorhanden)
data_dir <- here("data", "raw", "u2018_clc2018_v2020_20u1_raster100m", "DATA")

# 3) inneres ZIP im Results-Ordner finden
zipname <- dir(results_dir, pattern = "\\.zip$", full.names = TRUE)

# 4) inneres ZIP nach data/raw entpacken
unzip(zipname, exdir = here("data", "raw"))

# 5) TIF darin finden
clc_tif <- dir(data_dir, pattern = "\\.tif$", full.names = TRUE)[1]

# 6) Raster laden
clc_rast <- rast(clc_tif)

# 7) AOI croppen
aoi_clc   <- st_transform(aoi_burgwald_wgs, crs(clc_rast))
aoi_clc_v <- vect(aoi_clc)
clc_crop  <- crop(clc_rast, aoi_clc_v) |> mask(aoi_clc_v)
plot(clc_crop)

# 8) speichern
writeRaster(
  clc_crop,
  here("data/processed/clc5_2018_burgwald_from_zip.tif"),
  overwrite = TRUE
)

# -------------------------------------------------------
# 5) OSM Data – Hesse (Geofabrik)
# -------------------------------------------------------
# We use Geofabrik "Hessen" shapefile bundle (latest snapshot).
# This contains multiple layers (roads, landuse, water, etc.).
#
# Source:
#   https://download.geofabrik.de/europe/germany/hessen.html
#
# Stable "latest" URL for shapefiles:
url_osm_hessen_shp <- "https://download.geofabrik.de/europe/germany/hessen-latest-free.shp.zip"

osm_zip_file <- here::here("data", "raw", "osm_hessen_latest_free.shp.zip")
osm_dir      <- here::here("data", "raw", "osm_hessen")

download_if_missing(url_osm_hessen_shp, osm_zip_file, mode = "wb")

# Unzip if necessary
if (!dir_exists(osm_dir)) {
  dir_create(osm_dir)
  unzip(osm_zip_file, exdir = osm_dir)
}

# We will use:
# - Roads:      gis_osm_roads_free_1.shp
# - Land use:   gis_osm_landuse_a_free_1.shp
#
# NOTE: The exact filenames may vary slightly; we search for them.
roads_file <- dir(osm_dir, pattern = "gis_osm_roads_free_1\\.shp$", recursive = TRUE, full.names = TRUE)
landuse_file <- dir(osm_dir, pattern = "gis_osm_landuse_a_free_1\\.shp$", recursive = TRUE, full.names = TRUE)

if (length(roads_file) == 0) {
  stop("OSM roads shapefile 'gis_osm_roads_free_1.shp' not found in: ", osm_dir)
}
if (length(landuse_file) == 0) {
  stop("OSM landuse shapefile 'gis_osm_landuse_a_free_1.shp' not found in: ", osm_dir)
}

roads_file    <- roads_file[1]
landuse_file  <- landuse_file[1]

osm_roads   <- st_read(roads_file, quiet = TRUE)
osm_landuse <- st_read(landuse_file, quiet = TRUE)

# Clip to Burgwald AOI (WGS84)
osm_roads_burgwald <- osm_roads |>
  st_transform(4326) |>
  st_intersection(aoi_burgwald_wgs)

osm_landuse_burgwald <- osm_landuse |>
  st_transform(4326) |>
  st_intersection(aoi_burgwald_wgs)

# Save clipped OSM layers as GeoPackage
osm_roads_out_file <- here::here("data", "processed", "osm_roads_burgwald.gpkg")
osm_landuse_out_file <- here::here("data", "processed", "osm_landuse_burgwald.gpkg")

st_write(osm_roads_burgwald, osm_roads_out_file, delete_dsn = TRUE)
st_write(osm_landuse_burgwald, osm_landuse_out_file, delete_dsn = TRUE)

message("Saved clipped OSM roads to: ", osm_roads_out_file)
message("Saved clipped OSM land use to: ", osm_landuse_out_file)

# -------------------------------------------------------
# 6) Summary
# -------------------------------------------------------
# At this point, you should have in data/processed:
# - aoi_burgwald.gpkg                  (AOI polygon)
# - dem_dgm1_burgwald.tif              (clipped DEM 1 m)
# - clc5_2018_burgwald_demgrid.tif     (optional: rasterized CLC on DEM grid)
# - osm_roads_burgwald.gpkg            (clipped OSM roads)
# - osm_landuse_burgwald.gpkg          (clipped OSM landuse)
#
# These will be the inputs for further analysis in the next module.
