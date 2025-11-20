# -------------------------------------------------------
# SAGA-GIS Relief-Derivate mit Rsagacmd (sagacmd)
# -------------------------------------------------------
library(here)
library(terra)
library(Rsagacmd)

# einmal global initialisieren
saga <- saga_gis(
  raster_backend = "terra",
  vector_backend = "sf",
  all_outputs    = TRUE
)

saga_derive_relief <- function(dem,
                               prefix,
                               outdir     = here::here("data", "processed"),
                               minslope   = 0.1,
                               tpi_scales = c(100, 500, 1000),
                               saga       = saga) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ## 1) Sinks füllen (Wang & Liu) – Rsagacmd gibt direkt das gefüllte DEM zurück
  dem_filled <- saga$ta_preprocessor$fill_sinks_xxl_wang_liu(
    elev     = dem,
    minslope = minslope
  )
  
  terra::writeRaster(
    dem_filled,
    file.path(outdir, paste0(prefix, "_filled_dem.tif")),
    overwrite = TRUE
  )
  
  ## 2) Standard-Ableitungen (Slope, Aspect, Curvatures)
  lsps <- saga$ta_morphometry$slope_aspect_curvature(
    elevation   = dem_filled,  # nicht dem = ...
    method      = 6,           # z.B. "maximum slope" (siehst du in available_opts)
    unit_slope  = 0,           # 0 = radians
    unit_aspect = 0            # 0 = radians
  )
  
  # lsps ist eine LISTE von SpatRaster-Objekten -> Stack bauen:
  lsps_stack <- terra::rast(lsps)
  
  terra::writeRaster(
    lsps_stack,
    here::here("data", "processed", paste0(prefix, "_slope_aspect_curvature.tif")),
    overwrite = TRUE
  )
  
  ## 3) Multi-scale TPI
  tpi_list <- list()
  for (R in tpi_scales) {
    tpi_r <- saga$ta_morphometry$topographic_position_index_tpi(
      dem        = dem_filled,
      radius_min = 0,
      radius_max = R,
      standard   = TRUE   # standardisierte TPI-Werte
    )
    out_tpi <- file.path(outdir, paste0(prefix, "_tpi_r", R, "m.tif"))
    terra::writeRaster(tpi_r, out_tpi, overwrite = TRUE)
    tpi_list[[paste0("tpi_r", R)]] <- tpi_r
  }
  
  invisible(list(
    dem_filled = dem_filled,
    lsps       = lsps,
    tpi        = tpi_list
  ))
}

# -------------------------------------------------------
# 1) Originalauflösung (DGM1 Burgwald, ~1 m)
# -------------------------------------------------------
dem_1m <- dem = here::here("data", "processed", "dem_dgm1_burgwald.tif")


# Originalauflösung
relief_1m <- saga_derive_relief(
  dem    = dem_burgwald,
  prefix = "dem1m_burgwald",
  outdir = here::here("data", "processed", "relief_1m")
)

# 100 m DEM in R aggregieren (Beispiel)
dem_100m <- aggregate(rast(dem_1m), 100)

relief_100m <- saga_derive_relief(
  dem    = dem_100m,
  prefix = "dem100m_burgwald",
  outdir = here::here("data", "processed", "relief_100m")
)
