library(terra)
library(sf)
library(here)

# -------------------------------------------------------
# 0) Ausgangsdaten
# -------------------------------------------------------
dem_1m   <- dem_burgwald                        # DGM Burgwald (Originalauflösung)
clc_100m <- clc_crop                            # CLC (100 m, aus Copernicus)

# optional: 100-m-DEM (mittlerer Wert pro 100-m-Zelle)
fact_100 <- round(100 / res(dem_1m)[1])         # Faktor aus Zellgröße ableiten
dem_100m <- aggregate(dem_1m, fact = fact_100, fun = mean)

# -------------------------------------------------------
# 1) Watersheds (Basins) aus DEM
#    (einmal Original, einmal 100 m)
# -------------------------------------------------------

# Originalauflösung
fd_1m      <- terrain(dem_1m, v = "flowdir")    # Flow-Direction
basins_1m  <- basin(fd_1m)                      # IDs pro Einzugsgebiet
basins_1m_poly <- as.polygons(basins_1m, dissolve = TRUE)
basins_1m_sf   <- st_as_sf(basins_1m_poly)

st_write(
  basins_1m_sf,
  here("data/processed/basins_1m_burgwald.gpkg"),
  delete_dsn = TRUE
)

# 100-m-Version
fd_100m     <- terrain(dem_100m, v = "flowdir")
basins_100m <- basin(fd_100m)
basins_100m_poly <- as.polygons(basins_100m, dissolve = TRUE)
basins_100m_sf   <- st_as_sf(basins_100m_poly)

st_write(
  basins_100m_sf,
  here("data/processed/basins_100m_burgwald.gpkg"),
  delete_dsn = TRUE
)

# -------------------------------------------------------
# 2) Standard-DEM-Ableitungen
#    (slope, aspect, TPI, TRI, roughness, curvature)
# -------------------------------------------------------

derivs_1m <- terrain(
  dem_1m,
  v    = c("slope", "aspect", "TPI", "TRI", "roughness", "curvature"),
  unit = "radians"
)

writeRaster(
  derivs_1m,
  here("data/processed/dem_derivatives_1m_burgwald.tif"),
  overwrite = TRUE
)

derivs_100m <- terrain(
  dem_100m,
  v    = c("slope", "aspect", "TPI", "TRI", "roughness", "curvature"),
  unit = "radians"
)

writeRaster(
  derivs_100m,
  here("data/processed/dem_derivatives_100m_burgwald.tif"),
  overwrite = TRUE
)

# -------------------------------------------------------
# 3) Landformen (ridge / valley / flat / slope) – single scale
#    einfache TPI-basierte Klassifikation
# -------------------------------------------------------

classify_landforms <- function(dem, win_cells = 11,
                               tpi_thresh = 1,  # in Höhen-Einheiten (z.B. Meter)
                               slope_thresh_deg = 5) {
  # TPI: Zelle minus Mittelwert in Nachbarschaft
  w      <- matrix(1, nrow = win_cells, ncol = win_cells)
  mean_n <- focal(dem, w = w, fun = mean, na.policy = "omit", na.rm = TRUE)
  tpi    <- dem - mean_n
  
  slope  <- terrain(dem, v = "slope", unit = "radians")
  slope_deg <- slope * 180 / pi
  
  # Klassenraster initialisieren
  lf <- dem
  values(lf) <- NA_integer_
  
  # 1 = Valley, 2 = Ridge, 3 = Flat, 4 = Slope
  vals <- values(tpi, mat = FALSE)
  slp  <- values(slope_deg, mat = FALSE)
  
  lf_vals <- rep(NA_integer_, length(vals))
  lf_vals[ vals <= -tpi_thresh & slp >= slope_thresh_deg ] <- 1  # Valley
  lf_vals[ vals >=  tpi_thresh & slp >= slope_thresh_deg ] <- 2  # Ridge
  lf_vals[ abs(vals) < tpi_thresh & slp < slope_thresh_deg ] <- 3  # Flat
  lf_vals[ is.na(lf_vals) ] <- 4  # Rest: Slope / übrige Formen
  
  values(lf) <- lf_vals
  lf
}

# Single-scale Landformen auf Originalauflösung
landforms_1m <- classify_landforms(dem_1m, win_cells = 11,
                                   tpi_thresh = 1, slope_thresh_deg = 5)

writeRaster(
  landforms_1m,
  here("data/processed/landforms_1m_single.tif"),
  overwrite = TRUE
)

# Single-scale Landformen auf 100 m
landforms_100m <- classify_landforms(dem_100m, win_cells = 5,
                                     tpi_thresh = 1, slope_thresh_deg = 5)

writeRaster(
  landforms_100m,
  here("data/processed/landforms_100m_single.tif"),
  overwrite = TRUE
)

# -------------------------------------------------------
# 4) Landformen multiskalig (z.B. "lokal" vs. "regional")
# -------------------------------------------------------

# Originalauflösung: kleine und große Fenster
landforms_1m_small <- classify_landforms(dem_1m, win_cells = 7,
                                         tpi_thresh = 0.5, slope_thresh_deg = 3)

landforms_1m_large <- classify_landforms(dem_1m, win_cells = 25,
                                         tpi_thresh = 2, slope_thresh_deg = 5)

landforms_1m_stack <- c(landforms_1m_small, landforms_1m_large)
names(landforms_1m_stack) <- c("landform_small", "landform_large")

writeRaster(
  landforms_1m_stack,
  here("data/processed/landforms_1m_multiscale.tif"),
  overwrite = TRUE
)

# 100-m-DEM: Multi-Skala (z.B. lokales vs. überlokales Relief)
landforms_100m_small <- classify_landforms(dem_100m, win_cells = 3,
                                           tpi_thresh = 0.5, slope_thresh_deg = 3)

landforms_100m_large <- classify_landforms(dem_100m, win_cells = 9,
                                           tpi_thresh = 1.5, slope_thresh_deg = 5)

landforms_100m_stack <- c(landforms_100m_small, landforms_100m_large)
names(landforms_100m_stack) <- c("landform_small", "landform_large")

writeRaster(
  landforms_100m_stack,
  here("data/processed/landforms_100m_multiscale.tif"),
  overwrite = TRUE
)
