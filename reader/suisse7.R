library(sf)
library(sp)
library(gstat)
library(fields)
library(interp)
library(mgcv)
library(automap)
library(viridis)
library(spatstat)
library(sfheaders)
library(ggplot2)
library(terra)
library(raster)
library(tmap)
library(tmaptools)
library(dplyr)
library(openxlsx)

# Lade Vektordaten (z.B. Gebiet, Stationen)
area = sf::st_read("mc_session/plot.shp")
stations = sf::st_read("mc_session/stations_prelim.gpkg")
statdat = openxlsx::read.xlsx("~/Downloads/all_GW1000A-WIFIFC29(202308290000-202308291731).xlsx")

# Konvertiere Gebiet in Spatial für 'sp'-kompatible Funktionen
area_sp = as(area, "Spatial")
plotborder = sp::spTransform(area_sp, CRSobj = "EPSG:32633")

extent.plot <- raster::extent(plotborder)

# Erstelle Gitter (700m x 700m)
plot_grid <- expand.grid(
  x = seq(from = round(extent.plot@xmin), to = round(extent.plot@xmax), by = 700),
  y = seq(from = round(extent.plot@ymin), to = round(extent.plot@ymax), by = 700)
)

# Vektorisiere Gebiet für terra
plotborder_vect <- terra::vect(plotborder)

# Voronoi-Polygone (optional)
t_swiss_rain_voronoi <- terra::voronoi(plotborder_vect)
sf_swiss_rain_voronoi <- st_as_sf(t_swiss_rain_voronoi)

# Beispiel-Daten: Angenommen, 'rain_df' enthält Spalten x, y, rain (Werte)
rain_sf <- sf::st_as_sf(stations)  # Oder andere Punktdaten mit Variable 'rain'
rain_df <- sfheaders::sf_to_df(rain_sf, fill = TRUE) %>% 
  dplyr::select(rain = your_variable_column, x = geometry_X, y = geometry_Y) # passe hier Spaltennamen an

# Bounding Box erweitern (Buffer 30 km)
bbox <- extent(rain_df$x, rain_df$y)
bbox@ymin <- bbox@ymin - 30000
bbox@xmin <- bbox@xmin - 30000
bbox@ymax <- bbox@ymax + 30000
bbox@xmax <- bbox@xmax + 30000

# Raster-Gitter erzeugen (1 km Auflösung)
alt_grd_template_sf <- bbox %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_make_grid(cellsize = c(1000, 1000), what = "centers") %>%
  st_as_sf() %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  mutate(Z = 0)

alt_grd_template_raster <- raster::rasterFromXYZ(alt_grd_template_sf[, c("x", "y", "Z")], crs = crs)

# ANNAHME: stations_dfr enthält columns:
# x, y, rain, elevation, slope, aspect

stations_dfr <- rain_df %>%
  mutate(
    elevation = some_elevation_data,  # Hier tatsächliche Daten einsetzen
    slope = some_slope_data,
    aspect = some_aspect_data
  )

stations_sp <- as(stations_dfr, "Spatial")

# 1. Nearest Neighbor
fit_NN <- gstat::gstat(formula = rain ~ 1, data = stations_sp, nmax = 10, nmin = 3)

# 2. Inverse Distance Weighting
fit_IDW <- gstat::gstat(formula = rain ~ 1, data = stations_sp, nmax = 10, nmin = 3, set = list(idp = 0.5))

# 3. Thin Plate Spline
mtx <- cbind(stations_dfr$x, stations_dfr$y, stations_dfr$rain)
colnames(mtx) <- c("x", "y", "rain")
fit_TPS <- fields::Tps(x = as.matrix(mtx[, c("x", "y")]), Y = stations_dfr$rain, miles = FALSE)

# 4. Generalized Additive Model
fit_GAM <- mgcv::gam(rain ~ s(x, y), data = stations_dfr)

# 5. Triangular Irregular Surface
fit_TIN <- interp::interp(
  x = stations_dfr$x,
  y = stations_dfr$y,
  z = stations_dfr$rain,
  xo = alt_grd_template_sf$x,
  yo = alt_grd_template_sf$y,
  output = "points"
) %>% bind_cols()

# 6. Automatized Kriging (ordinary kriging)
fit_KRIG <- automap::autoKrige(rain ~ 1, input_data = stations_sp)

# 7. Kriging with External Drift (Elevation, Slope, Aspect)
fit_KED <- automap::autoKrige(rain ~ elevation + slope + aspect, input_data = stations_sp)

# Vorhersagen auf Gitter

interp_NN <- mask(raster::interpolate(alt_grd_template_raster, fit_NN), plotborder_vect)
interp_IDW <- mask(raster::interpolate(alt_grd_template_raster, fit_IDW), plotborder_vect)
interp_TPS <- mask(raster::interpolate(alt_grd_template_raster, fit_TPS), plotborder_vect)

interp_GAM <- alt_grd_template_sf %>%
  mutate(z = predict(fit_GAM, .)) %>%
  raster::rasterFromXYZ(crs = crs) %>%
  mask(plotborder_vect)

interp_TIN <- mask(raster::rasterFromXYZ(fit_TIN, crs = crs), plotborder_vect)

interp_KRIG <- mask(raster::rasterFromXYZ(as.data.frame(fit_KRIG$krige_output)[, c("x1", "x2", "var1.pred")], crs = crs), plotborder_vect)

interp_KED <- mask(raster::rasterFromXYZ(as.data.frame(fit_KED$krige_output)[, c("x1", "x2", "var1.pred")], crs = crs), plotborder_vect)

# Namen vergeben
names(interp_NN) <- "Nearest Neighbor"
names(interp_IDW) <- "Inverse Distance Weighted"
names(interp_KRIG) <- "Kriging"
names(interp_KED) <- "Kriging with External Drift"
names(interp_TPS) <- "Thin Plate Spline Regression"
names(interp_TIN) <- "Triangular Irregular Surface"
names(interp_GAM) <- "Generalized Additive Model"

# Funktion zum Plotten
plot_my_rasters <- function(raster_object) {
  tmap_mode("view")
  tm_shape(raster_object) +
    tm_raster(
      breaks = c(0, 10, 20, 30, 40, 50, 60),
      style = "fixed",
      title = "Precipitation [mm]",
      palette = "Blues"
    ) +
    tm_shape(st_as_sf(plotborder)) +
    tm_borders(col = "red") +
    tm_graticules() +
    tm_layout(
      main.title = names(raster_object),
      main.title.position = "center",
      main.title.size = 1.5,
      legend.outside = TRUE
    )
}

# Plots erzeugen und synchron anordnen
plotlist <- lapply(list(
  interp_NN,
  interp_IDW,
  interp_KRIG,
  interp_KED,
  interp_TPS,
  interp_TIN,
  interp_GAM
), plot_my_rasters)

m2 <- tmap_arrange(plotlist, ncol = 2, nrow = 4, sync = TRUE)
