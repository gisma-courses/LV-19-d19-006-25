# =========================
# Packages
# =========================
pkgs <- c(
  "sf","terra","gstat","fields","interp","mgcv","ggplot2","dplyr",
  "rnaturalearth","geostatsp","ranger","scales","mapview","htmlwidgets",
  "leaflet","RColorBrewer","viridisLite"
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) suppressMessages(install.packages(to_install, dependencies = TRUE))

suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))


options(sf_use_s2 = FALSE)
set.seed(42)

# =========================
# Data & base layers
# =========================
# Switzerland boundary (to LV95 / EPSG:2056)
ch_wgs84 <- rnaturalearth::ne_countries(country = "Switzerland", scale = 10, returnclass = "sf") |>
  st_make_valid() |>
  st_transform(4326)
ch_sf <- st_transform(ch_wgs84, 2056)
ch_v  <- terra::vect(ch_sf)

# Example data
data("swissRain",     package = "geostatsp")
data("swissAltitude", package = "geostatsp")
rain_v <- terra::unwrap(swissRain)      # SpatVector points
alt    <- terra::unwrap(swissAltitude)  # SpatRaster DEM

# Harmonize
rain_sf <- st_as_sf(rain_v)
if (is.na(st_crs(rain_sf))) st_crs(rain_sf) <- 2056 else rain_sf <- st_transform(rain_sf, 2056)
if (!grepl("2056", terra::crs(alt))) alt <- terra::project(alt, "EPSG:2056")
alt <- crop(alt, ch_v) |> mask(ch_v)
names(alt) <- "alt"
stopifnot("rain" %in% names(rain_sf))
rain_v <- terra::vect(rain_sf)

# Terrain derivatives (predictors)
slope_rad  <- terra::terrain(alt, v = "slope",  unit = "radians")
aspect_rad <- terra::terrain(alt, v = "aspect", unit = "radians")
slope_deg  <- terra::terrain(alt, v = "slope",  unit = "degrees"); names(slope_deg) <- "slope"
eastness   <- sin(aspect_rad); names(eastness) <- "eastness"
northness  <- cos(aspect_rad); names(northness) <- "northness"
tri        <- terra::terrain(alt, v = "TRI");       names(tri) <- "tri"
tpi        <- terra::terrain(alt, v = "TPI");       names(tpi) <- "tpi"
rough      <- terra::terrain(alt, v = "roughness"); names(rough) <- "roughness"

pred_stack <- c(alt, slope_deg, eastness, northness, tri, tpi, rough)

safe_extract <- function(stack, sfo) {
  as.data.frame(terra::extract(stack, terra::vect(sfo), ID = FALSE))
}

# =========================
# Common 1 km grid (template)
# =========================
border  <- 20000
res_m   <- 1000
ext0    <- terra::ext(rain_v)
ext_pad <- terra::ext(ext0$xmin - border, ext0$xmax + border,
                      ext0$ymin - border, ext0$ymax + border)

r_template <- terra::rast(ext = ext_pad, resolution = res_m, crs = terra::crs(rain_v))
xy <- terra::xyFromCell(r_template, 1:terra::ncell(r_template))
grid_sf <- sf::st_as_sf(data.frame(x = xy[,1], y = xy[,2]), coords = c("x","y"), crs = terra::crs(rain_v))
cc <- sf::st_coordinates(grid_sf)
grid_sf <- cbind(grid_sf, x = cc[,1], y = cc[,2])

# Points as plain df
rc <- st_coordinates(rain_sf)
rain_df <- data.frame(rain = rain_sf$rain, x = rc[,1], y = rc[,2])

# =========================
# Helper: XYZ -> Raster
# =========================
xyz_to_rast <- function(df_xyz, crs = "EPSG:2056") {
  if (!all(c("x","y") %in% names(df_xyz))) {
    if (all(c("X","Y") %in% names(df_xyz))) df_xyz <- dplyr::rename(df_xyz, x = X, y = Y)
    else stop("xyz_to_rast(): need x,y")
  }
  if (!"z" %in% names(df_xyz)) {
    if ("var1.pred" %in% names(df_xyz)) df_xyz$z <- df_xyz$var1.pred
    else if ("Z" %in% names(df_xyz))    df_xyz$z <- df_xyz$Z
    else stop("xyz_to_rast(): need z")
  }
  out <- df_xyz[, c("x","y","z")]
  out[] <- lapply(out, as.numeric)
  r <- terra::rast(out, type = "xyz", crs = crs)
  terra::crop(r, ch_v) |> terra::mask(ch_v)
}

# =========================
# Interpolation & Prediction (ALL models)
# =========================

## 1) Nearest neighbour (IDW nmax=1)
idw_NN <- gstat::idw(rain ~ 1, locations = rain_sf, newdata = grid_sf, nmax = 1)
r_NN <- xyz_to_rast(as.data.frame(cbind(st_coordinates(idw_NN), z = idw_NN$var1.pred)))
names(r_NN) <- "Nearest Neighbour"

## 2) IDW (idp = 0.5)
idw_IDW <- gstat::idw(rain ~ 1, locations = rain_sf, newdata = grid_sf, idp = 0.5, nmax = 10, nmin = 3)
r_IDW <- xyz_to_rast(as.data.frame(cbind(st_coordinates(idw_IDW), z = idw_IDW$var1.pred)))
names(r_IDW) <- "Inverse Distance Weighted"

## 3) Ordinary Kriging
vg_emp  <- gstat::variogram(rain ~ 1, locations = rain_sf)
vg_fit  <- gstat::fit.variogram(vg_emp, model = gstat::vgm(psill = NA, model = "Sph", range = NA, nugget = NA))
g_ok    <- gstat::gstat(formula = rain ~ 1, data = rain_sf, model = vg_fit)
kr_pred <- predict(g_ok, newdata = grid_sf)
r_KRIG  <- xyz_to_rast(as.data.frame(cbind(st_coordinates(kr_pred), z = kr_pred$var1.pred)))
names(r_KRIG) <- "Kriging"

## 4) Thin Plate Spline
coords_xy <- function(sfo) {
  m <- sf::st_coordinates(sfo); m <- m[,1:2, drop=FALSE]; storage.mode(m) <- "double"; m
}
Xy  <- coords_xy(rain_sf)
Xy0 <- coords_xy(grid_sf)
tps_fit  <- fields::Tps(Xy, rain_df$rain)
tps_pred <- predict(tps_fit, Xy0)
r_TPS    <- xyz_to_rast(data.frame(x = Xy0[,1], y = Xy0[,2], z = as.numeric(tps_pred)))
names(r_TPS) <- "Thin Plate Spline Regression"

## 5) TIN (triangulation) -> raster
tin_grid <- interp::interp(
  x = rain_df$x, y = rain_df$y, z = rain_df$rain,
  xo = sort(unique(grid_sf$x)),
  yo = sort(unique(grid_sf$y)),
  duplicate = "mean", extrap = FALSE, linear = TRUE, output = "grid"
)
r_TIN <- xyz_to_rast(interp::interp2xyz(tin_grid, data.frame = TRUE))
names(r_TIN) <- "Triangular Irregular Surface"

## 6) GAM s(x,y)
gam_fit  <- mgcv::gam(rain ~ s(x, y), data = rain_df)
gam_pred <- predict(gam_fit, newdata = data.frame(x = grid_sf$x, y = grid_sf$y))
r_GAM    <- xyz_to_rast(data.frame(x = grid_sf$x, y = grid_sf$y, z = as.numeric(gam_pred)))
names(r_GAM) <- "Generalized Additive Model"

## 7) KED (external drift: DEM + derivatives)
rain_ed <- cbind(rain_sf, safe_extract(pred_stack, rain_sf))
grid_ed <- cbind(grid_sf, safe_extract(pred_stack, grid_sf))
req <- c("alt","slope","eastness","northness")
ri  <- complete.cases(st_drop_geometry(rain_ed)[, c("rain", req), drop=FALSE])
gi  <- complete.cases(st_drop_geometry(grid_ed)[, req, drop=FALSE])
vg_ed_emp <- gstat::variogram(rain ~ alt + slope + eastness + northness, data = rain_ed[ri,])
vg_ed_fit <- gstat::fit.variogram(vg_ed_emp, model = gstat::vgm(model = "Sph"))
ked_pred  <- gstat::krige(rain ~ alt + slope + eastness + northness,
                          locations = rain_ed[ri,], newdata = grid_ed[gi,], model = vg_ed_fit)
r_KED <- xyz_to_rast(as.data.frame(cbind(st_coordinates(ked_pred), z = ked_pred$var1.pred)))
names(r_KED) <- "Kriging with External Drift"

## 8) Random Forests
features <- c("alt","slope","eastness","northness","tri","tpi","roughness")

# RF w/o XY
train_no_xy <- na.omit(data.frame(rain = rain_df$rain, safe_extract(pred_stack, rain_sf)[, features, drop=FALSE]))
rf_no_xy    <- ranger::ranger(rain ~ ., data = train_no_xy, num.trees = 500)
new_no_xy   <- safe_extract(pred_stack, grid_sf)[, features, drop=FALSE]
gi          <- complete.cases(new_no_xy)
rf_pred     <- predict(rf_no_xy, data = new_no_xy[gi, ])$predictions
r_RF_DEM    <- xyz_to_rast(data.frame(x = grid_sf$x[gi], y = grid_sf$y[gi], z = as.numeric(rf_pred)))
names(r_RF_DEM) <- "Random Forest (DEM only)"

# RF with XY
train_xy <- na.omit(data.frame(
  rain = rain_df$rain, x = rain_df$x, y = rain_df$y,
  safe_extract(pred_stack, rain_sf)[, features, drop=FALSE]
))
rf_xy   <- ranger::ranger(rain ~ ., data = train_xy, num.trees = 500)
new_xy  <- data.frame(x = grid_sf$x, y = grid_sf$y, safe_extract(pred_stack, grid_sf)[, features, drop=FALSE])
gi_xy   <- complete.cases(new_xy)
rf_pred_xy <- predict(rf_xy, data = new_xy[gi_xy, ])$predictions
r_RF_XY <- xyz_to_rast(data.frame(x = new_xy$x[gi_xy], y = new_xy$y[gi_xy], z = as.numeric(rf_pred_xy)))
names(r_RF_XY) <- "Random Forest (DEM + XY)"

# =========================
# Align (just in case) & build stacks
# =========================
# Use the OK grid as reference
ref <- r_KRIG
align <- function(x) terra::resample(x, ref, method = "bilinear")

pred_list <- list(
  r_TIN, r_IDW, r_TPS, r_KRIG, r_GAM, r_KED, r_RF_DEM, r_RF_XY
)
pred_list <- lapply(pred_list, align)

pred_stack_all <- do.call(c, pred_list)  # SpatRaster stack
# Save all predictions
terra::writeRaster(pred_stack_all, "predictions_all_models.tif", overwrite = TRUE)

# =========================
# Differences (aligned)
# =========================
r_DIFF_KED_OK   <- align(r_KED)   - align(r_KRIG)
names(r_DIFF_KED_OK) <- "KED − OK [mm]"
r_DIFF_RF_XY_NO <- align(r_RF_XY) - align(r_RF_DEM)
names(r_DIFF_RF_XY_NO) <- "RF(XY) − RF(no-XY) [mm]"

# =========================
# Build tidy data.frames for ggplot
# =========================
rast_to_df <- function(r, label = names(r)) {
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "value"
  d$model <- label
  d
}

models_df <- dplyr::bind_rows(lapply(1:nlyr(pred_stack_all), function(i) {
  rr <- pred_stack_all[[i]]
  rast_to_df(rr, names(rr))
}))

diff_df <- dplyr::bind_rows(
  rast_to_df(r_DIFF_KED_OK,   names(r_DIFF_KED_OK)),
  rast_to_df(r_DIFF_RF_XY_NO, names(r_DIFF_RF_XY_NO))
)

# Common color limits for model maps (robust range)
lims <- quantile(models_df$value, probs = c(0.02, 0.98), na.rm = TRUE)
lims <- c(floor(lims[1]), ceiling(lims[2]))

# Symmetric limits for differences
dlim <- quantile(abs(diff_df$value), probs = 0.98, na.rm = TRUE)
dlim <- as.numeric(dlim)

# =========================
# PLOTTING (ggplot2)
# =========================

# =========================
# PLOTTING (ggplot2) — fixed
# =========================

# Common limits (as before)
lims <- quantile(models_df$value, probs = c(0.02, 0.98), na.rm = TRUE)
lims <- c(floor(lims[1]), ceiling(lims[2]))

dlim <- quantile(abs(diff_df$value), probs = 0.98, na.rm = TRUE)
dlim <- as.numeric(dlim)

# 1) All model maps (faceted, 2 columns)
p_models <- ggplot(models_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = ch_sf, inherit.aes = FALSE,
          fill = NA, color = "grey15", linewidth = 0.2) +
  geom_point(data = rain_df, inherit.aes = FALSE,
             aes(x = x, y = y, size = rain),
             color = "dodgerblue", alpha = 0.9, stroke = 0) +
  scale_size(range = c(0.6, 3.0), guide = "none") +   # Größe aus Niederschlag, ohne Legende
  scale_fill_viridis_c(option = "inferno", limits = lims, oob = scales::squish,
                       name = "Rain [mm]") +
  coord_sf(crs = sf::st_crs(2056), expand = FALSE) +
  facet_wrap(~model, ncol = 2) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right")


# 2) Difference maps (faceted, 2 columns)
p_diffs <- ggplot(diff_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = ch_sf, inherit.aes = FALSE,  # <-- key fix
          fill = NA, color = "grey15", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c53030",
                       midpoint = 0, limits = c(-dlim, dlim), oob = scales::squish,
                       name = "Δ [mm]") +
  coord_sf(crs = sf::st_crs(2056), expand = FALSE) +
  facet_wrap(~model, ncol = 2) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text  = element_text(face = "bold"),
    legend.position = "right"
  )

# Show
print(p_models)
print(p_diffs)

# Save PNGs (same as before; adjust size if you want bigger maps / smaller points)
ggsave("models_faceted_2col.png", p_models, width = 10, height = 14, dpi = 300)
ggsave("diffs_faceted_2col.png",  p_diffs,  width = 10, height =  8, dpi = 300)

# ---- ONE-LEGEND MAPVIEW (predictions share one legend; rasters hidden by default) ----

# --- project to WGS84 for leaflet ---

to_wgs84 <- function(r) terra::project(r, "EPSG:4326")

pred_list <- list(
  "TIN"         = r_TIN,
  "IDW"         = r_IDW,
  "TPS"         = r_TPS,
  "OK"          = r_KRIG,
  "GAM"         = r_GAM,
  "KED"         = r_KED,
  "RF (DEM)"    = r_RF_DEM,
  "RF (DEM+XY)" = r_RF_XY
)
pred_wgs <- lapply(pred_list, to_wgs84)

# diffs (you already computed r_DIFF_KED_OK, r_DIFF_RF_XY_NO)
diff_wgs <- list(
  "KED − OK [mm]"            = to_wgs84(r_DIFF_KED_OK),
  "RF(XY) − RF(no-XY) [mm]"  = to_wgs84(r_DIFF_RF_XY_NO)
)

# --- palettes ---
# High-contrast BLUE -> WHITE for predictions

bw_vec   <- colorRampPalette(c("#ffffff", "#d8e7f7", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"))
pred_vals   <- unlist(lapply(pred_list, function(r) as.numeric(terra::values(r))))
pred_range  <- range(pred_vals, na.rm = TRUE)
pal_pred    <- colorNumeric(bw_vec(256), domain = pred_range, na.color = NA)
pred_breaks <- pretty(pred_range, n = 12)
if (tail(pred_breaks,1) <= pred_range[2]) pred_breaks <- c(pred_breaks, tail(pred_breaks,1) + diff(tail(pred_breaks,2)))
pred_cols   <- bw_vec(length(pred_breaks) - 1)

# Diverging BLUE–WHITE–RED for differences (no legend; keep hidden by default too)
dlim <- max(
  abs(stats::quantile(terra::values(r_DIFF_KED_OK),   c(.02,.98), na.rm=TRUE)),
  abs(stats::quantile(terra::values(r_DIFF_RF_XY_NO), c(.02,.98), na.rm=TRUE))
)
diff_breaks <- seq(-dlim, dlim, length.out = 11)
diff_cols   <- colorRampPalette(c("#08306b","#4292c6","#deebf7","#ffffff","#fdd0d0","#ef3b2c","#99000d"))(length(diff_breaks)-1)

# vectors
rain_sf_wgs <- sf::st_transform(rain_sf, 4326)
ch_sf_wgs   <- sf::st_transform(ch_sf,   4326)

# --- build mapview layers (suppress per-layer legends) ---
mv_pred <- mapply(function(r, nm)
  mapview(stars::st_as_stars(r),
          layer.name  = nm,
          at          = pred_breaks,
          col.regions = pred_cols,
          legend      = FALSE,
          maxpixels   = 8e5),
  pred_wgs, names(pred_wgs), SIMPLIFY = FALSE)

mv_diff <- mapply(function(r, nm)
  mapview(stars::st_as_stars(r),
          layer.name  = nm,
          at          = diff_breaks,
          col.regions = diff_cols,
          legend      = FALSE,
          maxpixels   = 8e5),
  diff_wgs, names(diff_wgs), SIMPLIFY = FALSE)

mv_pts    <- mapview(rain_sf_wgs, layer.name = "Stations",
                     cex = "rain" , zcol="rain", col.regions =bw_vec , alpha = 0.9)
mv_border <- mapview(ch_sf_wgs, layer.name = "CH border",
                     color = "black", lwd = 1.5, alpha.regions =0)  # transparent fill

# combine to one object
mv_all <- Reduce(`+`, c(mv_pred, mv_diff, list(mv_pts, mv_border)))

# convert to leaflet, hide all rasters by default
m <- mv_all@map
groups_to_hide <- c(names(pred_wgs), names(diff_wgs))
for (g in groups_to_hide) m <- leaflet::hideGroup(m, g)

# ONE shared legend for predictions (fix: position must be one of "topright" etc.)
m <- leaflet::addLegend(
  m, position = "topright",
  pal     = pal_pred,
  values  = pred_range,
  title   = "Rain [mm] (all models)",
  opacity = 1
)

# show / save
m
htmlwidgets::saveWidget(m, "swiss_models_one_legend_bluewhite.html", selfcontained = TRUE)

