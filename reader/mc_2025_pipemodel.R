
# Crisp figures
# EN: Crisp figures
knitr::opts_chunk$set(fig.width = 9, fig.height = 6, dpi = 150)
# Alle Chunk-Meldungen global weg
# EN: Silence messages/warnings in all chunks
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

# Hilfsfunktion: gstat & Co. ruhigstellen
# EN: Helper: silence gstat & friends
quiet <- function(expr) suppressWarnings(suppressMessages(force(expr)))


# Packages ---------------------------------------------------------------
req_pkgs <- c(
  "terra","sf","ggplot2","dplyr","tibble","tidyr",
  "suncalc","gstat","randomForest","mgcv","scales","patchwork",
  "knitr","kableExtra","RColorBrewer", "zoo"
)
inst <- rownames(installed.packages())
if (any(!req_pkgs %in% inst)) install.packages(setdiff(req_pkgs, inst), dependencies = TRUE)
invisible(lapply(req_pkgs, require, character.only = TRUE))
sf::sf_use_s2(FALSE)  # robust joins in small projected domains
set.seed(42)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Domain & grid ---------------------------------------------------------
crs_utm <- "EPSG:32632"
E0 <- 600000; N0 <- 5725000
len_x <- 16000; len_y <- 12000; res <- 10
ext <- terra::ext(E0 - len_x/2, E0 + len_x/2, N0 - len_y/2, N0 + len_y/2)
Rtemplate <- terra::rast(ext, resolution = res, crs = crs_utm)

xmin <- terra::xmin(ext); xmax <- terra::xmax(ext)
ymin <- terra::ymin(ext); ymax <- terra::ymax(ext)
x0 <- (xmin+xmax)/2;      y0 <- (ymin+ymax)/2

# Feature anchors (left/right thirds) ----------------------------------
x_hill_center <- xmin + len_x/3;   y_hill_center <- y0
x_lake_center <- xmin + 2*len_x/3; y_lake_center <- y0

# Scenario toggles ------------------------------------------------------
lake_choice <- "water"   # "none" | "water" | "hollow"
hill_choice <- "bump"    # "none" | "bump"
 
# Feature geometry
lake_diam_m  <- 1000; lake_depth_m <- 10; smooth_edges <- FALSE
hill_diam_m  <- 800; hill_height_m <- 90; hill_smooth  <- TRUE

# Night pooling reduction over hill (0..1)
pool_block_gain <- 0.4

no_random_hill = 300
micro_hill_diam_m = 40
micro_hill_height_m = 50
micro_hill_smooth = TRUE

# Stations --------------------------------------------------------------
station_mode      <- "random"     # "random" | "ns_transect" | "ew_transect"
n_st              <- 60
transect_margin_m <- 10
ns_offset_m <- 0   # + east / - west
ew_offset_m <- 0   # + north / - south

# CV & models -----------------------------------------------------------
# --- Global einheitlich ---
compute_block_size <- function(len_x, len_y, n_st,
                               target_st_per_block = 3,
                               min_blocks_axis = 3,
                               round_to = 50) {
  area <- len_x * len_y
  B_target <- max(min_blocks_axis^2, round(n_st / target_st_per_block))
  bs <- sqrt(area / B_target)                         # idealisierte Blockkante
  bs <- min(bs, len_x / min_blocks_axis, len_y / min_blocks_axis)  # mind. 3 pro Achse
  bs <- max(round_to, round(bs / round_to) * round_to)             # "schön" runden
  as.integer(bs)
}

# anwenden
block_size  <- compute_block_size(len_x, len_y, n_st)


.get_block_size <- function() {
  bs <- get0("block_size", ifnotfound = 100)
  if (!is.numeric(bs) || !is.finite(bs) || bs <= 0) 100 else as.numeric(bs)
}
models_use <- c("Voronoi","IDW","OK","KED","RF","GAM")

# Viz palettes ----------------------------------------------------------
temp_palette <- grDevices::colorRampPalette(c("#0000FF","#FF0000"))(256)  # blue->red
stretch_q    <- c(0.02, 0.98)



# Domain diagonal (used for variogram fallbacks)
dom_diag <- sqrt((xmax - xmin)^2 + (ymax - ymin)^2)

lat <- 51.8; lon <- 10.6
sun_pos_utc <- function(y, m, d, h, lat, lon) {
  t  <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", y, m, d, h), tz = "UTC")
  sp <- suncalc::getSunlightPosition(date = t, lat = lat, lon = lon)
  az_from_north <- (sp$azimuth + pi) %% (2*pi)
  list(alt = sp$altitude, az = az_from_north)
}
sun14 <- sun_pos_utc(2024, 6, 21, 14, lat, lon)
sun05 <- sun_pos_utc(2024, 6, 21,  5, lat, lon)

# Cosine of incidence (sun on slope/aspect)
cosi_fun <- function(alt, az, slp_r, asp_r) {
  zen <- (pi/2 - alt)
  ci  <- cos(slp_r)*cos(zen) + sin(slp_r)*sin(zen)*cos(az - asp_r)
  terra::ifel(ci < 0, 0, ci)
}

# 1) Topography ---------------------------------------------------------
# 1) Topography ---------------------------------------------------------
build_topography <- function(lake_mode = c("none","water","hollow"),
                             hill_mode = c("none","bump")) {
  lake_mode <- match.arg(lake_mode); hill_mode <- match.arg(hill_mode)
  XY <- as.data.frame(terra::xyFromCell(Rtemplate, 1:terra::ncell(Rtemplate))); names(XY) <- c("x","y")
  dy <- XY$y - y0
  a  <- 100 / ((len_y/2)^2)              # ~100 m rim height
  elev <- 500 + a * dy^2

  # Pond/hollow (right third)
  rl <- sqrt((XY$x - x_lake_center)^2 + (XY$y - y_lake_center)^2); lr <- lake_diam_m/2
  if (lake_mode %in% c("water","hollow")) {
    w_l <- if (smooth_edges) pmax(0, 1 - (rl/lr)^2) else as.numeric(rl <= lr)
    elev <- elev - lake_depth_m * w_l
  } else w_l <- 0

  # Hill (left third)
  if (hill_mode == "bump") {
    rh <- sqrt((XY$x - x_hill_center)^2 + (XY$y - y_hill_center)^2); hr <- max(1e-6, hill_diam_m/2)
    w_h <- if (hill_smooth) exp(- (rh/hr)^2) else as.numeric(rh <= hr)
    elev <- elev + hill_height_m * w_h
  } else w_h <- 0

  E <- Rtemplate; terra::values(E) <- elev; names(E) <- "elev"
  lakeR <- Rtemplate; terra::values(lakeR) <- if (lake_mode=="water") as.numeric(w_l>0) else 0; names(lakeR) <- "lake"
  hillW <- Rtemplate; terra::values(hillW) <- if (hill_mode=="bump") w_h else 0; names(hillW) <- "hillW"

  # Derivatives & sun
  slp  <- terra::terrain(E, v="slope",  unit="radians")
  asp  <- terra::terrain(E, v="aspect", unit="radians")
  slp0 <- terra::ifel(is.na(slp), 0, slp); asp0 <- terra::ifel(is.na(asp), 0, asp)
  I14  <- cosi_fun(sun14$alt, sun14$az, slp0, asp0)
  I05  <- cosi_fun(sun05$alt, sun05$az, slp0, asp0)

  list(E = E, lake = lakeR, hillW = hillW, slp = slp0, asp = asp0, I14 = I14, I05 = I05)
}


# 2) Land-cover ---------------------------------------------------------
# --- Land-Cover: zentrale Definition (1 forest, 2 water, 3 bare soil, 4 maize)
lc_levels_default <- c("forest","water","bare soil","maize")
lc_levels <- getOption("pipemodel.lc_levels", lc_levels_default)
lc_colors <- c("forest"="#2E8B57","water"="#5DADE2","bare soil"="#C49A6C","maize"="#F4D03F")
lc_colors_default <- c(
  "forest"          = "#2E8B57",
  "water"        = "#5DADE2",
  "bare soil"  = "#C49A6C",
  "maize"          = "#F4D03F"
)

# 2) Land-cover ---------------------------------------------------------
# --- Land-Cover: zentrale Definition (1 forest, 2 water, 3 bare soil, 4 maize)
lc_levels_default <- c("forest","water","bare soil","maize")
lc_levels <- getOption("pipemodel.lc_levels", lc_levels_default)
lc_colors <- c("forest"="#2E8B57","water"="#5DADE2","bare soil"="#C49A6C","maize"="#F4D03F")
lc_colors_default <- c(
  "forest"          = "#2E8B57",
  "water"        = "#5DADE2",
  "bare soil"  = "#C49A6C",
  "maize"          = "#F4D03F"
)
# 3) Physics: T14/T05 from topo + land-cover ----------------------------
# Day: solar sensitivity (alpha_I) and canopy shading factor for cos(i)
alpha_I_by_lc <- c("forest"=3.5, "water"=1.5, "bare soil"=6.0, "maize"=4.5)
shade_fac_by_lc <- c("forest"=0.6, "water"=1.0, "bare soil"=1.0, "maize"=0.9)
# Dawn: additive warm/cool biases (°C) + pooling modifiers (multiplier)
dawn_bias_by_lc <- c("forest"=0.3, "water"=1.2, "bare soil"=-0.5, "maize"=0.1)
pool_fac_by_lc  <- c("forest"=0.7, "water"=0.8, "bare soil"=1.1, "maize"=1.0)

# Koeffizienten (benannte Vektoren in Deutsch; nur falls nicht schon definiert)
if (!exists("alpha_I_by_lc")) {
  alpha_I_by_lc <- c("forest"=3.5,"water"=1.5,"bare soil"=6.0,"maize"=4.5)
}
if (!exists("shade_fac_by_lc")) {
  shade_fac_by_lc <- c("forest"=0.6,"water"=1.0,"bare soil"=1.0,"maize"=0.9)
}
if (!exists("dawn_bias_by_lc")) {
  dawn_bias_by_lc <- c("forest"=0.3,"water"=1.2,"bare soil"=-0.5,"maize"=0.1)
}
if (!exists("pool_fac_by_lc")) {
  pool_fac_by_lc  <- c("forest"=0.7,"water"=0.8,"bare soil"=1.1,"maize"=1.0)
}

# ------------------------------------------------------------------------------
# build_physics_fields(topography, landcover, noise14, noise05)
# Purpose:
#   Create the synthetic temperature "truth" fields for day (T14) and pre-dawn (T05)
#   using topography, land-cover-dependent coefficients, solar geometry, and small noise.
# Inputs:
#   - topography: list containing E, slp, I14, I05, hillW (from build_topography or scenario).
#   - landcover:  SpatRaster (or list with lc) of integer LC codes (1..4) per cell.
#   - noise14/noise05: SpatRasters with Gaussian noise to add small texture.
# Model (kept exactly as in code):
#   T14 = T0_14 + lapse_14*(E - mean(E)) + alpha_I(LC) * (I14 * shade_fac(LC)) + noise
#   T05 = T0_05 + inv_05 *(E - mean(E)) + eta_slope * slope
#                  - pool_base*(1 - pool_block_gain*hillW)*pool_fac(LC)
#                  + dawn_bias(LC) + noise
# Returns:
#   list(R14 = T14 raster, R05 = T05 raster), names preserved exactly.
# ------------------------------------------------------------------------------
build_physics_fields <- function(topography, landcover, noise14, noise05) {
  E    <- topography$E
  slp0 <- topography$slp
  I14  <- topography$I14
  I05  <- topography$I05
  hillW<- topography$hillW

  # accept either a list(list(lc=...)) or a SpatRaster directly
  lc <- if (inherits(landcover, "SpatRaster")) landcover else landcover$lc
  stopifnot(inherits(lc, "SpatRaster"))

  # lc (numerisch 1..4) -> Zeichenklassen gemäß lc_levels
  v  <- as.integer(terra::values(lc))
  v[!is.finite(v)] <- 1L
  v <- pmax(1L, pmin(v, length(lc_levels)))
  lc_char <- factor(lc_levels[v], levels = lc_levels)

  # Karten aus benannten Vektoren ableiten
  to_r <- function(x) terra::setValues(terra::rast(E), x)
  alpha_I <- to_r(as.numeric(alpha_I_by_lc[lc_char]))
  shade_f <- to_r(as.numeric(shade_fac_by_lc[lc_char]))
  dawn_b  <- to_r(as.numeric(dawn_bias_by_lc[lc_char]))
  pool_f  <- to_r(as.numeric(pool_fac_by_lc[lc_char]))

  # Effektive Strahlung (forest beschattet)
  I14_eff <- I14 * shade_f
  I05_eff <- I05 * shade_f  # der Vollständigkeit halber

  # Baselines
  E_mean <- terra::global(E, "mean", na.rm = TRUE)[1,1]
  Y <- terra::init(E, "y"); dist2ax <- abs(Y - (terra::ymax(E)+terra::ymin(E))/2); w_pool <- 70
  pool_base <- 4.0 * exp(- (dist2ax / w_pool)^2)
  pool_mod  <- pool_base * (1 - pool_block_gain * hillW) * pool_f

  # Tag (14 UTC)
  T0_14 <- 26.0; lapse_14 <- -0.0065
  R14 <- T0_14 + lapse_14 * (E - E_mean) + alpha_I * I14_eff + noise14; names(R14) <- "T14"

  # Nacht/Früh (05 UTC)
  T0_05 <- 8.5; inv_05 <- 0.003; eta_slope <- 0.6
  R05 <- T0_05 + inv_05 * (E - E_mean) + eta_slope * slp0 - pool_mod + dawn_b + noise05; names(R05) <- "T05"

  list(R14 = R14, R05 = R05)
}


# Noise rasters (generated once) --------------------------------------
set.seed(1001)
noise14_r <- terra::setValues(terra::rast(Rtemplate), rnorm(terra::ncell(Rtemplate), 0, 0.3))
set.seed(1002)
noise05_r <- terra::setValues(terra::rast(Rtemplate), rnorm(terra::ncell(Rtemplate), 0, 0.3))

# ------------------------------------------------------------------------------
# build_physics_fields(topography, landcover, noise14, noise05)
# Purpose:
#   Create the synthetic temperature "truth" fields for day (T14) and pre-dawn (T05)
#   using topography, land-cover-dependent coefficients, solar geometry, and small noise.
# Inputs:
#   - topography: list containing E, slp, I14, I05, hillW (from build_topography or scenario).
#   - landcover:  SpatRaster (or list with lc) of integer LC codes (1..4) per cell.
#   - noise14/noise05: SpatRasters with Gaussian noise to add small texture.
# Model (kept exactly as in code):
#   T14 = T0_14 + lapse_14*(E - mean(E)) + alpha_I(LC) * (I14 * shade_fac(LC)) + noise
#   T05 = T0_05 + inv_05 *(E - mean(E)) + eta_slope * slope
#                  - pool_base*(1 - pool_block_gain*hillW)*pool_fac(LC)
#                  + dawn_bias(LC) + noise
# Returns:
#   list(R14 = T14 raster, R05 = T05 raster), names preserved exactly.
# ------------------------------------------------------------------------------
build_physics_fields <- function(topography, landcover, noise14, noise05) {
  E    <- topography$E
  slp0 <- topography$slp
  I14  <- topography$I14
  I05  <- topography$I05
  hillW<- topography$hillW

  # accept either a list(list(lc=...)) or a SpatRaster directly
  lc <- if (inherits(landcover, "SpatRaster")) landcover else landcover$lc
  stopifnot(inherits(lc, "SpatRaster"))

  # lc (numerisch 1..4) -> Zeichenklassen gemäß lc_levels
  v  <- as.integer(terra::values(lc))
  v[!is.finite(v)] <- 1L
  v <- pmax(1L, pmin(v, length(lc_levels)))
  lc_char <- factor(lc_levels[v], levels = lc_levels)

  # Karten aus benannten Vektoren ableiten
  to_r <- function(x) terra::setValues(terra::rast(E), x)
  alpha_I <- to_r(as.numeric(alpha_I_by_lc[lc_char]))
  shade_f <- to_r(as.numeric(shade_fac_by_lc[lc_char]))
  dawn_b  <- to_r(as.numeric(dawn_bias_by_lc[lc_char]))
  pool_f  <- to_r(as.numeric(pool_fac_by_lc[lc_char]))

  # Effektive Strahlung (forest beschattet)
  I14_eff <- I14 * shade_f
  I05_eff <- I05 * shade_f  # der Vollständigkeit halber

  # Baselines
  E_mean <- terra::global(E, "mean", na.rm = TRUE)[1,1]
  Y <- terra::init(E, "y"); dist2ax <- abs(Y - (terra::ymax(E)+terra::ymin(E))/2); w_pool <- 70
  pool_base <- 4.0 * exp(- (dist2ax / w_pool)^2)
  pool_mod  <- pool_base * (1 - pool_block_gain * hillW) * pool_f

  # Tag (14 UTC)
  T0_14 <- 26.0; lapse_14 <- -0.0065
  R14 <- T0_14 + lapse_14 * (E - E_mean) + alpha_I * I14_eff + noise14; names(R14) <- "T14"

  # Nacht/Früh (05 UTC)
  T0_05 <- 8.5; inv_05 <- 0.003; eta_slope <- 0.6
  R05 <- T0_05 + inv_05 * (E - E_mean) + eta_slope * slp0 - pool_mod + dawn_b + noise05; names(R05) <- "T05"

  list(R14 = R14, R05 = R05)
}

# Noise rasters (generated once) --------------------------------------
set.seed(1001)
noise14_r <- terra::setValues(terra::rast(Rtemplate), rnorm(terra::ncell(Rtemplate), 0, 0.3))
set.seed(1002)
noise05_r <- terra::setValues(terra::rast(Rtemplate), rnorm(terra::ncell(Rtemplate), 0, 0.3))

# One-stop scenario (mit optionalen Mikro-Hügeln) -----------------------
# ------------------------------------------------------------------------------
# build_scenario(...)
# Purpose:
#   One-stop constructor for a complete synthetic scenario including:
#   - elevation (E), lake mask, hill weight
#   - slope/aspect and solar incidence rasters (I14, I05)
#   - land cover (lc) with fixed integer classes (do NOT rename class labels)
#   - physics-based temperatures (R14, R05)
#   This function can optionally add micro-hills for additional relief texture.
# Inputs:
#   - lake_mode, hill_mode: same semantics as build_topography.
#   - noise14, noise05: noise rasters injected into temperatures (unchanged here).
#   - random_hills, hills_xy, micro_*: optional micro-relief controls.
# Returns:
#   Named list with rasters and color/level metadata: (E, R14, R05, lake, hillW,
#   slp, asp, I14, I05, lc, lc_levels, lc_colors).
# Notes:
#   - Class labels like "forest", "water", "bare soil", "maize" remain unchanged
#     because other parts of the code rely on them as factor labels.
# ------------------------------------------------------------------------------
build_scenario <- function(lake_mode = c("none","water","hollow"),
                           hill_mode = c("none","bump"),
                           noise14 = noise14_r,
                           noise05 = noise05_r,
                           random_hills = NULL,
                           hills_xy = NULL,
                           micro_hill_diam_m =NULL,
                           micro_hill_height_m = NULL,
                           micro_hill_smooth = NULL,
                           micro_seed = NULL) {
  lake_mode <- match.arg(lake_mode); hill_mode <- match.arg(hill_mode)

  # Basis-Talform
  XY <- as.data.frame(terra::xyFromCell(Rtemplate, 1:terra::ncell(Rtemplate))); names(XY) <- c("x","y")
  dy <- XY$y - y0; a <- 100 / ((len_y/2)^2); elev <- 500 + a * dy^2

  # See/Grube
  rl <- sqrt((XY$x - x_lake_center)^2 + (XY$y - y_lake_center)^2); lr <- lake_diam_m/2
  if (lake_mode %in% c("water","hollow")) {
    w_l <- if (smooth_edges) pmax(0, 1 - (rl/lr)^2) else as.numeric(rl <= lr)
    elev <- elev - lake_depth_m * w_l
  } else w_l <- 0

  # Haupt-Hügel
  if (hill_mode == "bump") {
    rh <- sqrt((XY$x - x_hill_center)^2 + (XY$y - y_hill_center)^2); hr <- max(1e-6, hill_diam_m/2)
    w_h_main <- if (hill_smooth) exp(-(rh/hr)^2) else as.numeric(rh <= hr)
    elev <- elev + hill_height_m * w_h_main
  } else w_h_main <- 0

  # Mikro-Hügel (zufällig + manuell)
  centers <- NULL
  if (!is.null(hills_xy)) centers <- as.matrix(hills_xy[,1:2, drop = FALSE])
  if (random_hills > 0) {
    if (!is.null(micro_seed)) set.seed(micro_seed)
    margin <- micro_hill_diam_m/2 + 5
    cx <- runif(random_hills, xmin + margin, xmax - margin)
    cy <- runif(random_hills, ymin + margin, ymax - margin)
    centers <- rbind(centers, cbind(cx, cy))
  }
  w_h_micro <- rep(0, nrow(XY))
  if (!is.null(centers) && nrow(centers) > 0) {
    hr_m <- max(1e-6, micro_hill_diam_m/2)
    for (i in seq_len(nrow(centers))) {
      r  <- sqrt((XY$x - centers[i,1])^2 + (XY$y - centers[i,2])^2)
      wi <- if (micro_hill_smooth) exp(-(r/hr_m)^2) else as.numeric(r <= hr_m)
      w_h_micro <- w_h_micro + wi
    }
    w_h_micro <- pmin(1, w_h_micro)
    elev <- elev + micro_hill_height_m * w_h_micro
  }

  # Raster
  E     <- Rtemplate; terra::values(E) <- elev; names(E) <- "elev"
  lakeR <- Rtemplate; terra::values(lakeR) <- if (lake_mode=="water") as.numeric(w_l>0) else 0; names(lakeR) <- "lake"
  hillW_main  <- Rtemplate; terra::values(hillW_main)  <- w_h_main;  names(hillW_main)  <- "hillW"
  hillW_micro <- Rtemplate; terra::values(hillW_micro) <- w_h_micro; names(hillW_micro) <- "hillW"
  hillW <- terra::clamp(hillW_main + hillW_micro, lower = 0, upper = 1); names(hillW) <- "hillW"

  # Gelände & Sonne
  slp  <- terra::terrain(E, v="slope",  unit="radians")
  asp  <- terra::terrain(E, v="aspect", unit="radians")
  slp0 <- terra::ifel(is.na(slp), 0, slp); asp0 <- terra::ifel(is.na(asp), 0, asp)
  I14  <- cosi_fun(sun14$alt, sun14$az, slp0, asp0)
  I05  <- cosi_fun(sun05$alt, sun05$az, slp0, asp0)

  # Landnutzung (1 forest, 2 water, 3 bare soil, 4 maize)
  lc <- Rtemplate; terra::values(lc) <- 4L
  lc <- terra::ifel(lakeR > 0, 2L, lc)
  forest_mask <- (hillW > 0.2) | (slp0 > 0.15 & (terra::init(Rtemplate,"y") > y0))
  lc <- terra::ifel((forest_mask) & (lakeR <= 0), 1L, lc)
  v_slp   <- terra::values(slp0); thr_slp <- stats::quantile(v_slp[is.finite(v_slp)], 0.90, na.rm = TRUE)
  bare_mask <- (slp0 >= thr_slp) & (lakeR <= 0) & (!forest_mask)
  lc <- terra::ifel(bare_mask, 3L, lc); lc <- terra::clamp(lc, lower = 1, upper = 4); names(lc) <- "lc"

  # Physikfelder
  phys <- build_physics_fields(list(E=E, slp=slp0, I14=I14, I05=I05, hillW=hillW), lc, noise14, noise05)

  # NAs robust ersetzen (für metrische Raster)
  fix_nonfinite <- function(r) { v <- terra::values(r); m <- stats::median(v[is.finite(v)], na.rm = TRUE)
    v[!is.finite(v)] <- m; terra::values(r) <- v; r }

  # >>> FIX: Levels & Farben an Szenario anhängen (für Stationen & Plots)
  lc_levels <- lc_levels_default
  lc_colors <- lc_colors_default

  list(E = fix_nonfinite(E), R14 = fix_nonfinite(phys$R14), R05 = fix_nonfinite(phys$R05),
       lake = lakeR, hillW = hillW, slp = slp0, asp = asp0, I14 = I14, I05 = I05,
       lc = lc,
       lc_levels = lc_levels, lc_colors = lc_colors)
}

# Build the scenario ----------------------------------------------------
  scen <- build_scenario(lake_mode = lake_choice, hill_mode = hill_choice,random_hills = no_random_hill,micro_hill_diam_m = micro_hill_diam_m,micro_hill_height_m = micro_hill_height_m,micro_hill_smooth = micro_hill_smooth)


# ========== Landnutzungs- & Gelände-Visualizer ==========
# ------------------------------------------------------------------------------
# plot_landcover_terrain(scen, stations = NULL, show_contours = TRUE)
# Purpose:
#   Visualize land cover, elevation, and slope rasters side-by-side. Optionally
#   overlay station points and draw hill/lake contours for orientation.
# Inputs:
#   - scen: list returned by build_scenario(), must contain E, slp, lc, lake, hillW.
#   - stations: optional sf point layer to overlay.
#   - show_contours: toggle for drawing 0.5 level contours around lake/hill masks.
# Returns:
#   A patchwork-composed ggplot object (no side effects).
# ------------------------------------------------------------------------------
plot_landcover_terrain <- function(scen, stations = NULL, show_contours = TRUE,
                                   layout = c("grid","vertical")) {
  layout <- match.arg(layout)
  stopifnot(all(c("E","slp") %in% names(scen)))

  # Wenn LC fehlt (sollte jetzt nicht passieren), Fall  back bauen
  if (!("lc" %in% names(scen))) {
    E    <- scen$E; slp <- scen$slp
    lake <- if ("lake"  %in% names(scen)) scen$lake  else terra::rast(E)*0
    hill <- if ("hillW" %in% names(scen)) scen$hillW else terra::rast(E)*0
    y0loc <- (terra::ymax(E)+terra::ymin(E))/2
    slp_vals <- terra::values(slp)
    thr_slp  <- stats::quantile(slp_vals[is.finite(slp_vals)], 0.90, na.rm = TRUE)
    lc_fallback <- terra::rast(E); terra::values(lc_fallback) <- 4L  # maize
    lc_fallback <- terra::ifel(lake > 0, 2L, lc_fallback)            # water
    forest_mask <- (hill > 0.2) | (slp > 0.15 & (terra::init(E,"y") > y0loc))
    lc_fallback <- terra::ifel((forest_mask) & (lake <= 0), 1L, lc_fallback)  # forest
    bare_mask <- (slp >= thr_slp) & (lake <= 0) & (!forest_mask)
    lc_fallback <- terra::ifel(bare_mask, 3L, lc_fallback)           # bare soil
    names(lc_fallback) <- "lc"
    scen$lc <- lc_fallback
    lc_levels <- lc_levels_default
    scen$lc_colors <- lc_colors_default
  }

  lc_df  <- as.data.frame(scen$lc,  xy = TRUE); names(lc_df)  <- c("x","y","lc")
  E_df   <- as.data.frame(scen$E,   xy = TRUE); names(E_df)   <- c("x","y","elev")
  slp_df <- as.data.frame(scen$slp, xy = TRUE); names(slp_df) <- c("x","y","slp")

  lc_df$lc <- factor(lc_df$lc, levels = seq_along(lc_levels), labels = lc_levels)
  cols_lc  <- scen$lc_colors

  p_lc <- ggplot() +
    geom_raster(data = lc_df, aes(x, y, fill = lc)) +
    scale_fill_manual(values = cols_lc, na.value = "grey90", name = "Landuse") +
    coord_equal() + theme_minimal() +
    labs(title = "Landuse", x = "Easting", y = "Northing")

  p_elev <- ggplot() +
    geom_raster(data = E_df, aes(x, y, fill = elev)) +
    scale_fill_viridis_c(name = "Altitude [m]") +
    coord_equal() + theme_minimal() +
    labs(title = "Altitude", x = "Easting", y = "Northing")

  p_slp <- ggplot() +
    geom_raster(data = slp_df, aes(x, y, fill = slp)) +
    scale_fill_viridis_c(name = "Slope [rad]") +
    coord_equal() + theme_minimal() +
    labs(title = "Slope", x = "Easting", y = "Northing")

  if (isTRUE(show_contours)) {
    lake_df <- as.data.frame(scen$lake, xy = TRUE); names(lake_df) <- c("x","y","lake")
    hill_df <- as.data.frame(scen$hillW, xy = TRUE); names(hill_df) <- c("x","y","hillW")
    p_elev <- p_elev + geom_contour(data = E_df, aes(x, y, z = elev),
                                    bins = 10, colour = "black", alpha = 0.25,
                                    linewidth = 0.2, inherit.aes = FALSE)
    p_lc  <- p_lc  + geom_contour(data = lake_df, aes(x, y, z = lake),
                                  breaks = 0.5, colour = "black", linewidth = 0.35,
                                  inherit.aes = FALSE) +
                    geom_contour(data = hill_df, aes(x, y, z = hillW),
                                  breaks = 0.5, colour = "black", linetype = "22",
                                  linewidth = 0.3, inherit.aes = FALSE)
    p_slp <- p_slp + geom_contour(data = lake_df, aes(x, y, z = lake),
                                  breaks = 0.5, colour = "black", linewidth = 0.35,
                                  inherit.aes = FALSE) +
                    geom_contour(data = hill_df, aes(x, y, z = hillW),
                                  breaks = 0.5, colour = "black", linetype = "22",
                                  linewidth = 0.3, inherit.aes = FALSE)
  }

  if (!is.null(stations)) {
    add_st <- list(geom_sf(data = stations, colour = "black", fill = "white",
                           shape = 21, size = 1.6, stroke = 0.25, inherit.aes = FALSE))
    p_lc   <- p_lc   + add_st
    p_elev <- p_elev + add_st
    p_slp  <- p_slp  + add_st
  }
  if (layout == "vertical") {
    # stapeln, Legenden behalten
    return(
      (p_lc   + ggplot2::theme(plot.margin = ggplot2::margin(4,4,4,4, unit = "pt"))) /
        (p_elev + ggplot2::theme(plot.margin = ggplot2::margin(4,4,4,4, unit = "pt"))) /
        (p_slp  + ggplot2::theme(plot.margin = ggplot2::margin(4,4,4,4, unit = "pt"))) +
        patchwork::plot_layout(guides = "keep", heights = c(1,1,1))
    )
  } else {
    return((p_lc | (p_elev | p_slp)) + patchwork::plot_layout(guides = "keep"))
  }
}
# ---- Aufruf (direkt nach build_scenario):
plot_landcover_terrain(scen = scen,layout = "vertical")

# ---- 2x2 overview panel in English -----------------------------------------
plot_block_overview_2x2_en <- function(scen, pts_sf = NULL) {
  stopifnot(all(c("E","slp","I14","I05") %in% names(scen)))
  
  # Stack continuous rasters -> data.frame
  Rstack <- c(scen$E, scen$slp, scen$I14, scen$I05)
  df <- terra::as.data.frame(Rstack, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x","y","elev","slope","I14","I05")
  
  # Base theme & palettes (no extra packages)
  theme_base <- ggplot2::theme_minimal(base_size = 11)
  pal_terrain <- grDevices::hcl.colors(256, "Terrain")
  pal_slope   <- grDevices::hcl.colors(256, "Viridis")
  pal_hot     <- grDevices::hcl.colors(256, "YlOrRd")
  pal_cool    <- grDevices::hcl.colors(256, "PuBuGn")
  
  p_elev <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = elev)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_terrain, name = "m") +
    ggplot2::labs(title = "Terrain (Elevation)") + theme_base
  
  p_slope <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = slope)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_slope, name = "rad") +
    ggplot2::labs(title = "Slope (radians)") + theme_base
  
  p_I14 <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = I14)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_hot, name = "") +
    ggplot2::labs(title = "Insolation 14 UTC (cos i)") + theme_base
  
  p_I05 <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = I05)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_cool, name = "") +
    ggplot2::labs(title = "Insolation 05 UTC (cos i)") + theme_base
  
  # Land cover (if available); falls back to slope otherwise
  p_lc <- NULL
  if ("lc" %in% names(scen)) {
    lc_df <- terra::as.data.frame(scen$lc, xy = TRUE, na.rm = FALSE)
    names(lc_df) <- c("x","y","lc")
    lc_levels <- scen$lc_levels %||% scen$lc_lev
    if (!is.null(lc_levels)) {
      lc_df$lc <- factor(lc_df$lc, levels = seq_along(lc_levels), labels = lc_levels)
    }
    lc_cols <- scen$lc_colors %||% c("Forest"="#2E8B57","Water"="#5DADE2",
                                     "Bare soil"="#C49A6C","Maize"="#F4D03F")
    p_lc <- ggplot2::ggplot(lc_df, ggplot2::aes(x, y, fill = lc)) +
      ggplot2::geom_raster() + ggplot2::coord_equal() +
      ggplot2::scale_fill_manual(values = lc_cols, drop = FALSE, name = "Land cover") +
      ggplot2::labs(title = "Land cover") + theme_base
  }
  
  # Optional station overlay
  if (!is.null(pts_sf)) {
    pts_df <- sf::st_drop_geometry(pts_sf)
    add_pts <- function(p)
      p + ggplot2::geom_point(data = pts_df, ggplot2::aes(x = x, y = y),
                              inherit.aes = FALSE, size = 0.7, alpha = 0.7,
                              colour = "black")
    p_elev  <- add_pts(p_elev)
    p_slope <- add_pts(p_slope)
    p_I14   <- add_pts(p_I14)
    p_I05   <- add_pts(p_I05)
    if (!is.null(p_lc)) p_lc <- add_pts(p_lc)
  }
  
  # Choose 4th panel: land cover if present, otherwise slope
  p2 <- if (!is.null(p_lc)) p_lc else p_slope
  
  # 2x2 layout
  # 2x2 layout + rotate x-axis labels 45°
  (p_elev | p2) / (p_I14 | p_I05) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
      plot.margin = ggplot2::margin(5,5,5,5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
}




# (Optionale zweite Darstellung – Tippfehler korrigiert)
lc_df <- as.data.frame(scen$lc, xy = TRUE)
names(lc_df) <- c("x","y","lc")
lc_df$lc <- factor(lc_df$lc, levels = 1:length(lc_levels), labels = lc_levels)

ggplot(lc_df, aes(x, y, fill = lc)) +
  geom_raster() +
  scale_fill_manual(values = scen$lc_colors, drop = FALSE) +
  coord_equal() + theme_minimal() +
  labs(title = "Land-cover (classes)", x = "Easting", y = "Northing", fill = "LC")


# Station coordinates ---------------------------------------------------
if (station_mode == "random") {
  pts <- tibble::tibble(
    id = 1:n_st,
    x  = runif(n_st, xmin + transect_margin_m, xmax - transect_margin_m),
    y  = runif(n_st, ymin + transect_margin_m, ymax - transect_margin_m)
  )
} else if (station_mode == "ns_transect") {
  x_const <- min(max(x0 + ns_offset_m, xmin + transect_margin_m), xmax - transect_margin_m)
  y_seq   <- seq(ymin + transect_margin_m, ymax - transect_margin_m, length.out = n_st)
  pts <- tibble::tibble(id = 1:n_st, x = x_const, y = y_seq)
} else if (station_mode == "ew_transect") {
  y_const <- min(max(y0 + ew_offset_m, ymin + transect_margin_m), ymax - transect_margin_m)
  x_seq   <- seq(xmin + transect_margin_m, xmax - transect_margin_m, length.out = n_st)
  pts <- tibble::tibble(id = 1:n_st, x = x_seq, y = y_const)
} else stop("Unknown station_mode")

pts_sf <- sf::st_as_sf(pts, coords = c("x","y"), crs = crs_utm, remove = FALSE)
vpts   <- terra::vect(pts_sf)

# Extract covariates & targets -----------------------------------------
pts$z_surf <- as.numeric(terra::extract(scen$E,   vpts, ID = FALSE)[,1])
pts$slp    <- as.numeric(terra::extract(scen$slp, vpts, ID = FALSE)[,1])
pts$I14    <- as.numeric(terra::extract(scen$I14, vpts, ID = FALSE)[,1])
pts$I05    <- as.numeric(terra::extract(scen$I05, vpts, ID = FALSE)[,1])
pts$lc     <- as.integer(terra::extract(scen$lc,  vpts, ID = FALSE)[,1])
pts$T14    <- as.numeric(terra::extract(scen$R14, vpts, ID = FALSE)[,1])
pts$T05    <- as.numeric(terra::extract(scen$R05, vpts, ID = FALSE)[,1])

# Keep complete rows per time slot -------------------------------------
pts14 <- pts[stats::complete.cases(pts[, c("x","y","z_surf","slp","I14","lc","T14")]), ]
pts05 <- pts[stats::complete.cases(pts[, c("x","y","z_surf","slp","I05","lc","T05")]), ]

# Unify response name to 'temp' and carry LC as factor -----------------
# Unify response 'temp' and carry LC as factor with global levels ---------------
stn_sf_14 <- pts14 |>
  dplyr::transmute(
    id, x, y,
    z_surf = as.numeric(z_surf),
    slp    = as.numeric(slp),
    cosi   = as.numeric(I14),
    lc     = factor(lc_levels[pmax(1, pmin(lc, length(lc_levels)))], levels = lc_levels),
    temp   = as.numeric(T14)
  ) |>
  sf::st_as_sf(coords = c("x","y"), crs = crs_utm, remove = FALSE)

stn_sf_05 <- pts05 |>
  dplyr::transmute(
    id, x, y,
    z_surf = as.numeric(z_surf),
    slp    = as.numeric(slp),
    cosi   = as.numeric(I05),
    lc     = factor(lc_levels[pmax(1, pmin(lc, length(lc_levels)))], levels = lc_levels),
    temp   = as.numeric(T05)
  ) |>
  sf::st_as_sf(coords = c("x","y"), crs = crs_utm, remove = FALSE)

# With station overlay (e.g., T14)
plot_block_overview_2x2_en(scen, pts_sf = stn_sf_14)

#| label: lbo-helpers
#| echo: false
#| message: false
#| warning: false
#| results: 'asis'

# — Compute metrics from a CV table (fallback) —
lbo_compute_metrics <- function(cv_tbl) {
  stopifnot(all(c("model","obs","pred") %in% names(cv_tbl)))
  dplyr::group_by(cv_tbl, model) |>
    dplyr::summarise(
      n    = dplyr::n(),
      n_ok = sum(is.finite(obs) & is.finite(pred)),
      MAE  = mean(abs(pred - obs), na.rm = TRUE),
      RMSE = sqrt(mean((pred - obs)^2, na.rm = TRUE)),
      Bias = mean(pred - obs, na.rm = TRUE),
      R2   = {
        m <- mean(obs, na.rm=TRUE)
        1 - sum((pred-obs)^2, na.rm=TRUE) / sum((obs-m)^2, na.rm=TRUE)
      },
      .groups = "drop"
    ) |>
    dplyr::arrange(RMSE)
}

# — Show metrics from df or file —
lbo_show_metrics <- function(metrics_df = NULL, csv = NULL, html = NULL, caption = "LBO-CV metrics — baseline") {
  if (!is.null(metrics_df)) {
    return(knitr::kable(metrics_df, format="html", digits=3, caption = caption))
  }
  if (!is.null(html) && file.exists(html)) {
    return(htmltools::includeHTML(html))
  }
  if (!is.null(csv) && file.exists(csv)) {
    df <- read.csv(csv, check.names = FALSE)
    return(knitr::kable(df, format="html", digits=3, caption = caption))
  }
  cat("<em>metrics not found</em>")
}

# — Print a “panel” object robustly (patchwork or list-of-pages) —
lbo_print_panel <- function(panel_obj = NULL, page = 1, panel_png = NULL) {
  if (!is.null(panel_obj)) {
    # panel kann patchwork sein ODER Liste von Seiten (z.B. $pages)
    if (is.list(panel_obj) && !is.null(panel_obj$pages)) {
      print(panel_obj$pages[[page]])
      return(invisible(TRUE))
    }
    if (is.list(panel_obj) && all(vapply(panel_obj, inherits, logical(1), what="gg"))) {
      print(panel_obj[[page]])
      return(invisible(TRUE))
    }
    print(panel_obj)
    return(invisible(TRUE))
  }
  if (!is.null(panel_png) && file.exists(panel_png)) {
    cat(sprintf("![](%s){width=100%%}\n", panel_png))
    return(invisible(TRUE))
  }
  cat("<em>panel not available</em>\n")
  invisible(FALSE)
}



# ===================== SF-ONLY GEOSTATS PATCH =====================

# ---- Helpers ------------------------------------------------------
.align_factor_pair <- function(train_x, grid_x, fallback = NULL) {
  tl <- unique(na.omit(as.character(train_x)))
  if (!length(tl)) return(list(use = FALSE, train = NULL, grid = NULL, levels = character()))
  lev <- if (is.null(fallback)) tl else unique(c(tl, fallback))
  trc <- factor(as.character(train_x), levels = lev)
  gdc <- factor(as.character(grid_x),  levels = lev)
  list(use = TRUE, train = trc, grid = gdc, levels = lev)
}
.align_factor_to_model <- function(x, lev_model) {
  y <- factor(as.character(x), levels = lev_model)
  if (anyNA(y)) y[is.na(y)] <- lev_model[1]
  y
}
.fill_num_na_vec <- function(x, ref) {
  m <- stats::median(ref[is.finite(ref)], na.rm = TRUE)
  x[!is.finite(x)] <- m
  x
}
.default_vgm <- function(values, model = "Exp", range = 100) {
  psill <- stats::var(values, na.rm = TRUE); nug <- 0.1 * psill
  gstat::vgm(psill = psill, model = model, range = range, nugget = nug)
}

# ---- SF Learner: Voronoi / Next Neighbour -------------------------
pred_Voronoi <- function(train_sf, test_sf) {
  idx <- sf::st_nearest_feature(test_sf, train_sf)
  as.numeric(train_sf$temp)[idx]
}

# ---- SF Learner: IDW ----------------------------------------------
pred_IDW <- function(train_sf, test_sf, idp = 2) {
  # train_sf: sf POINTS mit Spalte 'temp'
  # test_sf : sf POINTS, beliebige Zusatzspalten
  pr <- suppressWarnings(gstat::idw(temp ~ 1, locations = train_sf, newdata = test_sf, idp = idp))
  as.numeric(pr$var1.pred)
}

# ---- SF Learner: Ordinary Kriging ---------------------------------
pred_OK <- function(train_sf, test_sf) {
  vg      <- suppressWarnings(gstat::variogram(temp ~ 1, data = train_sf))
  vgm_fit <- try(suppressWarnings(gstat::fit.variogram(vg, gstat::vgm("Exp"))), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) vgm_fit <- .default_vgm(train_sf$temp)
  pr <- suppressWarnings(gstat::krige(temp ~ 1, locations = train_sf, newdata = test_sf, model = vgm_fit))
  as.numeric(pr$var1.pred)
}

# ---- SF Learner: KED (UK mit externen Drifts) ---------------------
# --- sf-only KED (schluckt extra Args wie E=E) ----------------------
pred_KED <- function(train_sf, test_sf, ...) {
  stopifnot(inherits(train_sf, "sf"), inherits(test_sf, "sf"))
  need <- c("z_surf","slp","cosi")
  miss <- setdiff(need, names(train_sf))
  if (length(miss)) stop("pred_KED(): fehlende Drift-Spalten im Training: ",
                         paste(miss, collapse = ", "))
  
  # Optional LC als Faktor angleichen
  use_lc <- "lc" %in% names(train_sf) && "lc" %in% names(test_sf)
  tr <- train_sf
  te <- test_sf
  if (use_lc) {
    tr$lc <- droplevels(factor(tr$lc))
    te$lc <- factor(as.character(te$lc), levels = levels(tr$lc))
    te$lc[is.na(te$lc)] <- levels(tr$lc)[1]
  }
  
  # fehlende numerische Drifts im TEST mit Trainingsmedian auffüllen
  for (nm in need) {
    m <- stats::median(tr[[nm]][is.finite(tr[[nm]])], na.rm = TRUE)
    te[[nm]][!is.finite(te[[nm]])] <- m
  }
  
  # nur vollständige Trainingszeilen
  keep_tr <- c("temp", need, if (use_lc) "lc")
  dtr <- sf::st_drop_geometry(tr)[, keep_tr, drop = FALSE]
  ok  <- stats::complete.cases(dtr)
  tr  <- tr[ok, ]
  if (nrow(tr) < 5) return(rep(NA_real_, nrow(te)))
  
  # Formel: lineare Drifts + optional LC
  form <- stats::as.formula(paste("temp ~", paste(c(need, if (use_lc) "lc"), collapse = " + ")))
  
  # Variogramm + robuster Fit
  vg      <- suppressWarnings(gstat::variogram(form, data = tr))
  vgm_fit <- try(suppressWarnings(gstat::fit.variogram(vg, gstat::vgm("Exp"))), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) {
    ps <- stats::var(sf::st_drop_geometry(tr)$temp, na.rm = TRUE)
    vgm_fit <- gstat::vgm(psill = ps, model = "Exp", range = max(vg$dist, na.rm = TRUE)/3, nugget = 0.1*ps)
  }
  
  # Kriging mit externen Drifts (UK/KED)
  pr <- suppressWarnings(gstat::krige(form, locations = tr, newdata = te, model = vgm_fit))
  as.numeric(pr$var1.pred)
}


pred_RF <- function(train_sf, test_sf) {
  dtr <- sf::st_drop_geometry(train_sf)
  if (!("lc" %in% names(dtr))) dtr$lc <- factor("Wald", levels = lc_levels)
  dtr$lc <- droplevels(factor(as.character(dtr$lc), levels = lc_levels))
  dtr <- stats::na.omit(dtr)
  if (nrow(dtr) < 5) return(rep(NA_real_, nrow(test_sf)))
  rf  <- randomForest::randomForest(temp ~ x + y + z_surf + slp + cosi + lc, data = dtr, na.action = na.omit)
  
  dte <- sf::st_drop_geometry(test_sf)
  if (!("lc" %in% names(dte))) dte$lc <- factor("Wald", levels = lc_levels)
  lev <- levels(dtr$lc)
  dte$lc <- .align_factor_to_model(dte$lc, lev)
  
  good <- stats::complete.cases(dte[, c("x","y","z_surf","slp","cosi","lc")])
  out  <- rep(NA_real_, nrow(dte)); if (any(good)) out[good] <- stats::predict(rf, dte[good, ])
  out
}

# GAM-Fitter (nutzt safe_gam_formula)
fit_gam_safe <- function(stn_sf) {
  d <- stn_sf |> sf::st_drop_geometry()
  d <- d[stats::complete.cases(d[, c("x","y","temp","z_surf","slp","cosi")]), , drop = FALSE]
  if (nrow(d) < 10) stop("Too few stations for GAM: n=", nrow(d))
  mgcv::gam(formula = safe_gam_formula(d), data = d, method = "REML", select = TRUE)
}

pred_GAM <- function(train_sf, test_sf) {
  # training data (keep only columns we’ll use, drop NAs)
  dtr  <- sf::st_drop_geometry(train_sf)
  keep <- intersect(c("temp","x","y","z_surf","slp","cosi","lc"), names(dtr))
  dtr  <- dtr[stats::complete.cases(dtr[, keep, drop = FALSE]), keep, drop = FALSE]
  if (!nrow(dtr)) return(rep(NA_real_, nrow(test_sf)))
  
  # land cover: optionally include as factor main effect
  if ("lc" %in% names(dtr)) dtr$lc <- droplevels(factor(dtr$lc))
  inc_lc <- "lc" %in% names(dtr) && nlevels(dtr$lc) >= 2
  
  # fit GAM with guarded formula (dynamic k, optional terms)
  if (nrow(dtr) < 10) return(rep(NA_real_, nrow(test_sf)))
  gm <- mgcv::gam(
    formula = safe_gam_formula(dtr, include_lc = inc_lc),
    data    = dtr,
    method  = "REML",
    select  = TRUE
  )
  
  # prediction data (+ align LC levels to model, if used)
  dte <- sf::st_drop_geometry(test_sf)
  vars <- c("x","y","z_surf","slp","cosi", if (inc_lc) "lc")
  vars <- intersect(vars, names(dte))
  
  if (inc_lc) {
    # align factor levels; assumes you have .align_factor_to_model()
    lev <- levels(model.frame(gm)$lc)
    if (!("lc" %in% names(dte))) dte$lc <- lev[1]
    dte$lc <- .align_factor_to_model(dte$lc, lev)
  }
  
  good <- stats::complete.cases(dte[, vars, drop = FALSE])
  out  <- rep(NA_real_, nrow(dte))
  if (any(good)) {
    out[good] <- stats::predict(gm, dte[good, vars, drop = FALSE], type = "response")
  }
  out
}



# ---- SF-only predict_maps (kein sp, kein grid_sp) -------------------
# ---------- SF-only predict_maps ----------
predict_maps <- function(stn_sf, truth_raster,
                         which_time = c("T14","T05"),
                         scen, models = c("Voronoi","IDW","OK","KED","RF","GAM"),
                         lc_levels = NULL,
                         feature_rasters = NULL) {
  which_time <- match.arg(which_time)
  lc_levels  <- lc_levels %||% scen$lc_levels
  
  # Basis-/Feature-Raster setzen
  E      <- scen$E
  slp_r  <- scen$slp
  cosi_r <- if (which_time == "T14") scen$I14 else scen$I05
  if (!is.null(feature_rasters)) {
    if (!is.null(feature_rasters$E))     E     <- feature_rasters$E
    if (!is.null(feature_rasters$slp))   slp_r <- feature_rasters$slp
    if (!is.null(feature_rasters$cosi))  cosi_r<- feature_rasters$cosi
  }
  has_lc <- ("lc" %in% names(scen)) && !is.null(scen$lc)
  lc_r   <- if (has_lc) scen$lc else NULL
  
  # Trainingsdaten vervollständigen
  train_sf <- stn_sf
  if (!all(c("x","y") %in% names(train_sf))) {
    xy <- sf::st_coordinates(train_sf); train_sf$x <- xy[,1]; train_sf$y <- xy[,2]
  }
  if (!("z_surf" %in% names(train_sf)))
    train_sf$z_surf <- as.numeric(terra::extract(E,      sf::st_coordinates(train_sf))[,1])
  if (!("slp" %in% names(train_sf)))
    train_sf$slp    <- as.numeric(terra::extract(slp_r,  sf::st_coordinates(train_sf))[,1])
  if (!("cosi" %in% names(train_sf)))
    train_sf$cosi   <- as.numeric(terra::extract(cosi_r, sf::st_coordinates(train_sf))[,1])
  if (has_lc && !("lc" %in% names(train_sf))) {
    lc_codes <- as.integer(terra::extract(lc_r, sf::st_coordinates(train_sf))[,1])
    lc_codes[is.na(lc_codes)] <- 1L
    lc_codes <- pmax(1L, pmin(lc_codes, length(lc_levels)))
    train_sf$lc <- factor(lc_levels[lc_codes], levels = lc_levels)
  }
  
  # Vorhersage-Grid 1:1 zu Rasterzellen (nie filtern)
  xy <- as.data.frame(terra::xyFromCell(E, 1:terra::ncell(E))); names(xy) <- c("x","y")
  grid_df <- xy
  grid_df$z_surf <- as.numeric(terra::values(E))
  grid_df$slp    <- as.numeric(terra::values(slp_r))
  grid_df$cosi   <- as.numeric(terra::values(cosi_r))
  if (has_lc) {
    lc_codes <- as.integer(terra::values(lc_r))
    lc_codes[!is.finite(lc_codes)] <- 1L
    lc_codes <- pmax(1L, pmin(lc_codes, length(lc_levels)))
    grid_df$lc <- factor(lc_levels[lc_codes], levels = lc_levels)
  }
  grid_sf <- sf::st_as_sf(grid_df, coords = c("x","y"),
                          crs = sf::st_crs(train_sf), remove = FALSE)
  
  # LC-Levels train/grid angleichen (hart, ohne __OTHER__)
  use_lc <- has_lc && ("lc" %in% names(train_sf)) && ("lc" %in% names(grid_sf))
  if (use_lc) {
    lev <- levels(droplevels(factor(train_sf$lc)))
    train_sf$lc <- factor(as.character(train_sf$lc), levels = lev)
    grid_sf$lc  <- factor(as.character(grid_sf$lc),  levels = lev)
    if (anyNA(train_sf$lc) || anyNA(grid_sf$lc)) {
      use_lc <- FALSE; train_sf$lc <- NULL; grid_sf$lc <- NULL
    }
  }
  
  # Modelle laufen lassen
  pred_list <- list()
  if ("Voronoi" %in% models) pred_list$Voronoi <- pred_Voronoi(train_sf, grid_sf)
  if ("IDW"     %in% models) pred_list$IDW     <- pred_IDW    (train_sf, grid_sf, idp = 2)
  if ("OK"      %in% models) pred_list$OK      <- pred_OK     (train_sf, grid_sf)
  if ("KED"     %in% models) pred_list$KED     <- pred_KED    (train_sf, grid_sf)
  
  if ("RF" %in% models) {
    dtr <- sf::st_drop_geometry(train_sf)
    rf_vars <- c("x","y","z_surf","slp","cosi", if (use_lc) "lc")
    dtr <- stats::na.omit(dtr[, c("temp", rf_vars), drop = FALSE])
    pred_list$RF <- if (nrow(dtr) >= 5) {
      rf <- randomForest::randomForest(
        stats::as.formula(paste("temp ~", paste(rf_vars, collapse = " + "))),
        data = dtr, na.action = na.omit
      )
      as.numeric(stats::predict(rf, sf::st_drop_geometry(grid_sf)[, rf_vars, drop = FALSE]))
    } else rep(NA_real_, nrow(grid_sf))
  }
  
  if ("GAM" %in% models) {
    dtr <- sf::st_drop_geometry(train_sf)
    keep <- c("temp","x","y","z_surf","slp","cosi", if (use_lc) "lc")
    dtr  <- dtr[stats::complete.cases(dtr[, keep, drop = FALSE]), keep, drop = FALSE]
    if (nrow(dtr) >= 10) {
      form <- safe_gam_formula(dtr, include_lc = use_lc)
      gm   <- mgcv::gam(form, data = dtr, method = "REML", select = TRUE)
      vars_needed <- setdiff(all.vars(formula(gm)), "temp")
      nd <- sf::st_drop_geometry(grid_sf)[, vars_needed, drop = FALSE]
      mf <- try(model.frame(gm), silent = TRUE)
      if (!inherits(mf, "try-error")) {
        for (vn in vars_needed) if (is.factor(mf[[vn]])) {
          nd[[vn]] <- factor(as.character(nd[[vn]]), levels = levels(mf[[vn]]))
          na_idx <- is.na(nd[[vn]]); if (any(na_idx)) nd[[vn]][na_idx] <- levels(mf[[vn]])[1L]
        }
      }
      good <- stats::complete.cases(nd)
      tmp  <- rep(NA_real_, nrow(grid_sf))
      if (any(good)) tmp[good] <- stats::predict(gm, nd[good, , drop = FALSE], type = "response")
      pred_list$GAM <- tmp
    } else pred_list$GAM <- rep(NA_real_, nrow(grid_sf))
  }
  
  # Ausgaben
  pred_df <- dplyr::bind_rows(lapply(names(pred_list), function(nm) {
    tibble::tibble(model = nm, x = grid_df$x, y = grid_df$y, pred = pred_list[[nm]])
  }))
  
  make_r <- function(vals, template = E) {
    stopifnot(length(vals) == terra::ncell(template))
    r <- terra::rast(template); terra::values(r) <- as.numeric(vals); r
  }
  pred_rasters <- lapply(pred_list, make_r)
  
  truth_df <- as.data.frame(truth_raster, xy = TRUE, na.rm = FALSE)
  names(truth_df) <- c("x","y","truth")
  lims <- stats::quantile(truth_df$truth, probs = stretch_q, na.rm = TRUE)
  
  p_pred <- ggplot2::ggplot(pred_df, ggplot2::aes(x, y, fill = pred)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colors = temp_palette, limits = lims,
                                  oob = scales::squish, name = "Temp") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Predictions by model", which_time),
                  x = "Easting", y = "Northing") +
    ggplot2::facet_wrap(~ model, ncol = 3)
  
  p_truth <- ggplot2::ggplot(truth_df, ggplot2::aes(x, y, fill = truth)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colors = temp_palette, limits = lims,
                                  oob = scales::squish, name = "Temp") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Truth raster", which_time),
                  x = "Easting", y = "Northing")
  
  list(
    pred_df      = pred_df,
    pred_rasters = pred_rasters,
    p_pred       = p_pred,
    p_truth      = p_truth
  )
}


# ===================== Ende SF-Patch =====================



# Quick station table ---------------------------------------------------
pts |> dplyr::transmute(
  id, easting = round(x), northing = round(y),
  z_surf = round(z_surf,1), LC = lc_levels[lc],
  T14_C = round(T14,1), T05_C = round(T05,1)
) |> knitr::kable(caption = "Station sample (sanity check)", digits = 1)

 


assemble_simple_panel <- function(maps, cv_tbl, which_time = c("T14","T05")) {
  which_time <- match.arg(which_time)
  
  # Modelle in EINER Reihe
  p_models_row <- maps$p_pred +
    ggplot2::facet_wrap(~model, nrow = 1) +
    #ggplot2::labs(title = sprintf("%s — Predictions by model", which_time))
  
  # Eine Truth-Karte
  p_truth <- maps$p_truth + ggplot2::labs(title = sprintf("%s — truth", which_time))
  # 
  # Fehler-Panels
  p_block <- make_block_metric_box(cv_tbl, which_time)
  p_abs   <- make_abs_error_box(cv_tbl, which_time)
  p_dens  <- make_residual_density(cv_tbl, which_time)
  
  # Layout
  (      p_truth /
      p_models_row /
       p_block/
        p_abs /
      p_dens
  ) + patchwork::plot_layout(heights = c(1.1, 0.9, 1, 0.8))
}
#----------


# sorgt dafür, dass X_tr und X_te exakt die gleichen Spaltennamen/Anzahl haben
.reconcile_mm_cols <- function(X_tr, X_te) {
  miss_te <- setdiff(colnames(X_tr), colnames(X_te))
  if (length(miss_te)) {
    X_te <- cbind(
      X_te,
      matrix(0, nrow = nrow(X_te), ncol = length(miss_te),
             dimnames = list(NULL, miss_te))
    )
  }
  miss_tr <- setdiff(colnames(X_te), colnames(X_tr))
  if (length(miss_tr)) {
    X_tr <- cbind(
      X_tr,
      matrix(0, nrow = nrow(X_tr), ncol = length(miss_tr),
             dimnames = list(NULL, miss_tr))
    )
  }
  # identische Spaltenreihenfolge
  X_te <- X_te[, colnames(X_tr), drop = FALSE]
  list(X_tr = X_tr, X_te = X_te)
}

safe_r2 <- function(obs, pred) {
  idx <- is.finite(obs) & is.finite(pred)
  if (sum(idx) < 2) return(NA_real_)
  x <- obs[idx]; y <- pred[idx]
  sx <- stats::sd(x); sy <- stats::sd(y)
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
  stats::cor(x, y)^2
}

.k_for_xy <- function(n, n_xy) max(3, min(60, n_xy - 1L, floor(n * 0.8)))
.kcap_unique <- function(x, kmax) {
  ux <- unique(x[is.finite(x)])
  nu <- length(ux)
  if (nu <= 3) return(0L)                # treat as constant/near-constant
  max(4L, min(kmax, nu - 1L))
}

# Error/diagnostic panels ----------------------------------------------
block_metrics_long <- function(cv_tbl) {
  stopifnot(all(c("model","block_id","obs","pred") %in% names(cv_tbl)))
  cv_tbl |>
    dplyr::group_by(model, block_id) |>
    dplyr::summarise(RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)), MAE = mean(abs(obs - pred), na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_longer(c(RMSE, MAE), names_to = "Metric", values_to = "Value")
}

order_models_by_median_rmse <- function(cv_tbl) {
  bm <- block_metrics_long(cv_tbl)
  bm |>
    dplyr::filter(Metric == "RMSE") |>
    dplyr::group_by(model) |>
    dplyr::summarise(med = stats::median(Value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(med) |>
    dplyr::pull(model)
}

make_block_metric_box <- function(cv_tbl, which_time = "T14", tail_cap = 0.995) {
  bm <- block_metrics_long(cv_tbl) |>
    dplyr::filter(is.finite(Value))
  if (!is.null(tail_cap)) {
    ymax <- stats::quantile(bm$Value, tail_cap, na.rm = TRUE)
  }
  lev <- order_models_by_median_rmse(cv_tbl)
  bm$model <- factor(bm$model, levels = lev)
  
  ggplot2::ggplot(bm, ggplot2::aes(model, Value)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.35, width = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                          fill = "white", colour = "black", stroke = 0.5) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Block-wise errors (LBO-CV)", which_time),
                  subtitle = "Box = IQR · line = median · ◆ = mean",
                  x = "Model", y = "Error") +
    ggplot2::facet_wrap(~ Metric, scales = "free_y")
}

make_abs_error_box <- function(cv_tbl, which_time = "T14", tail_cap = 0.995) {
  df <- cv_tbl |>
    dplyr::mutate(abs_err = abs(pred - obs)) |>
    dplyr::filter(is.finite(abs_err))
  ymax <- if (!is.null(tail_cap)) stats::quantile(df$abs_err, tail_cap, na.rm = TRUE) else max(df$abs_err, na.rm = TRUE)
  lev <- df |>
    dplyr::group_by(model) |>
    dplyr::summarise(med = stats::median(abs_err, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(med) |>
    dplyr::pull(model)
  df$model <- factor(df$model, levels = lev)
  
  ggplot2::ggplot(df, ggplot2::aes(model, abs_err)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.3, width = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                          fill = "white", colour = "black", stroke = 0.5) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Absolute errors per station (LBO-CV)", which_time),
                  subtitle = "Box = IQR · line = median · ◆ = mean",
                  x = "Model", y = "|pred − obs|")
}


make_obs_pred_scatter <- function(cv_tbl, which_time = "T14") {
  lab <- .make_labeller(cv_tbl)
  ggplot(cv_tbl, aes(obs, pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.7, shape = 16) +
    coord_equal() + theme_minimal() +
    labs(title = sprintf("%s — Observed vs Predicted (LBO-CV)", which_time), x = "Observed", y = "Predicted") +
    facet_wrap(~ model, ncol = 3, labeller = ggplot2::as_labeller(lab))
}

make_residual_density <- function(cv_tbl, which_time = "T14") {
  cv_tbl |> dplyr::mutate(resid = pred - obs) |> ggplot2::ggplot(ggplot2::aes(resid, fill = model)) +
    ggplot2::geom_density(alpha = 0.4) + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Residual density", which_time), x = "Residual (°C)", y = "Density")
}

# Prediction maps & error panels ---------------------------------------
.make_labeller <- function(cv_tbl) {
  m <- cv_tbl |>
    dplyr::group_by(model) |>
    dplyr::summarise(RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)), MAE  = mean(abs(obs - pred), na.rm = TRUE), .groups = "drop")
  setNames(sprintf("%s  (RMSE=%.2f · MAE=%.2f)", m$model, m$RMSE, m$MAE), m$model)
}

# Helfer: Raster schön als ggplot
# --- Helfer: Raster hübsch als ggplot -----------------------------------
.plot_raster_gg <- function(r, title = "", palette = temp_palette, q = c(0.02,0.98), lims = NULL) {
  stopifnot(terra::nlyr(r) == 1)
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  nm <- names(df)[3]
  if (is.null(lims)) {
    vv <- terra::values(r, na.rm = TRUE)
    lims <- stats::quantile(vv, probs = q, na.rm = TRUE, names = FALSE)
  }
  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data[[nm]])) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = palette, limits = lims, oob = scales::squish) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "°C") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold"))
}

build_panels_with_errors <- function(
    maps, truth_raster, cv_tbl, stn_sf, which_time,
    temp_palette = temp_palette, stretch_q = c(0.02,0.98),
    layout = c("horizontal","vertical")
) {
  layout <- match.arg(layout)
  
  preds <- .get_preds_from_maps(maps)
  
  # Einheitliche Farbskala (Truth + alle Rastervorhersagen)
  all_vals <- c(terra::values(truth_raster, na.rm = TRUE))
  for (p in preds) if (inherits(p, "SpatRaster")) all_vals <- c(all_vals, terra::values(p, na.rm = TRUE))
  lims <- stats::quantile(all_vals, probs = stretch_q, na.rm = TRUE, names = FALSE)
  
  # Wahrheit
  p_truth <- .plot_raster_gg(truth_raster, title = paste0(which_time, " — truth"),
                             palette = temp_palette, q = stretch_q, lims = lims)
  
  # Vorhersagen (Raster -> ggplot; ggplot bleibt ggplot)
  p_pred_list <- Map(function(obj, nm) {
    if (inherits(obj, "SpatRaster")) .plot_raster_gg(obj, title = nm, palette = temp_palette, q = stretch_q, lims = lims)
    else if (inherits(obj, "ggplot")) obj + ggplot2::labs(title = nm)
    else stop("Vorhersage-Typ nicht unterstützt: ", class(obj)[1])
  }, preds, names(preds) %||% paste0("model_", seq_along(preds)))
  
  p_preds <- patchwork::wrap_plots(p_pred_list, ncol = if (layout == "vertical") 1 else 3)
  
  # Fehler/Diagnose
  p_box_rmse <- make_block_metric_box(cv_tbl, which_time = which_time, tail_cap = 0.995)
  p_box_ae   <- make_abs_error_box  (cv_tbl, which_time = which_time, tail_cap = 0.995)
  p_scatter  <- make_obs_pred_scatter(cv_tbl, which_time = which_time)
  p_dens     <- make_residual_density(cv_tbl, which_time = which_time)
  p_errs     <- (p_box_rmse | p_box_ae) / (p_scatter | p_dens)
  
  if (layout == "vertical") (p_truth / p_preds / p_errs) else ((p_truth | p_preds) / p_errs)
}
# --- Robust: Vorhersagen aus beliebigen 'maps'-Formen herausziehen ------
.get_preds_from_maps <- function(maps) {
  # 1) SpatRaster direkt
  if (inherits(maps, "SpatRaster")) {
    ul <- terra::unstack(maps)
    names(ul) <- names(maps)
    return(ul)
  }
  # 2) Liste mit typischen Feldern
  if (is.list(maps)) {
    if (!is.null(maps$preds))          return(maps$preds)
    if (!is.null(maps$pred_rasters))   return(maps$pred_rasters)
    if (!is.null(maps$pred_stack) && inherits(maps$pred_stack, "SpatRaster")) {
      ul <- terra::unstack(maps$pred_stack); names(ul) <- names(maps$pred_stack); return(ul)
    }
    if (!is.null(maps$stack) && inherits(maps$stack, "SpatRaster")) {
      ul <- terra::unstack(maps$stack); names(ul) <- names(maps$stack); return(ul)
    }
    if (!is.null(maps$maps) && inherits(maps$maps, "SpatRaster")) {
      ul <- terra::unstack(maps$maps); names(ul) <- names(maps$maps); return(ul)
    }
    # 3) Liste, die bereits einzelne SpatRaster oder ggplots enthält
    cand <- maps[ vapply(maps, function(x) inherits(x, "SpatRaster") || inherits(x, "ggplot"), logical(1)) ]
    if (length(cand) > 0) return(cand)
  }
  stop("build_panels_with_errors(): In 'maps' keine Vorhersagen gefunden.")
}

# --- Kartenplot mit optionalen Achsenticks/labels ----------------------
.plot_map_axes <- function(r, title, cols, lims, q = c(0.02,0.98),
                           base_size = 14, tick_n = 5,
                           show_axis_labels = FALSE, show_axis_ticks = TRUE) {
  stopifnot(terra::nlyr(r) == 1)
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  nm <- names(df)[3]

  if (is.null(lims) || !all(is.finite(lims)) || lims[1] >= lims[2]) {
    vv <- terra::values(r, na.rm = TRUE)
    lims <- stats::quantile(vv, probs = q, na.rm = TRUE, names = FALSE)
    if (!all(is.finite(lims)) || lims[1] == lims[2]) lims <- range(vv, na.rm = TRUE) + c(-1e-6, 1e-6)
  }
  if (is.function(cols)) cols <- cols(256)
  if (!is.atomic(cols) || length(cols) < 2) cols <- grDevices::hcl.colors(256, "YlOrRd")

  ggplot2::ggplot(df, ggplot2::aes(x, y, fill = .data[[nm]])) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0,0), breaks = scales::breaks_pretty(n = tick_n)) +
    ggplot2::scale_y_continuous(expand = c(0,0), breaks = scales::breaks_pretty(n = tick_n)) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = lims, oob = scales::squish) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "°C") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title   = ggplot2::element_blank(),
      axis.text    = if (show_axis_labels) ggplot2::element_text(size = base_size - 3) else ggplot2::element_blank(),
      axis.ticks   = if (show_axis_ticks)  ggplot2::element_line(linewidth = 0.25) else ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey40", linewidth = .4)
    )
}
build_panels_truth_preds_errors_paged <- function(
  maps, truth_raster, cv_tbl, which_time,
  models_per_page = 4,
  model_order      = NULL,
  temp_pal         = temp_palette,     # Vektor ODER Funktion -> wird zu Vektor
  stretch_q        = c(0.02, 0.98),
  errors_height    = 1.2,
  scatter_next_to_truth = TRUE,        # Scatter rechts von Truth?
  top_widths       = c(1.1, 0.9),      # Breitenverhältnis Truth | Scatter
  show_second_legend = FALSE           # zweite °C-Legende bei den Preds unterdrücken
) {
  stopifnot(length(stretch_q) == 2)
  if (is.function(temp_pal)) temp_pal <- temp_pal(256)
  stopifnot(is.atomic(temp_pal), length(temp_pal) >= 2)

  # Vorhersagen einsammeln / Reihenfolge
  preds_raw  <- .get_preds_from_maps(maps)
  pred_names <- names(preds_raw) %||% paste0("model_", seq_along(preds_raw))
  if (!is.null(model_order)) {
    keep <- intersect(model_order, pred_names)
    if (!length(keep)) stop("model_order enthält keine gültigen Modellnamen.")
    preds_raw  <- preds_raw[keep]
    pred_names <- keep
  }

  # Gemeinsame Farbskala
  all_vals <- c(terra::values(truth_raster, na.rm = TRUE))
  for (p in preds_raw) if (inherits(p, "SpatRaster")) all_vals <- c(all_vals, terra::values(p, na.rm = TRUE))
  lims <- stats::quantile(all_vals, probs = stretch_q, na.rm = TRUE, names = FALSE)
  if (all(is.finite(lims)) && lims[1] == lims[2]) {
    eps <- .Machine$double.eps * max(1, abs(lims[1])); lims <- lims + c(-eps, eps)
  }

  # Helfer: Raster -> ggplot
  make_tile <- function(obj, title_txt, show_legend = TRUE) {
    if (inherits(obj, "SpatRaster")) {
      g <- .plot_raster_gg(obj, title = title_txt, palette = temp_pal, q = stretch_q, lims = lims)
      if (!show_legend) g <- g + ggplot2::theme(legend.position = "none")
      g
    } else if (inherits(obj, "ggplot")) {
      obj + ggplot2::labs(title = title_txt)
    } else stop("Nicht unterstützter Prediction-Typ: ", class(obj)[1])
  }

  # Truth (+ optional Scatter daneben)
  p_truth   <- .plot_raster_gg(truth_raster, title = paste0(which_time, " — truth"),
                               palette = temp_pal, q = stretch_q, lims = lims)
  p_scatter <- make_obs_pred_scatter(cv_tbl, which_time = which_time)

  # Prediction-Kacheln bauen (nur erste mit °C-Legende falls gewünscht)
  pred_tiles <- lapply(seq_along(preds_raw), function(i) {
    show_leg <- if (isTRUE(show_second_legend)) TRUE else (i == 1L)
    make_tile(preds_raw[[i]], pred_names[i], show_legend = show_leg)
  })

  # Paginierung
  n <- length(pred_tiles)
  idx_split <- split(seq_len(n), ceiling(seq_len(n) / models_per_page))

  pages <- lapply(idx_split, function(idx) {
    preds_row <- patchwork::wrap_plots(pred_tiles[idx], nrow = 1, ncol = length(idx))

    top_row <- if (isTRUE(scatter_next_to_truth)) {
      (p_truth | (p_scatter + ggplot2::theme(legend.position = "none"))) +
        patchwork::plot_layout(widths = top_widths)
    } else {
      p_truth
    }

    # Errors unten: wenn Scatter schon oben, unten nur Dichte
    p_box_rmse <- make_block_metric_box(cv_tbl, which_time = which_time, tail_cap = 0.995)
    p_box_ae   <- make_abs_error_box  (cv_tbl, which_time = which_time, tail_cap = 0.995)
    p_dens     <- make_residual_density(cv_tbl, which_time = which_time)
    p_errors   <- (p_box_rmse | p_box_ae) / p_dens

    (top_row / preds_row / p_errors) +
      patchwork::plot_layout(heights = c(1, 1, errors_height), guides = "collect") &
      ggplot2::theme(legend.position = "right")
  })

  pages
}




# ------------------------------------------------------------------------------
# run_for_time(stn_sf, truth_r, label, scen_local, block_m, models)
# Purpose:
#   Convenience wrapper to run leave-block-out CV, build prediction maps,
#   and compose the "truth/prediction/error with CV residuals" panel.
# Inputs:
#   - stn_sf: station sf for one time slice (T14 or T05).
#   - truth_r: corresponding truth raster (R14 or R05).
#   - label: "T14" or "T05" for plot titles and selections.
#   - scen_local: scenario object with rasters.
#   - block_m: block size in meters for spatial CV.
#   - models: vector of model names to evaluate.
# Returns:
#   list(res = cv result, maps = predict_maps result, panel = patchwork plot).
# ------------------------------------------------------------------------------
run_for_time <- function(stn_sf, truth_r, label,
                         scen_local = scen,
                         block_m = block_size,
                         models = models_use,
                         layout = c("horizontal","vertical")) {
  layout <- match.arg(layout)
  res   <- run_lbo_cv(stn_sf, scen_local$E, block_size = block_m, models = models)
  maps <- predict_maps(stn_sf, truth_r, which_time = label,
                     scen = scen_local, models = models,
                     lc_levels = scen_local$lc_levels)
  # panel <- build_panels_with_errors(maps, truth_r, res$cv, stn_sf, label,
  #                                   temp_palette = temp_palette, stretch_q = stretch_q,
  #                                   layout = layout)
  # panel <- build_panels_truth_preds_errors_paged(maps =maps, truth_raster = truth_r, cv_tbl = res$cv,  which_time = label)
  panel <- build_panels_truth_preds_errors_paged(
  maps, truth_r, res$cv, label,
  models_per_page = 7,
  scatter_next_to_truth = TRUE
)
  list(res = res, maps = maps, panel = panel)
}

# --- Factor + Level helpers ---------------------------------------------------
.safe_levels <- function(x) unique(na.omit(as.character(x)))

# Train/Grid sauber auf gemeinsame Levels bringen (+ Fallback-Bucket)
.align_factor_pair <- function(train_x, grid_x, fallback = "__OTHER__") {
  tl <- .safe_levels(train_x)
  if (length(tl) == 0L) {
    return(list(use = FALSE, train = NULL, grid = NULL, levels = character()))
  }
  lev <- unique(c(tl, fallback))
  trc <- as.character(train_x); trc[is.na(trc)] <- fallback; trc[!(trc %in% lev)] <- fallback
  gdc <- as.character(grid_x);  gdc[is.na(gdc)]  <- fallback; gdc[!(gdc %in% lev)] <- fallback
  list(use = TRUE, train = factor(trc, levels = lev), grid = factor(gdc, levels = lev), levels = lev)
}

# Optional: numerische LC-Codes auf Labels mappen
.map_lc <- function(codes, lc_levels = NULL) {
  if (!is.null(lc_levels)) {
    ok <- !is.na(codes) & codes >= 1 & codes <= length(lc_levels)
    out <- rep(NA_character_, length(codes)); out[ok] <- lc_levels[codes[ok]]
    out
  } else as.character(codes)
}

# --- Kernel + Extract guards --------------------------------------------------
.mean_kernel_for_R <- function(r, R_m) {
  px <- mean(terra::res(r))
  half <- max(1L, ceiling(R_m / px))     # >= 1 px
  k <- 2L * half + 1L                    # odd size
  W <- matrix(1, nrow = k, ncol = k)
  W / sum(W)
}

smooth_mean_R <- function(r, R_m) {
  W <- .mean_kernel_for_R(r, R_m)
  terra::focal(r, w = W, fun = "mean", na.policy = "omit", pad = TRUE, normalize = FALSE)
}

.extract_to_pts <- function(r, pts_sf) {
  out <- try(terra::extract(r, terra::vect(pts_sf), ID = FALSE)[,1], silent = TRUE)
  if (inherits(out, "try-error") || length(out) == 0L) rep(NA_real_, nrow(pts_sf)) else out
}

# --- GAM guards ---------------------------------------------------------------
.k_for_xy <- function(n, n_xy) max(3, min(60, n_xy - 1L, floor(n * 0.8)))
.k_cap    <- function(x, kmax = 15) {
  ux <- dplyr::n_distinct(x[is.finite(x)])
  max(4, min(kmax, ux - 1L))
}

safe_gam_formula <- function(d, include_lc = FALSE) {
  stopifnot(all(c("temp","x","y") %in% names(d)))
  d <- d[stats::complete.cases(d[, c("temp","x","y")]), , drop = FALSE]
  n    <- nrow(d)
  n_xy <- dplyr::n_distinct(paste0(round(d$x,3), "_", round(d$y,3)))

  base <- if (n_xy >= 4) sprintf("temp ~ s(x,y,bs='tp',k=%d)", .k_for_xy(n, n_xy)) else "temp ~ x + y"

  add <- character(0)
  if ("z_surf" %in% names(d) && dplyr::n_distinct(d$z_surf) > 3)
    add <- c(add, sprintf("s(z_surf,bs='tp',k=%d)", .k_cap(d$z_surf, 20)))
  if ("slp" %in% names(d) && dplyr::n_distinct(d$slp) > 3)
    add <- c(add, sprintf("s(slp,bs='tp',k=%d)", .k_cap(d$slp, 12)))
  if ("cosi" %in% names(d) && dplyr::n_distinct(d$cosi) > 3)
    add <- c(add, sprintf("s(cosi,bs='tp',k=%d)", .k_cap(d$cosi, 12)))

  if (include_lc && "lc" %in% names(d)) {
    d$lc <- droplevels(factor(d$lc))
    if (nlevels(d$lc) >= 2) add <- c(add, "lc")
  }

  as.formula(paste(base, paste(add, collapse = " + "),
                   sep = if (length(add)) " + " else ""))
}


# --- Variogram/Scale utilities -----------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# plot_variogram_with_scales(vg, L50, L95, sill, title)
# Purpose:
#   Quick ggplot helper to display the empirical variogram with dotted sill and
#   dashed vertical markers at L50 and L95 for interpretation.
# Returns: a ggplot object.
# ------------------------------------------------------------------------------
plot_variogram_with_scales <- function(vg, L50, L95, sill, title = "Empirical variogram") {
  df <- as.data.frame(vg)
  ggplot2::ggplot(df, ggplot2::aes(dist, gamma)) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = sill, linetype = "dotted", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = L50, colour = "#2b8cbe", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = L95, colour = "#de2d26", linetype = "dashed") +
    ggplot2::annotate("text", x = L50, y = 0, vjust = -0.5,
                      label = sprintf("L50 = %.0f m", L50)) +
    ggplot2::annotate("text", x = L95, y = 0, vjust = -0.5,
                      label = sprintf("L95 = %.0f m", L95), colour = "#de2d26") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Distance (m)", y = "Semivariance")
}

# --- DEM smoothing + sun geometry --------------------------------------------
# ------------------------------------------------------------------------------
# gaussian_focal(r, radius_m, sigma_m = NULL)
# Purpose:
#   Build a separable, normalized 2D Gaussian kernel (in pixels) from a target
#   radius in meters, then return the kernel matrix to be used with terra::focal.
# Inputs:
#   - r: reference raster to read pixel size from.
#   - radius_m: target smoothing radius (meters).
#   - sigma_m: optional sigma; by default half the radius.
# Returns:
#   A normalized kernel matrix suitable for terra::focal() smoothing.
# ------------------------------------------------------------------------------
gaussian_focal <- function(r, radius_m, sigma_m = NULL) {
  resx <- terra::res(r)[1]
  if (is.null(sigma_m)) sigma_m <- radius_m / 2
  rad_px   <- max(1L, round(radius_m / resx))
  sigma_px <- max(0.5, sigma_m / resx)
  xs <- -rad_px:rad_px
  k1 <- exp(-0.5 * (xs / sigma_px)^2); k1 <- k1 / sum(k1)
  K  <- outer(k1, k1); K / sum(K)
}

# ------------------------------------------------------------------------------
# smooth_dem_and_derive(E, alt, az, radius_m)
# Purpose:
#   Smooth the DEM at a given metric radius and recompute slope and
#   cosine-of-incidence for a specified sun position (alt/az).
# Returns:
#   list(Es = smoothed DEM, slp = slope, cosi = cosine-of-incidence)
# ------------------------------------------------------------------------------
smooth_dem_and_derive <- function(E, alt, az, radius_m) {
  K   <- gaussian_focal(E, radius_m)
  Es  <- terra::focal(E, w = K, fun = mean, na.policy = "omit", pad = TRUE)
  slp <- terra::terrain(Es, v = "slope",  unit = "radians")
  asp <- terra::terrain(Es, v = "aspect", unit = "radians")
  zen <- (pi/2 - alt)
  ci  <- cos(slp)*cos(zen) + sin(slp)*sin(zen)*cos(az - asp)
  ci  <- terra::ifel(ci < 0, 0, ci)
  list(Es = Es, slp = slp, cosi = ci)
}

# --- helpers to cap k by available info --------------------------------
.k_for_xy <- function(n, n_xy) max(3, min(60, n_xy - 1L, floor(n * 0.8)))
.kcap_unique <- function(x, kmax) {
  ux <- unique(x[is.finite(x)])
  nu <- length(ux)
  if (nu <= 3) return(0L)                # treat as constant/near-constant
  max(4L, min(kmax, nu - 1L))
}

# --- CV of GAM with R-smoothed predictors (robust k) -------------------
# ------------------------------------------------------------------------------
# cv_gam_with_R(stn_sf, E, alt, az, R, block_size_m)
# Purpose:
#   Leave-block-out CV of a GAM whose predictors are computed from a DEM
#   smoothed at radius R (meters). This aligns the drift scale to the process
#   scale before fitting, then evaluates predictive skill via blocked holdouts.
# Notes:
#   - Contains guards for low sample size and dynamic k to avoid mgcv errors.
# Returns:
#   list(cv = per-point CV table, RMSE = numeric)
# ------------------------------------------------------------------------------
cv_gam_with_R <- function(stn_sf, E, alt = NULL, az = NULL, R, block_size_m = NULL) {
  
  # ---- 0) Blockgröße sauber auflösen (tuning-fähig)
  bs <- suppressWarnings(as.numeric(block_size_m)[1])             # bevorzugt: Tuning
  if (!is.finite(bs) || bs <= 0) {
    bs <- suppressWarnings(as.numeric(get0("block_size",          # Fallback: global
                                           ifnotfound = NA_real_)))
  }
  if (!is.finite(bs) || bs <= 0)
    stop("cv_gam_with_R(): keine gültige Blockgröße gefunden (Tuning oder global).")
  
  # ---- 1) R-gesmoothete Raster bauen (wie bei dir)
  zR   <- smooth_mean_R(E, R)
  slpR <- terra::terrain(zR, v = "slope",  unit = "radians")
  aspR <- terra::terrain(zR, v = "aspect", unit = "radians")
  cosiR <- if (!is.null(alt) && !is.null(az)) {
    ci <- cos(slpR)*cos(pi/2 - alt) + sin(slpR)*sin(pi/2 - alt)*cos(az - aspR)
    terra::ifel(ci < 0, 0, ci)
  } else NULL
  
  # ---- 2) Werte an Stationen extrahieren (wie bei dir)
  if (!all(c("x","y") %in% names(stn_sf))) {
    xy <- sf::st_coordinates(stn_sf); stn_sf$x <- xy[,1]; stn_sf$y <- xy[,2]
  }
  fill_med <- function(v) { m <- stats::median(v[is.finite(v)], na.rm = TRUE); ifelse(is.finite(v), v, m) }
  stn_sf$z_surf_R <- fill_med(.extract_to_pts(zR,   stn_sf))
  stn_sf$slp_R    <- fill_med(.extract_to_pts(slpR, stn_sf))
  stn_sf$cosi_R   <- if (is.null(cosiR)) rep(NA_real_, nrow(stn_sf)) else fill_med(.extract_to_pts(cosiR, stn_sf))
  
  # ---- 3) Blöcke bauen und zuordnen (einheitlich mit bs)
  bb_poly <- sf::st_as_sfc(sf::st_bbox(stn_sf), crs = sf::st_crs(stn_sf))
  blocks  <- sf::st_make_grid(bb_poly, cellsize = c(bs, bs), what = "polygons")
  blocks  <- sf::st_sf(block_id = seq_along(blocks), geometry = blocks)
  
  stn_blk <- sf::st_join(stn_sf, blocks, join = sf::st_intersects, left = TRUE)
  if (anyNA(stn_blk$block_id)) {
    i <- is.na(stn_blk$block_id)
    stn_blk$block_id[i] <- blocks$block_id[sf::st_nearest_feature(stn_blk[i,], blocks)]
  }
  
  if (!all(c("x","y") %in% names(stn_blk))) {
    xy <- sf::st_coordinates(stn_blk); stn_blk$x <- xy[,1]; stn_blk$y <- xy[,2]
  }
  
  # ---- 4) CV-Schleife (dein Code unverändert weiter)
  bids  <- sort(unique(stn_blk$block_id))
  preds <- vector("list", length(bids)); j <- 0L
  for (b in bids) {
    te <- stn_blk[stn_blk$block_id == b, ]
    tr <- stn_blk[stn_blk$block_id != b, ]
    
    dtr <- sf::st_drop_geometry(tr)
    need <- c("temp","x","y","z_surf_R","slp_R","cosi_R")
    dtr  <- dtr[stats::complete.cases(dtr[, intersect(need, names(dtr)), drop = FALSE]), , drop = FALSE]
    if (nrow(dtr) < 10) next
    
    n_xy <- dplyr::n_distinct(paste0(round(dtr$x,3), "_", round(dtr$y,3)))
    k_xy <- .k_for_xy(nrow(dtr), n_xy)
    k_z  <- .kcap_unique(dtr$z_surf_R, 20)
    k_sl <- .kcap_unique(dtr$slp_R,    12)
    k_ci <- .kcap_unique(dtr$cosi_R,   12)
    
    terms <- c()
    terms <- c(terms, if (n_xy >= 4) sprintf("s(x,y,bs='tp',k=%d)", k_xy) else "x + y")
    terms <- c(terms, if (k_z  >= 4) sprintf("s(z_surf_R,bs='tp',k=%d)", k_z)  else "z_surf_R")
    if (length(unique(dtr$slp_R[is.finite(dtr$slp_R)])) > 1)
      terms <- c(terms, if (k_sl >= 4) sprintf("s(slp_R,bs='tp',k=%d)", k_sl) else "slp_R")
    if (any(is.finite(dtr$cosi_R)) && length(unique(dtr$cosi_R[is.finite(dtr$cosi_R)])) > 1)
      terms <- c(terms, if (k_ci >= 4) sprintf("s(cosi_R,bs='tp',k=%d)", k_ci) else "cosi_R")
    
    form <- as.formula(paste("temp ~", paste(terms, collapse = " + ")))
    gm <- mgcv::gam(form, data = dtr, method = "REML", select = TRUE)
    
    dte <- sf::st_drop_geometry(te)
    ph  <- try(stats::predict(gm, newdata = dte, type = "response"), silent = TRUE)
    if (inherits(ph, "try-error")) ph <- rep(NA_real_, nrow(dte))
    
    j <- j + 1L
    preds[[j]] <- tibble::tibble(id = te$id, obs = te$temp, pred = as.numeric(ph), block_id = b)
  }
  
  preds <- preds[seq_len(j)]
  if (!length(preds)) {
    return(list(cv = tibble::tibble(id = integer(), obs = numeric(), pred = numeric(), block_id = integer()),
                RMSE = NA_real_))
  }
  out  <- dplyr::bind_rows(preds)
  rmse <- sqrt(mean((out$pred - out$obs)^2, na.rm = TRUE))
  list(cv = out, RMSE = rmse)
}

# --- U-curve tuning -----------------------------------------------------------
# ------------------------------------------------------------------------------
# tune_Rstar_ucurve(stn_sf, E, alt, az, L50, L95, block_fallback, n_grid, extra)
# Purpose:
#   Scan candidate R values (around the L50–L95 interval) and pick R* that
#   minimises blocked-CV RMSE. Returns the scan table and chosen R*.
# Returns:
#   list(grid = data.frame(R, RMSE), R_star, block_m)
# ------------------------------------------------------------------------------
tune_Rstar_ucurve <- function(stn_sf, E, alt, az, L50, L95, block_fallback = 120, n_grid = 6, extra = c(0.8, 1.2)) {
  L50 <- as.numeric(L50); L95 <- as.numeric(L95)
  if (!is.finite(L50) || !is.finite(L95) || L95 <= L50) {
    e <- terra::ext(E)
    dom_diag <- sqrt((terra::xmax(e)-terra::xmin(e))^2 + (terra::ymax(e)-terra::ymin(e))^2)
    L50 <- dom_diag/10; L95 <- dom_diag/4
  }
  block_m <- max(block_fallback, round(L50))
  R_min <- max(10, round(L50*extra[1])); R_max <- round(L95*extra[2])
  R_grid <- unique(round(seq(R_min, R_max, length.out = n_grid)))
  df <- do.call(rbind, lapply(R_grid, function(R) { z <- cv_gam_with_R(stn_sf, E, alt, az, R = R, block_size_m = block_m); c(R = R, RMSE = z$RMSE) })) |> as.data.frame()
  R_star <- df$R[which.min(df$RMSE)]
  list(grid = df, R_star = as.numeric(R_star), block_m = block_m)
}

# ------------------------------------------------------------------------------
# plot_ucurve(df, R_star, title)
# Purpose:
#   Visual helper to display the U-curve of RMSE vs. drift radius R with a
#   dashed marker at the selected R*.
# Returns: a ggplot object.
# ------------------------------------------------------------------------------
plot_ucurve <- function(df, R_star, title = "U-curve: tune R") {
  ggplot2::ggplot(df, ggplot2::aes(R, RMSE)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = R_star, linetype = "dashed", colour = "#de2d26") +
    ggplot2::theme_minimal() + ggplot2::labs(title = title, x = "Drift radius R (m)", y = "RMSE (block-CV)")
}

# --- Factor alignment (robust) -----------------------------------------------
.align_factor_to_model <- function(x, lev_model) {
  xs <- as.character(x)
  if (length(lev_model) == 0L) return(factor(rep(NA_character_, length(xs))))
  y <- factor(xs, levels = lev_model)
  if (anyNA(y)) {
    xs[is.na(y)] <- lev_model[1]
    y <- factor(xs, levels = lev_model)
  }
  y
}

# ---------- SF-only learners ----------
pred_Voronoi <- function(train_sf, test_sf) {
  idx <- sf::st_nearest_feature(test_sf, train_sf)
  as.numeric(train_sf$temp)[idx]
}

pred_IDW <- function(train_sf, test_sf, idp = 2) {
  pr <- suppressWarnings(gstat::idw(temp ~ 1, locations = train_sf, newdata = test_sf, idp = idp))
  as.numeric(pr$var1.pred)
}

.default_vgm <- function(values, model = "Exp", range = 100) {
  psill <- stats::var(values, na.rm = TRUE); nug <- 0.1 * psill
  gstat::vgm(psill = psill, model = model, range = range, nugget = nug)
}

pred_OK <- function(train_sf, test_sf) {
  vg      <- suppressWarnings(gstat::variogram(temp ~ 1, data = train_sf))
  vgm_fit <- try(suppressWarnings(gstat::fit.variogram(vg, gstat::vgm("Exp"))), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) vgm_fit <- .default_vgm(train_sf$temp)
  pr <- suppressWarnings(gstat::krige(temp ~ 1, locations = train_sf, newdata = test_sf, model = vgm_fit))
  as.numeric(pr$var1.pred)
}

.align_factor_to_model <- function(x, lev_model) {
  y <- factor(as.character(x), levels = lev_model)
  if (anyNA(y)) y[is.na(y)] <- lev_model[1]
  y
}
.fill_num_na_vec <- function(x, ref) {
  m <- stats::median(ref[is.finite(ref)], na.rm = TRUE)
  x[!is.finite(x)] <- m
  x
}

# --- sf-only KED (schluckt extra Args wie E=E) ----------------------
pred_KED <- function(train_sf, test_sf, ...) {
  stopifnot(inherits(train_sf, "sf"), inherits(test_sf, "sf"))
  need <- c("z_surf","slp","cosi")
  miss <- setdiff(need, names(train_sf))
  if (length(miss)) stop("pred_KED(): fehlende Drift-Spalten im Training: ",
                         paste(miss, collapse = ", "))
  
  # Optional LC als Faktor angleichen
  use_lc <- "lc" %in% names(train_sf) && "lc" %in% names(test_sf)
  tr <- train_sf
  te <- test_sf
  if (use_lc) {
    tr$lc <- droplevels(factor(tr$lc))
    te$lc <- factor(as.character(te$lc), levels = levels(tr$lc))
    te$lc[is.na(te$lc)] <- levels(tr$lc)[1]
  }
  
  # fehlende numerische Drifts im TEST mit Trainingsmedian auffüllen
  for (nm in need) {
    m <- stats::median(tr[[nm]][is.finite(tr[[nm]])], na.rm = TRUE)
    te[[nm]][!is.finite(te[[nm]])] <- m
  }
  
  # nur vollständige Trainingszeilen
  keep_tr <- c("temp", need, if (use_lc) "lc")
  dtr <- sf::st_drop_geometry(tr)[, keep_tr, drop = FALSE]
  ok  <- stats::complete.cases(dtr)
  tr  <- tr[ok, ]
  if (nrow(tr) < 5) return(rep(NA_real_, nrow(te)))
  
  # Formel: lineare Drifts + optional LC
  form <- stats::as.formula(paste("temp ~", paste(c(need, if (use_lc) "lc"), collapse = " + ")))
  
  # Variogramm + robuster Fit
  vg      <- suppressWarnings(gstat::variogram(form, data = tr))
  vgm_fit <- try(suppressWarnings(gstat::fit.variogram(vg, gstat::vgm("Exp"))), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) {
    ps <- stats::var(sf::st_drop_geometry(tr)$temp, na.rm = TRUE)
    vgm_fit <- gstat::vgm(psill = ps, model = "Exp", range = max(vg$dist, na.rm = TRUE)/3, nugget = 0.1*ps)
  }
  
  # Kriging mit externen Drifts (UK/KED)
  pr <- suppressWarnings(gstat::krige(form, locations = tr, newdata = te, model = vgm_fit))
  as.numeric(pr$var1.pred)
}


# --- Leave-Block-Out CV -------------------------------------------------------
# ------------------------------------------------------------------------------
# make_blocks_and_assign(pts_sf, E, block_size)
# Purpose:
#   Build a square grid of spatial blocks and assign each station to a block
#   (nearest if on edge). Used by leave-block-out cross-validation.
# Inputs:
#   - pts_sf: station sf with geometry.
#   - E: reference raster for domain extent/CRS.
#   - block_size: block edge length in meters.
# Returns:
#   list(blocks = sf polygons, pts = station sf with block_id).
# ------------------------------------------------------------------------------
make_blocks_and_assign <- function(pts_sf, E, block_size = 100) {
  bb <- sf::st_as_sfc(sf::st_bbox(c(xmin = terra::xmin(E), ymin = terra::ymin(E), xmax = terra::xmax(E), ymax = terra::ymax(E)), crs = sf::st_crs(pts_sf)))
  gr <- sf::st_make_grid(bb, cellsize = c(block_size, block_size), what = "polygons")
  blocks <- sf::st_sf(block_id = seq_along(gr), geometry = gr)
  pts_blk <- sf::st_join(pts_sf, blocks, join = sf::st_intersects, left = TRUE)
  if (any(is.na(pts_blk$block_id))) {
    nearest <- sf::st_nearest_feature(pts_blk[is.na(pts_blk$block_id), ], blocks)
    pts_blk$block_id[is.na(pts_blk$block_id)] <- blocks$block_id[nearest]
  }
  list(blocks = blocks, pts = pts_blk)
}

# saubere Farbskala für viele Blöcke
.discrete_cols <- function(n) scales::hue_pal()(n)

plot_blocks_grid <- function(blocks, pts_blk, title = "Blocks & stations") {
  # Ziel-CRS = CRS der Daten (UTM32N), Achsen in Metern
  crs_plot <- sf::st_crs(pts_blk)
  bb       <- sf::st_bbox(blocks)
  n_blocks <- dplyr::n_distinct(pts_blk$block_id)
  cols     <- .discrete_cols(max(1, n_blocks))

  ggplot() +
    geom_sf(data = blocks, fill = NA, color = "grey50", linewidth = 0.25) +
    geom_sf(data = pts_blk, aes(color = factor(block_id)), size = 2, alpha = 0.95) +
    scale_color_manual(values = cols, name = "Block") +
    coord_sf(
      crs  = crs_plot,    # <- erzwingt UTM32N als Plot-CRS (Meterachsen)
      datum = NA,         # keine Gradnetz-Beschriftung
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    theme_minimal() +
    labs(title = title, x = "Easting (m)", y = "Northing (m)")
}


# ------------------------------------------------------------------------------
# run_lbo_cv(stn_sf, E, block_size, models)
# Purpose:
#   Perform leave-block-out cross-validation across the requested set of models.
#   Each block is held out in turn; models are trained on the remainder and
#   predictions are collected for the held-out stations.
# Returns:
#   list(cv = long per-point table, metrics = summary table, diag_plot, blocks_plot)
# Notes:
#   - No model settings are changed; this wrapper only orchestrates the CV.
# ------------------------------------------------------------------------------
run_lbo_cv <- function(stn_sf, E, block_size = 100, models = models_use) {
  if (!all(c("x","y") %in% names(stn_sf))) { xy <- sf::st_coordinates(stn_sf); stn_sf$x <- xy[,1]; stn_sf$y <- xy[,2] }
  blk <- make_blocks_and_assign(stn_sf, E, block_size = block_size)
  blocks_sf <- blk$blocks; stn_blk <- blk$pts
  restore <- function(nm) if (!(nm %in% names(stn_blk))) stn_blk[[nm]] <<- stn_sf[[nm]][match(stn_blk$id, stn_sf$id)]
  for (nm in c("temp","z_surf","slp","cosi","lc","x","y")) restore(nm)

  block_ids <- sort(unique(stn_blk$block_id))
  out_list <- vector("list", length(block_ids))
  for (k in seq_along(block_ids)) {
    b <- block_ids[k]
    test_idx  <- which(stn_blk$block_id == b)
    train_idx <- which(stn_blk$block_id != b)
    train_sf <- stn_blk[train_idx, ]; test_sf <- stn_blk[test_idx, ]
    pred_tbl <- lapply(models, function(m) {
      p <- switch(m,
        "Voronoi" = pred_Voronoi(train_sf, test_sf),
        "IDW"     = pred_IDW(train_sf, test_sf),
        "OK"      = pred_OK(train_sf, test_sf),
        "KED"     = pred_KED(train_sf, test_sf, E = E),
        "RF"      = pred_RF(train_sf, test_sf),
        "GAM"     = pred_GAM(train_sf, test_sf),
        stop("Unknown model: ", m)
      )
      tibble::tibble(model = m, id = test_sf$id, obs = test_sf$temp, pred = p, block_id = b)
    })
    out_list[[k]] <- dplyr::bind_rows(pred_tbl)
  }

  cv_tbl <- dplyr::bind_rows(out_list)
  metrics <- cv_tbl %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      n_ok = sum(is.finite(obs) & is.finite(pred)),
      MAE  = {i <- is.finite(obs) & is.finite(pred); if (any(i)) mean(abs(pred[i]-obs[i])) else NA_real_},
      RMSE = {i <- is.finite(obs) & is.finite(pred); if (any(i)) sqrt(mean((pred[i]-obs[i])^2)) else NA_real_},
      Bias = {i <- is.finite(obs) & is.finite(pred); if (any(i)) mean(pred[i]-obs[i]) else NA_real_},
      R2   = safe_r2(obs, pred),
      .groups = "drop"
    ) |>
    dplyr::arrange(RMSE)

  diag_plot <- ggplot(cv_tbl, aes(obs, pred)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_point(alpha=0.7) +
    coord_equal() + theme_minimal() +
    labs(title = sprintf("LBO-CV (block = %dm) — Observed vs Predicted", block_size), x = "Observed", y = "Predicted") +
    facet_wrap(~ model)

  blocks_plot <- plot_blocks_grid(blocks_sf, stn_blk, title = sprintf("Blocks (%.0f m) & stations", block_size))
  list(cv = cv_tbl, metrics = metrics, diag_plot = diag_plot, blocks_plot = blocks_plot)
}

message("Running LBO-CV and building maps for T14 ...")
lbo_cv_14_result <- run_for_time(stn_sf_14, scen$R14, "T14")
message("Running LBO-CV and building maps for T05 ...")
lbo_cv_05_result <- run_for_time(stn_sf_05, scen$R05, "T05")




lbo_cv_14_result$res$blocks_plot 
lbo_cv_14_result$res$diag_plot
lbo_cv_05_result$res$blocks_plot  
lbo_cv_05_result$res$diag_plot
lbo_cv_14_result$panel
lbo_cv_05_result$panel
p_block_box14  <- make_block_metric_box(lbo_cv_14_result$res$cv, "T14")
p_abserr_box14 <- make_abs_error_box(lbo_cv_14_result$res$cv,  "T14")
p_block_box05  <- make_block_metric_box(lbo_cv_05_result$res$cv, "T05")
p_abserr_box05 <- make_abs_error_box(lbo_cv_05_result$res$cv, "T05")

(p_block_box14 | p_abserr_box14) / (p_block_box05 | p_abserr_box05)

knitr::kable(lbo_cv_14_result$res$metrics, digits = 3, caption = "LBO-CV metrics — T14")
knitr::kable(lbo_cv_05_result$res$metrics, digits = 3, caption = "LBO-CV metrics — T05")

make_obs_pred_scatter(lbo_cv_14_result$res$cv, "T14")
make_obs_pred_scatter(lbo_cv_05_result$res$cv, "T05")

make_residual_density(lbo_cv_14_result$res$cv, "T14")
make_residual_density(lbo_cv_05_result$res$cv, "T05")


# --- Scale inference (variogram -> L50/L95) ----------------------------------
# --- SCALE → TUNING → FEATURES @R* → CV → MAPS → PANELS ----------------------
# add smoothed predictors to both station sets
add_drifts_at_R <- function(stn_sf, E, alt, az, R, lc = NULL, lc_levels = NULL) {
  # R-geglättete Prädiktor-Raster bauen (E*, slp*, cosi*)
  dr <- smooth_dem_and_derive(E, alt, az, radius_m = R)
  XY <- sf::st_coordinates(stn_sf)
  
  # Extraktion der bei R geglätteten Prädiktoren
  stn_sf$z_surf <- as.numeric(terra::extract(dr$Es,   XY)[,1])
  stn_sf$slp    <- as.numeric(terra::extract(dr$slp,  XY)[,1])
  stn_sf$cosi   <- as.numeric(terra::extract(dr$cosi, XY)[,1])
  
  # Optional: Landnutzung als Faktor (nicht geglättet, aber konsistent gemappt)
  if (!is.null(lc)) {
    lc_codes <- as.integer(terra::extract(lc, XY)[,1])
    if (!is.null(lc_levels)) {
      lc_codes[is.na(lc_codes)] <- 1L
      lc_codes <- pmax(1L, pmin(lc_codes, length(lc_levels)))
      stn_sf$lc <- factor(lc_levels[lc_codes], levels = lc_levels)
    } else {
      stn_sf$lc <- factor(lc_codes)
    }
  }
  
  stn_sf
}

# ---- SF Variogram Scales (ohne sp::) --------------------------------
compute_Ls_from_points <- function(stn_sf, value_col = "temp",
                                   maxdist = NULL, nlag = 18, smooth_k = 3) {
  stopifnot(inherits(stn_sf, "sf"), value_col %in% names(stn_sf))
  pts <- stn_sf[is.finite(stn_sf[[value_col]]), ]
  if (is.null(maxdist)) {
    bb <- sf::st_bbox(pts)
    dom_diag <- sqrt((bb["xmax"]-bb["xmin"])^2 + (bb["ymax"]-bb["ymin"])^2)
    maxdist <- dom_diag / 2
  }
  form <- stats::as.formula(sprintf("%s ~ 1", value_col))
  vg  <- gstat::variogram(form, data = pts, cutoff = maxdist, width = maxdist/nlag)
  
  if (nrow(vg) >= smooth_k) {
    vg$gamma <- stats::filter(vg$gamma, rep(1/smooth_k, smooth_k), sides = 2)
    vg$gamma[!is.finite(vg$gamma)] <- zoo::na.approx(vg$gamma, na.rm = FALSE)
    vg$gamma <- zoo::na.locf(zoo::na.locf(vg$gamma, fromLast = TRUE))
  }
  sill <- max(vg$gamma, na.rm = TRUE)
  if (!is.finite(sill) || sill <= 0) sill <- stats::median(vg$gamma, na.rm = TRUE)
  
  L_at_q <- function(q) {
    thr <- q * sill
    i   <- which(vg$gamma >= thr)[1]
    if (is.na(i)) return(NA_real_)
    if (i == 1) return(vg$dist[1])
    d0 <- vg$dist[i-1]; d1 <- vg$dist[i]
    g0 <- vg$gamma[i-1]; g1 <- vg$gamma[i]
    if (!is.finite(d0) || !is.finite(d1) || g1 == g0) return(d1)
    d0 + (thr - g0) * (d1 - d0) / (g1 - g0)
  }
  list(vg = vg, sill = sill, L50 = L_at_q(0.5), L95 = L_at_q(0.95), cutoff = maxdist)
}


## 1) Skalen (Variogramm) aus den Stationspunkten
Ls14 <- compute_Ls_from_points(stn_sf_14, value_col = "temp")
Ls05 <- compute_Ls_from_points(stn_sf_05, value_col = "temp")

## (optional) Variogramm-Plots
p_vg14 <- plot_variogram_with_scales(Ls14$vg, Ls14$L50, Ls14$L95, Ls14$sill,
                                     title = "T14 — empirical variogram with L50/L95")
p_vg05 <- plot_variogram_with_scales(Ls05$vg, Ls05$L50, Ls05$L95, Ls05$sill,
                                     title = "T05 — empirical variogram with L50/L95")

## 2) R* via U-Kurve mit GAM@R (geblockte CV; Blockgröße an globalem block_size)
tune14 <- tune_Rstar_ucurve(
  stn_sf = stn_sf_14, E = scen$E, alt = sun14$alt, az = sun14$az,
  L50 = Ls14$L50, L95 = Ls14$L95, block_fallback = block_size, n_grid = 6
)
tune05 <- tune_Rstar_ucurve(
  stn_sf = stn_sf_05, E = scen$E, alt = sun05$alt, az = sun05$az,
  L50 = Ls05$L50, L95 = Ls05$L95, block_fallback = block_size, n_grid = 6
)


## (optional) U-Kurven plotten
print(plot_ucurve(tune14$grid, tune14$R_star, title = "T14 — U-curve"))
print(plot_ucurve(tune05$grid, tune05$R_star, title = "T05 — U-curve"))

message(sprintf("Chosen R* — T14: %d m | blocks ≈ %d m", tune14$R_star, tune14$block_m))
message(sprintf("Chosen R* — T05: %d m | blocks ≈ %d m", tune05$R_star, tune05$block_m))

## 3) Feature-Raster @R* (E*, slp*, cosi*)
fr14 <- smooth_dem_and_derive(scen$E, sun14$alt, sun14$az, radius_m = tune14$R_star)
fr05 <- smooth_dem_and_derive(scen$E, sun05$alt, sun05$az, radius_m = tune05$R_star)

## 4) Stations-Features @R* (inkl. LC, falls vorhanden)
stn14_R <- add_drifts_at_R(stn_sf_14, scen$E, sun14$alt, sun14$az, tune14$R_star,
                           lc = scen$lc, lc_levels = scen$lc_levels)
stn05_R <- add_drifts_at_R(stn_sf_05, scen$E, sun05$alt, sun05$az, tune05$R_star,
                           lc = scen$lc, lc_levels = scen$lc_levels)


## 5) Block-CV gegen E* (nicht gegen rohes E)
bench14 <- run_lbo_cv(stn14_R, E = fr14$Es, block_size = block_size, models = models_use)
bench05 <- run_lbo_cv(stn05_R, E = fr05$Es, block_size = block_size, models = models_use)



## 6) Karten mit feature_rasters = {E*, slp*, cosi*}
maps14_tuned <- predict_maps(
  stn_sf = stn14_R, truth_raster = scen$R14, which_time = "T14",
  scen = scen, models = models_use, lc_levels = scen$lc_levels,
  feature_rasters = list(E = fr14$Es, slp = fr14$slp, cosi = fr14$cosi)
)
maps05_tuned <- predict_maps(
  stn_sf = stn05_R, truth_raster = scen$R05, which_time = "T05",
  scen = scen, models = models_use, lc_levels = scen$lc_levels,
  feature_rasters = list(E = fr05$Es, slp = fr05$slp, cosi = fr05$cosi)
)

## 7) Panels (Truth | Predictions | Error/Residuals) – horizontal, gut lesbar
panel_pages_T14 <- build_panels_truth_preds_errors_paged(
  maps          = maps14_tuned,      # list with $pred_rasters etc.
  truth_raster  = scen$R14,
  cv_tbl        = bench14$cv,
  which_time    = "T14",
  models_per_page     = 7,            # all models on one page
  scatter_next_to_truth = TRUE,
  top_widths           = c(1.1, 0.9), # optional
  show_second_legend   = FALSE        # keep only one °C legend
)
# render the (only) page
panel_pages_T05 <- build_panels_truth_preds_errors_paged(
  maps          = maps05_tuned,      # list with $pred_rasters etc.
  truth_raster  = scen$R05,
  cv_tbl        = bench05$cv,
  which_time    = "T05",
  models_per_page     = 7,            # all models on one page
  scatter_next_to_truth = TRUE,
  top_widths           = c(1.1, 0.9), # optional
  show_second_legend   = FALSE        # keep only one °C legend
)
# render the (only) page
print(panel_pages_T14[[1]])
print(panel_pages_T05[[1]])

print(bench14)
print(bench05)
# --- Error budget from CV residuals -----------------------------------------

simple_error_budget <- function(res_cv, sigma_inst = 0.5, alpha = 0.6) {
  res <- res_cv$cv
  res <- res[is.finite(res$obs) & is.finite(res$pred), , drop = FALSE]
  RMSE <- sqrt(mean((res$pred - res$obs)^2))
  Bias <- mean(res$pred - res$obs)
  VarE <- var(res$pred - res$obs)
  # crude split of process vs. measurement noise
  meas <- sigma_inst^2
  proc <- max(0, VarE - meas)
  micro <- alpha * proc
  meso  <- (1 - alpha) * proc
  tibble::tibble(Component = c("RMSE","Bias","Total var","Instrument var","Microscale var","Mesoscale var"),
                 Value     = c(RMSE, Bias, VarE, meas, micro, meso))
}

# Use GAM@R* CV runs:
cv14 <- cv_gam_with_R(stn_sf_14, scen$E, sun14$alt, sun14$az, R = tune14$R_star, block_size_m = tune14$block_m)
cv05 <- cv_gam_with_R(stn_sf_05, scen$E, sun05$alt, sun05$az, R = tune05$R_star, block_size_m = tune05$block_m)

eb14 <- simple_error_budget(cv14, sigma_inst = 0.5, alpha = 0.6)
eb05 <- simple_error_budget(cv05, sigma_inst = 0.5, alpha = 0.6)

knitr::kable(eb14, digits = 3, caption = "Error budget — T14 (GAM @ R*)")
knitr::kable(eb05, digits = 3, caption = "Error budget — T05 (GAM @ R*)")



# ====================== EXPORT: Plots, Tabellen, Raster ======================
# Ordnerstruktur
out_dir <- "exports"
fig_dir <- file.path(out_dir, "figs")
tab_dir <- file.path(out_dir, "tables")
ras_dir <- file.path(out_dir, "rasters")
dat_dir <- file.path(out_dir, "data")
dir.create(out_dir, showWarnings = FALSE)
for (d in c(fig_dir, tab_dir, ras_dir, dat_dir)) dir.create(d, showWarnings = FALSE)

# Helfer
safe_save_plot <- function(p, file, w = 9, h = 6, dpi = 300) {
  if (!is.null(p) && inherits(p, c("gg","ggplot","patchwork"))) {
    try(ggplot2::ggsave(filename = file, plot = p, width = w, height = h,
                        dpi = dpi, bg = "white"), silent = TRUE)
  }
}
safe_write_csv <- function(x, file) {
  if (!is.null(x)) try(utils::write.csv(x, file, row.names = FALSE), silent = TRUE)
}
safe_save_kable <- function(tab, file_html, self_contained = TRUE) {
  # normalize and create the parent directory if needed
  out_path <- normalizePath(file_html, winslash = "/", mustWork = FALSE)
  out_dir  <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # try saving via kableExtra; on failure, write the raw HTML as a fallback
  tryCatch(
    {
      kableExtra::save_kable(tab, out_path, self_contained = self_contained)
      message("Saved table: ", out_path)
    },
    error = function(e) {
      warning("save_kable failed (", conditionMessage(e), "). Writing raw HTML fallback.")
      cat(as.character(tab), file = out_path)
      message("Fallback table written: ", out_path)
    }
  )

  invisible(out_path)
}




# ====================== EXPORT: Plots, Tabellen, Raster ======================
# Ordnerstruktur
base_dir=here::here()
out_dir <- file.path(base_dir,"block4_5/exports")
fig_dir <- file.path(out_dir, "figs")
tab_dir <- file.path(out_dir, "tables")
ras_dir <- file.path(out_dir, "rasters")
dat_dir <- file.path(out_dir, "data")
dir.create(out_dir, showWarnings = FALSE)
for (d in c(fig_dir, tab_dir, ras_dir, dat_dir)) dir.create(d, showWarnings = FALSE)

# Helfer
safe_save_plot <- function(p, file, w = 9, h = 6, dpi = 300) {
  if (!is.null(p) && inherits(p, c("gg","ggplot","patchwork"))) {
    try(ggplot2::ggsave(filename = file, plot = p, width = w, height = h,
                        dpi = dpi, bg = "white"), silent = TRUE)
  }
}
safe_write_csv <- function(x, file) {
  if (!is.null(x)) try(utils::write.csv(x, file, row.names = FALSE), silent = TRUE)
}
safe_save_kable <- function(df, file_html, caption = NULL) {
  if (!is.null(df) && requireNamespace("kableExtra", quietly = TRUE)) {
    tab <- knitr::kable(df, digits = 3, caption = caption, format = "html")
    kableExtra::save_kable(tab, file_html, self_contained = TRUE)
  }
}

# ---------- Plots sammeln & speichern ----------
# ---- Plots sicher speichern (fixte Referenzen) ----
plots <- list(
  # Panels (tuned)
  "T14_panel_tuned.png"  = if (exists("panel14_tuned")) panel14_tuned else NULL,
  "T05_panel_tuned.png"  = if (exists("panel05_tuned")) panel05_tuned else NULL,
  
  # LBO-CV Übersicht (BLOCKS + DIAG)
  "T14_blocks.png"       = if (exists("lbo_cv_14_result")) lbo_cv_14_result$res$blocks_plot else NULL,
  "T14_diag.png"         = if (exists("lbo_cv_14_result")) lbo_cv_14_result$res$diag_plot   else NULL,
  "T05_blocks.png"       = if (exists("lbo_cv_05_result")) lbo_cv_05_result$res$blocks_plot else NULL,
  "T05_diag.png"         = if (exists("lbo_cv_05_result")) lbo_cv_05_result$res$diag_plot   else NULL,
  
  # Truth/Pred (baseline)
  "T14_truth.png"        = if (exists("lbo_cv_14_result")) lbo_cv_14_result$maps$p_truth else NULL,
  "T14_pred.png"         = if (exists("lbo_cv_14_result")) lbo_cv_14_result$maps$p_pred  else NULL,
  "T05_truth.png"        = if (exists("lbo_cv_05_result")) lbo_cv_05_result$maps$p_truth else NULL,
  "T05_pred.png"         = if (exists("lbo_cv_05_result")) lbo_cv_05_result$maps$p_pred  else NULL,
  
  # Truth/Pred (tuned @ R*)
  "T14_truth_TUNED.png"  = if (exists("maps14_tuned")) maps14_tuned$p_truth else NULL,
  "T14_pred_TUNED.png"   = if (exists("maps14_tuned")) maps14_tuned$p_pred  else NULL,
  "T05_truth_TUNED.png"  = if (exists("maps05_tuned")) maps05_tuned$p_truth else NULL,
  "T05_pred_TUNED.png"   = if (exists("maps05_tuned")) maps05_tuned$p_pred  else NULL,
  
  # Boxplots, Scatter, Density (aus CV)
  "T14_block_box.png"    = if (exists("p_block_box14"))  p_block_box14  else NULL,
  "T14_abserr_box.png"   = if (exists("p_abserr_box14")) p_abserr_box14 else NULL,
  "T05_block_box.png"    = if (exists("p_block_box05"))  p_block_box05  else NULL,
  "T05_abserr_box.png"   = if (exists("p_abserr_box05")) p_abserr_box05 else NULL,
  "T14_obs_pred.png"     = if (exists("p_scatter14")) p_scatter14 else NULL,
  "T05_obs_pred.png"     = if (exists("p_scatter05")) p_scatter05 else NULL,
  "T14_resid_density.png"= if (exists("p_dens14"))    p_dens14    else NULL,
  "T05_resid_density.png"= if (exists("p_dens05"))    p_dens05    else NULL,
  
  # Variogramme & U-Kurven
  "T14_variogram.png"    = if (exists("p_vg14")) p_vg14 else NULL,
  "T05_variogram.png"    = if (exists("p_vg05")) p_vg05 else NULL,
  "T14_ucurve.png"       = if (exists("p_uc14")) p_uc14 else NULL,
  "T05_ucurve.png"       = if (exists("p_uc05")) p_uc05 else NULL,
  
  # Landuse / Terrain-Overview aus Teil 1
  "landcover_vertical.png" = if (exists("p_landcover_vert")) p_landcover_vert else NULL,
  "overview_2x2.png"       = if (exists("p_overview2x2"))    p_overview2x2    else NULL
)

for (nm in names(plots)) safe_save_plot(plots[[nm]], file.path(fig_dir, nm), w = 9, h = 6, dpi = 300)

# ---------- Tabellen & Daten exportieren ----------
# CV-Punktvorhersagen
# if (exists("lbo_cv_14_result")) safe_write_csv(lbo_cv_14_result$res$cv, file.path(dat_dir, "cv_points_T14.csv"))
# if (exists("lbo_cv_05_result")) safe_write_csv(lbo_cv_05_result$res$cv, file.path(dat_dir, "cv_points_T05.csv"))

# Grid-Vorhersagen (Truth/Pred) – “pred_df” der Map-Objekte
# if (exists("lbo_cv_14_result")) safe_write_csv(lbo_cv_14_result$maps$pred_df, file.path(dat_dir, "grid_pred_T14.csv"))
# if (exists("lbo_cv_05_result")) safe_write_csv(lbo_cv_05_result$maps$pred_df, file.path(dat_dir, "grid_pred_T05.csv"))
# if (exists("maps14_tuned")) safe_write_csv(maps14_tuned$pred_df, file.path(dat_dir, "grid_pred_T14_TUNED.csv"))
# if (exists("maps05_tuned")) safe_write_csv(maps05_tuned$pred_df, file.path(dat_dir, "grid_pred_T05_TUNED.csv"))

# Metriken
if (exists("lbo_cv_14_result")) {
  safe_write_csv(lbo_cv_14_result$res$metrics, file.path(tab_dir, "metrics_T14_base.csv"))
  safe_save_kable(lbo_cv_14_result$res$metrics, file.path(tab_dir, "metrics_T14_base.html"), "LBO-CV metrics — T14")
}
if (exists("lbo_cv_05_result")) {
  safe_write_csv(lbo_cv_05_result$res$metrics, file.path(tab_dir, "metrics_T05_base.csv"))
  safe_save_kable(lbo_cv_05_result$res$metrics, file.path(tab_dir, "metrics_T05_base.html"), "LBO-CV metrics — T05")
}
if (exists("bench14")) {
  safe_write_csv(bench14$metrics, file.path(tab_dir, "metrics_T14_tuned.csv"))
  safe_save_kable(bench14$metrics, file.path(tab_dir, "metrics_T14_tuned.html"), "Metrics — tuned @ R* (T14)")
}
if (exists("bench05")) {
  safe_write_csv(bench05$metrics, file.path(tab_dir, "metrics_T05_tuned.csv"))
  safe_save_kable(bench05$metrics, file.path(tab_dir, "metrics_T05_tuned.html"), "Metrics — tuned @ R* (T05)")
}

# Skalen & R*
if (exists("tune14") && exists("Ls14")) {
  safe_write_csv(tune14$grid, file.path(tab_dir, "Ucurve_T14.csv"))
  safe_write_csv(data.frame(L50 = Ls14$L50, L95 = Ls14$L95, R_star = tune14$R_star),
                 file.path(tab_dir, "scales_T14.csv"))
}
if (exists("tune05") && exists("Ls05")) {
  safe_write_csv(tune05$grid, file.path(tab_dir, "Ucurve_T05.csv"))
  safe_write_csv(data.frame(L50 = Ls05$L50, L95 = Ls05$L95, R_star = tune05$R_star),
                 file.path(tab_dir, "scales_T05.csv"))
}

# Optional: Error-Budget (falls berechnet)
# if (exists("eb14")) safe_write_csv(eb14, file.path(tab_dir, "error_budget_T14.csv"))
# if (exists("eb05")) safe_write_csv(eb05, file.path(tab_dir, "error_budget_T05.csv"))

# ---------- Raster exportieren ----------
if (exists("scen")) {
  try(terra::writeRaster(scen$E,   file.path(ras_dir, "E_dem.tif"),   overwrite = TRUE), silent = TRUE)
  try(terra::writeRaster(scen$R14, file.path(ras_dir, "R14_truth.tif"), overwrite = TRUE), silent = TRUE)
  try(terra::writeRaster(scen$R05, file.path(ras_dir, "R05_truth.tif"), overwrite = TRUE), silent = TRUE)
  if ("lc" %in% names(scen)) try(terra::writeRaster(scen$lc, file.path(ras_dir, "landcover.tif"),
                                                    overwrite = TRUE), silent = TRUE)
}
if (exists("bench14") && "E_star" %in% names(bench14))
  try(terra::writeRaster(bench14$E_star, file.path(ras_dir, "E_star_T14.tif"), overwrite = TRUE), silent = TRUE)
if (exists("bench05") && "E_star" %in% names(bench05))
  try(terra::writeRaster(bench05$E_star, file.path(ras_dir, "E_star_T05.tif"), overwrite = TRUE), silent = TRUE)

# ---------- Sessioninfo als Referenz ----------
try(saveRDS(sessionInfo(), file.path(out_dir, "sessionInfo.rds")), silent = TRUE)

message("✔ Export fertig. Siehe Ordner: ", normalizePath(out_dir))
# ========================================================================


