# =====================================================================
# block4_5/scenarios/scen_scaled_demo.R
# Works with your source_scenario() loader.
# Defines: SCEN_NAME, SCEN_DESC, make(...)
# Idea: Truth fields (T14/T05) are built from DEMs smoothed at known
#       metric radii R_true14 / R_true05 so that U-curve tuning
#       recovers these scales.
# =====================================================================

SCEN_NAME <- "scen_scaled_demo"
SCEN_DESC <- "Didactic scene with lake (right) and hill (left); T-fields built from DEM smoothed at known radii (R_true14/R_true05)."

# ------------------------- Defaults -----------------------------------
.defaults <- list(
  # Domain (UTM32N by default)
  center_E = 600000, center_N = 5725000,
  len_x = 1600, len_y = 1200, res = 10, crs = "EPSG:32632",
  
  # Geometric features
  lake_mode = "water",  # "none" | "water" | "hollow"
  hill_mode = "bump",   # "none" | "bump"
  lake_diam_m  = 220,   lake_depth_m = 10, smooth_edges = TRUE,
  hill_diam_m  = 600,   hill_height_m = 60, hill_smooth  = TRUE,
  
  # Micro-relief to inject sub-R* texture
  random_hills        = 120,
  micro_hill_diam_m   = 60,
  micro_hill_height_m = 12,
  micro_hill_smooth   = TRUE,
  micro_seed          = 1L,
  
  # >>> baked-in process scales (what the U-curve should find)
  R_true14 = 300,      # m (afternoon)
  R_true05 = 180,      # m (dawn)
  
  # Sun / geo
  lat = 51.8, lon = 10.6, sun_date = as.Date("2024-06-21"),
  
  # Stations
  station_mode = "random",  # "random" | "ns_transect" | "ew_transect"
  n_st = 60, transect_margin_m = 10, ns_offset_m = 0, ew_offset_m = 0,
  
  # Models & CV
  models = c("Voronoi","IDW","OK","KED","RF","GAM"),
  block_size = NA_real_    # NA -> compute automatically
)

# simple, side-effect free merge
.merge <- function(a, b) { a[names(b)] <- b; a }

# -------- local builder that bakes scales into the truth ----------------
.build_scenario_scaled <- function(
    domain,
    lake_mode, hill_mode,
    lake_diam_m, lake_depth_m, smooth_edges,
    hill_diam_m, hill_height_m, hill_smooth,
    random_hills, micro_hill_diam_m, micro_hill_height_m, micro_hill_smooth, micro_seed,
    R_true14, R_true05,
    lat, lon, sun_date
) {
  # Template raster from domain
  ext <- terra::ext(domain$xmin, domain$xmax, domain$ymin, domain$ymax)
  Rtemplate <- terra::rast(ext, resolution = domain$res, crs = domain$crs)
  
  xmin <- terra::xmin(ext); xmax <- terra::xmax(ext)
  ymin <- terra::ymin(ext); ymax <- terra::ymax(ext)
  len_x <- xmax - xmin;     len_y <- ymax - ymin
  x0 <- (xmin + xmax)/2;    y0 <- (ymin + ymax)/2
  
  # Raster coordinates (numeric)
  X <- terra::init(Rtemplate, "x")
  Y <- terra::init(Rtemplate, "y")
  
  # Base valley (parabolic)
  a  <- 120 / ((len_y/2)^2)     # ~120 m rim-to-center relief
  E  <- 500 + a * (Y - y0)^2
  names(E) <- "elev"
  
  # Lake / hollow (right third)
  x_lc <- xmin + 2*len_x/3; y_lc <- y0
  lr <- max(1e-6, lake_diam_m/2)
  rl <- sqrt((X - x_lc)^2 + (Y - y_lc)^2)
  w_l <- if (isTRUE(smooth_edges)) terra::ifel(rl <= lr, 1 - (rl/lr)^2, 0) else terra::ifel(rl <= lr, 1, 0)
  if (lake_mode %in% c("water","hollow")) E <- E - as.numeric(lake_depth_m) * w_l
  lakeR <- if (identical(lake_mode, "water")) terra::ifel(w_l > 0, 1L, 0L) else 0L * (0*X); names(lakeR) <- "lake"
  
  # Main hill (left third)
  x_hc <- xmin + len_x/3;  y_hc <- y0
  hr <- max(1e-6, hill_diam_m/2)
  rh <- sqrt((X - x_hc)^2 + (Y - y_hc)^2)
  w_h_main <- if (hill_mode == "bump") {
    if (isTRUE(hill_smooth)) exp(-(rh/hr)^2) else terra::ifel(rh <= hr, 1, 0)
  } else 0 * X
  E <- E + as.numeric(hill_height_m) * w_h_main
  
  # Micro-relief
  w_h_micro <- 0 * X
  if (random_hills > 0) {
    if (!is.null(micro_seed)) set.seed(micro_seed)
    margin <- micro_hill_diam_m/2 + 5
    hrm <- max(1e-6, micro_hill_diam_m/2)
    for (i in seq_len(random_hills)) {
      cx <- runif(1, xmin + margin, xmax - margin)
      cy <- runif(1, ymin + margin, ymax - margin)
      r  <- sqrt((X - cx)^2 + (Y - cy)^2)
      wi <- if (micro_hill_smooth) exp(-(r/hrm)^2) else terra::ifel(r <= hrm, 1, 0)
      w_h_micro <- terra::clamp(w_h_micro + wi, 0, 1)
    }
    E <- E + as.numeric(micro_hill_height_m) * w_h_micro
  }
  hillW <- terra::clamp(w_h_main + w_h_micro, 0, 1); names(hillW) <- "hillW"
  
  # Terrain derivatives
  slp <- terra::terrain(E, v="slope",  unit="radians")
  asp <- terra::terrain(E, v="aspect", unit="radians")
  
  # Sun geometry
  s14 <- sun_pos_utc(sun_date, 14L, lat, lon)
  s05 <- sun_pos_utc(sun_date,  5L, lat, lon)
  I14_raw <- cosi_fun(s14$alt, s14$az, slp, asp); names(I14_raw) <- "I14"
  I05_raw <- cosi_fun(s05$alt, s05$az, slp, asp); names(I05_raw) <- "I05"
  
  # Land cover: 1 forest, 2 water, 3 bare soil, 4 meadows
  lc_levels <- c("forest","water","bare soil","meadows")
  lc_colors <- c("forest"="#2E8B57","water"="#5DADE2","bare soil"="#C49A6C","meadows"="#7FBF7B")
  lc <- 4L + 0L*(0*X)
  lc <- terra::ifel(lakeR > 0, 2L, lc)
  forest_mask <- (hillW > 0.20) | (slp > 0.18 & (Y > y0))
  lc <- terra::ifel(forest_mask & (lakeR <= 0), 1L, lc)
  thr_slp <- stats::quantile(terra::values(slp)[,1], 0.90, na.rm = TRUE)
  bare_mask <- (slp >= thr_slp) & (lakeR <= 0) & (!forest_mask)
  lc <- terra::ifel(bare_mask, 3L, lc); lc <- terra::clamp(lc, 1L, 4L); names(lc) <- "lc"
  
  # -------- bake the process scales into truth via smoothed DEM ---------------
  fr14 <- smooth_dem_and_derive(E, s14$alt, s14$az, radius_m = R_true14)
  fr05 <- smooth_dem_and_derive(E, s05$alt, s05$az, radius_m = R_true05)
  
  # LC coefficients (simple but reproducible)
  alpha_I_by_lc <- c("forest"=3.5, "water"=1.5, "bare soil"=6.0, "meadows"=4.0)
  shade_fac_by_lc <- c("forest"=0.60, "water"=1.00, "bare soil"=1.00, "meadows"=0.95)
  dawn_bias_by_lc <- c("forest"=0.30, "water"=1.20, "bare soil"=-0.50, "meadows"=0.05)
  pool_fac_by_lc  <- c("forest"=0.70, "water"=0.80, "bare soil"=1.10, "meadows"=1.05)
  
  idx <- pmax(1L, pmin(4L, as.integer(terra::values(lc))))
  to_r <- function(x) terra::setValues(terra::rast(E), x)
  alpha_I <- to_r(as.numeric(alpha_I_by_lc[lc_levels[idx]]))
  shade_f <- to_r(as.numeric(shade_fac_by_lc[lc_levels[idx]]))
  dawn_b  <- to_r(as.numeric(dawn_bias_by_lc[lc_levels[idx]]))
  pool_f  <- to_r(as.numeric(pool_fac_by_lc [lc_levels[idx]]))
  
  # Pooling pattern across valley
  dist2ax <- abs(Y - (terra::ymax(E)+terra::ymin(E))/2)
  pool_base <- 4.0 * exp(-(dist2ax / 70)^2)
  pool_mod  <- pool_base * (1 - 0.4 * hillW) * pool_f
  
  # Small reproducible noise
  set.seed(2001); noise14 <- terra::setValues(terra::rast(E), rnorm(terra::ncell(E), 0, 0.25))
  set.seed(2002); noise05 <- terra::setValues(terra::rast(E), rnorm(terra::ncell(E), 0, 0.25))
  
  # Truth fields (from *smoothed* features)
  E_mean <- terra::global(E, "mean", na.rm = TRUE)[1,1]
  T0_14 <- 26.0; lapse_14 <- -0.0065
  T0_05 <-  8.5; inv_05   <-  0.003; eta_slope <- 0.6
  
  I14_eff <- fr14$cosi * shade_f
  R14 <- T0_14 + lapse_14 * (fr14$Es - E_mean) + alpha_I * I14_eff + noise14; names(R14) <- "T14"
  R05 <- T0_05 + inv_05   * (fr05$Es - E_mean) + eta_slope * fr05$slp - pool_mod + dawn_b + noise05; names(R05) <- "T05"
  
  # Return scenario list
  list(
    E = E, slp = slp, asp = asp,
    I14 = I14_raw, I05 = I05_raw,
    R14 = R14, R05 = R05,
    lake = lakeR, hillW = hillW,
    lc = lc, lc_levels = lc_levels, lc_colors = lc_colors,
    sun = list(T14 = list(alt = s14$alt, az = s14$az),
               T05 = list(alt = s05$alt, az = s05$az)),
    R_true = list(T14 = R_true14, T05 = R_true05)
  )
}

# ---------------------------------------------------------------------
# make(overrides = list(), do_cv = FALSE)
# Builds: domain -> scenario -> stations -> (optional) CV
# ---------------------------------------------------------------------
make <- function(overrides = list(), do_cv = FALSE) {
  p <- .merge(.defaults, overrides)
  
  # 1) Domain
  domain <- make_domain(
    center_E = p$center_E, center_N = p$center_N,
    len_x = p$len_x, len_y = p$len_y, res = p$res, crs = p$crs
  )
  
  # 2) Scenario (with baked scales)
  scen <- .build_scenario_scaled(
    domain       = domain,
    lake_mode    = p$lake_mode,
    hill_mode    = p$hill_mode,
    lake_diam_m  = p$lake_diam_m,
    lake_depth_m = p$lake_depth_m,
    smooth_edges = p$smooth_edges,
    hill_diam_m  = p$hill_diam_m,
    hill_height_m= p$hill_height_m,
    hill_smooth  = p$hill_smooth,
    random_hills        = p$random_hills,
    micro_hill_diam_m   = p$micro_hill_diam_m,
    micro_hill_height_m = p$micro_hill_height_m,
    micro_hill_smooth   = p$micro_hill_smooth,
    micro_seed          = p$micro_seed,
    R_true14 = p$R_true14, R_true05 = p$R_true05,
    lat = p$lat, lon = p$lon, sun_date = p$sun_date
  )
  
  # 3) Stations
  pts_sf <- make_stations(
    domain,
    n_st = p$n_st,
    station_mode = p$station_mode,
    transect_margin_m = p$transect_margin_m,
    ns_offset_m = p$ns_offset_m,
    ew_offset_m = p$ew_offset_m
  )
  
  # 4) Station features/targets
  stns <- stations_from_scenario(scen, pts_sf)
  stn_sf_14 <- stns$T14
  stn_sf_05 <- stns$T05
  
  # 5) Block size
  block_size <- if (is.finite(p$block_size)) {
    as.numeric(p$block_size)
  } else {
    compute_block_size(
      len_x = domain$xmax - domain$xmin,
      len_y = domain$ymax - domain$ymin,
      n_st  = p$n_st
    )
  }
  
  out <- list(
    name       = SCEN_NAME,
    desc       = SCEN_DESC,
    params     = p,
    domain     = domain,
    scen       = scen,
    pts_sf     = pts_sf,
    stn_sf_14  = stn_sf_14,
    stn_sf_05  = stn_sf_05,
    block_size = block_size
  )
  
  # 6) Optional CV
  if (isTRUE(do_cv)) {
    out$cv <- list(
      T14 = run_for_time(stn_sf_14, scen$R14, "T14",
                         scen_local = scen, block_m = block_size, models = p$models),
      T05 = run_for_time(stn_sf_05, scen$R05, "T05",
                         scen_local = scen, block_m = block_size, models = p$models)
    )
  }
  
  out
}
# =====================================================================
