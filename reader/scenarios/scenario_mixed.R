# block4_5/scenarios/mixed_scale.R
# Vollwertiges Szenario: Mixed-Scale (kleinräumige Insolation + mesoskaliges Pooling)
# Erwartet: Rtemplate, cosi_fun, sun14, sun05, lc_levels_default/lc_colors_default
# => kommen typischerweise aus block4_5/src/fun_pipemodel.R

# ---------- Kernbauer (dein bisheriger Code, minimal bereinigt) ----------
scenario_mixed_core <- function(
    template,                 # SpatRaster mit Ausdehnung/CRS/Res
    n_micro = 20,             # Anzahl Mikro-Hügel
    n_st    = 70,             # Stationsanzahl
    seed    = 42,             # Reproduzierbarkeit
    # LC-Koeffizienten (frei änderbar)
    alpha_I_by_lc = c(forest=3.0, water=1.5, "bare soil"=4.0, maize=2.5),
    shade_by_lc   = c(forest=0.6, water=1.0, "bare soil"=1.0, maize=0.9),
    dawn_by_lc    = c(forest=0.3, water=1.2, "bare soil"=-0.5, maize=0.1),
    pool_by_lc    = c(forest=0.7, water=0.8, "bare soil"=1.1, maize=1.0)
) {
  stopifnot(inherits(template, "SpatRaster"))
  R <- template
  ext <- terra::ext(R); x <- terra::init(R, "x"); y <- terra::init(R, "y")
  xmin <- terra::xmin(ext); xmax <- terra::xmax(ext)
  ymin <- terra::ymin(ext); ymax <- terra::ymax(ext)
  len_x <- xmax - xmin; y0 <- (ymin + ymax)/2; x0 <- (xmin + xmax)/2
  resx <- terra::res(R)[1]; crs_str <- terra::crs(R)
  
  gauss_bump <- function(x, y, cx, cy, radius_m, height_m, sign = +1) {
    sigma <- radius_m/2; r2 <- (x - cx)^2 + (y - cy)^2
    sign * height_m * exp(-r2 / (2 * sigma^2))
  }
  
  # DEM: Bowl (meso) + Hill (lokal) + mehrere Mikro-Hügel + schwache Ebene + Noise
  bowl <- gauss_bump(x, y, xmin + 2*len_x/3, y0, 1000, 60, -1)
  hill <- gauss_bump(x, y, xmin +   len_x/3, y0,  150, 50, +1)
  set.seed(seed+1); micro <- terra::rast(R); terra::values(micro) <- 0
  for (i in seq_len(n_micro)) {
    cx <- runif(1, xmin + 0.1*len_x, xmax - 0.1*len_x)
    cy <- runif(1, ymin + 0.1*len_x, ymax - 0.1*len_x)
    r  <- runif(1, 25, 60); h <- runif(1, 10, 30)
    micro <- micro + gauss_bump(x, y, cx, cy, r, h, +1)
  }
  slope <- 0.0005*(x - x0) - 0.0007*(y - y0)
  set.seed(seed+2); noise_dem <- terra::rast(R); terra::values(noise_dem) <- rnorm(terra::ncell(R), 0, 0.8)
  E <- slope + bowl + hill + micro + noise_dem; names(E) <- "E"
  
  # Landnutzung (1 forest, 2 water, 3 bare soil, 4 maize)
  lc <- terra::rast(R); terra::values(lc) <- 4L
  lc <- terra::ifel(((x - (xmin + 2*len_x/3))^2 + (y - y0)^2) <= (600^2), 2L, lc) # Wasser im Bowl
  lc <- terra::ifel(abs(x - (xmin + 0.45*len_x)) < 120, 1L, lc)                   # Wald-Streifen
  lc <- terra::ifel(abs(y - (ymin + 0.35*len_x)) <  90, 3L, lc)                   # Brach-Streifen
  lc <- terra::clamp(lc, 1, 4); names(lc) <- "lc"
  
  # Sonnengeometrie: hangabgeleitet
  slp <- terra::terrain(E, v="slope",  unit="radians")
  asp <- terra::terrain(E, v="aspect", unit="radians")
  
  # T14: Insolation (LC-spezifische Beschattung/Antwort) + leichte höhenabhängige Korrektur
  recode <- function(map, vec) terra::app(map, function(z) vec[pmax(1, pmin(4, z))])
  alphaI <- recode(lc, c(alpha_I_by_lc["forest"], alpha_I_by_lc["water"],
                         alpha_I_by_lc["bare soil"], alpha_I_by_lc["maize"]))
  shade  <- recode(lc, c(shade_by_lc["forest"], shade_by_lc["water"],
                         shade_by_lc["bare soil"], shade_by_lc["maize"]))
  set.seed(seed+3); noise14 <- terra::rast(R); terra::values(noise14) <- rnorm(terra::ncell(R), 0, 0.3)
  
  # Pooling (meso) via TPI ~600 m für T05
  w <- 2*round(600/resx) + 1
  E_mean600 <- terra::focal(E, w = w, fun = "mean", na.policy = "omit")
  TPI <- E - E_mean600
  pool_strength <- terra::clamp(-TPI, 0, Inf)
  pool_strength <- pool_strength / (terra::global(pool_strength, "max", na.rm=TRUE)[1] + 1e-6)
  poolF  <- recode(lc, c(pool_by_lc["forest"], pool_by_lc["water"],
                         pool_by_lc["bare soil"], pool_by_lc["maize"]))
  dawnB  <- recode(lc, c(dawn_by_lc["forest"],  dawn_by_lc["water"],
                         dawn_by_lc["bare soil"], dawn_by_lc["maize"]))
  set.seed(seed+4); noise05 <- terra::rast(R); terra::values(noise05) <- rnorm(terra::ncell(R), 0, 0.3)
  
  # Stationen (gleichmäßig zufällig)
  set.seed(seed+5)
  pts <- terra::as.points(terra::rast(R), values = FALSE)
  ids <- sample(nrow(pts), n_st)
  st  <- sf::st_as_sf(pts[ids,]) |> sf::st_set_crs(crs_str)
  
  list(E = E, slp = slp, asp = asp, lc = lc,
       alphaI = alphaI, shade = shade, poolF = poolF, dawnB = dawnB,
       noise14 = noise14, noise05 = noise05,
       stations = st)
}

# ---------- Factory: make() baut das komplette Szenario-Objekt ----------
make <- function(
    n_micro = as.integer(Sys.getenv("N_MICRO", 20)),
    n_st    = as.integer(Sys.getenv("N_ST",    70)),
    seed    = as.integer(Sys.getenv("SEED",    42)),
    alpha_I_by_lc = c(forest=3.0, water=1.5, "bare soil"=4.0, maize=2.5),
    shade_by_lc   = c(forest=0.6, water=1.0, "bare soil"=1.0, maize=0.9),
    dawn_by_lc    = c(forest=0.3, water=1.2, "bare soil"=-0.5, maize=0.1),
    pool_by_lc    = c(forest=0.7, water=0.8, "bare soil"=1.1, maize=1.0)
) {
  # Rtemplate & Utilities aus fun_pipemodel.R
  Rtemplate <- get0("Rtemplate", inherits = TRUE)
  if (is.null(Rtemplate)) stop("Rtemplate not found. Source block4_5/src/fun_pipemodel.R first.")
  cosi_fun  <- get0("cosi_fun", inherits = TRUE)
  if (is.null(cosi_fun)) stop("cosi_fun() not found (should come from fun_pipemodel.R).")
  sun14     <- get0("sun14", inherits = TRUE)
  sun05     <- get0("sun05", inherits = TRUE)
  if (is.null(sun14) || is.null(sun05)) stop("sun14/sun05 not found (fun_pipemodel.R).")
  
  core <- scenario_mixed_core(
    template = Rtemplate,
    n_micro = n_micro, n_st = n_st, seed = seed,
    alpha_I_by_lc = alpha_I_by_lc, shade_by_lc = shade_by_lc,
    dawn_by_lc = dawn_by_lc, pool_by_lc = pool_by_lc
  )
  
  # Insolation @ 14/05 aus Hang + Sonne
  I14 <- cosi_fun(sun14$alt, sun14$az, core$slp, core$asp)
  I05 <- cosi_fun(sun05$alt, sun05$az, core$slp, core$asp)
  
  # „Wahrheit“ berechnen (an core anknüpfen)
  E_mean <- terra::global(core$E, "mean", na.rm = TRUE)[1,1]
  R14 <- 20 + (core$alphaI * core$shade * I14) + 0.001*(core$E - E_mean) + core$noise14
  names(R14) <- "R14"
  R05 <- 10 + core$dawnB - 6*core$poolF*(terra::ifel(is.finite(core$poolF), 1, 1))* # skaliert über pool_strength
    (terra::ifel(is.finite(core$poolF), 1, 1))*0 + 0 # (Pool ist bereits in dawnB/poolF einkodiert)
  # -> Verwende explizit pool_strength:
  R05 <- 10 + core$dawnB - 6*core$poolF* (terra::mask(core$E, core$E, maskvalue=NA, updatevalue=1)*0 + 1) * 0 + 0
  # Korrigiert: nutze core$poolF * pool_strength:
  R05 <- 10 + core$dawnB - 6*core$poolF * {
    # pool_strength erneut berechnen (wie im core)
    resx <- terra::res(core$E)[1]
    w <- 2*round(600/resx) + 1
    E_mean600 <- terra::focal(core$E, w = w, fun = "mean", na.policy = "omit")
    TPI <- core$E - E_mean600
    ps <- terra::clamp(-TPI, 0, Inf)
    ps / (terra::global(ps, "max", na.rm=TRUE)[1] + 1e-6)
  } + 0.0005*(core$E - E_mean) + core$noise05
  names(R05) <- "R05"
  
  # LC-Metadaten
  lc_levels_default <- get0("lc_levels_default", ifnotfound = c("forest","water","bare soil","maize"), inherits = TRUE)
  lc_colors_default <- get0("lc_colors_default",
                            ifnotfound = c("forest"="#2E8B57","water"="#5DADE2","bare soil"="#C49A6C","maize"="#F4D03F"),
                            inherits = TRUE)
  
  # Szenario-Objekt (kompatibel zu predict_maps usw.)
  scen <- list(
    E = core$E, slp = core$slp, asp = core$asp,
    I14 = I14, I05 = I05,
    lc = core$lc,
    R14 = R14, R05 = R05,
    lc_levels = lc_levels_default,
    lc_colors = lc_colors_default
  )
  
  # Stationen -> stn_sf_14 / stn_sf_05 wie im Pipeline-Code
  st <- core$stations
  XY <- sf::st_coordinates(st)
  st$z_surf <- as.numeric(terra::extract(scen$E,   XY)[,1])
  st$slp_v  <- as.numeric(terra::extract(scen$slp, XY)[,1])
  st$I14_v  <- as.numeric(terra::extract(scen$I14, XY)[,1])
  st$I05_v  <- as.numeric(terra::extract(scen$I05, XY)[,1])
  st$lc_v   <- as.integer(terra::extract(scen$lc,  XY)[,1])
  st$T14_v  <- as.numeric(terra::extract(scen$R14, XY)[,1])
  st$T05_v  <- as.numeric(terra::extract(scen$R05, XY)[,1])
  
  # Einheiten wie im Original:
  to_stn <- function(df, which = c("T14","T05")) {
    which <- match.arg(which)
    cosi <- if (which == "T14") df$I14_v else df$I05_v
    temp <- if (which == "T14") df$T14_v  else df$T05_v
    lc   <- pmax(1L, pmin(df$lc_v, length(scen$lc_levels)))
    out <- tibble::tibble(
      id = seq_len(nrow(df)),
      x  = XY[,1], y = XY[,2],
      z_surf = as.numeric(df$z_surf),
      slp    = as.numeric(df$slp_v),
      cosi   = as.numeric(cosi),
      lc     = factor(scen$lc_levels[lc], levels = scen$lc_levels),
      temp   = as.numeric(temp)
    )
    sf::st_as_sf(out, coords = c("x","y"), crs = sf::st_crs(st), remove = FALSE)
  }
  
  stn_sf_14 <- to_stn(st, "T14")
  stn_sf_05 <- to_stn(st, "T05")
  
  list(
    scen      = scen,
    stn_sf_14 = stn_sf_14,
    stn_sf_05 = stn_sf_05,
    params    = list(n_micro = n_micro, n_st = n_st, seed = seed,
                     alpha_I_by_lc = alpha_I_by_lc, shade_by_lc = shade_by_lc,
                     dawn_by_lc = dawn_by_lc, pool_by_lc = pool_by_lc)
  )
}

# ---------- Optional: Preview-Helfer ----------
preview <- function(obj) {
  stopifnot(is.list(obj), !is.null(obj$scen))
  p <- plot_landcover_terrain(obj$scen, stations = obj$stn_sf_14, layout = "vertical")
  print(p)
  invisible(p)
}

# ---------- (Optional) Registry-Hook ----------
if (exists("register_scenario", mode = "function")) {
  register_scenario("mixed_scale", make)
}
