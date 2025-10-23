# =====================================================================
# block4_5/scenarios/lake_bump_dense.R
# Benötigt: source("block4_5/src/pipemodel_functions.R") davor
# Definiert: SCEN_NAME, SCEN_DESC, make(...)
# =====================================================================

SCEN_NAME <- "lake_bump_dense"
SCEN_DESC <- "Valley with lake (right), bump hill (left) + dense micro-hills; random stations."

# Defaults so wie bisher genutzt (nur Settings/Parameter)
.defaults <- list(
  # Domain
  center_E = 600000, center_N = 5725000,
  len_x = 600, len_y = 400, res = 10, crs = "EPSG:32632",
  
  # Topographie-Features
  lake_mode = "water",   # "none" | "water" | "hollow"
  hill_mode = "bump",    # "none" | "bump"
  lake_diam_m  = 100, lake_depth_m = 500, smooth_edges = FALSE,
  hill_diam_m  = 100, hill_height_m = 500, hill_smooth  = FALSE,
  
  # micro-relief
  random_hills        = 50,
  micro_hill_diam_m   = 30,
  micro_hill_height_m = 50,
  micro_hill_smooth   = FALSE,
  micro_seed          = NULL,
  

  # Sonne/Geo
  lat = 51.8, lon = 10.6, sun_date = as.Date("2024-06-21"),
  
  # Stationen
  station_mode = "random",  # "random" | "ns_transect" | "ew_transect"
  n_st = 60,
  transect_margin_m = 10,
  ns_offset_m = 0,
  ew_offset_m = 0,
  
  # Modelle + Block-CV
  models = c("Voronoi","IDW","OK","KED","RF","GAM"),
  block_size = NA_real_    # wenn NA -> automatisch berechnet
)

# einfaches Mergen (ohne Seiteneffekte)
.merge <- function(a, b) { a[names(b)] <- b; a }

# ---------------------------------------------------------------------
# make(overrides = list(), do_cv = FALSE)
# baut Domain -> Szenario -> Stationen -> (optional) CV
# nutzt ausschließlich Funktionen aus pipemodel_functions.R
# ---------------------------------------------------------------------
make <- function(overrides = list(), do_cv = FALSE) {
  p <- .merge(.defaults, overrides)
  
  # 1) Domain
  domain <- make_domain(
    center_E = p$center_E, center_N = p$center_N,
    len_x = p$len_x, len_y = p$len_y, res = p$res, crs = p$crs
  )
  
  # 2) Szenario (Topographie/Physikfelder)
  scen <- build_scenario(
    domain       = domain,
    lake_mode    = p$lake_mode,
    hill_mode    = p$hill_mode,
    
    # <<< diese Zeilen fehlten bisher
    lake_diam_m  = p$lake_diam_m,
    lake_depth_m = p$lake_depth_m,
    smooth_edges = p$smooth_edges,
    hill_diam_m  = p$hill_diam_m,
    hill_height_m= p$hill_height_m,
    hill_smooth  = p$hill_smooth,
    # >>>
    
    # micro-relief
    random_hills        = p$random_hills,
    micro_hill_diam_m   = p$micro_hill_diam_m,
    micro_hill_height_m = p$micro_hill_height_m,
    micro_hill_smooth   = p$micro_hill_smooth,
    micro_seed          = p$micro_seed,
    
    # Sonne/Geo
    lat = p$lat, lon = p$lon, sun_date = p$sun_date
  )
  
  # 3) Stationen
  pts_sf <- make_stations(
    domain,
    n_st = p$n_st,
    station_mode = p$station_mode,
    transect_margin_m = p$transect_margin_m,
    ns_offset_m = p$ns_offset_m,
    ew_offset_m = p$ew_offset_m
  )
  
  # 4) Station-Features/Targets extrahieren
  stns <- stations_from_scenario(scen, pts_sf)
  stn_sf_14 <- stns$T14
  stn_sf_05 <- stns$T05
  
  # 5) Blockgröße
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
  
  # 6) Optional: Block-CV (nur wenn gewünscht)
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
