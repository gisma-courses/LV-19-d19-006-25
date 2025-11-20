get_and_rasterize_clc5 <- function(aoi_wgs84,
                                   dem,
                                   field    = "clc18",
                                   vec_out  = here::here("data", "processed", "clc5_2018_burgwald_wfs.gpkg"),
                                   rast_out = here::here("data", "processed", "clc5_2018_burgwald_demgrid.tif")) {
  # aoi_wgs84 : sf oder sfc in EPSG:4326 (z.B. aoi_burgwald_wgs)
  # dem       : terra::SpatRaster (z.B. dem_burgwald, CRS = EPSG:25832)
  # field     : Attributname für Rasterisierung (bei dir: "clc18")
  
  # --- 0) AOI normalisieren ---
  if (inherits(aoi_wgs84, "sf")) {
    aoi_geom <- sf::st_geometry(aoi_wgs84)
  } else if (inherits(aoi_wgs84, "sfc")) {
    aoi_geom <- aoi_wgs84
  } else {
    stop("Argument 'aoi_wgs84' muss ein sf- oder sfc-Polygon sein (EPSG:4326).")
  }
  
  if (is.na(sf::st_crs(aoi_geom))) {
    stop("AOI hat keinen CRS. Bitte aoi_wgs84 mit st_set_crs(..., 4326) versehen.")
  }
  if (sf::st_crs(aoi_geom)$epsg != 4326) {
    aoi_geom <- sf::st_transform(aoi_geom, 4326)
  }
  
  if (!inherits(dem, "SpatRaster")) {
    stop("Argument 'dem' muss ein terra::SpatRaster sein.")
  }
  
  bb <- sf::st_bbox(aoi_geom)
  
  # --- 1) Zwei mögliche BBOX-Varianten bauen ---
  # Variante A: lat,lon,lat,lon (entspricht BKG-Beispiel)
  bbox_latlon <- paste0(
    bb$ymin, ",",
    bb$xmin, ",",
    bb$ymax, ",",
    bb$xmax, ",EPSG:4326"
  )
  
  # Variante B: lon,lat,lon,lat (klassisch x,y,x,y)
  bbox_lonlat <- paste0(
    bb$xmin, ",",
    bb$ymin, ",",
    bb$xmax, ",",
    bb$ymax, ",EPSG:4326"
  )
  
  wfs_base <- "https://sgx.geodatenzentrum.de/wfs_clc5_2018?"
  
  build_url <- function(bbox_str) {
    params <- list(
      SERVICE  = "WFS",
      VERSION  = "1.1.0",
      REQUEST  = "GetFeature",
      TYPENAME = "clc5_2018:clc5",
      BBOX     = bbox_str,
      SRSNAME  = "EPSG:25832"
    )
    paste0(
      wfs_base,
      paste(names(params), params, sep = "=", collapse = "&")
    )
  }
  
  # --- 2) Erster Versuch: BBOX als lat,lon,... ---
  wfs_url_A <- build_url(bbox_latlon)
  message("Hole CLC5-2018 per WFS (Variante A: lat,lon):\n  ", wfs_url_A)
  clc_vec <- try(sf::read_sf(wfs_url_A), silent = TRUE)
  
  if (inherits(clc_vec, "try-error") || nrow(clc_vec) == 0) {
    message("Variante A ergab keine Features – versuche Variante B (lon,lat).")
    
    # Zweiter Versuch: BBOX als lon,lat,...
    wfs_url_B <- build_url(bbox_lonlat)
    message("Hole CLC5-2018 per WFS (Variante B: lon,lat):\n  ", wfs_url_B)
    clc_vec <- sf::read_sf(wfs_url_B)
    
    if (nrow(clc_vec) == 0) {
      stop("WFS-Request mit beiden BBOX-Varianten ergab keine CLC-Features im AOI.")
    }
    
    used_url <- wfs_url_B
  } else {
    used_url <- wfs_url_A
  }
  
  # --- 3) CRS an DEM anpassen + Geometrien säubern ---
  clc_vec <- sf::st_transform(clc_vec, crs = terra::crs(dem))
  
  # Z/M-Koordinaten entfernen
  clc_vec <- sf::st_zm(clc_vec, drop = TRUE, what = "ZM")
  # Geometrien reparieren
  clc_vec <- sf::st_make_valid(clc_vec)
  # Alles auf MULTIPOLYGON casten (terra mag das am liebsten)
  clc_vec <- sf::st_cast(clc_vec, "MULTIPOLYGON", warn = FALSE)
  
  # Vektor speichern
  sf::st_write(clc_vec, vec_out, delete_dsn = TRUE)
  message("Saved CLC5-2018 vector to: ", vec_out)
  
  # --- 4) Rasterisierung ---
  clc_v <- terra::vect(clc_vec)
  
  if (!(field %in% names(clc_v))) {
    warning("Feld '", field, "' nicht gefunden. Verwende konstanten Wert 1 zur Rasterisierung.")
    clc_rast <- terra::rasterize(clc_v, dem, field = 1)
  } else {
    clc_rast <- terra::rasterize(clc_v, dem, field = field)
  }
  
  terra::writeRaster(clc_rast, rast_out, overwrite = TRUE)
  message("Saved CLC5-2018 raster (DEM grid) to: ", rast_out)
  
  invisible(list(
    clc_vector = clc_vec,
    clc_raster = clc_rast,
    wfs_url    = used_url
  ))
}


# ------------------------------------
# 2) Helper: download_if_missing()
# ------------------------------------
download_if_missing <- function(url, destfile, mode = "wb") {
  if (!file.exists(destfile)) {
    message("Downloading:\n  ", url, "\n  -> ", destfile)
    download.file(url, destfile = destfile, mode = mode, quiet = FALSE)
  } else {
    message("File already exists, skipping download:\n  ", destfile)
  }
}
