# Kopieren/umbenennen (z. B. scenario_meinsee.R) und anpassen
scenario_template <- function(template, ..., seed = 1) {
  stopifnot(inherits(template, "SpatRaster"))
  # Hier deine Parameter aus ... holen, Defaults setzen:
  p <- list(...)
  
  # Beispiel: p$see_radius <- p$see_radius %||% 400
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  p$see_radius <- p$see_radius %||% 400
  
  # --- baue dein DEM / LC / Felder ---
  # E <- ...
  # lc <- ...
  # R14 <- ...
  # R05 <- ...
  # st <- ...
  
  list(E = E, lc = lc, R14 = R14, R05 = R05, stations = st)
}

register_scenario("my_template", scenario_template)
