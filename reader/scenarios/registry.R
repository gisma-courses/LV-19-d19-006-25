# block4_5/scenarios/registry.R
if (!exists("scen_registry", inherits = FALSE)) scen_registry <- list()

register_scenario <- function(name, path) {
  if (!file.exists(path)) stop("Datei nicht gefunden: ", path)
  scen_registry[[name]] <<- normalizePath(path, winslash = "/", mustWork = TRUE)
}

register_all_scenarios <- function(dir = here::here("scenarios")) {
  files <- list.files(dir, pattern = "\\.R$", full.names = TRUE)
  for (f in files) {
    env <- new.env(parent = emptyenv())
    sys.source(f, envir = env)
    if (exists("SCEN_NAME", envir = env) && exists("make", envir = env, inherits = FALSE)) {
      register_scenario(get("SCEN_NAME", env), f)
    }
  }
}

source_scenario <- function(name) {
  path <- scen_registry[[name]]
  if (is.null(path)) stop("Scenario nicht registriert: ", name)
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  if (!exists("make", envir = env, inherits = FALSE))
    stop("Scenario-Datei '", basename(path), "' muss eine make()-Funktion definieren.")
  get("make", envir = env)
}

# Beispiel-Registrierungen:
register_scenario("lake_bump_dense",  here::here("scenarios/lake_bump_dense.R"))
register_scenario("scen_scaled_demo", here::here("scenarios/scen_scaled_demo.R"))
# oder: register_all_scenarios()
