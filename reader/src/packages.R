# --- Paketliste an EINER Stelle pflegen ---------------------------------
.req_pkgs <- list(
  core      = c("terra","sf","suncalc","gstat"),
  modeling  = c("randomForest","mgcv"),
  wrangling = c("dplyr","tibble","tidyr"),
  viz       = c("ggplot2","scales","patchwork","RColorBrewer"),
  report    = c("knitr","kableExtra","here","zoo", "gt", "openxlsx", "writexl")

)

ensure_packages <- function(pkgs = unlist(.req_pkgs)) {
  inst <- rownames(installed.packages())
  missing <- setdiff(pkgs, inst)
  if (length(missing)) install.packages(missing, dependencies = TRUE)
  invisible(lapply(pkgs, require, character.only = TRUE))
}

after_load <- function() {
  if (requireNamespace("sf", quietly = TRUE)) sf::sf_use_s2(FALSE)  # wie bisher
}

# Aufruf:
ensure_packages()
after_load()

# ---- Pfade ------------------------------------------------------------
base_dir <- tryCatch(here::here(), error = function(e) getwd())
src_dir  <- file.path(base_dir,  "src")
out_dir  <- file.path(base_dir,  "exports")
fig_dir  <- file.path(out_dir, "figs")
tab_dir  <- file.path(out_dir, "tables")
ras_dir  <- file.path(out_dir, "rasters")
dat_dir  <- file.path(out_dir, "data")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (d in c(fig_dir, tab_dir, ras_dir, dat_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

