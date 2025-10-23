# -----------------------------------------------------------------------------
# Control Script — Interpolate Ecowitt Air Temperature (Toolkit-driven)
# -----------------------------------------------------------------------------
# Purpose (mirrors the pipeline goals)
#   1) Data ingestion & cleaning (Ecowitt)
#   2) Spatial interpolation (per-timestep KED preview + method compare)
#   3) Variogram-based scale inference (L50/L95) & drift radius selection (R)
#   4) Scale-matched predictor construction (DEM + optional rasters)
#   5) Tuning R* via block-CV (U-curve) & method benchmarking
#   6) Diagnostics (panel plot; optional error budget)
#
# This script assumes the consolidated function set is saved as 'all_functions_1.R'.
# If you have the toolkit open in the Canvas, copy it to that file path first.
# -----------------------------------------------------------------------------

# 0 — Setup -------------------------------------------------------------------

pkgs <- c(
  "sf","terra","raster","dplyr","automap","gstat","mapview","stars",
  "readxl","stringr","tidyr","purrr","lubridate","rprojroot",
  "exactextractr","zoo","ggplot2","viridis","mgcv","randomForest","fields","sp","deldir"
)
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, dependencies = TRUE)
invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

# Project root and paths -------------------------------------------------------
wd <- tryCatch(find_rstudio_root_file(), error = function(...) getwd())
source(file.path(wd, "block4_5/all_functions_1.R"))   # <- loads the toolkit (previous Canvas file)

fn_DTM        <- file.path(wd, "block4_5/data_2024/copernicus_DEM.tif")
fn_stations   <- file.path(wd, "block4_5/data_2024/stations_prelim_modifiziert.gpkg")
fn_area       <- file.path(wd, "block4_5/data_2024/plot.shp")
fn_temp_FC29  <- file.path(wd, "block4_5/data_2024/all_GW1000A-WIFIFC29.xlsx")
fn_temp_DB2F  <- file.path(wd, "block4_5/data_2024/all_GW1000A-WIFIDB2F.xlsx")
cleandata_rds <- file.path(wd, "block4_5/data_2024/climdata.RDS")

out_dir   <- file.path(wd, "block4_5/interpolated")
fig_dir   <- file.path(out_dir, "fig")
method_dir<- file.path(out_dir, "methods_compare")
report_dir <- file.path(out_dir, "report")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

for (d in c(out_dir, fig_dir, method_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# CRS settings ----------------------------------------------------------------
epsg <- "EPSG:32633"  # UTM zone 33N
sf_crs_utm33 <- st_crs(epsg)

# 1 — Load & prepare base data -------------------------------------------------


# --- DEMs ---
DEM_scale  <- terra::rast("block4_5/data_2024/DEM.tif") |> terra::project(epsg)          # native
DEM_scale <- terra::aggregate(DEM_scale, c(20, 20))               # optional coarsen (e.g. ~20–25 m)
names(DEM_scale) <- "altitude"
DEM_render <- DEM_scale |> terra::aggregate(fact = c(10, 10))         # only for rendering products
# Stations and plot boundary → same CRS as DEM
stations_pos <- st_read(fn_stations) |> st_transform(sf_crs_utm33)
plot_area    <- st_read(fn_area)     |> st_transform(sf_crs_utm33)
plot_area    <- sf::st_make_valid(plot_area)  # geometry guard
  
# 2 — Ingest & clean Ecowitt data ---------------------------------------------

temp_FC29 <- extract_ecowitt_core_vars(fn_temp_FC29)
temp_DB2F <- extract_ecowitt_core_vars(fn_temp_DB2F)

t_rh_all <- merge_ecowitt_logger_vars(temp_FC29, temp_DB2F)

# Clean display names and map to verbose station names
for (meas in c("temperature","humidity")) {
  t_rh_all[[meas]] <- t_rh_all[[meas]] %>%
    rename_with(~ to_verbose(.x, ifelse(meas=="temperature","Temperature","Humidity")), -Time) %>%
    clean_names()
}

# 3 — Aggregate to 3-hour and reshape wide ------------------------------------

temp_agg <- t_rh_all$temperature %>%
  mutate(time = floor_date(Time, "3 hours")) %>%
  group_by(time) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

names(temp_agg) <- clean_ids(names(temp_agg))

# long → wide matrix: station rows, time columns
temp_matrix <- temp_agg %>%
  pivot_longer(cols = -time, names_to = "stationid", values_to = "value") %>%
  pivot_wider(names_from = time, values_from = value)

# 4 — Join altitude to stations, combine with values ---------------------------

stations_pos <- stations_pos %>%
  mutate(altitude = exactextractr::exact_extract(DEM_scale, st_buffer(stations_pos, 1), "mean")) %>%
  mutate(stationid = to_verbose(stationid))

m <- left_join(stations_pos, temp_matrix, by = "stationid")

# Name hygiene for station IDs and final column names
stations_pos$stationid <- gsub("\\(℃\\)|\\(％\\)|\\(\\%\\)", "", stations_pos$stationid)
m$stationid            <- gsub("\\(℃\\)|\\(％\\)|\\(\\%\\)", "", m$stationid)

names(m) <- fix_names(names(m))
saveRDS(m, cleandata_rds)

# 5 — Interpolation (per-timestep KED via convenience wrapper) -----------------
min_pts <- 5
vars <- as.list(grep("^A\\d{8,14}", names(m), value = TRUE))
kriged_list <- lapply(vars, function(v) interpolate_kriging(v, m, DEM_render, output_dir = out_dir, label = "pretty"))
names(kriged_list) <- vars  # <- WICHTIG: list names = Zeit-Keys (AYYYYMMDDHHMMSS)

# 6 — Multi-panel visualization (shared scale) ---------------------------------

panel <- timeseries_panel(
  kriged_list        = kriged_list,
  plot_area          = plot_area,
  stations_pos       = stations_pos,
  cells_target       = 150000,
  max_cols           = 4,
  label_pretty_time  = TRUE,
  out_png            = file.path(fig_dir, "timeseries_panel_grid.png"),
  out_pdf            = file.path(fig_dir, "timeseries_panel_grid.pdf"),
  fill_label         = "Temperature"
)
print(panel$plot)

# 8 run one ----------------------------------------------------------------------------

# Ein bestimmter time stamp (z.B. der dichteste)
pick_idx <- pick_densest_index(m, vars)
pick_ts = names(m)[pick_idx]
# Extra predictors (extend/replace as you like)
extra_list <- list(
  slope  = terra::terrain(DEM_scale, v = "slope",  unit = "degrees"),
  aspect = terra::terrain(DEM_scale, v = "aspect"),
  tri    = terra::terrain(DEM_scale, v = "TRI")
)
# Ein Timestamp (der dichteste), MIT Extras:
res_one <- run_one(
  v          = vars[[pick_idx]],
  m          = m,
  DEM_render = DEM_render,
  DEM_scale  = DEM_scale,
  method_dir = method_dir,
  fig_dir    = fig_dir,
  report_dir = report_dir,
  extra_preds = extra_list,   # oder NULL, wenn ohne Extras
  save_figs   = TRUE
)



# Alle Slots
# Alle Time-Steps rechnen + zusammenfassen
compute_all <- TRUE

if (isTRUE(compute_all)) {
  # --- Mindest-Vars prüfen (werden in run_one sonst via .need gezogen)
  req <- c("m","DEM_render","DEM_scale","method_dir","fig_dir","report_dir")
  miss <- req[!vapply(req, exists, logical(1), inherits = TRUE)]
  if (length(miss)) stop("Fehlende Objekte im Environment: ", paste(miss, collapse = ", "))
  
  # Kleiner Helfer: bestes Verfahren aus einer Bench-Tabelle ziehen
  .best_from_bench <- function(bench_obj) {
    if (is.null(bench_obj) || !is.data.frame(bench_obj$table) || nrow(bench_obj$table) < 1)
      return(NULL)
    b <- bench_obj$table
    # Nur endliche RMSE behalten
    b <- b[is.finite(b$RMSE), , drop = FALSE]
    if (!nrow(b)) return(NULL)
    b <- b[order(b$RMSE), , drop = FALSE]
    b[1, c("method","RMSE"), drop = FALSE]
  }
  
  
  # Fortschrittsausgabe
  message(sprintf("Starte compute_all für %d Zeitschritte …", length(vars)))
  
  res_all <- setNames(lapply(vars, function(vv) {
    message("→ run_one: ", pretty_time(vv))
    tryCatch(
      run_one(
        v          = vv,
        m          = m,
        DEM_render = DEM_render,
        DEM_scale  = DEM_scale,
        method_dir = method_dir,
        fig_dir    = fig_dir,
        report_dir = report_dir,
        extra_preds = extra_list,     # setze auf NULL, falls ohne Extras
        save_figs   = TRUE,
        save_tables = TRUE
      ),
      error = function(e) {
        warning("run_one fehlgeschlagen für ", vv, ": ", conditionMessage(e))
        NULL
      }
    )
  }), vars)
  
  # Komplettes Objekt sichern
  saveRDS(res_all, file.path(report_dir, "all_results.RDS"))
  
  # Zusammenfassung: R* + beste Methode (vergleicht bench vs bench_ex)
  # Zusammenfassung: R* + beste Methode (vergleicht bench vs bench_ex)
  summ <- do.call(rbind, lapply(names(res_all), function(k) {
    r <- res_all[[k]]
    if (is.null(r)) {
      return(data.frame(
        ts_key      = k, 
        stamp       = pretty_time(k),
        R_star      = NA_real_,
        best_source = NA_character_,
        best_method = NA_character_,
        best_RMSE   = NA_real_
      ))
    }
    rstar <- suppressWarnings(as.numeric(r$tune$R_star))
    if (!is.finite(rstar)) rstar <- NA_real_
    
    b0 <- .best_from_bench(r$bench)
    bE <- .best_from_bench(r$bench_ex)
    
    score0 <- if (!is.null(b0) && isTRUE(is.finite(b0$RMSE))) b0$RMSE else Inf
    scoreE <- if (!is.null(bE) && isTRUE(is.finite(bE$RMSE))) bE$RMSE else Inf
    
    if (is.infinite(score0) && is.infinite(scoreE)) {
      src <- NA_character_; bm <- NA_character_; br <- NA_real_
    } else if (score0 <= scoreE) {
      src <- "no_extras"; bm <- b0$method; br <- score0
    } else {
      src <- "with_extras"; bm <- bE$method; br <- scoreE
    }
    
    data.frame(
      ts_key      = k,
      stamp       = pretty_time(k),
      R_star      = rstar,
      best_source = src,
      best_method = bm,
      best_RMSE   = br
    )
  }))
  
  utils::write.csv(summ, file.path(report_dir, "summary_Rstar_bestmethod.csv"), row.names = FALSE)
  message("✔ Fertig: summary_Rstar_bestmethod.csv geschrieben.")
}

# ---- (Optional) CSV-Sicherungen für den EINEN pick_ts: aus res_one ziehen ----
ts_label <- pretty_time(pick_ts)
bench_base_csv <- file.path(report_dir, sprintf("benchmark_%s.csv",        slug(ts_label)))
bench_ex_csv   <- file.path(report_dir, sprintf("benchmark_extras_%s.csv", slug(ts_label)))
eb_base_csv    <- file.path(report_dir, sprintf("error_budget_%s.csv",     slug(ts_label)))
eb_ex_csv      <- file.path(report_dir, sprintf("error_budget_extras_%s.csv", slug(ts_label)))

if (is.list(res_one$bench)    && is.data.frame(res_one$bench$table))    write.csv(res_one$bench$table,    bench_base_csv, row.names = FALSE)
if (is.list(res_one$bench_ex) && is.data.frame(res_one$bench_ex$table)) write.csv(res_one$bench_ex$table, bench_ex_csv,   row.names = FALSE)
if (is.data.frame(res_one$errtab))    write.csv(res_one$errtab,    eb_base_csv, row.names = FALSE)
if (is.data.frame(res_one$errtab_ex)) write.csv(res_one$errtab_ex, eb_ex_csv,   row.names = FALSE)

# ---- Console summary: auf res_one verweisen ----
n_stations <- nrow(stations_pos)
n_pts_ts   <- sum(is.finite(m[[pick_ts]]))
Ls         <- get_Ls(res_one$wf$L)
Ls_e       <- if (!is.null(res_one$wf_ex)) get_Ls(res_one$wf_ex$L) else NULL
Rstar_base <- suppressWarnings(as.numeric(res_one$tune$R_star))
Rstar_ex   <- suppressWarnings(as.numeric(res_one$tune_ex$R_star))

cat(sprintf("Chosen R (micro/local): %s / %s m\n",
            ifelse(is.finite(res_one$wf$R['micro']), round(res_one$wf$R['micro']), "NA"),
            ifelse(is.finite(res_one$wf$R['local']), round(res_one$wf$R['local']), "NA")
))

best_row <- function(tab) if (is.data.frame(tab) && nrow(tab)) tab[which.min(tab$RMSE), , drop = FALSE] else NULL
best_b <- best_row(res_one$bench$table)
best_e <- if (!is.null(res_one$bench_ex)) best_row(res_one$bench_ex$table) else NULL

# ---- Viewer: Objekte aus res_one durchreichen ----
explanations= build_explanations(fig_dir = fig_dir,pick_ts =vars[[pick_idx]] )
run_mc_viewer(
  vars         = vars,
  method_dir   = method_dir,
  fig_dir      = fig_dir,
  stations_pos = stations_pos,
  plot_area    = plot_area,
  wf = res_one$wf,
  wf_ex = res_one$wf_ex,
  tune = res_one$tune,
  tune_ex = res_one$tune_ex,
  bench = res_one$bench,
  bench_ex = res_one$bench_ex,
  tab_err = res_one$errtab,
  tab_err_ex = res_one$errtab_ex,
  explanations = explanations
)


make_mc_report(
  res = res_one,
  out_pdf = file.path(report_dir, sprintf("MC_Report_%s.pdf", slug(res_one$stamp)))
)  
