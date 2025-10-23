  # =====================================================================
  # main_ultra.R — minimal: run + (optional) save-at-end
  #
  # Purpose:
  #   1) Source packages + your function library + scenario registry
  #   2) Build scenario + station sets from a chosen scenario "make()" factory
  #   3) Live preview: domain/land-cover/terrain + 2x2 overview + scenario preview
  #   4) Baseline: leave-block-out CV (LBO-CV) and prediction maps (T14, T05)
  #   5) Scale inference: empirical variogram -> (L50,L95) -> U-curve -> R*
  #   6) Tuned CV & maps at R*
  #   7) OPTIONAL: save all plots/tables/rasters at the end
  #
  # Design notes:
  #   - This script stays "thin": all heavy lifting lives in fun_pipemodel.R
  #     and the scenario files. This keeps the main pipe reproducible and
  #     testable.
  #   - Keep side effects (saving files) to the very end; set `export <- FALSE`
  #     if you just want to run and eyeball plots interactively.
  # =====================================================================
  
  message("=== pipe main (ultra) ===")
  
  # Toggle: set to FALSE to only run/plot without saving anything
  export <- TRUE
  
  # ---------------------------------------------------------------------
  # 0) Setup & packages (centralized in your packages.R)
  #    - Loads CRAN pkgs (terra, sf, ggplot2, mgcv, gstat, suncalc, ...)
  #    - Sets knitr options (if used in a notebook context)
  #    - We also disable spherical geometry in sf to keep planar ops robust
  #      for projected domains (UTM in our scenarios).
  # ---------------------------------------------------------------------
  source(here::here("block4_5/src/packages.R"))
  sf::sf_use_s2(FALSE)
  set.seed(42)  # one seed here (scenarios may set more where needed)
  
  # ---------------------------------------------------------------------
  # 1) Functions + Scenario registry
  #    - fun_pipemodel.R: your full function library (NO side effects)
  #    - registry.R: maps scenario names to files; exposes source_scenario()
  #      which returns a `make()` function to build the object.
  # ---------------------------------------------------------------------
  source(here::here("block4_5/src/fun_pipemodel.R"))
  source(here::here("block4_5/src/fun_learn_predict_core.R"))
  source(here::here("block4_5/scenarios/registry.R"))
  
  # ---------------------------------------------------------------------
  # 2) Pick a scenario
  #    - Choose via env var SCEN (e.g., export SCEN=scen_scaled_demo)
  #    - Defaults to "lake_bump_dense" which is a realistic didactic scene
  #      (valley, lake, bump hill, dense micro-relief, LC = forest/water/bare/baseline meadow).
  # ---------------------------------------------------------------------
  scen_name <- Sys.getenv("SCEN", "lake_bump_dense")
  #
  
  scen_name <- Sys.getenv("SCEN", "scen_scaled_demo")
  make_fun  <- source_scenario(scen_name)  # returns a function make(overrides=list(), do_cv=FALSE)
  
  # ---------------------------------------------------------------------
  # 3) Build the object (domain + scenario + stations + params)
  #    - `obj` is a list with a stable contract used downstream:
  #        scen: list of rasters (E, R14, R05, I14, I05, lc, sun, ...)
  #        stn_sf_14 / stn_sf_05: station sf at 14/05 UTC (features + temp)
  #        block_size: integer (meters) used by LBO-CV
  #        params$models: character vector of model names to run
  # ---------------------------------------------------------------------
  obj  <- make_fun()
  scen <- obj$scen
  st14 <- obj$stn_sf_14
  st05 <- obj$stn_sf_05
  bs   <- obj$block_size
  mods <- obj$params$models
  
  # --- Safety checks that catch common wiring issues early ----------------------
  stopifnot(inherits(st14, "sf"), inherits(st05, "sf"))
  stopifnot(all(c("E","R14","R05") %in% names(scen)))
  
  # Sun geometry must be present for the R*-tuning (cosine-of-incidence features)
  if (is.null(scen$sun) || is.null(scen$sun$T14) || is.null(scen$sun$T05)) {
    stop("Scenario did not populate scen$sun$T14 / scen$sun$T05 (alt/az). ",
         "Fix in the scenario builder (build_scenario) before tuning.")
  }
  if (any(is.null(c(scen$sun$T14$alt, scen$sun$T14$az,
                    scen$sun$T05$alt, scen$sun$T05$az)))) {
    stop("Sun angles (alt/az) are NULL. Scenario must supply numeric alt/az for T14/T05.")
  }
  
  # ---------------------------------------------------------------------
  # 4) Live preview: plots during the run (no side effects)
  #    Why plot first?
  #      - Instant sanity checks: station placement, LC, illumination maps
  #      - Early visual cues if something is off (e.g., CRS mismatch)
  # ---------------------------------------------------------------------
  print(plot_landcover_terrain(scen, stations = st14, layout = "vertical"))
  print(plot_block_overview_2x2_en(scen, pts_sf = st14))
  # preview_scenario() may show truth fields, histograms, thumbnails, etc.
  print(preview_scenario(obj))  # accepts obj or obj$scen in the robust version
  
    # ---------------------------------------------------------------------
  # 5) Baseline LBO-CV & prediction maps
  #    - For each time slice (T14 day / T05 pre-dawn):
  #      1) Run leave-block-out CV on station data
  #      2) Produce truth vs predicted raster maps (model ensemble)
  #    - Diagnostics printed/printed:
  #      * Metrics table (RMSE/MAE/R2 per model)
  #      * Blocks plot (spatial CV blocks)
  #      * Pred-vs-obs scatter plot
  #      * Truth and prediction maps
  # ---------------------------------------------------------------------
  run14 <- run_for_time(st14, scen$R14, "T14", scen_local = scen, block_m = bs, models = mods)
  run05 <- run_for_time(st05, scen$R05, "T05", scen_local = scen, block_m = bs, models = mods)
  
  cat("\n== Metrics T14 ==\n"); print(run14$res$metrics)
  cat("\n== Metrics T05 ==\n"); print(run05$res$metrics)
  
  print(run14$res$blocks_plot); print(run14$res$diag_plot)
  print(run05$res$blocks_plot); print(run05$res$diag_plot)
  print(run14$maps$p_truth);    print(run14$maps$p_pred)
  print(run05$maps$p_truth);    print(run05$maps$p_pred)
  
  # =====================================================================
  # 6) SCALE → R* tuning → tuned CV + maps
  #
  # Pipeline rationale:
  #   (a) Variogram reveals correlation ranges in the point field.
  #   (b) U-curve scans DEM-smoothing radii ~ [L50, L95] to find data-driven R*.
  #       Each radius implies a different "macro-signal" (E*, slope*, cosi*).
  #       We refit CV at each R and pick the RMSE-minimizer (R*).
  #   (c) With R* in hand, we derive feature rasters at that scale: E*, slp*, cosi*.
  #   (d) Re-extract station features at R* to ensure training/prediction consistency.
  #   (e) Run LBO-CV again (tuned) using E* as the reference raster for blocks/domain.
  #   (f) Predict tuned maps by injecting the feature rasters.
  #   (g) Build compact multi-model panels with residual diagnostics.
  #
  # Performance tips if tuning feels slow:
  #   - Reduce n_grid in the U-curve (e.g., 5 instead of 9).
  #     n_grid sets how many candidate smoothing radii are tested in the U-curve search for R*
  #   - Trim `mods` to a smaller set while teaching the concept.
  #   - Increase block size slightly (fewer blocks → fewer CV folds).
  # =====================================================================
  
  # --- (a) Variogram → L50/L95 -------------------------------------------------
  Ls14 <- compute_Ls_from_points(st14, value_col = "temp")
  Ls05 <- compute_Ls_from_points(st05, value_col = "temp")
  
  p_vg14 <- plot_variogram_with_scales(Ls14$vg, Ls14$L50, Ls14$L95, Ls14$sill,
                                       "T14 — empirical variogram")
  p_vg05 <- plot_variogram_with_scales(Ls05$vg, Ls05$L50, Ls05$L95, Ls05$sill,
                                       "T05 — empirical variogram")
  print(p_vg14); print(p_vg05)
  
  # --- (b) U-curve → R* --------------------------------------------------------
  # We pass *explicit* sun angles so tune_Rstar_ucurve() can build cosi@R
  # consistently with the scenario's solar geometry.
  tune14 <- tune_Rstar_ucurve(
    stn_sf = st14,
    E      = scen$E,
    alt    = scen$sun$T14$alt,
    az     = scen$sun$T14$az,
    L50    = Ls14$L50,
    L95    = Ls14$L95,
    block_fallback = bs,
    n_grid = 6
  )
  
  tune05 <- tune_Rstar_ucurve(
    stn_sf = st05,
    E      = scen$E,
    alt    = scen$sun$T05$alt,
    az     = scen$sun$T05$az,
    L50    = Ls05$L50,
    L95    = Ls05$L95,
    block_fallback = bs,
    n_grid = 6
  )
  
  # Plot the U-curves and report chosen R* (rounded for readability).
  p_uc14 <- plot_ucurve(tune14$grid, tune14$R_star, "T14 — U-curve")
  p_uc05 <- plot_ucurve(tune05$grid, tune05$R_star, "T05 — U-curve")
  print(p_uc14); print(p_uc05)
  
  # IMPORTANT: use %.0f (not %d) because R* is numeric (may be non-integer).
  message(sprintf("Chosen R* — T14: %.0f m | blocks ≈ %d m", tune14$R_star, tune14$block_m))
  message(sprintf("Chosen R* — T05: %.0f m | blocks ≈ %d m", tune05$R_star, tune05$block_m))
  
  # --- (c) Feature rasters @R* -------------------------------------------------
  # Smooth DEM at R* and derive slope/incident-cosine given the scenario sun angles.
  fr14 <- smooth_dem_and_derive(
    scen$E, scen$sun$T14$alt, scen$sun$T14$az, radius_m = tune14$R_star
  )
  fr05 <- smooth_dem_and_derive(
    scen$E, scen$sun$T05$alt, scen$sun$T05$az, radius_m = tune05$R_star
  )
  
  # --- (d) Station features @R* ------------------------------------------------
  # Re-extract E*, slope*, cosi* (plus consistent LC factors) at station points.
  # This keeps training features aligned with the tuned raster features.
  st14_R <- add_drifts_at_R(
    st14, scen$E, scen$sun$T14$alt, scen$sun$T14$az, tune14$R_star,
    lc = scen$lc, lc_levels = scen$lc_levels,
    na_action = "fill"   # or "drop" if you prefer to omit affected stations
  )
  st05_R <- add_drifts_at_R(
    st05, scen$E, scen$sun$T05$alt, scen$sun$T05$az, tune05$R_star,
    lc = scen$lc, lc_levels = scen$lc_levels,
    na_action = "fill"   # or "drop" if you prefer to omit affected stations
  )
  
  # --- (e) LBO-CV @R* ----------------------------------------------------------
  # Use the tuned smoothed DEM (E*) as the reference for CV blocks and domain
  # geometry so the CV respects the working resolution/scale of the model.
  bench14 <- run_lbo_cv(st14_R, E = scen$E, block_size = bs, models = mods)
  bench05 <- run_lbo_cv(st05_R, E = scen$E, block_size = bs, models = mods)
  print(bench14$metrics); print(bench05$metrics)
  
  # --- (f) Tuned maps ----------------------------------------------------------
  # Inject the tuned feature rasters so model predictions operate at R* scale.
  maps14_tuned <- predict_maps(
    stn_sf = st14_R, truth_raster = scen$R14, which_time = "T14",
    scen = scen, models = mods, lc_levels = scen$lc_levels,
    feature_rasters = list(E = fr14$Es, slp = fr14$slp, cosi = fr14$cosi)
  )
  maps05_tuned <- predict_maps(
    stn_sf = st05_R, truth_raster = scen$R05, which_time = "T05",
    scen = scen, models = mods, lc_levels = scen$lc_levels,
    feature_rasters = list(E = fr05$Es, slp = fr05$slp, cosi = fr05$cosi)
  )
  
  # --- (g) Panels: truth | predictions | residual diagnostics ------------------
  panel_T14 <- build_panels_truth_preds_errors_paged(
    maps = maps14_tuned, truth_raster = scen$R14, cv_tbl = bench14$cv,
    which_time = "T14", models_per_page = 7, scatter_next_to_truth = TRUE
  )
  panel_T05 <- build_panels_truth_preds_errors_paged(
    maps = maps05_tuned, truth_raster = scen$R05, cv_tbl = bench05$cv,
    which_time = "T05", models_per_page = 7, scatter_next_to_truth = TRUE
  )
  print(panel_T14[[1]]); print(panel_T05[[1]])
  
  
  # Sensor noise (°C) – from specs or co-location
  # 1) Choose the model you want to base alpha on (same as your main map model)
  model_for_alpha <- "RF"   # e.g., "RF", "KED", "GAM" — adjust to taste
  
  # 2) Sensor noise (°C). You can pull this from scenario params if you prefer.
  sigma_inst <- 0.43  # example; or 0.3–0.5 from specs/colocation
  
  # 3) Derive alpha (microscale share) from CV residual variograms (T14/T05)
  #    (joins station coords by 'id', averages residuals per station, fits residual variogram)
  alpha14 <- alpha_from_cv(bench14$cv, st_sf = st14, model = model_for_alpha,
                           id_col_st = "id", id_col_cv = "id")
  alpha05 <- alpha_from_cv(bench05$cv, st_sf = st05, model = model_for_alpha,
                           id_col_st = "id", id_col_cv = "id")
  
  # 4) Fallbacks if fit was unstable
  if (!is.finite(alpha14)) alpha14 <- 0.6
  if (!is.finite(alpha05)) alpha05 <- 0.6
  
  # 5) Optional: cap sigma_inst to the RMSE scale to avoid zero remainder cases
  rmse14_tuned <- rmse_from_cv(bench14$cv, model_for_alpha)
  rmse05_tuned <- rmse_from_cv(bench05$cv, model_for_alpha)
  sigma_inst14 <- cap_sigma_inst(sigma_inst, rmse14_tuned, frac = 0.95)
  sigma_inst05 <- cap_sigma_inst(sigma_inst, rmse05_tuned, frac = 0.95)
  
  # 6) Compute budgets (base & tuned). 
  #    If your simple_error_budget() accepts a scalar sigma_inst, pass sigma_inst14/05 accordingly.
  eb14_base  <- simple_error_budget(run14$res, sigma_inst = sigma_inst14, alpha = alpha14)
  eb05_base  <- simple_error_budget(run05$res, sigma_inst = sigma_inst05, alpha = alpha05)
  eb14_tuned <- simple_error_budget(bench14,   sigma_inst = sigma_inst14, alpha = alpha14)
  eb05_tuned <- simple_error_budget(bench05,   sigma_inst = sigma_inst05, alpha = alpha05)
  
  # 7) (Optional) also report SDs for readability
  variance_components <- c("Total var","Instrument var","Microscale var","Mesoscale var")
  
  eb_all <- dplyr::bind_rows(
    dplyr::mutate(eb14_base,  Time="T14", Mode="Base"),
    dplyr::mutate(eb05_base,  Time="T05", Mode="Base"),
    dplyr::mutate(eb14_tuned, Time="T14", Mode="Tuned"),
    dplyr::mutate(eb05_tuned, Time="T05", Mode="Tuned")
  ) |>
    dplyr::relocate(Time, Mode) |>
    dplyr::group_by(Time, Mode) |>
    dplyr::mutate(TotalVar = Value[Component == "Total var"]) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      # SD only for variance rows (units: °C)
      Value_SD = dplyr::case_when(
        Component %in% variance_components ~ sqrt(Value),
        TRUE ~ NA_real_
      ),
      # Share of total only for variance rows
      ShareOfTotal = dplyr::case_when(
        Component %in% variance_components & TotalVar > 0 ~ Value / TotalVar,
        TRUE ~ NA_real_
      ),
      Share_pct = dplyr::case_when(
        is.finite(ShareOfTotal) ~ scales::percent(ShareOfTotal, accuracy = 0.1),
        TRUE ~ NA_character_
      )
    )
  
  print(eb_all)
  
  # Fehlerbudgets berechnen (Base = runXX$res, Tuned = benchXX)
  eb14_base  <- simple_error_budget(run14$res, sigma_inst, alpha14) |>
    dplyr::mutate(Time = "T14", Mode = "Base")
  eb05_base  <- simple_error_budget(run05$res, sigma_inst, alpha05) |>
    dplyr::mutate(Time = "T05", Mode = "Base")
  eb14_tuned <- simple_error_budget(bench14,   sigma_inst, alpha14) |>
    dplyr::mutate(Time = "T14", Mode = "Tuned")
  eb05_tuned <- simple_error_budget(bench05,   sigma_inst, alpha05) |>
    dplyr::mutate(Time = "T05", Mode = "Tuned")
  
  eb_all <- dplyr::bind_rows(eb14_base, eb05_base, eb14_tuned, eb05_tuned) |>
    dplyr::relocate(Time, Mode)
  
  print(eb_all)
  
  order_levels <- c("RMSE","Bias","Total var","Instrument var","Microscale var","Mesoscale var")
  eb_all <- eb_all |>
    dplyr::mutate(Component = factor(Component, levels = order_levels)) |>
    dplyr::arrange(Time, Mode, Component)
  # einfache Stacked-Bar-Plot-Funktion
  plot_error_budget <- function(df) {
    d <- df |>
      dplyr::filter(Component %in% c("Instrument var","Microscale var","Mesoscale var"))
    ggplot2::ggplot(d,
                    ggplot2::aes(x = interaction(Time, Mode, sep = " "), y = Value, fill = Component)
    ) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = "Variance (°C²)", title = "Error budget by time & mode")
  }
  p_eb <- plot_error_budget(eb_all)
  print(p_eb)
  
  # =====================================================================
  # 7) Optional: save everything at the end (plots + tables + rasters)
  #    - Change `export <- FALSE` at the top to only run/plot interactively
  #    - We wrap saves in try() so a single failed save does not abort the run.
  # =====================================================================
  if (export) {
    # ---------- Ausgabe-Verzeichnis: results_<scen-name> ----------
    out_root <- here::here("block4_5")
    out_dir  <- file.path(out_root, sprintf("results_%s", scen_name))
    fig_dir  <- file.path(out_dir, "fig")
    tab_dir  <- file.path(out_dir, "tab")
    ras_dir  <- file.path(out_dir, "ras")
    # ohne Rückfrage, rekursiv, ohne Warnungen
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(ras_dir, recursive = TRUE, showWarnings = FALSE)
    
  
    
    # ---------- Baseline: Previews & CV-Plots ----------------------
    save_plot_min(plot_landcover_terrain(scen, stations = st14, layout = "vertical"),
                  fn_fig("landcover_terrain"))
    save_plot_min(plot_block_overview_2x2_en(scen, pts_sf = st14), fn_fig("overview_2x2"))
    save_plot_min(run14$res$blocks_plot, fn_fig("T14_blocks"))
    save_plot_min(run14$res$diag_plot,   fn_fig("T14_diag"))
    save_plot_min(run05$res$blocks_plot, fn_fig("T05_blocks"))
    save_plot_min(run05$res$diag_plot,   fn_fig("T05_diag"))
    save_plot_min(run14$maps$p_truth,    fn_fig("T14_truth"))
    save_plot_min(run14$maps$p_pred,     fn_fig("T14_pred"))
    save_plot_min(run05$maps$p_truth,    fn_fig("T05_truth"))
    save_plot_min(run05$maps$p_pred,     fn_fig("T05_pred"))
    save_plot_min(p_block_T14,   fn_fig("T14_block_box"))
    save_plot_min(p_abs_T14,     fn_fig("T14_abs_error_box"))
    save_plot_min(p_scatter_T14, fn_fig("T14_obs_pred_scatter"))
    save_plot_min(p_dens_T14,    fn_fig("T14_residual_density"))
    
    save_plot_min(p_block_T05,   fn_fig("T05_block_box"))
    save_plot_min(p_abs_T05,     fn_fig("T05_abs_error_box"))
    save_plot_min(p_scatter_T05, fn_fig("T05_obs_pred_scatter"))
    save_plot_min(p_dens_T05,    fn_fig("T05_residual_density"))
    # --- Tuned Panels ---
    save_plot_min(panel_T14[[1]], fn_fig("T14_panel_tuned"))
    save_plot_min(panel_T05[[1]], fn_fig("T05_panel_tuned"))
    
    # --- Raster ---
    save_raster_min(scen$E,   fn_ras("E_dem"))
    save_raster_min(scen$R14, fn_ras("R14_truth"))
    save_raster_min(scen$R05, fn_ras("R05_truth"))
    if ("lc" %in% names(scen)) save_raster_min(scen$lc, fn_ras("landcover"))
    # ---------- Scale inference + tuned panels ----------------------
    safe_save_plot(p_vg14, fn_fig("T14_variogram_scalet"))
    safe_save_plot(p_vg05, fn_fig("T05_variogram_scalet"))
    safe_save_plot(p_uc14, fn_fig("T14_ucurve"))
    safe_save_plot(p_uc05, fn_fig("T05_ucurve"))
    safe_save_plot([[1]], fn_fig("T14_panel_tuned"))
    safe_save_plot(panel_T05[[1]], fn_fig("T05_panel_tuned"))
    
    save_table_readable(bench14$metrics, "metrics_T14_tuned", "Tuned metrics — T14")
    save_table_readable(bench05$metrics, "metrics_T05_tuned", "Tuned metrics — T05")
    save_table_readable(tune14$grid,     "Ucurve_T14",       "U-curve grid — T14")
    save_table_readable(tune05$grid,     "Ucurve_T05",       "U-curve grid — T05")
    save_table_readable(data.frame(L50 = Ls14$L50, L95 = Ls14$L95, R_star = tune14$R_star),
                        "scales_T14", "Scales — T14 (L50/L95/R*)")
    save_table_readable(data.frame(L50 = Ls05$L50, L95 = Ls05$L95, R_star = tune05$R_star),
                        "scales_T05", "Scales — T05 (L50/L95/R*)")
    save_table_readable(run14$res$metrics, file.path(tab_dir, sprintf("metrics_T14_%s", scen_name)))
    save_table_readable(run05$res$metrics, file.path(tab_dir, sprintf("metrics_T05_%s", scen_name)))
    save_table_readable(bench14$metrics,   file.path(tab_dir, sprintf("metrics_T14_tuned_%s", scen_name)))
    save_table_readable(bench05$metrics,   file.path(tab_dir, sprintf("metrics_T05_tuned_%s", scen_name)))
    save_table_readable(eb_all,            file.path(tab_dir, sprintf("error_budget_%s", scen_name)))
    #
    # ---------- Raster mit Szenario-Präfix --------------------------
    try(terra::writeRaster(scen$E,   fn_ras("E_dem")),     silent = TRUE)
    try(terra::writeRaster(scen$R14, fn_ras("R14_truth")), silent = TRUE)
    try(terra::writeRaster(scen$R05, fn_ras("R05_truth")), silent = TRUE)
    if ("lc" %in% names(scen))
      try(terra::writeRaster(scen$lc, fn_ras("landcover")), silent = TRUE)
    
    # ---------- Sessioninfo -----------------------------------------
    try(saveRDS(sessionInfo(), file.path(out_dir, sprintf("%s_sessionInfo.rds", scen_name))),
        silent = TRUE)
    
    message("✔ Exports written to: ", normalizePath(out_dir, winslash = "/"))
  }
  
