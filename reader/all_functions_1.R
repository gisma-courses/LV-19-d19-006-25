  # -----------------------------------------------------------------------------
  # Microclimate Spatial Toolkit — Sorted & Commented (R)
  # -----------------------------------------------------------------------------
  # Purpose
  #   Cleaned, deduplicated, and logically ordered set of functions for:
  #     1) Ecowitt data ingestion and name hygiene
  #     2) Spatial interpolation methods (OK, IDW, TPS, GAM, RF, Voronoi, KED)
  #     3) Variogram-based scale inference (L50/L95) and drift radius selection (R)
  #     4) Scale-matched predictor construction (DEM features + generic rasters)
  #     5) Spatial CV blocks, R* tuning via U-curve, method benchmarking
  #     6) Diagnostics: error budget, multi-panel raster plotting
  #
  # Notes
  #   • All comments are in English. Internal helpers are prefixed with a dot.
  #   • Duplicates removed (e.g., get_Ls, make_block_folds, .safe_numeric_df, .square_raster).
  #   • Functions are grouped by topic and ready to split into separate files.
  #   • Minimal changes to your original logic; only consolidation and clarity edits.
  #
  # Suggested file layout (place under R/ in an R package or scripts/ in a project):
  #   00_utils_io.R               — Name cleaning, time formatting, registries
  #   01_ecowitt_ingest.R         — Ecowitt Excel parsing + logger merge
  #   10_spatial_utils.R          — Grid alignment, safe coercions, formulas, metrics
  #   20_interpolation_methods.R  — Voronoi, IDW, OK, TPS, GAM, RF, KED, single-step kriging
  #   30_variogram_scale.R        — estimate_L_all, get_Ls, choose_R_from_L, alpha
  #   40_predictors.R             — Smoothing kernels, DEM features, generic smoothing
  #   50_workflow.R               — End-to-end A→D workflow, summaries, default tables
  #   60_blockcv_tuning.R         — Block folds, R tuning via block-CV (U-curve)
  #   70_benchmark.R              — Benchmark at R* with spatial blocks
  #   80_panel_plot.R             — Faceted time-series raster plotting
  #   90_error_budget.R           — Error budget decomposition
  # -----------------------------------------------------------------------------
  
  # ============================= 00_utils_io.R ================================ #
  # --- small helpers for the viewer --------------------------------------------
  # Robust: nimmt bei extract() immer die erste Werte-Spalte (ohne ID/cell)
  .safe_extract_vals <- function(r, xy) {
    stopifnot(inherits(r, "SpatRaster"))
    df <- terra::extract(r, xy)        # kein ID-Argument verwenden
    if (!is.data.frame(df)) return(as.numeric(df))
    cn <- names(df)
    
    # Kandidaten: Spalt(en), die zum Layernamen passen
    vcol <- which(cn %in% names(r))
    # ID-ähnliche Spalten (je nach terra-Version: "ID", "cell")
    id_like <- which(tolower(cn) %in% c("id", "cell", "cells"))
    
    pick <- NA_integer_
    if (length(vcol) >= 1) {
      pick <- vcol[1]
    } else {
      # nimm die erste Spalte, die NICHT wie ID/cell aussieht
      pick <- setdiff(seq_along(cn), id_like)[1]
      if (is.na(pick)) pick <- 1L
    }
    as.numeric(df[[pick]])
  }
  
  # ---- CRS/Palette-Helper ------------------------------------------------------
  # Pick a reference raster (meters, square cells) for grid/blocks/smoothing
  .dem_ref <- function(dem, extra_preds = NULL) {
    if (inherits(dem, "SpatRaster")) return(.square_raster(dem))
    if (inherits(extra_preds, "SpatRaster")) return(.square_raster(extra_preds))
    if (is.list(extra_preds) && length(extra_preds)) return(.square_raster(extra_preds[[1]]))
    stop("Provide 'dem' or at least one raster in 'extra_preds'.")
  }
  
  
  # Project a SpatRaster to WGS84 if needed (for leaflet)
  .to_ll_sf <- function(x) {
    stopifnot(inherits(x, "sf"))
    # drop Z/M if present (Leaflet wants XY)
    x <- try(suppressWarnings(sf::st_zm(x, drop = TRUE, what = "ZM")), silent = TRUE)
    if (inherits(x, "try-error")) x <- x[[1]]  # if st_zm returned list
    if (is.na(sf::st_crs(x))) stop("sf has no CRS set — set it before calling the viewer.")
    sf::st_transform(x, 4326)
  }
  
  .to_ll_rast <- function(r) {
    stopifnot(inherits(r, "SpatRaster"))
    terra::project(r, "EPSG:4326")
  }
  
  
  # Plot-Polygon als Linien (vermeidet centroid-Warnungen komplett)
  .boundary_lines_ll <- function(poly_sf) {
    poly_ll <- .to_ll_sf(poly_sf)  # not .as_ll()
    sf::st_cast(sf::st_boundary(poly_ll), "MULTILINESTRING")
  }
  
  # Farbskala mit leichtem Puffer → keine "outside color scale"-Warnung beim Reproject
  # Palette + padded domain for stable legends
  .pal_with_domain <- function(r_ll, pad = 0.02) {
    v <- terra::values(r_ll)[, 1]
    v <- v[is.finite(v)]
    if (!length(v)) {
      dom <- c(0, 1)
    } else {
      rng  <- range(v)
      span <- diff(rng); if (!is.finite(span) || span == 0) span <- 1
      dom  <- seq(rng[1] - pad * span, rng[2] + pad * span, length.out = 256)
    }
    list(
      pal = leaflet::colorNumeric(viridisLite::viridis(length(dom)), domain = dom, na.color = "transparent"),
      dom = dom
    )
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  .pts_df_from_var <- function(sf_table, varname) {
    stopifnot(inherits(sf_table, "sf"))
    idx <- !is.na(sf_table[[varname]])
    xy  <- sf::st_coordinates(sf_table[idx, ])
    data.frame(x = xy[,1], y = xy[,2], altitude = sf_table$altitude[idx], value = sf_table[[varname]][idx])
  }
  # TODO: move to all_in.R — safe evaluator that returns NULL on error
  .safe <- function(expr, tag = NULL) {
    tryCatch(expr, error = function(e) { if (!is.null(tag)) message(sprintf("skip[%s]: %s", tag, conditionMessage(e))); NULL })
  }
  #' Clean column names in data frames
  #' Removes unit tokens, normalizes spaces and underscores.
  #' @param df data.frame or tibble
  #' @return same object with cleaned names
  clean_names <- function(df) {
    df |>
      dplyr::rename_with(~ stringr::str_replace_all(.x, "\\(％\\)|\\(℃\\)|\\(\\%\\)", "")) |>
      dplyr::rename_with(~ stringr::str_replace_all(.x, "\\s+", "_")) |>
      dplyr::rename_with(~ stringr::str_replace_all(.x, "_+", "_")) |>
      dplyr::rename_with(~ stringr::str_remove_all(.x, "_$"))
  }
  
  #' Clean column names in sf objects (defensive)
  #' @param sf_obj sf object
  #' @return sf object with cleaned names
  clean_names_sf <- function(sf_obj) {
    nm <- names(sf_obj)
    nm <- nm |>
      stringr::str_replace_all("\\((?:％|%|℃|°C)\\)", "") |>
      stringr::str_replace_all("[:\\s]+", "_") |>
      stringr::str_replace_all("[/\\\\]", "-") |>
      stringr::str_replace_all("[^0-9A-Za-z_.-]", "_") |>
      stringr::str_replace_all("_+", "_") |>
      stringr::str_replace_all("^_+|_+$", "")
    names(sf_obj) <- nm
    sf_obj
  }
  
  #' Fix variable names for R compatibility
  #' @param names_vec character vector
  #' @return cleaned, unique names
  fix_names <- function(names_vec) {
    names_vec <- gsub("[-: ]", "", names_vec)
    names_vec <- gsub("^([0-9])", "A\\1", names_vec)
    make.names(names_vec, unique = TRUE)
  }
  
  #' Format encoded timestamps to readable labels
  #' Accepts "AYYYYMMDDHHMMSS" → "YYYY-MM-DD HH:MM" or "AYYYYMMDD(_D)" → "YYYY-MM-DD".
  #' @param x character scalar
  #' @return formatted character
  pretty_time <- function(x) {
    vapply(as.character(x), function(s) {
      if (grepl("^A\\d{14}$", s)) {
        ts <- as.POSIXct(substr(s, 2, 15), format = "%Y%m%d%H%M%S", tz = "UTC")
        format(ts, "%Y-%m-%d %H:%M")
      } else if (grepl("^A\\d{8}(_D)?$", s)) {
        ts <- as.Date(substr(s, 2, 9), format = "%Y%m%d")
        format(ts, "%Y-%m-%d")
      } else {
        s
      }
    }, character(1))
  }
  
  
  #' Minimal ID cleaning (units/spaces → underscore)
  #' @param x character vector
  #' @return character
  clean_ids <- function(x) {
    x <- gsub("\\(℃\\)|\\(％\\)|\\(\\%\\)", "", x)
    x <- gsub("\\s+", "_", x)
    x
  }
  
  #' Method registry used by benchmark_at_Rstar()
  #' Maps human-friendly names to function symbols (strings) looked up via get().
  reg <- list(
    Voronoi = "voro_simple",
    IDW     = "pred_idw",
    OK      = "pred_ok",
    GAM     = "pred_gam",
    RF      = "pred_rf",
    KED     = "pred_ked",
    Trend   = "pred_trend_lm"   # <- neu
  )
  
  
  
  # ============================ 01_ecowitt_ingest.R =========================== #
  
  #' Extract core variables from an Ecowitt Excel export (two header rows)
  #' @param fn path to Ecowitt Excel export
  #' @return data.frame with Time + selected sensors and outdoor variables
  extract_ecowitt_core_vars <- function(fn) {
    library(readxl); library(zoo); library(stringr); library(dplyr)
    headers <- readxl::read_excel(fn, n_max = 2, col_names = FALSE)
    headers[1, 1] <- "Time"
    filled_row <- zoo::na.locf(as.character(headers[1, ]), fromLast = FALSE)
    headers[1, ] <- as.list(filled_row)
    raw_names <- paste0(headers[1, ], "_", headers[2, ])
    
    channels <- paste0("ch", 1:8)
    final_names <- sapply(raw_names, function(n) {
      ch_match <- stringr::str_extract(n, paste0(channels, collapse = "|"))
      group <- stringr::str_extract(n, "^[^_]+")
      meas <- stringr::str_replace_all(n, paste0(group, "_"), "")
      meas <- stringr::str_replace_all(meas, paste0(channels, collapse = "|"), "")
      meas <- stringr::str_trim(meas)
      if (!is.na(ch_match) && tolower(group) %in% c("temperature", "humidity", "soilmoisture")) {
        paste0(group, "_", ch_match, ifelse(meas == "", "", paste0("_", meas)))
      } else if (tolower(group) %in% c("solar", "uvi", "rainfall", "pressure", "wind", "outdoor")) {
        paste0(group, "_", meas)
      } else n
    })
    final_names[1] <- "Time"
    
    data <- readxl::read_excel(fn, skip = 2, col_names = FALSE)
    colnames(data) <- final_names
    data <- dplyr::mutate(data, Time = as.POSIXct(Time))
    
    core_pattern <- "ch[1-8]"
    core_types   <- c("temperature", "humidity", "soilmoisture")
    extra_groups <- c("outdoor", "solar", "uvi", "rainfall", "pressure", "wind")
    rainfall_keep <- "Rainfall_Rain_Rate"
    
    keep_cols <- names(data)[
      (grepl(core_pattern, names(data), ignore.case = TRUE) &
         grepl(paste(core_types, collapse = "|"), names(data), ignore.case = TRUE)) |
        grepl(paste(extra_groups, collapse = "|"), names(data), ignore.case = TRUE) |
        names(data) == "Time"
    ]
    data_filtered <- data[, keep_cols]
    
    drop_pattern <- paste(
      "Low\\(.*\\)",      # Low(...)
      "High\\(.*\\)",     # High(...)
      "Battery_",           # battery voltage
      "Feels\\s*Like\\(.*\\)",
      "Rainfall_",          # drop all Rainfall_* (we will keep one below)
      sep = "|"
    )
    to_drop <- grepl(drop_pattern, names(data_filtered), perl = TRUE)
    to_drop[names(data_filtered) == rainfall_keep] <- FALSE
    data_filtered <- data_filtered[, !to_drop]
    data_filtered
  }
  
  #' Merge Temperature/Humidity from two Ecowitt loggers
  #' @param df1,df2 two cleaned data.frames
  #' @param id1,id2 character tags (default "FC29","DB2F")
  #' @return list(temperature=..., humidity=...)
  merge_ecowitt_logger_vars <- function(df1, df2, id1 = "FC29", id2 = "DB2F") {
    library(dplyr); library(stringr)
    df1 <- dplyr::mutate(df1, Time = as.POSIXct(Time))
    df2 <- dplyr::mutate(df2, Time = as.POSIXct(Time))
    
    select_and_tag <- function(df, pattern, id) {
      df |>
        dplyr::select(Time, dplyr::contains(pattern)) |>
        dplyr::rename_with(~ paste0(., "_", id), -Time)
    }
    temp_df1 <- select_and_tag(df1, "Temperature(℃)", id1)
    temp_df2 <- select_and_tag(df2, "Temperature(℃)", id2)
    hum_df1  <- select_and_tag(df1, "Humidity(%)", id1)
    hum_df2  <- select_and_tag(df2, "Humidity(%)", id2)
    
    merge_clean <- function(df_a, df_b) {
      dplyr::full_join(df_a, df_b, by = "Time") |>
        dplyr::arrange(Time) |>
        dplyr::rename_with(~ stringr::str_remove(., "^Temp and Humidity "), -Time) |>
        dplyr::mutate(dplyr::across(where(is.character) & !dplyr::matches("Time"), ~ dplyr::na_if(., "-") |> as.numeric()))
    }
    list(temperature = merge_clean(temp_df1, temp_df2),
         humidity    = merge_clean(hum_df1, hum_df2))
  }
  
  #' Map short station IDs to verbose names
  #' @param id e.g., "base", "ch1", "ch1_r"
  #' @param measure "Temperature" or "Humidity"
  #' @param local_tag,remote_tag character
  #' @return character
  to_verbose <- function(id, measure = c("Temperature","Humidity"),
                         local_tag = "FC29", remote_tag = "DB2F") {
    measure <- match.arg(measure)
    id_low  <- tolower(id)
    is_remote <- stringr::str_ends(id_low, "_r")
    ch <- stringr::str_match(id_low, "^ch(\\d+)(?:_r)?$")[,2]
    dplyr::case_when(
      id_low == "base"    ~ paste0("Outdoor_", measure, "_", local_tag),
      id_low == "base_r"  ~ paste0("Outdoor_", measure, "_", remote_tag),
      !is.na(ch) & !is_remote ~ paste0("CH", ch, "_", measure, "_", local_tag),
      !is.na(ch) &  is_remote ~ paste0("CH", ch, "_", measure, "_", remote_tag),
      TRUE ~ id
    )
  }
  
  
  # ============================ 10_spatial_utils.R ============================ #
  
  #' Ensure cell size is square; resample if needed (bilinear)
  .square_raster <- function(r) {
    rs <- terra::res(r)
    if (abs(rs[1]-rs[2]) < 1e-6) return(r)
    r0 <- terra::rast(ext = terra::ext(r), resolution = mean(rs), crs = terra::crs(r))
    terra::resample(r, r0, method = "bilinear")
  }
  
  #' Align raster r to reference ref (CRS + grid), then mask to ref
  .align_to <- function(r, ref) {
    r2 <- .square_raster(r)
    if (!terra::compareGeom(r2, ref, stopOnError = FALSE)) {
      r2 <- terra::project(r2, terra::crs(ref))
      r2 <- terra::resample(r2, ref, method = "bilinear")
    }
    terra::mask(r2, ref[[1]], maskvalues = NA)
  }
  
  #' Safe numeric coercion for data.frames (drop ID col, remove all-NA cols)
  .safe_numeric_df <- function(x) {
    if (is.data.frame(x) && ncol(x) > 0 && names(x)[1] %in% c("ID","id"))
      x <- x[, -1, drop = FALSE]
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x[] <- lapply(x, function(col) {
      if (is.list(col)) col <- unlist(col, use.names = FALSE)
      suppressWarnings(as.numeric(col))
    })
    bad <- vapply(x, function(z) all(!is.finite(z)), logical(1))
    if (any(bad) && length(bad) < ncol(x)) x <- x[, !bad, drop = FALSE]
    x
  }
  
  #' Build model formula from a table
  .make_formula <- function(df, response = "value", drop_cols = c("fold_id","cell","x","y")) {
    xs <- setdiff(names(df), c(response, drop_cols))
    as.formula(paste(response, "~", paste(xs, collapse = " + ")))
  }
  
  #' Fit and predict (LM or GAM) with robust fallback
  .fit_and_predict <- function(train, test, response = "value", method = c("gam","lm")) {
    method <- match.arg(method)
    fml <- .make_formula(train, response = response)
    if (method == "gam" && requireNamespace("mgcv", quietly = TRUE)) {
      mod <- mgcv::gam(fml, data = train, method = "REML")
      preds <- as.numeric(predict(mod, newdata = test, type = "response"))
    } else {
      mod <- stats::lm(fml, data = train)
      preds <- as.numeric(predict(mod, newdata = test))
    }
    preds
  }
  
  #' Basic metrics
  .metrics <- function(obs, pred) {
    ok <- is.finite(obs) & is.finite(pred)
    if (!any(ok)) {
      return(data.frame(RMSE = NA_real_, MAE = NA_real_, Bias = NA_real_))
    }
    e <- pred[ok] - obs[ok]
    if (!length(e) || all(!is.finite(e))) {
      return(data.frame(RMSE = NA_real_, MAE = NA_real_, Bias = NA_real_))
    }
    data.frame(
      RMSE = sqrt(mean(e^2)),
      MAE  = mean(abs(e)),
      Bias = mean(e)
    )
  }
  
  
  
  # ======================== 20_interpolation_methods.R ======================== #
  
  #' Voronoi (Thiessen) interpolation on a raster grid
  #' @param m_sf sf points with the column `varname`
  #' @param varname name of the value column to rasterize via Voronoi
  #' @param DTM SpatRaster template (extent/CRS)
  #' @return SpatRaster with one layer "voronoi"
  voro_simple <- function(m_sf, varname, DTM) {
    stopifnot(inherits(m_sf, "sf"), inherits(DTM, "SpatRaster"))
    if (!varname %in% names(m_sf)) stop("Column '", varname, "' not found.")
    m_slice <- m_sf[!is.na(m_sf[[varname]]), c(varname, attr(m_sf, "sf_column"))]
    if (nrow(m_slice) < 2) stop("Need at least 2 points for Voronoi.")
    p <- terra::vect(sf::st_transform(m_slice, terra::crs(DTM)))
    voro <- terra::voronoi(p)
    voro <- terra::crop(voro, terra::as.polygons(terra::ext(DTM)))
    r <- terra::rasterize(voro, DTM, field = varname)
    names(r) <- "voronoi"
    r
  }
  
  #' Ordinary Kriging (OK) prediction using automap variogram
  pred_ok <- function(pts_df, DTM) {
    sp::coordinates(pts_df) <- ~ x + y
    sp::proj4string(pts_df) <- terra::crs(DTM, proj = TRUE)
    vgm_mod <- automap::autofitVariogram(value ~ 1, pts_df)$var_model
    newgrid <- stars::st_as_stars(DTM)
    kr <- gstat::krige(value ~ 1, locations = pts_df, newdata = newgrid, model = vgm_mod)
    r <- terra::rast(kr["var1.pred"]); names(r) <- "ok"; r
  }
  
  #' Inverse Distance Weighting (IDW)
  pred_idw <- function(pts_df, DTM, idp = 2) {
    sp::coordinates(pts_df) <- ~ x + y
    sp::proj4string(pts_df) <- terra::crs(DTM, proj = TRUE)
    newgrid <- stars::st_as_stars(DTM)
    idw <- gstat::idw(value ~ 1, locations = pts_df, newdata = newgrid, idp = idp)
    r <- terra::rast(idw["var1.pred"]); names(r) <- "idw"; r
  }
  
  #' Thin-Plate Spline (fields::Tps) with grid backfill
  .pred_to_rast <- function(pred_vec, template_rast) {
    r  <- template_rast[[1]]
    pv <- as.numeric(pred_vec)
    if (length(pv) != terra::ncell(r)) {
      stop("Prediction length (", length(pv), ") != ncell(template) (", terra::ncell(r), ").")
    }
    terra::values(r) <- pv
    names(r) <- "pred"
    r
  }
  
  # --- Simple trend predictor (LM) ---------------------------------------------
  # passt auf dein vorhandenes pts_df-Schema (x, y, altitude, value)
  pred_trend_lm <- function(pts_df, DTM, use_xy = TRUE, use_alt = TRUE) {
    stopifnot(inherits(DTM, "SpatRaster"))
    if (!all(c("x","y","value") %in% names(pts_df))) {
      stop("pred_trend_lm: pts_df braucht Spalten x, y, value (und optional altitude).")
    }
    df <- pts_df
    if (use_alt && !("altitude" %in% names(df))) use_alt <- FALSE
    
    # Formel zusammenbauen
    rhs <- c()
    if (use_xy)  rhs <- c(rhs, "x", "y")
    if (use_alt) rhs <- c(rhs, "altitude")
    if (!length(rhs)) rhs <- "1"  # nur Intercept
    fml <- stats::as.formula(paste("value ~", paste(rhs, collapse = " + ")))
    
    # Fit (LM) – robust gegen kleine n und ohne externe Abhängigkeiten
    fit <- stats::lm(fml, data = df)
    
    # Vorhersage auf Rastergrid
    grid_df <- terra::as.data.frame(DTM, xy = TRUE, cells = TRUE, na.rm = FALSE)
    names(grid_df) <- c("cell","x","y", "altitude")  # deine DTM-Layer heißt "altitude"
    pred <- as.numeric(predict(fit, newdata = grid_df))
    
    r <- .pred_to_rast(pred, DTM)
    names(r) <- "trend"
    r
  }
  
  pred_tps <- function(pts_df, DTM) {
    fit  <- fields::Tps(x = as.matrix(pts_df[, c("x","y")]), Y = pts_df$value)
    grid_df <- terra::as.data.frame(DTM, xy = TRUE, cells = TRUE, na.rm = FALSE)
    names(grid_df) <- c("cell","x","y","altitude")
    pred <- predict(fit, newdata = as.matrix(grid_df[, c("x","y")]))
    r <- .pred_to_rast(pred, DTM); names(r) <- "tps"; r
  }
  
  #' GAM with s(x,y) + s(altitude); robust fallbacks for small n
  pred_gam <- function(pts_df, DTM) {
    # 1) Daten vorbereiten
    pts <- pts_df[stats::complete.cases(pts_df[, c("x","y","altitude","value")]), ]
    if (!nrow(pts)) stop("No valid points for GAM.")
    pts <- stats::aggregate(value ~ x + y + altitude, data = pts, FUN = mean)
    
    n    <- nrow(pts)
    nxy  <- nrow(unique(pts[c("x","y")]))
    nalt <- length(unique(pts$altitude))
    
    # 2) Kleines-n-Fallback: einfache Höhen-Regression
    if (nxy < 4 || n < 6) {
      fit <- stats::lm(value ~ altitude, data = pts)
      grid_df <- terra::as.data.frame(DTM[[1]], xy = TRUE, cells = TRUE, na.rm = FALSE)
      names(grid_df)[1:4] <- c("cell","x","y","altitude")
      pred <- as.numeric(predict(fit, newdata = grid_df))
      r <- .pred_to_rast(pred, DTM); names(r) <- "gam_fallback_lm"; return(r)
    }
    
    # 3) GAM fit (einmal!), mit defensiven k-Werten
    k_xy  <- max(4, min(20, nxy - 1, floor(n/2)))
    k_alt <- max(3, min(8,  nalt - 1, floor(n/3)))
    
    fit <- tryCatch(
      mgcv::gam(
        value ~ s(x, y, k = k_xy,  bs = "tp") +
          s(altitude, k = k_alt, bs = "tp"),
        data = pts, method = "REML"
      ),
      error = function(e) {
        # robuster Fallback: nur Höhen-Spline
        mgcv::gam(value ~ s(altitude, k = k_alt, bs = "tp"),
                  data = pts, method = "REML")
      }
    )
    
    # 4) Vorhersage aufs Raster
    grid_df <- terra::as.data.frame(DTM[[1]], xy = TRUE, cells = TRUE, na.rm = FALSE)
    names(grid_df)[1:4] <- c("cell","x","y","altitude")
    pred <- as.numeric(predict(fit, newdata = grid_df, type = "response"))
    
    r <- .pred_to_rast(pred, DTM)
    names(r) <- "gam"
    r
  }

    #' Random Forest (x,y,altitude → value)
  pred_rf <- function(pts_df, DTM) {
    fit <- randomForest::randomForest(value ~ x + y + altitude, data = pts_df, ntree = 500)
    grid_df <- terra::as.data.frame(DTM, xy = TRUE, cells = TRUE, na.rm = FALSE)
    names(grid_df) <- c("cell","x","y","altitude")
    pred <- predict(fit, newdata = grid_df)
    r <- .pred_to_rast(as.numeric(pred), DTM); names(r) <- "rf"; r
  }
  
  #' Kriging with External Drift (KED) using altitude as drift
  # KED mit externer Drift aus einem Raster (z.B. DEM geglättet auf R_local)
  pred_ked <- function(pts_df, DTM_or_drift) {
    stopifnot(all(c("x","y","value") %in% names(pts_df)))
    stopifnot(inherits(DTM_or_drift, "SpatRaster"))
    
    drift <- DTM_or_drift[[1]]
    names(drift) <- "altitude"
    
    xy <- as.matrix(pts_df[, c("x","y")])
    pts_df$altitude <- .safe_extract_vals(drift, xy)
    
    ok <- is.finite(pts_df$value) & is.finite(pts_df$altitude)
    pts_df <- pts_df[ok, , drop = FALSE]
    if (nrow(pts_df) < 5) stop("pred_ked: <5 gültige Punkte nach Drift-Extraktion.")
    if (var(pts_df$altitude) < 1e-9) stop("pred_ked: Drift-Kontrast ~0 (fast konstant).")
    
    sp::coordinates(pts_df) <- ~ x + y
    sp::proj4string(pts_df) <- terra::crs(drift, proj = TRUE)
    
    vgm_mod <- automap::autofitVariogram(value ~ altitude, pts_df)$var_model
    
    nd <- stars::st_as_stars(drift); names(nd) <- "altitude"
    kr <- gstat::krige(value ~ altitude, locations = pts_df, newdata = nd, model = vgm_mod)
    
    r <- terra::rast(kr["var1.pred"]); names(r) <- "ked"; r
  }
  
  
  
  
  #' Single-time-step kriging convenience (universal with altitude drift if possible)
  #' @param varname column in data_sf to interpolate
  #' @param data_sf sf points with `altitude` and the `varname`
  #' @param dem_raster SpatRaster (prediction grid)
  #' @param output_dir output folder for GeoTIFF
  #' @param label 'raw' or 'pretty' (pretty formats layer name via pretty_time())
  interpolate_kriging <- function(varname, data_sf, dem_raster, output_dir = out_dir,
                                  label = c("raw","pretty")) {
    label <- match.arg(label)
    message("Interpolating: ", varname)
    if (!(varname %in% names(data_sf))) stop("Variable ", varname, " not found in data.")
    if (all(is.na(data_sf[[varname]]))) { warning("All values are NA for ", varname, " – skipping."); return(NULL) }
    f_drift <- stats::as.formula(paste(varname, "~ altitude"))
    vgm_model <- tryCatch(automap::autofitVariogram(f_drift, input_data = data_sf), error = function(e) NULL)
    if (is.null(vgm_model)) {
      vgm_model <- tryCatch({ list(var_model = automap::autofitVariogram(stats::as.formula(paste(varname, "~ 1")), input_data = data_sf)$var_model) }, error = function(e) NULL)
    }
    if (is.null(vgm_model)) { warning("Variogram failed for ", varname); return(NULL) }
    kriged_result <- tryCatch({
      gstat::krige(
        formula  = if (!is.null(vgm_model) && grepl("altitude", deparse(f_drift))) f_drift else stats::as.formula(paste(varname, "~ 1")),
        locations = data_sf,
        newdata   = stars::st_as_stars(dem_raster),
        model     = vgm_model$var_model
      )
    }, error = function(e) { warning("Kriging failed for ", varname, ": ", e$message); return(NULL) })
    if (is.null(kriged_result)) return(NULL)
    if ("var1.pred" %in% names(kriged_result)) kriged_result <- kriged_result["var1.pred"]
    new_name <- if (label == "pretty") pretty_time(varname) else varname
    names(kriged_result) <- new_name
    out_file <- file.path(output_dir, paste0(varname, "_interpolated.tif"))
    stars::write_stars(kriged_result, out_file, overwrite = TRUE)
    message("✔ Written: ", out_file)
    kriged_result
  }
  
  
  # ========================= 30_variogram_scale.R ============================ #
  
  #' Estimate practical correlation length(s) L from point data
  #' Returns model-based and empirical L_p, fitted vario, etc.
  estimate_L_all <- function(sf_pts,
                             p = c(0.95, 0.5),
                             crs_proj = 25832,
                             cutoff_frac = 0.7,
                             n_lags = 12,
                             plot = FALSE,
                             plot_type = c("base","ggplot2")) {
    plot_type <- match.arg(plot_type)
    stopifnot(inherits(sf_pts, "sf"))
    if (!"value" %in% names(sf_pts)) stop("Column 'value' not found in 'sf_pts'.")
    pts <- sf_pts[!is.na(sf_pts$value), ]
    if (nrow(pts) < 5) stop("Too few points for variogram fit (n < 5).")
    if (sf::st_is_longlat(pts)) pts <- sf::st_transform(pts, crs_proj)
    
    xy <- sf::st_coordinates(pts)
    dx <- diff(range(xy[,1])); dy <- diff(range(xy[,2]))
    maxd   <- sqrt(dx^2 + dy^2)
    cutoff <- cutoff_frac * maxd
    width  <- cutoff / n_lags
    
    vauto <- automap::autofitVariogram(value ~ 1, pts, cutoff = cutoff, width = width)
    vf <- vauto$var_model; vg <- vauto$exp_var
    j <- match(TRUE, as.character(vf$model) != "Nug", nomatch = NA_integer_)
    if (is.na(j)) stop("No structural variogram component found (only Nugget?).")
    model_name <- as.character(vf$model[j]); a <- vf$range[j]; kap <- vf$kappa[j]
    
    L_exp <- function(pp) -a * log(1 - pp)
    L_gau <- function(pp)  a * sqrt(-log(1 - pp))
    L_sph <- function(pp) { f <- function(x) 1.5*x - 0.5*x^3 - pp; if (pp >= 1) return(a); out <- try(uniroot(f, c(0, 1)), silent = TRUE); if (inherits(out, "try-error")) return(NA_real_); a * out$root }
    L_mat95 <- if (is.finite(kap) && kap > 0) a * sqrt(8 * kap) else 3 * a
    scale_matern <- function(pp) sqrt(-log(1 - pp) / -log(1 - 0.95))
    L_mat <- function(pp) L_mat95 * scale_matern(pp)
    
    L_model <- sapply(p, function(pi) {
      switch(model_name,
             "Exp" = L_exp(pi),
             "Gau" = L_gau(pi),
             "Sph" = L_sph(pi),
             "Mat" = L_mat(pi),
             L_exp(pi))
    })
    names(L_model) <- paste0("p", round(p*100))
    
    gmax <- max(vg$gamma, na.rm = TRUE)
    L_emp <- sapply(p, function(pi) {
      thr <- pi * gmax
      idx <- which(vg$gamma >= thr)
      if (length(idx) == 0) { i2 <- which.min(abs(vg$gamma - thr)); vg$dist[i2] } else vg$dist[min(idx)]
    })
    names(L_emp) <- paste0("p", round(p*100))
    
    explanation <- data.frame(
      type    = rep(c("model_based", "empirical"), each = length(p)),
      p       = rep(p, times = 2),
      L       = c(as.numeric(L_model), as.numeric(L_emp)),
      meaning = c(
        paste0("Model-based practical range: distance where variogram of model '", model_name, "' reaches ", round(100*p), "% of sill."),
        paste0("Empirical practical range: distance where experimental variogram reaches ", round(100*p), "% of max." )
      ),
      row.names = NULL
    )
    
    if (isTRUE(plot)) {
      if (plot_type == "base") {
        plot(vg$dist, vg$gamma, xlab = "Distance", ylab = "Semivariance",
             main = sprintf("Variogram (%s) with L markers", model_name), pch = 16)
        vline <- gstat::variogramLine(vf, maxdist = max(vg$dist, na.rm = TRUE), n = 200)
        lines(vline$dist, vline$gamma)
        if (length(L_model)) abline(v = as.numeric(L_model), lty = 2)
        if (length(L_emp))   abline(v = as.numeric(L_emp),   lty = 3)
        legend("bottomright",
               legend = c(paste0("model L_", names(L_model)), paste0("emp   L_", names(L_emp))),
               lty = c(rep(2, length(L_model)), rep(3, length(L_emp))), bty = "n", cex = 0.9)
      } else {
        if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_type='ggplot2' requires 'ggplot2'.")
        df_vgm <- vg
        vline <- gstat::variogramLine(vf, maxdist = max(df_vgm$dist, na.rm = TRUE), n = 200)
        df_lines <- data.frame(
          L     = c(as.numeric(L_model), as.numeric(L_emp)),
          label = c(paste0("model L_", names(L_model)), paste0("emp   L_", names(L_emp))),
          group = c(rep("model", length(L_model)), rep("emp", length(L_emp)))
        )
        p_g <- ggplot2::ggplot(df_vgm, ggplot2::aes(x = dist, y = gamma)) +
          ggplot2::geom_point() +
          ggplot2::geom_line(data = vline, ggplot2::aes(x = dist, y = gamma)) +
          ggplot2::geom_vline(data = df_lines, ggplot2::aes(xintercept = L, color = label, linetype = group)) +
          ggplot2::scale_linetype_manual(values = c(model = 2, emp = 3)) +
          ggplot2::labs(title = sprintf("Variogram (%s) with L markers", model_name),
                        x = "Distance", y = "Semivariance", color = "L values", linetype = "Type") +
          ggplot2::theme_minimal()
        print(p_g)
      }
    }
    
    list(L_model = L_model, L_emp = L_emp, explanation = explanation,
         var_model = vf, vg = vg, model_name = model_name, range_a = a, kappa = kap,
         cutoff = cutoff, width = width)
  }
  
  #' Extract model-based L50/L95 (fallback to empirical if needed)
  #' @return list(L50m, L95m, L50e, L95e)
  get_Ls <- function(Lres) {
    ex <- Lres$explanation
    pick <- function(kind, pp) { i <- which(ex$type == kind & abs(ex$p - pp) < 1e-8); if (length(i)) as.numeric(ex$L[i[1]]) else NA_real_ }
    list(L50m = pick("model_based", 0.50), L95m = pick("model_based", 0.95),
         L50e = pick("empirical",   0.50), L95e = pick("empirical",   0.95))
  }
  
  #' Choose radii R from L (micro=L50, local=L95); optional rounding
  choose_R_from_L <- function(Lres, round_to = 10) {
    Ls <- get_Ls(Lres)
    micro <- if (is.finite(Ls$L50m)) Ls$L50m else Ls$L50e
    local <- if (is.finite(Ls$L95m)) Ls$L95m else Ls$L95e
    if (!is.null(round_to)) {
      if (is.finite(micro)) micro <- round(micro/round_to)*round_to
      if (is.finite(local)) local <- round(local/round_to)*round_to
    }
    c(micro = micro, local = local)
  }
  
  #' Simple scale mismatch penalty g(R/L) with minimum at 1
  g_scale <- function(x) 0.5 * (x + 1/x) - 1
  
  #' Process standard deviation from variogram (sqrt structural sill)
  sigma_proc_from_vgm <- function(vf) {
    ps <- sum(vf$psill[vf$model != "Nug"], na.rm = TRUE)
    sqrt(ps)
  }
  
  #' Estimate alpha as drift strength proxy (altitude→value R^2)
  estimate_alpha <- function(sf_pts, value_col = "value", drift_col = "altitude", method = c("R2","sqrtR2")) {
    method <- match.arg(method)
    dat <- sf::st_drop_geometry(sf_pts)
    y <- dat[[value_col]]; x <- dat[[drift_col]]
    fit <- stats::lm(y ~ x)
    R2 <- summary(fit)$r.squared
    if (method == "sqrtR2") return(sqrt(R2))
    R2
  }
  
  
  # ============================= 40_predictors.R ============================= #
  
  # Smooth a raster at one or more metric radii (m)
  # Accepts a single numeric (e.g. 60) OR a vector (c(local=60, micro=20)).
  # Returns a *named list* of SpatRasters: names like "R60" or the names provided.
  smooth_raster_scales <- function(r, R, fun = mean) {
    stopifnot(inherits(r, "SpatRaster"))
    if (terra::crs(r, proj = TRUE) == "") stop("Raster must be projected (meters).")
    R <- as.numeric(R)
    if (length(R) < 1) stop("R must be a numeric radius (m) or vector.")
    
    # Create sensible names if none provided
    nm <- names(R)
    if (is.null(nm) || any(nm == "")) nm <- paste0("R", round(R))
    names(R) <- nm
    
    out <- vector("list", length(R))
    names(out) <- nm
    for (i in seq_along(R)) {
      rad <- R[i]
      w   <- .window_from_radius(rad, r)  # <-- pass numeric radius
      out[[i]] <- terra::focal(r, w = w, fun = fun, na.policy = "omit", na.rm = TRUE)
    }
    out
  }
  
  # DEM-derived, scale-matched features (relative height, slope, aspect)
  # R can be a single number or a named vector; layer names include that tag.
  derive_topo_features <- function(dem, R) {
    stopifnot(inherits(dem, "SpatRaster"), terra::nlyr(dem) == 1)
    sm <- smooth_raster_scales(dem, R)
    lay_list <- list()
    for (nm in names(sm)) {
      demR   <- sm[[nm]]
      altR   <- dem - demR
      slopeR <- terra::terrain(demR, v = "slope",  unit = "degrees")
      aspR   <- terra::terrain(demR, v = "aspect", unit = "degrees")
      lay_list[[paste0("altitude_",   nm)]] <- altR
      lay_list[[paste0("slope_deg_",  nm)]] <- slopeR
      lay_list[[paste0("aspect_sin_", nm)]] <- sin(aspR * pi/180)
      lay_list[[paste0("aspect_cos_", nm)]] <- cos(aspR * pi/180)
    }
    terra::rast(lay_list)
  }
  
  # Smooth arbitrary predictor rasters at R; returns a single SpatRaster stack.
  smooth_predictors <- function(preds, R) {
    if (inherits(preds, "SpatRaster")) preds <- list(extra = preds)
    stopifnot(is.list(preds), length(preds) > 0)
    
    out <- list()
    for (vname in names(preds)) {
      r  <- preds[[vname]]; stopifnot(inherits(r, "SpatRaster"))
      rs <- smooth_raster_scales(r, R)
      for (nm in names(rs)) {
        rr <- rs[[nm]]
        names(rr) <- paste0(vname, "_", nm, ".", names(rr))
        out[[paste0(vname, "_", nm)]] <- rr
      }
    }
    terra::rast(out)
  }
  
  
  # ============================== 50_workflow.R ============================== #
  
  # -----------------------------------------------------------------------------
  # run_one() – macht jetzt:
  #   • Raster (KED/OK/IDW/TPS/GAM/RF/Voronoi) + schreibt GeoTIFFs nach method_dir
  #   • A→D-Workflow (L50/L95, R), optional mit Extras
  #   • U-Kurve (R*-Tuning) + PNG
  #   • Benchmark bei R* + PNG + CSV
  #   • Error Budget + CSV
  #   • Rückgabe als gut strukturiertes List-Objekt pro Timestamp
  # -----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------- 
  # FIXED run_one(): keine self-referencing Defaults; robuster Env-Fallback
  # -----------------------------------------------------------------------------
  run_one <- function(
    v,
    m            = NULL,
    DEM_render   = NULL,   # nur fürs Anzeigen im Viewer; NICHT fürs Rechnen
    DEM_scale    = NULL,   # natives DEM (z.B. 2 m)
    method_dir   = NULL,
    fig_dir      = NULL,
    report_dir   = NULL,
    extra_preds  = NULL,   # optionale Zusatzprädiktoren (werden in Tuning/Benchmark genutzt)
    k_cv         = 5,
    model_tune   = "gam",
    use_blocks   = "L",    # "L" oder "R"
    save_rasters = TRUE,
    save_figs    = TRUE,
    save_tables  = TRUE
  ) {
    # --- kleine Helper -----------------------------------------------------------
    .need <- function(x, name) {
      if (!is.null(x)) return(x)
      y <- get0(name, envir = parent.frame(), inherits = TRUE, ifnotfound = NULL)
      if (is.null(y)) stop("run_one: object '", name, "' not found. Pass it explicitly or define it globally.")
      y
    }
    
    # DEM auf ~target_m Zellgröße bringen (integer-Faktor)
    dem_at <- function(dem_in, target_m, fun = mean) {
      stopifnot(inherits(dem_in, "SpatRaster"))
      resm <- mean(terra::res(dem_in))
      fact <- max(1L, round(as.numeric(target_m) / resm))
      out  <- terra::aggregate(dem_in, fact = c(fact, fact), fun = fun, na.rm = TRUE)
      names(out) <- "altitude"
      out
    }
    
    # Alle Methoden auf einem gegebenen Grid rechnen und schreiben (inkl. KED)
    .predict_all_on_grid <- function(m, v, grid_dem, out_dir, tag) {
      stopifnot(inherits(grid_dem, "SpatRaster"))
      pts <- .pts_df_from_var(m, v)
      r_vor   <- .safe(voro_simple(m, v, grid_dem),                              "Voronoi")
      r_idw   <- .safe(pred_idw(pts, grid_dem, idp = 2),                         "IDW")
      r_ok    <- .safe(pred_ok(pts, grid_dem),                                   "OK")
      r_trend <- .safe(pred_trend_lm(pts, grid_dem, use_xy = TRUE, use_alt = TRUE), "Trend")
      r_gam   <- .safe(pred_gam(pts, grid_dem),                                  "GAM")
      r_rf    <- .safe(pred_rf(pts, grid_dem),                                   "RF")
      
      # KED-Drift = exakt dieses Ziel-Grid
      driftR <- grid_dem; names(driftR) <- "altitude"
      r_ked  <- .safe(pred_ked(pts, driftR),                                     "KED")
      
      ts_tag <- slug(pretty_time(v))
      wr <- function(r, nm) if (isTRUE(save_rasters) && inherits(r, "SpatRaster"))
        terra::writeRaster(r, file.path(out_dir, sprintf("%s_%s_%s.tif", tolower(nm), ts_tag, tag)), overwrite = TRUE)
      
      wr(r_vor, "Voronoi"); wr(r_idw, "IDW"); wr(r_ok, "OK"); wr(r_trend, "Trend")
      wr(r_gam, "GAM");     wr(r_rf, "RF");   wr(r_ked,  "KED")
      
      invisible(list(Voronoi=r_vor, IDW=r_idw, OK=r_ok, Trend=r_trend, GAM=r_gam, RF=r_rf, KED=r_ked))
    }
    
    # --- Pflichtobjekte ziehen & Guarding ---------------------------------------
    m          <- .need(m,          "m")
    DEM_render <- .need(DEM_render, "DEM_render")   # bleibt für Viewer
    DEM_scale  <- .need(DEM_scale,  "DEM_scale")
    method_dir <- .need(method_dir, "method_dir")
    fig_dir    <- .need(fig_dir,    "fig_dir")
    report_dir <- .need(report_dir, "report_dir")
    
    if (is.null(v) || !nzchar(v)) return(invisible(NULL))
    dir.create(fig_dir,    recursive = TRUE, showWarnings = FALSE)
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    
    stamp  <- pretty_time(v)
    ts_tag <- slug(stamp)
    
    # Punkte prüfen
    pts <- .pts_df_from_var(m, v)
    if (!nrow(pts)) {
      warning("run_one: no points for ", v)
      return(invisible(NULL))
    }
    
    # ---------- 1) A→D Workflow (Basis + optional Extras) -----------------------
    wf <- workflow_L_R_predictors_blocks(
      sf_pts      = m,
      value_col   = v,
      dem         = DEM_scale,
      extra_preds = NULL,
      p           = c(0.95, 0.5),
      plot        = FALSE
    )
    
    wf_ex <- NULL; Ls_ex <- NULL; L95_ex <- NA_real_
    if (!is.null(extra_preds)) {
      wf_ex <- workflow_L_R_predictors_blocks(
        sf_pts      = m,
        value_col   = v,
        dem         = DEM_scale,
        extra_preds = extra_preds,
        p           = c(0.95, 0.5),
        plot        = FALSE
      )
      Ls_ex  <- get_Ls(wf_ex$L)
      L95_ex <- if (is.finite(Ls_ex$L95m)) Ls_ex$L95m else Ls_ex$L95e
    }
    
    # L95 (Fallbacks)
    Ls  <- get_Ls(wf$L)
    L95 <- if (is.finite(Ls$L95m)) Ls$L95m else if (is.finite(Ls$L95e)) Ls$L95e else 30
    
    # ---------- 2) Tuning (U-Kurve) --------------------------------------------
    tune <- tune_R_via_blockCV(
      wf, sf_pts = m, value_col = v,
      dem = DEM_scale, extra_preds = NULL,
      R_candidates = NULL, k = k_cv, model = model_tune, use_blocks = use_blocks
    )
    if (isTRUE(save_figs) && is.list(tune) && is.function(tune$plot)) {
      png(file.path(fig_dir, sprintf("u_curve_%s.png", ts_tag)), width = 900, height = 600)
      tune$plot(); dev.off()
    }
    
    tune_ex <- NULL
    if (!is.null(wf_ex)) {
      tune_ex <- tune_R_via_blockCV(
        wf_ex, sf_pts = m, value_col = v,
        dem = DEM_scale, extra_preds = extra_preds,
        R_candidates = NULL, k = k_cv, model = model_tune, use_blocks = use_blocks
      )
      if (isTRUE(save_figs) && is.list(tune_ex) && is.function(tune_ex$plot)) {
        png(file.path(fig_dir, sprintf("u_curve_extras_%s.png", ts_tag)), width = 900, height = 600)
        tune_ex$plot(); dev.off()
      }
    }
    
    # ---------- 3) Compute-Grids: DEM_Rstar & DEM_L95 ---------------------------
    Rstar <- suppressWarnings(as.numeric(tune$R_star))
    if (!is.finite(Rstar)) Rstar <- L95
    
    DEM_Rstar <- dem_at(DEM_scale, Rstar)  # Zellgröße ≈ R*
    DEM_L95   <- dem_at(DEM_scale, L95)    # Zellgröße ≈ L95
    
    # ---------- 4) ALLE Methoden auf beiden Grids rechnen -----------------------
    preds_R <- .predict_all_on_grid(m, v, DEM_Rstar, method_dir, tag = "Rstar")
    preds_L <- .predict_all_on_grid(m, v, DEM_L95,   method_dir, tag = "L95")
    
    # ---------- 5) Benchmark @ R* (einheitlich auf DEM_Rstar) -------------------
    bench <- benchmark_at_Rstar(
      wf, sf_pts = m, value_col = v,
      dem = DEM_Rstar, extra_preds = NULL,
      k = k_cv, use_blocks = use_blocks, model_for_tuning = model_tune,
      methods = c("OK","RF","Trend","GAM","IDW","Voronoi","KED"), make_plot = TRUE
    )
    if (isTRUE(save_figs) && !is.null(bench$plot)) {
      ggplot2::ggsave(file.path(fig_dir, sprintf("benchmark_%s.png", ts_tag)),
                      bench$plot, width = 8, height = 6, dpi = 200)
    }
    if (isTRUE(save_tables) && is.data.frame(bench$table)) {
      utils::write.csv(bench$table, file.path(report_dir, sprintf("benchmark_%s.csv", ts_tag)),
                       row.names = FALSE)
    }
    
    bench_ex <- NULL
    if (!is.null(wf_ex)) {
      bench_ex <- benchmark_at_Rstar(
        wf_ex, sf_pts = m, value_col = v,
        dem = DEM_Rstar, extra_preds = extra_preds,  # gleiche Grid-Logik auch mit Extras
        k = k_cv, use_blocks = use_blocks, model_for_tuning = model_tune,
        methods = c("OK","RF","Trend","GAM","IDW","Voronoi","KED"), make_plot = TRUE
      )
      if (isTRUE(save_figs) && !is.null(bench_ex$plot)) {
        ggplot2::ggsave(file.path(fig_dir, sprintf("benchmark_extras_%s.png", ts_tag)),
                        bench_ex$plot, width = 8, height = 6, dpi = 200)
      }
      if (isTRUE(save_tables) && is.data.frame(bench_ex$table)) {
        utils::write.csv(bench_ex$table, file.path(report_dir, sprintf("benchmark_extras_%s.csv", ts_tag)),
                         row.names = FALSE)
      }
    }
    
    # ---------- 6) Error Budget (OK-CV, R bzw. L95) -----------------------------
    cv <- cv_ok(
      sf_pts    = m,
      value_col = v,
      ref_r     = DEM_Rstar,      # gleiche Grid-Referenz
      fold_id   = wf$fold_id,     # bereits vorhandene Blockeinteilung
      vf        = wf$L$var_model
    )
    R_for_budget <- if (is.finite(tune$R_star)) tune$R_star else L95
    tab_err <- error_budget(
      y_obs = cv$y_obs, y_hat = cv$y_hat,
      vf = wf$L$var_model, R = R_for_budget, L95 = L95,
      sigma_inst = 0.5, alpha = 0.6
    )
    if (isTRUE(save_tables) && is.data.frame(tab_err)) {
      utils::write.csv(tab_err, file.path(report_dir, sprintf("error_budget_%s.csv", ts_tag)), row.names = FALSE)
    }
    
    tab_err_ex <- NULL
    if (!is.null(wf_ex)) {
      cv_ex <- cv_ok(
        sf_pts    = m,
        value_col = v,
        ref_r     = DEM_Rstar,
        fold_id   = wf_ex$fold_id,
        vf        = wf_ex$L$var_model
      )
      R_for_budget_ex <- if (!is.null(tune_ex) && is.finite(tune_ex$R_star)) tune_ex$R_star else L95_ex
      tab_err_ex <- error_budget(
        y_obs = cv_ex$y_obs, y_hat = cv_ex$y_hat,
        vf = wf_ex$L$var_model, R = R_for_budget_ex, L95 = L95_ex,
        sigma_inst = 0.5, alpha = 0.6
      )
      if (isTRUE(save_tables) && is.data.frame(tab_err_ex)) {
        utils::write.csv(tab_err_ex, file.path(report_dir, sprintf("error_budget_extras_%s.csv", ts_tag)), row.names = FALSE)
      }
    }
    
    # ---------- 7) Rückgabe + RDS -----------------------------------------------
    out <- list(
      ts_key     = v,
      stamp      = stamp,
      grids      = list(DEM_Rstar = DEM_Rstar, DEM_L95 = DEM_L95),
      files      = list(
        ucurve_png        = if (save_figs) file.path(fig_dir, sprintf("u_curve_%s.png", ts_tag)) else NA,
        ucurve_extras_png = if (!is.null(wf_ex) && save_figs) file.path(fig_dir, sprintf("u_curve_extras_%s.png", ts_tag)) else NA,
        bench_png         = if (!is.null(bench$plot) && save_figs) file.path(fig_dir, sprintf("benchmark_%s.png", ts_tag)) else NA,
        bench_ex_png      = if (!is.null(bench_ex) && !is.null(bench_ex$plot) && save_figs) file.path(fig_dir, sprintf("benchmark_extras_%s.png", ts_tag)) else NA,
        bench_csv         = if (save_tables) file.path(report_dir, sprintf("benchmark_%s.csv", ts_tag)) else NA,
        bench_ex_csv      = if (!is.null(wf_ex) && save_tables) file.path(report_dir, sprintf("benchmark_extras_%s.csv", ts_tag)) else NA,
        eb_csv            = if (save_tables) file.path(report_dir, sprintf("error_budget_%s.csv", ts_tag)) else NA,
        eb_ex_csv         = if (!is.null(wf_ex) && save_tables) file.path(report_dir, sprintf("error_budget_extras_%s.csv", ts_tag)) else NA
      ),
      wf        = wf,
      Ls        = Ls,
      tune      = tune,
      bench     = bench,
      errtab    = tab_err,
      wf_ex     = wf_ex,
      Ls_ex     = Ls_ex,
      tune_ex   = tune_ex,
      bench_ex  = bench_ex,
      errtab_ex = tab_err_ex
    )
    
    rds_path <- file.path(report_dir, sprintf("results_%s.RDS", ts_tag))
    saveRDS(out, rds_path)
    invisible(out)
  }
  
  
  
  
  #' Build spatial CV blocks and fold ids at ~L meters
  #' @param ref_r SpatRaster reference grid (projected meters)
  #' @param L target block size in meters
  #' @param pts sf POINTs for fold assignment
  #' @return list(block_raster, fold_id, cell_id)
  make_block_folds <- function(ref_r, L, pts) {
    ref_r <- .square_raster(ref_r)
    res_xy <- terra::res(ref_r)
    fact_xy <- pmax(1, round(L / res_xy))
    block_r <- terra::aggregate(ref_r, fact = fact_xy)
    pts_prj <- sf::st_transform(pts, terra::crs(block_r))
    cell_id <- terra::cellFromXY(block_r, sf::st_coordinates(pts_prj))
    fold_id <- as.integer(as.factor(cell_id))
    list(block_raster = block_r, fold_id = fold_id, cell_id = cell_id)
  }
  
  
  
  # --- Full L→R→predictors→blocks workflow (A–D) -------------------------------
  workflow_L_R_predictors_blocks <- function(sf_pts, value_col,
                                             dem = NULL,
                                             extra_preds = NULL,
                                             p = c(0.95, 0.5),
                                             micro = NULL,   # allow NULL to auto-choose
                                             local = NULL,   # allow NULL to auto-choose
                                             plot = FALSE, plot_type = c("base","ggplot2")) {
    
    pts <- sf_pts
    if (is.numeric(value_col)) {
      names(pts)[value_col] <- "value"
    } else {
      nm <- names(pts) == value_col
      if (!any(nm)) stop("value_col not found.")
      names(pts)[nm] <- "value"
    }
    pts <- pts[!is.na(pts$value), ]
    
    # A) Estimate L on this time-slice
    Lres <- estimate_L_all(pts, p = p, plot = plot, plot_type = plot_type)
    
    # Extract L50/L95 (prefer model-based, fallback to empirical)
    Ls  <- get_Ls(Lres)
    L50 <- if (is.finite(Ls$L50m)) Ls$L50m else Ls$L50e
    L95 <- if (is.finite(Ls$L95m)) Ls$L95m else Ls$L95e
    
    # B) Choose R (local + micro)
    R_local <- local %||% L95
    if (!is.finite(R_local) || R_local <= 0) R_local <- 30  # safe fallback
    R_micro <- micro %||% max(5, round(R_local/3))
    R       <- c(micro = as.numeric(R_micro), local = as.numeric(R_local))
    R_meters <- as.numeric(R_local)
    
    # Reference raster for grids/predictors
    ref_r <- NULL
    if (!is.null(dem)) ref_r <- dem
    if (is.null(ref_r) && !is.null(extra_preds)) ref_r <- extra_preds[[1]]
    if (is.null(ref_r)) stop("Provide 'dem' or at least one raster in 'extra_preds'.")
    
    # C) Build scale-matched predictors (smooth at R_meters)
    preds_list <- list(); res_checks <- list()
    
    if (!is.null(dem)) {
      topoR <- derive_topo_features(dem, R_meters)   # <-- numeric radius
      preds_list$topo <- topoR
      
      # resolution check for DEM (optional, keeps your table)
      res_xy <- terra::res(dem); res_m <- mean(res_xy)
      status <- if (!is.finite(L50) && !is.finite(L95)) "unknown"
      else if (res_m < L50) "too_fine"
      else if (res_m < L95) "ok_mid"
      else "too_coarse"
      todo <- switch(status,
                     too_fine = sprintf("Oversampled (%.2fm < L50). Smooth to ~R_local ≈ %.0fm.",
                                        res_m, round(R_local)),
                     ok_mid   = sprintf("Between L50 and L95. Consider smoothing to ~%.0fm.",
                                        round(R_local)),
                     too_coarse = sprintf("Coarser than L95 (%.2fm ≥ L95).", res_m),
                     "Fit failed; re-check variogram."
      )
      res_checks$DEM <- data.frame(
        section="Predictor resolution", item="DEM", value=sprintf("%.2f", res_m),
        units="m/cell", source="grid", status=status, todo=todo
      )
    }
    
    if (!is.null(extra_preds)) {
      # IMPORTANT: pass numeric radius, not the R vector
      smoothed <- smooth_predictors(extra_preds, R_meters)
      preds_list$extra <- smoothed
    }
    
    # --- BUILD & CLEAN predictor stack -----------------------------------------
    pred_stack <- terra::rast(Filter(Negate(is.null), preds_list))
    
    # 1) Auf DEM-Grid reprojizieren/resamplen (falls nötig)
    if (!is.null(pred_stack) && !terra::compareGeom(pred_stack, ref_r, stopOnError = FALSE)) {
      pred_stack <- terra::project(pred_stack, terra::crs(ref_r), method = "bilinear")
      pred_stack <- terra::resample(pred_stack, ref_r, method = "bilinear")
    }
    
    # 2) Alle Nicht-Finiten auf NA setzen (terra nutzt oft NaN als nodata)
    pred_stack <- terra::app(pred_stack, function(x) { x[!is.finite(x)] <- NA; x })
    
    # 3) Auf DEM-Maske maskieren (kein Off-Grid Gerümpel)
    pred_stack <- terra::mask(pred_stack, ref_r[[1]])
    
    # --- GRID DF ----------------------------------------------------------------
    coords  <- terra::xyFromCell(ref_r, 1:terra::ncell(ref_r))
    grid_df <- data.frame(cell = seq_len(nrow(coords)), x = coords[,1], y = coords[,2])
    
    if (!is.null(pred_stack)) {
      vals <- as.data.frame(terra::values(pred_stack))
      # Sicherheit: erneut Nicht-Finite raus
      for (j in seq_along(vals)) vals[[j]][!is.finite(vals[[j]])] <- NA
      names(vals) <- names(pred_stack)
      grid_df <- cbind(grid_df, vals)
    }
    
    # Spatial blocks/folds (default ≈ L95; you can switch to R_meters if you prefer)
    cv <- make_block_folds(ref_r, L = ifelse(is.finite(L95), L95, 30), pts = pts)
    
    list(
      L = Lres,
      R = R,                 # c(micro=..., local=...)
      R_meters = R_meters,   # single numeric radius for smoothing/blocks if needed
      predictors_stack = pred_stack,
      grid_df = grid_df,
      resolution_checks = res_checks,
      blocks = cv$block_raster,
      fold_id = cv$fold_id,
      stations_cell_id = cv$cell_id
    )
  }
  
  # --- Summaries for workflow A→D ----------------------------------------------
  summarize_workflow <- function(wf) {
    ex <- wf$L$explanation
    getL <- function(kind, pp) {
      i <- which(ex$type == kind & abs(ex$p - pp) < 1e-8)
      if (length(i)) as.numeric(ex$L[i[1]]) else NA_real_
    }
    L95m <- getL("model_based", 0.95); L50m <- getL("model_based", 0.50)
    L95e <- getL("empirical",   0.95); L50e <- getL("empirical",   0.50)
    
    rows_L <- rbind(
      data.frame(section="L (model)", item="L95_model", value=round(L95m,2), units="m",
                 source=wf$L$model_name, status="info",
                 todo="Use as practical range (R_local, CV blocks)."),
      data.frame(section="L (model)", item="L50_model", value=round(L50m,2), units="m",
                 source=wf$L$model_name, status="info",
                 todo="Half-range (sanity check)."),
      data.frame(section="L (emp)",   item="L95_emp",   value=round(L95e,2), units="m",
                 source="empirical", status="diag",
                 todo="Compare with model; large gap ⇒ recheck fit/bins."),
      data.frame(section="L (emp)",   item="L50_emp",   value=round(L50e,2), units="m",
                 source="empirical", status="diag",
                 todo="As above.")
    )
    
    R <- wf$R
    rows_R <- do.call(rbind, lapply(names(R), function(nm)
      data.frame(section="Chosen scales R", item=nm, value=R[[nm]], units="m",
                 source="rule R≈L",
                 status="action",
                 todo=if (nm=="local")
                   "Smooth predictors to this scale; set CV block size ≈ this."
                 else "Optional micro-scale for near-canopy effects.")
    ))
    
    block_m <- if (is.finite(L95m)) round(L95m) else if (is.finite(L95e)) round(L95e) else 30
    rows_cv   <- data.frame(section="CV", item="block_size", value=block_m, units="m",
                            source="L95 (pref. model)", status="action",
                            todo="Leave-block-out CV with this block size.")
    rows_fold <- data.frame(section="CV", item="n_folds",
                            value=length(unique(wf$fold_id)), units="count",
                            source="blocks", status="info",
                            todo="Inspect per-fold RMSE/MAE/Bias.")
    
    summary_tab <- rbind(rows_L, rows_R, rows_cv, rows_fold)
    rownames(summary_tab) <- NULL
    list(summary = summary_tab, defaults = default_process_table())
  }
  
  #' Default process → scale → predictor table (cm to ~500 m)
  default_process_table <- function() {
    data.frame(
      process = c(
        "Leaf boundary layer / needle-scale heat exchange",
        "Bark/leaf surface temperature & emissivity micro-contrast",
        "Soil crust / micro-roughness (cm–dm)",
        "Micro-relief (hummocks, root plates)",
        "Ground shading patterns (fine)",
        "Near-surface albedo patchiness (litter, moss)",
        "Local slope & aspect (micro)",
        "Sky View Factor (micro cells)",
        "Short canopy gaps / small clearings",
        "Wind shelter by small obstacles (H≈2–5 m)",
        "Edge effects (thin strips, small forest edges)",
        "SVF / shadow duration (neighborhood)",
        "Moisture patches / soil texture contrasts (local)",
        "Vegetation index patches (NDVI/NDMI) – local",
        "Advection from edges (leeward warming/cooling)",
        "Fetch / wind exposure (orographic channelling)",
        "Topographic position (hollow/ridge) – local",
        "Cold-air drainage tongues (gentle slopes)",
        "Terrain shading / horizon effects (low sun)",
        "Soil moisture gradients (drainage/accumulation)",
        "Edge plumes in strong wind (large clearings)",
        "Mesoscale topographic setting (plateau vs. valley)"
      ),
      scale_min_m = c(0.01, 0.01, 0.02, 1, 1, 1, 5, 5, 8, 5, 15, 15, 15, 20, 40, 50, 40, 100, 150, 120, 300, 300),
      scale_max_m = c(0.3, 0.5, 0.5, 5, 5, 5, 15, 15, 20, 20, 40, 40, 40, 60, 100, 120, 120, 300, 300, 250, 500, 500),
      typical_L_range = c("<1 m", "<1 m", "≈1 m", "2–10 m", "2–8 m", "2–8 m", "8–20 m", "8–20 m", "10–25 m", "10–25 m", "25–60 m", "25–60 m", "25–60 m", "30–70 m", "60–120 m", "70–150 m", "60–120 m", "120–300 m", "150–300 m", "120–250 m", "250–500 m", "250–500 m")
    )
  }
  
  #' Filter the default table by matched micro/local bands from wf
  #' @param wf result of workflow_L_R_predictors_blocks()
  #' @param show_all if TRUE, include match flags
  #' @return data.frame subset
  default_process_table_filtered <- function(wf, show_all = FALSE, micro_band = 0.5, local_band = 0.5, L95_band = c(0.5, 1.2)) {
    stopifnot(is.list(wf), !is.null(wf$L), !is.null(wf$R))
    tab <- default_process_table()
    ex  <- wf$L$explanation
    getL <- function(kind, pp) { i <- which(ex$type == kind & abs(ex$p - pp) < 1e-8); if (length(i)) as.numeric(ex$L[i[1]]) else NA_real_ }
    L50m <- getL("model_based", 0.50); L95m <- getL("model_based", 0.95)
    L50e <- getL("empirical",   0.50); L95e <- getL("empirical",   0.95)
    L50  <- if (is.finite(L50m)) L50m else L50e
    L95  <- if (is.finite(L95m)) L95m else L95e
    R_micro <- suppressWarnings(as.numeric(wf$R["micro"]))
    R_local <- suppressWarnings(as.numeric(wf$R["local"]))
    micro_win <- if (is.finite(R_micro)) c(R_micro*(1-micro_band), R_micro*(1+micro_band)) else c(NA_real_, NA_real_)
    local_win <- if (is.finite(R_local)) c(R_local*(1-local_band), R_local*(1+local_band)) else if (is.finite(L95)) c(L95*L95_band[1], L95*L95_band[2]) else c(NA_real_, NA_real_)
    overlaps <- function(a_min, a_max, b_min, b_max) { if (any(is.na(c(a_min,a_max,b_min,b_max)))) return(FALSE); (a_max >= b_min) && (a_min <= b_max) }
    tab$match_micro <- mapply(function(smin, smax) overlaps(smin, smax, micro_win[1], micro_win[2]), tab$scale_min_m, tab$scale_max_m)
    tab$match_local <- mapply(function(smin, smax) overlaps(smin, smax, local_win[1], local_win[2]), tab$scale_min_m, tab$scale_max_m)
    tab$overlap_with <- ifelse(tab$match_micro & tab$match_local, "micro+local", ifelse(tab$match_micro, "micro", ifelse(tab$match_local, "local", "none")))
    tab$R_micro <- if (is.finite(R_micro)) R_micro else NA_real_
    tab$R_local <- if (is.finite(R_local)) R_local else if (is.finite(L95)) L95 else NA_real_
    tab$L50 <- L50; tab$L95 <- L95
    if (isTRUE(show_all)) tab else subset(tab, overlap_with != "none")
  }
  
  
  # ============================ 60_blockcv_tuning.R =========================== #
  
  #' Internal: CV evaluation for a single R (build features @R and run spatial block-CV)
  .cv_for_R <- function(sf_pts, value_col, dem, extra_preds, R, L_for_blocks, k = 5, model = c("lm","gam")) {
    model <- match.arg(model)
    
    # Target column → 'value'
    pts <- sf_pts
    if (is.numeric(value_col)) names(pts)[value_col] <- "value" else names(pts)[names(pts) == value_col] <- "value"
    pts <- pts[!is.na(pts$value), ]
    if (nrow(pts) < max(5, k)) return(data.frame(R=R, RMSE=NA, MAE=NA, Bias=NA, n_folds=NA))
    
    # Reference raster (square, projected)
    ref_r <- .dem_ref(dem, extra_preds)
    
    # Build features at radius R (DEM topo + optional extras), aligned to ref_r
    topoR <- NULL
    if (inherits(dem, "SpatRaster")) {
      demA  <- .align_to(dem, ref_r)
      topoR <- try(derive_topo_features(demA, c(local = R)), silent = TRUE)
      if (inherits(topoR, "try-error")) topoR <- NULL
    }
    
    extraR <- NULL
    if (!is.null(extra_preds)) {
      preds_list    <- if (inherits(extra_preds, "SpatRaster")) list(extra = extra_preds) else extra_preds
      preds_aligned <- lapply(preds_list, .align_to, ref = ref_r)
      extraR <- try(smooth_predictors(preds_aligned, c(local = R)), silent = TRUE)
      if (inherits(extraR, "try-error")) extraR <- NULL
    }
    
    have <- Filter(Negate(is.null), list(topo = topoR, extra = extraR))
    if (!length(have)) return(data.frame(R=R, RMSE=NA, MAE=NA, Bias=NA, n_folds=NA))
    Xr <- terra::rast(have)
    
    # Extract predictors at station locations
    pts_prj <- sf::st_transform(pts, terra::crs(Xr))
    pred_df <- terra::extract(Xr, sf::st_coordinates(pts_prj))
    if (is.null(pred_df) || !nrow(pred_df)) return(data.frame(R=R, RMSE=NA, MAE=NA, Bias=NA, n_folds=NA))
    pred_df <- as.data.frame(pred_df)
    if ("ID" %in% names(pred_df)) pred_df$ID <- NULL
    for (j in seq_along(pred_df)) pred_df[[j]][!is.finite(pred_df[[j]])] <- NA
    
    dat <- data.frame(value = pts$value, pred_df, check.names = FALSE)
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < max(5, k)) return(data.frame(R=R, RMSE=NA, MAE=NA, Bias=NA, n_folds=NA))
    
    # Spatial folds on blocks ~ L_for_blocks
    folds <- make_block_folds(ref_r, L = L_for_blocks, pts = pts_prj)$fold_id
    folds <- folds[stats::complete.cases(pred_df)]
    tab <- table(folds); keep_folds <- as.integer(names(tab)[tab >= 2])
    use <- folds %in% keep_folds
    dat <- dat[use, ]; folds <- folds[use]
    if (length(unique(folds)) < 2) return(data.frame(R=R, RMSE=NA, MAE=NA, Bias=NA, n_folds=length(unique(folds))))
    
    uniq <- sort(unique(folds))
    mets <- lapply(uniq, function(fk) {
      tr <- dat[folds != fk, ]; te <- dat[folds == fk, ]
      pred <- try(.fit_and_predict(tr, te, response = "value", method = model), silent = TRUE)
      if (inherits(pred, "try-error")) return(NULL)
      .metrics(te$value, pred)
    })
    mets <- Filter(Negate(is.null), mets)
    if (!length(mets)) return(data.frame(R=R, RMSE=NA, MAE=NA, Bias=NA, n_folds=length(uniq)))
    M <- do.call(rbind, mets)
    for (nm in c("RMSE","MAE","Bias")) M[[nm]][!is.finite(M[[nm]])] <- NA_real_
    out <- colMeans(M, na.rm = TRUE)
    
    data.frame(R=R, RMSE=out["RMSE"], MAE=out["MAE"], Bias=out["Bias"], n_folds=length(uniq))
  }
  
  #' Tune R via block-CV (U-curve); optionally choose block size from R or L95
  #' @return list(table, R_star, plot=function())
  tune_R_via_blockCV <- function(wf, sf_pts, value_col, dem = NULL, extra_preds = NULL,
                                 R_candidates = NULL, k = 5, model = c("lm","gam"),
                                 use_blocks = c("R","L")) {
    model <- match.arg(model)
    use_blocks <- match.arg(use_blocks)
    
    # L50/L95 from wf (prefer model-based, fallback empirical)
    ex <- wf$L$explanation
    pickL <- function(kind,p) { i <- which(ex$type == kind & abs(ex$p - p) < 1e-8); if (length(i)) as.numeric(ex$L[i[1]]) else NA_real_ }
    L50 <- if (is.finite(pickL("model_based",0.50))) pickL("model_based",0.50) else pickL("empirical",0.50)
    L95 <- if (is.finite(pickL("model_based",0.95))) pickL("model_based",0.95) else pickL("empirical",0.95)
    
    # Reference raster
    ref_r <- .dem_ref(dem, extra_preds)
    
    # Candidate radii
    if (is.null(R_candidates)) {
      lo <- if (is.finite(L50)) L50 else max(5, min(20, L95/3))
      hi <- if (is.finite(L95)) L95 else max(30, lo*3)
      Rs <- unique(round(seq(lo, hi, length.out = 5)))
      px <- mean(terra::res(ref_r))
      R_candidates <- Rs[Rs >= 3*px]
    }
    
    block_size_for <- function(R) if (use_blocks == "R") R else if (is.finite(L95)) L95 else R
    
    rows <- lapply(
      R_candidates,
      function(R) {
        tryCatch(
          .cv_for_R(sf_pts, value_col, dem, extra_preds, R,
                    L_for_blocks = block_size_for(R), k = k, model = model),
          error = function(e) data.frame(R=R, RMSE=NA, MAE=NA, Bias=NA, n_folds=NA)
        )
      }
    )
    tab <- as.data.frame(do.call(rbind, rows))
    for (nm in c("RMSE","MAE","Bias")) if (nm %in% names(tab)) tab[[nm]][!is.finite(tab[[nm]])] <- NA_real_
    
    ok <- (("RMSE" %in% names(tab)) & is.finite(tab$RMSE)) & (("n_folds" %in% names(tab)) & tab$n_folds >= 2)
    R_star <- if (any(ok)) tab$R[which(ok)[ which.min(tab$RMSE[ok]) ]] else NA_real_
    
    plot_fun <- local({
      tt <- tab; m <- model; L50v <- L50; L95v <- L95; Rbest <- R_star
      function() {
        if (!requireNamespace("ggplot2", quietly = TRUE)) {
          plot(tt$R, tt$RMSE, type="b", pch=16, xlab="R (m)", ylab="RMSE", main="U-curve"); 
          if (is.finite(Rbest)) abline(v=Rbest, lty=2); if (is.finite(L50v)) abline(v=L50v, col="grey60", lty=3); if (is.finite(L95v)) abline(v=L95v, col="grey60", lty=3)
          return(invisible(NULL))
        }
        df <- subset(tt, is.finite(R) & is.finite(RMSE))
        p <- ggplot2::ggplot(df, ggplot2::aes(R, RMSE)) + ggplot2::geom_line() + ggplot2::geom_point() +
          ggplot2::labs(x="R (m)", y="RMSE (CV, drift-only)", title=sprintf("U-curve (%s)%s", m, if (is.finite(Rbest)) sprintf(" — R*=%dm", Rbest) else "")) +
          ggplot2::theme_minimal()
        if (is.finite(Rbest)) p <- p + ggplot2::geom_vline(xintercept = Rbest, linetype = 2)
        if (is.finite(L50v))  p <- p + ggplot2::geom_vline(xintercept = L50v,  color = "grey60", linetype = 3)
        if (is.finite(L95v))  p <- p + ggplot2::geom_vline(xintercept = L95v,  color = "grey60", linetype = 3)
        print(p); invisible(p)
      }
    })
    
    list(table = tab, R_star = R_star, plot = plot_fun)
  }
  
  
  
  
  # ============================== 70_benchmark.R ============================= #
  
  #' Benchmark multiple interpolation methods at tuned R* (spatial blocks)
  #' @return list(R_star, block_size, tuner, table, plot)
  benchmark_at_Rstar <- function(
    wf, sf_pts, value_col = "value", dem, extra_preds = NULL,
    k = 5, use_blocks = c("R","L95"), model_for_tuning = c("gam","lm"),
    methods = c("OK","RF","Trend","GAM","IDW","Voronoi","KED"),
    make_plot = TRUE
  ) {
    use_blocks <- match.arg(use_blocks, c("R","L","L95"))
    if (identical(use_blocks, "L")) use_blocks <- "L95"
    model_for_tuning  <- match.arg(model_for_tuning)
    
    .exists_fun <- function(nm) exists(nm, mode = "function", inherits = TRUE)
    
    # points
    pts <- sf_pts
    if (is.numeric(value_col)) names(pts)[value_col] <- "value" else names(pts)[names(pts) == value_col] <- "value"
    pts <- pts[!is.na(pts$value), ]
    if (nrow(pts) < 5) stop("Too few non-NA points.")
    
    # reference grid
    dem_sq <- .square_raster(dem)
    ref_r  <- dem_sq
    
    # tune R*
    res_tune <- tune_R_via_blockCV(
      wf, sf_pts = pts, value_col = "value",
      dem = dem_sq, extra_preds = extra_preds,
      R_candidates = NULL, k = k, model = model_for_tuning, use_blocks = "L"
    )
    R_star <- suppressWarnings(as.numeric(res_tune$R_star[1]))
    if (!is.finite(R_star)) {
      R_star <- suppressWarnings(as.numeric(wf$R["local"]))
      if (!is.finite(R_star)) {
        ex <- wf$L$explanation
        i  <- which(ex$type=="model_based" & abs(ex$p-0.95)<1e-8)
        R_star <- if (length(i)) as.numeric(ex$L[i[1]]) else 30
      }
      message(sprintf("R*: tuner failed; falling back to %g m.", R_star))
    }
    
    # predictors at R*
    Rvec       <- c(local = R_star)
    topo_star  <- derive_topo_features(dem_sq, Rvec)
    
    extra_star <- NULL
    if (!is.null(extra_preds)) {
      e_aligned <- if (inherits(extra_preds, "SpatRaster")) list(extra = extra_preds) else extra_preds
      e_aligned <- lapply(e_aligned, function(r){
        r2 <- .square_raster(r)
        if (!terra::compareGeom(r2, ref_r, stopOnError = FALSE)) {
          r2 <- terra::project(r2, terra::crs(ref_r))
          r2 <- terra::resample(r2, ref_r, method = "bilinear")
        }
        r2
      })
      extra_star <- smooth_predictors(e_aligned, Rvec)
    }
    
    # unified predictor stack (avoid NULLs)
    grid_stack <- terra::rast(Filter(Negate(is.null), list(topo_star, extra_star)))
    
    # drift raster for KED (same grid, smoothed to ~R*)
    R_for_drift <- if (is.finite(R_star)) R_star else mean(c(as.numeric(wf$R["local"]), wf$R_meters), na.rm = TRUE)
    if (!is.finite(R_for_drift)) R_for_drift <- 30
    drift_star <- try(safe_focal_mean(DEM_scale, R_for_drift, ref = dem_sq), silent = TRUE)
    if (!inherits(drift_star, "try-error") && inherits(drift_star, "SpatRaster")) {
      names(drift_star) <- "altitude"
    } else {
      drift_star <- ref_r
      names(drift_star) <- "altitude"
    }
    
    # (optional) grid_df if you need it later
    coords  <- terra::xyFromCell(ref_r, 1:terra::ncell(ref_r))
    vals    <- as.data.frame(terra::values(grid_stack))
    names(vals) <- names(grid_stack)
    grid_df <- cbind(data.frame(cell = seq_len(nrow(coords)), x = coords[,1], y = coords[,2]), vals)
    
    # helper to attach predictors to sf points (for non-KED algos that need xy+alt)
    enrich_pts <- function(m_sf) {
      xy <- sf::st_coordinates(m_sf)
      ext <- terra::extract(grid_stack, xy)
      if (!is.null(ext) && ncol(ext) > 0) {
        ext <- as.data.frame(ext)
        for (j in seq_along(ext)) m_sf[[names(ext)[j]]] <- ext[[j]]
      }
      m_sf
    }
    
    # block size and folds
    L95 <- {
      ex <- wf$L$explanation
      i  <- which(ex$type=="model_based" & abs(ex$p-0.95)<1e-8)
      v  <- if (length(i)) as.numeric(ex$L[i[1]]) else NA_real_
      if (!is.finite(v)) {
        j <- which(ex$type=="empirical" & abs(ex$p-0.95)<1e-8)
        v <- if (length(j)) as.numeric(ex$L[j[1]]) else NA_real_
      }
      v
    }
    block_size <- if (use_blocks == "R" && is.finite(R_star)) R_star else if (is.finite(L95)) L95 else 30
    fold_id    <- make_block_folds(ref_r, block_size, pts)$fold_id
    folds      <- sort(unique(fold_id))
    
    # methods to run
    keep <- intersect(methods, names(reg))
    keep <- keep[vapply(keep, function(m) .exists_fun(reg[[m]]), logical(1))]
    if (!length(keep)) stop("No available methods found among: ", paste(methods, collapse=", "))
    
    # CV loop
    results <- setNames(rep(list(NULL), length(keep)), keep)
    for (kfold in folds) {
      tr_idx <- (fold_id != kfold); te_idx <- (fold_id == kfold)
      pts_tr <- enrich_pts(pts[tr_idx, ])
      pts_te <- enrich_pts(pts[te_idx, ])
      
      .get_pts_df <- function(m_sf) {
        idx <- !is.na(m_sf$value)
        xy  <- sf::st_coordinates(m_sf[idx, ])
        alt <- m_sf$altitude[idx]
        if (!is.finite(mean(alt))) {
          # Fallback: nimm die erste Spalte, die wie ein Altitude-Feature aussieht
          cand <- grep("^altitude(_|$)", names(m_sf), value = TRUE)
          if (length(cand)) alt <- m_sf[[cand[1]]][idx]
        }
        data.frame(x = xy[,1], y = xy[,2], altitude = alt, value = m_sf$value[idx])
      }
      
      pts_df_tr <- .get_pts_df(pts_tr)
      if (!nrow(pts_df_tr)) next
      
      for (mname in keep) {
        fun <- get(reg[[mname]], mode = "function")
        
        r_pred <- tryCatch({
          if (identical(reg[[mname]], "voro_simple")) {
            # Voronoi needs sf + col name + template raster
            fun(pts_tr, "value", dem_sq)
            
          } else if (identical(mname, "KED")) {
            # Use same drift in train & predict
            # alt_tr statt ID holen
            alt_tr <- .safe_extract_vals(drift_star[[1]], sf::st_coordinates(pts_tr))
            pts_df_tr_ked <- pts_df_tr
            pts_df_tr_ked$altitude <- alt_tr
            pts_df_tr_ked <- pts_df_tr_ked[
              is.finite(pts_df_tr_ked$altitude) & is.finite(pts_df_tr_ked$value),
              , drop = FALSE
            ]
            r_pred <- if (!nrow(pts_df_tr_ked)) NULL else pred_ked(pts_df_tr_ked, drift_star)
            
          } else {
            fun(pts_df_tr, dem_sq)
          }
        }, error = function(e) NULL)
        
        if (!inherits(r_pred, "SpatRaster")) next
        xy_te <- sf::st_coordinates(pts_te)
        ext   <- terra::extract(r_pred, xy_te)
        if (is.null(ext) || ncol(ext) < 1) next
        pv    <- as.numeric(ext[,1]); ob <- as.numeric(pts_te$value)
        results[[mname]] <- rbind(results[[mname]], data.frame(method = mname, obs = ob, pred = pv))
      }
    }
    
    # metrics table
    non_empty <- Filter(Negate(is.null), results)
    tbl <- if (!length(non_empty)) {
      data.frame(method = keep, RMSE = NA_real_, MAE = NA_real_, Bias = NA_real_)
    } else {
      do.call(rbind, lapply(non_empty, function(df) {
        e <- df$pred - df$obs
        data.frame(
          method = df$method[1],
          RMSE   = sqrt(mean(e^2, na.rm = TRUE)),
          MAE    = mean(abs(e),   na.rm = TRUE),
          Bias   = mean(e,        na.rm = TRUE)
        )
      }))
    }
    bench <- unique(tbl); bench$R_star <- R_star; bench$block_size <- block_size; rownames(bench) <- NULL
    
    plt <- NULL
    if (isTRUE(make_plot) && requireNamespace("ggplot2", quietly = TRUE)) {
      plt <- ggplot2::ggplot(bench, ggplot2::aes(x = reorder(method, RMSE), y = RMSE)) +
        ggplot2::geom_col() + ggplot2::coord_flip() +
        ggplot2::labs(title = sprintf("Benchmark at R* = %dm (blocks: %dm)", round(R_star), round(block_size)),
                      x = "Method", y = "Block-CV RMSE") +
        ggplot2::theme_minimal()
    }
    
    list(R_star = R_star, block_size = block_size, tuner = res_tune, table = bench, plot = plt)
  }
  
  # ============================== 80_panel_plot.R ============================ #
  
  #' Build a multi-panel time-series raster plot (shared scale)
  #' Note: uses global pretty_time() for labels; no duplicate defined here.
  #' @return list(plot, r_stack, r_plot, df, png_path, pdf_path, layout, downsample_factor, value_range)
  timeseries_panel <- function(kriged_list, plot_area = NULL, stations_pos = NULL,
                               cells_target = 150000, max_cols = 4, label_pretty_time = TRUE,
                               out_png = "timeseries_panel_grid.png", out_pdf = "timeseries_panel_grid.pdf",
                               dpi = 200, fill_label = "Temperature") {
    require(terra); require(dplyr); require(tidyr); require(ggplot2); require(sf); require(stringr); require(lubridate)
    valid <- Filter(Negate(is.null), kriged_list)
    if (length(valid) < 1) stop("No valid rasters in `kriged_list`.")
    to_spatr <- function(s) { nm <- names(s); lyr <- if ("var1.pred" %in% nm) "var1.pred" else nm[1]; terra::rast(s[lyr]) }
    r_list  <- lapply(valid, to_spatr); r_stack <- terra::rast(r_list)
    nm <- names(valid)
    if (!is.null(nm) && all(nchar(nm) > 0)) {
      if (label_pretty_time) nm <- vapply(nm, pretty_time, character(1))
      names(r_stack) <- nm
    } else names(r_stack) <- sprintf("t%02d", seq_len(terra::nlyr(r_stack)))
    r_use <- tryCatch({ if (!is.null(plot_area)) terra::crop(r_stack, terra::vect(plot_area)) else r_stack }, error = function(e) r_stack)
    nc1  <- terra::ncell(r_use[[1]]); fact <- max(1L, round(sqrt(nc1 / cells_target)))
    r_plot <- if (fact > 1) terra::aggregate(r_use, fact = fact) else r_use
    df <- terra::as.data.frame(r_plot, xy = TRUE) |> tidyr::pivot_longer(cols = -(x:y), names_to = "layer", values_to = "value")
    df$layer <- factor(df$layer, levels = names(r_plot))
    val_range <- range(df$value, na.rm = TRUE)
    crs_r <- terra::crs(r_plot, proj = TRUE)
    stations_g <- tryCatch(if (!is.null(stations_pos)) sf::st_transform(stations_pos, crs_r) else NULL, error = function(...) stations_pos)
    area_g     <- tryCatch(if (!is.null(plot_area))   sf::st_transform(plot_area,  crs_r) else NULL, error = function(...) plot_area)
    nl <- nlevels(df$layer); ncol_facets <- min(max_cols, nl); nrow_facets <- ceiling(nl / ncol_facets)
    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::coord_equal() + ggplot2::scale_fill_viridis_c(limits = val_range, na.value = NA) +
      { if (!is.null(area_g))     ggplot2::geom_sf(data = area_g,     color = "red",  fill = NA, linewidth = 0.5, inherit.aes = FALSE) } +
      { if (!is.null(stations_g)) ggplot2::geom_sf(data = stations_g, color = "black", size = 0.8, inherit.aes = FALSE) } +
      ggplot2::facet_wrap(~ layer, ncol = ncol_facets, scales = "fixed") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), strip.text = element_text(face = "bold"), legend.position = "right") +
      ggplot2::labs(fill = fill_label)
    w_in <- 4.5 * ncol_facets; h_in <- 4.5 * nrow_facets
    ggplot2::ggsave(out_png, p, width = w_in, height = h_in, dpi = dpi)
    ggplot2::ggsave(out_pdf, p, width = w_in, height = h_in)
    list(plot = p, r_stack = r_stack, r_plot = r_plot, df = df, png_path = out_png, pdf_path = out_pdf, layout = c(nrow = nrow_facets, ncol = ncol_facets), downsample_factor = fact, value_range = val_range)
  }
  
  
  # ============================= 90_error_budget.R ============================ #
  
  # ============================= 90_error_budget.R ============================ #
  
  #' Block cross‑validation for Ordinary Kriging (OK)
  #'
  #' Runs leave-block-out CV using spatial blocks derived from a reference raster
  #' (same CRS/grid as your DEM). Optionally uses a provided variogram model;
  #' otherwise fits an OK variogram on each training fold with automap.
  #'
  #' @param sf_pts   sf POINT with a target column (name or index via `value_col`).
  #' @param value_col Name or index of the target column (default "value").
  #' @param ref_r    SpatRaster used to define CRS and the block grid.
  #' @param L        Block size in meters. Ignored if `fold_id` is supplied.
  #' @param vf       Optional gstat variogram model (from `estimate_L_all()` or your own).
  #'                 If NULL, a model is re-fitted in each training fold.
  #' @param fold_id  Optional integer vector of fold IDs (same length as `sf_pts`).
  #'
  #' @return list with `y_obs`, `y_hat`, `fold_id`, and `vf_used` (NULL or the
  #'         model passed in).
  #'
  #' @examples
  #' # Using L95 from workflow and a global variogram model:
  #' # cv <- cv_ok(m_points, value_col = "value", ref_r = DTM,
  #' #             L = get_Ls(wf$L)$L95m, vf = wf$L$var_model)
  #' # eb <- error_budget(cv$y_obs, cv$y_hat, vf = wf$L$var_model,
  #' #                    R = bench$R_star, L95 = get_Ls(wf$L)$L95m)
  cv_ok <- function(sf_pts, value_col = "value", ref_r, L = NULL, vf = NULL, fold_id = NULL) {
    stopifnot(inherits(sf_pts, "sf"), inherits(ref_r, "SpatRaster"))
    
    # 1) Normalize target column name → 'value'
    pts <- sf_pts
    if (is.numeric(value_col)) {
      names(pts)[value_col] <- "value"
    } else if (value_col != "value") {
      names(pts)[names(pts) == value_col] <- "value"
    }
    if (!"value" %in% names(pts)) stop("value_col not found.")
    pts <- pts[!is.na(pts$value), ]
    if (nrow(pts) < 4) stop("Too few points for CV (n < 4).")
    
    # 2) Project to ref_r CRS and make folds if needed
    pts <- sf::st_transform(pts, terra::crs(ref_r))
    n <- nrow(pts)
    if (is.null(fold_id)) {
      if (!is.finite(L)) stop("Provide block size L (meters) or a 'fold_id' vector.")
      fold_id <- make_block_folds(ref_r, L = L, pts = pts)$fold_id
    } else {
      if (length(fold_id) != n) stop("fold_id must have same length as rows in sf_pts.")
      fold_id <- as.integer(fold_id)
    }
    folds <- sort(unique(fold_id))
    if (length(folds) < 2) stop("Need at least 2 folds for CV.")
    
    # 3) Init outputs
    y_obs <- as.numeric(pts$value)
    y_hat <- rep(NA_real_, n)
    
    # 4) Loop folds: fit OK and predict at held-out points
    for (fk in folds) {
      idx_te <- which(fold_id == fk)
      idx_tr <- which(fold_id != fk)
      if (length(idx_tr) < 3 || length(idx_te) < 1) next
      
      train <- pts[idx_tr, ]
      test  <- pts[idx_te, ]
      
      train_sp <- methods::as(train, "Spatial")
      test_sp  <- methods::as(test,  "Spatial")
      
      # Per-fold variogram if none is supplied
      vgm_use <- vf
      if (is.null(vgm_use)) {
        # automap sometimes fails on very small/degenerate sets → try/catch
        vfit <- try(automap::autofitVariogram(value ~ 1, train), silent = TRUE)
        if (inherits(vfit, "try-error")) next
        vgm_use <- vfit$var_model
      }
      
      # gstat OK prediction at test locations
      kr <- try(gstat::krige(value ~ 1, locations = train_sp, newdata = test_sp, model = vgm_use), silent = TRUE)
      if (inherits(kr, "try-error")) next
      phat <- as.numeric(kr@data[["var1.pred"]])
      y_hat[idx_te] <- phat
    }
    
    # 5) Sanity: ensure we filled all test points at least once
    if (all(is.na(y_hat))) warning("cv_ok: no predictions produced (check variogram/folds).")
    list(y_obs = y_obs, y_hat = y_hat, fold_id = fold_id, vf_used = vf)
  }
  
  #' Decompose RMSE into bias, sensor noise, scale mismatch, and model rest
  #' @param vf gstat variogram model (from estimate_L_all)
  #' @param R chosen smoothing radius; @param L95 practical range
  #' @param sigma_inst instrument noise (SD); @param alpha weight (0..1) for scale term
  #' @return tibble with components and squared contributions
  error_budget <- function(y_obs, y_hat, vf, R, L95, sigma_inst = 0.25, alpha = 0.6) {
    stopifnot(length(y_obs) == length(y_hat))
    e     <- y_hat - y_obs
    rmse  <- sqrt(mean(e^2, na.rm = TRUE))
    bias  <- mean(e, na.rm = TRUE)
    sigma_proc <- sigma_proc_from_vgm(vf)
    x <- as.numeric(R) / as.numeric(L95)
    gx <- g_scale(x)
    term_bias2   <- bias^2
    term_inst2   <- sigma_inst^2
    term_scale2  <- (alpha * sigma_proc * gx)^2
    term_model2  <- max(0, rmse^2 - term_bias2 - term_inst2 - term_scale2)
    tibble::tibble(
      component = c("Bias", "Sensor (σ_inst)", "Scale mismatch", "Model rest", "Total (RMSE)"),
      formula   = c("Bias^2", "σ_inst^2", "(α σ_proc g(R/L))^2", "σ_model^2", "—"),
      value     = c(sqrt(term_bias2), sqrt(term_inst2), sqrt(term_scale2), sqrt(term_model2), rmse),
      squared   = c(term_bias2, term_inst2, term_scale2, term_model2, rmse^2)
    )
  }
  
  report_scales <- function(wf, label = "run") {
    Ls <- get_Ls(wf$L)
    cat("\n--- Scale report [", label, "] ---\n", sep = "")
    cat("Model:", wf$L$model_name, "\n")
    cat(sprintf("L50  model/emp: %.1f / %.1f m\n", Ls$L50m, Ls$L50e))
    cat(sprintf("L95  model/emp: %.1f / %.1f m\n", Ls$L95m, Ls$L95e))
    cat(sprintf("R (micro/local): %.1f / %.1f m\n", wf$R["micro"], wf$R["local"]))
    cat(sprintf("CV folds: %d\n", length(unique(wf$fold_id))))
  }
  
  # --- small helpers (move near the top, before first use) ---
  slug <- function(x) { x <- gsub("[^0-9A-Za-z_-]+","-", x); x <- gsub("-+","-", x); gsub("(^-|-$)","", x) }
  
  # pick densest time-slice instead of hard-coding 1L (optional but robust)
  pick_densest_index <- function(sf_wide, var_names) {
    nn <- sapply(var_names, function(v) sum(is.finite(sf_wide[[v]])))
    which.max(nn)
  }
  
  # --- add to all_functions_1.R -------------------------------------------------
  .window_from_radius <- function(r_m, rast) {
    resm <- mean(terra::res(rast))            # meters per cell
    rad_cells <- max(1L, round(as.numeric(r_m) / resm))
    w <- 2L * rad_cells + 1L                  # odd side length
    if (w < 3L) w <- 3L
    if (w %% 2L == 0L) w <- w + 1L
    w
  }
  
  
  # Ersetzt deine alte Version
  safe_focal_mean <- function(x, r_m, ref = NULL, shape = c("circle","square")) {
    stopifnot(inherits(x, "SpatRaster"))
    shape <- match.arg(shape)
    
    # 1) quadratische Zellen & ggf. auf 'ref' alignieren (CRS + Grid + Maske)
    xr <- .square_raster(x)
    if (!is.null(ref)) xr <- .align_to(xr, ref)
    
    # 2) Gewichtsmatrix in Meter-Radius
    w <- NULL
    if (shape == "circle") {
      # terra::focalMat nutzt Distanz in Karten­einheiten (hier Meter)
      w <- try(terra::focalMat(xr, r_m, type = "circle"), silent = TRUE)
    }
    if (inherits(w, "try-error") || is.null(w)) {
      # Fallback: quadratisches Fenster über Zellenradius
      k <- .window_from_radius(r_m, xr)
      w <- matrix(1, nrow = k, ncol = k)
    }
    
    terra::focal(xr, w = w, fun = mean, na.policy = "omit", na.rm = TRUE)
  }
  
  # ------------------------------------------------------------------------------
  
  quiet <- function(expr, log_file = file.path(report_dir, "run_warnings.log")) {
    dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
    warn_buf <- character()
    res <- withCallingHandlers(
      suppressMessages(suppressWarnings(capture.output({ val <- force(expr) }, file = NULL))),
      warning = function(w) { warn_buf <<- c(warn_buf, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    if (length(warn_buf)) writeLines(unique(warn_buf), log_file, useBytes = TRUE, sep = "\n")
    invisible(res)
  }
  
  
  
  #------------------------------------------------------------------------------------
  
  # (Remove the old global ".have"/present pre-check – it can mask the local one)
  # library(...) etc. can stay as you had them
  
  library(shiny)
  library(leaflet)
  library(DT)
  library(sf)
  library(terra)
  library(htmltools)
  
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  
  pretty_time <- function(x) {
    vapply(x, function(s) {
      if (grepl("^A\\d{14}$", s)) {
        ts <- as.POSIXct(substr(s, 2, 15), format = "%Y%m%d%H%M%S", tz = "UTC")
        format(ts, "%Y-%m-%d %H:%M")
      } else if (grepl("^A\\d{8}(_D)?$", s)) {
        ts <- as.Date(substr(s, 2, 9), format = "%Y%m%d")
        format(ts, "%Y-%m-%d")
      } else s
    }, character(1))
  }
  
  slug <- function(x) { x <- gsub("[^0-9A-Za-z_-]+","-", x); x <- gsub("-+","-", x); gsub("(^-|-$)","", x) }
  
  build_explanations <- function(fig_dir, pick_ts) {
    ts_label <- slug(pretty_time(pick_ts))
    files <- c(
      "timeseries_panel_grid.png",
      "timeseries_panel_grid.pdf",
      sprintf("u_curve_%s.png", ts_label),
      sprintf("u_curve_extras_%s.png", ts_label),
      sprintf("benchmark_%s.png", ts_label),
      sprintf("benchmark_extras_%s.png", ts_label)
    )
    paths <- file.path(fig_dir, files)
    desc  <- c(
      "Per-timestep KED previews; dots=stations; red=plot boundary.",
      "Same as PNG, vector PDF.",
      "U-curve for tuning R via block-CV (drift-only).",
      "U-curve with extra predictors.",
      "Method comparison at R* (lower RMSE is better).",
      "Benchmark with extras at R*."
    )
    keep <- file.exists(paths)
    out <- as.list(desc[keep]); names(out) <- basename(paths[keep])
    out
  }
  
  run_mc_viewer <- function(
      vars,
      method_dir,
      fig_dir,
      stations_pos,
      plot_area,
      wf = NULL, wf_ex = NULL,
      tune = NULL, tune_ex = NULL,
      bench = NULL, bench_ex = NULL,
      tab_err = NULL, tab_err_ex = NULL,
      explanations = NULL
  ) {
    stopifnot(length(vars) >= 1)
    
    if (is.null(explanations) && exists("build_explanations")) {
      explanations <- build_explanations(fig_dir, vars[1])
    } else if (is.null(explanations)) {
      explanations <- list()
    }
    
    # erwarteter TIFF-Pfad
    # erwartet: <method>_<stamp>.tif
    # erwartet: <method>_<stamp>.tif in method_dir/
    # Fallback: <ts>_interpolated.tif in der Elternmappe von method_dir (== out_dir)
    # --- Drop-in: robustes Datei-Matching mit _Rstar/_L95-Präferenz --------------
    raster_path <- function(method, ts) {
      stopifnot(length(method) == 1, length(ts) == 1)
      m <- tolower(method)
      
      # Tokens aus dem Timestamp (nutzt deine pretty_time-Logik)
      .ts_tokens <- function(ts_key) {
        raw <- tolower(as.character(ts_key))                       # "A20230829120000"
        pty <- tolower(pretty_time(ts_key))                        # "2023-08-29 12:00"
        slug_pt <- gsub("[^0-9A-Za-z_-]+","-", pty)                # "2023-08-29-12-00"
        d14 <- sub("^a", "", raw)                                  # "20230829120000" / "20230829"
        ymd  <- if (nchar(d14) >= 8) substr(d14,1,8) else NA_character_
        hhmm <- if (nchar(d14) >= 12) substr(d14,9,12) else NA_character_
        comp1 <- if (!is.na(ymd) && !is.na(hhmm))
          paste0(substr(ymd,1,4),"-",substr(ymd,5,6),"-",substr(ymd,7,8),"-",
                 substr(hhmm,1,2),"-",substr(hhmm,3,4)) else NA_character_
        comp2 <- gsub("-", "", comp1)                               # "202308291200"
        ymd_dash <- if (!is.na(ymd)) paste0(substr(ymd,1,4),"-",substr(ymd,5,6),"-",substr(ymd,7,8)) else NA_character_
        unique(na.omit(c(raw, slug_pt, comp1, comp2, ymd_dash, ymd)))
      }
      toks <- .ts_tokens(ts)
      tok_rx <- gsub("[-_]", "[-_]", toks)  # Trennertoleranz
      
      # 1) methodenspezifische Dateien in method_dir sammeln
      all_files <- list.files(method_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
      if (!length(all_files)) return(NA_character_)
      b <- tolower(basename(all_files))
      
      # nur Dateien mit Präfix "<method>_"
      keep_pref <- grepl(paste0("^", m, "_"), b)
      files_m <- all_files[keep_pref]; b_m <- b[keep_pref]
      if (length(files_m)) {
        # Scoring: längster passender Token im Dateinamen → höherer Score
        score <- vapply(seq_along(b_m), function(i) {
          max(c(0, vapply(tok_rx, function(rx) if (grepl(rx, b_m[i], perl = TRUE)) nchar(rx) else 0L, integer(1))))
        }, numeric(1))
        
        if (any(score > 0)) {
          best <- files_m[score == max(score)]
          # Präferenz: *_Rstar.tif vor *_L95.tif (bei gleichem Score)
          bbest <- tolower(basename(best))
          idxR <- grep("_rstar\\.tif$", bbest)
          if (length(idxR)) return(best[idxR[1]])
          idxL <- grep("_l95\\.tif$", bbest)
          if (length(idxL)) return(best[idxL[1]])
          return(best[1])
        }
      }
      
      # 2) Fallback NUR für KED: Preview <ts>_interpolated.tif in out_dir
      if (toupper(method) %in% c("KED","PREVIEW")) {
        out_dir <- dirname(method_dir)
        prev <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
        if (length(prev)) {
          bp <- tolower(basename(prev))
          rx_prev <- paste0("^(", paste0(tok_rx, collapse = "|"), ")_interpolated(_wgs84)?\\.tif$")
          hit <- grepl(rx_prev, bp, perl = TRUE)
          if (any(hit)) return(prev[which(hit)[1]])
        }
      }
      
      NA_character_
    }
    
    
    
    all_methods <- c("KED","OK","IDW","GAM","RF","Voronoi","Trend")
    # Prüfe über alle timestamps, ob es mind. 1 Datei pro Methode gibt
    present <- all_methods[vapply(
      all_methods,
      function(m) any(file.exists(sapply(vars, function(v) raster_path(m, v)))),
      logical(1)
    )]
    if (!length(present)) present <- all_methods
    
    # Vektoren in WGS84 für Leaflet (ohne Z/M)
    area_ll     <- sf::st_transform(sf::st_zm(plot_area,    drop = TRUE, what = "ZM"), 4326)
    stations_ll <- sf::st_transform(sf::st_zm(stations_pos, drop = TRUE, what = "ZM"), 4326)
    
    ts_choices <- stats::setNames(vars, pretty_time(vars))
    
    ui <- shiny::fluidPage(
      tags$head(tags$style(HTML("#map{height:70vh} .dt-caption{caption-side:top;font-weight:600;margin-bottom:.5em}"))),
      titlePanel("Microclimate Viewer"),
      sidebarLayout(
        sidebarPanel(width = 3,
                     radioButtons("method","Method",choices = present, selected = present[1]),
                     selectInput("ts","Timestamp", choices = ts_choices, selected = vars[1]),
                     fluidRow(
                       column(6, actionButton("prev","⟵ Prev", width="100%")),
                       column(6, actionButton("nextbtn","Next ⟶", width="100%"))
                     ),
                     # ---- Kontrastwahl -----------------------------------------------------
                     selectInput(
                       "stretch", "Contrast",
                       choices = c(
                         "Histogram equalized"           = "histeq",
                         "Quantiles (9 bins)"            = "quant9",
                         "Linear robust (2–98%)"         = "lin_robust",
                         "Linear (layer range)"          = "lin_local",
                         "Linear (global per timestamp)" = "lin_global"
                       ),
                       selected = "histeq"
                     ),
                     hr(),
                     h5("Figures"),
                     selectInput("which_fig","Choose figure",
                                 choices = names(explanations), selected = if (length(explanations)) names(explanations)[1] else NULL),
                     uiOutput("fig_link"),
                     verbatimTextOutput("fig_desc")
        ),
        mainPanel(width = 9,
                  tabsetPanel(
                    tabPanel("Map", leafletOutput("map"), br(), htmlOutput("map_explain")),
                    tabPanel("Tuning & Scales",
                             fluidRow(column(6, plotOutput("u_plot")), column(6, plotOutput("u_plot_ex"))),
                             hr(), h4("Scales & radii"), tableOutput("scale_tbl")
                    ),
                    tabPanel("Benchmarks",
                             h4("No extras"), DTOutput("bench_tbl"),
                             hr(), h4("With extras"), DTOutput("bench_tbl_ex")
                    ),
                    tabPanel("Error budget",
                             h4("No extras"), DTOutput("err_tbl"),
                             hr(), h4("With extras"), DTOutput("err_tbl_ex")
                    )
                  )
        )
      )
    )
    
    # ---- robuster Auto-Zoom auf Raster (3857), NA-sicher -----------------------
    zoom_to_raster <- function(lp, r3857) {
      e <- try(terra::ext(r3857), silent = TRUE)
      if (inherits(e, "try-error")) return(lp)
      box_num <- c(e$xmin, e$xmax, e$ymin, e$ymax)
      if (any(!is.finite(box_num))) return(lp)
      
      corners <- data.frame(x = c(e$xmin, e$xmax), y = c(e$ymin, e$ymax))
      pts_3857 <- sf::st_as_sf(corners, coords = c("x","y"), crs = 3857)
      pts_ll   <- try(sf::st_transform(pts_3857, 4326), silent = TRUE)
      if (inherits(pts_ll, "try-error")) return(lp)
      
      bb <- sf::st_bbox(pts_ll)
      if (any(!is.finite(as.numeric(bb)))) return(lp)
      
      padx <- as.numeric(bb["xmax"] - bb["xmin"]) * 0.02
      pady <- as.numeric(bb["ymax"] - bb["ymin"]) * 0.02
      
      leaflet::fitBounds(
        lp,
        lng1 = as.numeric(bb["xmin"]) - padx,
        lat1 = as.numeric(bb["ymin"]) - pady,
        lng2 = as.numeric(bb["xmax"]) + padx,
        lat2 = as.numeric(bb["ymax"]) + pady
      )
    }
    
    .safe_quantile_bin_pal <- function(vals, n_target = 9) {
      vals <- vals[is.finite(vals)]
      if (!length(vals)) return(list(
        pal        = leaflet::colorNumeric(viridisLite::viridis(256), c(0,1), na.color = NA),
        leg_values = numeric(0),
        note       = "linear_fallback_no_vals"
      ))
      ux <- sort(unique(vals))
      if (length(ux) < 3) {
        mm <- range(vals); span <- diff(mm); if (!is.finite(span) || span == 0) span <- (abs(mm[1]) + 1) * 1e-6
        dom <- c(mm[1]-0.5*span, mm[2]+0.5*span)
        return(list(
          pal        = leaflet::colorNumeric(viridisLite::viridis(256), dom, na.color = NA),
          leg_values = vals,
          note       = "linear_fallback_low_unique"
        ))
      }
      n   <- max(3, min(n_target, length(ux)))
      brw <- as.numeric(stats::quantile(vals, seq(0,1,length.out=n+1), na.rm=TRUE, names=FALSE))
      br  <- unique(brw)
      if (length(br) < 3) {
        mm <- range(vals); span <- diff(mm); if (!is.finite(span) || span == 0) span <- (abs(mm[1]) + 1) * 1e-6
        dom <- c(mm[1]-0.5*span, mm[2]+0.5*span)
        return(list(
          pal        = leaflet::colorNumeric(viridisLite::viridis(256), dom, na.color = NA),
          leg_values = vals,
          note       = "linear_fallback_flat_quantiles"
        ))
      }
      list(
        pal        = leaflet::colorBin(viridisLite::viridis(length(br)-1), domain = vals, bins = br, na.color = NA),
        leg_values = vals,
        note       = sprintf("bin_quantiles_%d", length(br)-1)
      )
    }
    
    .safe_addLegend <- function(lp, pal, values, title) {
      res <- try(leaflet::addLegend(lp, position="bottomright", pal=pal, values=values, title=title, opacity=1), silent=TRUE)
      if (inherits(res, "try-error")) {
        mm <- range(values, na.rm=TRUE); span <- diff(mm); if (!is.finite(span) || span==0) span <- (abs(mm[1])+1)*1e-6
        dom <- c(mm[1]-0.5*span, mm[2]+0.5*span)
        pal2 <- leaflet::colorNumeric(viridisLite::viridis(256), dom, na.color=NA)
        leaflet::addLegend(lp, position="bottomright", pal=pal2, values=values,
                           title=paste0(title," (linear fallback)"), opacity=1)
      } else res
    }
    
    
    server <- function(input, output, session) {
      # Globale Domäne pro Timestamp (falls mehrere Methoden existieren)
      domain_for_ts <- reactive({
        req(input$ts)
        rngs <- lapply(all_methods, function(m){
          fp <- raster_path(m, input$ts)
          if (!file.exists(fp)) return(NULL)
          r  <- try(terra::rast(fp), silent = TRUE)
          if (inherits(r, "try-error")) return(NULL)
          as.numeric(terra::minmax(r))
        })
        rngs <- Filter(Negate(is.null), rngs)
        if (!length(rngs)) return(NULL)
        c(
          min(vapply(rngs, `[`, numeric(1), 1), na.rm=TRUE),
          max(vapply(rngs, `[`, numeric(1), 2), na.rm=TRUE)
        )
      })
      
      # Static map shell
      output$map <- renderLeaflet({
        leaflet() |>
          addTiles(group = "OSM") |>
          leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron, group = "Carto") |>
          leaflet::addProviderTiles(leaflet::providers$Esri.WorldGrayCanvas, group = "Esri") |>
          addMapPane("raster",  zIndex = 410) |>
          addMapPane("vectors", zIndex = 420) |>
          addLayersControl(
            baseGroups = c("Carto","OSM","Esri"),
            overlayGroups = c("Layer","Boundary","Stations"),
            options = layersControlOptions(collapsed = FALSE)
          )
      })
      
      # Serve fig_dir als /figs/...
      shiny::addResourcePath("figs", normalizePath(fig_dir, winslash = "/", mustWork = FALSE))
      observe({
        figs <- basename(list.files(fig_dir, pattern = "\\.(png|pdf)$", full.names = FALSE))
        updateSelectInput(session, "which_fig",
                          choices = figs,
                          selected = if (length(figs)) figs[1] else character(0))
      })
      output$fig_link <- renderUI({
        req(input$which_fig, nzchar(input$which_fig))
        tags$p(tags$a(href = file.path("/figs", input$which_fig), target = "_blank", "Open figure in new tab"))
      })
      output$fig_desc <- renderText({ explanations[[input$which_fig]] %||% "" })
      
      # ---- Map refresh ----------------------------------------------------------
      
      # sichere Quantil-Palette: passt n an und fällt zurück auf linear, wenn nötig
      .safe_quantile_pal <- function(vals, n_target) {
        vals <- vals[is.finite(vals)]
        if (!length(vals)) return(list(pal = leaflet::colorNumeric(viridisLite::viridis(256), c(0,1), na.color=NA),
                                       leg_values = numeric(0), mode = "linear_fallback"))
        u <- length(unique(vals))
        if (u <= 2) {
          # zu wenig Spread → lineare Palette mit leichtem Puffer
          mm <- range(vals); span <- diff(mm); if (!is.finite(span) || span == 0) span <- (abs(mm[1])+1)*1e-6
          dom <- c(mm[1] - 0.5*span, mm[2] + 0.5*span)
          return(list(pal = leaflet::colorNumeric(viridisLite::viridis(256), dom, na.color=NA),
                      leg_values = vals, mode = "linear_fallback"))
        }
        n <- max(3, min(n_target, u))  # nie mehr Bins als unique Werte
        list(pal = leaflet::colorQuantile(viridisLite::viridis(n), domain = vals, n = n, na.color = NA),
             leg_values = vals, mode = sprintf("quantile_n=%d", n))
      }
      
      # sichere Legende: falls Quantil-Legende scheitert → linearer Fallback
      .safe_addLegend <- function(lp, pal, values, title) {
        res <- try(leaflet::addLegend(lp, position = "bottomright",
                                      pal = pal, values = values, title = title, opacity = 1),
                   silent = TRUE)
        if (inherits(res, "try-error")) {
          mm <- range(values, na.rm = TRUE)
          span <- diff(mm); if (!is.finite(span) || span == 0) span <- (abs(mm[1])+1)*1e-6
          dom <- c(mm[1] - 0.5*span, mm[2] + 0.5*span)
          pal2 <- leaflet::colorNumeric(viridisLite::viridis(256), dom, na.color = NA)
          leaflet::addLegend(lp, position = "bottomright", pal = pal2, values = values,
                             title = paste0(title, " (linear fallback)"), opacity = 1)
        } else res
      }
      
      refresh_map <- function(method, ts_key) {
        f  <- raster_path(method, ts_key)
        lp <- leafletProxy("map")
        
        # boundary + stations
        lp <- lp |>
          clearGroup("Boundary") |>
          clearGroup("Stations") |>
          addPolylines(data = area_ll, weight = 2, color = "#cc2b2b", group = "Boundary",
                       options = pathOptions(pane = "vectors")) |>
          addCircleMarkers(data = stations_ll, radius = 4, fillOpacity = 1,
                           color = "black", stroke = TRUE, weight = 1, group = "Stations",
                           options = pathOptions(pane = "vectors"))
        
        if (!file.exists(f)) {
          lp <- lp |> clearImages() |> clearControls()
          output$map_explain <- renderUI(HTML(
            sprintf("<b>%s @ %s</b><br><span style='color:#a00'>Missing file:</span> %s",
                    method, pretty_time(ts_key), htmltools::htmlEscape(basename(f)))
          ))
          bb <- sf::st_bbox(area_ll)
          leaflet::fitBounds(lp, bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])
          return(invisible(NULL))
        }
        
        sel_method <- method
        
        # ---- detail-schonendes Reproject → EPSG:3857 ---------------------------
        r_spat  <- terra::rast(f)
        res_src <- mean(terra::res(r_spat))
        
        ext0      <- terra::ext(r_spat)
        center_xy <- c((ext0$xmin + ext0$xmax)/2, (ext0$ymin + ext0$ymax)/2)
        pt_src    <- sf::st_sfc(sf::st_point(center_xy), crs = sf::st_crs(terra::crs(r_spat)))
        lat_deg   <- try(as.numeric(sf::st_coordinates(sf::st_transform(pt_src, 4326))[2]), silent = TRUE)
        if (inherits(lat_deg, "try-error") || !is.finite(lat_deg)) lat_deg <- 0
        
        scale_wm   <- 1 / max(0.25, cos(lat_deg * pi/180))  # clamp gegen extreme
        oversample <- 0.75                                   # <1 → feiner → Details bleiben
        res_3857   <- res_src * scale_wm * oversample
        
        # immer nearest → keine Glättung/Verwischung
        r_3857 <- terra::project(r_spat, "EPSG:3857", method = "near", res = res_3857)
        # ------------------------------------------------------------------------
        
        rr <- raster::raster(r_3857)
        raster::setMinMax(rr)
        # Werte
        vals <- raster::getValues(rr); vals <- vals[is.finite(vals)]
        if (!length(vals)) {# Werte lesen
          vals <- raster::getValues(rr)
          vals <- vals[is.finite(vals)]
          if (!length(vals)) {
            lp <- lp |> clearImages() |> clearControls()
            output$map_explain <- renderUI(HTML(
              sprintf("<b>%s @ %s</b><br><span style='color:#a00'>All values are NA.</span>",
                      sel_method, pretty_time(ts_key))
            ))
            # trotzdem auf Boundary zoomen, damit die Karte nicht „leer grau“ wirkt
            bb <- sf::st_bbox(area_ll)
            leaflet::fitBounds(lp, bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])
            return(invisible(NULL))
          }
          }
        
        mode <- isolate(input$stretch %||% "histeq")
        glob <- try(domain_for_ts(), silent = TRUE); if (inherits(glob, "try-error")) glob <- NULL
        
        make_lin_pal <- function(range_vec) leaflet::colorNumeric(viridisLite::viridis(256), domain=range_vec, na.color=NA)
        
        pal <- NULL; leg_values <- vals; legend_note <- NULL
        if (identical(mode, "histeq")) {
          sp <- .safe_quantile_bin_pal(vals, n_target = 256); pal <- sp$pal; leg_values <- sp$leg_values; legend_note <- sp$note
        } else if (identical(mode, "quant9")) {
          sp <- .safe_quantile_bin_pal(vals, n_target = 9);   pal <- sp$pal; leg_values <- sp$leg_values; legend_note <- sp$note
        } else if (identical(mode, "lin_robust")) {
          qs <- stats::quantile(vals, c(0.02,0.98), na.rm=TRUE, names=FALSE); if (!all(is.finite(qs)) || qs[1]==qs[2]) qs <- range(vals, na.rm=TRUE)
          pal <- make_lin_pal(qs); legend_note <- "linear_robust_2-98%"
        } else if (identical(mode, "lin_local")) {
          mm <- range(vals, na.rm=TRUE); if (mm[1]==mm[2]) { eps <- (abs(mm[1])+1)*1e-6; mm <- c(mm[1]-eps, mm[2]+eps) }
          pal <- make_lin_pal(mm); legend_note <- "linear_local"
        } else { # lin_global
          if (is.null(glob) || length(glob)!=2 || !all(is.finite(glob)) || glob[1]>=glob[2]) {
            glob <- range(vals, na.rm=TRUE); if (glob[1]==glob[2]) { eps <- (abs(glob[1])+1)*1e-6; glob <- c(glob[1]-eps, glob[2]+eps) }
          }
          pal <- make_lin_pal(glob); legend_note <- "linear_global"
        }
        
        lp <- lp |>
          clearImages() |>
          clearControls() |>
          leaflet::addRasterImage(rr, colors=pal, opacity=0.85, project=FALSE, group="Layer",
                                  options = leaflet::tileOptions(pane="raster"))
        
        lp <- .safe_addLegend(lp, pal, leg_values,
                              paste(sel_method, pretty_time(ts_key),
                                    if (!is.null(legend_note)) paste0(" • ", legend_note) else ""))
        
        lp <- lp |>
          leaflet::addLayersControl(
            baseGroups = c("Carto","OSM","Esri"),
            overlayGroups = c("Layer","Boundary","Stations"),
            options = leaflet::layersControlOptions(collapsed = FALSE)
          )
        
        # Auto-zoom
        lp <- zoom_to_raster(lp, r_3857)
        
        output$map_explain <- renderUI(HTML(paste0(
          "<b>Layer:</b> ", sel_method, " @ ", pretty_time(ts_key),
          "<br><b>File:</b> ", htmltools::htmlEscape(basename(f)),
          sprintf("<br><b>Legend stretch:</b> %s", mode)
        )))
        invisible(NULL)
      }
      
      # initial & reactive
      observe({ refresh_map(input$method %||% present[1], input$ts %||% vars[1]) })
      observeEvent(list(input$method, input$ts, input$stretch), {
        refresh_map(input$method, input$ts)
      }, ignoreInit = TRUE)
      
      # Prev/Next
      observeEvent(input$prev, {
        cur <- isolate(input$ts); i <- match(cur, vars); if (is.na(i)) i <- 1
        i2 <- if (i == 1) length(vars) else i - 1
        updateSelectInput(session, "ts", selected = vars[i2])
      }, ignoreInit = TRUE)
      
      observeEvent(input$nextbtn, {
        cur <- isolate(input$ts); i <- match(cur, vars); if (is.na(i)) i <- 1
        i2 <- if (i == length(vars)) 1 else i + 1
        updateSelectInput(session, "ts", selected = vars[i2])
      }, ignoreInit = TRUE)
      
      # Tuning-Plots
      output$u_plot <- renderPlot({
        req(tune, !is.null(tune$plot))
        tune$plot()     # gibt ggplot zurück und druckt ihn
      })
      
      output$u_plot_ex <- renderPlot({
        req(tune_ex, !is.null(tune_ex$plot))
        tune_ex$plot()
      })
      
      # Scales-Tabelle
      output$scale_tbl <- renderTable({
        rows <- list()
        addL <- function(wf_obj, tag) {
          if (is.null(wf_obj)) return()
          Ls <- try(get_Ls(wf_obj$L), silent = TRUE); if (inherits(Ls, "try-error")) return()
          rows[[length(rows) + 1]] <<- data.frame(
            Run = tag,
            L50m = round(Ls$L50m, 1), L95m = round(Ls$L95m, 1),
            L50e = round(Ls$L50e, 1), L95e = round(Ls$L95e, 1),
            R_micro = round(wf_obj$R["micro"], 1), R_local = round(wf_obj$R["local"], 1),
            check.names = FALSE
          )
        }
        addL(wf, "base"); addL(wf_ex, "extras")
        if (!length(rows)) return(NULL)
        do.call(rbind, rows)
      }, striped = TRUE, bordered = TRUE, spacing = "s")
      
      # Benchmarks
      output$bench_tbl <- renderDT({
        req(bench, is.data.frame(bench$table))
        datatable(bench$table, rownames = FALSE,
                  caption = tags$caption(class = "dt-caption",
                                         sprintf("Benchmark (no extras) — R*=%s, blocks=%s", bench$R_star %||% NA, bench$block_size %||% NA)),
                  options = list(pageLength = 10, dom = "tip")) |>
          formatRound(c("RMSE","MAE","Bias"), digits = 3)
      })
      output$bench_tbl_ex <- renderDT({
        req(bench_ex, is.data.frame(bench_ex$table))
        datatable(bench_ex$table, rownames = FALSE,
                  caption = tags$caption(class = "dt-caption",
                                         sprintf("Benchmark (with extras) — R*=%s, blocks=%s", bench_ex$R_star %||% NA, bench_ex$block_size %||% NA)),
                  options = list(pageLength = 10, dom = "tip")) |>
          formatRound(c("RMSE","MAE","Bias"), digits = 3)
      })
      
      # Error budgets
      output$err_tbl <- renderDT({
        req(tab_err, is.data.frame(tab_err))
        datatable(tab_err, rownames = FALSE,
                  caption = tags$caption(class = "dt-caption", "Error budget (no extras)"),
                  options = list(pageLength = 10, dom = "tip")) |>
          formatRound(c("value","squared"), digits = 3)
      })
      output$err_tbl_ex <- renderDT({
        req(tab_err_ex, is.data.frame(tab_err_ex))
        datatable(tab_err_ex, rownames = FALSE,
                  caption = tags$caption(class = "dt-caption", "Error budget (with extras)"),
                  options = list(pageLength = 10, dom = "tip")) |>
          formatRound(c("value","squared"), digits = 3)
      })
    }
    
    shinyApp(ui, server)
  }
  
  
  make_mc_report <- function(res, out_pdf = "Microclimate_Report.pdf") {
    stopifnot(is.list(res), length(out_pdf) == 1)
    
    # ---- leichte Abhängigkeiten
    have_png   <- requireNamespace("png",   quietly = TRUE)
    have_terra <- requireNamespace("terra", quietly = TRUE)
    stopifnot(requireNamespace("ggplot2", quietly = TRUE))
    library(ggplot2)
    library(grid)      # base grid
    library(grDevices) # pdf()
    
    # --------- kleine Helfer ----------------------------------------------------
    fmt <- function(x, d = 1) {
      if (length(x) < 1 || !is.finite(x)) return("NA")
      format(round(as.numeric(x), d), nsmall = d, trim = TRUE)
    }
    kv  <- function(...) {
      x <- list(...)
      val <- vapply(x, function(xx) {
        if (is.null(xx) || length(xx) < 1 || is.na(xx)) return("NA")
        as.character(xx[1])
      }, character(1))
      data.frame(key = names(x), value = val, row.names = NULL, check.names = FALSE)
    }
    grob_title <- function(txt, size = 18, face = "bold") grid::textGrob(txt, x=.02, y=.98, just=c("left","top"),
                                                                         gp=grid::gpar(fontsize=size, fontface=face))
    grob_sub <- function(text, size = 12) {
      grid::textGrob(
        label = text,
        x = grid::unit(1, "npc") - grid::unit(2, "mm"),
        y = grid::unit(1, "npc") - grid::unit(2, "mm"),
        just = c("right", "top"),
        gp = grid::gpar(col = "#555555", fontsize = size)
      )
    }
    
    gp=grid::gpar(col="#555", fontsize=size)
    draw_lines <- function(lines, y_start=.85, step=.045, size=11) {
      y <- y_start
      for (ln in lines) {
        grid::grid.text(ln, x=.04, y=y, just=c("left","center"),
                        gp=grid::gpar(fontsize=size))
        y <- y - step
      }
    }
    draw_kv_table <- function(df, y0=.8, step=.05, size=11) {
      if (nrow(df) == 0) return(invisible(NULL))
      y <- y0
      for (i in seq_len(nrow(df))) {
        grid::grid.text(df$key[i],   x=.04, y=y, just=c("left","center"),
                        gp=grid::gpar(fontsize=size, col="#333"))
        grid::grid.text(df$value[i], x=.48, y=y, just=c("left","center"),
                        gp=grid::gpar(fontsize=size, fontface="bold"))
        y <- y - step
      }
    }
    grob_image_file <- function(png_path, fallback_plot = NULL) {
      if (!is.null(png_path) && file.exists(png_path) && have_png) {
        img <- png::readPNG(png_path)
        return(rasterGrob(img, interpolate = TRUE))
      }
      if (inherits(fallback_plot, "ggplot")) return(ggplotGrob(fallback_plot))
      grid::textGrob("Kein Bild verfügbar", gp=gpar(col="#999"))
    }
    res_m <- function(r) {
      if (have_terra && inherits(r, "SpatRaster")) mean(terra::res(r)) else NA_real_
    }
    
    # --------- Daten aus res -----------------------------------------------------
    ts_label <- res$stamp %||% tryCatch(pretty_time(res$ts_key), error=function(...) res$ts_key)
    Ls       <- res$Ls %||% list(L50m=NA, L95m=NA, L50e=NA, L95e=NA)
    Rstar    <- suppressWarnings(as.numeric(res$tune$R_star))
    Rstar    <- if (is.finite(Rstar)) Rstar else NA_real_
    bench    <- res$bench; bench_ex <- res$bench_ex
    files    <- res$files %||% list()
    grids    <- res$grids %||% list()
    res_R    <- res_m(grids$DEM_Rstar)
    res_L    <- res_m(grids$DEM_L95)
    block_sz <- bench$block_size %||% NA_real_
    
    # Best rows
    best_row <- function(tab) if (is.data.frame(tab) && nrow(tab)) tab[which.min(tab$RMSE), , drop=FALSE] else NULL
    b0 <- best_row(bench$table)
    bE <- best_row(if (!is.null(bench_ex)) bench_ex$table else NULL)
    
    # Fehlerbudget Plots (falls vorhanden)
    eb   <- res$errtab
    ebE  <- res$errtab_ex
    eb_plot <- function(tab, title="Error budget") {
      if (!is.data.frame(tab) || !"component" %in% names(tab)) return(NULL)
      ggplot(tab, aes(x=reorder(component, value), y=value)) +
        geom_col() + coord_flip() +
        labs(title = title, x = NULL, y = "σ (Contribution)") +
        theme_minimal(base_size = 11)
    }
    p_eb  <- eb_plot(eb,  "Error budget (base)")
    p_ebE <- eb_plot(ebE, "Error budget (extras)")
    
    # ---------------------------------------------------------------------------
    #                               PDF START
    # ---------------------------------------------------------------------------
    pdf(out_pdf, width = 8.27, height = 11.69) # A4 hoch, in Zoll
    
    # --- Seite 1: Cover ---------------------------------------------------------
    grid.newpage()
    grid.draw(grob_title("Microclimate Report"))
    grid.draw(grob_sub(sprintf("Zeitstempel: %s", ts_label)))
    
    lines1 <- c(
      "Ziel: räumliche Interpolation je Zeitschritt, Skalenanalyse (Variogramm),",
      "Tuning des Drift-Radius R* via Block-CV (U-Kurve) und Methodenvergleich.",
      "KED nutzt ein auf die Zielauflösung geglättetes DEM als externe Drift."
    )
    draw_lines(lines1, y_start=.86)
    
    tbl_cover <- kv(
      `L50 (modell)`  = sprintf("%s m", fmt(Ls$L50m)),
      `L95 (modell)`  = sprintf("%s m", fmt(Ls$L95m)),
      `L50 (emp.)`    = sprintf("%s m", fmt(Ls$L50e)),
      `L95 (emp.)`    = sprintf("%s m", fmt(Ls$L95e)),
      `R* (getuned)`  = sprintf("%s m", fmt(Rstar)),
      `CV Blockgröße` = sprintf("%s m", fmt(block_sz)),
      `Rasterauflösung @R*`  = ifelse(is.finite(res_R), sprintf("%s m", fmt(res_R,2)), "NA"),
      `Rasterauflösung @L95` = ifelse(is.finite(res_L), sprintf("%s m", fmt(res_L,2)), "NA")
    )
    draw_kv_table(tbl_cover, y0=.60)
    
    grid.text("Dieser Report fasst die wichtigsten Ergebnisse zusammen und legt die Konzepte kurz dar.",
              x=.04, y=.20, just=c("left","center"), gp=gpar(col="#555"))
    
    # --- Seite 2: U-Kurve -------------------------------------------------------
    grid.newpage()
    grid.draw(grob_title("Tuning: U-Kurve (Block-CV)"))
    uc_png <- files$ucurve_png %||% NULL
    grob <- grob_image_file(uc_png, fallback_plot = res$tune$plot %||% NULL)
    grid.draw(editGrob(grob, vp = viewport(x=.5, y=.52, width=.9, height=.75)))
    draw_lines(c(
      "Vertikale Linie = R* (falls stabil gefunden).",
      "Gestrichelte Linien markieren L50 / L95 aus dem Variogramm."
    ), y_start=.12, step=.05)
    
    # --- Seite 3: Benchmarks (Base) ---------------------------------------------
    grid.newpage()
    grid.draw(grob_title("Benchmark der Methoden — Basisprädiktoren"))
    if (is.list(bench) && !is.null(bench$plot)) {
      grid.draw(editGrob(ggplotGrob(bench$plot), vp = viewport(x=.5, y=.55, width=.9, height=.8)))
    } else {
      grid.text("Kein Benchmark-Plot verfügbar", x=.5, y=.55)
    }
    if (!is.null(b0)) {
      draw_kv_table(kv(
        `Beste Methode` = as.character(b0$method[1]),
        `RMSE`          = fmt(b0$RMSE[1], 3),
        `MAE`           = fmt(b0$MAE[1], 3),
        `Bias`          = fmt(b0$Bias[1], 3)
      ), y0=.18, step=.06)
    }
    
    # --- Seite 4: Benchmarks (Extras) -------------------------------------------
    grid.newpage()
    grid.draw(grob_title("Benchmark der Methoden — mit Extras"))
    if (!is.null(bench_ex) && !is.null(bench_ex$plot)) {
      grid.draw(editGrob(ggplotGrob(bench_ex$plot), vp = viewport(x=.5, y=.55, width=.9, height=.8)))
    } else {
      grid.text("Kein Benchmark-Plot (Extras) verfügbar", x=.5, y=.55)
    }
    if (!is.null(bE)) {
      draw_kv_table(kv(
        `Beste Methode` = as.character(bE$method[1]),
        `RMSE`          = fmt(bE$RMSE[1], 3),
        `MAE`           = fmt(bE$MAE[1], 3),
        `Bias`          = fmt(bE$Bias[1], 3)
      ), y0=.18, step=.06)
    }
    
    # --- Seite 5: Fehlerbudget ---------------------------------------------------
    grid.newpage()
    grid.draw(grob_title("Fehlerbudget"))
    if (inherits(p_eb, "ggplot")) {
      grid.draw(editGrob(ggplotGrob(p_eb),  vp = viewport(x=.28, y=.60, width=.52, height=.70)))
    } else {
      grid.text("Kein Error-Budget (Base) vorhanden", x=.28, y=.60)
    }
    if (inherits(p_ebE, "ggplot")) {
      grid.draw(editGrob(ggplotGrob(p_ebE), vp = viewport(x=.72, y=.60, width=.52, height=.70)))
    } else {
      grid.text("Kein Error-Budget (Extras) vorhanden", x=.72, y=.60)
    }
    
    # Kurzer Text
    draw_lines(c(
      "Zerlegung von RMSE in Bias, Sensoraussteuerung (σ_inst),",
      "Skalen-Mismatch (α·σ_proc·g(R/L)) und Modellrest."
    ), y_start=.12)
    
    # --- Seite 6: Raster-Grids & KED Drift -------------------------------------
    grid.newpage()
    grid.draw(grob_title("Grids & Drift"))
    draw_kv_table(kv(
      `DEM @R* — mittlere Zellgröße`  = ifelse(is.finite(res_R), sprintf("%s m", fmt(res_R,2)), "NA"),
      `DEM @L95 — mittlere Zellgröße` = ifelse(is.finite(res_L), sprintf("%s m", fmt(res_L,2)), "NA"),
      `KED Drift`                     = "Glattes DEM (einzelner Layer: 'altitude')"
    ), y0=.78, step=.06)
    
    # versuche die gerenderten R* / L95 Rasterbilder (falls PNGs gespeichert) – ansonsten Hinweis
    p_R <- files$bench_png %||% NULL  # als Platzhalter; echte KED-Raster-Previews i.d.R. auf der Map
    grid.text("Hinweis: Kartenansichten (Leaflet) bieten die visuelle Kontrolle; dieser Report fasst die Kernergebnisse als Zahlen/Plots.",
              x=.04, y=.18, just=c("left","center"), gp=gpar(col="#555"))
    
    # --- Seite 7: Konzepte -------------------------------------------------------
    grid.newpage()
    grid.draw(grob_title("Begriffe & Konzepte"))
    bullets <- c(
      "L50 / L95: praktikable Korrelationslängen aus dem Variogramm (50% / 95% der Sill).",
      "R (micro/local): aus L abgeleitete Glättungsradien; hier wird R_local ≈ L95 genutzt.",
      "R*: aus Block-CV (U-Kurve) getunter Drift-Radius (Minimum der RMSE-Kurve).",
      "Block-CV: räumliche Faltung mit Blockkante ≈ L95 bzw. R*, um Leckage zu vermeiden.",
      "KED: Kriging mit externer Drift (hier geglättetes DEM). Drift wird auf dem Zielgrid beprobt.",
      "Fehlerbudget: RMSE wird in Bias, Sensor-Rauschen, Skalen-Mismatch und Modellrest zerlegt."
    )
    draw_lines(bullets, y_start=.82, step=.06)
    
    dev.off()
    message("✔ Report geschrieben: ", normalizePath(out_pdf))
    invisible(out_pdf)
  }
  
