## ======================================================================
## pipemodel_functions.R  —  nur Funktionen, keine Seiteneffekte
## ======================================================================
# ========================= I/O helpers (tables & plots) ======================


# ---- Export-Helper  ----
fn_fig <- function(stem, ext = "png") file.path(fig_dir, sprintf("%s.%s", stem, ext))
fn_ras <- function(stem, ext = "tif") file.path(ras_dir, sprintf("%s.%s", stem, ext))

save_plot_min <- function(p, file, width = 9, height = 6, dpi = 150, bg = "white") {
  # Speichert ggplot ODER "druckbare" Plot-Objekte
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (inherits(p, "ggplot")) {
    ggplot2::ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi, bg = bg)
  } else {
    grDevices::png(filename = file, width = width, height = height, units = "in", res = dpi, bg = bg)
    print(p)
    grDevices::dev.off()
  }
  invisible(normalizePath(file, winslash = "/"))
}

safe_save_plot <- function(p, file, width = 9, height = 6, dpi = 150, bg = "white") {
  try(save_plot_min(p, file, width, height, dpi, bg), silent = TRUE)
}

save_raster_min <- function(r, file, overwrite = TRUE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(r, file, overwrite = overwrite)
  invisible(normalizePath(file, winslash = "/"))
}


# Save a table in CSV (+ optional HTML via gt, XLSX via openxlsx/writexl)
# Robust to tibbles, list cols (ignored), and mistaken positional args.
# Save a table as CSV (always), HTML (if gt is installed), and XLSX
# file_stem: full path without extension, e.g. fn_tab("metrics_T14_base")
save_table_readable <- function(df, file_stem,
                                title = NULL,
                                digits = 3,
                                make_dirs = TRUE,
                                verbose = FALSE) {
  if (!inherits(df, "data.frame")) df <- as.data.frame(df)
  
  # Drop list-cols so write.csv/openxlsx don't choke
  is_listcol <- vapply(df, is.list, logical(1))
  if (any(is_listcol)) df <- df[ , !is_listcol, drop = FALSE]
  
  if (isTRUE(make_dirs)) dir.create(dirname(file_stem), showWarnings = FALSE, recursive = TRUE)
  
  # Round numeric columns safely
  numcols <- vapply(df, is.numeric, TRUE)
  if (any(numcols)) {
    for (nm in names(df)[numcols]) df[[nm]] <- round(df[[nm]], digits)
  }
  
  paths <- list()
  
  ## CSV
  csv_path <- paste0(file_stem, ".csv")
  utils::write.csv(df, csv_path, row.names = FALSE)
  paths$csv <- csv_path
  
  ## HTML via gt (optional)
  if (requireNamespace("gt", quietly = TRUE)) {
    gt_tbl <- gt::gt(df)
    if (!is.null(title)) gt_tbl <- gt::tab_header(gt_tbl, title = title)
    gt::gtsave(gt_tbl, paste0(file_stem, ".html"))
    paths$html <- paste0(file_stem, ".html")
  } else if (verbose) {
    message("[save_table_readable] Package 'gt' not installed → skipping HTML.")
  }
  
  ## XLSX via openxlsx (preferred) or writexl (fallback)
  xlsx_path <- paste0(file_stem, ".xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "table")
    
    # Optional title in A1, style it a bit
    start_row <- 1L
    if (!is.null(title)) {
      openxlsx::writeData(wb, "table", title, startRow = start_row, startCol = 1)
      # bold, bigger font for title
      st <- openxlsx::createStyle(textDecoration = "bold", fontSize = 14)
      openxlsx::addStyle(wb, "table", st, rows = start_row, cols = 1, gridExpand = TRUE, stack = TRUE)
      start_row <- start_row + 2L  # blank row after title
    }
    
    openxlsx::writeData(wb, "table", df, startRow = start_row, startCol = 1)
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    paths$xlsx <- xlsx_path
  } else if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(df, xlsx_path)
    paths$xlsx <- xlsx_path
  } else if (verbose) {
    message("[save_table_readable] Neither 'openxlsx' nor 'writexl' installed → skipping XLSX.")
  }
  
  invisible(paths)
}



#' Save a ggplot/patchwork safely (no-op if not a plot)
#'
#' @param p A ggplot/patchwork object.
#' @param file Output path (with extension, e.g. \code{.png}).
#' @param width,height Figure size in inches.
#' @param dpi Resolution in dots per inch.
#' @param bg Background color (default \code{"white"}).
#' @keywords io export plot
save_plot_safe <- function(p, file, width = 9, height = 6, dpi = 300, bg = "white") {
  if (inherits(p, c("gg", "ggplot", "patchwork"))) {
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    try(ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi, bg = bg),
        silent = TRUE)
  }
}
# =============================================================================

order_models_by_median_rmse <- function(cv_tbl) {
  bm <- block_metrics_long(cv_tbl)
  bm |>
    dplyr::filter(Metric == "RMSE") |>
    dplyr::group_by(model) |>
    dplyr::summarise(med = stats::median(Value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(med) |>
    dplyr::pull(model)
}

# Block-wise metrics (RMSE, MAE)
block_metrics_long <- function(cv_tbl) {
  stopifnot(all(c("model","block_id","obs","pred") %in% names(cv_tbl)))
  cv_tbl |>
    dplyr::group_by(model, block_id) |>
    dplyr::summarise(
      RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
      MAE  = mean(abs(obs - pred), na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(c(RMSE, MAE), names_to = "Metric", values_to = "Value")
}
make_block_metric_box <- function(cv_tbl, which_time = "T14", tail_cap = 0.995) {
  bm <- block_metrics_long(cv_tbl) |>
    dplyr::filter(is.finite(Value))
  if (!is.null(tail_cap)) {
    ymax <- stats::quantile(bm$Value, tail_cap, na.rm = TRUE)
  }
  lev <- order_models_by_median_rmse(cv_tbl)
  bm$model <- factor(bm$model, levels = lev)
  
  ggplot2::ggplot(bm, ggplot2::aes(model, Value)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.35, width = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                          fill = "white", colour = "black", stroke = 0.5) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Block-wise errors (LBO-CV)", which_time),
                  subtitle = "Box = IQR · line = median · ◆ = mean",
                  x = "Model", y = "Error") +
    ggplot2::facet_wrap(~ Metric, scales = "free_y")
}

make_abs_error_box <- function(cv_tbl, which_time = "T14", tail_cap = 0.995) {
  df <- cv_tbl |>
    dplyr::mutate(abs_err = abs(pred - obs)) |>
    dplyr::filter(is.finite(abs_err))
  ymax <- if (!is.null(tail_cap)) stats::quantile(df$abs_err, tail_cap, na.rm = TRUE) else max(df$abs_err, na.rm = TRUE)
  lev <- df |>
    dplyr::group_by(model) |>
    dplyr::summarise(med = stats::median(abs_err, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(med) |>
    dplyr::pull(model)
  df$model <- factor(df$model, levels = lev)
  
  ggplot2::ggplot(df, ggplot2::aes(model, abs_err)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.3, width = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                          fill = "white", colour = "black", stroke = 0.5) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Absolute errors per station (LBO-CV)", which_time),
                  subtitle = "Box = IQR · line = median · ◆ = mean",
                  x = "Model", y = "|pred − obs|")
}


make_obs_pred_scatter <- function(cv_tbl, which_time = "T14") {
  lab <- .make_labeller(cv_tbl)
  ggplot(cv_tbl, aes(obs, pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.7, shape = 16) +
    coord_equal() + theme_minimal() +
    labs(title = sprintf("%s — Observed vs Predicted (LBO-CV)", which_time), x = "Observed", y = "Predicted") +
    facet_wrap(~ model, ncol = 3, labeller = ggplot2::as_labeller(lab))
}

make_residual_density <- function(cv_tbl, which_time = "T14") {
  cv_tbl |> dplyr::mutate(resid = pred - obs) |> ggplot2::ggplot(ggplot2::aes(resid, fill = model)) +
    ggplot2::geom_density(alpha = 0.4) + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Residual density", which_time), x = "Residual (°C)", y = "Density")
}

# Prediction maps & error panels ---------------------------------------
.make_labeller <- function(cv_tbl) {
  m <- cv_tbl |>
    dplyr::group_by(model) |>
    dplyr::summarise(RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)), MAE  = mean(abs(obs - pred), na.rm = TRUE), .groups = "drop")
  setNames(sprintf("%s  (RMSE=%.2f · MAE=%.2f)", m$model, m$RMSE, m$MAE), m$model)
}
.plot_raster_gg <- function(r, title = "", palette = temp_palette, q = c(0.02,0.98), lims = NULL) {
  stopifnot(terra::nlyr(r) == 1)
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  nm <- names(df)[3]
  if (is.null(lims)) {
    vv <- terra::values(r, na.rm = TRUE)
    lims <- stats::quantile(vv, probs = q, na.rm = TRUE, names = FALSE)
  }
  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data[[nm]])) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = palette, limits = lims, oob = scales::squish) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "°C") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold"))
}

.get_preds_from_maps <- function(maps) {
  # 1) SpatRaster direkt
  if (inherits(maps, "SpatRaster")) {
    ul <- terra::unstack(maps)
    names(ul) <- names(maps)
    return(ul)
  }
  # 2) Liste mit typischen Feldern
  if (is.list(maps)) {
    if (!is.null(maps$preds))          return(maps$preds)
    if (!is.null(maps$pred_rasters))   return(maps$pred_rasters)
    if (!is.null(maps$pred_stack) && inherits(maps$pred_stack, "SpatRaster")) {
      ul <- terra::unstack(maps$pred_stack); names(ul) <- names(maps$pred_stack); return(ul)
    }
    if (!is.null(maps$stack) && inherits(maps$stack, "SpatRaster")) {
      ul <- terra::unstack(maps$stack); names(ul) <- names(maps$stack); return(ul)
    }
    if (!is.null(maps$maps) && inherits(maps$maps, "SpatRaster")) {
      ul <- terra::unstack(maps$maps); names(ul) <- names(maps$maps); return(ul)
    }
    # 3) Liste, die bereits einzelne SpatRaster oder ggplots enthält
    cand <- maps[ vapply(maps, function(x) inherits(x, "SpatRaster") || inherits(x, "ggplot"), logical(1)) ]
    if (length(cand) > 0) return(cand)
  }
  stop("build_panels_with_errors(): In 'maps' keine Vorhersagen gefunden.")
}
# --- Kartenplot mit optionalen Achsenticks/labels ----------------------
.plot_map_axes <- function(r, title, cols, lims, q = c(0.02,0.98),
                           base_size = 14, tick_n = 5,
                           show_axis_labels = FALSE, show_axis_ticks = TRUE) {
  stopifnot(terra::nlyr(r) == 1)
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  nm <- names(df)[3]
  
  if (is.null(lims) || !all(is.finite(lims)) || lims[1] >= lims[2]) {
    vv <- terra::values(r, na.rm = TRUE)
    lims <- stats::quantile(vv, probs = q, na.rm = TRUE, names = FALSE)
    if (!all(is.finite(lims)) || lims[1] == lims[2]) lims <- range(vv, na.rm = TRUE) + c(-1e-6, 1e-6)
  }
  if (is.function(cols)) cols <- cols(256)
  if (!is.atomic(cols) || length(cols) < 2) cols <- grDevices::hcl.colors(256, "YlOrRd")
  
  ggplot2::ggplot(df, ggplot2::aes(x, y, fill = .data[[nm]])) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0,0), breaks = scales::breaks_pretty(n = tick_n)) +
    ggplot2::scale_y_continuous(expand = c(0,0), breaks = scales::breaks_pretty(n = tick_n)) +
    ggplot2::scale_fill_gradientn(colours = cols, limits = lims, oob = scales::squish) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "°C") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title   = ggplot2::element_blank(),
      axis.text    = if (show_axis_labels) ggplot2::element_text(size = base_size - 3) else ggplot2::element_blank(),
      axis.ticks   = if (show_axis_ticks)  ggplot2::element_line(linewidth = 0.25) else ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey40", linewidth = .4)
    )
}
build_panels_truth_preds_errors_paged <- function(
    maps, truth_raster, cv_tbl, which_time,
    models_per_page = 4,
    model_order      = NULL,
    temp_pal         = temp_palette,     # Vektor ODER Funktion -> wird zu Vektor
    stretch_q        = c(0.02, 0.98),
    errors_height    = 1.2,
    scatter_next_to_truth = TRUE,        # Scatter rechts von Truth?
    top_widths       = c(1.1, 0.9),      # Breitenverhältnis Truth | Scatter
    show_second_legend = FALSE           # zweite °C-Legende bei den Preds unterdrücken
) {
  stopifnot(length(stretch_q) == 2)
  if (is.function(temp_pal)) temp_pal <- temp_pal(256)
  stopifnot(is.atomic(temp_pal), length(temp_pal) >= 2)
  
  # Vorhersagen einsammeln / Reihenfolge
  preds_raw  <- .get_preds_from_maps(maps)
  pred_names <- names(preds_raw) %||% paste0("model_", seq_along(preds_raw))
  if (!is.null(model_order)) {
    keep <- intersect(model_order, pred_names)
    if (!length(keep)) stop("model_order enthält keine gültigen Modellnamen.")
    preds_raw  <- preds_raw[keep]
    pred_names <- keep
  }
  
  # Gemeinsame Farbskala
  all_vals <- c(terra::values(truth_raster, na.rm = TRUE))
  for (p in preds_raw) if (inherits(p, "SpatRaster")) all_vals <- c(all_vals, terra::values(p, na.rm = TRUE))
  lims <- stats::quantile(all_vals, probs = stretch_q, na.rm = TRUE, names = FALSE)
  if (all(is.finite(lims)) && lims[1] == lims[2]) {
    eps <- .Machine$double.eps * max(1, abs(lims[1])); lims <- lims + c(-eps, eps)
  }
  
  # Helfer: Raster -> ggplot
  make_tile <- function(obj, title_txt, show_legend = TRUE) {
    if (inherits(obj, "SpatRaster")) {
      g <- .plot_raster_gg(obj, title = title_txt, palette = temp_pal, q = stretch_q, lims = lims)
      if (!show_legend) g <- g + ggplot2::theme(legend.position = "none")
      g
    } else if (inherits(obj, "ggplot")) {
      obj + ggplot2::labs(title = title_txt)
    } else stop("Nicht unterstützter Prediction-Typ: ", class(obj)[1])
  }
  
  # Truth (+ optional Scatter daneben)
  p_truth   <- .plot_raster_gg(truth_raster, title = paste0(which_time, " — truth"),
                               palette = temp_pal, q = stretch_q, lims = lims)
  p_scatter <- make_obs_pred_scatter(cv_tbl, which_time = which_time)
  
  # Prediction-Kacheln bauen (nur erste mit °C-Legende falls gewünscht)
  pred_tiles <- lapply(seq_along(preds_raw), function(i) {
    show_leg <- if (isTRUE(show_second_legend)) TRUE else (i == 1L)
    make_tile(preds_raw[[i]], pred_names[i], show_legend = show_leg)
  })
  
  # Paginierung
  n <- length(pred_tiles)
  idx_split <- split(seq_len(n), ceiling(seq_len(n) / models_per_page))
  
  pages <- lapply(idx_split, function(idx) {
    preds_row <- patchwork::wrap_plots(pred_tiles[idx], nrow = 1, ncol = length(idx))
    
    top_row <- if (isTRUE(scatter_next_to_truth)) {
      (p_truth | (p_scatter + ggplot2::theme(legend.position = "none"))) +
        patchwork::plot_layout(widths = top_widths)
    } else {
      p_truth
    }
    
    # Errors unten: wenn Scatter schon oben, unten nur Dichte
    p_box_rmse <- make_block_metric_box(cv_tbl, which_time = which_time, tail_cap = 0.995)
    p_box_ae   <- make_abs_error_box  (cv_tbl, which_time = which_time, tail_cap = 0.995)
    p_dens     <- make_residual_density(cv_tbl, which_time = which_time)
    p_errors   <- (p_box_rmse | p_box_ae) / p_dens
    
    (top_row / preds_row / p_errors) +
      patchwork::plot_layout(heights = c(1, 1, errors_height), guides = "collect") &
      ggplot2::theme(legend.position = "right")
  })
  
  pages
}
.pm_verbose <- function(v = NULL) {
  if (!is.null(v)) return(isTRUE(v))
  isTRUE(getOption("pipemodel.verbose", FALSE))
}
pm_say <- function(fmt, ..., v = NULL) {
  if (.pm_verbose(v)) message(sprintf(fmt, ...))
}
.k_for_xy <- function(n, n_xy) max(3, min(60, n_xy - 1L, floor(n * 0.8)))
.kcap_unique <- function(x, kmax) {
  ux <- unique(x[is.finite(x)])
  nu <- length(ux)
  if (nu <= 3) return(0L)                # treat as constant/near-constant
  max(4L, min(kmax, nu - 1L))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
# --- Sun helper (self-contained in the lib) --------------------------
sun_pos_utc <- function(y, m, d, h, lat, lon) {
  t  <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", y, m, d, h), tz = "UTC")
  sp <- suncalc::getSunlightPosition(date = t, lat = lat, lon = lon)
  list(
    alt = sp$altitude,
    az  = (sp$azimuth + pi) %% (2*pi)  # convert to [0, 2π) from north
  )
}

# Sun helper: pull sun14/sun05 from scen; else compute; else fallback
.get_sun <- function(scen, which = c("T14","T05")) {
  which <- match.arg(which)
  key   <- if (which == "T14") "sun14" else "sun05"
  
  # 1) stored in scen?
  s <- scen[[key]]
  if (is.list(s) && is.finite(s$alt) && is.finite(s$az)) return(s)
  
  # 2) compute from lat/lon/sun_date if available
  if (all(c("lat","lon","sun_date") %in% names(scen))) {
    hour <- if (which == "T14") 14L else 5L
    return(sun_pos_utc(scen$sun_date, hour, scen$lat, scen$lon))
  }
  
  # 3) hard fallback
  list(alt = if (which == "T14") 0.75 else 0.10, az = 0.0)
}

# -------------------------- Defaults -----------------------------------
lc_levels_default <- c("forest","water","bare soil","meadows")
lc_levels <- getOption("pipemodel.lc_levels", lc_levels_default)

lc_colors_default <- c(
  "forest"   = "#2E8B57",
  "water"    = "#5DADE2",
  "bare soil"= "#C49A6C",
  "meadows"  = "#7FBF7B"
)
temp_palette <- grDevices::colorRampPalette(c("#0000FF","#FF0000"))
stretch_q    <- c(0.02, 0.98)
# Normalize any CRS input to a non-empty character string
norm_crs <- function(crs, fallback = "EPSG:32632") {
  # allow sf::crs, numeric EPSG, character EPSG/WKT
  if (inherits(crs, "crs")) {
    out <- sf::st_crs(crs)$wkt
  } else if (is.numeric(crs) && length(crs) == 1 && is.finite(crs)) {
    out <- sprintf("EPSG:%d", as.integer(crs))
  } else if (is.character(crs) && length(crs) >= 1) {
    out <- crs[1]
  } else {
    out <- NA_character_
  }
  if (!length(out) || is.na(out) || identical(out, "")) out <- fallback
  out
}

# -------------------------- Domain/Template -----------------------------
make_domain <- function(center_E, center_N, len_x, len_y, res, crs = "EPSG:32632") {
  crs <- norm_crs(crs)
  ext <- terra::ext(center_E - len_x/2, center_E + len_x/2,
                    center_N - len_y/2, center_N + len_y/2)
  Rtemplate <- terra::rast(ext, resolution = res, crs = crs)
  list(
    xmin = terra::xmin(ext), xmax = terra::xmax(ext),
    ymin = terra::ymin(ext), ymax = terra::ymax(ext),
    x0 = center_E, y0 = center_N,
    res = as.numeric(res), crs = crs,
    Rtemplate = Rtemplate
  )
}


compute_block_size <- function(len_x, len_y, n_st,
                               target_st_per_block = 3,
                               min_blocks_axis = 3,
                               round_to = 50) {
  area <- len_x * len_y
  B_target <- max(min_blocks_axis^2, round(n_st / target_st_per_block))
  bs <- sqrt(area / B_target)
  bs <- min(bs, len_x / min_blocks_axis, len_y / min_blocks_axis)
  bs <- max(round_to, round(bs / round_to) * round_to)
  as.integer(bs)
}

# -------------------------- Sonne/Geometrie -----------------------------

# Cosine of incidence on a slope/aspect for a given sun (radians)
cosi_fun <- function(alt, az, slp_r, asp_r) {
  zen <- (pi/2 - alt)
  ci  <- cos(slp_r)*cos(zen) + sin(slp_r)*sin(zen)*cos(az - asp_r)
  terra::ifel(ci < 0, 0, ci)
}

# Sun position at a given UTC date + hour (numeric hour), return radians
sun_pos_utc <- function(sun_date, hour_utc, lat, lon) {
  stopifnot(inherits(sun_date, "Date"), length(hour_utc) == 1)
  t  <- as.POSIXct(sprintf("%s %02d:00:00", format(sun_date, "%Y-%m-%d"), as.integer(hour_utc)), tz = "UTC")
  sp <- suncalc::getSunlightPosition(date = t, lat = lat, lon = lon)
  list(
    alt = as.numeric(sp$altitude),                   # radians
    az  = as.numeric((sp$azimuth + pi) %% (2*pi))    # convert to [0..2π) from north
  )
}


# -------------------------- Rauschen ------------------------------------
make_noise_pair <- function(template, sd14 = 0.3, sd05 = 0.3,
                            seed14 = 1001, seed05 = 1002) {
  set.seed(seed14)
  n14 <- terra::setValues(terra::rast(template), rnorm(terra::ncell(template), 0, sd14))
  set.seed(seed05)
  n05 <- terra::setValues(terra::rast(template), rnorm(terra::ncell(template), 0, sd05))
  list(noise14 = n14, noise05 = n05)
}

# -------------------------- Topographie ---------------------------------
build_topography <- function(domain,
                             lake_mode = c("none","water","hollow"),
                             hill_mode = c("none","bump"),
                             lake_diam_m  = 50,  lake_depth_m = 10, smooth_edges = FALSE,
                             hill_diam_m  = 80,  hill_height_m = 50, hill_smooth  = FALSE) {
  lake_mode <- match.arg(lake_mode); hill_mode <- match.arg(hill_mode)
  Rtemplate <- domain$Rtemplate
  x0 <- domain$x0; y0 <- domain$y0
  xmin <- domain$xmin; xmax <- domain$xmax
  len_x <- xmax - xmin; y_hc <- y0
  x_hc <- xmin + len_x/3; x_lc <- xmin + 2*len_x/3; y_lc <- y0
  
  XY <- as.data.frame(terra::xyFromCell(Rtemplate, 1:terra::ncell(Rtemplate)))
  names(XY) <- c("x","y")
  dy <- XY$y - y0
  a  <- 100 / (( (domain$ymax - domain$ymin)/2 )^2)
  elev <- 500 + a * dy^2
  
  # See/Grube
  rl <- sqrt((XY$x - x_lc)^2 + (XY$y - y_lc)^2); lr <- lake_diam_m/2
  if (lake_mode %in% c("water","hollow")) {
    w_l <- if (smooth_edges) pmax(0, 1 - (rl/lr)^2) else as.numeric(rl <= lr)
    elev <- elev - lake_depth_m * w_l
  } else w_l <- 0
  
  # Haupt-Hügel
  if (hill_mode == "bump") {
    rh <- sqrt((XY$x - x_hc)^2 + (XY$y - y_hc)^2); hr <- max(1e-6, hill_diam_m/2)
    w_h <- if (hill_smooth) exp(-(rh/hr)^2) else as.numeric(rh <= hr)
    elev <- elev + hill_height_m * w_h
  } else w_h <- 0
  
  E <- Rtemplate; terra::values(E) <- elev; names(E) <- "elev"
  lakeR <- Rtemplate; terra::values(lakeR) <- if (lake_mode=="water") as.numeric(w_l>0) else 0; names(lakeR) <- "lake"
  hillW <- Rtemplate; terra::values(hillW) <- w_h; names(hillW) <- "hillW"
  
  slp  <- terra::terrain(E, v="slope",  unit="radians")
  asp  <- terra::terrain(E, v="aspect", unit="radians")
  
  list(E = E, lake = lakeR, hillW = hillW,
       slp = terra::ifel(is.na(slp), 0, slp),
       asp = terra::ifel(is.na(asp), 0, asp))
}
# --- Sun helpers (UTC) -------------------------------------------------
sun_pos_utc <- function(date, hour, lat, lon) {
  t  <- as.POSIXct(sprintf("%s %02d:00:00", as.Date(date), hour), tz = "UTC")
  sp <- suncalc::getSunlightPosition(date = t, lat = lat, lon = lon)
  # Azimut: 0 = Nord, positiv im Uhrzeigersinn
  list(alt = sp$altitude, az = (sp$azimuth + pi) %% (2*pi))
}

# -------------------------- Physikfelder --------------------------------
build_physics_fields <- function(topography, landcover,
                                 noise14, noise05,
                                 alpha_I_by_lc = c("forest" = 3.5, "water" = 1.5, "bare soil" = 6.0, "meadows" = 4.0),
                                 shade_fac_by_lc = c("forest" = 0.60, "water" = 1.00, "bare soil" = 1.00, "meadows" = 0.95),
                                 dawn_bias_by_lc = c("forest" = 0.30, "water" = 1.20, "bare soil" = -0.50, "meadows" = 0.05),
                                 pool_fac_by_lc  = c("forest" = 0.70, "water" = 0.80, "bare soil" = 1.10, "meadows" = 1.05),
                                 pool_block_gain = 0.4,
                                 sun14 = list(alt = 0.75, az = 0.0),
                                 sun05 = list(alt = 0.10, az = 0.0))
                                 {
  E    <- topography$E
  slp0 <- topography$slp
  asp0 <- topography$asp
  hillW<- topography$hillW
  
  # Sonnen-Inzidenz
  I14 <- cosi_fun(sun14$alt, sun14$az, slp0, asp0)
  I05 <- cosi_fun(sun05$alt, sun05$az, slp0, asp0)
  
  lc <- if (inherits(landcover, "SpatRaster")) landcover else landcover$lc
  stopifnot(inherits(lc, "SpatRaster"))
  
  v <- as.integer(terra::values(lc))
  v[!is.finite(v)] <- 1L
  v <- pmax(1L, pmin(v, length(lc_levels_default)))
  lc_char <- factor(lc_levels_default[v], levels = lc_levels_default)
  
  to_r <- function(x) terra::setValues(terra::rast(E), x)
  alpha_I <- to_r(as.numeric(alpha_I_by_lc[lc_char]))
  shade_f <- to_r(as.numeric(shade_fac_by_lc[lc_char]))
  dawn_b  <- to_r(as.numeric(dawn_bias_by_lc[lc_char]))
  pool_f  <- to_r(as.numeric(pool_fac_by_lc[lc_char]))
  
  I14_eff <- I14 * shade_f
  
  E_mean <- terra::global(E, "mean", na.rm = TRUE)[1,1]
  Y <- terra::init(E, "y")
  dist2ax <- abs(Y - (terra::ymax(E)+terra::ymin(E))/2)
  w_pool <- 70
  pool_base <- 4.0 * exp(- (dist2ax / w_pool)^2)
  pool_mod  <- pool_base * (1 - pool_block_gain * hillW) * pool_f
  
  T0_14 <- 26.0; lapse_14 <- -0.0065
  R14 <- T0_14 + lapse_14 * (E - E_mean) + alpha_I * I14_eff + noise14; names(R14) <- "T14"
  
  T0_05 <- 8.5; inv_05 <- 0.003; eta_slope <- 0.6
  R05 <- T0_05 + inv_05 * (E - E_mean) + eta_slope * slp0 - pool_mod + dawn_b + noise05; names(R05) <- "T05"
  
  list(R14 = R14, R05 = R05, I14 = I14, I05 = I05)
}

# --- Sun + cos(i) helpers (safe to keep once in your lib) --------------------
cosi_fun <- function(alt, az, slp_r, asp_r) {
  zen <- (pi/2 - alt)
  ci  <- cos(slp_r)*cos(zen) + sin(slp_r)*sin(zen)*cos(az - asp_r)
  terra::ifel(ci < 0, 0, ci)
}

sun_pos_utc <- function(sun_date, hour_utc, lat, lon) {
  stopifnot(inherits(sun_date, "Date"))
  t  <- as.POSIXct(sprintf("%s %02d:00:00",
                           format(sun_date, "%Y-%m-%d"),
                           as.integer(hour_utc)), tz = "UTC")
  sp <- suncalc::getSunlightPosition(date = t, lat = lat, lon = lon)
  list(
    alt = as.numeric(sp$altitude),                  # radians
    az  = as.numeric((sp$azimuth + pi) %% (2*pi))   # [0..2π) from north
  )
}

build_scenario <- function(
    domain,
    lake_mode = c("none","water","hollow"),
    hill_mode = c("none","bump"),
    # main hill / lake geometry (meters)
    lake_diam_m  = 50,  lake_depth_m = 10, smooth_edges = FALSE,
    hill_diam_m  = 80,  hill_height_m = 50, hill_smooth  = FALSE,
    # micro-relief (meters)
    random_hills        = 0,
    micro_hill_diam_m   = 30,
    micro_hill_height_m = 50,
    micro_hill_smooth   = TRUE,
    micro_seed          = NULL,
    # sun / geo
    lat = 51.8, lon = 10.6, sun_date = as.Date("2024-06-21"),
    # optional noise
    noise14 = NULL, noise05 = NULL
) {
  lake_mode <- match.arg(lake_mode)
  hill_mode <- match.arg(hill_mode)
  
  # --- 0) Template & CRS guard (must be meters) -----------------------
  ext <- terra::ext(domain$xmin, domain$xmax, domain$ymin, domain$ymax)
  Rtemplate <- terra::rast(ext, resolution = domain$res, crs = domain$crs)
  
  crs_sf <- sf::st_crs(terra::crs(Rtemplate, proj=TRUE))
  if (isTRUE(sf::st_is_longlat(crs_sf))) {
    stop(
      "build_scenario(): Domain CRS is geographic (degrees). ",
      "All geometry is in meters. Use a projected CRS (e.g. UTM / EPSG:32632)."
    )
  }
  
  xmin <- terra::xmin(ext); xmax <- terra::xmax(ext)
  ymin <- terra::ymin(ext); ymax <- terra::ymax(ext)
  len_x <- xmax - xmin;     len_y <- ymax - ymin
  x0 <- (xmin + xmax)/2;    y0 <- (ymin + ymax)/2
  
  # coordinate rasters
  X <- terra::init(Rtemplate, "x")
  Y <- terra::init(Rtemplate, "y")
  
  # quick sanity for scale
  px <- mean(terra::res(Rtemplate))
  lr_px <- (lake_diam_m/2) / px
  hr_px <- (hill_diam_m/2) / px
  message(sprintf("[build_scenario] pixel=%.2f m; lake r=%.1f px; hill r=%.1f px", px, lr_px, hr_px))
  
  # --- 1) Base valley --------------------------------------------------
  a  <- 100 / ((len_y/2)^2)
  E  <- 500 + a * (Y - y0)^2
  names(E) <- "elev"
  
  # --- 2) Lake (mirror of hill, negative) ------------------------------
  x_lc <- xmin + 2*len_x/3;  y_lc <- y0
  lr   <- max(1e-6, lake_diam_m/2)
  rl   <- sqrt((X - x_lc)^2 + (Y - y_lc)^2)
  
  w_l <- if (isTRUE(smooth_edges)) {
    exp(-(rl/lr)^2)            # Gaussian "bump"
  } else {
    terra::ifel(rl <= lr, 1, 0) # hard disc
  }
  
  if (lake_mode %in% c("water","hollow")) {
    E <- E - as.numeric(lake_depth_m) * w_l
  }
  lakeR <- if (identical(lake_mode, "water")) terra::ifel(w_l > 1e-6, 1L, 0L)
  else terra::setValues(terra::rast(Rtemplate), 0L)
  names(lakeR) <- "lake"
  
  # --- 3) Main hill ----------------------------------------------------
  x_hc <- xmin + len_x/3;  y_hc <- y0
  hr   <- max(1e-6, hill_diam_m/2)
  rh   <- sqrt((X - x_hc)^2 + (Y - y_hc)^2)
  
  w_h_main <- if (hill_mode == "bump") {
    if (isTRUE(hill_smooth)) exp(-(rh/hr)^2) else terra::ifel(rh <= hr, 1, 0)
  } else {
    0 * X
  }
  E <- E + as.numeric(hill_height_m) * w_h_main
  
  # --- 4) Micro hills (additive, clamped to 1) ------------------------
  w_h_micro <- 0 * X
  if (random_hills > 0) {
    if (!is.null(micro_seed)) set.seed(micro_seed)
    margin <- micro_hill_diam_m/2 + 5
    hrm <- max(1e-6, micro_hill_diam_m/2)
    for (i in seq_len(random_hills)) {
      cx <- runif(1, xmin + margin, xmax - margin)
      cy <- runif(1, ymin + margin, ymax - margin)
      r  <- sqrt((X - cx)^2 + (Y - cy)^2)
      wi <- if (isTRUE(micro_hill_smooth)) exp(-(r/hrm)^2) else terra::ifel(r <= hrm, 1, 0)
      sum_i <- w_h_micro + wi
      w_h_micro <- terra::ifel(sum_i > 1, 1, sum_i)  # clamp without pmin()
    }
    E <- E + as.numeric(micro_hill_height_m) * w_h_micro
  }
  
  hillW <- w_h_main + w_h_micro
  hillW <- terra::ifel(hillW > 1, 1, hillW); names(hillW) <- "hillW"
  
  # --- 5) Derivatives --------------------------------------------------
  slp <- terra::terrain(E, v = "slope",  unit = "radians")
  asp <- terra::terrain(E, v = "aspect", unit = "radians")
  
  # --- 6) Sun & cos(i) -------------------------------------------------
  sun14 <- sun_pos_utc(sun_date, 14L, lat, lon)
  sun05 <- sun_pos_utc(sun_date,  5L, lat, lon)
  I14   <- cosi_fun(sun14$alt, sun14$az, slp, asp); names(I14) <- "I14"
  I05   <- cosi_fun(sun05$alt, sun05$az, slp, asp); names(I05) <- "I05"
  
  # --- 7) Land cover (1 forest, 2 water, 3 bare, 4 meadows) -----------
  lc <- terra::setValues(terra::rast(Rtemplate), 4L)  # meadows
  lc <- terra::ifel(lakeR > 0, 2L, lc)                # water overrides
  forest_mask <- terra::ifel((hillW > 0.2) | ((slp > 0.15) & (Y > y0)), 1, 0)
  lc <- terra::ifel((forest_mask == 1) & (lakeR <= 0), 1L, lc)
  v_slp   <- terra::values(slp)
  thr_slp <- stats::quantile(v_slp[is.finite(v_slp)], 0.90, na.rm = TRUE)
  bare_mask <- terra::ifel((slp >= thr_slp) & (lakeR <= 0) & (forest_mask == 0), 1, 0)
  lc <- terra::ifel(bare_mask == 1, 3L, lc)
  lc <- terra::clamp(lc, 1L, 4L); names(lc) <- "lc"
  
  lc_levels <- c("forest","water","bare soil","meadows")
  lc_colors <- c("forest"="#2E8B57","water"="#5DADE2","bare soil"="#C49A6C","meadows"="#7FBF7B")
  
  # --- 8) Noise --------------------------------------------------------
  if (is.null(noise14)) {
    set.seed(1001)
    noise14 <- terra::setValues(terra::rast(E), rnorm(terra::ncell(E), 0, 0.3))
  }
  if (is.null(noise05)) {
    set.seed(1002)
    noise05 <- terra::setValues(terra::rast(E), rnorm(terra::ncell(E), 0, 0.3))
  }
  
  # --- 9) Physics fields ----------------------------------------------
  topo <- list(E = E, slp = slp, asp = asp, hillW = hillW)
  phys <- build_physics_fields(
    topography = topo, landcover = lc,
    noise14 = noise14, noise05 = noise05,
    sun14 = sun14, sun05 = sun05
  )
  R14 <- phys$R14; R05 <- phys$R05
  
  # --- 10) Return ------------------------------------------------------
  list(
    E = E, slp = slp, asp = asp,
    I14 = I14, I05 = I05,
    R14 = R14, R05 = R05,
    lake = lakeR, hillW = hillW,
    lc = lc, lc_levels = lc_levels, lc_colors = lc_colors,
    sun = list(T14 = sun14, T05 = sun05)
  )
}




# -------------------------- Stationen/Features -------------------------
make_stations <- function(domain, n_st = 60,
                          station_mode = c("random","ns_transect","ew_transect"),
                          transect_margin_m = 10, ns_offset_m = 0, ew_offset_m = 0,
                          crs = sf::st_crs(domain$Rtemplate)) {
  station_mode <- match.arg(station_mode)
  with(domain, {
    if (station_mode == "random") {
      pts <- tibble::tibble(
        id = 1:n_st,
        x  = runif(n_st, xmin + transect_margin_m, xmax - transect_margin_m),
        y  = runif(n_st, ymin + transect_margin_m, ymax - transect_margin_m)
      )
    } else if (station_mode == "ns_transect") {
      x_const <- min(max(x0 + ns_offset_m, xmin + transect_margin_m), xmax - transect_margin_m)
      y_seq   <- seq(ymin + transect_margin_m, ymax - transect_margin_m, length.out = n_st)
      pts <- tibble::tibble(id = 1:n_st, x = x_const, y = y_seq)
    } else {
      y_const <- min(max(y0 + ew_offset_m, ymin + transect_margin_m), ymax - transect_margin_m)
      x_seq   <- seq(xmin + transect_margin_m, xmax - transect_margin_m, length.out = n_st)
      pts <- tibble::tibble(id = 1:n_st, x = x_seq, y = y_const)
    }
    sf::st_as_sf(pts, coords = c("x","y"), crs = crs, remove = FALSE)
  })
}

stations_from_scenario <- function(scen, pts_sf) {
  vpts <- terra::vect(pts_sf)
  df <- tibble::as_tibble(pts_sf) %>%
    dplyr::mutate(
      z_surf = as.numeric(terra::extract(scen$E,   vpts, ID = FALSE)[,1]),
      slp    = as.numeric(terra::extract(scen$slp, vpts, ID = FALSE)[,1]),
      I14    = as.numeric(terra::extract(scen$I14, vpts, ID = FALSE)[,1]),
      I05    = as.numeric(terra::extract(scen$I05, vpts, ID = FALSE)[,1]),
      lc     = as.integer(terra::extract(scen$lc,  vpts, ID = FALSE)[,1]),
      T14    = as.numeric(terra::extract(scen$R14, vpts, ID = FALSE)[,1]),
      T05    = as.numeric(terra::extract(scen$R05, vpts, ID = FALSE)[,1])
    )
  lc_levels <- scen$lc_levels
  pts14 <- df[stats::complete.cases(df[, c("x","y","z_surf","slp","I14","lc","T14")]), ]
  pts05 <- df[stats::complete.cases(df[, c("x","y","z_surf","slp","I05","lc","T05")]), ]
  stn_sf_14 <- pts14 %>%
    dplyr::transmute(id, x, y,
                     z_surf = as.numeric(z_surf), slp = as.numeric(slp), cosi = as.numeric(I14),
                     lc = factor(lc_levels[pmax(1, pmin(lc, length(lc_levels)))], levels = lc_levels),
                     temp = as.numeric(T14)) %>%
    sf::st_as_sf(coords = c("x","y"), crs = sf::st_crs(pts_sf), remove = FALSE)
  stn_sf_05 <- pts05 %>%
    dplyr::transmute(id, x, y,
                     z_surf = as.numeric(z_surf), slp = as.numeric(slp), cosi = as.numeric(I05),
                     lc = factor(lc_levels[pmax(1, pmin(lc, length(lc_levels)))], levels = lc_levels),
                     temp = as.numeric(T05)) %>%
    sf::st_as_sf(coords = c("x","y"), crs = sf::st_crs(pts_sf), remove = FALSE)
  list(T14 = stn_sf_14, T05 = stn_sf_05)
}

# -------------------------- Plots: Übersicht ---------------------------

# -------------------------- Preview: Domain ----------------------------
# Zeigt Extent, optional ein Grid, und annotiert Kern-Parameter.
preview_domain <- function(domain, grid = TRUE, grid_step = NULL, annotate = TRUE) {
  stopifnot(is.list(domain), !is.null(domain$Rtemplate))
  crs <- sf::st_crs(domain$Rtemplate)
  bb  <- sf::st_as_sfc(sf::st_bbox(c(
    xmin = domain$xmin, ymin = domain$ymin,
    xmax = domain$xmax, ymax = domain$ymax
  ), crs = crs))

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = bb, fill = NA, color = "black", linewidth = 0.7) +
    ggplot2::coord_sf(crs = crs, datum = NA, expand = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Domain preview", x = "Easting (m)", y = "Northing (m)")

  if (isTRUE(grid)) {
    len_x <- domain$xmax - domain$xmin
    len_y <- domain$ymax - domain$ymin
    step  <- if (is.null(grid_step)) max(100, round(min(len_x, len_y)/6, -2)) else grid_step
    gr    <- sf::st_make_grid(bb, cellsize = c(step, step), what = "polygons")
    gr_ln <- sf::st_as_sf(sf::st_boundary(gr))
    p <- p + ggplot2::geom_sf(data = gr_ln, color = "grey70", linewidth = 0.2)
  }

  if (isTRUE(annotate)) {
    r   <- terra::res(domain$Rtemplate)
    txt <- sprintf(
      "Center: (%.0f, %.0f)\nSize: %.0f × %.0f m\nRes: %.0f × %.0f m\nCRS: %s",
      domain$x0, domain$y0,
      domain$xmax - domain$xmin, domain$ymax - domain$ymin,
      r[1], r[2], as.character(terra::crs(domain$Rtemplate))
    )
    p <- p + ggplot2::annotate(
      "label",
      x = domain$xmin + 0.02 * (domain$xmax - domain$xmin),
      y = domain$ymin + 0.06 * (domain$ymax - domain$ymin),
      label = txt, hjust = 0, vjust = 0, size = 3
    )
  }

  print(p)
  invisible(p)
}

# -------------------------- Preview: Szenario-Raster -------------------
# Visualisiert die vorhandenen Raster im Szenario (ohne Modelle/CV).
preview_scenario <- function(x,
                             which = c("lc","E","slp","R14","R05"),
                             stations = NULL,
                             show_contours = TRUE,
                             layout = c("grid","vertical")) {
  layout <- match.arg(layout)
  
  # ---- Szenario + Stationen erkennen ---------------------------------
  scen <- if (is.list(x) && !is.null(x$scen)) x$scen else x
  if (is.null(stations) && is.list(x)) {
    stations <- x$pts_sf %||% x$stn_sf_14 %||% x$stn_sf_05
  }
  
  # ---- Verfügbare Ebenen sammeln -------------------------------------
  layers <- list(
    lc  = scen$lc,
    E   = scen$E,
    slp = scen$slp,
    R14 = scen$R14,
    R05 = scen$R05
  )
  keep <- names(layers) %in% which & vapply(layers, function(r) inherits(r,"SpatRaster"), TRUE)
  layers <- layers[keep]
  if (!length(layers)) stop("preview_scenario(): Keine der angefragten Ebenen im Szenario vorhanden.")
  
  # ---- optionale Konturen vorbereiten --------------------------------
  add_contours <- function(p) p
  if (isTRUE(show_contours)) {
    adders <- list()
    
    has_lake <- !is.null(scen$lake) && inherits(scen$lake, "SpatRaster")
    if (has_lake) {
      lake_df <- as.data.frame(scen$lake, xy = TRUE); names(lake_df) <- c("x","y","lake")
      adders[[length(adders)+1]] <- ggplot2::geom_contour(
        data = lake_df,
        mapping = ggplot2::aes(x = x, y = y, z = lake),
        breaks = 0.5, colour = "black", linewidth = 0.35,
        inherit.aes = FALSE
      )
    }
    
    has_hill <- !is.null(scen$hillW) && inherits(scen$hillW, "SpatRaster")
    if (has_hill) {
      hill_df <- as.data.frame(scen$hillW, xy = TRUE); names(hill_df) <- c("x","y","hillW")
      adders[[length(adders)+1]] <- ggplot2::geom_contour(
        data = hill_df,
        mapping = ggplot2::aes(x = x, y = y, z = hillW),
        breaks = 0.5, colour = "black", linetype = "22", linewidth = 0.3,
        inherit.aes = FALSE
      )
    }
    
    add_contours <- function(p) {
      if (length(adders)) for (a in adders) p <- p + a
      p
    }
  }
  
  # ---- Einzelplots bauen ----------------------------------------------
  make_plot <- function(name, r) {
    df <- as.data.frame(r, xy = TRUE); names(df) <- c("x","y","val")
    
    if (name == "lc") {
      levs <- scen$lc_levels %||% c("1","2","3","4")
      pal  <- scen$lc_colors %||% setNames(scales::hue_pal()(length(levs)), levs)
      df$val <- factor(df$val, levels = seq_along(levs), labels = levs)
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = val)) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_manual(
          values = pal[levels(df$val)],  # <-- sicher auf Levels abbilden
          na.value = "grey90", name = "Landuse"
        ) +
        ggplot2::coord_equal() + ggplot2::theme_minimal() +
        ggplot2::labs(title = "Landuse", x = "Easting", y = "Northing")
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = val)) +
        ggplot2::geom_raster() +
        ggplot2::coord_equal() + ggplot2::theme_minimal() +
        ggplot2::labs(
          title = switch(name,
                         E   = "Elevation (m)",
                         slp = "Slope (rad)",
                         R14 = "T 14 UTC",
                         R05 = "T 05 UTC",
                         name),
          x = "Easting", y = "Northing"
        )
      if (name %in% c("E","slp")) {
        p <- p + ggplot2::scale_fill_viridis_c()
      } else {
        p <- p + ggplot2::scale_fill_gradientn(colors = temp_palette(256), name = "Temp")
      }
    }
    
    # Konturen (nur falls vorhanden)
    p <- add_contours(p)
    
    # Stationen optional
    if (!is.null(stations) && inherits(stations, "sf")) {
      # stilles CRS-Align (Fehlermeldungen unterdrücken)
      suppressWarnings({
        stations_plot <- try(sf::st_transform(stations, sf::st_crs(scen$lc %||% scen$E)), silent = TRUE)
        if (inherits(stations_plot, "try-error")) stations_plot <- stations
        p <- p + ggplot2::geom_sf(
          data = stations_plot, colour = "black", fill = "white",
          shape = 21, size = 1.6, stroke = 0.25, inherit.aes = FALSE
        )
      })
    }
    p
  }
  
  plots <- Map(make_plot, names(layers), layers)
  
  # ---- kombinieren ----------------------------------------------------
  if (length(plots) == 1) {
    p_out <- plots[[1]]
  } else if (layout == "vertical") {
    p_out <- patchwork::wrap_plots(plots, ncol = 1)
  } else {
    p_out <- patchwork::wrap_plots(plots, ncol = min(3, length(plots)))
  }
  
  print(p_out)
  invisible(p_out)
}




plot_landcover_terrain <- function(scen, stations = NULL, show_contours = TRUE,
                                   layout = c("grid","vertical")) {
  layout <- match.arg(layout)
  lc_df  <- as.data.frame(scen$lc,  xy = TRUE); names(lc_df)  <- c("x","y","lc")
  E_df   <- as.data.frame(scen$E,   xy = TRUE); names(E_df)   <- c("x","y","elev")
  slp_df <- as.data.frame(scen$slp, xy = TRUE); names(slp_df) <- c("x","y","slp")
  lc_df$lc <- factor(lc_df$lc, levels = seq_along(scen$lc_levels), labels = scen$lc_levels)
  
  p_lc <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = lc_df, ggplot2::aes(x, y, fill = lc)) +
    ggplot2::scale_fill_manual(values = scen$lc_colors, na.value = "grey90", name = "Landuse") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = "Landuse", x = "Easting", y = "Northing")
  
  p_elev <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = E_df, ggplot2::aes(x, y, fill = elev)) +
    ggplot2::scale_fill_viridis_c(name = "Altitude [m]") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = "Altitude", x = "Easting", y = "Northing")
  
  p_slp <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = slp_df, ggplot2::aes(x, y, fill = slp)) +
    ggplot2::scale_fill_viridis_c(name = "Slope [rad]") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = "Slope", x = "Easting", y = "Northing")
  
  if (isTRUE(show_contours)) {
    lake_df <- as.data.frame(scen$lake, xy = TRUE); names(lake_df) <- c("x","y","lake")
    hill_df <- as.data.frame(scen$hillW, xy = TRUE); names(hill_df) <- c("x","y","hillW")
    p_lc  <- p_lc  +
      ggplot2::geom_contour(data = lake_df, ggplot2::aes(x, y, z = lake),
                            breaks = 0.5, colour = "black", linewidth = 0.35) +
      ggplot2::geom_contour(data = hill_df, ggplot2::aes(x, y, z = hillW),
                            breaks = 0.5, colour = "black", linetype = "22", linewidth = 0.3)
    p_slp <- p_slp +
      ggplot2::geom_contour(data = lake_df, ggplot2::aes(x, y, z = lake),
                            breaks = 0.5, colour = "black", linewidth = 0.35) +
      ggplot2::geom_contour(data = hill_df, ggplot2::aes(x, y, z = hillW),
                            breaks = 0.5, colour = "black", linetype = "22", linewidth = 0.3)
  }
  if (!is.null(stations)) {
    add_st <- list(ggplot2::geom_sf(data = stations, colour = "black", fill = "white",
                                    shape = 21, size = 1.6, stroke = 0.25, inherit.aes = FALSE))
    p_lc   <- p_lc   + add_st
    p_elev <- p_elev + add_st
    p_slp  <- p_slp  + add_st
  }
  if (layout == "vertical") {
    (p_lc / p_elev / p_slp) + patchwork::plot_layout(guides = "keep")
  } else {
    (p_lc | (p_elev | p_slp)) + patchwork::plot_layout(guides = "keep")
  }
}

plot_block_overview_2x2_en <- function(scen, pts_sf = NULL) {
  Rstack <- c(scen$E, scen$slp, scen$I14, scen$I05)
  df <- terra::as.data.frame(Rstack, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x","y","elev","slope","I14","I05")
  theme_base <- ggplot2::theme_minimal(base_size = 11)
  pal_terrain <- grDevices::hcl.colors(256, "Terrain")
  pal_slope   <- grDevices::hcl.colors(256, "Viridis")
  pal_hot     <- grDevices::hcl.colors(256, "YlOrRd")
  pal_cool    <- grDevices::hcl.colors(256, "PuBuGn")
  p_elev <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = elev)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_terrain, name = "m") +
    ggplot2::labs(title = "Terrain (Elevation)") + theme_base
  p_slope <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = slope)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_slope, name = "rad") +
    ggplot2::labs(title = "Slope (radians)") + theme_base
  p_I14 <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = I14)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_hot, name = "") +
    ggplot2::labs(title = "Insolation 14 UTC (cos i)") + theme_base
  p_I05 <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = I05)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(colours = pal_cool, name = "") +
    ggplot2::labs(title = "Insolation 05 UTC (cos i)") + theme_base
  if (!is.null(pts_sf)) {
    pts_df <- sf::st_drop_geometry(pts_sf)
    add_pts <- function(p)
      p + ggplot2::geom_point(data = pts_df, ggplot2::aes(x = x, y = y),
                              inherit.aes = FALSE, size = 0.7, alpha = 0.7,
                              colour = "black")
    p_elev  <- add_pts(p_elev); p_slope <- add_pts(p_slope)
    p_I14   <- add_pts(p_I14);  p_I05   <- add_pts(p_I05)
  }
  (p_elev | (p_slope)) / (p_I14 | p_I05) + patchwork::plot_layout(guides = "collect")
}

# -------------------------- Geostat/Models -----------------------------
# helpers
.align_factor_to_model <- function(x, lev_model) {
  xs <- as.character(x)
  if (!length(lev_model)) return(factor(rep(NA_character_, length(xs))))
  y <- factor(xs, levels = lev_model)
  if (anyNA(y)) { xs[is.na(y)] <- lev_model[1]; y <- factor(xs, levels = lev_model) }
  y
}
.default_vgm <- function(values, model = "Exp", range = 100) {
  psill <- stats::var(values, na.rm = TRUE); nug <- 0.1 * psill
  gstat::vgm(psill = psill, model = model, range = range, nugget = nug)
}
safe_r2 <- function(obs, pred) {
  idx <- is.finite(obs) & is.finite(pred)
  if (sum(idx) < 2) return(NA_real_)
  x <- obs[idx]; y <- pred[idx]
  sx <- stats::sd(x); sy <- stats::sd(y)
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(NA_real_)
  stats::cor(x, y)^2
}
safe_gam_formula <- function(d, include_lc = FALSE) {
  stopifnot(all(c("temp","x","y") %in% names(d)))
  d <- d[stats::complete.cases(d[, c("temp","x","y")]), , drop = FALSE]
  n    <- nrow(d)
  n_xy <- dplyr::n_distinct(paste0(round(d$x,3), "_", round(d$y,3)))
  k_xy <- max(3, min(60, n_xy - 1L, floor(n * 0.8)))
  base <- if (n_xy >= 4) sprintf("temp ~ s(x,y,bs='tp',k=%d)", k_xy) else "temp ~ x + y"
  add <- character(0)
  kcap <- function(x, kmax) {
    ux <- unique(x[is.finite(x)]); nu <- length(ux)
    if (nu <= 3) return(0L); max(4L, min(kmax, nu - 1L))
  }
  if ("z_surf" %in% names(d) && dplyr::n_distinct(d$z_surf) > 3) add <- c(add, sprintf("s(z_surf,bs='tp',k=%d)", kcap(d$z_surf, 20)))
  if ("slp"    %in% names(d) && dplyr::n_distinct(d$slp)    > 3) add <- c(add, sprintf("s(slp,bs='tp',k=%d)",    kcap(d$slp, 12)))
  if ("cosi"   %in% names(d) && dplyr::n_distinct(d$cosi)   > 3) add <- c(add, sprintf("s(cosi,bs='tp',k=%d)",   kcap(d$cosi, 12)))
  if (include_lc && "lc" %in% names(d)) { d$lc <- droplevels(factor(d$lc)); if (nlevels(d$lc) >= 2) add <- c(add, "lc") }
  stats::as.formula(paste(base, paste(add, collapse = " + "), sep = if (length(add)) " + " else ""))
}
# learners
# NOTE:
# The following learner functions have been moved to a dedicated file
# (e.g., learners_geostat_core.R):
#   - pred_Voronoi
#   - pred_IDW
#   - pred_OK
#   - pred_KED
#   - pred_RF
#   - pred_GAM
#
# Source that file alongside your helpers BEFORE any code that calls them.
# -------------------------- Block-CV -----------------------------------
make_blocks_and_assign <- function(pts_sf, E, block_size = 100) {
  bb <- sf::st_as_sfc(sf::st_bbox(c(xmin = terra::xmin(E), ymin = terra::ymin(E),
                                    xmax = terra::xmax(E), ymax = terra::ymax(E)),
                                  crs = sf::st_crs(pts_sf)))
  gr <- sf::st_make_grid(bb, cellsize = c(block_size, block_size), what = "polygons")
  blocks <- sf::st_sf(block_id = seq_along(gr), geometry = gr)
  pts_blk <- sf::st_join(pts_sf, blocks, join = sf::st_intersects, left = TRUE)
  if (any(is.na(pts_blk$block_id))) {
    nearest <- sf::st_nearest_feature(pts_blk[is.na(pts_blk$block_id), ], blocks)
    pts_blk$block_id[is.na(pts_blk$block_id)] <- blocks$block_id[nearest]
  }
  list(blocks = blocks, pts = pts_blk)
}
plot_blocks_grid <- function(blocks, pts_blk, title = "Blocks & stations") {
  crs_plot <- sf::st_crs(pts_blk)
  bb       <- sf::st_bbox(blocks)
  n_blocks <- dplyr::n_distinct(pts_blk$block_id)
  cols     <- scales::hue_pal()(max(1, n_blocks))
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = blocks, fill = NA, color = "grey50", linewidth = 0.25) +
    ggplot2::geom_sf(data = pts_blk, ggplot2::aes(color = factor(block_id)), size = 2, alpha = 0.95) +
    ggplot2::scale_color_manual(values = cols, name = "Block") +
    ggplot2::coord_sf(crs  = crs_plot, datum = NA,
                      xlim = c(bb["xmin"], bb["xmax"]),
                      ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
    ggplot2::theme_minimal() + ggplot2::labs(title = title, x = "Easting (m)", y = "Northing (m)")
}
run_lbo_cv <- function(stn_sf, E, block_size = 100, models = c("Voronoi","IDW","OK","KED","RF","GAM")) {
  if (!all(c("x","y") %in% names(stn_sf))) { xy <- sf::st_coordinates(stn_sf); stn_sf$x <- xy[,1]; stn_sf$y <- xy[,2] }
  blk <- make_blocks_and_assign(stn_sf, E, block_size = block_size)
  blocks_sf <- blk$blocks; stn_blk <- blk$pts
  for (nm in c("temp","z_surf","slp","cosi","lc","x","y")) if (!(nm %in% names(stn_blk))) stn_blk[[nm]] <- stn_sf[[nm]][match(stn_blk$id, stn_sf$id)]
  
  block_ids <- sort(unique(stn_blk$block_id))
  out_list <- vector("list", length(block_ids))
  for (k in seq_along(block_ids)) {
    b <- block_ids[k]
    test_idx  <- which(stn_blk$block_id == b)
    train_idx <- which(stn_blk$block_id != b)
    train_sf <- stn_blk[train_idx, ]; test_sf <- stn_blk[test_idx, ]
    pred_tbl <- lapply(models, function(m) {
      p <- switch(m,
                  "Voronoi" = pred_Voronoi(train_sf, test_sf),
                  "IDW"     = pred_IDW    (train_sf, test_sf),
                  "OK"      = pred_OK     (train_sf, test_sf),
                  "KED"     = pred_KED    (train_sf, test_sf, E = E),
                  "RF"      = pred_RF     (train_sf, test_sf),
                  "GAM"     = pred_GAM    (train_sf, test_sf),
                  stop("Unknown model: ", m)
      )
      tibble::tibble(model = m, id = test_sf$id, obs = test_sf$temp, pred = p, block_id = b)
    })
    out_list[[k]] <- dplyr::bind_rows(pred_tbl)
  }
  
  cv_tbl <- dplyr::bind_rows(out_list)
  metrics <- cv_tbl %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      n_ok = sum(is.finite(obs) & is.finite(pred)),
      MAE  = {i <- is.finite(obs) & is.finite(pred); if (any(i)) mean(abs(pred[i]-obs[i])) else NA_real_},
      RMSE = {i <- is.finite(obs) & is.finite(pred); if (any(i)) sqrt(mean((pred[i]-obs[i])^2)) else NA_real_},
      Bias = {i <- is.finite(obs) & is.finite(pred); if (any(i)) mean(pred[i]-obs[i]) else NA_real_},
      R2   = safe_r2(obs, pred),
      .groups = "drop"
    ) |>
    dplyr::arrange(RMSE)
  
  diag_plot <- ggplot2::ggplot(cv_tbl, ggplot2::aes(obs, pred)) +
    ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed") +
    ggplot2::geom_point(alpha=0.7) +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("LBO-CV (block = %dm) — Observed vs Predicted", block_size), x = "Observed", y = "Predicted") +
    ggplot2::facet_wrap(~ model)
  
  blocks_plot <- plot_blocks_grid(blocks_sf, stn_blk, title = sprintf("Blocks (%.0f m) & stations", block_size))
  list(cv = cv_tbl, metrics = metrics, diag_plot = diag_plot, blocks_plot = blocks_plot)
}

# -------------------------- „run_for_time“ Wrapper ---------------------
run_for_time <- function(stn_sf, truth_r, label,
                         scen_local,
                         block_m,
                         models = c("Voronoi","IDW","OK","KED","RF","GAM"),
                         layout = c("horizontal","vertical")) {
  layout <- match.arg(layout)
  res   <- run_lbo_cv(stn_sf, scen_local$E, block_size = block_m, models = models)
  maps  <- predict_maps(stn_sf, truth_r, which_time = label,
                        scen = scen_local, models = models,
                        lc_levels = scen_local$lc_levels)
  list(res = res, maps = maps)
}

# -------------------------- Skalen & Tuning ----------------------------
plot_variogram_with_scales <- function(vg, L50, L95, sill, title = "Empirical variogram") {
  df <- as.data.frame(vg)
  ggplot2::ggplot(df, ggplot2::aes(dist, gamma)) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = sill, linetype = "dotted", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = L50, colour = "#2b8cbe", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = L95, colour = "#de2d26", linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Distance (m)", y = "Semivariance")
}
.mean_kernel_for_R <- function(r, R_m) {
  px <- mean(terra::res(r))
  half <- max(1L, ceiling(R_m / px))
  k <- 2L * half + 1L
  W <- matrix(1, nrow = k, ncol = k)
  W / sum(W)
}
smooth_mean_R <- function(r, R_m) {
  W <- .mean_kernel_for_R(r, R_m)
  terra::focal(r, w = W, fun = "mean", na.policy = "omit", pad = TRUE, normalize = FALSE)
}
gaussian_focal <- function(r, radius_m, sigma_m = NULL) {
  resx <- terra::res(r)[1]
  if (is.null(sigma_m)) sigma_m <- radius_m / 2
  rad_px   <- max(1L, round(radius_m / resx))
  sigma_px <- max(0.5, sigma_m / resx)
  xs <- -rad_px:rad_px
  k1 <- exp(-0.5 * (xs / sigma_px)^2); k1 <- k1 / sum(k1)
  K  <- outer(k1, k1); K / sum(K)
}
smooth_dem_and_derive <- function(E, alt, az, radius_m) {
  resx <- terra::res(E)[1]
  pad_cells <- ceiling(radius_m / resx) + 2L
  
  E_pad <- terra::extend(E, pad_cells)
  
  K <- gaussian_focal(E_pad, radius_m)
  Es_pad  <- terra::focal(E_pad, w = K, fun = mean, na.policy = "omit", pad = TRUE)
  
  slp_pad <- terra::terrain(Es_pad, v = "slope",  unit = "radians")
  asp_pad <- terra::terrain(Es_pad, v = "aspect", unit = "radians")
  
  ci_pad  <- cosi_fun(alt, az, slp_pad, asp_pad)
  
  list(
    Es   = terra::crop(Es_pad,  E),
    slp  = terra::crop(slp_pad, E),
    cosi = terra::crop(ci_pad,  E)
  )
}

 .extract_to_pts <- function(r, pts_sf) {
  out <- try(terra::extract(r, terra::vect(pts_sf), ID = FALSE)[,1], silent = TRUE)
  if (inherits(out, "try-error") || length(out) == 0L) rep(NA_real_, nrow(pts_sf)) else out
}
 cv_gam_with_R <- function(stn_sf, E, alt = NULL, az = NULL, R, block_size_m = NULL, verbose = TRUE,...) 
   {
   t0 <- proc.time()
   # robust check whether to compute cos(i)
   use_cosi <- isTRUE(!is.null(alt) && !is.null(az) &&
                        is.finite(alt) && is.finite(az)) 
   
   # --- 0) Block size guard
   bs <- suppressWarnings(as.numeric(block_size_m)[1])
   if (!is.finite(bs) || bs <= 0) {
     bs <- suppressWarnings(as.numeric(get0("block_size", ifnotfound = NA_real_)))
   }
   if (!is.finite(bs) || bs <= 0)
     stop("cv_gam_with_R(): no valid block size (block_size_m or global block_size).")
   
   # --- 1) DEM smoothing + derived fields
   zR   <- smooth_mean_R(E, R)
   slpR <- terra::terrain(zR, v = "slope",  unit = "radians")
   aspR <- terra::terrain(zR, v = "aspect", unit = "radians")
   
   # Only compute cos(i) if BOTH angles are clean, finite scalars
   use_cosi <- isTRUE(is.numeric(alt) && is.numeric(az) &&
                        length(alt) == 1L && length(az) == 1L &&
                        is.finite(alt) && is.finite(az))
   if (use_cosi) {
     zen  <- (pi/2 - alt)
     ci   <- cos(slpR)*cos(zen) + sin(slpR)*sin(zen)*cos(az - aspR)
     cosiR <- terra::ifel(ci < 0, 0, ci)
   } else {
     cosiR <- NULL
   }
   
   # --- 2) Extract to stations (fill missing with medians)
   if (!all(c("x","y") %in% names(stn_sf))) {
     xy <- sf::st_coordinates(stn_sf); stn_sf$x <- xy[,1]; stn_sf$y <- xy[,2]
   }
   fill_med <- function(v) {
     m <- stats::median(v[is.finite(v)], na.rm = TRUE)
     v[!is.finite(v)] <- m
     v
   }
   stn_sf$z_surf_R <- fill_med(.extract_to_pts(zR,   stn_sf))
   stn_sf$slp_R    <- fill_med(.extract_to_pts(slpR, stn_sf))
   if (use_cosi) {
     stn_sf$cosi_R <- fill_med(.extract_to_pts(cosiR, stn_sf))
   } else {
     stn_sf$cosi_R <- NA_real_
   }
   
   # --- 3) Build blocks
   bb_poly <- sf::st_as_sfc(sf::st_bbox(stn_sf), crs = sf::st_crs(stn_sf))
   blocks  <- sf::st_make_grid(bb_poly, cellsize = c(bs, bs), what = "polygons")
   blocks  <- sf::st_sf(block_id = seq_along(blocks), geometry = blocks)
   
   stn_blk <- sf::st_join(stn_sf, blocks, join = sf::st_intersects, left = TRUE)
   if (anyNA(stn_blk$block_id)) {
     i <- is.na(stn_blk$block_id)
     stn_blk$block_id[i] <- blocks$block_id[sf::st_nearest_feature(stn_blk[i,], blocks)]
   }
   if (!all(c("x","y") %in% names(stn_blk))) {
     xy <- sf::st_coordinates(stn_blk); stn_blk$x <- xy[,1]; stn_blk$y <- xy[,2]
   }
   
   # --- 4) LBO-CV
   bids  <- sort(unique(stn_blk$block_id))
   pm_say("[cv_gam_with_R] R=%.0f m | block=%.0f m | stations=%d | blocks=%d | cos(i)=%s",
          R, block_size_m, nrow(stn_sf), length(bids),
          if (use_cosi) "yes" else "no", v = verbose)
   preds <- vector("list", length(bids)); j <- 0L
   
   for (b in bids) {
     te <- stn_blk[stn_blk$block_id == b, ]
     tr <- stn_blk[stn_blk$block_id != b, ]
     pm_say("  - block %d: train=%d, test=%d", b, nrow(tr), nrow(te), v = verbose)
     
     
     dtr  <- sf::st_drop_geometry(tr)
     # include cosi_R only if it’s present with finite values
     need <- c("temp","x","y","z_surf_R","slp_R")
     inc_cosi <- ("cosi_R" %in% names(dtr)) && any(is.finite(dtr$cosi_R))
     if (inc_cosi) need <- c(need, "cosi_R")
     
     dtr  <- dtr[stats::complete.cases(dtr[, need, drop = FALSE]), need, drop = FALSE]
     if (nrow(dtr) < 10) next
     
     # dynamic k guards
     n_xy <- dplyr::n_distinct(paste0(round(dtr$x,3), "_", round(dtr$y,3)))
     k_xy <- .k_for_xy(nrow(dtr), n_xy)
     k_z  <- .kcap_unique(dtr$z_surf_R, 20)
     k_sl <- .kcap_unique(dtr$slp_R,    12)
     if (inc_cosi) k_ci <- .kcap_unique(dtr$cosi_R, 12)
     
     # assemble formula with only informative terms
     terms <- c()
     terms <- c(terms, if (n_xy >= 4) sprintf("s(x,y,bs='tp',k=%d)", k_xy) else "x + y")
     terms <- c(terms, if (k_z  >= 4) sprintf("s(z_surf_R,bs='tp',k=%d)", k_z)  else "z_surf_R")
     if (length(unique(dtr$slp_R[is.finite(dtr$slp_R)])) > 1)
       terms <- c(terms, if (k_sl >= 4) sprintf("s(slp_R,bs='tp',k=%d)", k_sl) else "slp_R")
     if (inc_cosi && any(is.finite(dtr$cosi_R)) &&
         length(unique(dtr$cosi_R[is.finite(dtr$cosi_R)])) > 1)
       terms <- c(terms, if (k_ci >= 4) sprintf("s(cosi_R,bs='tp',k=%d)", k_ci) else "cosi_R")
     
     form <- stats::as.formula(paste("temp ~", paste(terms, collapse = " + ")))
     gm <- mgcv::gam(form, data = dtr, method = "REML", select = TRUE)
     
     dte <- sf::st_drop_geometry(te)
     # restrict to variables actually in the model
     vars_needed <- setdiff(all.vars(form), "temp")
     dte <- dte[, vars_needed, drop = FALSE]
     ph  <- try(stats::predict(gm, newdata = dte, type = "response"), silent = TRUE)
     if (inherits(ph, "try-error")) ph <- rep(NA_real_, nrow(dte))
     
     j <- j + 1L
     preds[[j]] <- tibble::tibble(id = te$id, obs = te$temp, pred = as.numeric(ph), block_id = b)
   }
   
   preds <- preds[seq_len(j)]
   if (!length(preds)) {
     return(list(cv = tibble::tibble(id = integer(), obs = numeric(), pred = numeric(), block_id = integer()),
                 RMSE = NA_real_))
   }
   out  <- dplyr::bind_rows(preds)
   rmse <- sqrt(mean((out$pred - out$obs)^2, na.rm = TRUE))
   list(cv = out, RMSE = rmse)
 }
 
 
 

 tune_Rstar_ucurve <- function(stn_sf, E, alt = NULL, az = NULL,
                               L50, L95, block_fallback = 120,
                               n_grid = 6, extra = c(0.8, 1.2),
                               scen = NULL, which_time = c("T14","T05")) {
   
   which_time <- match.arg(which_time)
   
   # fallback L50/L95 if broken
   e <- terra::ext(E)
   dom_diag <- sqrt((terra::xmax(e)-terra::xmin(e))^2 + (terra::ymax(e)-terra::ymin(e))^2)
   if (!is.finite(L50) || !is.finite(L95) || L95 <= L50) {
     L50 <- dom_diag/10; L95 <- dom_diag/4
   }
   block_m <- max(block_fallback, round(L50))
   
   # sun from scen if not given
   if ((is.null(alt) || is.null(az)) && !is.null(scen)) {
     s <- .get_sun(scen, which_time)
     alt <- s$alt; az <- s$az
   }
   
   R_min <- max(10, round(L50*extra[1])); R_max <- round(L95*extra[2])
   R_grid <- unique(round(seq(R_min, R_max, length.out = n_grid)))
   
   rows <- lapply(R_grid, function(R) {
     z <- cv_gam_with_R(stn_sf, E, alt = alt, az = az, R = R,
                        block_size_m = block_m, scen = NULL, which_time = which_time)
     data.frame(R = R, RMSE = z$RMSE)
   })
   df <- do.call(rbind, rows)
   R_star <- df$R[which.min(df$RMSE)]
   list(grid = df, R_star = as.numeric(R_star), block_m = block_m)
 }
 
 

plot_ucurve <- function(df, R_star, title = "U-curve: tune R") {
  ggplot2::ggplot(df, ggplot2::aes(R, RMSE)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = R_star, linetype = "dashed", colour = "#de2d26") +
    ggplot2::theme_minimal() + ggplot2::labs(title = title, x = "Drift radius R (m)", y = "RMSE (block-CV)")
}

# Drop-in replacement
add_drifts_at_R <- function(stn_sf, E, alt, az, R,
                            lc = NULL, lc_levels = NULL,
                            na_action = c("error","fill","drop")) {
  na_action <- match.arg(na_action)
  
  # 0) Align CRS (key cause of NA extractions)
  crs_r <- sf::st_crs(E)
  if (!isTRUE(sf::st_crs(stn_sf) == crs_r)) {
    stn_sf <- sf::st_transform(stn_sf, crs_r)
  }
  
  # 1) Build @R* features (Es, slope, cosi) — your existing function
  fr <- smooth_dem_and_derive(E, alt, az, radius_m = R)
  
  # 2) Extract to points
  v <- terra::vect(stn_sf)
  stn_sf$E_R    <- as.numeric(terra::extract(fr$Es,   v, ID = FALSE)[, 1])
  stn_sf$slp_R  <- as.numeric(terra::extract(fr$slp,  v, ID = FALSE)[, 1])
  stn_sf$cosi_R <- as.numeric(terra::extract(fr$cosi, v, ID = FALSE)[, 1])
  
  # Optional LC (factor) — unchanged logic
  if (!is.null(lc)) {
    if (is.null(lc_levels)) lc_levels <- lc_levels_default
    lc_idx <- as.integer(terra::extract(lc, v, ID = FALSE)[, 1])
    lc_idx[!is.finite(lc_idx)] <- 1L
    lc_idx <- pmax(1L, pmin(lc_idx, length(lc_levels)))
    stn_sf$lc <- factor(lc_levels[lc_idx], levels = lc_levels)
  }
  
  # 3) Handle NA per policy
  d <- sf::st_drop_geometry(stn_sf)
  miss <- !stats::complete.cases(d[, c("E_R","slp_R","cosi_R"), drop = FALSE])
  
  if (any(miss)) {
    if (na_action == "error") {
      stop("Station features at R* contain NA. Increase padding in smooth_dem_and_derive(), ",
           "reduce R*, or call add_drifts_at_R(..., na_action='fill'/'drop').")
    }
    if (na_action == "fill") {
      fill_med <- function(x) { m <- stats::median(x[is.finite(x)], na.rm = TRUE); x[!is.finite(x)] <- m; x }
      stn_sf$E_R    <- fill_med(stn_sf$E_R)
      stn_sf$slp_R  <- fill_med(stn_sf$slp_R)
      stn_sf$cosi_R <- fill_med(stn_sf$cosi_R)
    }
    if (na_action == "drop") {
      stn_sf <- stn_sf[!miss, ]
    }
  }
  
  stn_sf
}

compute_Ls_from_points <- function(stn_sf, value_col = "temp",
                                   maxdist = NULL, nlag = 18, smooth_k = 3) {
  stopifnot(inherits(stn_sf, "sf"), value_col %in% names(stn_sf))
  pts <- stn_sf[is.finite(stn_sf[[value_col]]), ]
  if (is.null(maxdist)) {
    bb <- sf::st_bbox(pts)
    dom_diag <- sqrt((bb["xmax"]-bb["xmin"])^2 + (bb["ymax"]-bb["ymin"])^2)
    maxdist <- dom_diag / 2
  }
  form <- stats::as.formula(sprintf("%s ~ 1", value_col))
  vg  <- gstat::variogram(form, data = pts, cutoff = maxdist, width = maxdist/nlag)
  if (nrow(vg) >= smooth_k) {
    vg$gamma <- stats::filter(vg$gamma, rep(1/smooth_k, smooth_k), sides = 2)
    vg$gamma[!is.finite(vg$gamma)] <- zoo::na.approx(vg$gamma, na.rm = FALSE)
    vg$gamma <- zoo::na.locf(zoo::na.locf(vg$gamma, fromLast = TRUE))
  }
  sill <- max(vg$gamma, na.rm = TRUE)
  if (!is.finite(sill) || sill <= 0) sill <- stats::median(vg$gamma, na.rm = TRUE)
  L_at_q <- function(q) {
    thr <- q * sill
    i   <- which(vg$gamma >= thr)[1]
    if (is.na(i)) return(NA_real_)
    if (i == 1) return(vg$dist[1])
    d0 <- vg$dist[i-1]; d1 <- vg$dist[i]
    g0 <- vg$gamma[i-1]; g1 <- vg$gamma[i]
    if (!is.finite(d0) || !is.finite(d1) || g1 == g0) return(d1)
    d0 + (thr - g0) * (d1 - d0) / (g1 - g0)
  }
  list(vg = vg, sill = sill, L50 = L_at_q(0.5), L95 = L_at_q(0.95), cutoff = maxdist)
}

# ===== Alpha from CV residuals (nugget fraction) =====

# Build per-station coordinates from an sf station layer
build_station_xy <- function(st_sf, id_col = "id") {
  stopifnot(inherits(st_sf, "sf"), id_col %in% names(st_sf))
  xy <- sf::st_coordinates(st_sf)
  out <- st_sf |>
    sf::st_drop_geometry() |>
    dplyr::transmute(id = .data[[id_col]], X = xy[,1], Y = xy[,2])
  out
}

# Aggregate CV residuals per station for a chosen model
cv_residuals_by_station <- function(cv_tbl, model, id_col = "id") {
  stopifnot(all(c("model", id_col, "obs", "pred") %in% names(cv_tbl)))
  cv_tbl |>
    dplyr::filter(.data$model == !!model) |>
    dplyr::mutate(resid = .data$obs - .data$pred) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      id    = dplyr::first(.data[[id_col]]),
      resid = mean(.data$resid, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::select(id, resid)
}

# Empirical residual variogram from X/Y/resid table
residual_variogram_xy <- function(df_xy, cutoff = NULL, width = NULL, cressie = TRUE) {
  stopifnot(all(c("X", "Y", "resid") %in% names(df_xy)))
  df_xy <- dplyr::filter(df_xy, is.finite(.data$X), is.finite(.data$Y), is.finite(.data$resid))
  if (nrow(df_xy) < 5) return(NULL)
  if (is.null(cutoff)) {
    dmax <- max(stats::dist(as.matrix(df_xy[, c("X","Y")])))
    cutoff <- as.numeric(dmax) * 0.5
  }
  if (is.null(width)) width <- cutoff / 12
  gstat::variogram(resid ~ 1, data = df_xy, locations = ~ X + Y,
                   cutoff = cutoff, width = width, cressie = cressie)
}

# Robust nugget fraction (tries multiple models; falls back to non-parametric ratio)
alpha_from_vg_safe <- function(vg) {
  if (is.null(vg) || nrow(vg) < 3) return(NA_real_)
  sill0  <- max(vg$gamma, na.rm=TRUE)
  nug0   <- max(min(vg$gamma, na.rm=TRUE), 1e-6)
  range0 <- stats::quantile(vg$dist, 0.33, na.rm=TRUE)
  
  cands <- list(
    gstat::vgm(psill = sill0 - nug0, model="Exp", range=range0, nugget=nug0),
    gstat::vgm(psill = sill0 - nug0, model="Sph", range=range0, nugget=nug0),
    gstat::vgm(psill = sill0 - nug0, model="Gau", range=range0, nugget=nug0),
    gstat::vgm(psill = sill0 - nug0, model="Mat", range=range0, nugget=nug0, kappa=0.5)
  )
  fits <- lapply(cands, function(m) try(gstat::fit.variogram(vg, m), silent=TRUE))
  ok   <- vapply(fits, function(f) inherits(f, "variogramModel"), logical(1))
  if (any(ok)) {
    sse <- vapply(fits[ok], function(f) attr(f, "SSErr"), numeric(1))
    fit <- fits[ok][[which.min(sse)]]
    nug <- fit$psill[fit$model == "Nug"]; sill <- sum(fit$psill)
    if (length(nug) && is.finite(sill) && sill > 0) return(as.numeric(nug/sill))
  }
  # Non-parametric fallback: early/late semivariance ratio (clamped)
  g_first <- min(vg$gamma[vg$np > 0], na.rm=TRUE)
  g_last  <- max(vg$gamma[vg$np > 0], na.rm=TRUE)
  nf <- g_first / g_last
  max(0, min(1, nf))
}

# One-call wrapper: get alpha from CV table + stations (by model)
alpha_from_cv <- function(cv_tbl, st_sf, model, id_col_st = "id", id_col_cv = "id",
                          cutoff = NULL, width = NULL, cressie = TRUE) {
  stopifnot(inherits(st_sf, "sf"))
  # per-station residuals for that model
  res <- cv_residuals_by_station(cv_tbl, model = model, id_col = id_col_cv)
  # station coordinates
  st_xy <- build_station_xy(st_sf, id_col = id_col_st)
  dat   <- dplyr::left_join(res, st_xy, by = dplyr::join_by(id == id))
  vg    <- residual_variogram_xy(dat, cutoff = cutoff, width = width, cressie = cressie)
  alpha_from_vg_safe(vg)
}

# ===== Optional helpers for the budget =====

# Cap sigma_inst so instrument variance never exceeds CV variance
cap_sigma_inst <- function(sigma_inst, rmse, frac = 0.95) {
  pmin(sigma_inst, frac * rmse)
}

# Convenience RMSE from a CV table for a given model
rmse_from_cv <- function(cv_tbl, model) {
  d <- dplyr::filter(cv_tbl, .data$model == !!model)
  sqrt(mean((d$obs - d$pred)^2, na.rm = TRUE))
}



simple_error_budget <- function(res_cv, sigma_inst = 0.5, alpha = 0.6) {
  res <- res_cv$cv
  res <- res[is.finite(res$obs) & is.finite(res$pred), , drop = FALSE]
  RMSE <- sqrt(mean((res$pred - res$obs)^2))
  Bias <- mean(res$pred - res$obs)
  VarE <- stats::var(res$pred - res$obs, na.rm = TRUE)
  meas <- sigma_inst^2
  proc <- max(0, VarE - meas)
  micro <- alpha * proc
  meso  <- (1 - alpha) * proc
  tibble::tibble(Component = c("RMSE","Bias","Total var","Instrument var","Microscale var","Mesoscale var"),
                 Value     = c(RMSE, Bias, VarE, meas, micro, meso))
}

# --- Save a tuned maps bundle (rasters + plots + table) ---
save_maps_bundle <- function(maps, which_time,
                             save_rasters = TRUE,
                             save_plots   = TRUE,
                             save_df      = TRUE) {
  stopifnot(!missing(which_time))
  
  out_files <- character(0)
  
  # 1) SpatRaster predictions per model  → GeoTIFF
  if (save_rasters && !is.null(maps$pred_rasters)) {
    for (mdl in names(maps$pred_rasters)) {
      r <- maps$pred_rasters[[mdl]]
      # optional: trim NA borders if you prefer tight extents
      # r <- terra::trim(r)
      f <- fn_ras(sprintf("%s_pred_%s", which_time, mdl))
      save_raster_min(r, f)
      out_files <- c(out_files, f)
    }
  }
  
  # 2) Prediction/Truth panels  → PNG
  if (save_plots) {
    if (!is.null(maps$p_pred)) {
      f <- fn_fig(sprintf("%s_pred_panel", which_time))
      save_plot_min(maps$p_pred, f); out_files <- c(out_files, f)
    }
    if (!is.null(maps$p_truth)) {
      f <- fn_fig(sprintf("%s_truth_panel", which_time))
      save_plot_min(maps$p_truth, f); out_files <- c(out_files, f)
    }
  }
  
  # 3) pred_df (x,y,pred per model)  → CSV
  if (save_df && !is.null(maps$pred_df)) {
    # ensure tab_dir exists in your setup; if not, use fn for tables you already have
    f <- file.path(tab_dir, sprintf("%s_pred_df.csv", which_time))
    readr::write_csv(maps$pred_df, f)
    out_files <- c(out_files, f)
  }
  
  message(sprintf("✔ Saved %s bundle: %s", which_time,
                  paste(basename(out_files), collapse = ", ")))
  invisible(out_files)
}


## ======================================================================
## Ende der Bibliothek
## ======================================================================

## ---------- Mini-Beispiel (nicht Teil der Bibliothek) -----------------
## domain  <- make_domain()
## scen    <- build_scenario(domain, lake_mode="water", hill_mode="bump",
##                           random_hills = 100, micro_seed = 1)
## pts_sf  <- make_stations(domain, n_st = 60, station_mode = "random")
## stns    <- stations_from_scenario(scen, pts_sf)
## bs      <- compute_block_size(len_x = domain$xmax-domain$xmin,
##                               len_y = domain$ymax-domain$ymin, n_st = nrow(pts_sf))
## out14   <- run_for_time(stns$T14, scen$R14, "T14", scen, bs)
## out05   <- run_for_time(stns$T05, scen$R05, "T05", scen, bs)
## # Plot-Beispiele:
## # print(plot_landcover_terrain(scen, stations = stns$T14))
## # print(out14$res$blocks_plot); print(out14$res$diag_plot)
