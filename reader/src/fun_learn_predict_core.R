#' Geostatistical Learners & Map Predictor (Core Only)
#'
#' @title Learners and Raster Predictor (no helpers inside)
#' @description
#' A compact set of **model-specific predictors** used in your teaching/
#' pipeline code, plus a high-level `predict_maps()` convenience that
#' evaluates multiple learners on a full grid.  
#'
#' This file intentionally contains **no helpers**. It assumes that common
#' utilities and constants are sourced from your *helpers* module, including:
#' - `%||%` — null-coalescing helper
#' - `.default_vgm()` — conservative variogram fallback
#' - `.align_factor_to_model()` — align factor levels at predict time
#' - `safe_gam_formula()` — guarded GAM formula constructor
#' - `lc_levels_default` — global land-cover levels
#' - `temp_palette`, `stretch_q` — visualization defaults
#'
#' @details
#' **Contract expected by all learners**:
#' - `train_sf`, `test_sf` are `sf` objects with at least:
#'   - `temp` (numeric): the response variable to be learned
#'   - `x`, `y` (numeric): planar coordinates (will be derived from geometry
#'     if absent)
#'   - Drift/covariate columns depending on the learner (see each function)
#' - Each learner returns a numeric vector of predictions aligned with
#'   `nrow(test_sf)`.
#'
#' **Coordinate Reference System**: all learners assume that `x` and `y`
#' are in a **projected CRS** with meter-like units (e.g., UTM).
#'
#' **Error handling**:
#' - Learners are defensive; if inputs are insufficient (e.g., too few rows,
#'   missing drift columns), they return `NA_real_` predictions of the correct
#'   length instead of failing hard (except where a *hard requirement* is unmet
#'   such as missing KED drifts in training).
#'
#' @section Dependencies:
#' - **Packages**: `sf`, `gstat`, `mgcv`, `randomForest`, `terra`, `ggplot2`,
#'   `tibble`, `dplyr`, `stats`, `scales`
#' - **Helpers (sourced elsewhere)**: `%||%`, `.default_vgm`, `.align_factor_to_model`,
#'   `safe_gam_formula`, `lc_levels_default`, `temp_palette`, `stretch_q`
#'
#' @seealso
#' - Your helpers/utilities module for the functions noted above.
#' - `gstat::krige`, `gstat::idw`, `gstat::variogram`, `gstat::fit.variogram`
#' - `mgcv::gam`, `randomForest::randomForest`
#'
#' @keywords geostatistics interpolation kriging regression GAM randomForest
#' @family learners


#' Voronoi / Nearest-Station Predictor
#'
#' @description
#' Assigns each prediction point the observed value from the **nearest**
#' training station (a fast proxy for Voronoi interpolation).
#'
#' @param train_sf `sf` with at least `temp` and geometry.
#' @param test_sf  `sf` with geometry to predict for.
#'
#' @return Numeric vector `length(nrow(test_sf))` with nearest-neighbor temps.
#' @examples
#' # y_hat <- pred_Voronoi(train_sf, grid_sf)
pred_Voronoi <- function(train_sf, test_sf) {
  idx <- sf::st_nearest_feature(test_sf, train_sf)
  as.numeric(train_sf$temp)[idx]
}


#' Inverse Distance Weighting (IDW)
#'
#' @description
#' Classic **IDW** using `gstat::idw`, predicting from training points to
#' the test geometry.
#'
#' @param train_sf `sf` with `temp` and geometry.
#' @param test_sf  `sf` with geometry.
#' @param idp      Inverse distance power (default `2`).
#'
#' @return Numeric vector of predictions for `test_sf`.
#' @examples
#' # y_hat <- pred_IDW(train_sf, grid_sf, idp = 2)
pred_IDW <- function(train_sf, test_sf, idp = 2) {
  pr <- suppressWarnings(gstat::idw(temp ~ 1, locations = train_sf, newdata = test_sf, idp = idp))
  as.numeric(pr$var1.pred)
}


#' Ordinary Kriging (OK)
#'
#' @description
#' Univariate **OK** with an automatically fitted **exponential** variogram.
#' Falls back to `.default_vgm()` if fitting fails (e.g., too few points).
#'
#' @param train_sf `sf` with `temp` and geometry.
#' @param test_sf  `sf` with geometry.
#'
#' @return Numeric vector of kriged predictions.
#' @examples
#' # y_hat <- pred_OK(train_sf, grid_sf)
pred_OK <- function(train_sf, test_sf) {
  vg      <- suppressWarnings(gstat::variogram(temp ~ 1, data = train_sf))
  vgm_fit <- try(suppressWarnings(gstat::fit.variogram(vg, gstat::vgm("Exp"))), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) vgm_fit <- .default_vgm(train_sf$temp)
  pr <- suppressWarnings(gstat::krige(temp ~ 1, locations = train_sf, newdata = test_sf, model = vgm_fit))
  as.numeric(pr$var1.pred)
}


#' Kriging with External Drift (KED)
#'
#' @description
#' **KED** with additive drift terms. Requires drifts in *training*, fills
#' non-finite values in *test* by median of training. If `lc` is present in
#' both sets, it is included as a categorical drift with aligned levels.
#'
#' @details
#' **Required drift columns** in `train_sf`: `z_surf`, `slp`, `cosi`.  
#' If any are missing in training, this function errors (by design).
#'
#' @param train_sf `sf`; must contain `temp`, `z_surf`, `slp`, `cosi`, geometry,
#'   and optionally `lc`.
#' @param test_sf  `sf` with geometry and preferably the same drift columns
#'   (non-finite values are median-filled).
#' @param ...      Unused (placeholder for compatibility).
#'
#' @return Numeric vector of KED predictions, `length(nrow(test_sf))`.
#' @examples
#' # y_hat <- pred_KED(train_sf, grid_sf)
pred_KED <- function(train_sf, test_sf, ...) {
  need <- c("z_surf","slp","cosi")
  miss <- setdiff(need, names(train_sf))
  if (length(miss)) stop("pred_KED(): missing drifts in training: ", paste(miss, collapse = ", "))
  use_lc <- "lc" %in% names(train_sf) && "lc" %in% names(test_sf)
  tr <- train_sf; te <- test_sf
  if (use_lc) {
    tr$lc <- droplevels(factor(tr$lc))
    te$lc <- factor(as.character(te$lc), levels = levels(tr$lc))
    te$lc[is.na(te$lc)] <- levels(tr$lc)[1]
  }
  for (nm in need) {
    m <- stats::median(tr[[nm]][is.finite(tr[[nm]])], na.rm = TRUE)
    te[[nm]][!is.finite(te[[nm]])] <- m
  }
  keep_tr <- c("temp", need, if (use_lc) "lc")
  dtr <- sf::st_drop_geometry(tr)[, keep_tr, drop = FALSE]
  ok  <- stats::complete.cases(dtr); tr <- tr[ok, ]
  if (nrow(tr) < 5) return(rep(NA_real_, nrow(te)))
  form <- stats::as.formula(paste("temp ~", paste(c(need, if (use_lc) "lc"), collapse = " + ")))
  vg      <- suppressWarnings(gstat::variogram(form, data = tr))
  vgm_fit <- try(suppressWarnings(gstat::fit.variogram(vg, gstat::vgm("Exp"))), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) {
    ps <- stats::var(sf::st_drop_geometry(tr)$temp, na.rm = TRUE)
    vgm_fit <- gstat::vgm(psill = ps, model = "Exp", range = max(vg$dist, na.rm = TRUE)/3, nugget = 0.1*ps)
  }
  pr <- suppressWarnings(gstat::krige(form, locations = tr, newdata = te, model = vgm_fit))
  as.numeric(pr$var1.pred)
}


#' Random Forest Regressor (RF)
#'
#' @description
#' A **RandomForest** on spatial and drift features. If `lc` is absent, a
#' harmless single-level factor is injected (levels provided by
#' `lc_levels_default`). At prediction, factor levels are aligned using
#' `.align_factor_to_model()`.
#'
#' @param train_sf `sf` with `temp`, `x`, `y`, `z_surf`, `slp`, `cosi`,
#'   optionally `lc` (factor), and geometry.
#' @param test_sf  `sf` with the same predictors (geometry required).
#'
#' @return Numeric vector of RF predictions.
#' @examples
#' # y_hat <- pred_RF(train_sf, grid_sf)
pred_RF <- function(train_sf, test_sf) {
  dtr <- sf::st_drop_geometry(train_sf)
  if (!("lc" %in% names(dtr))) dtr$lc <- factor(lc_levels_default[1], levels = lc_levels_default)
  dtr$lc <- droplevels(factor(as.character(dtr$lc), levels = lc_levels_default))
  dtr <- stats::na.omit(dtr)
  if (nrow(dtr) < 5) return(rep(NA_real_, nrow(test_sf)))
  rf  <- randomForest::randomForest(temp ~ x + y + z_surf + slp + cosi + lc, data = dtr, na.action = na.omit)
  dte <- sf::st_drop_geometry(test_sf)
  if (!("lc" %in% names(dte))) dte$lc <- factor(lc_levels_default[1], levels = lc_levels_default)
  lev <- levels(dtr$lc)
  dte$lc <- .align_factor_to_model(dte$lc, lev)
  good <- stats::complete.cases(dte[, c("x","y","z_surf","slp","cosi","lc")])
  out  <- rep(NA_real_, nrow(dte)); if (any(good)) out[good] <- stats::predict(rf, dte[good, ])
  out
}


#' Generalized Additive Model (GAM)
#'
#' @description
#' A **GAM** (thin-plate splines) built with a protective formula from
#' `safe_gam_formula()` that caps basis sizes and includes `lc` only if
#' useful. Requires a minimal number of complete rows.
#'
#' @param train_sf `sf` with `temp`, `x`, `y`, `z_surf`, `slp`, `cosi`,
#'   optionally `lc` (factor).
#' @param test_sf  `sf` with matching predictors.
#'
#' @return Numeric vector of GAM predictions; `NA_real_` if the model could
#'   not be trained.
#' @examples
#' # y_hat <- pred_GAM(train_sf, grid_sf)
pred_GAM <- function(train_sf, test_sf) {
  dtr  <- sf::st_drop_geometry(train_sf)
  keep <- intersect(c("temp","x","y","z_surf","slp","cosi","lc"), names(dtr))
  dtr  <- dtr[stats::complete.cases(dtr[, keep, drop = FALSE]), keep, drop = FALSE]
  if (!nrow(dtr)) return(rep(NA_real_, nrow(test_sf)))
  if ("lc" %in% names(dtr)) dtr$lc <- droplevels(factor(dtr$lc))
  inc_lc <- "lc" %in% names(dtr) && nlevels(dtr$lc) >= 2
  if (nrow(dtr) < 10) return(rep(NA_real_, nrow(test_sf)))
  gm <- mgcv::gam(formula = safe_gam_formula(dtr, include_lc = inc_lc), data = dtr, method = "REML", select = TRUE)
  dte <- sf::st_drop_geometry(test_sf)
  vars <- c("x","y","z_surf","slp","cosi", if (inc_lc) "lc"); vars <- intersect(vars, names(dte))
  if (inc_lc) {
    lev <- levels(model.frame(gm)$lc)
    if (!("lc" %in% names(dte))) dte$lc <- lev[1]
    dte$lc <- .align_factor_to_model(dte$lc, lev)
  }
  good <- stats::complete.cases(dte[, vars, drop = FALSE])
  out  <- rep(NA_real_, nrow(dte)); if (any(good)) out[good] <- stats::predict(gm, dte[good, vars, drop = FALSE], type = "response")
  out
}


#' Predict on a Raster Grid with Multiple Learners + Pretty Plots
#'
#' @description
#' High-level utility that:
#' 1. Ensures station covariates exist (E, slope, cos(i), optional LC).
#' 2. Builds a **full-grid** data frame of covariates from rasters.
#' 3. Runs selected learners (`Voronoi`, `IDW`, `OK`, `KED`, `RF`, `GAM`).
#' 4. Returns both **prediction rasters** and **ggplot** panels.
#'
#' @param stn_sf `sf` training stations; must have `temp` and (if missing)
#'   this function will derive `x`, `y` and extract missing covariates from
#'   rasters.
#' @param truth_raster `SpatRaster` (single-layer) used only for common
#'   color scaling in plots (and optional “truth” visualization).
#' @param which_time Character; `"T14"` or `"T05"` (plot titles only).
#' @param scen A scenario list containing at least: `E`, `slp`, and either
#'   `I14` or `I05` (for cos(i)) and optionally `lc` + `lc_levels`.
#' @param models Character vector of learners to run.
#' @param lc_levels Optional character vector of LC levels (defaults to
#'   `scen$lc_levels`).
#' @param feature_rasters Optional list with named rasters `E`, `slp`, `cosi`
#'   to **override** the scenario’s baseline (e.g., when using tuned R*).
#'
#' @return A list with:
#' \describe{
#'   \item{pred_df}{Tidy `tibble` of predictions for all models & grid cells}
#'   \item{pred_rasters}{`list` of `SpatRaster` predictions, one per model}
#'   \item{p_pred}{`ggplot` facet showing all model maps}
#'   \item{p_truth}{`ggplot` of the truth raster (for reference)}
#' }
#'
#' @note
#' Requires helpers/constants: `%||%`, `temp_palette`, `stretch_q`, plus
#' land-cover level alignment utilities.
#'
#' @examples
#' # out <- predict_maps(stn_sf, scen$R14, which_time = "T14", scen = scen)
#' # print(out$p_truth); print(out$p_pred)
predict_maps <- function(stn_sf, truth_raster,
                         which_time = c("T14","T05"),
                         scen, models = c("Voronoi","IDW","OK","KED","RF","GAM"),
                         lc_levels = NULL,
                         feature_rasters = NULL) {
  which_time <- match.arg(which_time)
  lc_levels  <- lc_levels %||% scen$lc_levels
  E      <- feature_rasters$E   %||% scen$E
  slp_r  <- feature_rasters$slp %||% scen$slp
  cosi_r <- feature_rasters$cosi %||% if (which_time == "T14") scen$I14 else scen$I05
  has_lc <- ("lc" %in% names(scen)) && !is.null(scen$lc)
  lc_r   <- if (has_lc) scen$lc else NULL
  
  train_sf <- stn_sf
  if (!all(c("x","y") %in% names(train_sf))) {
    xy <- sf::st_coordinates(train_sf); train_sf$x <- xy[,1]; train_sf$y <- xy[,2]
  }
  if (!("z_surf" %in% names(train_sf)))
    train_sf$z_surf <- as.numeric(terra::extract(E,      sf::st_coordinates(train_sf))[,1])
  if (!("slp" %in% names(train_sf)))
    train_sf$slp    <- as.numeric(terra::extract(slp_r,  sf::st_coordinates(train_sf))[,1])
  if (!("cosi" %in% names(train_sf)))
    train_sf$cosi   <- as.numeric(terra::extract(cosi_r, sf::st_coordinates(train_sf))[,1])
  if (has_lc && !("lc" %in% names(train_sf))) {
    lc_codes <- as.integer(terra::extract(lc_r, sf::st_coordinates(train_sf))[,1])
    lc_codes[is.na(lc_codes)] <- 1L
    lc_codes <- pmax(1L, pmin(lc_codes, length(lc_levels)))
    train_sf$lc <- factor(lc_levels[lc_codes], levels = lc_levels)
  }
  
  xy <- as.data.frame(terra::xyFromCell(E, 1:terra::ncell(E))); names(xy) <- c("x","y")
  grid_df <- xy
  grid_df$z_surf <- as.numeric(terra::values(E))
  grid_df$slp    <- as.numeric(terra::values(slp_r))
  grid_df$cosi   <- as.numeric(terra::values(cosi_r))
  if (has_lc) {
    lc_codes <- as.integer(terra::values(lc_r))
    lc_codes[!is.finite(lc_codes)] <- 1L
    lc_codes <- pmax(1L, pmin(lc_codes, length(lc_levels)))
    grid_df$lc <- factor(lc_levels[lc_codes], levels = lc_levels)
  }
  grid_sf <- sf::st_as_sf(grid_df, coords = c("x","y"),
                          crs = sf::st_crs(train_sf), remove = FALSE)
  
  use_lc <- has_lc && ("lc" %in% names(train_sf)) && ("lc" %in% names(grid_sf))
  if (use_lc) {
    lev <- levels(droplevels(factor(train_sf$lc)))
    train_sf$lc <- factor(as.character(train_sf$lc), levels = lev)
    grid_sf$lc  <- factor(as.character(grid_sf$lc),  levels = lev)
    if (anyNA(train_sf$lc) || anyNA(grid_sf$lc)) {
      use_lc <- FALSE; train_sf$lc <- NULL; grid_sf$lc <- NULL
    }
  }
  
  pred_list <- list()
  if ("Voronoi" %in% models) pred_list$Voronoi <- pred_Voronoi(train_sf, grid_sf)
  if ("IDW"     %in% models) pred_list$IDW     <- pred_IDW    (train_sf, grid_sf, idp = 2)
  if ("OK"      %in% models) pred_list$OK      <- pred_OK     (train_sf, grid_sf)
  if ("KED"     %in% models) pred_list$KED     <- pred_KED    (train_sf, grid_sf)
  if ("RF"      %in% models) {
    dtr <- sf::st_drop_geometry(train_sf)
    rf_vars <- c("x","y","z_surf","slp","cosi", if (use_lc) "lc")
    dtr <- stats::na.omit(dtr[, c("temp", rf_vars), drop = FALSE])
    pred_list$RF <- if (nrow(dtr) >= 5) {
      rf <- randomForest::randomForest(stats::as.formula(paste("temp ~", paste(rf_vars, collapse = " + "))),
                                       data = dtr, na.action = na.omit)
      as.numeric(stats::predict(rf, sf::st_drop_geometry(grid_sf)[, rf_vars, drop = FALSE]))
    } else rep(NA_real_, nrow(grid_sf))
  }
  if ("GAM"     %in% models) pred_list$GAM     <- pred_GAM    (train_sf, grid_sf)
  
  pred_df <- dplyr::bind_rows(lapply(names(pred_list), function(nm) {
    tibble::tibble(model = nm, x = grid_df$x, y = grid_df$y, pred = pred_list[[nm]])
  }))
  
  make_r <- function(vals, template = E) { r <- terra::rast(template); terra::values(r) <- as.numeric(vals); r }
  pred_rasters <- lapply(pred_list, make_r)
  
  truth_df <- as.data.frame(truth_raster, xy = TRUE, na.rm = FALSE)
  names(truth_df) <- c("x","y","truth")
  lims <- stats::quantile(truth_df$truth, probs = stretch_q, na.rm = TRUE)
  
  p_pred <- ggplot2::ggplot(pred_df, ggplot2::aes(x, y, fill = pred)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colors = temp_palette(256), limits = lims,
                                  oob = scales::squish, name = "Temp") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Predictions by model", which_time),
                  x = "Easting", y = "Northing") +
    ggplot2::facet_wrap(~ model, ncol = 3)
  
  p_truth <- ggplot2::ggplot(truth_df, ggplot2::aes(x, y, fill = truth)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colors = temp_palette(256), limits = lims,
                                  oob = scales::squish, name = "Temp") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Truth raster", which_time),
                  x = "Easting", y = "Northing")
  
  list(pred_df = pred_df, pred_rasters = pred_rasters, p_pred = p_pred, p_truth = p_truth)
}
