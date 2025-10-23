# =========================================================
# PipeModel end-to-end (clean & robust, unified "temp")
# - Half-pipe with left-side hill & right-side pond/hollow
# - Two time slots: 14 UTC (T14), 05 UTC (T05)
# - Station sampling (random or transect)
# - Models: Voronoi (NN), IDW, OK, KED, RF, GAM
# - Leave-Block-Out CV, metrics, truth & prediction maps
# - Enhanced outputs: overlays, obs–pred scatter, block-wise errors,
#   residual distributions, pretty tables
# =========================================================

# --------------------------- 0) Packages & options ---------------------------
req_pkgs <- c(
  "terra","sf","sp","ggplot2","dplyr","tibble","tidyr","scales","patchwork",
  "suncalc","gstat","randomForest","mgcv","kableExtra","knitr"
)
inst <- rownames(installed.packages())
if (any(!req_pkgs %in% inst)) install.packages(setdiff(req_pkgs, inst), dependencies = TRUE)
invisible(lapply(req_pkgs, require, character.only = TRUE))

sf::sf_use_s2(FALSE)   # robust nearest-feature/joins on small projected domains
set.seed(42)

# --------------------------- 1) Global knobs ---------------------------
# Domain & grid
crs_utm <- "EPSG:32632"
E0 <- 600000; N0 <- 5725000
len_x <- 500; len_y <- 300; res <- 5

# Scenario controls
lake_choice <- "water"   # "none" | "water" | "hollow"
hill_choice <- "bump"    # "none" | "bump"
lake_diam_m  <- 80; lake_depth_m <- 10; smooth_edges <- FALSE
hill_diam_m  <- 80; hill_height_m <- 50; hill_smooth <- FALSE
pool_block_gain <- 0.4   # weaken pooling over the hill (0..1)

# Temperature palette (blue -> red)
temp_palette <- grDevices::colorRampPalette(c("#0000FF","#FF0000"))(256)
stretch_q    <- c(0.02, 0.98)

# Stations
station_mode      <- "random"     # "random" | "ns_transect" | "ew_transect"
n_st              <- 60
transect_margin_m <- 10
ns_offset_m <- 0     # + east / - west
ew_offset_m <- 0     # + north / - south

# LBO-CV
block_size <- 100  # meters
models_use <- c("Voronoi","IDW","OK","KED","RF","GAM")

# --------------------------- 2) Domain helpers ---------------------------
ext <- terra::ext(E0 - len_x/2, E0 + len_x/2, N0 - len_y/2, N0 + len_y/2)
Rtemplate <- terra::rast(ext, resolution = res, crs = crs_utm)
xmin <- terra::xmin(ext); xmax <- terra::xmax(ext)
ymin <- terra::ymin(ext); ymax <- terra::ymax(ext)
x0 <- (xmin+xmax)/2; y0 <- (ymin+ymax)/2
x_hill_center <- xmin + len_x/3      # left third
y_hill_center <- y0
x_lake_center <- xmin + 2*len_x/3    # right third
y_lake_center <- y0

# Sun geometry
lat <- 51.8; lon <- 10.6
sun_pos_utc <- function(y, m, d, h, lat, lon) {
  t  <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00", y, m, d, h), tz = "UTC")
  sp <- suncalc::getSunlightPosition(date = t, lat = lat, lon = lon)
  az_from_north <- (sp$azimuth + pi) %% (2*pi)
  list(alt = sp$altitude, az = az_from_north)
}
sun14 <- sun_pos_utc(2024, 6, 21, 14, lat, lon)
sun05 <- sun_pos_utc(2024, 6, 21,  5, lat, lon)
cosi_fun <- function(alt, az, slp_r, asp_r) {
  zen <- (pi/2 - alt)
  ci  <- cos(slp_r)*cos(zen) + sin(slp_r)*sin(zen)*cos(az - asp_r)
  terra::ifel(ci < 0, 0, ci)
}

# --------------------------- 3) Scenario builder ---------------------------
build_scenario <- function(lake_mode = c("none","water","hollow"),
                           hill_mode = c("none","bump")) {
  lake_mode <- match.arg(lake_mode); hill_mode <- match.arg(hill_mode)
  XY <- as.data.frame(terra::xyFromCell(Rtemplate, 1:terra::ncell(Rtemplate))); names(XY) <- c("x","y")
  dy <- XY$y - y0
  a  <- 100 / ((len_y/2)^2)      # ~100 m rim height
  elev <- 500 + a * dy^2
  
  # Lake / hollow (right third)
  rl <- sqrt((XY$x - x_lake_center)^2 + (XY$y - y_lake_center)^2); lr <- lake_diam_m/2
  if (lake_mode %in% c("water","hollow")) {
    w_l <- if (smooth_edges) pmax(0, 1 - (rl/lr)^2) else as.numeric(rl <= lr)
    elev <- elev - lake_depth_m * w_l
  } else w_l <- 0
  
  # Hill (left third)
  if (hill_mode == "bump") {
    rh <- sqrt((XY$x - x_hill_center)^2 + (XY$y - y_hill_center)^2); hr <- max(1e-6, hill_diam_m/2)
    w_h <- if (hill_smooth) exp(- (rh/hr)^2) else as.numeric(rh <= hr)
    elev <- elev + hill_height_m * w_h
  } else w_h <- 0
  
  E <- Rtemplate; terra::values(E) <- elev; names(E) <- "elev"
  lakeR <- Rtemplate; terra::values(lakeR) <- if (lake_mode=="water") as.numeric(w_l>0) else 0; names(lakeR) <- "lake"
  hillW <- Rtemplate; terra::values(hillW) <- if (hill_mode=="bump") w_h else 0; names(hillW) <- "hillW"
  
  # Derivatives & sun
  slp  <- terra::terrain(E, v="slope",  unit="radians")
  asp  <- terra::terrain(E, v="aspect", unit="radians")
  slp0 <- terra::ifel(is.na(slp), 0, slp); asp0 <- terra::ifel(is.na(asp), 0, asp)
  I14  <- cosi_fun(sun14$alt, sun14$az, slp0, asp0)
  I05  <- cosi_fun(sun05$alt, sun05$az, slp0, asp0)
  
  # Physics
  set.seed(1001); noise14 <- terra::setValues(terra::rast(E), rnorm(terra::ncell(E), 0, 0.3))
  set.seed(1002); noise05 <- terra::setValues(terra::rast(E), rnorm(terra::ncell(E), 0, 0.3))
  E_mean <- terra::global(E, "mean", na.rm = TRUE)[1,1]
  Y <- terra::init(E, "y"); dist2ax <- abs(Y - y0); w_pool <- 70
  pool_base <- 4.0 * exp(- (dist2ax / w_pool)^2)
  pool      <- pool_base * (1 - pool_block_gain * hillW)
  alpha_I_grass <- 5.0; alpha_I_water <- 1.5
  alpha_map <- alpha_I_grass * (1 - lakeR) + alpha_I_water * lakeR
  warm_bias_water_dawn <- if (lake_mode=="water") 1.5 else 0
  
  T0_14 <- 26.0; lapse_14 <- -0.0065
  R14 <- T0_14 + lapse_14 * (E - E_mean) + alpha_map * I14 + noise14; names(R14) <- "T14"
  T0_05 <- 8.5; inv_05 <- 0.003; eta_slope <- 0.6
  R05 <- T0_05 + inv_05 * (E - E_mean) + eta_slope * slp0 - pool + warm_bias_water_dawn * lakeR + noise05; names(R05) <- "T05"
  
  # Replace non-finite with median
  fix_nonfinite <- function(r) {
    v <- terra::values(r); m <- stats::median(v[is.finite(v)], na.rm = TRUE); v[!is.finite(v)] <- m; terra::values(r) <- v; r
  }
  list(E = fix_nonfinite(E), R14 = fix_nonfinite(R14), R05 = fix_nonfinite(R05),
       lake = lakeR, hillW = hillW, slp = slp0, asp = asp0, I14 = I14, I05 = I05)
}
scen <- build_scenario(lake_choice, hill_choice)

# --------------------------- 4) Stations & extraction ---------------------------
if (station_mode == "random") {
  pts <- tibble::tibble(
    id = 1:n_st,
    x  = runif(n_st, xmin + transect_margin_m, xmax - transect_margin_m),
    y  = runif(n_st, ymin + transect_margin_m, ymax - transect_margin_m)
  )
} else if (station_mode == "ns_transect") {
  x_const <- min(max(x0 + ns_offset_m, xmin + transect_margin_m), xmax - transect_margin_m)
  y_seq <- seq(ymin + transect_margin_m, ymax - transect_margin_m, length.out = n_st)
  pts <- tibble::tibble(id = 1:n_st, x = x_const, y = y_seq)
} else if (station_mode == "ew_transect") {
  y_const <- min(max(y0 + ew_offset_m, ymin + transect_margin_m), ymax - transect_margin_m)
  x_seq <- seq(xmin + transect_margin_m, xmax - transect_margin_m, length.out = n_st)
  pts <- tibble::tibble(id = 1:n_st, x = x_seq, y = y_const)
} else stop("Unknown station_mode")

pts_sf <- sf::st_as_sf(pts, coords = c("x","y"), crs = crs_utm, remove = FALSE)
vpts   <- terra::vect(pts_sf)

# Extract station covariates & targets
pts$z_surf <- as.numeric(terra::extract(scen$E,   vpts, ID = FALSE)[,1])
pts$slp    <- as.numeric(terra::extract(scen$slp, vpts, ID = FALSE)[,1])
pts$I14    <- as.numeric(terra::extract(scen$I14, vpts, ID = FALSE)[,1])
pts$I05    <- as.numeric(terra::extract(scen$I05, vpts, ID = FALSE)[,1])
pts$T14    <- as.numeric(terra::extract(scen$R14, vpts, ID = FALSE)[,1])
pts$T05    <- as.numeric(terra::extract(scen$R05, vpts, ID = FALSE)[,1])

# Keep complete rows per time slot and build sf with unified response 'temp'
pts14 <- pts[stats::complete.cases(pts[, c("x","y","z_surf","slp","I14","T14")]), ]
pts05 <- pts[stats::complete.cases(pts[, c("x","y","z_surf","slp","I05","T05")]), ]

stn_sf_14 <- pts14 |>
  dplyr::transmute(id, x, y,
                   z_surf = as.numeric(z_surf),
                   slp    = as.numeric(slp),
                   cosi   = as.numeric(I14),
                   temp   = as.numeric(T14)) |>
  sf::st_as_sf(coords = c("x","y"), crs = crs_utm, remove = FALSE)

stn_sf_05 <- pts05 |>
  dplyr::transmute(id, x, y,
                   z_surf = as.numeric(z_surf),
                   slp    = as.numeric(slp),
                   cosi   = as.numeric(I05),
                   temp   = as.numeric(T05)) |>
  sf::st_as_sf(coords = c("x","y"), crs = crs_utm, remove = FALSE)

# --------------------------- 5) Blocking for spatial CV ---------------------------
make_blocks_and_assign <- function(pts_sf, E, block_size = 100) {
  bb <- sf::st_as_sfc(
    sf::st_bbox(
      c(xmin = terra::xmin(E), ymin = terra::ymin(E),
        xmax = terra::xmax(E), ymax = terra::ymax(E)),
      crs = sf::st_crs(pts_sf)
    )
  )
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
  ggplot() +
    geom_sf(data = blocks, fill = NA, color = "grey50", linewidth = 0.2) +
    geom_sf(data = pts_blk, aes(color = factor(block_id)), size = 2) +
    scale_color_brewer(palette = "Set2", name = "Block") +
    coord_sf(expand = FALSE) + theme_minimal() +
    labs(title = title, x = "Easting", y = "Northing")
}

# --------------------------- 6) Interpolation algorithms ---------------------------
# Voronoi / Nearest neighbour
pred_Voronoi <- function(train_sf, test_sf) {
  idx <- sf::st_nearest_feature(test_sf, train_sf)
  as.numeric(train_sf$temp)[idx]
}
# IDW
pred_IDW <- function(train_sf, test_sf, idp = 2) {
  pr <- gstat::idw(temp ~ 1, as(train_sf["temp"], "Spatial"), newdata = as(test_sf, "Spatial"), idp = idp)
  as.numeric(pr$var1.pred)
}
# Ordinary Kriging
pred_OK <- function(train_sf, test_sf) {
  tr_sp <- as(train_sf["temp"], "Spatial")
  vg <- gstat::variogram(temp ~ 1, tr_sp)
  vgm_fit <- try(gstat::fit.variogram(vg, gstat::vgm("Exp")), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) vgm_fit <- gstat::vgm(variance = stats::var(train_sf$temp, na.rm=TRUE), model="Exp", range=100)
  kr <- gstat::krige(temp ~ 1, locations = tr_sp, newdata = as(test_sf, "Spatial"), model = vgm_fit)
  as.numeric(kr$var1.pred)
}
# KED / Universal Kriging with external drift (z_surf)
pred_KED <- function(train_sf, test_sf, E = NULL) {
  if (is.null(E)) stop("pred_KED: provide raster E for fallback elevation extraction.")
  add_z <- function(s) {
    if (!("z_surf" %in% names(s)) || any(!is.finite(s$z_surf))) {
      z <- terra::extract(E, sf::st_coordinates(s))[,1]; s$z_surf <- as.numeric(z)
    }; s
  }
  train_sf <- add_z(train_sf); test_sf <- add_z(test_sf)
  tr_sp <- sp::SpatialPointsDataFrame(
    coords = sf::st_coordinates(train_sf),
    data   = data.frame(temp = train_sf$temp, z_surf = as.numeric(train_sf$z_surf)),
    proj4string = sp::CRS(sf::st_crs(train_sf)$wkt)
  )
  tr_sp <- tr_sp[stats::complete.cases(tr_sp@data), ]
  te_sp <- sp::SpatialPointsDataFrame(
    coords = sf::st_coordinates(test_sf),
    data   = data.frame(z_surf = as.numeric(test_sf$z_surf)),
    proj4string = sp::CRS(sf::st_crs(test_sf)$wkt)
  )
  vg <- gstat::variogram(temp ~ z_surf, tr_sp)
  vgm_fit <- try(gstat::fit.variogram(vg, gstat::vgm("Exp")), silent = TRUE)
  if (inherits(vgm_fit, "try-error")) vgm_fit <- gstat::vgm(variance = stats::var(tr_sp$temp, na.rm=TRUE), model="Exp", range=100)
  pr <- gstat::krige(temp ~ z_surf, locations = tr_sp, newdata = te_sp, model = vgm_fit)
  as.numeric(pr$var1.pred)
}
# Random Forest
pred_RF <- function(train_sf, test_sf) {
  dtr <- sf::st_drop_geometry(train_sf); dtr <- stats::na.omit(dtr)
  if (nrow(dtr) < 5) return(rep(NA_real_, nrow(test_sf)))
  rf  <- randomForest::randomForest(temp ~ x + y + z_surf + slp + cosi, data = dtr, na.action = na.omit)
  dte <- sf::st_drop_geometry(test_sf)
  good <- stats::complete.cases(dte[, c("x","y","z_surf","slp","cosi")])
  out <- rep(NA_real_, nrow(dte)); if (any(good)) out[good] <- stats::predict(rf, dte[good, ])
  out
}
# GAM (safe)
fit_gam_safe <- function(stn_sf) {
  d <- stn_sf |> sf::st_drop_geometry()
  d <- d[stats::complete.cases(d), , drop = FALSE]
  n <- nrow(d); if (n < 10) stop("Too few stations for GAM: n=", n)
  n_xy <- dplyr::n_distinct(paste0(round(d$x,3), "_", round(d$y,3)))
  k_xy <- max(10, min(60, n_xy - 1, floor(n * 0.8)))
  k1   <- function(v, kmax=15) { ku <- length(unique(d[[v]])); max(4, min(ku - 1, kmax)) }
  mgcv::gam(temp ~ s(x,y,bs='tp',k=k_xy) + s(z_surf,bs='tp',k=k1("z_surf",20)) +
              s(slp,bs='tp',k=k1("slp",12)) + s(cosi,bs='tp',k=k1("cosi",12)),
            data = d, method = "REML", select = TRUE)
}
pred_GAM <- function(train_sf, test_sf) {
  dtr <- sf::st_drop_geometry(train_sf); dtr <- dtr[stats::complete.cases(dtr), , drop = FALSE]
  if (nrow(dtr) < 10) return(rep(NA_real_, nrow(test_sf)))
  gm <- fit_gam_safe(train_sf)
  dte <- sf::st_drop_geometry(test_sf)
  good <- stats::complete.cases(dte[, c("x","y","z_surf","slp","cosi")])
  out <- rep(NA_real_, nrow(dte)); if (any(good)) out[good] <- stats::predict(gm, dte[good, ], type="response")
  out
}

# --------------------------- 7) Leave-Block-Out CV ---------------------------
run_lbo_cv <- function(stn_sf, E, block_size = 100, models = models_use) {
  if (!all(c("x","y") %in% names(stn_sf))) {
    xy <- sf::st_coordinates(stn_sf); stn_sf$x <- xy[,1]; stn_sf$y <- xy[,2]
  }
  blk <- make_blocks_and_assign(stn_sf, E, block_size = block_size)
  blocks_sf <- blk$blocks; stn_blk <- blk$pts
  restore <- function(nm) if (!(nm %in% names(stn_blk))) stn_blk[[nm]] <<- stn_sf[[nm]][match(stn_blk$id, stn_sf$id)]
  for (nm in c("temp","z_surf","slp","cosi","x","y")) restore(nm)
  
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
                  "IDW"     = pred_IDW(train_sf, test_sf),
                  "OK"      = pred_OK(train_sf, test_sf),
                  "KED"     = pred_KED(train_sf, test_sf, E = E),
                  "RF"      = pred_RF(train_sf, test_sf),
                  "GAM"     = pred_GAM(train_sf, test_sf),
                  stop("Unknown model: ", m))
      tibble::tibble(model = m, id = test_sf$id, obs = test_sf$temp, pred = p, block_id = b)
    })
    out_list[[k]] <- dplyr::bind_rows(pred_tbl)
  }
  cv_tbl <- dplyr::bind_rows(out_list)
  metrics <- cv_tbl |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      n    = dplyr::n(),
      MAE  = mean(abs(obs - pred), na.rm = TRUE),
      RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
      Bias = mean(pred - obs, na.rm = TRUE),
      R2   = cor(obs, pred, use = "complete.obs")^2,
      .groups = "drop"
    ) |>
    dplyr::arrange(RMSE)
  
  diag_plot <- ggplot(cv_tbl, aes(obs, pred)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_point(alpha=0.7) +
    coord_equal() + theme_minimal() +
    labs(title = sprintf("LBO-CV (block = %dm) — Observed vs Predicted", block_size),
         x = "Observed", y = "Predicted") +
    facet_wrap(~ model)
  
  blocks_plot <- plot_blocks_grid(blocks_sf, stn_blk,
                                  title = sprintf("Blocks (%.0f m) & stations", block_size))
  list(cv = cv_tbl, metrics = metrics, diag_plot = diag_plot, blocks_plot = blocks_plot)
}

# --------------------------- 8) Prediction maps (grid) ---------------------------
predict_maps <- function(stn_sf, truth_raster, which_time = c("T14","T05"),
                         scen, models = c("Voronoi","IDW","OK","KED","RF","GAM")) {
  which_time <- match.arg(which_time)
  E      <- scen$E
  slp_r  <- scen$slp
  cosi_r <- if (which_time == "T14") scen$I14 else scen$I05
  
  # ensure training has all features & coords
  train_sf <- stn_sf
  if (!all(c("x","y") %in% names(train_sf))) {
    xy <- sf::st_coordinates(train_sf); train_sf$x <- xy[,1]; train_sf$y <- xy[,2]
  }
  if (!("z_surf" %in% names(train_sf))) train_sf$z_surf <- as.numeric(terra::extract(E,      sf::st_coordinates(train_sf))[,1])
  if (!("slp"    %in% names(train_sf))) train_sf$slp    <- as.numeric(terra::extract(slp_r,  sf::st_coordinates(train_sf))[,1])
  if (!("cosi"   %in% names(train_sf))) train_sf$cosi   <- as.numeric(terra::extract(cosi_r, sf::st_coordinates(train_sf))[,1])
  
  # prediction grid with z_surf
  grid_df <- as.data.frame(c(E, slp_r, cosi_r), xy = TRUE, na.rm = FALSE)
  names(grid_df) <- c("x","y","elev","slp","cosi")
  grid_df$z_surf <- grid_df$elev
  grid_sf <- sf::st_as_sf(grid_df, coords = c("x","y"),
                          crs = sf::st_crs(train_sf), remove = FALSE)
  good_grid <- stats::complete.cases(grid_df[, c("x","y","z_surf","slp","cosi")])
  grid_sp <- as(grid_sf, "Spatial")
  
  pred_list <- list()
  # Voronoi
  if ("Voronoi" %in% models) {
    idx <- sf::st_nearest_feature(grid_sf, train_sf)
    tmp <- rep(NA_real_, nrow(grid_df))
    tmp[good_grid] <- as.numeric(train_sf$temp)[idx[good_grid]]
    pred_list$Voronoi <- tmp
  }
  # IDW
  if ("IDW" %in% models) {
    tr_sp <- as(train_sf["temp"], "Spatial")
    pr    <- gstat::idw(temp ~ 1, locations = tr_sp, newdata = grid_sp, idp = 2)
    tmp   <- rep(NA_real_, nrow(grid_df)); tmp[good_grid] <- as.numeric(pr$var1.pred)[good_grid]
    pred_list$IDW <- tmp
  }
  # OK
  if ("OK" %in% models) {
    tr_sp <- as(train_sf["temp"], "Spatial")
    vg    <- gstat::variogram(temp ~ 1, tr_sp)
    vgm_fit <- try(gstat::fit.variogram(vg, gstat::vgm("Exp")), silent = TRUE)
    if (inherits(vgm_fit, "try-error"))
      vgm_fit <- gstat::vgm(variance = stats::var(train_sf$temp, na.rm = TRUE), model = "Exp", range = 100)
    pr  <- gstat::krige(temp ~ 1, locations = tr_sp, newdata = grid_sp, model = vgm_fit)
    tmp <- rep(NA_real_, nrow(grid_df)); tmp[good_grid] <- as.numeric(pr$var1.pred)[good_grid]
    pred_list$OK <- tmp
  }
  # KED
  if ("KED" %in% models) {
    tr_sp <- sp::SpatialPointsDataFrame(
      coords = sf::st_coordinates(train_sf),
      data   = data.frame(temp = train_sf$temp, z_surf = as.numeric(train_sf$z_surf)),
      proj4string = sp::CRS(sf::st_crs(train_sf)$wkt)
    )
    tr_sp <- tr_sp[stats::complete.cases(tr_sp@data), ]
    stopifnot("z_surf" %in% names(grid_sp@data))
    vg    <- gstat::variogram(temp ~ z_surf, tr_sp)
    vgm_fit <- try(gstat::fit.variogram(vg, gstat::vgm("Exp")), silent = TRUE)
    if (inherits(vgm_fit, "try-error"))
      vgm_fit <- gstat::vgm(variance = stats::var(tr_sp$temp, na.rm = TRUE), model = "Exp", range = 100)
    pr  <- gstat::krige(temp ~ z_surf, locations = tr_sp, newdata = grid_sp, model = vgm_fit)
    tmp <- rep(NA_real_, nrow(grid_df)); tmp[good_grid] <- as.numeric(pr$var1.pred)[good_grid]
    pred_list$KED <- tmp
  }
  # RF
  if ("RF" %in% models) {
    dtr <- sf::st_drop_geometry(train_sf); dtr <- stats::na.omit(dtr)
    rf  <- randomForest::randomForest(temp ~ x + y + z_surf + slp + cosi, data = dtr, na.action = na.omit)
    tmp <- rep(NA_real_, nrow(grid_df)); tmp[good_grid] <- stats::predict(rf, grid_df[good_grid, c("x","y","z_surf","slp","cosi")])
    pred_list$RF <- tmp
  }
  # GAM
  if ("GAM" %in% models) {
    dtr <- sf::st_drop_geometry(train_sf)
    n <- nrow(dtr); k_xy <- max(10, min(60, n - 5))
    gm  <- mgcv::gam(temp ~ s(x,y,k=k_xy) + s(z_surf,k=15) + s(slp,k=12) + s(cosi,k=12),
                     data = dtr, method = "REML")
    tmp <- rep(NA_real_, nrow(grid_df))
    tmp[good_grid] <- stats::predict(gm, grid_df[good_grid, c("x","y","z_surf","slp","cosi")], type="response")
    pred_list$GAM <- tmp
  }
  
  pred_df <- dplyr::bind_rows(lapply(names(pred_list), function(nm) {
    tibble::tibble(model = nm, x = grid_df$x, y = grid_df$y, pred = pred_list[[nm]])
  }))
  truth_df <- as.data.frame(truth_raster, xy = TRUE, na.rm = FALSE); names(truth_df) <- c("x","y","truth")
  lims <- stats::quantile(truth_df$truth, probs = c(0.02, 0.98), na.rm = TRUE)
  
  # (plain) maps: Predictions & Truth (used by panel builder below)
  p_pred <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = pred_df, ggplot2::aes(x, y, fill = pred)) +
    ggplot2::scale_fill_gradientn(colors = temp_palette, limits = lims,
                                  oob = scales::squish, name = "Temp") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Predictions by model", which_time),
                  x = "Easting", y = "Northing") +
    ggplot2::facet_wrap(~ model, ncol = 3)
  
  p_truth <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = truth_df, ggplot2::aes(x, y, fill = truth)) +
    ggplot2::scale_fill_gradientn(colors = temp_palette, limits = lims,
                                  oob = scales::squish, name = "Temp") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Truth raster", which_time),
                  x = "Easting", y = "Northing")
  
  list(pred_df = pred_df, p_pred = p_pred, p_truth = p_truth)
}

# --------------------------- 9) Visual helpers & tables ---------------------------
# Facet labels with RMSE/MAE
.make_labeller <- function(cv_tbl) {
  m <- cv_tbl |>
    dplyr::group_by(model) |>
    dplyr::summarise(RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
                     MAE  = mean(abs(obs - pred),    na.rm = TRUE),
                     .groups = "drop")
  setNames(sprintf("%s  (RMSE=%.2f · MAE=%.2f)", m$model, m$RMSE, m$MAE), m$model)
}

# Obs vs Pred scatter (per model)
make_obs_pred_scatter <- function(cv_tbl, which_time = "T14") {
  lab <- .make_labeller(cv_tbl)
  ggplot(cv_tbl, aes(obs, pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.7, shape = 16) +
    coord_equal() + theme_minimal() +
    labs(title = sprintf("%s — Observed vs Predicted (LBO-CV)", which_time),
         x = "Observed", y = "Predicted") +
    facet_wrap(~ model, ncol = 3, labeller = ggplot2::as_labeller(lab))
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
order_models_by_median_rmse <- function(cv_tbl) {
  bm <- block_metrics_long(cv_tbl)
  bm |>
    dplyr::filter(Metric == "RMSE") |>
    dplyr::group_by(model) |>
    dplyr::summarise(med = stats::median(Value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(med) |>
    dplyr::pull(model)
}
make_block_metric_box <- function(cv_tbl, which_time = "T14") {
  bm <- block_metrics_long(cv_tbl)
  lev <- order_models_by_median_rmse(cv_tbl)
  bm$model <- factor(bm$model, levels = lev)
  ggplot2::ggplot(bm, ggplot2::aes(model, Value)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.35, width = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                          fill = "white", colour = "black", stroke = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Block-wise errors (LBO-CV)", which_time),
                  subtitle = "Box = IQR · line = median · ◆ = mean",
                  x = "Model", y = "Error") +
    ggplot2::facet_wrap(~ Metric, scales = "free_y")
}
make_block_metric_scatter <- function(cv_tbl, which_time = "T14") {
  bm <- block_metrics_long(cv_tbl)
  order_rmse <- bm |>
    dplyr::filter(Metric == "RMSE") |>
    dplyr::group_by(model) |>
    dplyr::summarise(mu = mean(Value, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(mu) |>
    dplyr::pull(model)
  bm$model <- factor(bm$model, levels = order_rmse)
  ggplot(bm, aes(model, Value)) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 1.8) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "white", colour = "black", stroke = 0.5) +
    theme_minimal() +
    labs(title = sprintf("%s — Block-wise errors (LBO-CV)", which_time),
         subtitle = "Points = blocks · white diamond = mean per model",
         x = "Model", y = "Error") +
    facet_wrap(~ Metric, scales = "free_y")
}
make_abs_error_box <- function(cv_tbl, which_time = "T14") {
  df <- cv_tbl |> dplyr::mutate(abs_err = abs(pred - obs))
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
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("%s — Absolute errors per station (LBO-CV)", which_time),
                  subtitle = "Box = IQR · line = median · ◆ = mean",
                  x = "Model", y = "|pred − obs|")
}
make_residual_density <- function(cv_tbl, which_time = "T14") {
  cv_tbl2 <- cv_tbl |> dplyr::mutate(res = pred - obs)
  ggplot(cv_tbl2, aes(res)) +
    geom_histogram(bins = 25, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = sprintf("%s — Residual distribution per model (LBO-CV)", which_time),
         x = "Residual (pred − obs)", y = "Count") +
    facet_wrap(~ model, scales = "free_y")
}

# Pretty tables
pretty_table <- function(df, caption) {
  df |>
    knitr::kable(digits = 3, caption = caption, format = "html") |>
    kableExtra::kable_styling(full_width = FALSE, bootstrap_options = c("striped","hover","condensed")) |>
    kableExtra::row_spec(0, bold = TRUE)
}

# --------------------------- 10) Panels with overlays ---------------------------
# Prepare overlay contours (lake ring, hill disk)
lake_df <- as.data.frame(scen$lake,  xy = TRUE); names(lake_df) <- c("x","y","lake")
hill_df <- as.data.frame(scen$hillW, xy = TRUE); names(hill_df) <- c("x","y","hillW")
overlay_layers <- list(
  geom_contour(data = lake_df, aes(x, y, z = lake),
               breaks = 0.5, colour = "black", linewidth = 0.35),
  geom_contour(data = hill_df, aes(x, y, z = hillW),
               breaks = 0.5, colour = "black", linetype = "22", linewidth = 0.3)
)
theme_small_map <- theme(
  axis.title = element_text(size = 9),
  axis.text  = element_text(size = 8),
  strip.text = element_text(size = 9, face = "bold"),
  legend.title = element_text(size = 9),
  legend.text  = element_text(size = 8)
)

build_panels_with_errors <- function(maps, truth_raster, cv_tbl, stn_sf, which_time,
                                     temp_palette = temp_palette, stretch_q = stretch_q) {
  truth_df <- as.data.frame(truth_raster, xy = TRUE, na.rm = FALSE); names(truth_df) <- c("x","y","truth")
  pred_df  <- maps$pred_df
  err_df   <- dplyr::inner_join(pred_df, truth_df, by = c("x","y")) |>
    dplyr::mutate(err = pred - truth)
  emax <- max(abs(err_df$err), na.rm = TRUE); emax <- if (is.finite(emax)) emax else 1
  lims_T <- stats::quantile(truth_df$truth, probs = stretch_q, na.rm = TRUE)
  coords  <- sf::st_coordinates(stn_sf)
  st_xy   <- tibble::tibble(id = stn_sf$id, x = coords[,1], y = coords[,2])
  st_res  <- cv_tbl |> dplyr::mutate(resid = pred - obs) |> dplyr::select(model, id, resid) |> dplyr::left_join(st_xy, by = "id")
  lab <- .make_labeller(cv_tbl)
  
  p_truth <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = truth_df, ggplot2::aes(x, y, fill = truth)) +
    ggplot2::scale_fill_gradientn(colors = temp_palette, limits = lims_T, oob = scales::squish, name = "Temp") +
    overlay_layers + ggplot2::coord_equal() + ggplot2::theme_minimal() + theme_small_map +
    ggplot2::labs(title = sprintf("%s — Truth raster", which_time), x = "Easting", y = "Northing")
  
  p_pred <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = pred_df, ggplot2::aes(x, y, fill = pred)) +
    ggplot2::scale_fill_gradientn(colors = temp_palette, limits = lims_T, oob = scales::squish, name = "Temp") +
    overlay_layers + ggplot2::coord_equal() + ggplot2::theme_minimal() + theme_small_map +
    ggplot2::labs(title = sprintf("%s — Predictions by model", which_time), x = "Easting", y = "Northing") +
    ggplot2::facet_wrap(~ model, ncol = 3, labeller = ggplot2::as_labeller(lab))
  
  p_err <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = err_df, ggplot2::aes(x, y, fill = err)) +
    ggplot2::scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#de2d26",
                                  midpoint = 0, limits = c(-emax, emax), name = "Error") +
    overlay_layers +
    ggplot2::geom_point(data = st_res, ggplot2::aes(x, y, fill = resid),
                        shape = 21, colour = "black", size = 2, stroke = 0.25) +
    ggplot2::scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#de2d26",
                                  midpoint = 0, limits = c(-emax, emax), name = "Error") +
    ggplot2::coord_equal() + ggplot2::theme_minimal() + theme_small_map +
    ggplot2::labs(title = sprintf("%s — Error (pred − truth) with CV residuals", which_time),
                  x = "Easting", y = "Northing") +
    ggplot2::facet_wrap(~ model, ncol = 3, labeller = ggplot2::as_labeller(lab))
  
  (p_truth | p_pred) / p_err
}

# --------------------------- 11) RUN: CV + maps + visuals ---------------------------
message("Running LBO-CV and building maps for T14 ...")
res14  <- run_lbo_cv(stn_sf_14, scen$E, block_size = block_size, models = models_use)
maps14 <- predict_maps(stn_sf_14, scen$R14, which_time = "T14", scen = scen, models = models_use)

message("Running LBO-CV and building maps for T05 ...")
res05  <- run_lbo_cv(stn_sf_05, scen$E, block_size = block_size, models = models_use)
maps05 <- predict_maps(stn_sf_05, scen$R05, which_time = "T05", scen = scen, models = models_use)

# Panels (truth | predictions / errors) with overlays and nicer strips
panel14 <- build_panels_with_errors(maps14, scen$R14, res14$cv, stn_sf_14, "T14",
                                    temp_palette = temp_palette, stretch_q = stretch_q)
panel05 <- build_panels_with_errors(maps05, scen$R05, res05$cv, stn_sf_05, "T05",
                                    temp_palette = temp_palette, stretch_q = stretch_q)
panel14; panel05
# ggsave("panel_T14.png", panel14, width = 12, height = 8, dpi = 300)
# ggsave("panel_T05.png", panel05, width = 12, height = 8, dpi = 300)

# Scatterplots (Obs vs Pred)
p_scatter14 <- make_obs_pred_scatter(res14$cv, "T14")
p_scatter05 <- make_obs_pred_scatter(res05$cv, "T05")
p_scatter14; p_scatter05

# Block-wise error summaries (box & scatter) + absolute error box
p_block_box14 <- make_block_metric_box(res14$cv, "T14")
p_block_box05 <- make_block_metric_box(res05$cv, "T05")
p_block_sc14  <- make_block_metric_box(res14$cv, "T14")
p_block_sc05  <- make_block_metric_box(res05$cv, "T05")
p_abserr_box14 <- make_abs_error_box(res14$cv, "T14")
p_abserr_box05 <- make_abs_error_box(res05$cv, "T05")
p_block_box14; p_block_sc14; p_abserr_box14
p_block_box05; p_block_sc05; p_abserr_box05

# Residual distributions
make_residual_density(res14$cv, "T14")
make_residual_density(res05$cv, "T05")

# Blocks, diagnostics
print(res14$blocks_plot); print(res14$diag_plot)
print(res05$blocks_plot); print(res05$diag_plot)

# Pretty metrics tables
pretty_table(res14$metrics, "LBO-CV metrics — T14")
pretty_table(res05$metrics, "LBO-CV metrics — T05")

# Station table (sanity check)
pts |>
  dplyr::transmute(
    id,
    easting  = round(pts$x),
    northing = round(pts$y),
    z_surf   = round(z_surf,1),
    T14_C    = round(T14,1),
    T05_C    = round(T05,1)
  ) |>
  pretty_table("Station sample (sanity check)")
