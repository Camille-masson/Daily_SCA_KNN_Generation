download_hydroweb_snow <- function(
    AREA,
    START_DATE,
    END_DATE,
    COLLECTION_ID,
    API_KEY,
    raw_dir,
    downloads_dir = dirname(raw_dir),  # where to store the temporary ZIP
    ONLY_FSC = FALSE,
    verbose = TRUE
) {
  
  ### Description:
  # download_hydroweb_snow()
  # This function downloads Theia Sentinel-2 Snow Cover (L2B) tiles from hydroweb.next
  # for a given bbox (AREA) and time range (START_DATE–END_DATE) using the official
  # Python client py_hydroweb via reticulate. Results are extracted and copied into
  # raw_dir, while the ZIP archive is removed afterwards. A marker file (dates + bbox)
  # prevents re-running an identical request, and existing tiles are not duplicated.
  # If ONLY_FSC = TRUE, only FSC rasters are kept. Set verbose = TRUE for console logs.
  
  ## Check required packages ----
  stopifnot(requireNamespace("reticulate", quietly = TRUE))
  stopifnot(requireNamespace("fs", quietly = TRUE))
  stopifnot(requireNamespace("stringr", quietly = TRUE))
  library(reticulate)
  library(fs)
  library(stringr)
  
  ## Validate input parameters ----
  if (API_KEY == "" || is.null(API_KEY)) stop("Missing API_KEY.")
  if (length(AREA) != 4) stop("AREA must be c(xmin, ymin, xmax, ymax).")
  
  # Helpers
  to_iso <- function(x) format(x, "%Y-%m-%dT%H:%M:%SZ")
  stamp  <- function(x) gsub("[^0-9]", "", to_iso(x))
  
  # Dates
  start_dt <- as.POSIXct(START_DATE, tz = "UTC")
  end_dt   <- as.POSIXct(END_DATE, tz = "UTC")
  if (is.na(start_dt) || is.na(end_dt)) {
    stop("Invalid START_DATE / END_DATE. Use ISO UTC, e.g. '2023-09-01T00:00:00Z'.")
  }
  if (start_dt > end_dt) stop("START_DATE must be <= END_DATE.")
  
  ## Load py_hydroweb (Python) ----
  if (!py_module_available("py_hydroweb")) {
    if (verbose) message("py_hydroweb not found -> installing in the reticulate Python environment...")
    py_install("py-hydroweb", pip = TRUE)
  }
  py_require("py_hydroweb")
  py_hydroweb <- import("py_hydroweb", delay_load = FALSE)
  
  ## Create Hydroweb client ----
  client <- py_hydroweb$Client(
    "https://hydroweb.next.theia-land.fr/api",
    api_key = API_KEY
  )
  
  ## Skip duplicate requests (marker file) ----
  area_tag <- paste0(round(AREA, 5), collapse = "_")  # short signature of bbox
  marker <- file.path(
    raw_dir,
    paste0("done_", stamp(start_dt), "_", stamp(end_dt), "_", area_tag, ".txt")
  )
  
  if (file.exists(marker)) {
    if (verbose) message("OK: this exact request (dates + extent) was already downloaded -> skip.")
    return(invisible(list(status = "skipped", raw_dir = raw_dir, marker = marker)))
  }
  
  if (verbose) {
    message("==== Download ", COLLECTION_ID, " (",
            to_iso(start_dt), " -> ", to_iso(end_dt),
            ", bbox: ", paste(AREA, collapse = ","), ") ====")
  }
  
  ## Build download request (basket + filters) ----
  basket <- py_hydroweb$DownloadBasket("my_download_basket")
  q <- dict(
    start_datetime = dict(lte = to_iso(end_dt)),
    end_datetime   = dict(gte = to_iso(start_dt))
  )
  
  basket$add_collection(
    COLLECTION_ID,
    bbox  = as.list(AREA),
    query = q
  )
  
  ## Submit basket and download ZIP ----
  zip_name <- paste0(
    "hydroweb_", COLLECTION_ID, "_",
    stamp(start_dt), "_", stamp(end_dt), ".zip"
  )
  
  if (verbose) message("Downloading ZIP…")
  zip_path <- client$submit_and_download_zip(
    basket,
    zip_filename  = zip_name,
    output_folder = downloads_dir
  )
  
  ## Extract ZIP to a temporary folder ----
  tmp_unzip <- tempfile(pattern = "unz_")
  dir_create(tmp_unzip)
  unzip(zip_path, exdir = tmp_unzip)
  
  files <- dir_ls(tmp_unzip, recurse = TRUE, type = "file")
  
  # Optional FSC filter
  if (ONLY_FSC && length(files) > 0) {
    keep <- str_detect(basename(files), "FSC") & !str_detect(basename(files), "QC")
    files <- files[keep]
  }
  
  ## Copy files (skip duplicates) ----
  n_new  <- 0
  n_skip <- 0
  
  if (length(files) == 0) {
    if (verbose) message("No files found in the ZIP for this time range.")
  } else {
    for (f in files) {
      dest <- file.path(raw_dir, path_file(f))
      
      if (file_exists(dest)) {
        n_skip <- n_skip + 1
      } else {
        file_copy(f, dest)
        n_new <- n_new + 1
      }
    }
  }
  
  ## Cleanup temporary files ----
  dir_delete(tmp_unzip)
  file_delete(zip_path)
  
  ## Write marker file ----
  writeLines(
    c(
      paste0("collection=", COLLECTION_ID),
      paste0("bbox=", paste(AREA, collapse=",")),
      paste0("start=", to_iso(start_dt)),
      paste0("end=", to_iso(end_dt)),
      paste0("new_files=", n_new),
      paste0("skipped_existing=", n_skip),
      paste0("done_at_utc=", to_iso(Sys.time()))
    ),
    marker
  )
  
  if (verbose) {
    message("Done: +", n_new, " new files, ", n_skip, " already present.")
    message("Raw data in: ", normalizePath(raw_dir))
  }
  
  invisible(list(
    status = "downloaded",
    raw_dir = raw_dir,
    marker = marker,
    new_files = n_new,
    skipped_existing = n_skip
  ))
}






mosaic_crop_save_snow_binary_daily <- function(
    raw_dir,
    AOI_file,
    YEAR,
    output_dir,
    STUDY_AREA,
    cloud_thr = 5,     # % max de cloud (205) accepté
    snow_thr  = 70,    # seuil FSC -> neige
    verbose   = TRUE
) {
  
  
  START_DATE <- sprintf("%d-01-01T00:00:00", YEAR)
  END_DATE   <- sprintf("%d-12-31T23:59:59", YEAR)
  
  start_dt <- as.POSIXct(START_DATE, tz = "UTC")
  end_dt   <- as.POSIXct(END_DATE,   tz = "UTC")
  
  aoi <- terra::vect(AOI_file)
  
  # ---- OUTPUT folders (dans output_dir / 1. Snow_cover_by_site / STUDY_AREA / YYYY) ----
  safe_name <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)
  site_tag <- safe_name(STUDY_AREA)
  
  
  # ---- List FSC files ----
  f_all <- fs::dir_ls(raw_dir, recurse = FALSE, type = "file")
  f_fsc <- f_all[
    str_detect(basename(f_all), "FSC") &
      !str_detect(basename(f_all), "QCFLAGS") &
      str_detect(basename(f_all), "\\.tif$")
  ]
  if (length(f_fsc) == 0) stop("No FSC tif found in raw_dir.")
  
  # ---- Extract timestamps ----
  get_stamp <- function(x) str_extract(basename(x), "\\d{8}T\\d{6}")
  stamps <- vapply(f_fsc, get_stamp, character(1))
  ok <- !is.na(stamps)
  f_fsc <- f_fsc[ok]
  stamps <- stamps[ok]
  
  dt <- as.POSIXct(stamps, format = "%Y%m%dT%H%M%S", tz = "UTC")
  
  # ---- Temporal filter ----
  keep_time <- dt >= start_dt & dt <= end_dt
  f_fsc <- f_fsc[keep_time]
  dt    <- dt[keep_time]
  if (length(f_fsc) == 0) stop("No FSC files within the requested time range.")
  
  # ---- Spatial prefilter ----
  keep_space <- logical(length(f_fsc))
  for (i in seq_along(f_fsc)) {
    r_meta <- terra::rast(f_fsc[i])
    aoi_r  <- terra::project(aoi, terra::crs(r_meta))
    keep_space[i] <- !is.null(terra::intersect(terra::ext(r_meta), terra::ext(aoi_r)))
  }
  f_fsc <- f_fsc[keep_space]
  dt    <- dt[keep_space]
  if (length(f_fsc) == 0) stop("No FSC tiles intersecting the AOI.")
  
  # ---- Group by UTC day ----
  day_tag <- format(dt, "%Y%m%d")
  files_by_day <- split(f_fsc, day_tag)
  days_unique <- names(files_by_day)
  
  # ---- Mosaic rules ----
  mosaic_rules <- function(v) {
    if (all(is.na(v))) return(NA)
    if (any(v == 205, na.rm = TRUE)) return(205)   # cloud dominates
    vv <- v[v >= 0 & v <= 100]
    if (length(vv) > 0) return(max(vv))
    return(255)                                    # nodata code
  }
  
  template_2154 <- NULL
  n_saved <- 0L
  
  for (d in days_unique) {
    
    if (verbose) message("Processing day: ", d, " (", length(files_by_day[[d]]), " tile(s))")
    
    ras_list <- lapply(files_by_day[[d]], function(f) {
      tryCatch(terra::rast(f), error = function(e) NULL)
    })
    ras_list <- Filter(Negate(is.null), ras_list)
    if (length(ras_list) == 0) next
    
    # --- reference tile + AOI in same CRS ---
    ref <- ras_list[[1]]
    ref_crs <- terra::crs(ref)
    aoi_ref <- terra::project(aoi, ref_crs)
    aoi_bbox <- terra::ext(aoi_ref)
    
    # template aligned to ref grid covering AOI bbox
    template_native <- terra::crop(ref, aoi_bbox, snap = "out")
    
    # --- align each tile to native grid over AOI bbox ---
    ras_list <- lapply(ras_list, function(r) {
      
      if (!terra::same.crs(r, ref)) {
        r <- terra::project(r, ref_crs, method = "near")
      }
      
      if (is.null(terra::intersect(terra::ext(r), aoi_bbox))) return(NULL)
      
      r2 <- terra::crop(r, aoi_bbox, snap = "out")
      r2 <- terra::resample(r2, template_native, method = "near")
      r2
    })
    ras_list <- Filter(Negate(is.null), ras_list)
    if (length(ras_list) == 0) {
      if (verbose) message(" -> skipped (no tile overlaps AOI after alignment).")
      next
    }
    
    # --- mosaic ---
    mos <- if (length(ras_list) == 1) {
      ras_list[[1]]
    } else {
      r_stack <- terra::rast(ras_list)
      terra::app(r_stack, mosaic_rules)
    }
    
    # --- IMPORTANT: crop + mask AOI (sinon tu calcules sur la bbox, pas sur l’AOI) ---
    mos <- terra::crop(mos, aoi_bbox, snap = "out")
    mos <- terra::mask(mos, aoi_ref)
    
    # --- AOI checks (NoData + cloud%) ---
    vals <- terra::values(mos, mat = FALSE)
    vals <- vals[!is.na(vals)]
    
    if (length(vals) == 0 || all(vals == 255)) {
      if (verbose) message(" -> skipped (AOI is 100% NoData for this day).")
      next
    }
    
    vals_valid <- vals[vals != 255]
    if (length(vals_valid) == 0) {
      if (verbose) message(" -> skipped (AOI has no valid pixels for this day).")
      next
    }
    
    cloud_pct <- 100 * sum(vals_valid == 205) / length(vals_valid)
    if (cloud_pct > cloud_thr) {
      if (verbose) message(sprintf(" -> skipped (AOI cloud cover = %.2f%% > %.2f%%).", cloud_pct, cloud_thr))
      next
    }
    
    # --- binarize (keep cloud 205; keep nodata 255) ---
    mos_bin <- mos
    mos_bin[mos <  snow_thr] <- 0
    mos_bin[mos >= snow_thr & mos <= 100] <- 1
    # 205 stays 205; 255 stays 255; outside AOI stays NA
    
    # --- reproject to 2154 + align to common grid ---
    mos_2154 <- terra::project(mos_bin, "EPSG:2154", method = "near")
    if (is.null(template_2154)) {
      template_2154 <- mos_2154
    } else {
      mos_2154 <- terra::resample(mos_2154, template_2154, method = "near")
    }
    
    # --- output folder per CALENDAR YEAR (année classique) ---
    year_tag <- substr(d, 1, 4)
    
    
    out_name <- paste0(site_tag, "_", d, "_SnowCoverBin_2154.tif")
    out_path <- file.path(output_case, out_name)
    
    # (optionnel mais propre) : écrire nodata en 255 via NAflag
    # on convertit 255 -> NA juste avant écriture, et on dit que NA = 255 dans le GeoTIFF
    out_r <- mos_2154
    out_r[out_r == 255] <- NA
    
    terra::writeRaster(
      out_r, out_path, overwrite = TRUE,
      NAflag = 255,
      wopt = list(datatype = "INT1U")
    )
    
    n_saved <- n_saved + 1L
    if (verbose) message(" -> saved: ", normalizePath(out_path))
  }
  
  if (verbose) {
    message("Done. Days saved: ", n_saved, " / ", length(days_unique))
    message("Output folder: ", normalizePath(output_case))
  }
  
  invisible(output_case)
}









#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#
build_dataset_snow_modeling <- function(
    input_snow_image_file,
    topo_file,
    output_file
) {
  ### Description:
  # build_dataset_snow_modeling()
  # This function builds a training dataset for snow-cover modelling from:
  #   - a multi-layer snow cover raster time series (one band per date),
  #   - a stack of topographic variables.
  # Steps:
  #   1) Read the snow stack and use its grid as template.
  #   2) Reproject/resample the topographic stack onto the snow grid.
  #   3) Stack snow and topographic rasters into a single SpatRaster.
  #   4) Convert the stack to a data.frame with x, y, cell and all bands.
  #   5) Remove pixels that have no valid snow information (only NA/255 across all dates).
  #   6) Reshape the snow stack to long format (one row per cell × date).
  #   7) Create a binary snow target:
  #        - 0 = no snow
  #        - 1 = snow
  #        - NA = cloud/shadow/no data (e.g. 205, 255, NA)
  #   8) Save the resulting long-format dataset as an .rds file.
  
  
  ## Check required packages ----
  stopifnot(requireNamespace("terra", quietly = TRUE))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  library(terra)
  library(dplyr)
  library(tidyr)
  
  ## Define snow cover template to align topography ----
  snow_stack    <- terra::rast(input_snow_image_file) 
  snow_template <- snow_stack[[1]]
  
  topo_stack_raw <- terra::rast(topo_file)
  if (!terra::same.crs(topo_stack_raw, snow_template)) {
    topo_on_snow <- terra::project(topo_stack_raw, snow_template, method = "bilinear")
  } else {
    topo_on_snow <- terra::resample(topo_stack_raw, snow_template, method = "bilinear")
  }
  
  # Check 
  stopifnot(
    nrow(snow_stack) == nrow(topo_on_snow),
    ncol(snow_stack) == ncol(topo_on_snow),
    terra::same.crs(snow_stack, topo_on_snow),
    all(terra::ext(snow_stack) == terra::ext(topo_on_snow))
  )
  
  ## Stack snow cover and topographic rasters ----
  full_stack <- c(snow_stack, topo_on_snow)
  snow_names <- names(snow_stack)      # ex. "20221001", "20221006", ...
  topo_names <- names(topo_on_snow)    # ex. "elev", "slope_deg", "tpi", "dah", ...
  
  
  ## 4. Build long-format dataset for modelling ----
  # - extract all pixels (including NA, 205, 255)
  # - remove pixels with no valid snow info (only NA/255 over all dates)
  # - reshape snow stack to long format and create binary snow target
  
  df <- as.data.frame(
    full_stack,
    xy    = TRUE,
    cells = TRUE,
    na.rm = FALSE
  )
  
  no_snow_all <- apply(
    df[, snow_names, drop = FALSE],
    1,
    function(z) all(is.na(z) | z == 255)
  )
  df <- df[!no_snow_all, ]
  
  dataset <- df |>
    tidyr::pivot_longer(
      cols      = all_of(snow_names),
      names_to  = "date_str",
      values_to = "snow_raw"
    ) |>
    dplyr::mutate(
      date = as.Date(date_str, format = "%Y%m%d"),
      snow = dplyr::case_when(
        snow_raw == 0 ~ 0L,
        snow_raw == 1 ~ 1L,
        TRUE          ~ NA_integer_
      )
    ) |>
    dplyr::select(
      cell, x, y,
      date_str, date,
      snow,
      dplyr::all_of(topo_names)
    )
  
  saveRDS(dataset, file = output_file)
  
}