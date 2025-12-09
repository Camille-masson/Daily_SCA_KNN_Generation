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







#------------------------------------------------------------------------------#
################################################################################
#------------------------------------------------------------------------------#
mosaic_crop_stack_snow_binary <- function(
    raw_dir,
    AOI_file,
    START_DATE,
    END_DATE,
    output_file,
    verbose = TRUE
) {
  
  ### Description:
  # mosaic_crop_stack_snow_binary()
  # This function post-processes raw Theia Sentinel-2 Snow Cover (L2B) FSC tiles
  # downloaded from hydroweb.next. For a given AOI shapefile and time range:
  #   1) It selects FSC GeoTIFFs in raw_dir, filters them by acquisition date
  #      (START_DATE–END_DATE) and by spatial overlap with the AOI.
  #   2) Files are grouped by UTC day so that only one output layer is produced per day.
  #   3) For each day, all intersecting tiles are aligned on a common native grid,
  #      then mosaicked pixel-wise using priority rules:
  #        - cloud (205) dominates,
  #        - otherwise the highest FSC value (0–100) is kept,
  #        - 255 (NoData) is kept only when no valid data exists.
  #   4) The daily mosaic is cropped/masked to the AOI; days fully 255 over the AOI
  #      are skipped to avoid useless layers.
  #   5) FSC is binarized (<70→0, 70–100→1) while keeping 205 and 255 unchanged.
  #   6) Each daily layer is reprojected to EPSG:2154 (nearest neighbour),
  #      aligned to a common 2154 grid, stacked into a time series, and saved to output_file.
  
  ## Check required packages ----
  stopifnot(requireNamespace("terra", quietly = TRUE))
  stopifnot(requireNamespace("fs", quietly = TRUE))
  stopifnot(requireNamespace("stringr", quietly = TRUE))
  library(terra)
  library(fs)
  library(stringr)
  
  ## Validate input parameters ----
  if (!file.exists(AOI_file)) stop("AOI_file not found.")
  start_dt <- as.POSIXct(START_DATE, tz = "UTC")
  end_dt   <- as.POSIXct(END_DATE,   tz = "UTC")
  if (is.na(start_dt) || is.na(end_dt)) stop("Invalid START_DATE / END_DATE (ISO UTC expected)")
  aoi <- terra::vect(AOI_file)
  if (is.na(terra::crs(aoi))) stop("AOI shapefile has no CRS")
  
  ## List FSC files ----
  # List all raw files in raw_dir, then keep only FSC GeoTIFFs
  # (exclude QCFLAGS and XML/other formats). Stop if none are found.
  f_all <- fs::dir_ls(raw_dir, recurse = FALSE, type = "file")
  f_fsc <- f_all[
    str_detect(basename(f_all), "FSC") &
      !str_detect(basename(f_all), "QCFLAGS") &
      str_detect(basename(f_all), "\\.tif$")
  ]
  if (length(f_fsc) == 0) stop("No FSC tif found in raw_dir.")
  
  ## Extract acquisition timestamps ----
  # Extract Sentinel-2 acquisition timestamps (YYYYMMDDTHHMMSS) from filenames,
  # drop files without a valid timestamp, then convert to UTC datetimes for filtering/grouping.
  get_stamp <- function(x) str_extract(basename(x), "\\d{8}T\\d{6}")
  stamps <- vapply(f_fsc, get_stamp, character(1))
  ok <- !is.na(stamps)
  f_fsc <- f_fsc[ok]
  stamps <- stamps[ok]
  
  dt <- as.POSIXct(stamps, format = "%Y%m%dT%H%M%S", tz = "UTC")
  
  ## Temporal filter ----
  # Keep only FSC files whose acquisition datetime falls within START_DATE–END_DATE; stop if none remain.
  keep_time <- dt >= start_dt & dt <= end_dt
  f_fsc <- f_fsc[keep_time]
  dt    <- dt[keep_time]
  if (length(f_fsc) == 0) stop("No FSC files within the requested time range.")
  
  ## Spatial filter (quick) ----
  # Quick spatial pre-filter: keep only FSC tiles whose extent overlaps the AOI
  # (after reprojecting AOI to each tile CRS), and stop if none intersect.
  keep_space <- logical(length(f_fsc))
  for (i in seq_along(f_fsc)) {
    r_meta <- terra::rast(f_fsc[i])
    aoi_r  <- terra::project(aoi, terra::crs(r_meta))
    inter  <- terra::intersect(terra::ext(r_meta), terra::ext(aoi_r))
    keep_space[i] <- !is.null(inter)
  }
  f_fsc <- f_fsc[keep_space]
  dt    <- dt[keep_space]
  if (length(f_fsc) == 0) stop("No FSC tiles intersecting the AOI.")
  
  ## Prepare inputs for daily mosaics ----
  # Group FSC files by UTC day so there is only one mosaic/stack layer per day.
  # Create a temporary folder to store daily processed rasters before stacking,
  # and initialize trackers (kept days, temp files, and a 2154 template grid).
  # Define pixel-wise mosaic rules: 205 (cloud) dominates, else max of 0–100 wins,
  # and 255 is kept only when no valid data is available.
  day_tag <- format(dt, "%Y%m%d")
  files_by_day <- split(f_fsc, day_tag)
  days_unique <- names(files_by_day)
  
  
  tmp_dir <- file.path(dirname(output_file), "tmp_layers")
  fs::dir_create(tmp_dir)
  on.exit(fs::dir_delete(tmp_dir), add = TRUE)
  
  tmp_layers <- character(0)
  kept_days  <- character(0)
  template_2154 <- NULL
  
  ## Rules for mosaic
  mosaic_rules <- function(v) {
    if (all(is.na(v))) return(255)
    if (any(v == 205, na.rm = TRUE)) return(205)
    vv <- v[v >= 0 & v <= 100]
    if (length(vv) > 0) return(max(vv))
    return(255)
  }
  
  ## Process each DAY ----
  # Loop over each UTC day to build one final layer per day:
  # - read daily FSC tiles safely
  # - choose first tile as reference, reproject AOI to its CRS
  # - create an AOI-aligned reference grid (template)
  # - reproject/crop/resample/mask each tile to the same grid (skip non-overlapping tiles)
  # - apply custom pixel-wise mosaic rules (205 > max 0–100 > 255)
  # - skip the day if AOI is fully 255 (NoData)
  # - binarize FSC (<70=0, 70–100=1, keep 205 & 255)
  # - reproject to EPSG:2154 and resample to a common 2154 template for stacking
  # - write the daily result as a temporary layer for the final stack
  for (d in days_unique) {
    
    if (verbose) message("Processing day: ", d, " (", length(files_by_day[[d]]), " tile(s))")
    
    ras_list <- lapply(files_by_day[[d]], function(f) {
      tryCatch(terra::rast(f), error = function(e) NULL)
    })
    ras_list <- Filter(Negate(is.null), ras_list)
    if (length(ras_list) == 0) next
    
    ## --- REF + AOI in REF CRS ---
    ref   <- ras_list[[1]]
    ref_crs <- terra::crs(ref)
    aoi_ref <- terra::project(aoi, ref_crs)
    aoi_bbox <- terra::ext(aoi_ref)
    
    ## template aligned on ref grid, covering AOI
    template_native <- terra::crop(ref, aoi_bbox, snap = "out")
    
    ## --- FORCE SAME GRID BEFORE MOSAIC (with overlap check) ---
    ras_list <- lapply(ras_list, function(r) {
      
      # reproject tile to ref CRS if needed
      if (!terra::same.crs(r, ref)) {
        r <- terra::project(r, ref_crs, method = "near")
      }
      
      # if still no overlap -> skip safely
      if (is.null(terra::intersect(terra::ext(r), aoi_bbox))) {
        return(NULL)
      }
      
      # crop/resample/mask on ref grid
      r2 <- terra::crop(r, aoi_bbox, snap = "out")
      r2 <- terra::resample(r2, template_native, method = "near")
      
      r2
    })
    
    ras_list <- Filter(Negate(is.null), ras_list)
    if (length(ras_list) == 0) {
      if (verbose) message(" -> skipped (no tile overlaps AOI after alignment).")
      next
    }
    
    ## Mosaic per day with rules ----
    mos <- if (length(ras_list) == 1) {
      ras_list[[1]]
    } else {
      r_stack <- terra::rast(ras_list)
      terra::app(r_stack, mosaic_rules)
    }
    
    mos[is.na(mos)] <- 255
    
    ## Skip day if AOI is 100% NoData ----
    vals <- terra::values(mos, mat = FALSE)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0 || all(vals == 255)) {
      if (verbose) message(" -> skipped (AOI is 100% NoData for this day).")
      next
    }
    
    ## Binarize FSC (keep 205 & 255) ----
    mos_bin <- mos
    mos_bin[mos < 70] <- 0
    mos_bin[mos >= 70 & mos <= 100] <- 1
    mos_bin[is.na(mos_bin)] <- 255
    
    ## Reproject to EPSG:2154 and align stack grid ----
    mos_2154 <- terra::project(mos_bin, "EPSG:2154", method = "near")
    
    if (is.null(template_2154)) {
      template_2154 <- mos_2154
    } else {
      mos_2154 <- terra::resample(mos_2154, template_2154, method = "near")
    }
    
    ## Write temp layer ----
    tmp_file <- file.path(tmp_dir, paste0("SNOWBIN_", d, ".tif"))
    terra::writeRaster(mos_2154, tmp_file, overwrite = TRUE,
                       wopt = list(datatype = "INT1U"))
    
    tmp_layers <- c(tmp_layers, tmp_file)
    kept_days  <- c(kept_days, d)
  }
  
  ## Stack and save ----
  if (length(tmp_layers) == 0) stop("All days were 100% NoData over the AOI. Nothing to stack.")
  
  if (verbose) message("Stacking layers...")
  stack_r <- terra::rast(tmp_layers)
  names(stack_r) <- kept_days
  
  terra::writeRaster(stack_r, output_file, overwrite = TRUE,
                     wopt = list(datatype = "INT1U"))
  
  if (verbose) {
    message("Output saved to: ", normalizePath(output_file))
    message("Layers kept: ", terra::nlyr(stack_r),
            " / ", length(days_unique))
  }
  
  invisible(output_file)
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