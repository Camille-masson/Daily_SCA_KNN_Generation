






#### 0. LIBRARIES AND CONSTANTS ####
#----------------------------------#

gc()
# Loading configuration
source("config.R")


# This first part is for the product hydroweb
if(TRUE){
#### 1. DOWNLOAD RAWDATA ####
#---------------------------#

if(TRUE){
  # Description:
  # This script automatically downloads the “Theia Sentinel-2 Snow Cover (Level-2B)”
  # product from hydroweb.next.
  #
  # Before running:
  #   1) Create an account on https://hydroweb.next.theia-land.fr
  #   2) Generate an API key in your user settings
  #   3) Copy/paste this key into the API_KEY parameter in R.
  #
  # Parameters to set:
  #   - AREA: study extent (xmin, ymin, xmax, ymax in WGS84 lon/lat)
  #   - START_DATE / END_DATE: time range in UTC
  #
  # All Sentinel-2 tiles intersecting the bbox and dates are downloaded, unzipped,
  # and stored in:
  #   downloads/raw_data/
  # ZIP archives are deleted after extraction. If the same period is requested
  # again, the script skips the download and prints “OK”; tiles already present
  # are not duplicated.
  #
  # Product summary:
  # Theia Sentinel-2 Snow Cover L2B provides fractional snow cover maps derived
  # from Sentinel-2 optical imagery at 20 m spatial resolution, with an effective
  # revisit of ~5 days. Data are delivered on the native Sentinel-2 UTM tiling grid
  # (~110 × 110 km tiles). Each product includes:
  #   - LIS_S2-SNOW-FSC*.tif: fractional snow cover (0 = no snow, 1–100 = % snow,
  #     205 = cloud/shadow, 255 = no data)
  #   - LIS_S2-SNOW-FSC-QCFLAGS*.tif: quality flags
  #   - LIS_S2-SNOW-FSC_*.xml: metadata.
  # Filenames contain the acquisition date/time in UTC.
  
  ## Function et package ----
  library(reticulate)
  library(fs)
  library(stringr)
  source(file.path(functions_dir , "Functions_utility.R"))
  
  
  ## PARAMETERS ----
  AREA <- c(6.359784,
            45.033859,
            6.425705,
            45.086778 )  # xmin,ymin,xmax,ymax lon/lat WGS84
  START_DATE <- "2024-01-01T00:00:00Z"
  END_DATE   <- "2024-12-31T23:59:59Z"
  COLLECTION_ID <- "LIS_FSC_PREOP"
  API_KEY <- "3btR8r5Q6nuV4pqANmBcPUucYcPfRmcw2BuKR79sQSEjSeWloe"  # <-- replace with your hydroweb.next API key
  ONLY_FSC <- FALSE  # TRUE = keep only FSC rasters (ignore QCFLAGS and XML metadata)
  
  ## OUTPUT ----
  # Create the output folder for raw products
  raw_dir <- file.path(downloads_dir, "raw_data")
  if (!dir.exists(downloads_dir)) {
    dir.create(downloads_dir, recursive = TRUE)
  }
  
  ## CODE ----
  
  # Function to download the snow data
  download_hydroweb_snow(
    AREA = AREA, 
    START_DATE = START_DATE,
    END_DATE = END_DATE,
    COLLECTION_ID = COLLECTION_ID,
    API_KEY = API_KEY,
    raw_dir = raw_dir,
    ONLY_FSC = ONLY_FSC
  )
}

  
  
  
  
  
  
  
#### 2. MOSAIC AND CROP TO STUDY AREA ####
#----------------------------------------#

if(TRUE){
  # Description:
  # This part post-processes the raw Theia Sentinel-2 Snow Cover (L2B) tiles downloaded in Part 1.
  # For the selected water year (START_DATE–END_DATE) and AOI shapefile:
  #   1) FSC tiles are filtered by date and by spatial intersection with the AOI.
  #   2) For each acquisition timestamp, all intersecting tiles are mosaicked using
  #      priority rules:
  #        - Values 0–100: keep the highest snow fraction where tiles overlap.
  #        - Cloud (205): dominates any other value.
  #        - No data (255): converted to NA so valid data is kept when available.
  #   3) The daily mosaic is cropped and masked to the AOI.
  #   4) The FSC mosaic is binarized:
  #        - FSC < 70  -> 0 (no/low snow)
  #        - FSC 70–100 -> 1 (snow)
  #        - Cloud (205) and NA remain NA.
  #   5) All binary mosaics are stacked into a multi-layer raster time series.
  #   6) The final stack is reprojected to EPSG:2154 (Lambert-93) with nearest-neighbour
  #      resampling to preserve class values, then saved as a GeoTIFF in outputs/.
  
  
  
  ## Function et package ----
  source(file.path(functions_dir , "Functions_utility.R"))
  
  ## PARAMETERS ----
  STUDY_AREA <- "Roche-Noire" # <- Replace with the name of the study area
  YEAR <- 2024 # WARNING: the water year is defined here from 1 Oct of WATER_YEAR to 30 Sep of WATER_YEAR + 1.
  
  
  ## INPUT ----
  # Folder for raw data downloaded in Part 1
  raw_dir <- file.path(downloads_dir, "raw_data")
  
  # Create folder to store the Area Of Interest (AOI) shapefile.
  # Place in this folder a shapefile named: AOI_<STUDY_AREA>.shp
  # Example: AOI_Cayolle.shp
  # The AOI can be provided in EPSG:2154 (Lambert-93) or any other CRS.
  # It will be reprojected automatically to match the raster CRS during cropping.
  AOI_case <- file.path(inputs_dir, "Area_Of_Interest")
  if (!dir.exists(AOI_case)) {
    dir.create(AOI_case, recursive = TRUE)
  }
  AOI_file <- file.path (AOI_case, paste0("AOI_",STUDY_AREA,".shp"))
  
  
  ## OUTPUT ----
  # Create output folder for the processed Snow Cover product
  output_snow_case <- file.path(output_dir, "1. Snow_cover_by_site")
  if (!dir.exists(output_snow_case)) {dir.create(output_snow_case, recursive = TRUE)}
  
  output_case <- file.path(output_snow_case, paste0(STUDY_AREA))
  if (!dir.exists(output_case)) {dir.create(output_case, recursive = TRUE)}
  
  
  ## CODE ----
  mosaic_crop_save_snow_binary_daily(
    raw_dir     = raw_dir,
    AOI_file    = AOI_file,
    YEAR        = YEAR,
    output_dir  = output_dir,
    STUDY_AREA  = STUDY_AREA,
    cloud_thr   = 10,
    snow_thr    = 70,
    verbose     = TRUE
  )
}


}



#---------------------------#
#### 1. DOWNLOAD RAWDATA ####
#---------------------------#








#---------------------------#
#### 1. DOWNLOAD RAWDATA ####
#---------------------------#

if (TRUE) {
  
  ## Packages ----
  library(httr2)
  library(jsonlite)
  library(fs)
  
  ## PARAMETERS ----
  AREA <- c(6.304012148350035,
            45.02484125162051,
            6.43194986658607,
            45.11260606895114)  # xmin,ymin,xmax,ymax lon/lat WGS84
  
  # ✅ Tu modifies ça à la main selon tes besoins
  START_DATE <- "2025-11-01T00:00:00Z"
  END_DATE   <- "2025-12-31T23:59:59Z"
  
  # Sentinel Hub classic OAuth (realm main)
  CLIENT_ID     <- "900acc8d-fd83-4c7b-9f09-bc0d8dbaf358"
  CLIENT_SECRET <- "8muCNdPQ5EhAXEO0w7yrv0RwUiPUuDaw"
  
  AUTH_URL    <- "https://services.sentinel-hub.com/auth/realms/main/protocol/openid-connect/token"
  CATALOG_URL <- "https://services-uswest2.sentinel-hub.com/api/v1/catalog/1.0.0/search"
  PROCESS_URL <- "https://services-uswest2.sentinel-hub.com/api/v1/process"
  
  # Output
  downloads_dir <- "downloads"
  raw_dir <- file.path(downloads_dir, "raw_data")
  dir_create(raw_dir, recurse = TRUE)
  
  resolution_m <- 30
  
  ## Helpers ----
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  get_sh_token <- function(client_id, client_secret, auth_url) {
    resp <- request(auth_url) |>
      req_body_form(
        grant_type = "client_credentials",
        client_id = client_id,
        client_secret = client_secret
      ) |>
      req_perform()
    
    if (resp_status(resp) >= 400) {
      msg <- tryCatch(resp_body_string(resp), error = function(e) "")
      stop("Auth error ", resp_status(resp), "\n", msg)
    }
    
    resp_body_json(resp)$access_token
  }
  
  # Approx bbox (WGS84) -> pixel size for ~30 m
  bbox_to_size_wgs84 <- function(bbox, resolution_m = 30) {
    xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]
    mid_lat <- (ymin + ymax) / 2
    
    m_per_deg_lat <- 111320
    m_per_deg_lon <- 111320 * cos(mid_lat * pi/180)
    
    width_m  <- (xmax - xmin) * m_per_deg_lon
    height_m <- (ymax - ymin) * m_per_deg_lat
    
    width_px  <- max(1, ceiling(width_m  / resolution_m))
    height_px <- max(1, ceiling(height_m / resolution_m))
    
    c(width_px, height_px)
  }
  
  # ✅ Evalscript NDSI-ready + QA
  # Bands order:
  # 1) Green
  # 2) SWIR1
  # 3) QA  (HLS quality mask)
  # 4) dataMask
  build_green_swir1_qa_evalscript <- function() {
    lines <- c(
      "//VERSION=3",
      "function setup() {",
      "  return {",
      '    input: ["Green", "SWIR1", "QA", "dataMask"],',
      '    output: { bands: 4, sampleType: "FLOAT32" }',
      "  };",
      "}",
      "",
      "function evaluatePixel(sample) {",
      "  return [",
      "    sample.Green,",
      "    sample.SWIR1,",
      "    sample.QA,",
      "    sample.dataMask",
      "  ];",
      "}"
    )
    paste(lines, collapse = "\n")
  }
  
  # Catalog simple (sans distinct)
  get_hls_dates_manual <- function(token, catalog_url, area, start_date, end_date) {
    
    body <- list(
      bbox = as.list(area),
      datetime = paste0(start_date, "/", end_date),
      collections = list("hls"),
      limit = 100
    )
    
    resp <- request(catalog_url) |>
      req_headers(
        Authorization = paste("Bearer", token),
        `Content-Type` = "application/json"
      ) |>
      req_body_json(body, auto_unbox = TRUE) |>
      req_perform()
    
    if (resp_status(resp) >= 400) {
      msg <- tryCatch(resp_body_string(resp), error = function(e) "")
      stop("Catalog error ", resp_status(resp), "\n", msg)
    }
    
    out <- resp_body_json(resp)
    feats <- out$features %||% list()
    
    if (length(feats) == 0) return(character(0))
    
    dts <- vapply(feats, function(f) f$properties$datetime %||% NA_character_, character(1))
    sort(unique(substr(dts[!is.na(dts)], 1, 10)))
  }
  
  # Process one date
  request_hls_green_swir1_qa_one_date <- function(token, process_url, area, date_str,
                                                  resolution_m = 30,
                                                  constellation = NULL) {
    
    size <- bbox_to_size_wgs84(area, resolution_m)
    width <- size[1]; height <- size[2]
    
    evalscript <- build_green_swir1_qa_evalscript()
    
    dataFilter <- list(
      timeRange = list(
        from = paste0(date_str, "T00:00:00Z"),
        to   = paste0(date_str, "T23:59:59Z")
      )
    )
    
    # ✅ Option utile : en 2013, tu peux forcer LANDSAT
    if (!is.null(constellation)) {
      dataFilter$constellation <- constellation
    }
    
    body <- list(
      input = list(
        bounds = list(
          properties = list(
            crs = "http://www.opengis.net/def/crs/OGC/1.3/CRS84"
          ),
          bbox = as.list(area)
        ),
        data = list(list(
          type = "hls",
          dataFilter = dataFilter
        ))
      ),
      output = list(
        width = width,
        height = height,
        responses = list(list(
          identifier = "default",
          format = list(type = "image/tiff")
        ))
      ),
      evalscript = evalscript
    )
    
    resp <- request(process_url) |>
      req_headers(
        Authorization = paste("Bearer", token),
        `Content-Type` = "application/json"
      ) |>
      req_body_json(body, auto_unbox = TRUE) |>
      req_perform()
    
    if (resp_status(resp) >= 400) {
      msg <- tryCatch(resp_body_string(resp), error = function(e) "")
      stop("Process error ", resp_status(resp), " for date ", date_str, "\n", msg)
    }
    
    resp_body_raw(resp)
  }
  
  ## =========================
  ## MAIN (manuel)
  ## =========================
  
  token <- get_sh_token(CLIENT_ID, CLIENT_SECRET, AUTH_URL)
  
  dates <- get_hls_dates_manual(
    token = token,
    catalog_url = CATALOG_URL,
    area = AREA,
    start_date = START_DATE,
    end_date = END_DATE
  )
  
  message("Nombre de dates HLS trouvées: ", length(dates))
  
  if (length(dates) == 0) {
    message("Aucune date trouvée.")
  } else {
    
    for (d in dates) {
      
      out_name <- paste0("HLS_GREEN_SWIR1_QA_", gsub("-", "", d), "_30m_approx_CRS84.tif")
      out_path <- file.path(raw_dir, out_name)
      
      if (file_exists(out_path)) next
      
      message("Téléchargement: ", d)
      
      # ✅ Recommandé pour 2013:
      # constellation = "LANDSAT"
      # Pour 2018+ tu peux mettre NULL pour laisser tout HLS
      tif_raw <- request_hls_green_swir1_qa_one_date(
        token = token,
        process_url = PROCESS_URL,
        area = AREA,
        date_str = d,
        resolution_m = resolution_m,
        constellation = if (substr(d, 1, 4) == "2013") "LANDSAT" else NULL
      )
      
      writeBin(tif_raw, out_path)
    }
    
    message("Terminé. Données dans: ", normalizePath(raw_dir))
  }
}

















#---------------------------#
#### 2. BUILD SNOW BINARY ###
#---------------------------#

if (TRUE) {
  
  ## Packages ----
  library(terra)
  library(fs)
  library(stringr)
  
  ## INPUT / OUTPUT ----
  raw_dir <- file.path("downloads", "raw_data")
  out_dir <- file.path("downloads", "snow_binary")
  dir_create(out_dir, recurse = TRUE)
  
  # On cible uniquement tes stacks utiles
  files <- dir_ls(raw_dir, glob = "*.tif")
  files <- files[str_detect(basename(files), "HLS_GREEN_SWIR1_QA")]
  
  message("Nb de fichiers trouvés: ", length(files))
  
  ## Parameters ----
  ndsi_thresh   <- 0.4
  mask_adjacent <- TRUE
  mask_cirrus   <- TRUE
  invalid_value <- 205   # ✅ ton code invalid
  
  ## Bit test helper for SpatRaster ----
  bit_is_set_raster <- function(QA_r, bit) {
    terra::app(QA_r, fun = function(x) {
      x <- as.integer(round(x))
      bitwAnd(x, bitwShiftL(1L, bit)) != 0L
    })
  }
  
  ## Core function ----
  make_snow_binary_one <- function(f) {
    
    r <- rast(f)
    n <- nlyr(r)
    
    # Expected order:
    # 1) Green
    # 2) SWIR1
    # 3) QA
    # 4) dataMask (optional)
    Green <- r[[1]]
    SWIR1 <- r[[2]]
    QA    <- r[[3]]
    
    # dataMask optional
    if (n >= 4) {
      dataMask <- r[[4]]
    } else {
      dataMask <- Green * 0 + 1
    }
    
    # QA bits
    cloud  <- bit_is_set_raster(QA, 1)
    shadow <- bit_is_set_raster(QA, 3)
    adj    <- bit_is_set_raster(QA, 2)
    cir    <- bit_is_set_raster(QA, 0)
    
    invalid <- cloud | shadow
    if (mask_adjacent) invalid <- invalid | adj
    if (mask_cirrus)   invalid <- invalid | cir
    
    # add no-data
    invalid <- invalid | (dataMask == 0)
    
    # NDSI
    ndsi <- (Green - SWIR1) / (Green + SWIR1)
    
    # ✅ Output coding:
    # 205 = invalid
    # 1   = snow
    # 0   = no snow
    snowbin <- ifel(invalid, invalid_value,
                    ifel(ndsi >= ndsi_thresh, 1, 0))
    
    names(snowbin) <- "snow_binary"
    snowbin
  }
  
  ## Loop ----
  for (f in files) {
    
    base <- basename(f)
    out_name <- str_replace(base, "HLS_GREEN_SWIR1_QA", "HLS_SNOWBIN_NDSI")
    out_path <- file.path(out_dir, out_name)
    
    if (file_exists(out_path)) {
      message("Skip (exists): ", out_name)
      next
    }
    
    message("Processing: ", base)
    
    snow <- make_snow_binary_one(f)
    
    writeRaster(
      snow,
      out_path,
      overwrite = TRUE,
      wopt = list(datatype = "INT1U") # 0..255
    )
  }
  
  message("Terminé. Sorties dans: ", normalizePath(out_dir))
}





#-----------------------------------------------#
#### 3. FILTER SCENES BY % INVALID (205)     ####
#-----------------------------------------------#

if (TRUE) {
  
  library(terra)
  library(fs)
  library(stringr)
  
  ## PARAMETERS ----
  invalid_value <- 205
  max_invalid_ratio <- 0.2  # 5%
  
  ## INPUT / OUTPUT ----
  in_dir  <- file.path("downloads", "snow_binary")
  out_dir <- file.path("downloads", "snow_binary_filtered")
  dir_create(out_dir, recurse = TRUE)
  
  # Fichiers snow binaires
  files <- dir_ls(in_dir, glob = "*.tif")
  files <- files[str_detect(basename(files), "HLS_SNOWBIN_NDSI")]
  
  message("Nb de scènes candidates: ", length(files))
  
  # Helper: extraire la date depuis le nom
  extract_date <- function(x) {
    m <- str_match(basename(x), "(\\d{8})")
    if (is.na(m[1,2])) return(NA_character_)
    m[1,2]
  }
  
  # ✅ Helper: calcul robuste du ratio de 205
  compute_invalid_ratio <- function(f) {
    r <- rast(f)
    
    total <- terra::global(!is.na(r), "sum", na.rm = TRUE)[1,1]
    if (is.na(total) || total == 0) return(1)
    
    inv <- terra::global(r == invalid_value, "sum", na.rm = TRUE)[1,1]
    if (is.na(inv)) inv <- 0
    
    as.numeric(inv / total)
  }
  
  ## LOOP ANALYSE ----
  report <- data.frame(
    file = files,
    date_yyyymmdd = vapply(files, extract_date, character(1)),
    invalid_ratio = NA_real_,
    keep = FALSE,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(files)) {
    f <- files[i]
    message("[", i, "/", length(files), "] Analyse: ", basename(f))
    
    ratio <- compute_invalid_ratio(f)
    
    report$invalid_ratio[i] <- ratio
    report$keep[i] <- is.finite(ratio) && ratio <= max_invalid_ratio
  }
  
  # CSV de contrôle
  write.csv(report,
            file = file.path(out_dir, "snow_invalid_report.csv"),
            row.names = FALSE)
  
  ## KEEP FILES ----
  keep_files <- report$file[report$keep]
  message("Scènes conservées (<=5% 205): ", length(keep_files))
  message("Scènes rejetées: ", sum(!report$keep))
  
  for (f in keep_files) {
    file_copy(f, file.path(out_dir, basename(f)), overwrite = TRUE)
  }
  
  ## learning_date ----
  learning_date <- report$date_yyyymmdd[report$keep]
  learning_date <- learning_date[!is.na(learning_date)]
  learning_date <- paste0(substr(learning_date, 1, 4), "-",
                          substr(learning_date, 5, 6), "-",
                          substr(learning_date, 7, 8))
  learning_date <- sort(unique(learning_date))
  
  message("Nb de learning_date: ", length(learning_date))
  
  ## STACK léger via VRT ----
  filtered_files <- dir_ls(out_dir, glob = "*.tif")
  filtered_files <- filtered_files[str_detect(basename(filtered_files), "HLS_SNOWBIN_NDSI")]
  
  if (length(filtered_files) > 0) {
    snow_stack <- rast(filtered_files)
    
    # Option VRT
    writeRaster(snow_stack,
                filename = file.path(out_dir, "snow_stack.vrt"),
                overwrite = TRUE)
    
    message("Stack prêt: ", nlyr(snow_stack), " couches.")
  } else {
    message("Aucun fichier filtré à stacker.")
  }
  
  ## Export objets en session
  assign("learning_date", learning_date, envir = .GlobalEnv)
  if (exists("snow_stack")) assign("snow_stack", snow_stack, envir = .GlobalEnv)
  
  message("Terminé. Dossier: ", normalizePath(out_dir))
}
