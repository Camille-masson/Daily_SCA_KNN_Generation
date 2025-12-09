#### 0. LIBRARIES AND CONSTANTS ####
#----------------------------------#

gc()
# Loading configuration
source("config.R")


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
  AREA <- c(5.230422275471649,
            43.94091874168251,
            7.441074966219298,
            45.534980058544996)  # xmin,ymin,xmax,ymax lon/lat WGS84
  START_DATE <- "2020-02-01T00:00:00Z"
  END_DATE   <- "2020-02-29T23:59:59Z"
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
  STUDY_AREA <- "Cayolle" # <- Replace with the name of the study area
  WATER_YEAR <- 2021 # WARNING: the water year is defined here from 1 Oct of WATER_YEAR to 30 Sep of WATER_YEAR + 1.
  
  WY_LABEL   <- sprintf("%d-%d", WATER_YEAR, WATER_YEAR + 1)
  START_DATE <- sprintf("%d-10-01T00:00:00Z", WATER_YEAR)
  END_DATE   <- sprintf("%d-09-30T23:59:59Z", WATER_YEAR + 1)
  
  
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
  output_case <- file.path(output_dir, "1. Snow_cover_by_site")
  if (!dir.exists(output_case)) {
    dir.create(output_case, recursive = TRUE)
  }
  
  date_tag <- paste0(gsub("[^0-9]", "", substr(START_DATE, 1, 10)), "_",
                     gsub("[^0-9]", "", substr(END_DATE,   1, 10)))
  
  output_file <- file.path(output_case,paste0("Theia_Sentinel-2_Snow_Cover_Level-2B_",STUDY_AREA,
                                              "_WY", WY_LABEL, "_", date_tag,".tif"))
  
  ## CODE ----
  mosaic_crop_stack_snow_binary(
    raw_dir     = raw_dir,
    AOI_file    = AOI_file,
    START_DATE  = START_DATE,
    END_DATE    = END_DATE,
    output_file = output_file,
    verbose     = TRUE
  )
}



