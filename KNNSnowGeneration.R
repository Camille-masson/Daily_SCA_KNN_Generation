#### 0. LIBRARIES & CONSTANTS ####
#--------------------------------#

gc()
# Loading configuration
source("config.R")

site <- "Roche-Noire"



#### 1. Optimisation of parameters and alpha for KNN ####
#-------------------------------------------------------#
if (TRUE){
  
  # Description: 
  # Une partie provisoire qui lit directement les données téléchargé par Philippe 
  # via anciennement Théïa (maintetant hydroweb) qui cumule les images satellite 
  # traiter par Simon Gascouin en binaire (snow/no snow) regroupant pour la plage 
  # temporelle 2013-2023 les satellite Landsat 7-8, Sentinelle 2a. Formant un pool 
  # de 137 images viable sur cette fenètre (viable étant définis comme avec moins 
  # de 5% de nuage sur la tuile)
  #
  # A therme cette partie devra etre complètement remodelé pour fusionner avec la 
  # le script de ce meme projet : SatelliteSnowDataAcquisition.R qui téléchargé via
  # l'API hydroweb les images satellites de n'importe quelle zones d'interets ! 
  #
  # De plus pour les données climatiques d'entrée (Tmin; Tmax; Precipitation) elles
  # sont pour le moment basé sur la station météo fluxalpes dans le site d'entrainement 
  # de Roche noir et devra a therme peut etre utiliser les données température issue
  # de la réanalyse S2M (crocus)
  
  
  
  ## Function et package ----
  
  library(terra)
  library(zoo)
  library(dplyr)
  library(rgenoud)
  library(parallel)
  
  source(file.path(functions_dir , "Functions_KNN.R"))
  
  ## INPUT ----
  # Create folder to store the Area Of Interest (AOI) shapefile.
  KNN_case <- file.path(inputs_dir, "KNN_Training_Data")
  if (!dir.exists(KNN_case)) {
    dir.create(KNN_case, recursive = TRUE)
  }
  
  # Stack of binary snow cover maps (0/1) (Sentienlle / Landsat)
  SCA_file <- file.path (KNN_case, "SCA13to23.L1s.tif")
  
  # Acquisition metadata
  fSCA_file <- file.path (KNN_case, "THE.fSCA.csv")
  
  # Daily meteorological data (FLUXALP)
  Fluxalp_file <- file.path (KNN_case, "FLUXALP.daily.csv")
  
  
  ## OUTPUT ----
  output_case <- file.path(output_dir, "3. KNN_Snow_Generation")
  if (!dir.exists(output_case)) {
    dir.create(output_case, recursive = TRUE)
  }
  
  output_equation_case <- file.path(output_case, "equation_parameters")
  if (!dir.exists(output_equation_case)) {
    dir.create(output_equation_case, recursive = TRUE)
  }
  
  parameters_equation <- file.path(output_equation_case, paste0("Equation_parameters_",site,".csv"))
  
  
  ## CODE ----
  ## 1.1 Data Formatting ----
  if (TRUE){
    
    # Description :
    # Prepare all inputs needed for the KNN climate-distance workflow. This block
    # turns raw satellite snow maps + acquisition metadata + daily meteorology into
    # standardized, lagged climate predictors and pairwise distance components between
    # learning dates.
    # We build a 0–90 day climate cube (short/long temporal neighborhoods) for each
    # acquisition date, rescale variables to [0,1], then compute Manhattan distances
    # per lag and per variable to feed the wshort / wlong / K optimisation steps.
    
    # Load snow cover
    img <- terra::rast(SCA_file)   # SpatRaster : nPixels x nDates
    
    # Load fSCA metadata, filter WS == "TOT" and build acquisition dates
    df <- read.csv(fSCA_file)
    df <- df[df$WS == "TOT", ]                                                                              # Modifier plus tard pour faire directement le match avec le stack img plutot
    df$DATE <- as.Date(paste(df$YEAR, df$DOY, sep = " "), format = "%Y %j")                                 # que de ce baser sur WS
    
    # Load daily FLUXALP meteorological data and convert DATE to Date class
    FXA <- read.csv(Fluxalp_file)
    FXA$DATE <- as.Date(FXA$DATE)
    
    # Linearly interpolate missing FLUXALP values (TN, TX, SNOW, NDVI, etc.)
    for (i in 2:6) {
      FXA[, i] <- zoo::na.approx(FXA[, i])
    }
    
    # Join satellite acquisition dates with corresponding daily FLUXALP meteorology
    df1 <- dplyr::left_join(df,FXA,by=c("DATE"))
    
    # Build 0–90 day meteorological lag cube for each acquisition date
    TMP <- list()
    for(i in 0:90){
      print(i)
      tmp          <- df1$DATE-i
      TMP[[i+1]]   <- FXA[match(tmp,FXA$DATE),c("TN","TX","SNOW","NDVI")]
    }
    
    # Rescale input variable from 0 to 1 (x-tmin/max-min)
    tmp   <- do.call(rbind.data.frame,lapply(TMP,f.MIN))
    MIN   <- t(apply(tmp,2,min,na.rm=T))
    colnames(MIN) <- c("TN","TX","SNOW","NDVI")
    tmp   <- do.call(rbind.data.frame,lapply(TMP,f.MAX))
    MAX   <- t(apply(tmp,2,max,na.rm=T))
    colnames(MAX) <- c("TN","TX","SNOW","NDVI")
    
    INPUT <- lapply(TMP,f.SCALE)
    names(INPUT) <- paste0("DAY_minus",0:90)
    
    # apply this to each lag; list of 91 matrices
    DMfull <- lapply(INPUT, f.MAN)
    # DMfull[[lag]] : n_pairs x n_vars  (col1 = TN, col2 = TX, col3 = SNOW, col4 = NDVI)
    dim(DMfull[[1]])
    
    # Snow Cover maps as matrix
    # img: SpatRaster (nPixels x nDates)
    SCA_mat <- terra::values(img, mat = TRUE)  # nPixels x nDates
    n_dates <- ncol(SCA_mat)
    
    # check
    stopifnot(n_dates == nrow(df1))
  }
  
  ## 1.2 Optimize climate temporal windows (wshort & wlong) ----
  if (TRUE){
  
    # Description:
    # This block tunes the temporal windows used to compare climate histories
    # between satellite acquisition dates.
    #
    # A grid of candidate window lengths is created:
    # - wshort from 1 to 7 days
    # - wlong  from 1 to 90 days
    #
    # For each (wshort, wlong) pair:
    # - build_climate_distance() aggregates the precomputed Manhattan distances
    #   (DMfull) over the selected short and long lag ranges,
    #   averages by window length, and returns a full n_dates x n_dates
    #   climate distance matrix.
    # - eval_w() runs a leave-one-out test with K = 1:
    #   the target date is removed from training, the closest learning date
    #   in climate space is selected, and its snow map is used as prediction.
    #   TP/TN/FP/FN are accumulated across all dates.
    #
    # Global OA and Kappa are computed from the aggregated counts,
    # and epsilon = OA + Kappa is stored for each parameter pair.
    #
    # The best (wshort, wlong) is the combination that maximizes epsilon.
    
    wshort_range <- 1:7
    wlong_range  <- 1:90   
    
    grid <- expand.grid(
      wshort = wshort_range,
      wlong  = wlong_range
    )
    
    grid$OA      <- NA_real_
    grid$kappa   <- NA_real_
    grid$epsilon <- NA_real_
    
    K0 <- 1  # K fixé pour cette première étape
    
    for (r in seq_len(nrow(grid))) {
      ws <- grid$wshort[r]
      wl <- grid$wlong[r]
      cat("Combo", r, "/", nrow(grid), " : wshort =", ws, " wlong =", wl, "\n")
      
      res <- eval_w(ws, wl)
      grid$OA[r]      <- res$OA
      grid$kappa[r]   <- res$kappa
      grid$epsilon[r] <- res$eps
    }
    
    # Meilleur couple selon OA + Kappa
    best <- grid[which.max(grid$epsilon), ]
    print(best)
  }
  
  ## 1.3 Optimize K (wshort & wlong fixed) ----
  if (TRUE){
    # Description :
    # This block selects the optimal number of neighbours K once the best
    # temporal windows (wshort, wlong) have been fixed.
    #
    # A single climate distance matrix (D_best) is built using the chosen
    # wshort_best and wlong_best.
    #
    # K is then tested over a predefined range (here 1:50).
    # For each K, eval_K() performs a leave-one-out reconstruction:
    # - the target date is removed from the learning set,
    # - the K closest learning dates in climate space are selected,
    # - their snow maps are averaged pixel-wise,
    # - a majority vote (>= 0.5) produces the binary prediction.
    #
    # TP/TN/FP/FN are accumulated across all dates, then global OA, Kappa
    # and epsilon = OA + Kappa are computed.
    #
    # The best K is the value that maximizes epsilon.
    
    wshort_best <- best$wshort   
    wlong_best  <- best$wlong    
    cat("Using wshort =", wshort_best, " wlong =", wlong_best, "\n")
    
    D_best <- build_climate_distance(wshort_best, wlong_best)
    
    K_range <- 1:50
    K_res <- data.frame(
      K      = K_range,
      OA     = NA_real_,
      kappa  = NA_real_,
      epsilon = NA_real_
    )
    
    for (kk in K_range) {
      resK <- eval_K(kk, D_best, SCA_mat)
      K_res$OA[K_res$K == kk]      <- resK$OA
      K_res$kappa[K_res$K == kk]   <- resK$kappa
      K_res$epsilon[K_res$K == kk] <- resK$eps
    }
    
    K_best <- K_res[which.max(K_res$epsilon), ]
    print(K_best)
 
  
  }
  
  ## 1.4 Optimize alpha (wshort, wlong, K connus) ----
  if(TRUE){
    # Description :
    # Optimise the climate-variable weights (alpha) once the temporal windows
    # (wshort, wlong) and the number of neighbours K are fixed.
    #
    # The goal is to find separate weights for the short and long windows
    # for the four climate variables (TN, TX, SNOW, NDVI), i.e. 8 parameters in total.
    #
    # The function f.ALPHA() builds a weighted climate-distance matrix from DMfull,
    # runs a leave-one-out KNN reconstruction using the fixed K,
    # and computes a global performance score epsilon = OA + Kappa.
    #
    # genoud searches the alpha vector (integers 0–100, constrained to sum to 100)
    # that minimizes the objective 2 - epsilon.
    #
    # The optimal alpha values are saved together with wshort, wlong and K
    # into a single CSV file for the reconstruction step.
    # Fixed parameters from the previous steps
    
    wshort <- wshort_best
    wlong  <- wlong_best
    K      <- K_best$K
    
    # Climate variables used in the alpha optimisation (TN, TX, SNOW, NDVI)
    vars_idx <- 1:4
    nvars    <- 8  
    
    # Distance components restricted to the selected short and long windows
    myDMshort <- DMfull[1:wshort]
    myDMlong  <- DMfull[1:wlong]
    
    # Homogeneous starting alpha vector for rgenoud (integers summing to 100)
    START <- rep(floor(100 / nvars), nvars)
    reste <- 100 - sum(START)
    if (reste > 0) START[1:reste] <- START[1:reste] + 1
    
    # Local parallel cluster 
    ncores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makePSOCKcluster(ncores)
    
    # Make sure all required objects/functions are available on the workers
    parallel::clusterExport(
      cl,
      varlist = c(
        "g", "f.ALPHA", "compute_confusion",
        "myDMshort", "myDMlong",
        "wshort", "wlong", "K",
        "SCA_mat", "n_dates", "vars_idx", "nvars"
      ),
      envir = environment()
    )
    
    # Reproducible random streams (master + workers)
    set.seed(123)
    parallel::clusterSetRNGStream(cl, 123)
    
    # Genetic optimisation of alpha (constraint handled inside g)
    OPT <- genoud(
      fn = g,                                                       # constrained objective: checks sum(x)=100 + returns cost
      nvars = nvars,                                                
      max = FALSE,                                                  # we MINIMIZE the objective
      pop.size = 2000,                                              
      max.generations = 15,                                         
      lexical = 2,                                                  # standard robust settings
      starting.values = START,                                      # valid starting vector (integers, sum=100)
      Domains = matrix(c(rep(0, nvars), rep(100, nvars)), nvars, 2),# bounds for each alpha
      boundary.enforcement = 2,                                     # strict handling of bounds
      data.type.int = TRUE,                                         # optimise integer alphas (0..100)
      cluster = cl                                                  
    )
    
    parallel::stopCluster(cl)
    
    # Best alpha solution returned by genoud
    alpha_opt_raw <- OPT$par
    alpha_opt     <- alpha_opt_raw / 100
    
    # Name the alpha weights clearly
    names(alpha_opt) <- c("alpha_TN_short", "alpha_TX_short", "alpha_SNOW_short",
                          "alpha_NDVI_short", "alpha_TN_long", "alpha_TX_long",
                          "alpha_SNOW_long", "alpha_NDVI_long")
    
    alpha_opt
    OPT$value  
    
    # Save parameters and optimised alphas for the reconstruction step
    df_equation <- data.frame(
      wshort = wshort_best,
      wlong  = wlong_best,
      K      = K,
      t(alpha_opt),
      check.names = FALSE
    )
    write.csv(df_equation, parameters_equation, row.names = FALSE)
 
    
    
   }
}

#### 2. Reconstruction of Query Date ####
#---------------------------------------#
if (TRUE){
  # Description :
  # This block reconstructs daily snow cover within one selected water year.
  # It uses the optimal parameters estimated in Part 1 (wshort, wlong, K and alpha weights).
  #
  # A daily calendar is built for the water year, then split into:
  # - observed dates (existing Landsat/Sentinel acquisitions),
  # - query dates (all missing days).
  #
  # Learning dates are NOT restricted to the water year:
  # the reconstruction always searches neighbours within the full multi-year
  # learning pool (the ~137 viable satellite images).
  #
  # For each query day, the code computes the alpha-weighted climate distance
  # to all learning dates using the fixed short/long windows, selects the K
  # closest neighbours, and applies a majority vote to generate a binary map.
  #
  # The final output is a single raster stack ordered by day, containing both
  # observed layers (with an "_obs" tag) and reconstructed layers (with a "_knn" tag),
  # plus a metadata table to track the source of each layer for plotting.
  
  
  ## Packages ----
  library(terra)
  library(zoo)
  library(dplyr)
  
  source(file.path(functions_dir , "Functions_KNN.R"))
  
  ## PARAMETERS ----
  site <- site
  WATER_YEAR <- 2022  
  
  WY_LABEL   <- sprintf("%d-%d", WATER_YEAR, WATER_YEAR + 1)
  START_DATE <- as.Date(sprintf("%d-10-01", WATER_YEAR))
  END_DATE   <- as.Date(sprintf("%d-09-30", WATER_YEAR + 1))
  
  date_tag <- paste0(format(START_DATE, "%Y%m%d"), "_", format(END_DATE, "%Y%m%d"))
  
  ## INPUT ----
  KNN_case <- file.path(inputs_dir, "KNN_Training_Data")
  
  SCA_file     <- file.path(KNN_case, "SCA13to23.L1s.tif")
  fSCA_file    <- file.path(KNN_case, "THE.fSCA.csv")
  Fluxalp_file <- file.path(KNN_case, "FLUXALP.daily.csv")
  
  gen_case      <- file.path(output_dir, "3. KNN_Snow_Generation")
  equation_case <- file.path(gen_case, "equation_parameters")
  parameters_equation <- file.path(equation_case, paste0("Equation_parameters_", site, ".csv"))
  
  ## OUTPUT ----
  output_case <- file.path(gen_case, paste0("WY", WY_LABEL))
  if (!dir.exists(output_case)) dir.create(output_case, recursive = TRUE)
  
  output_stack_file <- file.path(
    output_case,
    paste0("SCA_KNN_", site, "_WY", WY_LABEL, "_", date_tag, ".tif")
  )
  
  output_meta_file <- file.path(
    output_case,
    paste0("SCA_KNN_", site, "_WY", WY_LABEL, "_", date_tag, "_meta.csv")
  )
  
  
  ## CODE ----
  # Read parameters found in part 1
  params <- read.csv(parameters_equation)
  
  wshort <- params$wshort[1]
  wlong  <- params$wlong[1]
  K <- params$K[1]
  
  vars     <- c("TN","TX","SNOW","NDVI")
  vars_idx <- 1:4
  
  alpha_short <- as.numeric(params[1, c("alpha_TN_short","alpha_TX_short",
                                        "alpha_SNOW_short","alpha_NDVI_short")])
  alpha_long  <- as.numeric(params[1, c("alpha_TN_long","alpha_TX_long",
                                        "alpha_SNOW_long","alpha_NDVI_long")])
  
  
  ## Load Learning data (same objects as part 1)
  img <- terra::rast(SCA_file)
  SCA_mat <- terra::values(img, mat = TRUE)
  n_dates <- ncol(SCA_mat)
  
  df <- read.csv(fSCA_file)
  df <- df[df$WS == "TOT", ]
  df$DATE <- as.Date(paste(df$YEAR, df$DOY, sep = " "), format = "%Y %j")
  
  FXA <- read.csv(Fluxalp_file)
  FXA$DATE <- as.Date(FXA$DATE)
  for (i in 2:6) FXA[, i] <- zoo::na.approx(FXA[, i])
  
  df1 <- dplyr::left_join(df, FXA, by = "DATE")
  stopifnot(nrow(df1) == n_dates)
  
  learning_dates <- df1$DATE
  
  
  ## Build Water Year calendar 
  all_dates <- seq(START_DATE, END_DATE, by = "day")
  
  # Keep only dates covered by the meteorological series
  all_dates <- all_dates[all_dates %in% FXA$DATE]
  
  # Identify observed vs query dates
  obs_mask    <- all_dates %in% learning_dates
  obs_dates   <- all_dates[obs_mask]
  query_dates <- all_dates[!obs_mask]
  
  
  ## Build climate cube for Learning (for MIN/MAX) 
  TMP_L <- list()
  for (i in 0:90){
    tmp_date    <- learning_dates - i
    TMP_L[[i+1]] <- FXA[match(tmp_date, FXA$DATE), vars]
  }
  
  # Same MIN/MAX logic as in part 1
  tmp <- do.call(rbind.data.frame, lapply(TMP_L, f.MIN))
  MIN <- t(apply(tmp, 2, min, na.rm = TRUE))
  colnames(MIN) <- vars
  
  tmp <- do.call(rbind.data.frame, lapply(TMP_L, f.MAX))
  MAX <- t(apply(tmp, 2, max, na.rm = TRUE))
  colnames(MAX) <- vars
  
  # Scale Learning cube using the existing f.SCALE
  INPUT_L <- lapply(TMP_L, f.SCALE)
  names(INPUT_L) <- paste0("DAY_minus", 0:90)
  
  
  ## Build climate cube for Query 
  TMP_Q <- list()
  for (i in 0:90){
    tmp_date    <- query_dates - i
    TMP_Q[[i+1]] <- FXA[match(tmp_date, FXA$DATE), vars]
  }
  
  INPUT_Q <- lapply(TMP_Q, f.SCALE)
  names(INPUT_Q) <- paste0("DAY_minus", 0:90)
  
  # Short window
  Qlag <- Reduce("+", lapply(1:wshort, function(lag) {
    as.matrix(INPUT_Q[[lag]][, vars_idx, drop = FALSE])
  })) / wshort
  
  Llag <- Reduce("+", lapply(1:wshort, function(lag) {
    as.matrix(INPUT_L[[lag]][, vars_idx, drop = FALSE])
  })) / wshort
  
  D_short <- dist_one_lag(Qlag, Llag, alpha_short)
  
  # Long window
  Qlag <- Reduce("+", lapply(1:wlong, function(lag) {
    as.matrix(INPUT_Q[[lag]][, vars_idx, drop = FALSE])
  })) / wlong
  
  Llag <- Reduce("+", lapply(1:wlong, function(lag) {
    as.matrix(INPUT_L[[lag]][, vars_idx, drop = FALSE])
  })) / wlong
  
  D_long <- dist_one_lag(Qlag, Llag, alpha_long)
  
  # Final Query -> Learning distance matrix
  D_QL <- D_short + D_long   # n_query x n_learning
  
  ## Build final WY stack (observed + reconstructed) ----
  template <- img[[1]]
  
  layers      <- vector("list", length(all_dates))
  layer_names <- character(length(all_dates))
  
  q_counter <- 0
  
  for (k in seq_along(all_dates)) {
    
    d <- all_dates[k]
    
    if (d %in% learning_dates) {
      
      idx <- match(d, learning_dates)
      layers[[k]] <- img[[idx]]
      layer_names[k] <- paste0(format(d, "%Y%m%d"), "_obs")
      
    } else {
      
      q_counter <- q_counter + 1
      
      di <- D_QL[q_counter, ]
      nn <- order(di)[1:K]
      
      neigh <- SCA_mat[, nn, drop = FALSE]
      pred_prob <- rowMeans(neigh, na.rm = TRUE)
      pred <- ifelse(pred_prob >= 0.5, 1, 0)
      
      layers[[k]] <- terra::setValues(template, pred)
      layer_names[k] <- paste0(format(d, "%Y%m%d"), "_knn")
    }
  }
  
  # Assemble stack
  stack_r <- layers[[1]]
  if (length(layers) > 1){
    for (k in 2:length(layers)) stack_r <- c(stack_r, layers[[k]])
  }
  
  names(stack_r) <- layer_names
  
  terra::writeRaster(stack_r, output_stack_file, overwrite = TRUE)
  
  
  ## Save metadata marker for plots ----
  df_meta <- data.frame(
    DATE = all_dates,
    layer = layer_names,
    source = ifelse(obs_mask, "observed", "reconstructed"),
    WATER_YEAR = WY_LABEL
  )
  
  write.csv(df_meta, output_meta_file, row.names = FALSE)
  
  cat("WY stack written to:\n", output_stack_file, "\n")
  cat("Metadata written to:\n", output_meta_file, "\n")
  
}



























#### 3. Génération of dayly snow cover for visualisation ####
#-----------------------------------------------------------#
if (TRUE){
  
  ## Packages ----
  library(terra)
  
  ## PARAMETERS ----
  site <- site
  WATER_YEAR <- 2022
  
  WY_LABEL   <- sprintf("%d-%d", WATER_YEAR, WATER_YEAR + 1)
  START_DATE <- as.Date(sprintf("%d-10-01", WATER_YEAR))
  END_DATE   <- as.Date(sprintf("%d-09-30", WATER_YEAR + 1))
  
  date_tag <- paste0(format(START_DATE, "%Y%m%d"), "_", format(END_DATE, "%Y%m%d"))
  
  # Visual window inside the water year
  START_VIS <- as.Date(sprintf("%d-02-01", WATER_YEAR + 1))
  END_VIS   <- as.Date(sprintf("%d-07-30", WATER_YEAR + 1))
  
  ## INPUT ----
  gen_case   <- file.path(output_dir, "3. KNN_Snow_Generation")
  input_case <- file.path(gen_case, paste0("WY", WY_LABEL))
  
  input_stack_file <- file.path(
    input_case,
    paste0("SCA_KNN_", site, "_WY", WY_LABEL, "_", date_tag, ".tif")
  )
  
  input_meta_file <- file.path(
    input_case,
    paste0("SCA_KNN_", site, "_WY", WY_LABEL, "_", date_tag, "_meta.csv")
  )
  
  ## OUTPUT ----
  output_pdf_file <- file.path(
    input_case,
    paste0("SCA_daily_", site, "_WY", WY_LABEL, "_",
           format(START_VIS, "%Y%m%d"), "_", format(END_VIS, "%Y%m%d"), ".pdf")
  )
  
  ## ---- Read data ----
  stack_r <- terra::rast(input_stack_file)
  
  if (file.exists(input_meta_file)) {
    meta <- read.csv(input_meta_file)
    meta$DATE <- as.Date(meta$DATE)
  } else {
    # fallback if meta missing
    layer_names <- names(stack_r)
    dates_txt <- sub("(_obs|_knn)$", "", layer_names)
    meta <- data.frame(
      DATE = as.Date(dates_txt, format = "%Y%m%d"),
      layer = layer_names,
      source = ifelse(grepl("_obs$", layer_names), "observed", "reconstructed"),
      WATER_YEAR = WY_LABEL
    )
  }
  
  ## ---- Filter dates ----
  meta_vis <- meta[meta$DATE >= START_VIS & meta$DATE <= END_VIS, ]
  if (nrow(meta_vis) == 0) stop("No layers found in the requested date window.")
  
  meta_vis <- meta_vis[order(meta_vis$DATE), ]
  
  ## ---- Helpers ----
  snow_fraction_pct <- function(r){
    m <- terra::global(r, "mean", na.rm = TRUE)[1,1]
    as.numeric(m) * 100
  }
  
  draw_binary_map <- function(r, is_obs, date_txt){
    
    # clean palette
    col_nosnow <- "#F5F6F7"
    col_snow   <- "#2C7FB8"
    
    terra::plot(r,
                col = c(col_nosnow, col_snow),
                legend = FALSE,
                axes = FALSE,
                mar = 0)
    
    # soft thin border + highlight if observed
    box(col = if (is_obs) "#D7191C" else "#444444", lwd = if (is_obs) 3 else 1.5)
    
    # title inside left panel
    mtext(date_txt, side = 3, line = -1.2, adj = 0.02, cex = 1.2, font = 2)
    
    # small clean legend
    legend("topright",
           legend = c("No snow", "Snow"),
           fill = c(col_nosnow, col_snow),
           border = NA,
           bty = "n",
           cex = 0.95)
  }
  
  draw_fraction_bar <- function(frac){
    
    # canvas
    plot(NA, xlim = c(0,100), ylim = c(0,1),
         xlab = "Snow-covered area (%)",
         ylab = "",
         yaxt = "n",
         bty = "n",
         cex.lab = 1.05,
         cex.axis = 1.0)
    
    # background bar
    rect(0, 0.25, 100, 0.75, col = "#EEF1F4", border = "#D0D6DC", lwd = 1)
    
    # filled bar
    rect(0, 0.25, frac, 0.75, col = "#2C7FB8", border = NA)
    
    # ticks
    axis(1, at = seq(0,100,20))
    
    # label value
    text(frac, 0.5, labels = sprintf("%.1f%%", frac),
         pos = 4, cex = 1.2, font = 2, col = "#1B1B1B")
    
    title("Daily snow fraction", cex.main = 1.25, font.main = 2)
  }
  
  ## ---- PDF ----
  pdf(output_pdf_file, width = 11.5, height = 6)
  
  for (k in seq_len(nrow(meta_vis))) {
    
    lay_name <- meta_vis$layer[k]
    d        <- meta_vis$DATE[k]
    src      <- meta_vis$source[k]
    is_obs   <- src == "observed"
    
    r <- stack_r[[lay_name]]
    frac <- snow_fraction_pct(r)
    
    # stable layout
    par(
      mfrow = c(1,2),
      mar = c(3, 3, 3, 2),
      oma = c(0, 0, 2.5, 0),
      bg = "white"
    )
    
    # LEFT : map
    draw_binary_map(r, is_obs, format(d, "%Y-%m-%d"))
    
    # RIGHT : fraction bar
    draw_fraction_bar(frac)
    
    # global header
    mtext(
      if (is_obs) "OBSERVED (satellite)" else "RECONSTRUCTED (KNN)",
      side = 3, outer = TRUE, line = 0.3,
      cex = 1.1, font = 2,
      col = if (is_obs) "#D7191C" else "#2C7FB8"
    )
  }
  
  dev.off()
  
  cat("PDF written to:\n", output_pdf_file, "\n")
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
























