# Generates a daily snow visualisation PDF from a multi-layer SCA stack over a 
# chosen date window. Optionally prepares a cached hillshade background at a target
# resolution and overlays binary snow maps with a clean cartographic layout. Computes
# daily snow-covered area (%) and adds a compact vertical gauge alongside each map,
# returning the filtered metadata and outputs invisibly.

create_daily_snow_pdf <- function(
    input_stack_file,
    input_meta_file = NULL,
    hill_file = NULL,
    map_case = ".",
    output_pdf_file,
    START_VIS,
    END_VIS,
    WY_LABEL = NA_character_,
    FAST = FALSE,
    HILL_TARGET_RES = 5
){
  
  # Read data
  stack_r <- terra::rast(input_stack_file)
  
  # Meta 
  if (!is.null(input_meta_file) && file.exists(input_meta_file)) {
    
    meta <- read.csv(input_meta_file)
    meta$DATE <- as.Date(meta$DATE)
    
  } else {
    
    layer_names <- names(stack_r)
    dates_txt   <- sub("(_obs|_knn)$", "", layer_names)
    
    meta <- data.frame(
      DATE       = as.Date(dates_txt, format = "%Y%m%d"),
      layer      = layer_names,
      source     = ifelse(grepl("_obs$", layer_names), "observed", "reconstructed"),
      WATER_YEAR = WY_LABEL
    )
  }
  
  # Filter dates
  meta_vis <- meta[meta$DATE >= START_VIS & meta$DATE <= END_VIS, ]
  if (nrow(meta_vis) == 0) stop("No layers found in the requested date window.")
  meta_vis <- meta_vis[order(meta_vis$DATE), ]
  
  # Hill handling (basemap-like, single prep + cache)
  hill_coarse <- NULL
  
  if (!is.null(hill_file) && file.exists(hill_file)) {
    
    template_r <- stack_r[[ meta_vis$layer[1] ]]
    
    hill_cache <- file.path(
      map_case,
      paste0(
        "Hill_bg_",
        tools::file_path_sans_ext(basename(hill_file)),
        "_", HILL_TARGET_RES, "m.tif"
      )
    )
    
    if (file.exists(hill_cache) && !FAST) {
      
      hill_coarse <- terra::rast(hill_cache)
      
    } else {
      
      hill_raw <- terra::rast(hill_file)
      
      # aggregate FIRST on raw hill resolution
      res_raw <- mean(terra::res(hill_raw))
      fact <- max(1, round(HILL_TARGET_RES / res_raw))
      
      hill_agg <- if (fact > 1) {
        terra::aggregate(hill_raw, fact = fact, fun = "mean", na.rm = TRUE)
      } else {
        hill_raw
      }
      
      if (FAST) {
        hill_agg <- terra::aggregate(hill_agg, fact = 2, fun = "mean", na.rm = TRUE)
      }
      
      # project ONLY CRS (do not force snow grid resolution)
      hill_proj <- terra::project(hill_agg, terra::crs(template_r))
      
      # crop to snow extent
      hill_coarse <- terra::crop(hill_proj, terra::ext(template_r))
      
      if (!FAST) {
        terra::writeRaster(hill_coarse, hill_cache, overwrite = TRUE)
      }
    }
  }
  
  # Helpers
  snow_fraction_pct <- function(r){
    m <- terra::global(r, "mean", na.rm = TRUE)[1,1]
    as.numeric(m) * 100
  }
  
  draw_extent_border <- function(r, col = "#333333", lwd = 1.5){
    e_poly <- terra::as.polygons(terra::ext(r))
    terra::lines(e_poly, col = col, lwd = lwd)
  }
  
  draw_binary_map <- function(r, is_obs, hill = NULL){
    
    col_nosnow <- "#F5F6F7"
    col_snow   <- "#2C7FB8"
    
    if (!is.null(hill)) {
      
      terra::plot(
        hill,
        col    = gray.colors(100, start = 0.15, end = 0.95),
        legend = FALSE,
        axes   = TRUE,
        xlab   = "Longitude",
        ylab   = "Latitude",
        xaxs   = "i",
        yaxs   = "i"
      )
      
      r_snow <- r
      r_snow[r_snow == 0] <- NA
      
      terra::plot(
        r_snow,
        col    = col_snow,
        legend = FALSE,
        axes   = FALSE,
        add    = TRUE
      )
      
    } else {
      
      terra::plot(
        r,
        col    = c(col_nosnow, col_snow),
        legend = FALSE,
        axes   = TRUE,
        xlab   = "Longitude",
        ylab   = "Latitude",
        xaxs   = "i",
        yaxs   = "i"
      )
    }
    
    draw_extent_border(
      r,
      col = if (is_obs) "#D7191C" else "#333333",
      lwd = if (is_obs) 3 else 1.5
    )
    
    par(xpd = NA)
    legend(
      "topright",
      inset  = c(-0.04, 0.02),
      legend = c("No snow", "Snow"),
      fill   = c(col_nosnow, col_snow),
      border = NA,
      bty    = "n",
      cex    = 1.05,
      x.intersp = 0.6,
      y.intersp = 0.9
    )
    par(xpd = FALSE)
  }
  
  draw_fraction_gauge <- function(frac){
    
    col_snow <- "#2C7FB8"
    
    plot(
      NA,
      xlim = c(0, 1),
      ylim = c(0, 100),
      xlab = "",
      ylab = "",
      xaxt = "n",
      yaxt = "n",
      bty  = "n"
    )
    
    mtext("Daily snow fraction", side = 3, line = 0.6, cex = 1.1, font = 2)
    
    axis(
      side = 2,
      at   = seq(0, 100, 20),
      las  = 1,
      pos  = 0.48,
      cex.axis = 0.9
    )
    
    rect(0.58, 0, 0.82, 100, col = "#EEF1F4", border = "#D0D6DC", lwd = 1)
    rect(0.58, 0, 0.82, frac, col = col_snow, border = NA)
    
    text(
      x = 0.85, y = frac,
      labels = sprintf("%.1f%%", frac),
      pos = 4, cex = 1.2, font = 2
    )
  }
  
  # PDF
  pdf(output_pdf_file, width = 12, height = 6.5, pointsize = 11)
  
  for (k in seq_len(nrow(meta_vis))) {
    
    lay_name <- meta_vis$layer[k]
    d        <- meta_vis$DATE[k]
    src      <- meta_vis$source[k]
    is_obs   <- src == "observed"
    
    r    <- stack_r[[lay_name]]
    frac <- snow_fraction_pct(r)
    
    layout(
      matrix(
        c(1, 1,
          2, 3),
        nrow = 2, byrow = TRUE
      ),
      heights = c(0.12, 0.88),
      widths  = c(1.35, 0.65)
    )
    
    ## HEADER
    par(mar = c(0, 0, 0, 0), bg = "white")
    plot.new()
    header_txt <- sprintf(
      "%s — %s",
      if (is_obs) "OBSERVED (satellite)" else "RECONSTRUCTED (KNN)",
      format(d, "%Y-%m-%d")
    )
    text(
      0.5, 0.5,
      header_txt,
      cex = 1.2, font = 2,
      col = if (is_obs) "#D7191C" else "#2C7FB8"
    )
    
    ## MAP
    par(mar = c(2.4, 3.0, 1.2, 3.2), bg = "white")
    draw_binary_map(r, is_obs, hill = hill_coarse)
    
    ## GAUGE
    par(mar = c(2.0, 1.6, 1.2, 1.0), bg = "white")
    draw_fraction_gauge(frac)
  }
  
  dev.off()
  
  invisible(list(
    stack = stack_r,
    meta_vis = meta_vis,
    hill_bg = hill_coarse,
    output_pdf = output_pdf_file
  ))
}



