#### 0. LIBRARIES AND CONSTANTS ####
#----------------------------------#

gc()
# Loading configuration
source("config.R")



#### 1. PREP DATA  ####
#---------------------#

if(TRUE){
  
  # Description:
  # This section prepares the input and output paths for one study area
  # and one water year, then builds a long-format dataset for snow
  # modelling by combining the Sentinel-2 snow stack with topographic
  # variables (Lidar) and saving the result as an .rds file in : outputs/2. Dataset_Snow_Modeling
  
  
  
  ## PARAMETERS ----
  STUDY_AREA <- "Cayolle" # <- Replace with the name of the study area
  WATER_YEAR <- 2021 # WARNING: the water year is defined here from 1 Oct of WATER_YEAR to 30 Sep of WATER_YEAR + 1.
  
  WY_LABEL   <- sprintf("%d-%d", WATER_YEAR, WATER_YEAR + 1)
  START_DATE <- sprintf("%d-10-01T00:00:00Z", WATER_YEAR)
  END_DATE   <- sprintf("%d-09-30T23:59:59Z", WATER_YEAR + 1)
  
  
  ## INPUT ----
  # Input folder for the processed Snow Cover product
  input_case <- file.path(output_dir, "1. Snow_cover_by_site")
  
  date_tag <- paste0(gsub("[^0-9]", "", substr(START_DATE, 1, 10)), "_",
                     gsub("[^0-9]", "", substr(END_DATE,   1, 10)))
  
  input_snow_image_file <- file.path(input_case,paste0("Theia_Sentinel-2_Snow_Cover_Level-2B_",STUDY_AREA,"_WY", WY_LABEL, "_", date_tag,".tif"))
  
  # Create folder to store the topographic variables.
  topo_case <- file.path(inputs_dir, "Topography")
  if (!dir.exists(topo_case)) {
    dir.create(topo_case, recursive = TRUE)
  }
  topo_file <- file.path (topo_case, paste0("TOPO_1_",STUDY_AREA,".tif"))
  
  
  ## OUTPUT ----
  # Create output folder for the processed Snow Cover product
  output_case <- file.path(output_dir, "2. Dataset_Snow_Modeling")
  if (!dir.exists(output_case)) {
    dir.create(output_case, recursive = TRUE)
  }
  
  output_file <- file.path(output_case,paste0("Dataset_Snow_Modeling_",STUDY_AREA,
                                              "_WY", WY_LABEL, "_", date_tag,".rds"))
  
  ## CODE ----
  build_dataset_snow_modeling(input_snow_image_file,topo_file,output_file)
  
  
}


#### 2. Cloud-gap filling  ####
#-----------------------------#

if(TRUE){
  
  # Description:


  ## PARAMETERS ----
  STUDY_AREA <- "Cayolle" # <- Replace with the name of the study area
  WATER_YEAR <- 2021 # WARNING: the water year is defined here from 1 Oct of WATER_YEAR to 30 Sep of WATER_YEAR + 1.
  
  WY_LABEL   <- sprintf("%d-%d", WATER_YEAR, WATER_YEAR + 1)
  START_DATE <- sprintf("%d-10-01T00:00:00Z", WATER_YEAR)
  END_DATE   <- sprintf("%d-09-30T23:59:59Z", WATER_YEAR + 1)
  
  
  ## INPUT ----
  # Input folder for the processed Snow Cover product
  input_case <- file.path(output_dir, "2. Dataset_Snow_Modeling")
  
  date_tag <- paste0(gsub("[^0-9]", "", substr(START_DATE, 1, 10)), "_",
                     gsub("[^0-9]", "", substr(END_DATE,   1, 10)))
  
  snow_dataset_file <- file.path(input_case,paste0("Dataset_Snow_Modeling_",STUDY_AREA,
                                              "_WY", WY_LABEL, "_", date_tag,".rds"))
  
  
  
  ## CODE ----
  
  snow <- readRDS(snow_dataset_file)
  
  summary(snow)
  
  library(dplyr)
  library(ggplot2)
  
  gap_threshold <- 5L
  
  ## 1. Unique acquisition dates
  acq_dates <- snow |>
    distinct(date) |>
    arrange(date) |>
    mutate(
      y     = 0,
      label = format(date, "%Y-%m-%d")
    )
  
  ## 2. Gaps between acquisitions
  date_gaps <- acq_dates |>
    mutate(
      prev_date  = lag(date),
      delta_days = as.integer(date - prev_date)
    ) |>
    filter(!is.na(delta_days) & delta_days > gap_threshold) |>
    mutate(
      mid_date = prev_date + (date - prev_date) / 2,
      gap_lab  = paste0(delta_days, " d")
    )
  
  ## 3. Plot
ggplot() +
  # red bands = gaps > threshold, full panel height
  geom_rect(
    data = date_gaps,
    aes(xmin = prev_date, xmax = date),
    ymin = -Inf, ymax = Inf,
    fill = "red",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  # points = acquisition dates
  geom_point(
    data = acq_dates,
    aes(x = date, y = y),
    size = 2.5
  ) +
  # date labels, larger, diagonal, just below the point
  geom_text(
    data = acq_dates,
    aes(x = date, y = y, label = label),
    angle = 55,
    vjust = 0.5,
    hjust = 1.1,
    size  = 4
  ) +
  # labels for gap length (lower inside the red band)
  geom_text(
    data = date_gaps,
    aes(x = mid_date, y = 0.07, label = gap_lab),  # <- plus bas qu'avant
    size = 3.5,
    vjust = 0.5
  ) +
  scale_y_continuous(
    NULL,
    breaks = NULL,
    limits = c(-0.05, 0.12),  # <- on rogne le haut du graph
    expand = c(0, 0)
  ) +
  scale_x_date(
    name        = "Acquisition date",
    date_breaks = "1 month",
    date_labels = "%Y-%m-%d"
  ) +
  labs(
    title    = "Sentinel-2 snow cover acquisitions for Water year 2021/2022",
    subtitle = paste0(
      "Red bands = gaps > ", gap_threshold,
      " days between acquisitions"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1)
  )
  
  
  
library(dplyr)
library(lubridate)
library(ranger)

## ------------------------------------------------------------------
## 0. Préparation : date au bon format
## ------------------------------------------------------------------

snow <- snow %>%
  mutate(date = as.Date(date))   # au cas où

## ------------------------------------------------------------------
## 1. Règle T-5 / T+5 : comblement simple des NA
## ------------------------------------------------------------------

# Décaler les valeurs de snow de +/- 5 jours pour chaque pixel
snow_minus5 <- snow %>%
  transmute(
    cell,
    date = date + 5,           # ce qui était à T-5 se recale sur T
    snow_minus5 = snow
  )

snow_plus5 <- snow %>%
  transmute(
    cell,
    date = date - 5,           # ce qui était à T+5 se recale sur T
    snow_plus5 = snow
  )

# Jointure + application de la règle
snow_filled <- snow %>%
  left_join(snow_minus5, by = c("cell", "date")) %>%
  left_join(snow_plus5,  by = c("cell", "date")) %>%
  mutate(
    can_fill = is.na(snow) &
      !is.na(snow_minus5) &
      !is.na(snow_plus5) &
      (snow_minus5 == snow_plus5),
    
    snow_filled = if_else(can_fill, snow_minus5, snow)
  )

# Stats
n_na_before <- sum(is.na(snow$snow))
n_filled    <- sum(snow_filled$can_fill)
n_na_after  <- sum(is.na(snow_filled$snow_filled))

cat("Règle T-5/T+5\n")
cat("  NA avant   :", n_na_before, "\n")
cat("  NA comblés :", n_filled,    "\n")
cat("  NA restants:", n_na_after,  "\n\n")

# On écrase la colonne snow et on nettoie les colonnes intermédiaires
snow <- snow_filled %>%
  mutate(snow = snow_filled) %>%
  select(-snow_filled, -snow_minus5, -snow_plus5, -can_fill)


  

library(dplyr)
library(lubridate)

library(ranger)

## ------------------------------------------------------------------
## 2. Construction des features pour la RF
## ------------------------------------------------------------------

snow_feat <- snow %>%
  mutate(
    # variables temporelles
    doy = yday(date),
    doy_sin = sin(2 * pi * doy / 365),
    doy_cos = cos(2 * pi * doy / 365),
    
    # aspect en radians + sin/cos (patterns spatiaux liés au soleil/vent)
    asp_rad  = aspect_deg * pi / 180,
    asp_sin  = sin(asp_rad),
    asp_cos  = cos(asp_rad),
    
    # coord normalisées (pour capturer un gradient spatial large)
    x_centered = scale(x)[, 1],
    y_centered = scale(y)[, 1]
  )

# On garde une version de travail
snow_rf <- snow_feat

## ------------------------------------------------------------------
## 3. Préparation des données d'apprentissage (pixels clairs)
## ------------------------------------------------------------------

# Lignes avec neige connue (0/1)
obs <- snow_rf %>%
  filter(!is.na(snow))

n_obs <- nrow(obs)
cat("Nombre total d'observations claires (0/1) pour l'entraînement :", n_obs, "\n")

# Pour éviter d'avoir 10x plus de 0 que de 1 : échantillonnage équilibré
set.seed(123)

obs_1 <- obs %>% filter(snow == 1)
obs_0 <- obs %>% filter(snow == 0)

n1 <- nrow(obs_1)
n0 <- nrow(obs_0)

cat("  1 (neige) :", n1, "pixels\n")
cat("  0 (pas neige) :", n0, "pixels\n")

# On équilibre : on prend min(n1, n0) de chaque
n_train_each <- min(n1, n0, 200000)  # max 200k par classe pour rester raisonnable

train_1 <- obs_1 %>% sample_n(n_train_each)
train_0 <- obs_0 %>% sample_n(n_train_each)

train_rf <- bind_rows(train_0, train_1)

cat("Taille du jeu d'entraînement équilibré :", nrow(train_rf), "\n")

# Réponse en facteur pour classification
train_rf <- train_rf %>%
  mutate(snow_factor = factor(snow, levels = c(0, 1)))

# Choix des prédicteurs topo + temps + espace
predictors <- c(
  "elev",
  "slope_deg",
  "asp_sin", "asp_cos",
  "tpi", "tri", "roughness", "curvature",
  "dah",
  "doy", "doy_sin", "doy_cos",
  "x_centered", "y_centered"
)

## ------------------------------------------------------------------
## 4. Entraînement de la Random Forest (ranger)
## ------------------------------------------------------------------

set.seed(123)

rf_model <- ranger(
  formula         = snow_factor ~ .,
  data            = train_rf[, c("snow_factor", predictors)],
  num.trees       = 500,
  mtry            = floor(sqrt(length(predictors))),  # classique : sqrt(p)
  min.node.size   = 50,        # un peu de régularisation (évite l'overfit)
  importance      = "impurity",
  probability     = TRUE,      # on veut des proba
  classification  = TRUE
)

cat("\nOOB error rate RF (%) :", rf_model$prediction.error * 100, "\n")

# Importance des variables (optionnel, mais intéressant)
print(sort(rf_model$variable.importance, decreasing = TRUE))

## ------------------------------------------------------------------
## 5. Prédiction sur les pixels NA (gap filling)
## ------------------------------------------------------------------

idx_na <- which(is.na(snow_rf$snow))

cat("\nNombre de pixels NA à prédire :", length(idx_na), "\n")

if (length(idx_na) > 0) {
  pred_rf <- predict(
    rf_model,
    data = snow_rf[idx_na, predictors]
  )
  
  # Matrice de proba (colonnes "0" et "1")
  p1 <- pred_rf$predictions[, "1"]
  
  # Classification avec seuil 0.5 (on pourra optimiser plus tard)
  snow_pred_class <- ifelse(p1 >= 0.5, 1, 0)
  
  # Remplissage des NA
  snow_rf$snow[idx_na] <- snow_pred_class
}

cat("NA restants dans snow après RF :", sum(is.na(snow_rf$snow)), "\n")

# On remplace l'objet snow par la version complétée, en gardant les colonnes utiles
snow <- snow_rf
summary(snow)

library(dplyr)
library(tidyr)
library(terra)

# On part de ton objet 'snow' actuel (déjà gap-fillé, avec la colonne snow)
# Colonnes importantes : x, y, date, snow

snow <- snow %>%
  mutate(date = as.Date(date))

## 1) Passer en large : 1 ligne = 1 pixel (x,y), 1 colonne = 1 date ----

snow_wide <- snow %>%
  select(x, y, date, snow) %>%
  mutate(
    # nom de couche propre (et dans l'ordre chrono)
    date_str = format(date, "%Y%m%d")
  ) %>%
  select(-date) %>%
  tidyr::pivot_wider(
    id_cols   = c(x, y),       # un point par coordonnée
    names_from  = date_str,    # une couche par date
    values_from = snow
  )

# (optionnel) trier par y, x pour être sûr de la cohérence spatiale
snow_wide <- snow_wide %>%
  arrange(y, x)

## 2) Construire le stack raster (SpatRaster terra) ----

# La première colonne = x, la deuxième = y, le reste = couches de neige
r_snow <- terra::rast(
  snow_wide,
  type = "xyz",              # 1: x, 2: y, 3+ : valeurs
  crs  = "EPSG:2154"         # Lambert 93, comme tu veux
)

# Vérifier vite fait
print(r_snow)
# plot(r_snow[[1]])  # pour voir la première date par exemple

## 3) Écrire en GeoTIFF multi-bandes (une bande = une date) ----

writeRaster(
  r_snow,
  filename  = "snow_gapfill_2021_2022_stack.tif",
  overwrite = TRUE
)

# Si tu veux aussi un fichier par jour (optionnel) :
# dates_layers <- names(r_snow)
# for (nm in dates_layers) {
#   writeRaster(r_snow[[nm]],
#               filename = paste0("snow_", nm, ".tif"),
#               overwrite = TRUE)
# }



}





































































