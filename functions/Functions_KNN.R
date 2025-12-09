# Computes column-wise minima for a climate lag table.
# Used to build global scaling bounds across all lags/dates.
# Input: data.frame/matrix of climate variables. Output: numeric vector.
f.MIN <- function(x) apply(x,2,min,na.rm=T)

# Computes column-wise maxima for a climate lag table.
# Used to build global scaling bounds across all lags/dates.
# Input: data.frame/matrix of climate variables. Output: numeric vector.
f.MAX <- function(x) apply(x,2,max,na.rm=T)


# Computes Manhattan distances between all pairs of dates for each variable.
# Input: one lag table (n_dates x n_vars). 
# Output: matrix (n_pairs x n_vars), one distance vector per variable.
f.MAN <- function(mydf) {
  # mydf : n_dates x n_vars  (TN, TX, SNOW, NDVI)
  # output : matrix (n_pairs x n_vars), one column per variable
  sapply(seq_len(ncol(mydf)), function(j) {
    dist(mydf[, j], method = "manhattan")
  })
}


# Scales climate variables to [0,1] using global MIN/MAX computed beforehand.
# Input: one lag table with raw values.
# Output: scaled table with the same structure.
f.SCALE <- function(x) {
  tmp1 <- x %>% sweep(2, MIN)
  tmp2 <- tmp1 %>% sweep(2,MAX-MIN,"/")
}


# Helper to verify the distance logic on TX across a chosen lag window.
# Input: two date indices + lag set. Output: scalar TX distance.
# Can return a sum or a mean depending on scale_by_length.
compute_TX_window_dist <- function(i, j, lags = 1:10, scale_by_length = TRUE) {
  v_i <- unlist(lapply(lags, function(l) INPUT[[l]][i, "TX"]))
  v_j <- unlist(lapply(lags, function(l) INPUT[[l]][j, "TX"]))
  
  d <- sum(abs(v_i - v_j))   # somme
  if (scale_by_length) d <- d / length(lags)  # moyenne si TRUE
  
  return(d)
}


# Builds the full climate distance matrix between all learning dates.
# Uses precomputed per-lag/per-variable distances (DMfull).
# Output: n_dates x n_dates matrix combining short and long windows.
build_climate_distance <- function(wshort, wlong, vars_idx = 1:3) {
  
  # --- Short window ---
  d_short_vec <- Reduce("+", lapply(seq_len(wshort), function(lag) {
    rowSums(DMfull[[lag]][, vars_idx, drop = FALSE])
  })) / wshort
  
  D_short <- as.matrix(structure(d_short_vec,
                                 class = "dist",
                                 Size  = n_dates))
  # --- Long window ---
  d_long_vec <- Reduce("+", lapply(seq_len(wlong), function(lag) {
    rowSums(DMfull[[lag]][, vars_idx, drop = FALSE])
  })) / wlong
  
  D_long <- as.matrix(structure(d_long_vec,
                                class = "dist",
                                Size  = n_dates))
  # --- Final matrice ---
  D_short + D_long
}


# Computes confusion counts between observed and predicted binary maps.
# Removes NA pairs before counting.
# Output: TP/TN/FP/FN for one comparison.
compute_confusion <- function(obs, pred) {
  ok   <- !is.na(obs) & !is.na(pred)
  obs  <- obs[ok]
  pred <- pred[ok]
  
  TP <- sum(obs == 1 & pred == 1)
  TN <- sum(obs == 0 & pred == 0)
  FP <- sum(obs == 0 & pred == 1)
  FN <- sum(obs == 1 & pred == 0)
  
  list(TP = TP, TN = TN, FP = FP, FN = FN)
}


# Evaluates one (wshort, wlong) pair with LOOCV and K=1 behaviour.
# For each target date, selects the single closest climate neighbour.
# Returns global OA, Kappa and epsilon = OA + Kappa.
eval_w <- function(wshort, wlong, vars_idx = 1:3) {
  D <- build_climate_distance(wshort, wlong, vars_idx = vars_idx)
  
  allTP <- allTN <- allFP <- allFN <- 0
  
  for (i in seq_len(n_dates)) {
    
    idx_train <- setdiff(seq_len(n_dates), i)
    
    # K = 1 → closest neighbour only
    di <- D[i, idx_train]
    nn <- idx_train[which.min(di)]
    
    pred <- SCA_mat[, nn]     # predicted snow map
    obs  <- SCA_mat[, i]      # true snow map
    
    m <- compute_confusion(obs, pred)
    allTP <- allTP + m$TP
    allTN <- allTN + m$TN
    allFP <- allFP + m$FP
    allFN <- allFN + m$FN
  }
  
  N  <- allTP + allTN + allFP + allFN
  OA <- (allTP + allTN) / N
  
  pe <- ((allTP + allFP) * (allTP + allFN) +
           (allFN + allTN) * (allFP + allTN)) / (N^2)
  kappa <- if (pe != 1) (OA - pe) / (1 - pe) else NA_real_
  
  list(OA = OA, kappa = kappa, eps = OA + kappa)
}


# Evaluates a candidate number of neighbours K with LOOCV.
# Uses a fixed climate distance matrix and majority vote on K maps.
# Returns global OA, Kappa and epsilon for this K.
eval_K <- function(K, D, SCA_mat) {
  cat("   eval K =", K, "\n")
  D <- as.matrix(D)
  n_dates <- ncol(SCA_mat)
  
  allTP <- allTN <- allFP <- allFN <- 0
  
  for (i in seq_len(n_dates)) {
    idx_train <- setdiff(seq_len(n_dates), i)
    di <- D[i, idx_train]
    nn <- idx_train[order(di)[1:K]]
    
    # KNN : mod
    pred_prob <- rowMeans(SCA_mat[, nn, drop = FALSE])
    pred      <- ifelse(pred_prob >= 0.5, 1, 0)
    
    obs <- SCA_mat[, i]
    # Compute confusion on one date
    m <- compute_confusion(obs, pred)
    
    allTP <- allTP + m$TP
    allTN <- allTN + m$TN
    allFP <- allFP + m$FP
    allFN <- allFN + m$FN
  }
  
  N  <- allTP + allTN + allFP + allFN
  OA <- (allTP + allTN) / N
  
  pe <- ((allTP + allFP) * (allTP + allFN) +
           (allFN + allTN) * (allFP + allTN)) / (N^2)
  
  kappa <- if (pe != 1) (OA - pe) / (1 - pe) else NA_real_
  
  list(OA = OA, kappa = kappa, eps = OA + kappa)
}


# Objective function for alpha optimisation with fixed wshort/wlong/K.
# Builds a weighted climate distance matrix and runs LOOCV KNN.
# Returns the cost to minimize: 2 - (OA + Kappa).
f.ALPHA <- function(x) {
  
  xrshort <- x[1:4] / 100
  xrlong  <- x[5:8] / 100
  
  SHshort <- lapply(myDMshort, function(X) {
    Xuse <- as.matrix(X[, vars_idx, drop = FALSE])
    sweep(Xuse, 2, xrshort, "*")
  })
  
  SHlong <- lapply(myDMlong, function(X) {
    Xuse <- as.matrix(X[, vars_idx, drop = FALSE])
    sweep(Xuse, 2, xrlong, "*")
  })
  
  # moyenne sur les fenêtres
  SHORT <- Reduce("+", SHshort) / wshort
  LONG  <- Reduce("+", SHlong)  / wlong
  
  SHmean <- rowSums(SHORT) + rowSums(LONG)
  
  DISTmat <- as.matrix(structure(SHmean, class = "dist", Size = n_dates))
  
  allTP <- allTN <- allFP <- allFN <- 0
  
  for (i in seq_len(n_dates)) {
    idx_train <- setdiff(seq_len(n_dates), i)
    di  <- DISTmat[i, idx_train]
    nn  <- idx_train[order(di)[1:K]]
    
    neigh <- SCA_mat[, nn, drop = FALSE]
    pred  <- ifelse(rowMeans(neigh) >= 0.5, 1, 0)
    obs   <- SCA_mat[, i]
    
    m <- compute_confusion(obs, pred)
    allTP <- allTP + m$TP
    allTN <- allTN + m$TN
    allFP <- allFP + m$FP
    allFN <- allFN + m$FN
  }
  
  N  <- allTP + allTN + allFP + allFN
  OA <- (allTP + allTN) / N
  
  pe <- ((allTP + allFP) * (allTP + allFN) +
           (allFN + allTN) * (allFP + allTN)) / (N^2)
  
  kappa <- if (pe != 1) (OA - pe) / (1 - pe) else 0
  
  eps <- OA + kappa
  eps <- max(min(eps, 2), 0)
  
  2 - eps
}


# Constraint wrapper for genoud.
# Enforces sum(x) = 100 for alpha candidates.
# Returns a two-part vector: constraint status + objective value.
g <- function(x) {
  c(
    ifelse(sum(x) == 100, 0, 1),
    f.ALPHA(x)
  )
}


## Distance Query -> Learning (alpha weighted) 
dist_one_lag <- function(Qlag, Llag, alpha_vec){
  D <- 0
  for (v in 1:4){
    D <- D + alpha_vec[v] * abs(outer(Qlag[, v], Llag[, v], "-"))
  }
  D
}

