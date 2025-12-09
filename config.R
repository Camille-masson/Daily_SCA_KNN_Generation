## CONFIGURATION DU PROJET ###

# Définition des chemins 
root_dir <- getwd()  # Récupère le chemin du projet automatiquement
downloads_dir <- file.path(root_dir, "downloads")
output_dir <- file.path(root_dir, "outputs")
functions_dir <- file.path(root_dir, "functions")
inputs_dir <- file.path(root_dir, "inputs")

if (! dir.exists(downloads_dir)) dir.create(downloads_dir)
if (! dir.exists(output_dir)) dir.create(output_dir)
if (! dir.exists(functions_dir)) dir.create(functions_dir)
if (! dir.exists(inputs_dir)) dir.create(inputs_dir)

# --- renv ignore big data folders ---
renvignore_path <- file.path(getwd(), ".renvignore")
if (!file.exists(renvignore_path)) {
  writeLines(c(
    "downloads/",
    "outputs/",
    "inputs/"
  ), renvignore_path)
}

#Installed.Package
if (!require("renv")) install.packages("renv")
library(renv)
deps <- renv::dependencies()
packages <- unique(deps$Package)

installed_pkgs <- installed.packages()[, "Package"]
missing_pkgs <- setdiff(packages, installed_pkgs)

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE, repos = "https://cran.rstudio.com/")
}

#sapply(packages, require, character.only = TRUE)




