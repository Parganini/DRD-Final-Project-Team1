# DRD Final Project - Team 1
# Setup script

# Install packages if needed:
# BiocManager::install(c("minfi", "IlluminaHumanMethylation450kmanifest"))
# install.packages(c("tidyverse", "ggplot2", "pheatmap", "qqman"))

library(minfi)
library(tidyverse)
library(ggplot2)
library(pheatmap)

# Team 1 parameters
TEAM_ID <- 1
TEAM_ADDRESS <- "42796479"
DETP_THRESHOLD <- 0.05
NORMALIZATION_METHOD <- "preprocessNoob"
DIFFERENTIAL_TEST <- "t-test"

# Paths
RAW_DATA_DIR <- "data/raw"
PROCESSED_DATA_DIR <- "data/processed"
FIGURES_DIR <- "figures"
RESULTS_DIR <- "results"
