rm(list = ls())

source("scripts/00_setup.R")

################################################################################
## SCRIPT 01: Load raw methylation data
################################################################################

## Goal:
## - Inspect raw data folder
## - Load the SampleSheet
## - Create RGset using minfi
## - Save RGset for the next scripts

################################################################################
## 1. Inspect raw data folder
################################################################################

list.files(RAW_DATA_DIR)

sample_sheet_file <- file.path(RAW_DATA_DIR, "SampleSheet_Report_II.csv")
SampleSheet <- read.csv(sample_sheet_file, header = TRUE)
SampleSheet


################################################################################
## 2. Load sample sheet with minfi
################################################################################

baseDir <- RAW_DATA_DIR

targets <- read.metharray.sheet(baseDir)
print(targets)
str(targets)

################################################################################
## 3. Create RGChannelSet object
################################################################################

RGset <- read.metharray.exp(targets = targets)

print(RGset)
str(RGset)

################################################################################
## 4. Save RGset
################################################################################

save(RGset, file = file.path(PROCESSED_DATA_DIR, "RGset.RData"))

message("RGset saved to: ", file.path(PROCESSED_DATA_DIR, "RGset.RData"))
