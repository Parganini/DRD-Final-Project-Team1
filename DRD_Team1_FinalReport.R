###############################################################################
### Team 1 methylation array workflow
###############################################################################
##
##
## Team 1 assignments:
## - Step 3 address: 42796479
## - Step 5 detection p-value threshold: 0.05
## - Step 7 normalization method: preprocessNoob
## - Step 9 differential methylation test: t-test
##

# Step 00. Setup

### Clean the workspace before starting the analysis.
rm(list=ls())

### Load all packages used by the complete workflow.
library(minfi)
library(qqman)
library(gplots)

### Define input and output folders relative to the repository root.
baseDir <- "data/raw"
rdsDir <- "results/rds"
tablesDir <- "results/tables"
figuresDir <- "results/figures"

### Create output folders for intermediate objects, tables, and figures.
dir.create(rdsDir, recursive=TRUE, showWarnings=FALSE)
dir.create(tablesDir, recursive=TRUE, showWarnings=FALSE)
dir.create(figuresDir, recursive=TRUE, showWarnings=FALSE)

###############################################################################
### STEP 01: load raw data with minfi
###############################################################################



###############################################################################
### STEP 02: extract Red and Green fluorescence values
###############################################################################



###############################################################################
### STEP 03: fluorescence values for the Team 1 address
###############################################################################
## Team 1 address: 42796479.



###############################################################################
### STEP 04: create the raw MethylSet object
###############################################################################



###############################################################################
### STEP 05: quality checks
###############################################################################
## Team 1 detection p-value threshold: 0.05.



###############################################################################
### STEP 06: calculate raw beta and M values
###############################################################################



###############################################################################
### STEP 07: normalization using preprocessNoob
###############################################################################
## Team 1 normalization method: preprocessNoob.



###############################################################################
### STEP 08: PCA on normalized beta values
###############################################################################



###############################################################################
### STEP 09: differential methylation analysis using t-test
###############################################################################
## Team 1 differential methylation test: t-test.



###############################################################################
### STEP 10: multiple testing correction
###############################################################################



###############################################################################
### STEP 11: volcano and Manhattan plots
###############################################################################



###############################################################################
### STEP 12: heatmap of the top 100 CpG probes
###############################################################################



