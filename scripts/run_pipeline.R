################################################################################
## Run full Team 1 pipeline
################################################################################

source("scripts/01_load_raw_data.R")
source("scripts/02_red_green_fluorescence.R")
source("scripts/03_mset_raw_and_qc.R")
source("scripts/04_raw_beta_m_values.R")
source("scripts/05_normalization_noob.R")
source("scripts/06_pca_normalized_beta.R")
source("scripts/07_differential_methylation_ttest.R")
source("scripts/08_multiple_testing_and_plots.R")
source("scripts/09_heatmap_top100.R")
