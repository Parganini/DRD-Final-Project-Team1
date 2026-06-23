## Team 1 assignments:
## - Step 3 address: 42796479
## - Step 5 detection p-value threshold: 0.05
## - Step 7 normalization method: preprocessNoob
## - Step 9 differential methylation test: t-test

###############################################################################
### STEP 00: SETUP
###############################################################################

cat("\n--- Running Step 00: Setup ---\n")

### Clean the workspace before starting the analysis.
rm(list = ls())

### Load all packages used by the complete workflow.
### Install if need it
# BiocManager::install("IlluminaHumanMethylation450kmanifest")
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(minfi)
library(qqman)
library(gplots)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

### Define input and output folders relative to the repository root.
baseDir <- "data/raw"
rdsDir <- "results/rds"
tablesDir <- "results/tables"
figuresDir <- "results/figures"

### Create output folders for intermediate objects, tables, and figures.
dir.create(rdsDir, recursive = TRUE, showWarnings = FALSE)
dir.create(tablesDir, recursive = TRUE, showWarnings = FALSE)
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)

### Report saved outputs in a consistent way when the script is sourced.
report_saved <- function(output_type, path) {
  cat(sprintf("%s saved: %s\n", output_type, path))
}

###############################################################################
### STEP 01: LOAD RAW DATA WITH MINFI
###############################################################################

cat("\n--- Running Step 01: Load Raw Data with minfi ---\n")

### Check the raw data folder and read the sample sheet.
list.files(baseDir)

SampleSheet <- read.csv(file.path(baseDir, "SampleSheet_Report_II.csv"),
                        header = TRUE,
                        stringsAsFactors = FALSE)
head(SampleSheet)

### Read the sample sheet with minfi and load the IDAT files.
targets <- read.metharray.sheet(baseDir)
head(targets)

RGset <- read.metharray.exp(targets = targets)

### Save the RGChannelSet object for the following analysis steps.
save(RGset, file = file.path(rdsDir, "RGset_Report.RData"))
report_saved("RData object", file.path(rdsDir, "RGset_Report.RData"))

### Explore the RGChannelSet object.
RGset
class(RGset)
str(RGset)


###############################################################################
### STEP 02: EXTRACT RED AND GREEN FLUORESCENCE VALUES
###############################################################################

cat("\n--- Running Step 02: Extract Red and Green Fluorescence Values ---\n")

### Extract the Red fluorescence channel from RGset.
Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)

### Extract the Green fluorescence channel from RGset.
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

###############################################################################
### STEP 03: PROBE INFORMATION AND ADDRESS LOOKUP
### Group 1 address: 42796479
###############################################################################

cat("\n--- Running Step 03: Probe Information and Address Lookup ---\n")

TARGET_ADDRESS <- "42796479"

### Check the manifest and extract Type I and Type II probe information.
getManifest(RGset)
ProbeInfo_I  <- data.frame(getProbeInfo(RGset, type = "I"))
head(ProbeInfo_I)

ProbeInfo_II <- data.frame(getProbeInfo(RGset, type = "II"))
head(ProbeInfo_II)

### Look up the Team 1 address in the probe information tables.
ProbeInfo_I[ProbeInfo_I$AddressA == TARGET_ADDRESS, ]
ProbeInfo_I[ProbeInfo_I$AddressB == TARGET_ADDRESS, ]
ProbeInfo_II[ProbeInfo_II$AddressA == TARGET_ADDRESS, ]

# Result: Type I, AddressA → UNMETHYLATED probe of cg10602367, Color = Red

### Extract the Red and Green fluorescence values for the Team 1 address.
Red_fluorescences <- Red[rownames(Red) == TARGET_ADDRESS, , drop = FALSE]
Red_fluorescences

Green_fluorescences <- Green[rownames(Green) == TARGET_ADDRESS, , drop = FALSE]
Green_fluorescences

### Build a summary table with the fluorescence values and probe annotation.
df_summary_probes <- data.frame(
  Sample             = colnames(Red_fluorescences),
  Red_Fluorescence   = as.numeric(Red_fluorescences[1, ]),
  Green_Fluorescence = as.numeric(Green_fluorescences[1, ]),
  Type               = "I",
  Color              = "Red"
)

cat("Team 1 address fluorescence summary:\n")
print(df_summary_probes)
write.csv(
  df_summary_probes,
  file.path(tablesDir, "step03_team1_address_fluorescence.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step03_team1_address_fluorescence.csv"))

### Save the Team 1 probe annotation identified from the manifest lookup.
team1_probe_match <- ProbeInfo_I[ProbeInfo_I$AddressA == TARGET_ADDRESS, , drop = FALSE]
team1_probe_annotation <- data.frame(
  ProbeID = team1_probe_match$Name,
  Address = TARGET_ADDRESS,
  Address_Match = "AddressA",
  Probe_Type = "I",
  Signal = "Unmethylated",
  Color = if ("Color" %in% colnames(team1_probe_match)) team1_probe_match$Color else NA
)

write.csv(
  team1_probe_annotation,
  file.path(tablesDir, "step03_team1_probe_annotation.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step03_team1_probe_annotation.csv"))


###############################################################################
### STEP 04: CREATE THE MSet.raw OBJECT
###############################################################################

cat("\n--- Running Step 04: Create the MSet.raw Object ---\n")

### Convert RGset to a raw MethylSet object without normalization.
MSet.raw <- preprocessRaw(RGset)
MSet.raw

### Save MSet.raw for later workflow steps.
save(MSet.raw, file = file.path(rdsDir, "MSet_raw.RData"))
report_saved("RData object", file.path(rdsDir, "MSet_raw.RData"))

### Extract methylated and unmethylated signal matrices.
Meth <- as.matrix(getMeth(MSet.raw))
str(Meth)
head(Meth)

Unmeth <- as.matrix(getUnmeth(MSet.raw))
str(Unmeth)
head(Unmeth)

### Check probe cg10602367, identified from the Team 1 address in Step 03.
Meth[rownames(Meth) == "cg10602367", , drop = FALSE]
Unmeth[rownames(Unmeth) == "cg10602367", , drop = FALSE]


###############################################################################
### STEP 05: QUALITY CHECKS
###############################################################################

cat("\n--- Running Step 05: Quality Checks ---\n")

### Team 1 detection p-value threshold: 0.05.
DETECTION_P_THRESHOLD <- 0.05

### 5a. QC plot
### Calculate median methylated and unmethylated signal intensities.
qc <- getQC(MSet.raw)
qc
qc_summary <- data.frame(
  Sample = rownames(qc),
  qc,
  row.names = NULL
)

write.csv(
  qc_summary,
  file.path(tablesDir, "step05_qc_summary.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step05_qc_summary.csv"))

### Plot the QC results.
png(file.path(figuresDir, "step05_qcplot.png"),
    width = 7, height = 6, units = "in", res = 300, pointsize = 11)
plotQC(qc)
dev.off()
report_saved("Figure", file.path(figuresDir, "step05_qcplot.png"))

### 5b. Negative control intensities
### Extract control probes from the array and count the control types.
NegativeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
table(NegativeControl$Type)
head(NegativeControl)

### Plot negative control intensities.
png(file.path(figuresDir, "step05_controlstripplot.png"),
    width = 7, height = 6, units = "in", res = 300, pointsize = 11)
controlStripPlot(RGset, controls = "NEGATIVE")
dev.off()
report_saved("Figure", file.path(figuresDir, "step05_controlstripplot.png"))


### 5c. Detection p-values and failed positions per sample
### Calculate detection p-values and save them for later checks.
detP <- detectionP(RGset)
save(detP, file = file.path(rdsDir, "detP.RData"))
report_saved("RData object", file.path(rdsDir, "detP.RData"))

### Identify failed probe-sample tests using the Team 1 threshold.
Over_Threshold <- detP > DETECTION_P_THRESHOLD
head(Over_Threshold)
table(Over_Threshold)
sum(Over_Threshold)

Total_Failed_Positions <- sum(Over_Threshold)
Total_Detection_Tests <- length(Over_Threshold)
Total_Failed_Percentage <- Total_Failed_Positions / Total_Detection_Tests * 100

failed_probe_message <- sprintf(
  "There are %s failed detection tests with a detection p-value higher than %.2f",
  Total_Failed_Positions,
  DETECTION_P_THRESHOLD
)

### Count the number and percentage of failed positions per sample.
Over_Threshold_Per_Sample <- colSums(Over_Threshold)
Over_Threshold_Per_Sample

Over_Threshold_Fraction_Per_Sample <- colMeans(Over_Threshold)
Over_Threshold_Fraction_Per_Sample

failed <- data.frame(
  Sample = names(Over_Threshold_Per_Sample),
  n_Failed_Positions = Over_Threshold_Per_Sample,
  Fraction_Failed_Positions = Over_Threshold_Fraction_Per_Sample,
  Percentage_Failed_Positions = Over_Threshold_Fraction_Per_Sample * 100
)
cat("Failed detection tests summary:\n")
print(failed)
print(failed_probe_message)
write.csv(
  failed,
  file.path(tablesDir, "step05_failed_detection_by_sample.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step05_failed_detection_by_sample.csv"))




###############################################################################
### STEP 06: CALCULATE RAW BETA AND M VALUES
###############################################################################

cat("\n--- Running Step 06: Calculate Raw Beta and M Values ---\n")

### getBeta() and getM() retrieve beta and M value matrices from MSet.raw.

### Extract raw beta values and inspect their dimensions and distribution.
beta <- getBeta(MSet.raw)
class(beta)
dim(beta)
head(beta)
summary(beta)

### Extract raw M values and inspect their dimensions and distribution.
M <- getM(MSet.raw)
dim(M)
head(M)
summary(M)

### Save raw beta and M matrices.
saveRDS(beta, file = file.path(rdsDir, "beta_raw.rds"))
saveRDS(M, file = file.path(rdsDir, "M_raw.rds"))
report_saved("RDS object", file.path(rdsDir, "beta_raw.rds"))
report_saved("RDS object", file.path(rdsDir, "M_raw.rds"))

### Recover phenotype data to split samples by group.
pheno <- data.frame(pData(RGset))
str(pheno)
table(pheno$Group)

### Split beta values into CTRL and DIS groups.
beta_CTRL <- beta[, pheno$Group == "CTRL"]
beta_DIS <- beta[, pheno$Group == "DIS"]

### Split M values into CTRL and DIS groups.
M_CTRL <- M[, pheno$Group == "CTRL"]
M_DIS <- M[, pheno$Group == "DIS"]

dim(beta_CTRL)
dim(beta_DIS)
dim(M_CTRL)
dim(M_DIS)
cat(sprintf("Raw beta matrix dimensions: %s probes x %s samples\n",
            nrow(beta), ncol(beta)))
cat(sprintf("Raw M matrix dimensions: %s probes x %s samples\n",
            nrow(M), ncol(M)))

### Calculate mean beta and M values per CpG within each group.
mean_of_beta_CTRL <- apply(beta_CTRL, 1, mean, na.rm = TRUE)
mean_of_beta_DIS <- apply(beta_DIS, 1, mean, na.rm = TRUE)

mean_of_M_CTRL <- apply(M_CTRL, 1, mean, na.rm = TRUE)
mean_of_M_DIS <- apply(M_DIS, 1, mean, na.rm = TRUE)

### Remove infinite M values before density estimation.
mean_of_M_CTRL <- mean_of_M_CTRL[is.finite(mean_of_M_CTRL)]
mean_of_M_DIS <- mean_of_M_DIS[is.finite(mean_of_M_DIS)]

### Estimate density distributions for raw beta and M means.
d_mean_of_beta_CTRL <- density(mean_of_beta_CTRL, na.rm = TRUE)
d_mean_of_beta_DIS <- density(mean_of_beta_DIS, na.rm = TRUE)

d_mean_of_M_CTRL <- density(mean_of_M_CTRL, na.rm = TRUE)
d_mean_of_M_DIS <- density(mean_of_M_DIS, na.rm = TRUE)

### Plot raw beta density distributions by group.
png(file.path(figuresDir, "step06_raw_beta_density_CTRL_DIS.png"),
    width = 7, height = 6, units = "in", res = 300, pointsize = 11)

plot(d_mean_of_beta_CTRL, col = "blue",
     main = "Raw beta values", xlab = "Mean beta")

lines(d_mean_of_beta_DIS, col = "red")

legend("topright", legend = c("CTRL", "DIS"),
       col = c("blue", "red"), lty = 1, bty = "n")

dev.off()
report_saved("Figure", file.path(figuresDir, "step06_raw_beta_density_CTRL_DIS.png"))

### Plot raw M value density distributions by group.
png(file.path(figuresDir, "step06_raw_M_density_CTRL_DIS.png"),
    width = 7, height = 6, units = "in", res = 300, pointsize = 11)

plot(d_mean_of_M_CTRL, col = "blue",
     main = "Raw M values", xlab = "Mean M")

lines(d_mean_of_M_DIS, col = "red")

legend("topright", legend = c("CTRL", "DIS"),
       col = c("blue", "red"), lty = 1, bty = "n")

dev.off()
report_saved("Figure", file.path(figuresDir, "step06_raw_M_density_CTRL_DIS.png"))

### Save summary vectors used for the density plots.
saveRDS(mean_of_beta_CTRL,
        file = file.path(rdsDir, "mean_of_beta_raw_CTRL.rds"))
report_saved("RDS object", file.path(rdsDir, "mean_of_beta_raw_CTRL.rds"))

saveRDS(mean_of_beta_DIS,
        file = file.path(rdsDir, "mean_of_beta_raw_DIS.rds"))
report_saved("RDS object", file.path(rdsDir, "mean_of_beta_raw_DIS.rds"))

saveRDS(mean_of_M_CTRL,
        file = file.path(rdsDir, "mean_of_M_raw_CTRL.rds"))
report_saved("RDS object", file.path(rdsDir, "mean_of_M_raw_CTRL.rds"))

saveRDS(mean_of_M_DIS,
        file = file.path(rdsDir, "mean_of_M_raw_DIS.rds"))
report_saved("RDS object", file.path(rdsDir, "mean_of_M_raw_DIS.rds"))

###############################################################################
### STEP 07: NORMALIZATION USING preprocessNoob
###############################################################################

cat("\n--- Running Step 07: Normalization Using preprocessNoob ---\n")

### Team 1 normalization method: preprocessNoob.

### 7a. Raw beta distributions by probe chemistry
### Retrieve Type I and Type II probe annotations from the manifest.
dfI <- data.frame(getProbeInfo(RGset, type = "I"))
dfII <- data.frame(getProbeInfo(RGset, type = "II"))

dim(dfI)
dim(dfII)
str(dfI)
str(dfII)
head(dfI$Name)
head(dfII$Name)

### Split raw beta values by probe chemistry.
beta_I <- beta[rownames(beta) %in% dfI$Name, ]
beta_II <- beta[rownames(beta) %in% dfII$Name, ]

dim(beta_I)
dim(beta_II)

### Calculate mean raw beta by probe chemistry.
mean_of_beta_I <- apply(beta_I, 1, mean, na.rm = TRUE)
mean_of_beta_II <- apply(beta_II, 1, mean, na.rm = TRUE)

### Estimate density of raw beta means by probe chemistry.
d_mean_of_beta_I <- density(mean_of_beta_I, na.rm = TRUE)
d_mean_of_beta_II <- density(mean_of_beta_II, na.rm = TRUE)

### Calculate and estimate density of raw beta standard deviations.
sd_of_beta_I <- apply(beta_I, 1, sd, na.rm = TRUE)
sd_of_beta_II <- apply(beta_II, 1, sd, na.rm = TRUE)

d_sd_of_beta_I <- density(sd_of_beta_I, na.rm = TRUE)
d_sd_of_beta_II <- density(sd_of_beta_II, na.rm = TRUE)

### 7b. Normalization with preprocessNoob
### Normalize the raw RGChannelSet with preprocessNoob.
preprocessNoob_results <- preprocessNoob(RGset)

### Inspect the normalized object.
str(preprocessNoob_results)
class(preprocessNoob_results)
preprocessNoob_results

### Extract normalized beta and M values.
beta_preprocessNoob <- getBeta(preprocessNoob_results)
M_preprocessNoob <- getM(preprocessNoob_results)

head(beta_preprocessNoob)
dim(beta_preprocessNoob)

### Save normalized objects for reproducibility.
saveRDS(preprocessNoob_results, file = file.path(rdsDir, "preprocessNoob_results.rds"))
saveRDS(beta_preprocessNoob, file = file.path(rdsDir, "beta_preprocessNoob.rds"))
saveRDS(M_preprocessNoob, file = file.path(rdsDir, "M_preprocessNoob.rds"))
report_saved("RDS object", file.path(rdsDir, "preprocessNoob_results.rds"))
report_saved("RDS object", file.path(rdsDir, "beta_preprocessNoob.rds"))
report_saved("RDS object", file.path(rdsDir, "M_preprocessNoob.rds"))

### 7c. Normalized beta distributions by probe chemistry
### Split normalized beta values by probe chemistry.
beta_preprocessNoob_I <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfI$Name, ]
beta_preprocessNoob_II <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfII$Name, ]

### Calculate mean normalized beta values by probe chemistry.
mean_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, mean, na.rm = TRUE)
mean_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, mean, na.rm = TRUE)

### Estimate density of normalized beta means.
d_mean_of_beta_preprocessNoob_I <- density(mean_of_beta_preprocessNoob_I, na.rm = TRUE)
d_mean_of_beta_preprocessNoob_II <- density(mean_of_beta_preprocessNoob_II, na.rm = TRUE)

### Calculate and estimate density of normalized beta standard deviations.
sd_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, sd, na.rm = TRUE)
sd_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, sd, na.rm = TRUE)

d_sd_of_beta_preprocessNoob_I <- density(sd_of_beta_preprocessNoob_I, na.rm = TRUE)
d_sd_of_beta_preprocessNoob_II <- density(sd_of_beta_preprocessNoob_II, na.rm = TRUE)

### 7d. Six-panel raw vs normalized comparison
### Define colors for the raw and normalized beta boxplots.
pheno <- data.frame(pData(RGset))
pheno$Group

type_col <- c("Type I" = "blue", "Type II" = "red")
group_col <- c(CTRL = "violet", DIS = "orange")
boxplot_col <- group_col[as.character(pheno$Group)]
boxplot_col

### Compare raw and normalized beta values using plots.
png(file.path(figuresDir, "step07_raw_vs_preprocessNoob_6panel.png"),
    width = 14, height = 8, units = "in", res = 300, pointsize = 11)
par(mfrow = c(2, 3), mar = c(7, 4.5, 3, 1), oma = c(1, 1, 1, 1))

### Panel 1: Raw beta mean density by probe type.
plot(d_mean_of_beta_I, col = type_col["Type I"], main = "Raw beta")
lines(d_mean_of_beta_II, col = type_col["Type II"])
legend("topright", legend = names(type_col), col = type_col,
       lty = 1, bty = "n", cex = 0.8)

### Panel 2: Raw beta standard deviation density by probe type.
plot(d_sd_of_beta_I, col = type_col["Type I"], main = "Raw SD")
lines(d_sd_of_beta_II, col = type_col["Type II"])
legend("topright", legend = names(type_col), col = type_col,
       lty = 1, bty = "n", cex = 0.8)

### Panel 3: Raw beta boxplot by sample group.
boxplot(beta, ylim = c(0, 1), col = boxplot_col, main = "Raw beta",
        names = FALSE)
axis(1, at = seq_len(ncol(beta)), labels = seq_len(ncol(beta)), cex.axis = 0.8)
usr <- par("usr")
legend(x = mean(usr[1:2]), y = usr[3] - 0.16 * diff(usr[3:4]),
       legend = names(group_col), xjust = 0.5, yjust = 1,
       fill = group_col, title = "Samples", horiz = TRUE, bty = "n",
       cex = 0.85, xpd = NA)

### Panel 4: Noob-normalized beta mean density by probe type.
plot(d_mean_of_beta_preprocessNoob_I, col = type_col["Type I"],
     main = "Noob beta")
lines(d_mean_of_beta_preprocessNoob_II, col = type_col["Type II"])
legend("topright", legend = names(type_col), col = type_col,
       lty = 1, bty = "n", cex = 0.8)

### Panel 5: Noob-normalized beta standard deviation density by probe type.
plot(d_sd_of_beta_preprocessNoob_I, col = type_col["Type I"],
     main = "Noob SD")
lines(d_sd_of_beta_preprocessNoob_II, col = type_col["Type II"])
legend("topright", legend = names(type_col), col = type_col,
       lty = 1, bty = "n", cex = 0.8)

### Panel 6: Noob-normalized beta boxplot by sample group.
boxplot(beta_preprocessNoob, ylim = c(0, 1), col = boxplot_col,
        main = "Noob beta", names = FALSE)
axis(1, at = seq_len(ncol(beta_preprocessNoob)), labels = seq_len(ncol(beta_preprocessNoob)), cex.axis = 0.8)
usr <- par("usr")
legend(x = mean(usr[1:2]), y = usr[3] - 0.16 * diff(usr[3:4]),
       legend = names(group_col), xjust = 0.5, yjust = 1,
       fill = group_col, title = "Samples", horiz = TRUE, bty = "n",
       cex = 0.85, xpd = NA)
dev.off()
report_saved("Figure", file.path(figuresDir, "step07_raw_vs_preprocessNoob_6panel.png"))


###############################################################################
### STEP 08: PCA ON NORMALIZED BETA VALUES
###############################################################################

cat("\n--- Running Step 08: PCA on Normalized Beta Values ---\n")

### Extract sample information from RGset.
sample_info <- data.frame(pData(RGset))
sample_info$Group <- factor(sample_info$Group)
sample_info$Sex <- factor(sample_info$Sex)

### In pData(RGset), Sentrix_ID is stored as Slide.
sample_info$Sentrix_ID <- factor(sample_info$Slide)

### 8a. Run PCA on normalized beta values.
pca_results <- prcomp(t(beta_preprocessNoob), scale. = TRUE)

### Calculate the proportion of variance explained by each component.
pca_var <- pca_results$sdev^2
pca_var_explained <- pca_var / sum(pca_var)
pca_var_cumulative <- cumsum(pca_var_explained)
cat(sprintf("PCA variance explained: PC1 = %.1f%%, PC2 = %.1f%%\n",
            pca_var_explained[1] * 100, pca_var_explained[2] * 100))

pca_variance_table <- data.frame(
  Principal_Component = paste0("PC", seq_along(pca_var_explained)),
  Variance_Explained = pca_var_explained,
  Percentage_Explained = pca_var_explained * 100,
  Cumulative_Variance = pca_var_cumulative,
  Cumulative_Percentage = pca_var_cumulative * 100
)

write.csv(
  pca_variance_table,
  file.path(tablesDir, "step08_pca_variance_explained.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step08_pca_variance_explained.csv"))

pca_sample_scores <- data.frame(
  Sample = rownames(pca_results$x),
  PC1 = pca_results$x[, "PC1"],
  PC2 = pca_results$x[, "PC2"],
  Group = sample_info$Group,
  Sex = sample_info$Sex,
  Sentrix_ID = sample_info$Sentrix_ID,
  row.names = NULL
)

write.csv(
  pca_sample_scores,
  file.path(tablesDir, "step08_pca_sample_scores.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step08_pca_sample_scores.csv"))

### Define shared axis limits for the PCA scatter plots.
x_padding <- diff(range(pca_results$x[, "PC1"])) * 0.12
y_padding <- diff(range(pca_results$x[, "PC2"])) * 0.12
xlim_pca <- range(pca_results$x[, "PC1"]) + c(-x_padding, x_padding)
ylim_pca <- range(pca_results$x[, "PC2"]) + c(-y_padding, y_padding)

### 8b. Create scree plot.
png(file.path(figuresDir, "step08_PCA_ScreePlot.png"),
    width = 8, height = 6, units = "in", res = 300, pointsize = 11)
bp <- barplot(
  pca_var_explained[1:8],
  names.arg = paste0("PC", 1:8),
  ylim = c(0, 1),
  main = "Variance Explained by Principal Components",
  xlab = "Principal Components",
  ylab = "Proportion of Variance"
)

text(
  bp,
  pca_var_explained[1:8],
  labels = paste0(round(pca_var_explained[1:8] * 100, 1), "%"),
  pos = 3,
  cex = 0.8
)

lines(
  x = bp,
  y = pca_var_cumulative[1:8],
  type = "b",
  pch = 16,
  col = "red"
)

text(
  bp,
  pca_var_cumulative[1:8],
  labels = paste0(round(pca_var_cumulative[1:8] * 100, 1), "%"),
  pos = 1,
  cex = 0.7,
  col = "red"
)

legend(
  x = bp[6],
  y = 0.78,
  legend = c("Variance explained", "Cumulative variance"),
  fill = c("grey", NA),
  border = c("black", NA),
  lty = c(NA, 1),
  pch = c(NA, 16),
  col = c("black", "red"),
  bty = "n"
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step08_PCA_ScreePlot.png"))

### 8c. Create PCA plot colored by disease group.
png(file.path(figuresDir, "step08_PCA_Group.png"),
    width = 8, height = 6, units = "in", res = 300, pointsize = 11)

plot(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  col = sample_info$Group,
  pch = 17,
  cex = 1.5,
  xlim = xlim_pca,
  ylim = ylim_pca,
  xlab = paste0(
    "PC1 (",
    round(pca_var_explained[1] * 100, 1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(pca_var_explained[2] * 100, 1),
    "%)"
  ),
  main = "PCA of Normalized Beta Values by Group"
)

text(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  labels = sample_info$SampleID,
  pos = 3,
  cex = 0.8
)

legend(
  "topright",
  legend = levels(sample_info$Group),
  col = seq_along(levels(sample_info$Group)),
  pch = 17
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step08_PCA_Group.png"))

### 8d. Create PCA plot colored by sex.
png(file.path(figuresDir, "step08_PCA_by_Sex.png"),
    width = 8, height = 6, units = "in", res = 300, pointsize = 11)

plot(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  col = sample_info$Sex,
  pch = 17,
  cex = 1.5,
  xlim = xlim_pca,
  ylim = ylim_pca,
  xlab = paste0(
    "PC1 (",
    round(pca_var_explained[1] * 100, digits = 1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(pca_var_explained[2] * 100, digits = 1),
    "%)"
  ),
  main = "PCA of Normalized Beta Values by Sex"
)

text(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  labels = sample_info$SampleID,
  pos = 3,
  cex = 0.8
)

legend(
  "topright",
  legend = levels(sample_info$Sex),
  col = seq_along(levels(sample_info$Sex)),
  pch = 17
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step08_PCA_by_Sex.png"))

### 8e. Create PCA plot colored by Sentrix_ID / slide batch.
png(file.path(figuresDir, "step08_PCA_by_SentrixID.png"),
    width = 8, height = 6, units = "in", res = 300, pointsize = 11)

plot(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  col = sample_info$Sentrix_ID,
  pch = 17,
  cex = 1.5,
  xlim = xlim_pca,
  ylim = ylim_pca,
  xlab = paste0("PC1 (", round(pca_var_explained[1] * 100, 1), "%)"),
  ylab = paste0("PC2 (", round(pca_var_explained[2] * 100, 1), "%)"),
  main = "PCA of Normalized Beta Values by Sentrix_ID"
)

text(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  labels = sample_info$SampleID,
  pos = 3,
  cex = 0.8
)

legend(
  "topright",
  legend = levels(sample_info$Sentrix_ID),
  col = seq_along(levels(sample_info$Sentrix_ID)),
  pch = 17
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step08_PCA_by_SentrixID.png"))

###############################################################################
### STEP 09: DIFFERENTIAL METHYLATION ANALYSIS USING T-TEST
###############################################################################

cat("\n--- Running Step 09: Differential Methylation Analysis Using t-test ---\n")

### Team 1 differential methylation test: t-test.

### Use simple object names for the t-test analysis.
beta_matrix <- beta_preprocessNoob
sample_data <- sample_info

### Identify CTRL and DIS samples.
control_samples <- sample_data$Group == "CTRL"
disease_samples <- sample_data$Group == "DIS"

### Calculate the mean beta value for each probe in the CTRL group.
mean_ctrl <- rowMeans(
  beta_matrix[, control_samples],
  na.rm = TRUE
)

### Calculate the mean beta value for each probe in the DIS group.
mean_dis <- rowMeans(
  beta_matrix[, disease_samples],
  na.rm = TRUE
)

### Calculate the methylation difference between the two groups.
delta_beta <- mean_dis - mean_ctrl

### Define a function to perform a t-test for one CpG probe.
My_ttest_function <- function(probe_values) {
  t_test <- t.test(
    probe_values[control_samples],
    probe_values[disease_samples]
  )
  return(t_test$p.value)
}

### Apply the t-test to each CpG probe and collect the p-values.
p_values <- apply(
  beta_matrix,
  1,
  My_ttest_function
)

### Create a table containing the differential methylation results.
t_test_results <- data.frame(
  ProbeID = rownames(beta_matrix),
  Mean_CTRL = mean_ctrl,
  Mean_DIS = mean_dis,
  Delta_Beta = delta_beta,
  P_Value = p_values
)

### Count probes significant at the nominal p-value threshold.
Nominal_Significant_Probes <- sum(t_test_results$P_Value <= 0.05, na.rm = TRUE)
cat(sprintf("Nominally significant probes (P_Value <= 0.05): %s\n",
            Nominal_Significant_Probes))

### Save the results for the next step of the pipeline.
write.csv(
  t_test_results,
  file.path(tablesDir, "step09_t_test_results.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step09_t_test_results.csv"))

###############################################################################
### STEP 10: MULTIPLE TESTING CORRECTION
###############################################################################

cat("\n--- Running Step 10: Multiple Testing Correction ---\n")

### Team 1 significance threshold for multiple testing correction.
alpha <- 0.05

### Apply Bonferroni and Benjamini-Hochberg corrections to the t-test p-values.
t_test_results$P_Bonferroni <- p.adjust(t_test_results$P_Value, method = "bonferroni")
t_test_results$P_BH <- p.adjust(t_test_results$P_Value, method = "BH")

head(t_test_results[, c("ProbeID", "P_Value", "P_Bonferroni", "P_BH")])

### Visualize the distribution of nominal and corrected p-values.
png(
  file.path(figuresDir, "step10_pvalue_correction_boxplot.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)

boxplot(
  t_test_results[, c("P_Value", "P_Bonferroni", "P_BH")],
  names = c("Nominal", "Bonferroni", "BH"),
  col = c("grey80", "lightblue", "lightgreen"),
  ylab = "P-value",
  main = "Nominal and Corrected P-values"
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step10_pvalue_correction_boxplot.png"))

### Count significant probes before and after multiple testing correction.
sig_nominal <- sum(t_test_results$P_Value <= alpha, na.rm = TRUE)
sig_bonferroni <- sum(t_test_results$P_Bonferroni <= alpha, na.rm = TRUE)
sig_bh <- sum(t_test_results$P_BH <= alpha, na.rm = TRUE)

correction_summary <- data.frame(
  Method = c("Nominal p-value", "Bonferroni", "Benjamini-Hochberg"),
  Threshold = alpha,
  Significant_Probes = c(sig_nominal, sig_bonferroni, sig_bh)
)

cat("Multiple testing correction summary:\n")
print(correction_summary)
write.csv(
  correction_summary,
  file.path(tablesDir, "step10_correction_summary.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step10_correction_summary.csv"))

### Save the table containing nominal and corrected p-values.
write.csv(
  t_test_results,
  file.path(tablesDir, "step10_t_test_results_adjusted.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step10_t_test_results_adjusted.csv"))


###############################################################################
### STEP 11: VOLCANO AND MANHATTAN PLOTS
###############################################################################

cat("\n--- Running Step 11: Volcano and Manhattan Plots ---\n")

### 11a. Volcano plot.
delta_beta_threshold <- 0.1

hyper <- which(
  t_test_results$Delta_Beta > delta_beta_threshold &
    t_test_results$P_Value <= alpha
)

hypo <- which(
  t_test_results$Delta_Beta < -delta_beta_threshold &
    t_test_results$P_Value <= alpha
)

volcano_summary <- data.frame(
  Category = c("Hypermethylated in DIS", "Hypomethylated in DIS"),
  Delta_Beta_Threshold = c(delta_beta_threshold, -delta_beta_threshold),
  P_Value_Threshold = alpha,
  Candidate_Probes = c(length(hyper), length(hypo))
)

cat("Volcano plot exploratory candidate summary:\n")
print(volcano_summary)
volcano_candidates <- rbind(
  if (length(hyper) > 0) {
    data.frame(
      Category = rep("Hypermethylated in DIS", length(hyper)),
      t_test_results[hyper, ],
      row.names = NULL
    )
  },
  if (length(hypo) > 0) {
    data.frame(
      Category = rep("Hypomethylated in DIS", length(hypo)),
      t_test_results[hypo, ],
      row.names = NULL
    )
  }
)

write.csv(
  volcano_candidates,
  file.path(tablesDir, "step11_volcano_candidates.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step11_volcano_candidates.csv"))

### Produce the volcano plot using nominal p-values.
png(
  file.path(figuresDir, "step11_volcano_plot.png"),
  width = 8,
  height = 6,
  units = "in",
  res = 300
)

plot(
  x = t_test_results$Delta_Beta,
  y = -log10(t_test_results$P_Value),
  pch = 16,
  cex = 0.5,
  col = "darkgray",
  xlab = "Delta Beta (DIS - CTRL)",
  ylab = "-log10(Nominal P-value)",
  main = "Volcano Plot of DNA Methylation Differences",
  xlim = c(
    -max(abs(t_test_results$Delta_Beta), na.rm = TRUE),
    max(abs(t_test_results$Delta_Beta), na.rm = TRUE)
  )
)

points(
  t_test_results$Delta_Beta[hyper],
  -log10(t_test_results$P_Value)[hyper],
  col = "red",
  pch = 16,
  cex = 0.6
)

points(
  t_test_results$Delta_Beta[hypo],
  -log10(t_test_results$P_Value)[hypo],
  col = "blue",
  pch = 16,
  cex = 0.6
)

abline(h = -log10(alpha), col = "black", lty = 2)
abline(v = c(-delta_beta_threshold, delta_beta_threshold), col = "black", lty = 2)

legend(
  "topright",
  legend = c("Hypermethylated in DIS", "Hypomethylated in DIS"),
  col = c("red", "blue"),
  pch = 16,
  cex = 0.8,
  bty = "n"
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step11_volcano_plot.png"))

### 11b. Manhattan plot.
ann <- getAnnotation(preprocessNoob_results)
locations <- data.frame(
  ProbeID = rownames(ann),
  CHR = ann$chr,
  MAPINFO = ann$pos
)

manhattan_data <- merge(t_test_results, locations, by = "ProbeID")

order_chr <- c(as.character(1:22), "X", "Y")
manhattan_data$CHR <- gsub("chr", "", manhattan_data$CHR)
manhattan_data$CHR <- factor(manhattan_data$CHR, levels = order_chr)
manhattan_data$CHR <- as.numeric(manhattan_data$CHR)

manhattan_data <- manhattan_data[
  !is.na(manhattan_data$CHR) &
    !is.na(manhattan_data$MAPINFO) &
    !is.na(manhattan_data$P_Value),
]

png(
  file.path(figuresDir, "step11_manhattan_plot.png"),
  width = 10,
  height = 6,
  units = "in",
  res = 300
)

manhattan(
  manhattan_data,
  chr = "CHR",
  bp = "MAPINFO",
  snp = "ProbeID",
  p = "P_Value",
  col = rainbow(length(order_chr)),
  suggestiveline = -log10(0.00001),
  genomewideline = FALSE,
  annotatePval = 0.00001,
  main = "Manhattan Plot of Epigenome-Wide Association",
  ylab = "-log10(Nominal P-value)"
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step11_manhattan_plot.png"))


###############################################################################
### STEP 12: HEATMAP OF THE TOP 100 CpG PROBES
###############################################################################

cat("\n--- Running Step 12: Heatmap of Top 100 CpGs ---\n")

### Select the top 100 CpG probes ranked by nominal t-test p-value.
top_100_results <- t_test_results[order(t_test_results$P_Value), ][1:100, ]
top_100_probes <- top_100_results$ProbeID

cat("Top 100 CpGs selected by nominal p-value for heatmaps:\n")
print(head(top_100_results))
write.csv(
  top_100_results,
  file.path(tablesDir, "step12_top100_heatmap_probes.csv"),
  row.names = FALSE
)
report_saved("Table", file.path(tablesDir, "step12_top100_heatmap_probes.csv"))

### Extract normalized beta values for the selected probes.
beta_top100 <- beta_preprocessNoob[top_100_probes, ]
input_heatmap <- as.matrix(beta_top100)

dim(input_heatmap)
cat(sprintf("Heatmap matrix dimensions: %s probes x %s samples\n",
            nrow(input_heatmap), ncol(input_heatmap)))

### Define sample group colors for the heatmap column annotation.
heatmap_group_colors <- c(CTRL = "royalblue", DIS = "orange")
colorbar <- heatmap_group_colors[as.character(sample_info$Group)]
heatmap_palette <- colorRampPalette(c("green", "black", "red"))(100)

### Produce the heatmap using complete linkage clustering.
png(
  file.path(figuresDir, "step12_heatmap_top100_complete.png"),
  width = 8,
  height = 8,
  units = "in",
  res = 300
)

heatmap.2(
  input_heatmap,
  col = heatmap_palette,
  Rowv = TRUE,
  Colv = TRUE,
  dendrogram = "both",
  key = TRUE,
  ColSideColors = colorbar,
  density.info = "none",
  trace = "none",
  scale = "none",
  symm = FALSE,
  margins = c(10, 8),
  cexCol = 0.8,
  cexRow = 0.35,
  labRow = rownames(input_heatmap),
  main = "Top 100 CpGs - Complete Linkage"
)

legend(
  "topright",
  legend = names(heatmap_group_colors),
  fill = heatmap_group_colors,
  border = FALSE,
  bty = "n",
  cex = 0.8
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step12_heatmap_top100_complete.png"))

### Produce the heatmap using single linkage clustering.
png(
  file.path(figuresDir, "step12_heatmap_top100_single.png"),
  width = 8,
  height = 8,
  units = "in",
  res = 300
)

heatmap.2(
  input_heatmap,
  col = heatmap_palette,
  Rowv = TRUE,
  Colv = TRUE,
  hclustfun = function(x) hclust(x, method = "single"),
  dendrogram = "both",
  key = TRUE,
  ColSideColors = colorbar,
  density.info = "none",
  trace = "none",
  scale = "none",
  symm = FALSE,
  margins = c(10, 8),
  cexCol = 0.8,
  cexRow = 0.35,
  labRow = rownames(input_heatmap),
  main = "Top 100 CpGs - Single Linkage"
)

legend(
  "topright",
  legend = names(heatmap_group_colors),
  fill = heatmap_group_colors,
  border = FALSE,
  bty = "n",
  cex = 0.8
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step12_heatmap_top100_single.png"))

### Produce the heatmap using average linkage clustering.
png(
  file.path(figuresDir, "step12_heatmap_top100_average.png"),
  width = 8,
  height = 8,
  units = "in",
  res = 300
)

heatmap.2(
  input_heatmap,
  col = heatmap_palette,
  Rowv = TRUE,
  Colv = TRUE,
  hclustfun = function(x) hclust(x, method = "average"),
  dendrogram = "both",
  key = TRUE,
  ColSideColors = colorbar,
  density.info = "none",
  trace = "none",
  scale = "none",
  symm = FALSE,
  margins = c(10, 8),
  cexCol = 0.8,
  cexRow = 0.35,
  labRow = rownames(input_heatmap),
  main = "Top 100 CpGs - Average Linkage"
)

legend(
  "topright",
  legend = names(heatmap_group_colors),
  fill = heatmap_group_colors,
  border = FALSE,
  bty = "n",
  cex = 0.8
)

dev.off()
report_saved("Figure", file.path(figuresDir, "step12_heatmap_top100_average.png"))

cat("\n--- Pipeline Execution Complete ---\n")
