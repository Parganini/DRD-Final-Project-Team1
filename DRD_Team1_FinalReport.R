## Team 1 assignments:
## - Step 3 address: 42796479
## - Step 5 detection p-value threshold: 0.05
## - Step 7 normalization method: preprocessNoob
## - Step 9 differential methylation test: t-test

# Step 00. Setup

### Clean the workspace before starting the analysis.
rm(list = ls())


### Load all packages used by the complete workflow.

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

###############################################################################
### STEP 01: load raw data with minfi
###############################################################################

list.files("./data/raw/")
SampleSheet <- read.csv("./data/raw/SampleSheet_Report_II.csv", header = TRUE)
head(SampleSheet)

targets <- read.metharray.sheet("./data/raw")
head(targets)

RGset <- read.metharray.exp(targets = targets)

save(RGset, file = file.path(rdsDir, "RGset_Report.RData"))
load(file.path(rdsDir, "RGset_Report.RData"))

RGset
str(RGset)


###############################################################################
### STEP 02: extract Red and Green fluorescence values
###############################################################################

Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)

Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

###############################################################################
### STEP 3. PROBE INFORMATION AND ADDRESS LOOKUP
### Group 1 address: 42796479
###############################################################################

TARGET_ADDRESS <- "42796479"

# --- Probe type lookup ---
getManifest(RGset)
ProbeInfo_I  <- data.frame(getProbeInfo(RGset, type = "I"))
head(ProbeInfo_I)

ProbeInfo_II <- data.frame(getProbeInfo(RGset, type = "II"))
head(ProbeInfo_II)

ProbeInfo_I[ProbeInfo_I$AddressA == TARGET_ADDRESS, ]
ProbeInfo_I[ProbeInfo_I$AddressB == TARGET_ADDRESS, ]
ProbeInfo_II[ProbeInfo_II$AddressA == TARGET_ADDRESS, ]
# Result: Type I, AddressA → UNMETHYLATED probe of cg10602367, Color = Red

# --- Extract fluorescence and build summary ---
Red_fluorescences   <- Red[rownames(Red) == TARGET_ADDRESS, ]
Red_fluorescences

Green_fluorescences <- Green[rownames(Green) == TARGET_ADDRESS, ]
Green_fluorescences

df_summary_probes <- data.frame(
  Sample             = colnames(Red_fluorescences),
  Red_Fluorescence   = as.numeric(Red_fluorescences[1, ]),
  Green_Fluorescence = as.numeric(Green_fluorescences[1, ]),
  Type               = "I",
  Color              = "Red"
)

print(df_summary_probes)


# ============================================================
# STEP 4. Creation of the MSet.raw object
# ============================================================

# Convert RGset to MethylSet (raw, no normalization)
MSet.raw <- preprocessRaw(RGset)
MSet.raw


# Save for later use
save(MSet.raw, file = file.path(rdsDir, "MSet_raw.RData"))

# Extract methylated and unmethylated signal matrices
Meth <- as.matrix(getMeth(MSet.raw))
str(Meth)
head(Meth)

Unmeth <- as.matrix(getUnmeth(MSet.raw))
str(Unmeth)
head(Unmeth)

# ============================================================
# Check probe cg10602367 (our Group 1 probe from Step 3)
# ============================================================
Meth[rownames(Meth) == "cg10602367",]
Unmeth[rownames(Unmeth) == "cg10602367",]


###############################################################################
### STEP 05: quality checks
###############################################################################
## Team 1 detection p-value threshold: 0.05.

# --- a. QCplot ---
qc <- getQC(MSet.raw)
qc

# High median values for both signals indicate high-quality data.
# Some samples show comparatively lower methylated and unmethylated signals.
# Overall, the data can be kept for the next steps, but sample quality should
# be considered when interpreting downstream results.
png(file.path(figuresDir, "step05_qcplot.png"), width=7, height=6, units="in", res=300, pointsize=11)
plotQC(qc)
dev.off()

# Extracts all control probes from the array
# Count Control Types

NegativeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
table(NegativeControl$Type)
head(NegativeControl)
png(file.path(figuresDir, "step05_controlstripplot.png"), width=7, height=6, units="in", res=300, pointsize=11)
controlStripPlot(RGset, controls = "NEGATIVE")
dev.off()


# --- c. Detection p-values (Group 1 threshold = 0.05) ---
Detection_pvalue <- detectionP(RGset)
save(Detection_pvalue, file = file.path(rdsDir, "Detection_pvalue.RData"))

Over_Threshold <- Detection_pvalue > 0.05
head(Over_Threshold)
table(Over_Threshold)
sum(Over_Threshold)

Total_Failed_Positions <- sum(Over_Threshold)
Total_Detection_Tests <- length(Over_Threshold)
Total_Failed_Percentage <- Total_Failed_Positions / Total_Detection_Tests * 100

failed_probe_message <- sprintf(
  "There are %s failed detection tests with a detection p-value higher than 0.05",
  Total_Failed_Positions
)

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
print(failed)
print(failed_probe_message)




###############################################################################
### STEP 06: calculate raw beta and M values
###############################################################################

## getBeta and getM retrieve beta and M value matrices from the raw MethylSet.

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
saveRDS(beta, file=file.path(rdsDir, "beta_raw.rds"))
saveRDS(M, file=file.path(rdsDir, "M_raw.rds"))

### Recover phenotype data to split samples by group.
pheno <- data.frame(pData(RGset))
str(pheno)
table(pheno$Group)

### Split beta values into CTRL and DIS groups.
beta_CTRL <- beta[, pheno$Group=="CTRL"]
beta_DIS <- beta[, pheno$Group=="DIS"]

### Split M values into CTRL and DIS groups.
M_CTRL <- M[, pheno$Group=="CTRL"]
M_DIS <- M[, pheno$Group=="DIS"]

dim(beta_CTRL)
dim(beta_DIS)
dim(M_CTRL)
dim(M_DIS)

### Calculate mean beta and M values per CpG within each group.
mean_of_beta_CTRL <- apply(beta_CTRL, 1, mean, na.rm=TRUE)
mean_of_beta_DIS <- apply(beta_DIS, 1, mean, na.rm=TRUE)

mean_of_M_CTRL <- apply(M_CTRL, 1, mean, na.rm=TRUE)
mean_of_M_DIS <- apply(M_DIS, 1, mean, na.rm=TRUE)

### Remove infinite M values before density estimation.
mean_of_M_CTRL <- mean_of_M_CTRL[is.finite(mean_of_M_CTRL)]
mean_of_M_DIS <- mean_of_M_DIS[is.finite(mean_of_M_DIS)]

### Estimate density distributions for raw beta and M means.
d_mean_of_beta_CTRL <- density(mean_of_beta_CTRL, na.rm=TRUE)
d_mean_of_beta_DIS <- density(mean_of_beta_DIS, na.rm=TRUE)

d_mean_of_M_CTRL <- density(mean_of_M_CTRL, na.rm=TRUE)
d_mean_of_M_DIS <- density(mean_of_M_DIS, na.rm=TRUE)

### Plot raw beta density distributions by group.
png(file.path(figuresDir, "step06_raw_beta_density_CTRL_DIS.png"),
    width=7, height=6, units="in", res=300, pointsize=11)

plot(d_mean_of_beta_CTRL, col="blue",
     main="Raw beta values", xlab="Mean beta")

lines(d_mean_of_beta_DIS, col="red")

legend("topright", legend=c("CTRL", "DIS"),
       col=c("blue", "red"), lty=1, bty="n")

dev.off()

### Plot raw M value density distributions by group.
png(file.path(figuresDir, "step06_raw_M_density_CTRL_DIS.png"),
    width=7, height=6, units="in", res=300, pointsize=11)

plot(d_mean_of_M_CTRL, col="blue",
     main="Raw M values", xlab="Mean M")

lines(d_mean_of_M_DIS, col="red")

legend("topright", legend=c("CTRL", "DIS"),
       col=c("blue", "red"), lty=1, bty="n")

dev.off()

### Save summary vectors used for the density plots.
saveRDS(mean_of_beta_CTRL,
        file=file.path(rdsDir, "mean_of_beta_raw_CTRL.rds"))

saveRDS(mean_of_beta_DIS,
        file=file.path(rdsDir, "mean_of_beta_raw_DIS.rds"))

saveRDS(mean_of_M_CTRL,
        file=file.path(rdsDir, "mean_of_M_raw_CTRL.rds"))

saveRDS(mean_of_M_DIS,
        file=file.path(rdsDir, "mean_of_M_raw_DIS.rds"))

###############################################################################
### STEP 07: normalization using preprocessNoob
###############################################################################
## Team 1 normalization method: preprocessNoob.

### Open the raw RGChannelSet and raw beta matrix saved in previous steps.
load(file.path(rdsDir, "RGset_Report.RData"))
beta <- readRDS(file.path(rdsDir, "beta_raw.rds"))

### Inspect the raw RGChannelSet and beta matrix.
RGset
head(beta)
dim(beta)

### Retrieve Type I and Type II probe annotations from the manifest.
dfI <- data.frame(getProbeInfo(RGset, type="I"))
dfII <- data.frame(getProbeInfo(RGset, type="II"))

dim(dfI)
dim(dfII)
str(dfI)
str(dfII)
head(dfI$Name)
head(dfII$Name)

### Split raw beta values by probe chemistry.
beta_I <- beta[rownames(beta) %in% dfI$Name,]
beta_II <- beta[rownames(beta) %in% dfII$Name,]

dim(beta_I)
dim(beta_II)

### Calculate mean raw beta by probe chemistry.
mean_of_beta_I <- apply(beta_I, 1, mean, na.rm=T)
mean_of_beta_II <- apply(beta_II, 1, mean, na.rm=T)

### Estimate density of raw beta means by probe chemistry.
d_mean_of_beta_I <- density(mean_of_beta_I, na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II, na.rm=T)

### Calculate and estimate density of raw beta standard deviations.
sd_of_beta_I <- apply(beta_I, 1, sd, na.rm=T)
sd_of_beta_II <- apply(beta_II, 1, sd, na.rm=T)

d_sd_of_beta_I <- density(sd_of_beta_I, na.rm=T)
d_sd_of_beta_II <- density(sd_of_beta_II, na.rm=T)

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
saveRDS(preprocessNoob_results, file=file.path(rdsDir, "preprocessNoob_results.rds"))
saveRDS(beta_preprocessNoob, file=file.path(rdsDir, "beta_preprocessNoob.rds"))
saveRDS(M_preprocessNoob, file=file.path(rdsDir, "M_preprocessNoob.rds"))

### Split normalized beta values by probe chemistry.
beta_preprocessNoob_I <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfI$Name,]
beta_preprocessNoob_II <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfII$Name,]

### Calculate mean normalized beta values by probe chemistry.
mean_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, mean, na.rm=T)
mean_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, mean, na.rm=T)

### Estimate density of normalized beta means.
d_mean_of_beta_preprocessNoob_I <- density(mean_of_beta_preprocessNoob_I, na.rm=T)
d_mean_of_beta_preprocessNoob_II <- density(mean_of_beta_preprocessNoob_II, na.rm=T)

### Calculate and estimate density of normalized beta standard deviations.
sd_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, sd, na.rm=T)
sd_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, sd, na.rm=T)

d_sd_of_beta_preprocessNoob_I <- density(sd_of_beta_preprocessNoob_I, na.rm=T)
d_sd_of_beta_preprocessNoob_II <- density(sd_of_beta_preprocessNoob_II, na.rm=T)

### Define colors for the raw and normalized beta boxplots.
pheno <- data.frame(pData(RGset))
pheno$Group

boxplot_col <- c(CTRL = "blue", DIS = "red")[as.character(pheno$Group)]
boxplot_col

### Compare raw and normalized beta values using plots.
png(file.path(figuresDir, "step07_raw_vs_preprocessNoob_6panel.png"), width=13, height=7, units="in", res=300, pointsize=10)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), oma = c(1, 1, 1, 5))

plot(d_mean_of_beta_I, col = "blue", main = "raw beta", xlim = c(0, 1), ylim = c(0, 5))
lines(d_mean_of_beta_II, col="red")
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, bty = "n", cex = 0.8)

plot(d_sd_of_beta_I, col="blue", main="raw sd", xlim=c(0,0.6), ylim=c(0,60))
lines(d_sd_of_beta_II, col = "red")
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, bty = "n", cex = 0.8)

boxplot(beta, ylim = c(0, 1), col = boxplot_col, main = "raw beta", names=FALSE)
axis( 1, at = seq_len(ncol(beta)), labels = seq_len(ncol(beta)), cex.axis = 0.8)
legend( "right", inset = c(-0.4, 0), legend = c("CTRL", "DIS"), fill = c("blue", "red"), bty = "n", cex = 0.8, xpd = NA)

plot(d_mean_of_beta_preprocessNoob_I, col="blue", main="preprocessNoob beta", xlim=c(0,1), ylim=c(0,5))
lines(d_mean_of_beta_preprocessNoob_II, col = "red")
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, bty = "n", cex = 0.8)

plot(d_sd_of_beta_preprocessNoob_I, col="blue", main="preprocessNoob sd", xlim=c(0,0.6), ylim=c(0,60))
lines(d_sd_of_beta_preprocessNoob_II, col = "red")
legend("topright", legend = c("Type I", "Type II"), col = c("blue", "red"), lty = 1, bty = "n", cex = 0.8)

boxplot(beta_preprocessNoob, ylim = c(0, 1), col = boxplot_col, main = "preprocessNoob beta", names=FALSE)
axis( 1, at = seq_len(ncol(beta_preprocessNoob)), labels = seq_len(ncol(beta_preprocessNoob)), cex.axis = 0.8)
legend( "right", inset = c(-0.4, 0), legend = c("CTRL", "DIS"), fill = c("blue", "red"), bty = "n", cex = 0.8, xpd = NA)
dev.off()


###############################################################################
### STEP 08: PCA on normalized beta values
###############################################################################

# Load normalized beta values if they are not already in memory
beta_preprocessNoob <- readRDS("results/rds/beta_preprocessNoob.rds")
load("results/rds/RGset_Report.RData")

# Extract sample information from RGset
sample_info <- data.frame(pData(RGset))

# Run PCA on normalized beta values
pca_results <- prcomp(t(beta_preprocessNoob), scale. = TRUE)

# Calculate the proportion of variance explained by each component
pca_var <- pca_results$sdev^2
pca_var_explained <- pca_var / sum(pca_var)

## Create Scree Plot
png(
  "results/figures/step08_PCA_ScreePlot.png",
  width = 1200,
  height = 800
)

bp <- barplot(
  pca_var_explained[1:8],
  names.arg = paste0("PC", 1:8),
  ylim = c(0, 0.5),
  main = "Variance Explained by Principal Components",
  xlab = "Principal Components",
  ylab = "Proportion of Variance Explained"
)

text(
  bp,
  pca_var_explained[1:8],
  labels = paste0(round(pca_var_explained[1:8] * 100, 1), "%"),
  pos = 3,
  cex = 0.8
)

dev.off()

## Create PCA plot colored by disease group
png(
  "results/figures/step08_PCA_Group.png",
  width = 1200,
  height = 800
)

plot(
  pca_results$x[,1],
  pca_results$x[,2],
  col = as.factor(sample_info$Group),
  pch = 17,
  cex = 1.5,
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
  pca_results$x[,1],
  pca_results$x[,2],
  labels = sample_info$SampleID,
  pos = 3,
  cex = 0.8
)

legend(
  "topright",
  legend = levels(as.factor(sample_info$Group)),
  col = 1:length(levels(as.factor(sample_info$Group))),
  pch = 17
)

dev.off()

## PCA colored by sex to check whether samples cluster according to gender
sample_info$Sex <- as.factor(sample_info$Sex)

png(
  "results/figures/step08_PCA_by_Sex.png",
  width = 1200,
  height = 800
)

plot(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  col = sample_info$Sex,
  pch = 17,
  cex = 1.5,
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

## PCA colored by Sentrix_ID (batch)

sample_info$Slide <- as.factor(sample_info$Slide)

png(
  "results/figures/step08_PCA_by_SentrixID.png",
  width = 1200,
  height = 800
)

plot(
  pca_results$x[, "PC1"],
  pca_results$x[, "PC2"],
  col = sample_info$Slide,
  pch = 17,
  cex = 1.5,
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
  legend = levels(sample_info$Slide),
  col = seq_along(levels(sample_info$Slide)),
  pch = 17
)

dev.off()

###############################################################################
### STEP 09: differential methylation analysis using t-test
###############################################################################
## Team 1 differential methylation test: t-test.

# Load the RGset object generated in Step 1
load(
  "results/rds/RGset_Report.RData"
)

# Load the normalized beta values generated in Step 7
beta_preprocessNoob <- readRDS(
  "results/rds/beta_preprocessNoob.rds"
)

# Extract sample information from RGset
sample_info <- data.frame(
  pData(RGset)
)

# Use simpler names for the analysis
beta_matrix <- beta_preprocessNoob
sample_data <- sample_info

# Identify CTRL and DIS samples
control_samples <- sample_data$Group == "CTRL"
disease_samples <- sample_data$Group == "DIS"

# Calculate the mean beta value for each probe in the CTRL group
mean_ctrl <- rowMeans(
  beta_matrix[, control_samples]
)

# Calculate the mean beta value for each probe in the DIS group
mean_dis <- rowMeans(
  beta_matrix[, disease_samples]
)

# Calculate the methylation difference between the two groups
delta_beta <- mean_dis - mean_ctrl

# Perform a t-test for each probe and collect the p-values
p_values <- apply(
  beta_matrix,
  1,
  function(probe_values) {
    t.test(
      probe_values[control_samples],
      probe_values[disease_samples])$p.value}
)

# Create a table containing the differential methylation results
t_test_results <- data.frame(
  ProbeID = rownames(beta_matrix),
  Mean_CTRL = mean_ctrl,
  Mean_DIS = mean_dis,
  Delta_Beta = delta_beta,
  P_Value = p_values
)

# Save the results for the next step of the pipeline
write.csv(
  t_test_results,
  "results/tables/step09_t_test_results.csv",
  row.names = FALSE
)

###############################################################################
### STEP 10: Multiple testing correction
###############################################################################

cat("\n--- Running Step 10: Multiple Testing Correction ---\n")

# Load the t-test results generated in Step 09
t_test_results <- read.csv("results/tables/step09_t_test_results.csv")

# Apply Bonferroni correction (Family-Wise Error Rate)
t_test_results$P_Bonferroni <- p.adjust(t_test_results$P_Value, method = "bonferroni")

# Apply Benjamini-Hochberg correction (False Discovery Rate)
t_test_results$P_FDR <- p.adjust(t_test_results$P_Value, method = "BH")

# Set the significance threshold
alpha <- 0.05

# Count significant probes
sig_nominal <- sum(t_test_results$P_Value <= alpha, na.rm = TRUE)
sig_bonferroni <- sum(t_test_results$P_Bonferroni <= alpha, na.rm = TRUE)
sig_fdr <- sum(t_test_results$P_FDR <= alpha, na.rm = TRUE)

cat("Significant probes (Nominal P <= 0.05):", sig_nominal, "\n")
cat("Significant probes (Bonferroni <= 0.05):", sig_bonferroni, "\n")
cat("Significant probes (BH/FDR <= 0.05):", sig_fdr, "\n")

# Save the adjusted table
write.csv(t_test_results, "results/tables/step10_t_test_results_adjusted.csv", row.names = FALSE)


###############################################################################
### STEP 11: Volcano and Manhattan plots
###############################################################################

cat("\n--- Running Step 11: Volcano and Manhattan Plots ---\n")

# --- A. Volcano Plot ---
png("results/figures/step11_Volcano_Plot.png", width=8, height=6, units="in", res=300)

plot(
  x = t_test_results$Delta_Beta,
  y = -log10(t_test_results$P_Value),
  pch = 16,
  cex = 0.5,
  col = "darkgray",
  xlab = "Delta Beta (DIS - CTRL)",
  ylab = "-log10(Nominal P-value)",
  main = "Volcano Plot of DNA Methylation Differences",
  xlim = c(-max(abs(t_test_results$Delta_Beta), na.rm=TRUE), max(abs(t_test_results$Delta_Beta), na.rm=TRUE))
)

# Highlight significant targets (|Delta Beta| > 0.1 & p <= 0.05)
hyper <- which(t_test_results$Delta_Beta > 0.1 & t_test_results$P_Value <= alpha)
hypo <- which(t_test_results$Delta_Beta < -0.1 & t_test_results$P_Value <= alpha)

points(t_test_results$Delta_Beta[hyper], -log10(t_test_results$P_Value)[hyper], col = "red", pch = 16, cex = 0.6)
points(t_test_results$Delta_Beta[hypo], -log10(t_test_results$P_Value)[hypo], col = "blue", pch = 16, cex = 0.6)

abline(h = -log10(alpha), col = "black", lty = 2)
abline(v = c(-0.1, 0.1), col = "black", lty = 2)

legend("topright", legend = c("Hypermethylated", "Hypomethylated", "Not Significant"), 
       col = c("red", "blue", "darkgray"), pch = 16, cex = 0.8, bty="n")
dev.off()


# --- B. Manhattan Plot ---
library(qqman)
preprocessNoob_results <- readRDS("results/rds/preprocessNoob_results.rds")
ann <- getAnnotation(preprocessNoob_results)
locations <- data.frame(ProbeID = rownames(ann), CHR = ann$chr, MAPINFO = ann$pos)

# Merge coordinates with t-test results
manhattan_data <- merge(t_test_results, locations, by = "ProbeID")

# Format chromosomes for qqman
manhattan_data$CHR <- gsub("chr", "", manhattan_data$CHR)
manhattan_data$CHR[manhattan_data$CHR == "X"] <- 23
manhattan_data$CHR[manhattan_data$CHR == "Y"] <- 24
manhattan_data$CHR <- as.numeric(manhattan_data$CHR)
manhattan_data <- manhattan_data[!is.na(manhattan_data$CHR) & !is.na(manhattan_data$P_Value), ]

png("results/figures/step11_Manhattan_Plot.png", width=10, height=6, units="in", res=300)
manhattan(
  manhattan_data,
  chr = "CHR",
  bp = "MAPINFO",
  snp = "ProbeID",
  p = "P_Value",
  col = c("gray50", "gray80"),
  suggestiveline = FALSE, 
  genomewideline = FALSE, 
  main = "Manhattan Plot of Epigenome-Wide Association",
  ylab = "-log10(Nominal P-value)"
)
dev.off()


###############################################################################
### STEP 12: Heatmap of the top 100 CpG probes
###############################################################################

cat("\n--- Running Step 12: Heatmap of Top 100 CpGs ---\n")
library(gplots)

# Extract top 100 probes by nominal P-value
top_100_results <- t_test_results[order(t_test_results$P_Value), ][1:100, ]
top_100_probes <- top_100_results$ProbeID

# Subset normalized matrix
beta_preprocessNoob <- readRDS("results/rds/beta_preprocessNoob.rds")
beta_top100 <- beta_preprocessNoob[top_100_probes, ]

# Setup phenotype colors
load("results/rds/RGset_Report.RData")
pheno <- data.frame(pData(RGset))
group_colors <- ifelse(pheno$Group == "CTRL", "blue", "red")

png("results/figures/step12_Heatmap_Top100.png", width=8, height=8, units="in", res=300)
heatmap.2(
  as.matrix(beta_top100),
  main = "Top 100 Differentially Methylated CpGs",
  trace = "none",              
  col = bluered(100),          
  scale = "none",              
  dendrogram = "both",         
  ColSideColors = group_colors,
  margins = c(10, 5),          
  cexCol = 0.8,
  labRow = FALSE               
)

legend("topright", 
       legend = levels(as.factor(pheno$Group)), 
       fill = c("blue", "red"), 
       border = FALSE, 
       bty = "n", 
       cex = 0.8)
dev.off()

cat("\n--- Pipeline Execution Complete ---\n")