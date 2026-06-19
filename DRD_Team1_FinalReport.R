## Team 1 assignments:
## - Step 3 address: 42796479
## - Step 5 detection p-value threshold: 0.05
## - Step 7 normalization method: preprocessNoob
## - Step 9 differential methylation test: t-test

# Step 00. Setup

### Clean the workspace before starting the analysis.
rm(list=ls())


### Load all packages used by the complete workflow.

BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")


library(minfi)
library(qqman)
library(gplots)
library(IlluminaHumanMethylation450kmanifest)

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

list.files("./data/raw/")
SampleSheet <- read.csv("./data/raw/SampleSheet_Report_II.csv", header = TRUE)
View(SampleSheet)

targets <- read.metharray.sheet("./data/raw")
View(targets)

RGset <- read.metharray.exp(targets = targets)

save(RGset, file = "RGset_Report.RData")
load("RGset_Report.RData")

RGset
str(RGset)


###############################################################################
### STEP 02: extract Red and Green fluorescence values
###############################################################################

Red <- data.frame(getRed(RGset))
dim(Red)
View(Red)
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
View(ProbeInfo_I)

ProbeInfo_II <- data.frame(getProbeInfo(RGset, type = "II"))
View(ProbeInfo_II)

ProbeInfo_I[ProbeInfo_I$AddressA == TARGET_ADDRESS, ]
ProbeInfo_I[ProbeInfo_I$AddressB == TARGET_ADDRESS, ]
ProbeInfo_II[ProbeInfo_II$AddressA == TARGET_ADDRESS, ]
# Result: Type I, AddressA → UNMETHYLATED probe of cg10602367, Color = Red

# --- Extract fluorescence and build summary ---
Red_fluorescences   <- Red[rownames(Red) == TARGET_ADDRESS, ]
View(Red_fluorescences)

Green_fluorescences <- Green[rownames(Green) == TARGET_ADDRESS, ]
View((Green_fluorescences))

df_summary_probes <- data.frame(
  Sample            = colnames(Red_fluorescences),
  Red_Fluorescence  = as.numeric(Red_fluorescences[1, ]),
  Green_Fluorescence= as.numeric(Green_fluorescences[1, ]),
  Type              = "I",
  Color             = "Red"
)

print(df_summary_probes)


# ============================================================
# STEP 4. Creation of the MSet.raw object
# ============================================================

# Convert RGset to MethylSet (raw, no normalization)
MSet.raw <- preprocessRaw(RGset)
MSet.raw


# Save for later use
save(MSet.raw, file = "MSet_raw.RData")

# Extract methylated and unmethylated signal matrices
Meth <- as.matrix(getMeth(MSet.raw))
str(Meth)
head(Meth)
View(Meth)

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
# In this case, the samples are spread apart, and most show low 
# median values for both signals, suggesting poor quality.
# Out of the eight samples, only three are considered to be of good quality
plotQC(qc)

# Extracts all control probes from the array
# Count Control Types
# 
NegativeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
table(NegativeControl$Type)
head(NegativeControl)
controlStripPlot(RGset, controls = "NEGATIVE")


# --- c. Detection p-values (Group 1 threshold = 0.05) ---
Detection_pvalue <- detectionP(RGset)
save(Detection_pvalue, file = "Detection_pvalue.RData")

Over_Threshold <- Detection_pvalue > 0.05
head(Over_Threshold)
table(Over_Threshold)
sum(Over_Threshold)

failed_probe_message <- sprintf(
  "There are %s failed probes with a detection p-value higher than 0.05",
  sum(Over_Threshold)
)

Over_Threshold_Per_Sample <- colSums(Over_Threshold)
Over_Threshold_Per_Sample

summary(Over_Threshold)

failed <- data.frame(
  Sample = names(Over_Threshold_Per_Sample),
  n_Failed_Positions = Over_Threshold_Per_Sample
)
print(failed)
print(failed_probe_message)




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



