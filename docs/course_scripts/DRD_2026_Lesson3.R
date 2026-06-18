rm(list=ls())
setwd("D:/Google Drive/UNIBO/Materials/DNA & RNA Dynamics/Module 2/DRD_2026_Module2_Lesson3")
library(minfi)
library(minfiData)
vignette("minfi")

## CHUNK 1: Import raw data
# List the files included in the folder
list.files("Input_Data/")
# Let's have a look at the csv
SampleSheet <- read.table("Input_Data/SampleSheet.csv",sep=",",header=T)
SampleSheet
# or
SampleSheet <- read.csv("Input_Data/SampleSheet.csv",header=T)
SampleSheet

# Set the directory in which the raw data are stored and load the samplesheet using the function read.metharray.sheet
baseDir <- ("Input_data")
targets <- read.metharray.sheet(baseDir)
targets

# Create an object of class RGChannelSet using the function read.metharray.exp
?read.metharray.exp
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")

# Let's explore the RGset object:
RGset
str(RGset)
?RGChannelSet

# We extract the Green and Red Channels using the functions getGreen and getRed
Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)

Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

## CHUNK 2: Addresses, probes and out of band probes
load('../DRD_2026_Module2_Lesson2/Illumina450Manifest_clean.RData')
ls()
head(Illumina450Manifest_clean)

# From Addresses to probes: 

# I want to check the probe having address 10600313 (The first rowname in Red and Green objects)
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="10600313",]
# there is not an AddressA_ID with code 10600313...
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="10600313",]
# ...but there is an AddressB_ID with code 10600313, and it is associated to the probe cg25192902, type I, Red; the Address_A of this probe is 10679328

# I want to check the probe having address 10600322 (The second rowname in Red and Green objects)
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="10600322",]
# There is an AddressA_ID with code 10600322, and it is associated to the probe cg00226849, type II. Note that the column Color_Channel is empty, as all the Type II probes will emit in Red or Green depending on the methylation status of the target CpG. In addition, also Address_B is empty, as type II probes never have an Address_B. Indeed, if you digit
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="10600322",]
# the resulting dataframe will be empty. 


# From probes to Addresses
head(Illumina450Manifest_clean)
# cg00035864, type II (the first probe in the Illumina450Manifest_clean object)
Red[rownames(Red)=="31729416",]
Green[rownames(Green)=="31729416",]

# cg00050873, type I Red (the second probe in the Illumina450Manifest_clean object)
Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg00050873",]
Red[rownames(Red)=="32735311",]
Red[rownames(Red)=="31717405",]
# But the two addresses are present also when I look at the Green object: these are out of band signals!
Green[rownames(Green)=="32735311",]
Green[rownames(Green)=="31717405",]

# Exercise 1: Let's consider the probe cg03360716. Is it a type I or type II probe? What is/are its addresses? What are the associated red and green fluorescences?
#Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg03360716",c('Infinium_Design_Type')]

## CHUNK 3: Class "IlluminaMethylationManifest"
# The RGChannelSet stores also a manifest object that contains the probe design information of the array (IlluminaMethylationManifest).
?IlluminaMethylationManifest

# From the help page, you see that getManifest and getManifestInfo are two accessor functions for this class. Let's check their results:
getManifest(RGset)
getManifestInfo(RGset)
# getManifestInfo does not seem really useful! :-D

# What about the getProbeInfo() function? It returns a data.frame giving the type I, type II or control probes
getProbeInfo(RGset)
# Only the first and the last rows of the object are printed. We will create an object df containing the result of this function:
df <- data.frame(getProbeInfo(RGset))
dim(df)
# Strange...this df has 135476 rows, not 485512...why?
head(getProbeInfo(RGset, type = "II"))
df_II <- data.frame(getProbeInfo(RGset, type = "II"))
dim(df_II)
350036+135476 #485512
# This means that in the default settings, the getProbeInfo function returns only Type I probes. However, this is not really specified in the help page. Take home message: when you work with a package, try before you trust!


## CHUNK 4: Extract methylated and unmethylated signals
# We will use the function MSet.raw
MSet.raw <- preprocessRaw(RGset)
MSet.raw
?MethylSet
# Note that now the number of rows is 485512, exactly the number of probes according to the Manifest!
save(MSet.raw,file="MSet_raw.RData")
Meth <- as.matrix(getMeth(MSet.raw))
str(Meth)
head(Meth)
Unmeth <- as.matrix(getUnmeth(MSet.raw))
str(Unmeth)
head(Unmeth)

# Let's check what happens to the probes that we considered before when we move from RGset to MethylSet
# cg00035864, type II (the first probe in the Illumina450Manifest_clean object)
Red[rownames(Red)=="31729416",]
Green[rownames(Green)=="31729416",]
Unmeth[rownames(Unmeth)=="cg00035864",]
Meth[rownames(Meth)=="cg00035864",]

# cg00050873, type I Red (the second probe in the Illumina450Manifest_clean object)
Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg00050873",]
Red[rownames(Red)=="32735311",]
Red[rownames(Red)=="31717405",]
# But the two addresses are present also when I look at the Green object: these are out of band signals!
Green[rownames(Green)=="32735311",]
Green[rownames(Green)=="31717405",]
Unmeth[rownames(Unmeth)=="cg00050873",]
Meth[rownames(Meth)=="cg00050873",]

## CHUNK 5: QC plot and control probes
# First of all, we consider the median of methylation and unmethylation channels for each sample, using the function getQC 

qc <- getQC(MSet.raw)
qc
plotQC(qc)

# All the samples have high median methylation and unmethylation signals, i.e. they have good quality!

# Then, we check the control probes. The Illumina guide (GenomeStudio_control_probes.pdf) tells us the expected intensities for each type of control probes.

# To know the names and the numbers of the different types of control probes, we apply the getProbeInfo function to the RGset object.
getProbeInfo(RGset, type = "Control")
df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
str(df_TypeControl)
table(df_TypeControl$Type)

# Then, we can use the controlSripPlot function to plot the intensity values of each type of controls probes in our samples

controlStripPlot(RGset, controls="NEGATIVE")
# Negative controls are all fine, as all below 1000 (log2(1000)=10)
# Note that there is a minor error in the package...green and red labels are swapped (only in R4.0.0 or higher)
# You can compare the results of these controls with the plots in Lesson 3 slides...everything looks fine!
controlStripPlot(RGset, controls="EXTENSION")
controlStripPlot(RGset, controls="STAINING")
controlStripPlot(RGset, controls="BISULFITE CONVERSION I")
controlStripPlot(RGset, controls="BISULFITE CONVERSION II")
controlStripPlot(RGset, controls="SPECIFICITY I")
controlStripPlot(RGset, controls="SPECIFICITY II")

## CHUNK 6: Detection p-value and filtering
?detectionP
# The detectionP() functions takes a bit. You can try to jump the next 2 lines and load the object
detP <- detectionP(RGset) # do not launch it now, it takes some minutes
save(detP,file="detP.RData")

load("detP.RData") # or detP <- readRDS("detP.rds")
str(detP)
dim(detP)
head(detP)
# Note that in the detP object, rownames are the CpG probes names

# Now we want to consider a detectionP threshold of 0.05
failed <- detP>0.05
head(failed)
dim(failed)
table(failed)
# summary(failed) is particularly useful, as it returns the number of failed (TRUE) and not failed (FALSE) positions for each sample
summary(failed) 

# Unfortunately, minfi does not provide any function that performs filtering of samples/probes on the basis of detection p-value! Other packages, like wateRmelon, allow to perform this action using convenient functions.

# In any case, we can manually operate on the detP object to filter out failed samples and probes

# First of all, we want to calculate the percentage of failed positions in each sample, and filter out those samples that have a % of failed positions higher than a certain threshold, for example 5%. To calculate the fraction of failed postions per sample, we can use the function colMeans() and apply it to the "failed" object that we have previoulsy created
?colMeans
means_of_columns <- colMeans(failed) 
head(means_of_columns)
# Each value corresponds to the ratio (number of TRUE)/(number of TRUE+FALSE) for each column (that is, for each sample)

# For example, for sample 9344737127_R01C02 we had 499 TRUE and 485013 FALSE --> 499/(499+485013)= 0.001027781

# We can decide, for example, to retain only the samples in which <5 % of sites has a bad detection p-value (greater than 0.05, the threshold that we used to calculate the object "failed")
samples_to_be_retained <- means_of_columns<0.05
samples_to_be_retained
# Then I subset the RGset object in order to retain only the selected samples:
RGset <- RGset[samples_to_be_retained]
RGset
# In this example, all the samples have < 5% of probes with a bad detection pValue, and therefore I can retain all the samples


# The next step is to identify the probes that have a bad detectionP. To calculate the fraction of failed samples per probes, we can use the function rowMeans() and apply it to the "failed" object that we have previoulsy created
means_of_rows <- rowMeans(failed)
head(means_of_rows)
# Each value corresponds to the ratio (number of TRUE)/(number of TRUE+FALSE) for each row (that is, for each probe)
head(failed)
# For example, for the probe cg00050873 we had 6 TRUE and 2 FALSE --> 6/(6+2)= 0.750

# We can decide, for example, to retain only the probes in which <1 % of samples has a bad detection p-value (greater than 0.05, the threshold that we used to calculate the object "failed")
probes_to_be_retained <- means_of_rows<0.01
head(probes_to_be_retained)
str(probes_to_be_retained)
# Note that the "probes_to_be_retained" object is a vector of logical vectors with names (that is, each element of the vector has a name). We can access to these names using the function names()
names(probes_to_be_retained[1:10])
# returns the names of the first 10 probes
table(probes_to_be_retained)
# counts the number of TRUE and FALSE in the vector.

# Now I want to select the probes that I want to remove: these are the probes with "FALSE" in the probes_to_be_retained vector
probes_to_be_removed <- probes_to_be_retained[probes_to_be_retained==FALSE]
length(probes_to_be_removed)
names_probes_to_be_removed <- names(probes_to_be_removed)
head(names_probes_to_be_removed)

# Unfortunately, as we will see in the next lesson, most of the normalization procedures implemented in minfi work on the RGset object and do not like missing probes. If the number of bad probes is low, they will not affect the normalization procedures. Therefore, we can retain the bad probes at this step, but save them in order to remove them after the normalizaion step.
save(names_probes_to_be_removed,file="names_probes_to_be_removed.RData")

# Exercise 2: load the RGset_CRC.RData or RGset_CRC.rds object (this is from an IlluminaHumanMethylationEPIC experiment and you will nee to install the IlluminaHumanMethylationEPICmanifest from Bioconductor). How many samples does it contain? Use plotQC the check the quality of the samples and check NEGATIVE and EXTENSION controls...do you see something wrong with these data? 
