## The Illumina 450k manifest
# First of all, let's clean the workspace and set our the working directory:
rm(list=ls())
setwd("D:/Google Drive/UNIBO/Materials/DNA & RNA Dynamics/Module 2/DRD_2026_Module2_Lesson2")

# CHUNK1: upload the manifest and explore it ####
# The Illumina 450k manifest has been downloaded from this link:
# http://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html

# The csv file is quite large (200 MB). It is a good idea to import in R just the first 20 rows of the file to check that everything is ok. To this aim, we set the nrows argument=20. As this file is a csv, we will use the comma as separator.
Illumina450Manifest <- read.table("HumanMethylation450_15017482_v1-2.csv",sep=",",header=T, nrows=20)
# There are some problems :-( 
# Let's try to use the readLines() function to check the first lines of the file. You will see that the problem is that the first 7 lines of the file are a sort of header, and that the spreadsheet actually starts from the 8th line. Now we try to import again the file skipping the first 7 lines (argument skip=7). Furthermore, we will set stringsAsFactors = T
readLines("HumanMethylation450_15017482_v1-2.csv",n=20)
Illumina450Manifest <- read.table("HumanMethylation450_15017482_v1-2.csv",sep=",",header=T, nrows=20, skip=7, stringsAsFactors = T) 
dim(Illumina450Manifest)
head(Illumina450Manifest)
str(Illumina450Manifest)

# You can also use the read.csv function, which uses the comma as default separator
?read.csv
Illumina450Manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv",header=T, nrows=20,skip=7, stringsAsFactors = T)
head(Illumina450Manifest)
str(Illumina450Manifest)

# The descriptions of the Manifest Column Headings can be found here: https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf. In summary: 
# IlmnID: Unique CpG locus identifier from the Illumina CG database
# Name: Unique CpG locus identifier from the Illumina CG database
# AddressA_ID: Address of probe A
# AlleleA_ProbeSeq: Sequence for probe A
# AddressB_ID: Address of probe  B - Infinium I assays only
# AlleleB_ProbeSeq: Sequence for probe B - Infinium I assays only
# Infinium_Design_Type: Defines Assay type - Infinium I or Infinium II
# Next_Base: Base added at SBE step - Infinium I assays only
# Color_Channel: Color of the incorporated base  (Red or Green) - Infinium I assays only
# Forward_Sequence: Sequence (in 5'-3' orientation) flanking query site
# Genome_Build: Genome build on which forward sequence is based
# CHR: Chromosome - genome build 37
# MAPINFO: Coordinates - genome build 37
# SourceSeq: Unconverted design sequence
# Chromosome_36: Chromosome - genome build 36
# Coordinate_36: Coordinates - genome build 36
# Strand: Design strand
# Probe_SNPs: Assays with SNPs present within probe >10bp from query site
# Probe_SNPs_10: Assays with SNPs present within probe ≤10bp from query site (HM27 carryover or recently discovered)
# Random_Loci: Loci which were chosen randomly in the design proccess
# Methyl27_Loci: Present or absent on HumanMethylation27 array
# UCSC_RefGene_Name: Gene name (UCSC)
# UCSC_RefGene_Accession: Accession number (UCSC)
# UCSC_RefGene_Group: Gene region feature category (UCSC)
# UCSC_CpG_Islands_Name: CpG island name (UCSC)
# Relation_to_UCSC_CpG_Island: Relationship to Canonical CpG Island: Shores - 0-2 kb from CpG island; Shelves - 2-4 kb from CpG island.
# Phantom: FANTOM-derived promoter
# DMR: Differentially methylated region (experimentally determined)
# Enhancer: Enhancer element (informatically-determined)
# HMM_Island: Hidden Markov Model Island
# Regulatory_Feature_Name: Regulatory feature (informatically determined)
# Regulatory_Feature_Group: Regulatory feature category
# DHS: DNAse hypersensitive site (experimentally determined)

# Now we can upload the whole manifest. This step can take some minutes, if you encounter problems you can skip the next two command lines and directly go to "Load the Illumina450Manifest.RData":
Illumina450Manifest <- read.csv("HumanMethylation450_15017482_v1-2.csv",header=T, skip=7, stringsAsFactors = T) 
save(Illumina450Manifest,file="Illumina450Manifest.RData") 
# Remember that if you save your data with save(), it cannot be restored under different name. The original object name will be automatically used when you will load the data.

# Load the Illumina450Manifest.RData
load("Illumina450Manifest.RData")
# To know the name of the object that we have just uploaded, we can list the objects in the workspace:
ls()
# Let's apply some exploratory functions:
dim(Illumina450Manifest)
head(Illumina450Manifest)
str(Illumina450Manifest)
colnames(Illumina450Manifest)
head(rownames(Illumina450Manifest))

# Let's check the distribution of probes across the chromosomes
table(Illumina450Manifest$CHR)
# Remember that by default,factor levels are ordered according to alphabetical order (1,10,11, etc).

# Remember that, to select the column CHR, you can write also:
table(Illumina450Manifest[,"CHR"])
# or
table(Illumina450Manifest[,]$CHR)
# or index it. From the result of the colnames function, we know that the column CHR is the 12th column, so you can write:
table(Illumina450Manifest[,12])


# CHUNK 2: create a "clean" version of the manifest ####
# There are 916 probes that do not map in any chromosome. Let's have a look at them: we will subset the dataframe in order to retain only the rows for which the column CHR is empty (nothing between the quotes)
notMappedToCHR <- Illumina450Manifest[Illumina450Manifest$CHR=="",]
dim(notMappedToCHR)
str(notMappedToCHR)
# You will see that R remembers all the factor levels from the original dataframe. Use the droplevels function to remove the unused levels!
notMappedToCHR <- droplevels(notMappedToCHR)
dim(notMappedToCHR)
str(notMappedToCHR)

# 916 lines are not too much...let's try to print them, considering only the first 2 columns:
notMappedToCHR[,c(1,2)]
# We can use the table function to have some hints:
table(notMappedToCHR$Name) #Try the other possibilities to select the column "Name"
# There are 65 probes whose name starts with "rs". These probes do not measure DNA methylation, but genetic polymorphisms (SNPs). They are included in the 450k design to allow sample tracking and check for sample mixups. The remaining probes are control probes that we will describe in the next lesson

# Now we want to create a new dataframe in which we remove the control and the rs probes. We can perform this task in different ways:
# we can create an object that stores the result of the condition:
only_in_a_chromosome <- Illumina450Manifest$CHR!=""
head(only_in_a_chromosome)
table(only_in_a_chromosome)
Illumina450Manifest_clean <- Illumina450Manifest[only_in_a_chromosome,]
dim(Illumina450Manifest_clean)

# or we can apply directly tha condition to subset the dataframe
Illumina450Manifest_clean <- Illumina450Manifest[Illumina450Manifest$CHR!="",]
dim(Illumina450Manifest_clean)

# or we can select the probes which are not included in the notMappedToCHR object
Illumina450Manifest_clean <- Illumina450Manifest[!Illumina450Manifest$IlmnID %in% notMappedToCHR$IlmnID,]
dim(Illumina450Manifest_clean)

# In any case, always remember to drop the unused levels!
str(Illumina450Manifest_clean)
Illumina450Manifest_clean <- droplevels(Illumina450Manifest_clean)
str(Illumina450Manifest_clean)
# We will save the clean manifest in a new file:
save(Illumina450Manifest_clean,file="Illumina450Manifest_clean.RData")

# CHUNK3: Islands, shores and shelves ####
# Let's use the exploratory function table: 
table(Illumina450Manifest_clean$Relation_to_UCSC_CpG_Island)

# How many probes map in Islands?
ProbesIslands <- Illumina450Manifest_clean[Illumina450Manifest_clean$Relation_to_UCSC_CpG_Island=="Island",]
dim(ProbesIslands)
table(ProbesIslands$Relation_to_UCSC_CpG_Island)
ProbesIslands <- droplevels(ProbesIslands)
table(ProbesIslands$Relation_to_UCSC_CpG_Island)

# How many probes in Islands map on chromosome 1? I will use a double condition:
ProbesIslands_chromosome1 <- Illumina450Manifest_clean[Illumina450Manifest_clean$Relation_to_UCSC_CpG_Island=="Island" & Illumina450Manifest_clean$CHR=="1",]
dim(ProbesIslands_chromosome1)
# The table function allows also to investigate the combination of 2 columns: 
table(Illumina450Manifest_clean$Relation_to_UCSC_CpG_Island,Illumina450Manifest_clean$CHR) 
table(ProbesIslands_chromosome1$Relation_to_UCSC_CpG_Island, ProbesIslands_chromosome1$CHR) 
ProbesIslands_chromosome1 <- droplevels(ProbesIslands_chromosome1)
table(ProbesIslands_chromosome1$Relation_to_UCSC_CpG_Island,ProbesIslands_chromosome1$CHR) 

# How many probes in Islands map also in a gene? Again, I will use a double condition
ProbesIslands_Gene <- Illumina450Manifest_clean[Illumina450Manifest_clean$Relation_to_UCSC_CpG_Island=="Island" & Illumina450Manifest_clean$UCSC_RefGene_Name!="",]
dim(ProbesIslands_Gene)
ProbesIslands_Gene <- droplevels(ProbesIslands_Gene)
str(ProbesIslands_Gene)

# Now we want to identifiy the CpG island (considering also its shores and shelves) with the highest number of probes
table(Illumina450Manifest_clean$UCSC_CpG_Islands_Name)
# Uh, that was a bad idea, it printed a lot of things in the Console!
nlevels(Illumina450Manifest_clean$UCSC_CpG_Islands_Name)
head(table(Illumina450Manifest_clean$UCSC_CpG_Islands_Name))
# I will convert the output of the table function in a data.frame
island_table <- data.frame(table(Illumina450Manifest_clean$UCSC_CpG_Islands_Name))
head(island_table)
# I will order the rows of the dataframe using the order function as follows:
island_table <- island_table[order(island_table[,2]),]
# or, in alternative, I can write:
island_table <- island_table[order(island_table$Freq),]
head(island_table)
tail(island_table)
# I can also order from the largest to the smallest:
island_table <- island_table[order(-island_table[,2]),]
head(island_table)
# I will select two CpG islands with 117 probes each
highest_probes <- Illumina450Manifest_clean[Illumina450Manifest_clean$UCSC_CpG_Islands_Name=="chr6:31830299-31830948"|Illumina450Manifest_clean$UCSC_CpG_Islands_Name=="chr6:31939730-31940559",]
dim(highest_probes) 
# or, in alternative, I can write:
highest_probes <- Illumina450Manifest_clean[Illumina450Manifest_clean$UCSC_CpG_Islands_Name %in% c("chr6:31830299-31830948","chr6:31939730-31940559"),]
dim(highest_probes)
highest_probes <- droplevels(highest_probes)
str(highest_probes)

# Exercise 1: Try to have a look at the probe having the IlmnID "cg19560547", which is in the promoter of WRN gene...does it map to a CpG island or not according to the Illumina Infinium450k manifest?
# Exercise 2: How many probes in the array are in the open sea?
# Exercise 3: How many probes map in genic sequences but not in CpG island, shores and shelves?
# Exercise 4: What is the coordinate (defined according to the MAPINFO column) of the 10th probe on chromosome 7?

## CHUNK4: Genes
levels(Illumina450Manifest_clean$UCSC_RefGene_Name)[1:20]

# Have a look at this:
Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg03817621",]
# See this example: A1CF http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr10%3A52472902-52731702&hgsid=213618992_FWAq8vB3Gv9iJTJYY7xyDqECVmFA


# Have a look at this:
Illumina450Manifest_clean[Illumina450Manifest_clean$Name=="cg00048759",]
#http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A99773884-99776433&hgsid=228211450_jK6pt7GD6xw8kdT4sfPebiyQiI2f


## CHUNK 5: Infinium I and Infinium II probes
# Let's investigate Infinium I and Infinium II probes in the manifest
str(Illumina450Manifest_clean)
table(Illumina450Manifest_clean$Infinium_Design_Type)
table(Illumina450Manifest_clean$Color_Channel)

# What is the number of Infinium I and Infinium II probes in Islands, shores, shelves and not CpG rich regions??
TypeI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dim(TypeI)
TypeI <- droplevels(TypeI)
TypeII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dim(TypeII)
TypeII <- droplevels(TypeII)
table(TypeI$Relation_to_UCSC_CpG_Island)
table(TypeII$Relation_to_UCSC_CpG_Island)

# Alternatively (and much faster):
table(Illumina450Manifest_clean$Relation_to_UCSC_CpG_Island, Illumina450Manifest_clean$Infinium_Design_Type)
table(Illumina450Manifest_clean$Infinium_Design_Type,Illumina450Manifest_clean$Relation_to_UCSC_CpG_Island)

## CHUNK 6: Barplots
# A nice and short tutorial can be found here: http://www.statmethods.net/graphs/bar.html
# I want to generate a barplot of the probes mapping in each chromosome:
counts <- table(Illumina450Manifest_clean$CHR)
barplot(counts)
levels(Illumina450Manifest_clean$CHR)

# R uses alphabetical order by default. I can reorder the levels of the chromosomes:
Illumina450Manifest_clean$CHR <- factor(Illumina450Manifest_clean$CHR,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
levels(Illumina450Manifest_clean$CHR)
barplot(table(Illumina450Manifest_clean$CHR))
barplot(table(Illumina450Manifest_clean$CHR),main="Probes per chromosome",xlab="Chromosome",ylab="Counts")
barplot(table(Illumina450Manifest_clean$CHR),main="Probes per chromosome",xlab="Chromosome",ylab="Counts",col="red")
barplot(table(Illumina450Manifest_clean$CHR),main="Probes per chromosome",xlab="Chromosome",ylab="Counts",col=rainbow(24))
barplot(table(Illumina450Manifest_clean$CHR),main="Probes per chromosome",xlab="Chromosome",ylab="Counts",col=rainbow(12))
barplot(table(Illumina450Manifest_clean$CHR, Illumina450Manifest_clean$Infinium_Design_Type))
barplot(table(Illumina450Manifest_clean$Infinium_Design_Type, Illumina450Manifest_clean$CHR))
barplot(table(Illumina450Manifest_clean$Infinium_Design_Type, Illumina450Manifest_clean$CHR), legend=levels(Illumina450Manifest_clean$Infinium_Design_Type))
barplot(table(Illumina450Manifest_clean$Infinium_Design_Type, Illumina450Manifest_clean$CHR), legend=levels(Illumina450Manifest_clean$Infinium_Design_Type),col=c("blue","red"))
barplot(table(Illumina450Manifest_clean$Infinium_Design_Type, Illumina450Manifest_clean$CHR),beside=T,legend=levels(Illumina450Manifest_clean$Infinium_Design_Type))
barplot(table(Illumina450Manifest_clean$Infinium_Design_Type, Illumina450Manifest_clean$CHR),beside=T,legend=levels(Illumina450Manifest_clean$Infinium_Design_Type),col=c("blue","red"))

# Exercise 5: Create a barplot representing the number of type I and II probes in each region (islands, shores, shelves and open sea)
# Exercise 6 Load the Infinium27k manifest: how many probes are in common between the 27k and 450k? And where do the common probes map (Island, shore, shelves, genic or not genic)

# Supplementary
# Is the ratio between typeI and typeII constant between the chromosomes?
table(Illumina450Manifest_clean$Infinium_Design_Type, Illumina450Manifest_clean$CHR)
df_table <- data.frame(table(Illumina450Manifest_clean$Infinium_Design_Type, Illumina450Manifest_clean$CHR))
df_table
df_table_I <- df_table[df_table$Var1=="I",]
df_table_I <- droplevels(df_table_I)
df_table_II <- df_table[df_table$Var1=="II",]
df_table_II <- droplevels(df_table_II)
head(df_table_I)
head(df_table_II)
# Let's check that the order of the chromosomes is the same in df_table_I and df_table_II
df_table_I[,2]==df_table_II[,2]
table(df_table_I[,2]==df_table_II[,2])
ratio_I_II <- df_table_I[,3]/df_table_II[,3]
barplot(ratio_I_II)
ratio_I_II
# ratio_I_II is jus a vector, without names. If I want to add the chromosome names to the plots
barplot(ratio_I_II,names.arg=levels(df_table_I$Var2))
# or, in alternative:
names(ratio_I_II) <- levels(df_table_I$Var2)
barplot(ratio_I_II)
