rm(list=ls())
setwd("D:/Google Drive/UNIBO/Materials/DNA & RNA Dynamics/Module 2/Lesson4")
suppressMessages(library(minfi))

## CHUNK 1: Beta and M values
# We load the MSet object, which contains the methylation and unmethylation signals
load("../Lesson3/MSet_raw.RData")
ls()
MSet.raw
# Remember that it is a MethylSet class object. If we want to know the accessor functions for this class of objects, we can go to the help page:
?MethylSet
# The functions getBeta and getM allow to retrieve beta and M values matrices:
beta <- getBeta(MSet.raw)
class(beta)
dim(beta)
head(beta)
summary(beta)
# In our dataset, minimum for beta values is 0, maximum is 1. We have also some NAs: these are the positions for which both the Methylation and Unmethylation values in the MSet.raw are equal to 0 (by default, in minfi the offset is set to 0)

M <- getM(MSet.raw)
dim(M)
head(M)
summary(M)
# In our dataset, minimum for M values is -Inf (Methylation value=0, Unmethylation value>0), maximum is +Inf (Methylation value>0, Unmethylation value=0). We have also some NAs: these are the positions for which both the Methylation and Unmethylation values in the MSet.raw are equal to 0 (by default, in minfi the offset is set to 0)

# Note! The output of getBeta and getM functions is a matrix. We can also convert the output of these functions to a dataframe (as we have done with the output of getRed and getGreen functions in the last lesson):
beta_df <- data.frame(getBeta(MSet.raw))
class(beta_df)
# Do you note some differences between beta and beta_df?
dim(beta)
dim(beta_df)
head(beta)
head(beta_df)

# Now we want to plot the density distribution of mean (across the 8 samples) beta and M values for each CpG

# First of all, we have to calculate the mean of each row of the matrix across the 8 samples.

# To calculate the mean, we will use the mean() function. As we know that in our input there are some missing values, we will specify that NA values should be stripped (see this page https://www.statmethods.net/input/missingdata.html) using the na.rm=T argument.

# To calculate the mean of each row of the matrix, we will use the apply() function, which takes as input a dataframe or a matrix and returns a vector (see here for example: http://petewerner.blogspot.it/2012/12/using-apply-sapply-lapply-in-r.html). Note that the "MARGIN" argument in the function allows to specify that the function mean() should be applied to each row (set "MARGIN" to 1) or to each column (set "MARGIN" to 2) of the matrix/dataframe.

mean_of_beta <- apply(beta,1,mean,na.rm=T)
head(beta)
head(mean_of_beta)
dim(beta)
length(mean_of_beta)

# Now we can calculate the density distribution
?density
d_mean_of_beta <- density(mean_of_beta)
d_mean_of_beta
plot(d_mean_of_beta,main="Density of Beta Values",col="orange")

# We apply the same steps on M values:
mean_of_M <- apply(M,1,mean,na.rm=T)
d_mean_of_M <- density(mean_of_M)
plot(d_mean_of_M)
plot(d_mean_of_M,main="Density of M Values",col="purple")

# We can put several plots in the same panel by setting the option "mfrow" in the par() function (http://www.statmethods.net/advgraphs/layout.html). Yu can find a quick guide to plots here: http://www.statmethods.net/graphs/index.html
par(mfrow=c(1,2))
plot(d_mean_of_beta,main="Density of Beta Values",col="orange")
plot(d_mean_of_M,main="Density of M Values",col="purple")
# or
par(mfrow=c(2,1))
plot(d_mean_of_beta,main="Density of Beta Values",col="orange")
plot(d_mean_of_M,main="Density of M Values",col="purple")

# Another useful plot to represent our data is boxplots:
boxplot(beta,main="beta values")

## CHUNK 2: Save plots
# You can save the graph in a variety of formats from the menu File -> Save As.
# You can also save the graph via code using one of the following functions:

# Function 						Output to

# pdf("mygraph.pdf") 			pdf file

# win.metafile("mygraph.wmf") 	windows metafile

# png("mygraph.png") 			png file

# jpeg("mygraph.jpg") 			jpeg file

# bmp("mygraph.bmp") 			bmp file

# postscript("mygraph.ps") 		postscript file

# For example (the file will be saved in the working directory):
pdf("Picture1.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(d_mean_of_beta,main="Density of Beta Values",col="blue")
boxplot(beta)
dev.off() #shuts down the graphical device


## CHUNK 3: plot the density distribution of beta values for Infinium I and II probes
# We want to subset the dataframe of beta values according to Type I and II
# First of all, we need to know what are the Type I and Type II probes. To this aim, we need our manifest file.
load('../Lesson2/Illumina450Manifest_clean.RData')
ls()
str(Illumina450Manifest_clean)

# We subset the Illumina450Manifest_clean in two dataframes, containing only type I (dfI) or type II (dfII) probes
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dim(dfI)
str(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

# Remember that in the beta matrix, which contains the beta values, the names of the probes are stored in the rownames of the matrix:
head(beta)
# Now we subset the beta matrix in order to retain only the rows whose name is in the first column of dfI...
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
dim(beta_I)
# ... or in the first column of dfII
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]
dim(beta_II)

# For each probe in the mean_of_beta_I and mean_of_beta_II matrices, we calculate the mean of beta values across the 8 samples...
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)

# ... and then we calculate the density distribution of the 2 vectors of mean values:
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

# Finally, we can plot the densities.
plot(d_mean_of_beta_I,col="blue")
plot(d_mean_of_beta_II,col="red")

# But in this case it is more meaningful to overlay the two density distributions. To this aim, we first open a plot using the plot() function, and then we add a new line, coloured with a different color, using the lines() function
plot(d_mean_of_beta_I,col="blue")
lines(d_mean_of_beta_II,col="red")

# Note that the density distribution of beta values drom Infinium type II probes is shifted towards the centre. Moreover,type I probes have higher density at low methylation levels.


## CHUNK 4: Normalization
# We will apply the preprocessQuantile function. To evaluate the effect of the normalization, we will compare raw and normalized data according to 3 types of plots:

# - densities of mean beta values for type I and type II probes

# - densities of standard deviation of beta values for type I and type II probes

# - boxplots of beta values


# For Raw data, we already have the matrix of beta values, the d_mean_of_beta_I and d_mean_of_beta_II objects; we need to calculate the densities of the standard deviations, which can be calculated using the function sd():
sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,)
d_sd_of_beta_II <- density(sd_of_beta_II)

# Now we can perform the normalization using the preprocessQuantile function.
?preprocessQuantile
# According to the help page, the input can be an RGset  or MehylSet object. Let's load the RGset object:
load("../Lesson3/RGset.RData")
RGset

preprocessQuantile_results <- preprocessQuantile(RGset)
# You can have the following error: "there is no package called ‘IlluminaHumanMethylation450kanno.ilmn12.hg19"; in this case you have to install the missing package (it will take a bit):
# if (!requireNamespace("BiocManager", quietly = TRUE))

#     install.packages("BiocManager")

# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# preprocessQuantile_results <- preprocessQuantile(RGset)

str(preprocessQuantile_results)
class(preprocessQuantile_results)
preprocessQuantile_results
?GenomicRatioSet

# You see that getBeta is among the accessor functions of a GenomicRatioSet class object
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)
head(beta_preprocessQuantile)

save(beta_preprocessQuantile, file="beta_preprocessQuantile.RData")

# Now, as we have done before for raw data, we divide the beta_preprocessQuantile matrix according to type I and type II probes, calculate the mean and the standartd deviation for each probe across the 8 samples and calculate the density distributions
beta_preprocessQuantile_I <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfI$IlmnID,]
beta_preprocessQuantile_II <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfII$IlmnID,]
mean_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,mean)
mean_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,mean)
d_mean_of_beta_preprocessQuantile_I <- density(mean_of_beta_preprocessQuantile_I,na.rm=T)
d_mean_of_beta_preprocessQuantile_II <- density(mean_of_beta_preprocessQuantile_II,na.rm=T)
sd_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,sd)
sd_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,sd)
d_sd_of_beta_preprocessQuantile_I <- density(sd_of_beta_preprocessQuantile_I,na.rm=T)
d_sd_of_beta_preprocessQuantile_II <- density(sd_of_beta_preprocessQuantile_II,na.rm=T)

# And here it is our beatiful plot!
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta")
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd")
lines(d_sd_of_beta_II,col="red")
boxplot(beta)
plot(d_mean_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile beta")
lines(d_mean_of_beta_preprocessQuantile_II,col="red")
plot(d_sd_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile sd")
lines(d_sd_of_beta_preprocessQuantile_II,col="red")
boxplot(beta_preprocessQuantile)

# Note that it is easier to compare the plots if we have the same scales on x and y axes; to this aim, we can specify the xlim and ylim arguments:
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_II,col="red")
boxplot(beta,ylim=c(0,1))
plot(d_mean_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_preprocessQuantile_II,col="red")
plot(d_sd_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_preprocessQuantile_II,col="red")
boxplot(beta_preprocessQuantile,ylim=c(0,1))

# We can save this plot in a pdf
pdf("Plot_comparison_raw_preprocessQuantile.pdf",height=7,width=15)
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_II,col="red")
boxplot(beta,ylim=c(0,1))
plot(d_mean_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_preprocessQuantile_II,col="red")
plot(d_sd_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile sd",xlim=c(0,0.6),ylim=c(0,60))
lines(d_sd_of_beta_preprocessQuantile_II,col="red")
boxplot(beta_preprocessQuantile,ylim=c(0,1))
dev.off()

# Try by yourself preprocessFunnorm, preprocessNoob and preprocessSWAN
