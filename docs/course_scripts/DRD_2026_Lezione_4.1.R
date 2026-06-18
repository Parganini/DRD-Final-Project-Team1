## CHUNK 1: PCA
# PCA can be used as a diagnostic plot for the detection of outliers and batch effects.
# We will use the function prcomp() to calculate the PCA on our matrix of normalized beta values.

#Let's load our quantile-normalized DNA methylation data (Lesson 4) and the samplesheet with phenotypic information (Lesson 3):
load("../Lesson4/beta_preprocessQuantile.RData")
str(beta_preprocessQuantile)
pheno <- read.csv("../Lesson3/Input_data/SampleSheet.csv",header=T, stringsAsFactors=T)
str(pheno)

# Our beta matrix has samples in columns and CpG probes in rows: the prcomp() function should be applied to the transposed matrix, which is achieved using the t() function
?prcomp
pca_results <- prcomp(t(beta_preprocessQuantile),scale=T)

# We can print and plot the variance accounted for each component
print(summary(pca_results))
plot(pca_results)
# We ask the str() of pca_object
str(pca_results)

# The principal components of interest are stored in the element named "x" of the list.
pca_results$x

# We can plot PC1 and PC2 and check if samples cluster according to some variable; the argument "cex" defines the size of the dots, the argument "pch" the dot type:
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2)

# We can label the dots using the text() function and setting the "labels" argument to the rownames of pca_results$x; the argument "pos" defines the position of the text:
?text
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),pos=1)
# Sample 9376538140_R02C01 seems to be an outlier;

# In this dataset we have few samples, so it is ok to plot their names. When you work with larger datasets, on the contrary, it is easier to colour the dots according to some column of the pheno file of interest, for example Group:
pheno$Group

# We define a palette of colours; each colour will be assigned to each levels of the factor according to the order of levels. For example, in this case Group A will be in orange, Group B in purple. We will also adjust the margins and add a legend to the plot.

levels(pheno$Group)
palette(c("orange","purple"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=pheno$Group,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=2)

# Another example coloring according to the sex of the subjects
levels(pheno$Gender)
palette(c("pink","blue"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=pheno$Gender,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Gender),col=c(1:nlevels(pheno$Group)),pch=2)

## CHUNK 2: Batch Effect Removal
# Omics experiments are sometimes affected by batch effects, that is unwanted variations in the data that arise due to external factors like experimental conditions or equipment, rather than biological differences.

#Let's preview our data to check for batch effects
levels(pheno$Slide)
#pheno$Slide has no levels because it is read as a numeric vector. Convert it into a factor for correct pca plot coloring. 
pheno$Slide <- factor(pheno$Slide)

palette(c("grey","orange"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=2,col=pheno$Slide,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Slide),col=c(1:nlevels(pheno$Group)),pch=2, cex=0.5)

boxplot(beta_preprocessQuantile, col=pheno$Slide)

# We do not see a clear batch effect in our dataset, good!
# Although our data seems to be devoid of batch effects let's check what happens if we apply batch effect correction.

# Batch effect correction on DNA methylation data is performed using the *ComBat()* function, which is part of the sva package. 
# ComBat() adjusts for known batches using an empirical Bayesian framework, helping to remove these batch effects and make the data more reliable for downstream analysis.

# Package sva can be downloaded and installed by running:
# BiocManager::install("sva")

library(sva)

# *ComBat()* works on a matrix with features (in this case, methylation values) in the rows and samples in the columns. We will create a first model, in which we will include the variable in which we are interested in (in this case, Sample_Group) and an object named batch in which we will include the known source of batch (in this case, the column Slide).
mod.null <- model.matrix(~as.factor(Group), data = pheno)
batch <- as.factor(pheno$Slide)
combat_edata <- ComBat(dat = beta_preprocessQuantile, batch = batch, mod = mod.null)  
str(combat_edata)
Quantile_Combat_beta <- as.data.frame(combat_edata)

#And now let's see how ComBat has modified our data! We will calculate again PCA and plot PC1 and PC2 according to Group, Gender and Slide.

Tbeta_ComBat <- t(Quantile_Combat_beta)
pca_result_ComBat <- prcomp(Tbeta_ComBat, scale = TRUE)

levels(pheno$Group)
palette(c("orange","purple"))
plot(pca_result_ComBat$x[,1], pca_result_ComBat$x[,2],cex=2,pch=2,col=pheno$Group,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_result_ComBat$x[,1], pca_result_ComBat$x[,2],labels=rownames(pca_result_ComBat$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=2)

levels(pheno$Gender)
palette(c("pink","blue"))
plot(pca_result_ComBat$x[,1], pca_result_ComBat$x[,2],cex=2,pch=2,col=pheno$Gender,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_result_ComBat$x[,1], pca_result_ComBat$x[,2],labels=rownames(pca_result_ComBat$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=2)

levels(pheno$Slide)
palette(c("grey","orange"))
plot(pca_result_ComBat$x[,1], pca_result_ComBat$x[,2],cex=2,pch=2,col=pheno$Slide,xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_result_ComBat$x[,1], pca_result_ComBat$x[,2],labels=rownames(pca_result_ComBat$x),cex=0.5,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=2)

boxplot(Quantile_Combat_beta,col=pheno$Slide)

# Batch has minimal effect on data variability. 
# Batch effect correction did not result in large changes in beta value distribution nor in data variability.

# Exercise 1: Load beta_CRC.RData and Samplesheet_CRC.txt. Plot beta value distribution (using boxplot) and Principal Components (using PCA plot) and check for batch effects.
# Do you notice a batch effect on data? Remove batch effect from data using ComBat() function. Check the dataset again using the plotting functions above. Are there any differences compared to uncorrected data?
