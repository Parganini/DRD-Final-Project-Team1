rm(list=ls())
setwd("../Lesson6")

# For these last analyses, we will use the object final_ttest_first20k, that we previously prepared by applying the ttest function to the first 20000 rows od the beta_preprocessQuantile dataframe. If it is not present in your workspace or in Lesson 5 folder, you can download it from Lesson 6 folder:
load("../Lesson5/final_ttest_first20k.RData")
str(final_ttest_first20k)
ls()
# If the pheno pbject is not present in your workspace, you can load it again from Lesson 3 folder:
pheno <- read.csv("../Lesson3/Input_data/SampleSheet.csv",header=T, stringsAsFactors=T)
str(pheno)

## CHUNK 1: multiple test correction
# We will use the function p.adjust to perform the multiple test correction. For the sake of semplicity, we will focus on the results of the ttest, but you can apply the same procedure also to the p-values resulting from the non parametric test.
?p.adjust
# We will apply the Benjamini & Hochberg and the Bonferroni corrections:
corrected_pValues_BH <- p.adjust(final_ttest_first20k$pValues_ttest_first20k,"BH")
corrected_pValues_Bonf <- p.adjust(final_ttest_first20k$pValues_ttest_first20k,"bonferroni")
final_ttest_first20k_corrected <- data.frame(final_ttest_first20k, corrected_pValues_BH, corrected_pValues_Bonf)
head(final_ttest_first20k_corrected)
# We can visualize the distributions of the p-values and of the corrected p-values by boxplots:
colnames(final_ttest_first20k_corrected)
boxplot(final_ttest_first20k_corrected[,9:11])

# How many probes survive the multiple test correction?
dim(final_ttest_first20k_corrected[final_ttest_first20k_corrected$pValues_ttest<=0.05,])
dim(final_ttest_first20k_corrected[final_ttest_first20k_corrected$corrected_pValues_BH<=0.05,])
dim(final_ttest_first20k_corrected[final_ttest_first20k_corrected$corrected_pValues_Bonf<=0.05,])


## CHUNK 2: Volcano plots
# First of all, we have to calculate the difference between the averge of group A values and the average of group B values. To this aim, we will first create two matrixes containing the beta-values of group A and group B samples, and then we will calculate the mean within each group for each row.
pheno$Group

beta_first20k <- final_ttest_first20k_corrected[,1:8]

beta_first20k_groupA <- beta_first20k[,pheno$Group=="A"]
mean_beta_first20k_groupA <- apply(beta_first20k_groupA,1,mean)
beta_first20k_groupB <- beta_first20k[,pheno$Group=="B"]
mean_beta_first20k_groupB <- apply(beta_first20k_groupB,1,mean)

# Now we can calculate the difference between average values:
delta_first20k <- mean_beta_first20k_groupB-mean_beta_first20k_groupA
head(delta_first20k)

# Now we create a dataframe with two columns, one containing the delta values and the other with the -log10 of p-values
toVolcPlot <- data.frame(delta_first20k, -log10(final_ttest_first20k_corrected$pValues_ttest))
head(toVolcPlot)
plot(toVolcPlot[,1], toVolcPlot[,2])

# Let's improve a bit the plot: http://www.statmethods.net/advgraphs/parameters.html
# We can also add a threshold for pvalue significance (for example, 0.01) by using the abline function

plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
-log10(0.01)
?abline
abline(h=-log10(0.01),col="red")

# Now I want to highlight the probes (that is, the points), that have a nominal pValue<0.01 and a delta > 0.1)
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
toHighlight <- toVolcPlot[toVolcPlot[,1]>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="yellow")

# If I want to highlight the points with an absolute delta > 0.01
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5)
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="red")

## CHUNK 3: Manhattan plots
# We will use the "qqman" package, from the CRAN repository:
# to install it type: 
# install.packages("qqman")
library(qqman)

# To calculate the Manhattan plot, during the lesson we will use the object final_ttest_reduced, that I prepared by applying the ttest function on a subset of the probes of the dataset. The object is in the Lesson 6 folder:
load("final_ttest_reduced.RData")
# First we have to annotate our dataframe, that is add genome annotation information for each cpg probe. We will use the Illumina450Manifest_clean object that we previously created:
load('../Lesson2/Illumina450Manifest_clean.RData')

# We will use the merge() function to merge the final_ttest_corrected with the Illumina450Manifest_clean object
?merge
# The merge function performs the merging by using a column which is common to two dataframes and which has the same name in the two dataframes
head(Illumina450Manifest_clean)
head(final_ttest_reduced)
# We want to merge on the basis of the CpG probes, but unfortunately in the final_ttest_corrected object the CpG probes are stored in the rownames, not in a column. We can overcome this issue as follows:
final_ttest_reduced <- data.frame(rownames(final_ttest_reduced),final_ttest_reduced)
head(final_ttest_reduced)
colnames(final_ttest_reduced)
colnames(final_ttest_reduced)[1] <- "IlmnID"
colnames(final_ttest_reduced)

final_ttest_reduced_annotated <- merge(final_ttest_reduced, Illumina450Manifest_clean,by="IlmnID")
dim(final_ttest_reduced)
dim(Illumina450Manifest_clean)
dim(final_ttest_reduced_annotated)
str(final_ttest_reduced_annotated)

# Note that the dataframe is automathically reordered on the basis of alphabetical order of the column used for merging
head(final_ttest_reduced_annotated)

# Now we can create the input for the Manhattan plot analysis. The input object should contain 4 info: probe, chromosome, position on the chromosome and p-value. We will select these columns in the final_ttest_corrected_annotated object
input_Manhattan <- final_ttest_reduced_annotated[colnames(final_ttest_reduced_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_ttest")]
dim(input_Manhattan)
head(input_Manhattan)
str(input_Manhattan$CHR)
levels(input_Manhattan$CHR)
# It is better to reorder the levels of the CHR
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )
levels(input_Manhattan$CHR)

# The function that we will use is "manhattan"
?manhattan
# As you see, the column "CHR" should be numeric --> we will convert factors to numbers
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)

# and finally we can produce our Manhattan plot
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest" )
# what is the blue line?
-log10(0.00001)
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest",annotatePval = 0.00001,col=rainbow(24) )


## CHUNK 4: heatmaps with hierarchical clustering
# We need the gplots package from CRAN repository
#install.packages("gplots")
library(gplots)

# We will use the function heatmap.2. 
?heatmap.2
# It wants as input a matrix. As an heatmap on several thousands of probes is computationally demanding, for our analysis we will use just the top 100 most significant CpG probes. We will extract only beta values for these top 100 probes and convert them in a matrix:
input_heatmap=as.matrix(final_ttest_first20k[1:100,1:8])

#The default method for distance is "euclidean" (click on "dist" in distfun paragraph), while the default function for hierarchical clustering is "complete" (click on "hclust" in hclustfun paragraph)

# Finally, we we create a bar of colors for Group A or Group B membership
pheno$Group
colorbar <- c("green","green","orange","orange","green","green","orange","orange")

# In the following lines we will compare the results of hierchical clustering using different methods.

# Complete (default options)
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")

# You can see that hierarchical clustering divides well group A and group B samples; in addition, you see that some probes are hypermethylated in Group A compared to group B, others are hypomethylated.

# Single
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Single linkage")

# Average
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Average linkage")

# You can set your palette of colors
col2=colorRampPalette(c("green","black","red"))(100)
heatmap.2(input_heatmap,col=col2,Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
