# In this lesson we will use preprocessQuantile as normalization method for our data. Let's recalculate beta values after preprocessQuantile normalization and read the file containing the samplesheet:
# We will recalculate the preprocessQuantile_results object
rm(list=ls())
library(minfi)
setwd("../Lesson5")
suppressMessages(library(minfi))

load("../Lesson3/RGset.RData")
preprocessQuantile_results <- preprocessQuantile(RGset)

## CHUNK 1: homoscedasticity and heteroscedasticity

# We want to check for homoscedasticity and heteroscedasticity in beta and M values

beta_preprocessQuantile <- getBeta(preprocessQuantile_results)
M_preprocessQuantile <- getM(preprocessQuantile_results)

beta_preprocessQuantile_mean <- apply(beta_preprocessQuantile,1, mean,na.rm=T)
beta_preprocessQuantile_stdev <- apply(beta_preprocessQuantile,1, sd,na.rm=T)

M_preprocessQuantile_mean <- apply(M_preprocessQuantile,1, mean,na.rm=T)
M_preprocessQuantile_stdev <- apply(M_preprocessQuantile,1, sd,na.rm=T)

# Now we want to plot, for each probe, sd against mean. The plot will contain 485512 dots, they are really a lot. For this reason in this case it is more convenient to use the function smoothScatter() instead of plot()

? smoothScatter
smoothScatter(beta_preprocessQuantile_mean, beta_preprocessQuantile_stdev)
?lowess
lines(lowess(beta_preprocessQuantile_mean, beta_preprocessQuantile_stdev), col="red") # lowess carries out a locally weighted regression of y on x

pdf("Beta_M_heteoschedasticity_omoscedasticity.pdf")
par(mfrow=c(1,2))
smoothScatter(beta_preprocessQuantile_mean, beta_preprocessQuantile_stdev)
lines(lowess(beta_preprocessQuantile_mean, beta_preprocessQuantile_stdev), col="red") # lowess carries out a locally weighted regression of y on x
smoothScatter(M_preprocessQuantile_mean, M_preprocessQuantile_stdev)
lines(lowess(M_preprocessQuantile_mean, M_preprocessQuantile_stdev), col="red") # lowess line (x,y)
dev.off()

# For this dataset, we do not observe large differences between Beta and M values!


## CHUNK 2: Identification of Differentially methylated positions
# We will perform the differential methylation analysis on beta values from data normalized using preprocessQuantile. The analytical steps performed in this tutorial can be applied to any normalization method and to both beta values and M values.  

# Let's read the file containing the samplesheet:
pheno <- read.csv("../Lesson3/Input_data/SampleSheet.csv",header=T, stringsAsFactors=T)
str(pheno)

# Our aim is to identify the CpG probes that are differentially methyled between samples belonging to groups A and B. We will use a parametric test (t-test) and a non parametric test (Mann-Whitney test). Please note that these statistical appraoches are very basic, and usually researchers adopt much more sophisticated statistical approaches, which allows also to take into account the effect of potential confounding factors. However, in our examples we will keep the things simple, and we will learn how to extract the info in which we are interested in (in this case, the p-value) from the output of a function. 

# Let's use the t-test to compare the first CpG of beta_preprocessQuantile between groups A and B:
? t.test
t_test <- t.test(beta_preprocessQuantile[1,] ~ pheno$Group)
t_test


# Let's use the Mann-Whitney test to compare the first CpG of beta_preprocessQuantile between groups A and B:
wilcox <- wilcox.test(beta_preprocessQuantile[1,]~ pheno$Group)
wilcox

# Note that both the t.test and the wicox.test functions return a list:
str(t_test)
str(wilcox)

# We are interested in the "p.value" element of these lists.
t_test$p.value
wilcox$p.value

# For 1 row (1 CpG probe) it is easy, but...how can I apply the 2 tests to each row of the dataframe (that is, to each probe of the microarray) and extract the p-value?
# We will use the apply() function, that we have already met in Lesson 4 to calculate the mean of each row of the matrix. However, unlike the mean() function, t.test() and wilcox.test() functions do not return a value, but a list; therefore, we have to create an ad hoc function:

My_ttest_function <- function(x) {
  t_test <- t.test(x~ pheno$Group)
  return(t_test$p.value)
} 

# Let's apply the function to few rows of beta_preprocessQuantile
firstbeta_preprocessQuantile <- beta_preprocessQuantile[1:10,]
pValues_ttest <- apply(firstbeta_preprocessQuantile,1, My_ttest_function)
pValues_ttest

# Great, everything works fine!

# Let's do the same with the non-parametric test
My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x~ pheno$Group)
  return(wilcox$p.value)
} 
pValues_wilcox <- apply(firstbeta_preprocessQuantile,1, My_mannwhitney_function)
pValues_wilcox

# You can apply your brand new functions to all the rows of the beta_preprocessQuantile matrix. As this step can take several minutes, in this turial we will apply the parametric and non parametric functions only to the first 20000 probes

first20k_beta_preprocessQuantile <- beta_preprocessQuantile[1:20000,]
pValues_ttest_first20k <- apply(first20k_beta_preprocessQuantile,1, My_ttest_function)
length(pValues_ttest_first20k)
pValues_wilcox_first20k <- apply(first20k_beta_preprocessQuantile,1, My_mannwhitney_function)
length(pValues_wilcox_first20k)

# We can create a data.frame with all the beta values and the pValue column
final_ttest_first20k <- data.frame(first20k_beta_preprocessQuantile, pValues_ttest_first20k)
head(final_ttest_first20k)
dim(final_ttest_first20k)
final_wilcox_first20k <- data.frame(first20k_beta_preprocessQuantile, pValues_wilcox_first20k)
head(final_wilcox_first20k)
dim(final_wilcox_first20k)

# We can order the probes on the basis of the pValues column (from the smallest to the largest value)
final_ttest_first20k <- final_ttest_first20k[order(final_ttest_first20k$pValues_ttest_first20k),]
head(final_ttest_first20k)

final_wilcox_first20k <- final_wilcox_first20k[order(final_wilcox_first20k$pValues_wilcox_first20k),]
head(final_wilcox_first20k) #Note that in these settings the pvalues from Mann-Whitney tests are discrete

# How many probes have a pValue<=0.05?
final_ttest_first20k_0.05 <- final_ttest_first20k[final_ttest_first20k$pValues_ttest_first20k<=0.05,]
dim(final_ttest_first20k_0.05)
final_wilcox_first20k_0.05 <- final_wilcox_first20k[final_wilcox_first20k$pValues_wilcox_first20k<=0.05,]
dim(final_wilcox_first20k_0.05)

intersection <- intersect(rownames(final_ttest_first20k_0.05),rownames(final_wilcox_first20k_0.05))
length(intersection)

## CHUNK 3: dmpFinder
?dmpFinder
# As you can see, dmpFinder is based on limma

pValues_dmpFinder_first20k <- dmpFinder(first20k_beta_preprocessQuantile,pheno=pheno$Group,type="categorical")
length(pValues_dmpFinder_first20k)
head(pValues_dmpFinder_first20k)

table(rownames(final_ttest_first20k)==rownames(pValues_dmpFinder_first20k))
table(rownames(final_ttest_first20k)%in%rownames(pValues_dmpFinder_first20k))

# We will use the merge() function to merge the final_ttest_first20k with the pValues_dmpFinder_first20k dataframe
final_dmpFinder_first20k <- merge(first20k_beta_preprocessQuantile,pValues_dmpFinder_first20k,by="row.names")
head(final_dmpFinder_first20k)
rownames(final_dmpFinder_first20k) <- final_dmpFinder_first20k[,1]
head(final_dmpFinder_first20k)

final_dmpFinder_first20k <- final_dmpFinder_first20k[order(final_dmpFinder_first20k$pval),]
head(final_dmpFinder_first20k)
# How many probes have a pValue<=0.05?
final_dmpFinder_first20k_0.05 <- final_dmpFinder_first20k[final_dmpFinder_first20k$pval<=0.05,]
dim(final_dmpFinder_first20k_0.05)

# Lest's visualize the intersection between the 3 approaches that we used; we will prepare a Venn Diagram using the VennDiagram function
install.packages("VennDiagram")
library(VennDiagram)
?venn.diagram
venn.diagram(list(A=rownames(final_dmpFinder_first20k_0.05),B=rownames(final_ttest_first20k_0.05),C=rownames(final_wilcox_first20k_0.05)),filename = "Venn_differentially_methylated_positions.png")

## CHUNK 4: analysis of a large dataset
# Now let's try to perform parametric and non parametric tests on a larger dataset (chr18_large_example.RData). I have already prepared an R object containing chr18 beta values for 64 subjects, 32 normal and 32 with cancer. The pheno vector is:
pheno_large <- factor(c(rep("normal",32),rep("cancer",32)))
pheno_large

load("../Lesson5/chr18_large_example.RData")
ls()
head(chr18_large_example)

save(final_ttest_first20k,file="final_ttest_first20k.RData")

# Build a function to perform ttest and a function to perform Mann-Whitney test
