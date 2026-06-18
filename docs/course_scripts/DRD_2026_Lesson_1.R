#############################################################################################
## CHUNK 1: Let's start
# Let's print version information about R:
sessionInfo()

# To check the objects that are present in your workspace:
ls()

# It is good practice to remove all the variables from the workspace (unless you need them).
# You can do this in 2 ways:
# 1) Go to Workspace --> Erase workspace
# or 
# 2) digit
rm(list=ls())

#To know the current working directory

getwd()

#To set the new working directory

setwd("~/Dropbox/DRD_2026/1")
getwd()

#############################################################################################
## CHUNK 2: Vectors
#Let's create a vector named V which contains the integer number 1

V1 <- 1
V1

#Let's create a vector named V1 which contains the integer numbers from 1 to 5

V1 <- c(1,2,3,4,5,NA)
V1

#or you can use

V1 <- c(1:5)
V1

#Another example

V1 <- c(1:500)
V1

#Remember: you cannot create an object called "1"
#Try to digit: 1 <- c(1,2,3,4,5)	

#Let's explore the V1 object

length(V1)
class(V1)
str(V1)
summary(V1)
head(V1)
table(V1)

#Now we want to create a vector named V2 which contains decimal numbers from 0 to 2. We will use the function "seq"

?seq
V2 <- seq(from=0,to=10,by=0.1)
V2
length(V2)
class(V2)
str(V2)
summary(V2)
head(V2)
table(V2)

#Now we want to create a vector named V3 which contains the names of the numbers from "one" to "ten". 

V3 <- c("one","two","three","four","five","six","seven","eight","nine","ten")
V3
length(V3)
class(V3)
str(V3)
summary(V3)
head(V3)
table(V3)

#What does happen if we do not use the quotes?

#Try to digit: V3 <- c(one,two,three)	

#Let's try in this way: create 3 objects and concatenate them:

one <- 1
two <- 2
three <- 3
V3 <- c(one,two,three)
V3
length(V3)
class(V3)
str(V3)
summary(V3)
head(V3)
table(V3)

#Now we play a bit with datatypes. First of all, we create a vector of characters, named V4

V4 <- c("one","two","three","one","one","two")
V4
length(V4)
class(V4)
str(V4)
summary(V4)
head(V4)
table(V4)

#An now we create a vector named V5 by converting V4 from datatype=character to datatype=factor

V4
V5 <- factor(V4)
V5
#Note that the levels of the factor are ordered according to the alphabetical order.
#We can reorder the levels:
V5 <- factor(V5, levels=c("one","two","three"))
V5

length(V5)
class(V5)
str(V5)
summary(V5)
head(V5)
table(V5)
levels(V5)
length(levels(V5))
nlevels(V5)

#We can go from factors to numbers...

V5
V6 <- as.numeric(V5)
V6
V5

#...or from factors to characters

V5
V7 <- as.character(V5)
V7
V4

## Indexing a vector
#Let's go back to our object V1, which is a vector including integers from 1 to 500

V1

#Now we want to consider only the 20th element of this vector

V1[20]

#I can index a vector of any datatype:

V2
V2[30]
V5
V5[6]
V5[20]

#############################################################################################
## CHUNK 3: Operators
#I want to subset V1 in order to retain only the elements that satisfy the condition "higher that 31"

V1
#Try to digit: V1_>31 <- V1[V1>31]

#Ops...I tried to create an object with a not correct name. Let's try again:

V1_greater31 <- V1[V1>31]
V1_greater31
summary(V1_greater31)

#Now we want to check if two vectors, V8 and V9, are identical:

V8 <- c(1:10)
V9 <- c(1:10)
V8==V9

#I can store the result of the operation in a new object, V10

V10 <- V8==V9 
V10
str(V10)
summary(V10)
table(V10)

#Let's compare two character vectors:

V11 <- c("A","B","D","F","L")
V12 <- c("Q","B","D","A","L")
V11==V12
V13 <- V11==V12
V13
table(V13)

#Note that if the two vectors do not have the same length, we achieve a warning message:

V11 <- c("A","B","D","F","L")
V12_long <- c("Q","B","D","A","L","Z")
V11==V12_long

#Now we want to know if two vectors are different:

V11!=V12
V14 <- V11!=V12
table(V14)

#And now we want to subet the vector V11 on the basis of the condition that V11 is equal to V12. Note that the condition is checked for each position of the two vectors

V11
V12
V15<- V11[V11==V12]
V15

#Now we want to know the elements of V11 which are included in V12, independently of their position

V11 %in% V12
V15 <- V11[V11 %in% V12]
V15

#Now we want to know the elements of V11 which are not included in V12, independently of their position

#Try to digit: V15 <- V11[V11 !%in% V12]
# Ops, something is wrong! Let's try in this way:
V15 <- V11[!V11 %in% V12]

#############################################################################################
## CHUNK 4: Dataframes
#We use the function "read.table" to load a small dataframe from a txt file:

SampleSheet <- read.table("~/Dropbox/DRD_2026/1/SampleSheet.txt",sep="\t",header=T)
#NB: Replace "~/Dropbox/DRD_2026/1/SampleSheet.txt" with the current file path on your PC

SampleSheet
class(SampleSheet)
dim(SampleSheet)
length(SampleSheet)
str(SampleSheet)
head(SampleSheet)
summary(SampleSheet)

# Important note:
# BEFORE R4.0.0, when using read.table() or data.frame(), R converted strings to factors by default.
# STARTING FROM R4.0.0, when using read.table() or data.frame(), R no longer converts strings to factors by default.
# https://cran.r-project.org/doc/manuals/r-devel/NEWS.html
# https://developer.r-project.org/Blog/public/2022/02/16/stringsasfactors/
# We can set stringsAsFactors = TRUE to convert strings to factors when reading a dataframe

SampleSheet2 <- read.table("~/Dropbox/DRD_2026/1/SampleSheet.txt",sep="\t",header=T, stringsAsFactors = T)
str(SampleSheet2)

#Let's check again the str of SampleSheet...we can specifically convert the column "Sex" to factor:

str(SampleSheet)
SampleSheet$Sex <- factor(SampleSheet$Sex)
str(SampleSheet)

#We can also use exploratory functions on specific columns:

summary(SampleSheet$CD8.naive)

table(SampleSheet$Group)
class(SampleSheet$Group)
summary(SampleSheet$Group)

SampleSheet$Group <- factor(SampleSheet$Group)
table(SampleSheet$Group)
class(SampleSheet$Group)
summary(SampleSheet$Group)
levels(SampleSheet$Group)


## Indexing a dataframe
#We want to select the first 2 columns of the dataframe. Remember: on the left of the comma you operate on rows, on the right of the comma you operate on columns

SampleSheet[,1:2]

#If we select only 1 column we do not have anymore a dataframe...

SampleSheet[,1]

#...unless we set drop=F

SampleSheet[,1,drop=F]

#If we know the name of the column of interest, we can directly indicate it; we can do this in different ways

SampleSheet$Sample_Name
SampleSheet[,"Sample_Name"]
SampleSheet[1:2,]$Sample_Name

#We can operate also on rows...

SampleSheet[1,]

#... or at the same time on rows and columns:

SampleSheet[1:3,1:2]
SampleSheet[1:3,"Sample_Name"]
SampleSheet[1:3,]$Sample_Name
SampleSheet[1:3,3:5]

## Order a dataframe
#We can order the rows of a dataframe according a certain column using the order() function. For example, to order according to the Age:
head(SampleSheet)
#...from the smallest to the largest:
SampleSheet <- SampleSheet[order(SampleSheet$Age),]
head(SampleSheet)
SampleSheet$Age

#...from the largest to the smallest:
SampleSheet <- SampleSheet[order(-SampleSheet$Age),]
head(SampleSheet)
SampleSheet$Age

#############################################################################################
## CHUNK 5: Subsetting a dataframe
#We want retrieve only the rows for which the column "Group" is equal to "W" (remember that now this column is a factor)

class(SampleSheet$Group)
levels(SampleSheet$Group)
#Try to digit: SampleSheet[SampleSheet$Group==W,]
#Ops! Without quotes, W is read by R as an object, and in the Console there is not an object called W. Let's try in this way:

SampleSheet[SampleSheet$Group=="W",]
SampleSheet_W <- SampleSheet[SampleSheet$Group=="W",]
SampleSheet_W
SampleSheet
dim(SampleSheet_W)

#Alternatively:

condition <- SampleSheet$Group=="W"
condition
SampleSheet_W <- SampleSheet[condition,]
SampleSheet_W
dim(SampleSheet_W)

#I can subset on both rows and columns at the same time: see this example...

SampleSheet_W <- SampleSheet[SampleSheet$Group=="W",c("Sample_Name","Group")]
SampleSheet_W
dim(SampleSheet_W)

#...and this example

SampleSheet_W <- SampleSheet[SampleSheet$Group=="W",c(1,2,5)]
SampleSheet_W
dim(SampleSheet_W)

#Let's check the str of SampleSheet_W:

str(SampleSheet_W)

#There is something strange...I have only "W", but R remembers also the "C"! I can avoid this by using the droplevels function:

SampleSheet_W <- droplevels(SampleSheet_W)
dim(SampleSheet_W)
str(SampleSheet_W)

#Another example of subsetting rows according to a condition:

SampleSheet[SampleSheet$CD8.naive>300,]
SampleSheet_CD8.naive_300 <- SampleSheet[SampleSheet$CD8.naive>300,]
str(SampleSheet_CD8.naive_300)
SampleSheet_CD8.naive_300 <- droplevels(SampleSheet_CD8.naive_300)
str(SampleSheet_CD8.naive_300)

#Now we want to subset the dataframe to retain only the samples whose name i stored in the vector v_A_C_D:

v_A_C_D <- c("A","C","D")
SampleSheet
SampleSheet$Sample_Name %in% v_A_C_D
SampleSheet_A_C_D <- SampleSheet[SampleSheet$Sample_Name %in% v_A_C_D,]
SampleSheet_A_C_D
str(SampleSheet_A_C_D)
SampleSheet_A_C_D <- droplevels(SampleSheet_A_C_D)
str(SampleSheet_A_C_D)

#############################################################################################
## CHUNK 6: Save options
# To write a txt file (or a csv file, specifying the comma as separator)
write.table(SampleSheet_A_C_D,file="SampleSheet_A_C_D.txt",sep="\t")
# Let's try to open the file that you have created and check how it looks like. Then, try to save the same file setting row.names=F. What is the difference?
write.table(SampleSheet_A_C_D,file="SampleSheet_A_C_D_rownamesF.txt",sep="\t",row.names=F)

# Saving data into R data formats can reduce considerably the size of large files. The objects can be read back from the file at a later date by using specific functions.

# Saving in Rds format:
saveRDS(SampleSheet_A_C_D, file = "SampleSheet_A_C_D.rds")
# Note that if you save your data as rds, it will be possible to restore the object under a different name .

rm(list=ls())
readRDS("SampleSheet_A_C_D.rds")
SampleSheet_A_C_D <- readRDS("SampleSheet_A_C_D.rds")
SampleSheet_A_C_D
new <- readRDS("SampleSheet_A_C_D.rds")
new
ls()

# Saving in Rdata format: 
save(SampleSheet_A_C_D,file="SampleSheet_A_C_D.RData")
# Note that if you save your data as RData, when you will load it in the workspace the object cannot be restored under different name, but the original object name will be automatically used.
rm(list=ls())
load("SampleSheet_A_C_D.RData")
ls()
rm(list=ls())
new <- load("SampleSheet_A_C_D.RData")

# The function save.image() is just a short-cut for 'save my current workspace', i.e., save(list = ls(all.names = TRUE), file = ".RData", envir = .GlobalEnv). It is also what happens with q("yes").
save.image("Lesson_1_final_part.RData")

#############################################################################################
## CHUNK 7: Exercises:
# Load the SampleSheet file
# Order the rows according to the CD4T column (from the smallest to the largest)
# How many males and how many females do you have?
# How many rows have CD8T <0.05?
# What is the sex of the subjects with CD8T <0.05?
# What is the sex of the subjects with CD8.naive <300?
# What are the mean and the median age of males and females?

#############################################################################################
## CHUNK 8: A larger dataframe: the Infinium 27k manifest
#The file is in a spreadsheet format . Let's import in R the first 20 rows of the file

Infinium27 <- read.table("~/Dropbox/DRD_2026/1/Infinium27.txt",sep="\t",header=T,nrow=20)
dim(Infinium27)
head(Infinium27)

#Now we import the entire dataframe:

Infinium27  <- read.table("~/Dropbox/DRD_2026/1/Infinium27.txt",sep="\t",header=T)
dim(Infinium27)
head(Infinium27)
str(Infinium27)

#Let's set stringsAsFactors = TRUE

Infinium27  <- read.table("~/Dropbox/DRD_2026/1/Infinium27.txt",sep="\t",header=T,stringsAsFactors = T)
dim(Infinium27)
head(Infinium27)
str(Infinium27)

#How many probes are there in each chromosome?
#We can proceed in two ways:
#1) We subset the dataframe for each chromosome

chr1 <- Infinium27[Infinium27$Chr=="1",]
dim(chr1)
chr2 <- Infinium27[Infinium27$Chr=="2",]
dim(chr2)

#etc
#2) We use the table function

table(Infinium27$Chr)

#Are there missing values?

str(Infinium27)
TSS_NA <- Infinium27[is.na(Infinium27$TSS_Coordinate),]
dim(TSS_NA)
head(TSS_NA)
TSS_NA_2 <- Infinium27[Infinium27$TSS_Coordinate=="NA",]
dim(TSS_NA_2)

#How many probes are not associated to a gene?

Infinium27_No_gene <- Infinium27[Infinium27$Symbol=="",]
dim(Infinium27_No_gene)
str(Infinium27_No_gene)
Infinium27_No_gene <- droplevels(Infinium27_No_gene)
str(Infinium27_No_gene)

#How many genes there are on each chromosome?

str(Infinium27$Gene_ID)
Infinium27_gene <- Infinium27[Infinium27$Gene_ID!="",]
str(Infinium27_gene)
Infinium27_gene <- droplevels(Infinium27_gene)
dim(Infinium27_gene)
str(Infinium27_gene)
summary(Infinium27_gene)
table(Infinium27_gene$Chr)
table(Infinium27_gene$Chr, Infinium27_gene$Gene_Strand)
df_table <- data.frame(table(Infinium27_gene$Chr))
df_table
colnames(df_table) <- c("chr","number of genes")
df_table

#What is the gene with the highest number of probes?

gene_table <- table(Infinium27$Symbol)
head(gene_table)
str(gene_table)
gene_table <- data.frame(gene_table)
str(gene_table)
colnames(gene_table) <- c("Gene","Freq")
head(gene_table)
gene_table <- gene_table[order(gene_table$Freq),]
head(gene_table)
tail(gene_table)
gene_table <- gene_table[order(-gene_table$Freq),]
head(gene_table)
gene_table[gene_table$Freq == max(gene_table$Freq),]
gene_table <- gene_table[!is.na(gene_table$Gene),]
head(gene_table)
gene_table <- gene_table[gene_table$Gene !="",]
head(gene_table)
gene_table <- droplevels(gene_table)
gene_table[gene_table$Freq == max(gene_table$Freq),]

#############################################################################################
## CHUNK 9: Exercises:
#What is the chromosome with the highest number of probes?

#What is the chromosome with the highest number of probes with SourceStrand equal to TOP?

#On which chromosomes are the probes not associated to a gene?

#What is the maximum distance of a probe from TSS? And the smallest?

#How many CpG islands are on chromosome 7? 	

  
