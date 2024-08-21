library(tidyverse)
remotes::install_version('Matrix', version = '1.6.4') #Required for installing Seurat
install.packages("Seurat")
library(Seurat)
library(SeuratObject)
library(patchwork)
library(Matrix)

#Read in UMI count date
UMIcounts <- read_delim('GSE181919_UMI_counts.txt') #Didn't work, too large I think. 
UMI_counts <- data.table::fread('GSE181919_UMI_counts.txt')

head(UMI.mat[,1:10])
rnames <- as.vector(UMI_counts$V1) #Saves gene names as character vector
UMI_count.df <- as.data.frame(UMI_counts)
UMI.mat <- data.matrix(UMI_count.df[,2:ncol(UMI_count.df)]) #Creates matrix and drops gene names as column
rownames(UMI.mat) <- rnames #Assigns character vector of gene names as row names
UMI.mat <- as.sparse(UMI.mat) #Necessary for
#Remember to free unused memory
#Create Seurat object
?Seurat
Sobject <- CreateSeuratObject(counts = UMI.mat, project = "Bish", min.features = 200, min.cells = 3)

#A few stats on matrix
dim(UMI.m)
dimnames(UMI.m)
class(UMI.mat)
sum(duplicated(rownames(UMI.mat)))

#Stats on seurat object
dim(Sobject)
class(Sobject)               
head(Sobject[,1:10])

##QC and filtering -----

Sobject[["percent.mt"]] <- PercentageFeatureSet(Sobject, pattern = "^MT-") #Adds column with % of mitochondrial genes
head(Sobject, 5)

#Visualize QC parameters
VlnPlot(Sobject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Sobject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Filtering

Sobject <- subset(Sobject, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)

#Normalization

Sobject <- NormalizeData(Sobject, normalization.method = "LogNormalize", scale.factor = 10000)

