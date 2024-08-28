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

#ID variable features

Sobject <- FindVariableFeatures(Sobject, selection.method = "vst")

top10 <- head(VariableFeatures(Sobject), 10)

vfPlot1 <- VariableFeaturePlot(Sobject)
vfPLot2 <- LabelPoints(plot = vfPlot1, points = top10, repel = TRUE)
vfPlot1 + vfPLot2

#Scaling

all.genes <- rownames(Sobject)
Sobject <- ScaleData(Sobject, features = all.genes)

#PCA

Sobject <- RunPCA(Sobject, features = VariableFeatures(object = Sobject))
VizDimLoadings(Sobject, dims = 1:2, reduction = "pca") #1D plots with gene names
DimPlot(Sobject, reduction = "pca") + NoLegend() #2D plot
DimHeatmap(Sobject, dims = 1:15, cells = 500, balanced = TRUE) #15 dimensions, heat maps
ElbowPlot(Sobject, ndims = 50) #50 dimensions, dips around PC25?

#Clustering (25 dimensions)

Sobject <- FindNeighbors(Sobject, dims = 1:25)
Sobject <- FindClusters(Sobject, resolution = 0.1)
head(Idents(Sobject), 10)

#UMAP 

Sobject <- RunUMAP(Sobject, dims = 1:25)
DimPlot(Sobject, reduction = "umap")
saveRDS(Sobject, file = "/fp/homes01/u01/ec-evanmv/scRNAseq_GSE181919/Sobject.Rds")

#Id Cluster Biomarkers

Sobject.markers <- FindAllMarkers(Sobject, only.pos = TRUE)
saveRDS(Sobject.markers, file = "../scRNAseq_GSE181919/Sobject.markers.Rds")

FeaturePlot(Sobject, features = c("DCN", "COL1A2", "C1S", "COL6A2", "C1R"))

Sobject.markers %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% 
  ungroup() -> top10

ClusterHeatMap <- DoHeatmap(subset(Sobject, downsample = 100) , features = top10$gene) + NoLegend()

dev.off()

write_csv(Sobject.markers, "markers.csv")
write_csv(top10, "top10Markers.csv")
