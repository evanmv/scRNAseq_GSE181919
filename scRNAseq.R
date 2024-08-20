library(tidyverse)
remotes::install_version('Matrix', version = '1.6.4') #Required for installing Seurat
install.packages("Seurat")
library(SeuratObject)
library(dplyr)
library(patchwork)

#Read in UMI count date
UMIcounts <- read_delim('GSE181919_UMI_counts.txt') #Didn't work, too large I think. 
UMI_counts <- data.table::fread('GSE181919_UMI_counts.txt')

head(UMI_counts[, 1:10])
rownames(UMI_counts) <- UMI_counts$V1 #Set gene names from column V1 as row names
UMI.t <- as.tibble(UMI_counts) #Tidy data

#Create Seurat object
?Seurat
Sobject <- CreateSeuratObject(counts = UMI_counts, project = "Bish", min.features = 200, min.cells = 3)
