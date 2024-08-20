library(tidyverse)
install.packages("Seurat")
library(Seurat)

#Read in UMI count date
UMIcounts <- read_delim('GSE181919_UMI_counts.txt') #Didn't work, too large I think. 
UMI_counts <- data.table::fread('GSE181919_UMI_counts.txt')

head(UMI_counts[, 1:10])
rownames(UMI_counts) <- UMI_counts$V1 #Set gene names from column V1 as row names
UMI.t <- as.tibble(UMI_counts) #Tidy data

#Create Seurat object