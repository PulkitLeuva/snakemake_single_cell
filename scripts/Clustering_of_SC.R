########### Loading Libraries #####################
library(rhdf5)
library(hdf5r)
library(Seurat)
library(SeuratObject)
library(stringr)
library(ggplot2)
library(dplyr)

# Capture command-line arguments
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_meta_data <- args[2]
output_clustering <- args[3]
dims <- as.integer(args[4])
Cluster_resolutions = args[5]
# Read the RDS file
data_ <- readRDS(input_file)

# Perform clustering using the user-provided number of dimensions
Scaled_pbmc_PCA <- FindNeighbors(data_, dims = 1:dims)
Scaled_pbmc_PCA <- FindClusters(Scaled_pbmc_PCA, resolution = seq(0.05, 2, 0.05))

# Save the meta data and clustering results
Meta_Data <- Scaled_pbmc_PCA@meta.data
write.csv(Meta_Data, output_meta_data)
Idents(Scaled_pbmc_PCA)
saveRDS(Scaled_pbmc_PCA, output_clustering)


cluster_counts <- data.frame(Resolution = numeric(), Num_Clusters = integer())
Cluster_count = Scaled_pbmc_PCA@meta.data
# Loop through each resolution and count the number of clusters
for (res in seq(0.05,2,0.05)) {
  cluster_col <- paste0("RNA_snn_res.", res)
  num_clusters <- length(unique(Cluster_count[[cluster_col]]))
  cluster_counts <- rbind(cluster_counts, data.frame(Resolution = res, Num_Clusters = num_clusters))
}

write.csv(cluster_counts, Cluster_resolutions)
