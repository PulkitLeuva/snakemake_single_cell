########### Loading Libraries #####################
library(rhdf5)
library(hdf5r)
library(Seurat)
library(SeuratObject)
library(stringr)
library(ggplot2)
library(dplyr)

data= readRDS((commandArgs(trailingOnly=TRUE)[1]))

data <- subset(data, subset = nFeature_RNA > (commandArgs(trailingOnly=TRUE)[3]) & nFeature_RNA <(commandArgs(trailingOnly=TRUE)[4]) & percentage.mt < (commandArgs(trailingOnly=TRUE)[5]))

Normalized_pbmc=NormalizeData(data)

Normalized_pbmc <- FindVariableFeatures(Normalized_pbmc, selection.method = "vst", nfeatures = as.numeric(commandArgs(trailingOnly=TRUE)[6]))

#Normalized_pbmc@assays$RNA@counts
# Normalized_pbmc@assays$RNA@data
top10 <- head(VariableFeatures(Normalized_pbmc), as.numeric((commandArgs(trailingOnly=TRUE)[7])))


plot1 <- VariableFeaturePlot(Normalized_pbmc,cols=c("black","blue"))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)

pdf((commandArgs(trailingOnly=TRUE)[8]),height=6,width=5)
plot2
dev.off()

################# sacling data
all_genes=rownames(Normalized_pbmc)
Scaled_pbmc=ScaleData(Normalized_pbmc,features = all_genes)

pdf((commandArgs(trailingOnly=TRUE)[9]),height=6.5,width=6)
VlnPlot(Scaled_pbmc, features = c("nFeature_RNA", "nCount_RNA", "percentage.mt"), ncol = 3,cols = "violet")
dev.off()

plot1 <- VariableFeaturePlot(Scaled_pbmc,cols=c("black","blue"))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)

pdf((commandArgs(trailingOnly=TRUE)[10]),height=8,width=5)
plot2
dev.off()


##################### Performing LInear Dimensionality of The Dataset #############

Scaled_pbmc_PCA=RunPCA(Scaled_pbmc,features = VariableFeatures(object = Scaled_pbmc))
Scaled_pbmc_PCA@reductions$pca

pdf((commandArgs(trailingOnly=TRUE)[11]),height=8,width=5)
DimPlot(Scaled_pbmc_PCA, reduction = "pca", cols = "magenta") + NoLegend()
dev.off()

########################3 Heatmap of PCA #######################################
pdf((commandArgs(trailingOnly=TRUE)[12]),height=7,width=5)
DimHeatmap(Scaled_pbmc_PCA, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()
############### Elbow olot for determining the PC's
pdf((commandArgs(trailingOnly=TRUE)[13]),height=5.5,width=5)
ElbowPlot(Scaled_pbmc_PCA)
dev.off()

saveRDS(Scaled_pbmc_PCA,commandArgs(trailingOnly=TRUE)[2])
############################ Finding Clusters of cells ##########################

# Scaled_pbmc_PCA <- FindNeighbors(Scaled_pbmc_PCA, dims = 1:10)
# Scaled_pbmc_PCA@tools
# Scaled_pbmc_PCA <- FindClusters(Scaled_pbmc_PCA, resolution = 0.25)

# Idents(Scaled_pbmc_PCA)

