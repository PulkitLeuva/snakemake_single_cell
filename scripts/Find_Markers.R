
############### Loading Libraries #####################
library(rhdf5)
library(hdf5r)
library(Seurat)
library(SeuratObject)
library(stringr)
library(ggplot2)
library(dplyr)



# Extract command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_file <- args[1]
output_file <- args[2]
output_file2 <- args[3]
Resolutions <- as.character(args[4])
Onlypositive_FC_ <- as.logical(args[5])
min_pct <- as.numeric(args[6])
Top_listed_Genes <- as.integer(args[7])
pdf1 <- args[8]
pdf2 <- args[9]
dims <- as.numeric(args[10])

# Load the Seurat object
data_ <- readRDS(input_file)

################ Set the resolution for clusters ##################
Idents(object = data_) <- as.character(paste0("RNA_snn_res.",Resolutions))

data_ <- RunUMAP(data_, reduction = "pca", dims = 1:dims)
pdf(pdf2, height = 6, width = 5)

DimPlot(data_,
        reduction = "umap",
        label = TRUE,
        label.size = 6, pt.size = 1) + labs(title = "Uniform Manifold Approximation and Projection (UMAP)") +
  theme(plot.title = element_text(color = "steelblue4", size = 13, face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

################ Finding markers of all Genes #####################
combined_markers <- FindAllMarkers(object = data_, only.pos = Onlypositive_FC_, min.pct = min_pct, logfc.threshold = min_pct, min.cells.group = 3)
combined_markers
write.table(combined_markers, output_file, sep = "\t", quote = FALSE, col.names = NA)

################   Top Markers   ##################################

top <- combined_markers %>% group_by(cluster) %>% top_n(n = Top_listed_Genes, wt = avg_log2FC)
top
write.table(top, output_file2, sep = "\t", quote = FALSE, col.names = NA)

###################  Heatmap #####################################

pdf(pdf1, height = 6, width = 5)
DoHeatmap(data_, features = top$gene, label = TRUE, group.bar = TRUE, size = 5) +
  theme(text = element_text(size = 2), legend.position = "right", plot.margin = margin(2, 2, 2, 2, "cm"),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10), axis.text = element_text(size = 4))

dev.off()
