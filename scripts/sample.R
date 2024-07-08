
library(Seurat)
library(rhdf5)
library(hdf5r)
library(SeuratObject)
library(stringr)
library(ggplot2)
library(dplyr)

############### Loading Files ###################################################
scRNAseq_Data_10x=Read10X_h5((commandArgs(trailingOnly=TRUE)[1]),use.names = TRUE)


####### Create Seurta Object 
data <- CreateSeuratObject(counts = scRNAseq_Data_10x, project = "scRNA_Seq", min.cells = commandArgs(trailingOnly=TRUE)[3], min.features = commandArgs(trailingOnly=TRUE)[4])


###############################  barplot for the number of Genes ###############
pdf((commandArgs(trailingOnly=TRUE)[5]),height=5.5,width=4)
data@meta.data %>% 
  ggplot(aes(x=orig.ident),fill=orig.ident) + 
  geom_bar() +
  theme_classic() +
  theme(text = element_text(size=15),axis.title.x = element_blank()) +
  theme(plot.title = element_text(hjust=0.5, face="bold"))+
  labs(title = "Number of cells per group",
       subtitle = "Determined by number of unique cellular barcodes detected")+
  theme(plot.title = element_text(color = "steelblue4", size = 13, face = "bold"),
  plot.subtitle = element_text(color = "orange", size = 9),
  plot.caption = element_text(color = "green", size = 6, face = "italic")) + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), plot.caption = element_text(hjust = 1))


############### Mitochondrial percentage of Cells ##################################
data[["percentage.mt"]]=PercentageFeatureSet(data,pattern = "^MT-")
#data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

pdf((commandArgs(trailingOnly=TRUE)[6]),height=6.5,width=6)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percentage.mt"), ncol = 3,cols = "violet")
dev.off()


##############   Scatter Plot ############################
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percentage.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf((commandArgs(trailingOnly=TRUE)[7]),height=8,width=5)
plot1 + plot2
dev.off()

saveRDS(data,commandArgs(trailingOnly=TRUE)[2])
############## Filtering Cells ##############################


#pbmc <- subset(pbmc, subset = nFeature_RNA > number & nFeature_RNA <2000 & percent.mt < 25)
