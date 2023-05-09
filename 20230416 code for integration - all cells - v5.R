library(Seurat)
library(ggplot2)
library(patchwork)

#Load data
Heart <- readRDS(file = "20210518 R Code SCtransform filtered umap tsne Heart.rds")
Liver <- readRDS(file = "20210517 R Code SCtransform filtered umap tsne Liver.rds")
Brain <- readRDS(file = "20210511 R Code SCtransform filtered umap tsne Brain.rds")
Lung <- readRDS(file = "20210517 R Code SCtransform filtered umap tsne Lung.rds")

#save old cluster idents
Heart[["original.cluster"]] <- Idents(Heart)
Brain[["original.cluster"]] <- Idents(Brain)
Liver[["original.cluster"]] <- Idents(Liver)
Lung[["original.cluster"]] <- Idents(Lung)


#renaming clusters
{Heart <- RenameIdents(Heart, "0" = "Monocyte Derived", "4" = "Other", "6" = "Monocyte Derived", "7" = "Monocyte Derived", 
                      "1" = "Cycling", "9" = "Other", "10" = "Other", "11" = "Other", "14" = "Other", 
                      "2" = "Dendritic Cells", "3" = "Dendritic Cells", 
                      "5" = "Cycling", "12" = "Pro-inflammatory", "8" = "Heart Tissue Resident", "13" = "Anti-inflammatory")
}

{Brain <- RenameIdents(Brain, "0" = "Microglia", "1" = "Microglia", "2" = "Microglia", "3" = "Other", "4" = "Cycling", "5" = "Microglia", 
                       "6" = "Cycling", "7" = "Border-Associated", "8" = "Microglia", "9" = "Border-Associated",  "10" = "Microglia", 
                       "11" = "Microglia", "12" = "Other", "13" = "Border-Associated", "14" = "Other", "15" = "Other")
}

{Liver <- RenameIdents(Liver, "0" = "Kupffer Cells", "1" = "Kupffer Cells", "2" = "Monocyte Derived", "3" = "Monocyte Derived", "4" = "Other", 
                       "5" = "Pro-inflammatory", "6" = "Anti-inflammatory", "7" = "Lymphocytes", "8" = "Other", "9" = "Dendritic Cells",  "10" = "Lymphocytes", "11" = "Peritoneal", "12" = "Other")

}  

{Lung <- RenameIdents(Lung, "0" = "Alveolar", "1" = "Monocyte Derived", "2" = "Cycling", "3" = "Pro-inflammatory", "4" = "Interstitial", "5" = "Interstitial", 
                        "6" = "Other", "7" = "Monocyte Derived", "8" = "Pro-inflammatory", "9" = "Monocyte Derived",  "10" = "Alveolar", "11" = "Cycling", "12" = "Lymphocytes", "13" = "Other", 
                        "14" = "Alveolar", "15" = "Interstitial", "16" = "Other", "17" = "Lymphocytes", "18" = "Other", "19" = "Interstitial")
}   

#integration
#Heart <- subset(Heart, idents = "Other", invert = TRUE)
#Brain <- subset(Brain, idents = "Other", invert = TRUE)
#Liver <- subset(Liver, idents = c("Other", "Lymphocytes"), invert = TRUE)
#Lung <- subset(Lung, idents = c("Other", "Lymphocytes"), invert = TRUE)

macrophage.list <- c(Heart, Brain, Liver, Lung)
macrophage.features <- SelectIntegrationFeatures(object.list = macrophage.list, nfeatures = 3000)
options(future.globals.maxSize = 8000 * 1024^2)
macrophage.list <- PrepSCTIntegration(object.list = macrophage.list, anchor.features = macrophage.features, verbose = FALSE)
macrophage.anchors <- FindIntegrationAnchors(object.list = macrophage.list, normalization.method = "SCT", anchor.features = macrophage.features, verbose = FALSE)
macrophage.integrated <- IntegrateData(anchorset = macrophage.anchors, normalization.method = "SCT", verbose = TRUE)

#PCA, tsne, umap
macrophage.integrated <- RunPCA(macrophage.integrated, verbose = FALSE)
ElbowPlot(macrophage.integrated, ndims = 50, reduction = "pca")
macrophage.integrated <- RunUMAP(macrophage.integrated, dims = 1:35)
DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE)
macrophage.integrated <- RunTSNE(object = macrophage.integrated, dims = 1:35, verbose = FALSE)
DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE)

#save old idents
macrophage.integrated[["cell.ident"]] <- Idents(macrophage.integrated)

#clustering 
macrophage.integrated <- FindNeighbors(macrophage.integrated, dims = 1:35)
macrophage.integrated <- FindClusters(macrophage.integrated, resolution = 0.43)
DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE)

#save new idents
macrophage.integrated[["agg.ident"]] <- Idents(macrophage.integrated)

#SAVE RDS!!!
saveRDS(macrophage.integrated, file = "20230416 macrophage integrated V5 - ALL cells.rds")

#Start here!!!

macrophage.integrated <- readRDS(file = "20230416 macrophage integrated V5 - ALL cells.rds")

#clusterplot


source("custom_seurat_functions.R")
plot_integrated_clusters(macrophage.integrated)

source("custom_seurat_functions V2.R")
plot_integrated_clusters2(macrophage.integrated)

#plots into PDF, highlight each organ
pdf("20230321 labeled aggr plots - ALL cells.pdf", width = 8.5, height = 6)
plot_integrated_clusters(macrophage.integrated)
plot_integrated_clusters2(macrophage.integrated)
DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE)
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", label = TRUE, repel = TRUE) + ggtitle(" ") + theme(legend.position = c(0.8, 0.8), legend.text = element_text(size = 10), legend.key.size = unit(1,'mm'))
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", label = TRUE, repel = TRUE) + ggtitle(" ") + theme(legend.position = c(0.82, 0.9), legend.text = element_text(size = 9), legend.key.size = unit(1,'mm'))
DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident") + ggtitle("Organ")
DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE, group.by = "orig.ident") + ggtitle("Organ")
DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE, group.by = "original.cluster") + ggtitle("Original cluster")
DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE, group.by = "original.cluster") + ggtitle("Original cluster")
DimPlot(macrophage.integrated, reduction = "umap", cols = c("red", "gray", "gray", "gray"), group.by = "orig.ident") + ggtitle("Brain Cells")
DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", cols = c("gray", "red", "gray", "gray")) + ggtitle("Heart Cells")
DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", cols = c("gray",  "gray", "red", "gray")) + ggtitle("Liver Cells")
DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", cols = c("gray", "gray", "gray", "red")) + ggtitle("Lung Cells")
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "orig.ident", cols = c("red", "gray", "gray", "gray")) + ggtitle("Brain Cells")
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "orig.ident", cols = c("gray", "red", "gray", "gray")) + ggtitle("Heart Cells")
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "orig.ident", cols = c("gray",  "gray", "red", "gray")) + ggtitle("Liver Cells")
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "orig.ident", cols = c("gray", "gray", "gray", "red")) + ggtitle("Lung Cells")
dev.off()

#for 20230416 pdf
pdf("20230416 labeled aggr plots - ALL cells - new v5.pdf", width = 8.5, height = 6)
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", label = TRUE, repel = TRUE) + ggtitle(" ") + theme(legend.position = c(0.06, 0.24), legend.text = element_text(size = 10), legend.key.size = unit(1,'mm'))
DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE) + ggtitle(" ") + theme(legend.position = c(0.9, 0.15))
DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE)
plot_integrated_clusters(macrophage.integrated)
plot_integrated_clusters2(macrophage.integrated)
dev.off()


count.table1 <- table(macrophage.integrated@meta.data$seurat_clusters, macrophage.integrated@meta.data$cell.ident)
count.table1

count.table2 <- table(macrophage.integrated@meta.data$orig.ident, macrophage.integrated@meta.data$cell.ident)
count.table2

count.table3 <- table(macrophage.integrated@meta.data$seurat_clusters, macrophage.integrated@meta.data$orig.ident)
count.table3

count.total <- table(macrophage.integrated@meta.data$cell.ident)
count.total
