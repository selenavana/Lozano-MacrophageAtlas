library(Seurat)
library(ggplot2)
library(patchwork)

#Load data
Brain <- readRDS(file = "20210511 R Code SCtransform filtered umap tsne Brain.rds")

#Subsetting BAM clusters
Brain7 <- subset(Brain, idents = "7")
Brain9 <- subset(Brain, idents = "9")
Brain13 <- subset(Brain, idents = "13")

#Run PCA on each subset
Brain7 <- RunPCA(Brain7, verbose = FALSE)
ElbowPlot(Brain7, ndims = 50, reduction = "pca")

Brain9 <- RunPCA(Brain9, verbose = FALSE)
ElbowPlot(Brain9, ndims = 50, reduction = "pca")

Brain13 <- RunPCA(Brain13, verbose = FALSE)
ElbowPlot(Brain13, ndims = 50, reduction = "pca")

#Run UMAP and TSNE for each subset
Brain7 <- RunUMAP(Brain7, dims = 1:35)
DimPlot(Brain7, reduction = "umap", label = TRUE, repel = TRUE)
Brain7 <- RunTSNE(object = Brain7, dims = 1:35, verbose = FALSE)
DimPlot(Brain7, reduction = "tsne", label = TRUE, repel = TRUE)

Brain9 <- RunUMAP(Brain9, dims = 1:35)
DimPlot(Brain9, reduction = "umap", label = TRUE, repel = TRUE)
Brain9 <- RunTSNE(object = Brain9, dims = 1:35, verbose = FALSE)
DimPlot(Brain9, reduction = "tsne", label = TRUE, repel = TRUE)

Brain13 <- RunUMAP(Brain13, dims = 1:35)
DimPlot(Brain13, reduction = "umap", label = TRUE, repel = TRUE)
Brain13 <- RunTSNE(object = Brain13, dims = 1:35, verbose = FALSE)
DimPlot(Brain13, reduction = "tsne", label = TRUE, repel = TRUE)

#Run FindNeighbors and FindClusters for each subset
Brain7 <- FindNeighbors(Brain7, dims = 1:35)
Brain7 <- FindClusters(Brain7, resolution = 0.5)
DimPlot(Brain7, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(Brain7, reduction = "tsne", label = TRUE, repel = TRUE)

Brain9 <- FindNeighbors(Brain9, dims = 1:35)
Brain9 <- FindClusters(Brain9, resolution = 0.5)
DimPlot(Brain9, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(Brain9, reduction = "tsne", label = TRUE, repel = TRUE)

Brain13 <- FindNeighbors(Brain13, dims = 1:35)
Brain13 <- FindClusters(Brain13, resolution = 0.5)
DimPlot(Brain13, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(Brain13, reduction = "tsne", label = TRUE, repel = TRUE)

#finding markers of differentially expressed genes in all clusters, for each subset
all.markers7 <- FindAllMarkers(Brain7)
write.csv(all.markers7, "20220529 DEG Brain cluster 7 subclusters.csv")

all.markers9 <- FindAllMarkers(Brain9)
write.csv(all.markers9, "20220529 DEG Brain cluster 9 subclusters.csv")

all.markers13 <- FindAllMarkers(Brain13)
write.csv(all.markers13, "20220529 DEG Brain cluster 13 subclusters.csv")

#Rename idents

Brain7 <- RenameIdents(Brain7, "0" = "7_0", "1" = "7_1", "2" = "7_2")
Brain9 <- RenameIdents(Brain9, "0" = "9_0", "1" = "9_1")
Brain13 <- RenameIdents(Brain13, "0" = "13_0", "1" = "13_1")

pdf("20220529 brain BAM subcluster graphs.pdf", width = 8.5, height = 6)
DimPlot(Brain7, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Brain7, reduction = "tsne", label = FALSE, repel = TRUE)
DimPlot(Brain9, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Brain9, reduction = "tsne", label = FALSE, repel = TRUE)
DimPlot(Brain13, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Brain13, reduction = "tsne", label = FALSE, repel = TRUE)
dev.off()

#BAM identification with dotplots

BAM <- c(Brain7, Brain9, Brain13)
BAM_Identification_Dotplot_allclusters <- DotPlot(Brain, assay = "SCT", features = c("S100a4", "Csf1r", "Cx3cr1", "Csf1", "Apoe", "Cd74", "Klra2", "Mrc1", "Lyve1", "Clec4n", "H2-Aa", "H2-Ab1", "H2-Eb1", "Clec12a", "Axl", "Maf", "Folr2", "Igf1", "Fcgr2b")) +  RotatedAxis() + xlab("Microglia/BAM Markers") + ylab("Brain Cluster")
BAM_Identification_Dotplot_allclusters

p7 <- DotPlot(Brain7, assay = "SCT", features = c("S100a4", "Csf1r", "Cx3cr1", "Csf1", "Apoe", "Cd74", "Klra2", "Mrc1", "Lyve1", "Clec4n", "H2-Aa", "H2-Ab1", "H2-Eb1", "Clec12a", "Axl", "Maf", "Folr2", "Igf1", "Fcgr2b")) +  RotatedAxis() + xlab(" ") + ylab(" ") + theme(legend.key.size = unit(0.3, "cm"))
p9 <- DotPlot(Brain9, assay = "SCT", features = c("S100a4", "Csf1r", "Cx3cr1", "Csf1", "Apoe", "Cd74", "Klra2", "Mrc1", "Lyve1", "Clec4n", "H2-Aa", "H2-Ab1", "H2-Eb1", "Clec12a", "Axl", "Maf", "Folr2", "Igf1", "Fcgr2b")) +  RotatedAxis() + ylab("BAM Cluster") + xlab(" ") + theme(legend.key.size = unit(0.3, "cm"))
p13 <- DotPlot(Brain13, assay = "SCT", features = c("S100a4", "Csf1r", "Cx3cr1", "Csf1", "Apoe", "Cd74", "Klra2", "Mrc1", "Lyve1", "Clec4n", "H2-Aa", "H2-Ab1", "H2-Eb1", "Clec12a", "Axl", "Maf", "Folr2", "Igf1", "Fcgr2b")) +  RotatedAxis() + xlab("Microglia/BAM Markers") + ylab(" ") + theme(legend.key.size = unit(0.3, "cm"))


p7 / p9 / p13

#BAM Marker Violin plot
features <- c("Apoe", "Axl", "Cd74", "Cd163", "Clec12a", "Igf1", "Lyve1", "Lyz2", "Mrc1", "Ms4a7", "Cybb", "Ifi27l2a", "Wfdc17")
plot1 <- VlnPlot(Brain, features, stack = TRUE, flip = TRUE, assay = "SCT") + theme(legend.position = "none")
plot1


#BAM vs other brain cells DEG
all.markers.bam <- FindMarkers(Brain, ident.1 = c("7", "9", "13"), ident.2 = NULL)
write.csv(all.markers.bam, "20220713 DEG for BAM.csv")
