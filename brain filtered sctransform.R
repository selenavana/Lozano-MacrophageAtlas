#clustering with the sctranform normalization
#setwd("Macintosh HD://Users/boisvertlab/Desktop/HHWNJBGXH/Brain")

library(Seurat)
Brain.data <- Read10X(data.dir = "~/Desktop/HHWNJBGXH/Brain/Brain/outs/raw_feature_bc_matrix")
Brain <- CreateSeuratObject(counts = Brain.data, project = "Brain", min.cells = 3, min.features = 200)
Brain

#filter cells
Brain[["percent.mt"]] <- PercentageFeatureSet(Brain, pattern = "^mt-")
VlnPlot(Brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
Brain <- subset(Brain, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(Brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)


Brain <- SCTransform(object = Brain, verbose = TRUE)
Brain <- RunPCA (object = Brain, verbose = TRUE)
DimPlot(Brain, reduction = "pca")

#did not visualize PCAs
#VizDimLoadings(Brain, dims = 1:2, reduction = "pca")
#DimHeatmap(Brain, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(Brain, dims = 1:12, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
#JackStraw cannot be used on SCTransformed data

ElbowPlot(Brain, ndims = 50, reduction = "pca")

Brain <- FindNeighbors(Brain, dims = 1:40)

#using the resolution value that outputs 15 clusters
Brain <- FindClusters(Brain, resolution = 0.43)
head(Idents(Brain),5)
Brain <- RunUMAP(object = Brain, dims = 1:40, verbose = TRUE)
DimPlot(Brain, reduction = "umap", label = TRUE)

Brain <- RunTSNE(object = Brain, dims = 1:40, verbose = TRUE)
DimPlot(Brain, reduction = "tsne", label = TRUE)

saveRDS(Brain, file = "20210511 Brain filtered SCtransform.rds")

#projecting each PC onto a UMAP - did PCs 1:40
#FeaturePlot(Brain, features = (paste0('PC_', 1:9)), reduction = "umap")
#FeaturePlot(Brain, features = (paste0('PC_', 10:18)), reduction = "umap")
#FeaturePlot(Brain, features = (paste0('PC_', 19:27)), reduction = "umap")
#print(Brain[["pca"]], dims = 1:30, nfeatures = 30)

#finding markers of differentially expressed genes in all clusters
all.markers <- FindAllMarkers(Brain)
write.csv(all.markers, "20210511 scTransform filtered all genes expression brain.csv")

#building dendrogram
#install.packages('ape')
Brain <- BuildClusterTree(Brain)
#PlotClusterTree(Brain)
data.tree <- Tool(Brain, slot = "BuildClusterTree")

#save plots as PDF
pdf("20210511 Brain plots sctransform filtered.pdf", width = 6.5, height = 5)
ElbowPlot(Brain, ndims = 50, reduction = "pca")
DimPlot(Brain, reduction = "umap", label = TRUE)
DimPlot(Brain, reduction = "tsne", label = TRUE)
ape::plot.phylo(x = data.tree, direction = "downwards", srt = 90, no.margin = TRUE, label.offset = 0.3, adj = 0.5, use.edge.length = FALSE)
FeaturePlot(Brain, features = (paste0('PC_', 1:4)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 5:8)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 9:12)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 13:16)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 17:20)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 21:24)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 25:28)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 29:32)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 33:36)), reduction = "umap")
FeaturePlot(Brain, features = (paste0('PC_', 37:40)), reduction = "umap")
dev.off()

saveRDS(Brain, file = "20210511 R Code SCtransform filtered umap tsne brain.rds")
readRDS(Brain, file = "20210511 R Code SCtransform filtered umap tsne brain.rds")


library(patchwork)
library(ggplot2)

M1dotplot <- DotPlot(Brain, assay = "SCT", features = c("Cd38","Ccr7","Cd80","Cd86","Cxcl9","Cybb","Fcgr1","Fpr2","Gpr18","H-2", "Il12a", "Il12b","Il1b", "Il23a","Il4ra","Il6","Marco", "Nos2", "Slamf1","Socs3","Stat1","Tnf")) +  RotatedAxis() + xlab("M1 Macrophage Markers") + ylab("Brain Cluster")
M1dotplot

M2dotplot <- DotPlot(Brain, assay = "SCT", features = c("Arg1", "Chil3", "Egr2", "Il27ra", "Mrc1", "Myc", "Pparg", "Retnla", "Stat6", "CCL17", "Clec7a", "Fn1", "Igf1", "Il10", "Il1r1", "Tgfb1", "Ccl1", "Cd86", "Il12a", "Il12b", "Il1b", "Il6", "Sphk1", "Tnf", "Tnfsf14", "Cd163", "Cxcl13", "Mertk", "Nos2", "Vegfa")) +RotatedAxis() + xlab("M2 Macrophage Markers") + ylab("Brain Cluster")
M2dotplot

Macmonodotplot <- DotPlot(Brain, assay = "SCT", features = c("Apoe","Ccr2","Cd209a","Cd36","Cd74","Cebpa","Cebpd","Csf1r","Cx3cr1","Fos","H2-Aa","H2-Ab1","Itgam","Itgax","Lgals3","Ly6c1","Ly6c2", "Lyz1", "Mef2a","Mmp8","Nr4a1","Pparg", "Akr1b10", "Blvrb", "Camk1", "Glul", "Myo7a", "Nln", "Pcyox1", "Pla2g15", "Pon3", "Slc48a1")) +  RotatedAxis() + xlab("Monocyte / Macrophage Markers") + ylab("Brain Cluster")
Macmonodotplot

DCCytTBCelldotplot <- DotPlot(Brain, assay = "SCT", features = c("Ciita", "Clec9a", "Cst3", "Fcer1a", "Itgae", "Itgam", "Itgax", "Ptprc", "Xcr1", "Il7r", "S100a4", "Cd3d", "Cd3e", "Cd3g", "Cd69", "Cd8a", "Cd8b1", "Itga2", "Cd19", "Cd79a", "Cd79b", "Cr2", "Ebf1", "Ighm", "Igkc", "Iglc1", "Iglc2", "Iglc3", "Ly6d", "Ms4a1","Gnly", "Nkg7",  "Ppbp")) +  RotatedAxis() + xlab("Dendritic Cell / Cytotoxic T Cell / B Cell Markers") + ylab("Brain Cluster")
DCCytTBCelldotplot


pdf("20211210 Brain dot plots 2.pdf", width = 7.5, height = 5)
M1dotplot
M2dotplot
Macmonodotplot
DCCytTBCelldotplot
dev.off()
