#clustering with the sctranform normalization
#setwd("Macintosh HD://Users/boisvertlab/Desktop/HHWNJBGXH/Lung")

library(Seurat)
Lung.data <- Read10X(data.dir = "~/Desktop/HHWNJBGXH/Lung/Lung/outs/raw_feature_bc_matrix")
Lung <- CreateSeuratObject(counts = Lung.data, project = "Lung", min.cells = 3, min.features = 200)
Lung

#filter cells
Lung[["percent.mt"]] <- PercentageFeatureSet(Lung, pattern = "^mt-")
VlnPlot(Lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
Lung <- subset(Lung, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(Lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

Lung <- SCTransform(object = Lung, verbose = TRUE)
Lung <- RunPCA (object = Lung, verbose = TRUE)
#DimPlot(Lung, reduction = "pca")

#did not visualize PCAs
#VizDimLoadings(Lung, dims = 1:2, reduction = "pca")
#DimHeatmap(Lung, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(Lung, dims = 1:12, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
#JackStraw cannot be used on SCTransformed data

ElbowPlot(Lung, ndims = 50, reduction = "pca")

Lung <- FindNeighbors(Lung, dims = 1:40)

#using the resolution value that outputs 15 clusters
Lung <- FindClusters(Lung, resolution = 0.43)
head(Idents(Lung),5)
Lung <- RunUMAP(object = Lung, dims = 1:40, verbose = TRUE)
DimPlot(Lung, reduction = "umap", label = TRUE)

Lung <- RunTSNE(object = Lung, dims = 1:40, verbose = TRUE)
DimPlot(Lung, reduction = "tsne", label = TRUE)

saveRDS(Lung, file = "20210517 Lung filtered SCtransform.rds")

#projecting each PC onto a UMAP - did PCs 1:40
#FeaturePlot(Lung, features = (paste0('PC_', 1:9)), reduction = "umap")
#FeaturePlot(Lung, features = (paste0('PC_', 10:18)), reduction = "umap")
#FeaturePlot(Lung, features = (paste0('PC_', 19:27)), reduction = "umap")
#print(Lung[["pca"]], dims = 1:30, nfeatures = 30)

#finding markers of differentially expressed genes in all clusters
all.markers <- FindAllMarkers(Lung)
write.csv(all.markers, "20210517 scTransform filtered all genes expression Lung.csv")

#building dendrogram
#install.packages('ape')
Lung <- BuildClusterTree(Lung)
#PlotClusterTree(Lung)
data.tree <- Tool(Lung, slot = "BuildClusterTree")

#save plots as PDF
pdf("20210517 Lung plots sctransform filtered.pdf", width = 6.5, height = 5)
ElbowPlot(Lung, ndims = 50, reduction = "pca")
DimPlot(Lung, reduction = "umap", label = TRUE)
DimPlot(Lung, reduction = "tsne", label = TRUE)
ape::plot.phylo(x = data.tree, direction = "downwards", srt = 90, no.margin = TRUE, label.offset = 0.3, adj = 0.5, use.edge.length = FALSE)
FeaturePlot(Lung, features = (paste0('PC_', 1:4)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 5:8)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 9:12)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 13:16)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 17:20)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 21:24)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 25:28)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 29:32)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 33:36)), reduction = "umap")
FeaturePlot(Lung, features = (paste0('PC_', 37:40)), reduction = "umap")
dev.off()

saveRDS(Lung, file = "20210517 R Code SCtransform filtered umap tsne Lung.rds")
# readRDS(Lung, file = "20210517 R Code SCtransform filtered umap tsne Lung.rds")

library(patchwork)
library(ggplot2)

lung.plot <- DotPlot(Lung, assay = "SCT", features = c("Adgre1","Atp10a","Blnk","C1qa","C1qb","C1qc","C4b","Car4","Ccl12","Ccl2","Ccl24","Ccl3","Ccl4","Ccl6","Ccl7","Ccl8","Ccl9","Ccr1","Ccr2","Ccr5","Ccr7","Cd14","Cd163","Cd36","Cd38","Cd4","Cd68","Cd86","Chil3","Colec12","Csf1r","Csf2ra","Csf2rb","Cx3cr1","Cxcl13","Cxcl14","Cxcl16","Cxcl2","Dpp4","Enpep","Fcgr1","Fcgr3","Fcna","Flvcr2","Folr2","Gas7","H2-Ab1","Hpgd","Ifnar1","Ifnar2","Ifngr1","Ifngr2","Il10","Il10ra","Il10rb","Il17ra","Il18","Il1b","Il21r","Il4ra","Il6ra","Itga9","Itgae","Itgam","Itgax","Junb","Krt79","Lacc1","Lamp1","Ly6g","Lyve1","Lyz2","Maf","Mafb","Marco","Mertk","Mrc1","Olfml3","Retnla","Scimp","Siglecf","Socs3","Tgm2","Tlr2","Tlr4","Tnf","Ttc9","Zbtb46")) +  RotatedAxis()

lung.plot

lungvlnplots <- VlnPlot(Lung, features = c("Itgam", "Siglecf", "Olfml3", "Gas7"), pt.size = 0, combine = FALSE, assay = "SCT")


M1dotplot <- DotPlot(Lung, assay = "SCT", features = c("Cd38","Ccr7","Cd80","Cd86","Cxcl9","Cybb","Fcgr1","Fpr2","Gpr18","H-2", "Il12a", "Il12b","Il1b", "Il23a","Il4ra","Il6","Marco", "Nos2", "Slamf1","Socs3","Stat1","Tnf")) +  RotatedAxis() + xlab("M1 Macrophage Markers") + ylab("Lung Cluster")
M1dotplot

M2dotplot <- DotPlot(Lung, assay = "SCT", features = c("Arg1", "Chil3", "Egr2", "Il27ra", "Mrc1", "Myc", "Pparg", "Retnla", "Stat6", "CCL17", "Clec7a", "Fn1", "Igf1", "Il10", "Il1r1", "Tgfb1", "Ccl1", "Cd86", "Il12a", "Il12b", "Il1b", "Il6", "Sphk1", "Tnf", "Tnfsf14", "Cd163", "Cxcl13", "Mertk", "Nos2", "Vegfa")) +RotatedAxis() + xlab("M2 Macrophage Markers") + ylab("Lung Cluster")
M2dotplot

Macmonodotplot <- DotPlot(Lung, assay = "SCT", features = c("Apoe","Ccr2","Cd209a","Cd36","Cd74","Cebpa","Cebpd","Csf1r","Cx3cr1","Fos","H2-Aa","H2-Ab1","Itgam","Itgax","Lgals3","Ly6c1","Ly6c2", "Lyz1", "Mef2a","Mmp8","Nr4a1","Pparg", "Akr1b10", "Blvrb", "Camk1", "Glul", "Myo7a", "Nln", "Pcyox1", "Pla2g15", "Pon3", "Slc48a1")) +  RotatedAxis() + xlab("Monocyte / Macrophage Markers") + ylab("Lung Cluster")
Macmonodotplot

DCCytTBCelldotplot <- DotPlot(Lung, assay = "SCT", features = c("Ciita", "Clec9a", "Cst3", "Fcer1a", "Itgae", "Itgam", "Itgax", "Ptprc", "Xcr1", "Il7r", "S100a4", "Cd3d", "Cd3e", "Cd3g", "Cd69", "Cd8a", "Cd8b1", "Itga2", "Cd19", "Cd79a", "Cd79b", "Cr2", "Ebf1", "Ighm", "Igkc", "Iglc1", "Iglc2", "Iglc3", "Ly6d", "Ms4a1","Gnly", "Nkg7",  "Ppbp")) +  RotatedAxis() + xlab("Dendritic Cell / Cytotoxic T Cell / B Cell Markers") + ylab("Lung Cluster")
DCCytTBCelldotplot


pdf("20211228 Lung dot plots 2.pdf", width = 9, height = 5)
M1dotplot
M2dotplot
Macmonodotplot
DCCytTBCelldotplot
dev.off()
