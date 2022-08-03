#clustering with the sctranform normalization
#setwd("Macintosh HD://Users/boisvertlab/Desktop/HHWNJBGXH/Liver")

library(Seurat)
Liver.data <- Read10X(data.dir = "~/Desktop/HHWNJBGXH/Liver/Liver/outs/raw_feature_bc_matrix")
Liver <- CreateSeuratObject(counts = Liver.data, project = "Liver", min.cells = 3, min.features = 200)
Liver

#filter cells
Liver[["percent.mt"]] <- PercentageFeatureSet(Liver, pattern = "^mt-")
VlnPlot(Liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
Liver <- subset(Liver, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(Liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

Liver <- SCTransform(object = Liver, verbose = TRUE)
Liver <- RunPCA (object = Liver, verbose = TRUE)
#DimPlot(Liver, reduction = "pca")

#did not visualize PCAs
#VizDimLoadings(Liver, dims = 1:2, reduction = "pca")
#DimHeatmap(Liver, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(Liver, dims = 1:12, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
#JackStraw cannot be used on SCTransformed data

ElbowPlot(Liver, ndims = 50, reduction = "pca")

Liver <- FindNeighbors(Liver, dims = 1:40)

#using the resolution value that outputs 15 clusters
Liver <- FindClusters(Liver, resolution = 0.43)
head(Idents(Liver),5)
Liver <- RunUMAP(object = Liver, dims = 1:40, verbose = TRUE)
DimPlot(Liver, reduction = "umap", label = TRUE)

Liver <- RunTSNE(object = Liver, dims = 1:40, verbose = TRUE)
DimPlot(Liver, reduction = "tsne", label = TRUE)

saveRDS(Liver, file = "20210517 Liver filtered SCtransform.rds")

#projecting each PC onto a UMAP - did PCs 1:40
#FeaturePlot(Liver, features = (paste0('PC_', 1:9)), reduction = "umap")
#FeaturePlot(Liver, features = (paste0('PC_', 10:18)), reduction = "umap")
#FeaturePlot(Liver, features = (paste0('PC_', 19:27)), reduction = "umap")
#print(Liver[["pca"]], dims = 1:30, nfeatures = 30)

#finding markers of differentially expressed genes in all clusters
all.markers <- FindAllMarkers(Liver)
write.csv(all.markers, "20210517 scTransform filtered all genes expression Liver.csv")

#building dendrogram
#install.packages('ape')
Liver <- BuildClusterTree(Liver)
#PlotClusterTree(Liver)
data.tree <- Tool(Liver, slot = "BuildClusterTree")

#save plots as PDF
pdf("20210511 Liver plots sctransform filtered.pdf", width = 6.5, height = 5)
ElbowPlot(Liver, ndims = 50, reduction = "pca")
DimPlot(Liver, reduction = "umap", label = TRUE)
DimPlot(Liver, reduction = "tsne", label = TRUE)
ape::plot.phylo(x = data.tree, direction = "downwards", srt = 90, no.margin = TRUE, label.offset = 0.3, adj = 0.5, use.edge.length = FALSE)
FeaturePlot(Liver, features = (paste0('PC_', 1:4)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 5:8)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 9:12)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 13:16)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 17:20)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 21:24)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 25:28)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 29:32)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 33:36)), reduction = "umap")
FeaturePlot(Liver, features = (paste0('PC_', 37:40)), reduction = "umap")
dev.off()

saveRDS(Liver, file = "20210517 R Code SCtransform filtered umap tsne Liver.rds")
# readRDS(Liver, file = "20210517 R Code SCtransform filtered umap tsne Liver.rds")

library(patchwork)
library(ggplot2)

M1dotplot <- DotPlot(Liver, assay = "SCT", features = c("Cd38","Ccr7","Cd80","Cd86","Cxcl9","Cybb","Fcgr1","Fpr2","Gpr18","Il12a", "Il12b","Il1b", "Il23a","Il4ra","Il6","Marco","Mhc2", "Nos2", "Slamf1","Socs3","Stat1","Tnf")) +  RotatedAxis() + xlab("M1 Macrophage Markers") + ylab("Liver Cluster")
M1dotplot

M2dotplot <- DotPlot(Liver, assay = "SCT", features = c("Arg1","Ccl1", "Ccl17","Cd163","Cd86","Chil3", "Clec7a", "Cxcl13","Egr2","Fn1","Igf1","Il10","Il12a","Il12b","Il1b","Il1r1","Il27ra","Il6","Mertk","Mrc1","Myc", "Nos2", "Pparg", "Retnla","Sphk1","Stat6","Tgfb1","Tnf","Tnfsf14","Vegfa")) +RotatedAxis() + xlab("M2 Macrophage Markers") + ylab("Liver Cluster")
M2dotplot

Macmonodotplot <- DotPlot(Liver, assay = "SCT", features = c("A930039a15Rik", "Apoe","Ccr2","Cd209a","Cd36","Cd74","Cebpa","Cebpd","Csf1r","Cx3cr1","Fos","H2-Aa","H2-Ab1","Itgam","Itgax","Lgals3","Ly6c1","Mef2a","Mmp8","Nr4a1","Pparg", "Akr1b10", "Blvrb", "Camk1", "Glul", "Myo7a", "Nln", "Pcyox1", "Pla2g15", "Pon3", "Slc48a1")) +  RotatedAxis() + xlab("Monocyte / Macrophage Markers") + ylab("Liver Cluster")
Macmonodotplot

DCCytTBCelldotplot <- DotPlot(Liver, assay = "SCT", features = c("Clec9a", "Itgae", "Itgax", "Ptprc", "Xcr1", "Cd3d", "Cd3e", "Cd3g", "Cd69", "Cd8a", "Cd8b1", "Cd19", "Cd79a", "Cd79b", "Ebf1", "Ighm", "Igkc", "Iglc1", "Iglc2", "Iglc3", "Ly6d", "Ms4a1")) +  RotatedAxis() + xlab("Dendritic Cell / Cytotoxic T Cell / B Cell Markers") + ylab("Liver Cluster")
DCCytTBCelldotplot

Kupfferdotplot <- DotPlot(Liver, assay = "SCT", features = c("Clec4f", "Vsig4", "Cd5l", "C1qa", "C1qb", "C1qc", "Mrc1", "Timd4", "Adgre1", "Cx3cr1", "Itgax", "Ly6g")) +  RotatedAxis() + xlab("Kupffer Cell Markers") + ylab("Liver Cluster")
Kupfferdotplot

MDMdotplot <- DotPlot(Liver, assay = "SCT", features = c("Ly6c2", "Clec4f", "Csf1r", "Clec4d", "Clec4e", "Clec5a", "Havcr2", "Tlr2")) +  RotatedAxis() + xlab("Monocyte-derived Macrophage Markers") + ylab("Liver Cluster")
MDMdotplot

DCdotplot <- DotPlot(Liver, assay = "SCT", features = c("H2", "Bst2", "Itgax", "Clec9a", "Cd209a", "Siglech", "Sirpa", "Cd24a", "Btla", "Ly75", "Cst3")) +  RotatedAxis() + xlab("Dendritic Cell Markers") + ylab("Liver Cluster")
DCdotplot

EosinNeutrodotplot <- DotPlot(Liver, assay = "SCT", features = c("Il5ra", "Ccr3", "Siglecf", "Fcgr2b", "Csf2ra", "Itgam", "Fcgr3", "Ly6g", "Itgb2")) +  RotatedAxis() + xlab("Eosinophil / Neutrophil Markers") + ylab("Liver Cluster")
EosinNeutrodotplot

pdf("20211202 liver dotplots elia.pdf", width = 6.5, height = 5)
Kupfferdotplot
MDMdotplot
DCdotplot
EosinNeutrodotplot
dev.off()

SunsMoMFdotplot <- DotPlot(Liver, assay = "SCT", features = c("Eng", "Cd36", "Nr3c2", "Fcgr2b", "Stab1", "Stab2", "Scarb1", "Scarf1", "Lyve1", "Lrp1", "Aoc3", "Prox1", "Cxcl9", "Cxcl11", "Cxcl1", "Cxcl16", "Clec4g", "Spn", "Fcgr1", "Ccr5", "Tlr2", "H2")) +  RotatedAxis() + xlab("Sinusoidal / MoMF Markers") + ylab("Liver Cluster")
SunsMoMFdotplot


#Common dot plots
library(Seurat)
library(patchwork)
library(ggplot2)

M1dotplot <- DotPlot(Liver, assay = "SCT", features = c("Cd38","Ccr7","Cd80","Cd86","Cxcl9","Cybb","Fcgr1","Fpr2","Gpr18","H-2", "Il12a", "Il12b","Il1b", "Il23a","Il4ra","Il6","Marco", "Nos2", "Slamf1","Socs3","Stat1","Tnf")) +  RotatedAxis() + xlab("M1 Macrophage Markers") + ylab("Liver Cluster")
M1dotplot

M2dotplot <- DotPlot(Liver, assay = "SCT", features = c("Arg1", "Chil3", "Egr2", "Il27ra", "Mrc1", "Myc", "Pparg", "Retnla", "Stat6", "CCL17", "Clec7a", "Fn1", "Igf1", "Il10", "Il1r1", "Tgfb1", "Ccl1", "Cd86", "Il12a", "Il12b", "Il1b", "Il6", "Sphk1", "Tnf", "Tnfsf14", "Cd163", "Cxcl13", "Mertk", "Nos2", "Vegfa")) +RotatedAxis() + xlab("M2 Macrophage Markers") + ylab("Liver Cluster")
M2dotplot

Macmonodotplot <- DotPlot(Liver, assay = "SCT", features = c("Apoe","Ccr2","Cd209a","Cd36","Cd74","Cebpa","Cebpd","Csf1r","Cx3cr1","Fos","H2-Aa","H2-Ab1","Itgam","Itgax","Lgals3","Ly6c1","Ly6c2", "Lyz1", "Mef2a","Mmp8","Nr4a1","Pparg", "Akr1b10", "Blvrb", "Camk1", "Glul", "Myo7a", "Nln", "Pcyox1", "Pla2g15", "Pon3", "Slc48a1")) +  RotatedAxis() + xlab("Monocyte / Macrophage Markers") + ylab("Liver Cluster")
Macmonodotplot

DCCytTBCelldotplot <- DotPlot(Liver, assay = "SCT", features = c("Ciita", "Clec9a", "Cst3", "Fcer1a", "Itgae", "Itgam", "Itgax", "Ptprc", "Xcr1", "Il7r", "S100a4", "Cd3d", "Cd3e", "Cd3g", "Cd69", "Cd8a", "Cd8b1", "Itga2", "Cd19", "Cd79a", "Cd79b", "Cr2", "Ebf1", "Ighm", "Igkc", "Iglc1", "Iglc2", "Iglc3", "Ly6d", "Ms4a1","Gnly", "Nkg7",  "Ppbp")) +  RotatedAxis() + xlab("Dendritic Cell / Cytotoxic T Cell / B Cell Markers") + ylab("Liver Cluster")
DCCytTBCelldotplot


pdf("20211228 Liver dot plots 2.pdf", width = 9, height = 5)
M1dotplot
M2dotplot
Macmonodotplot
DCCytTBCelldotplot
dev.off()
