library(Seurat)
library(ggplot2)
library(patchwork)


macrophage.integrated <- readRDS(file = "20230412 macrophage integrated V5 - other cells removed.rds")

#finding markers of differentially expressed genes in all clusters
all.markers <- FindAllMarkers(macrophage.integrated)
write.csv(all.markers, "20230412 all genes expression - aggr no other.csv")

#building dendrogram
#install.packages('ape')
macrophage.integrated <- BuildClusterTree(macrophage.integrated)
#PlotClusterTree(macrophage.integrated)
data.tree <- Tool(macrophage.integrated, slot = "BuildClusterTree")
ape::plot.phylo(x = data.tree, direction = "downwards", srt = 90, no.margin = TRUE, label.offset = 0.3, adj = 0.5, use.edge.length = FALSE)

#Find markers between organs
macrophage.integrated <- SetIdent(macrophage.integrated, value = macrophage.integrated@meta.data$orig.ident)
all.markers.organ <- FindAllMarkers(macrophage.integrated)
write.csv(all.markers.organ, "20230412 all genes expression ORGAN - macs only.csv")

#find markers between idents
macrophage.integrated <- SetIdent(macrophage.integrated, value = macrophage.integrated@meta.data$cell.ident)
all.markers.cellident <- FindAllMarkers(macrophage.integrated)
write.csv(all.markers.cellident, "20230412 all genes expression CELL IDENT - macs only.csv")


#highlight each cluster, cell type
{pdf("20230412 aggregate cell identity - macs only - umap.pdf", width = 17, height = 18)
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Cell Identity")  
up1 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("red", "gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Alveolar")+ NoLegend()
up2 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88","red",  "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Anti-inflammatory")+ NoLegend()
up3 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "red", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Border-Associated")+ NoLegend()
up4 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c( "gray88", "gray88", "gray88", "red","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Cycling")+ NoLegend()
up5 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "red", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Dendritic Cells")+ NoLegend()
up6 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Heart Tissue Resident")+ NoLegend()
up7 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","red", "gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Interstitial")+ NoLegend()
up8 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88","red",  "gray88", "gray88","gray88", "gray88")) + ggtitle("Kupffer Cells")+ NoLegend()
up9 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "red", "gray88","gray88", "gray88")) + ggtitle("Microglia")+ NoLegend()
up10 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","red", "gray88", "gray88")) + ggtitle("Monocyte Derived")+ NoLegend()
up11 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "red", "gray88")) + ggtitle("Peritoneal")+ NoLegend()
up12 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "red")) + ggtitle("Pro-inflammatory")+ NoLegend()
ucomp1 <- up1 + up2 + up3 + up4 + up5 + up6 + plot_layout(ncol = 2)
ucomp2 <- up7 + up8 + up9 + up10 + up11 + up12 + plot_layout(ncol = 2)
ucomp1
ucomp2
dev.off()

pdf("20230412 aggregate cell identity - macs only - tsne.pdf", width = 17, height = 18)
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Cell Identity")  
tp1 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("red", "gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Alveolar")+ NoLegend()
tp2 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88","red",  "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Anti-inflammatory")+ NoLegend()
tp3 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "red", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Border-Associated")+ NoLegend()
tp4 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c( "gray88", "gray88", "gray88", "red","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Cycling")+ NoLegend()
tp5 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "red", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Dendritic Cells")+ NoLegend()
tp6 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Heart Tissue Resident")+ NoLegend()
tp7 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","red", "gray88", "gray88", "gray88","gray88", "gray88")) + ggtitle("Interstitial")+ NoLegend()
tp8 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88","red",  "gray88", "gray88","gray88", "gray88")) + ggtitle("Kupffer Cells")+ NoLegend()
tp9 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "red", "gray88","gray88", "gray88")) + ggtitle("Microglia")+ NoLegend()
tp10 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","red", "gray88", "gray88")) + ggtitle("Monocyte Derived")+ NoLegend()
tp11 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "red", "gray88")) + ggtitle("Peritoneal")+ NoLegend()
tp12 <- DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "red")) + ggtitle("Pro-inflammatory")+ NoLegend()
tcomp1 <- tp1 + tp2 + tp3 + tp4 + tp5 + tp6 + plot_layout(ncol = 2)
tcomp1
tcomp2 <- tp7 + tp8 + tp9 + tp10 + tp11 + tp12 + plot_layout(ncol = 2)
tcomp2
dev.off()
}

#highlight each cluster, number - didn't do this 3/16/23
{pdf("20230412 aggregate cluster - macs only - umap.pdf", width = 8.5, height = 6)
  DimPlot(macrophage.integrated, reduction = "umap",  label = TRUE, repel = TRUE)
  DimPlot(macrophage.integrated, reduction = "umap", cols = c("red", "gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88","red",  "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "red", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c( "gray88", "gray88", "gray88", "red","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "red", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","red", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88","red",  "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","red", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "red", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88", "red","gray88", "gray88", "gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88",  "red","gray88", "gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88",  "gray88","gray88", "gray88", "red","gray88","gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","red", "gray88", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "red", "gray88")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88",  "gray88","gray88", "gray88", "gray88","gray88", "gray88", "red")) 
dev.off()

pdf("20220303 aggregate cluster - macs only - tsne.pdf", width = 8.5, height = 6)
DimPlot(macrophage.integrated, reduction = "tsne",  label = TRUE, repel = TRUE)
DimPlot(macrophage.integrated, reduction = "tsne", cols = c("red", "gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88","red",  "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "red", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c( "gray88", "gray88", "gray88", "red","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "red", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","red", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88","red",  "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","red", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "red", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "red", "gray88","gray88", "gray88", "gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88", "red","gray88", "gray88", "gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88",  "red","gray88", "gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88",  "gray88","gray88", "gray88", "red","gray88","gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","red", "gray88", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "red", "gray88")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray88", "gray88", "gray88", "gray88", "gray88", "gray88","gray88", "gray88", "gray88","gray88", "gray88",  "gray88","gray88", "gray88", "gray88","gray88", "gray88", "red")) 
dev.off()
}

#just monocyte-derived umap
pdf("20230412 aggregate cell identity - macs only - umap - figure 3.pdf", width = 8.5, height = 6)
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", label = TRUE, repel = TRUE) + ggtitle(" ") + theme(legend.position = c(0.78, 0.9), legend.text = element_text(size = 10), legend.key.size = unit(1,'mm'))
up10
dev.off()

#violin plots for each cluster/ident  - didnt do yet
cluster.features <- c("Btnl9", "Cyp1b1", "Hfe", "H3f3b", "Cda", "Rgs2", "Cobll1", "H2-DMb2", "C5ar1", "Gstm1", "Prss23", "F11r", "Cxcl16", "Tuba1b", "Csf3r", "Lama4", "Pltp", "Mbd1")
cluster.vlnplot <- VlnPlot(macrophage.integrated, cluster.features, stack = TRUE, flip = TRUE, assay = "SCT")
cluster.vlnplot
  

macrophage.integrated <- SetIdent(macrophage.integrated, value = macrophage.integrated@meta.data$cell.ident)
ident.features <- c("Blvrb", "Cxcl3", "Sirpb1a", "Lrrc8c", "Cd3eap", "Tsc22d1", "Cadm1", "Rerg", "Cda", "Plac8", "Retnla", "Csf3r")
ident.vlnplot <- VlnPlot(macrophage.integrated, ident.features, stack = TRUE, flip = TRUE, assay = "SCT")
ident.vlnplot

#figure 4 updates

pdf("20230412 aggregate organ identity - macs only - umap - figure 4.pdf", width = 8.5, height = 6)
DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE) + ggtitle(" ")  + theme(legend.position = c(0.84, 0.15)) 
DimPlot(macrophage.integrated, reduction = "umap", group.by = "agg.ident", label = TRUE, repel = TRUE) + ggtitle(" ")  
dev.off()

pdf("20230321 aggregate organ identity - macs only - tsne - figure 4 9x6.pdf", width = 8.5, height = 6)
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "orig.ident", label = FALSE, repel = TRUE) + ggtitle(" ")  + theme(legend.position = c(0.84, 0.15)) 
DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle(" ")  
dev.off()

#clusterplot 

pdf("20230417 aggregate organ identity - macs only - fraction of cells.pdf", width = 8.5, height = 6)
source("custom_seurat_functions no other.R")
plot_integrated_clusters(macrophage.integrated)

source("custom_seurat_functions V2.R")
plot_integrated_clusters2(macrophage.integrated)

ape::plot.phylo(x = data.tree, direction = "downwards", srt = 90, no.margin = TRUE, label.offset = 0.3, adj = 0.5, use.edge.length = FALSE)

dev.off()
