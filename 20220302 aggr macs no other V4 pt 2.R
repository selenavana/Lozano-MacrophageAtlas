library(Seurat)
library(ggplot2)
library(patchwork)


macrophage.integrated <- readRDS(file = "~/Selena/20220226 macrophage integrated V4 - other cells removed.rds")

#finding markers of differentially expressed genes in all clusters
all.markers <- FindAllMarkers(macrophage.integrated)
write.csv(all.markers, "20220302 all genes expression - aggr no other.csv")

#building dendrogram
#install.packages('ape')
macrophage.integrated <- BuildClusterTree(macrophage.integrated)
#PlotClusterTree(macrophage.integrated)
data.tree <- Tool(macrophage.integrated, slot = "BuildClusterTree")
ape::plot.phylo(x = data.tree, direction = "downwards", srt = 90, no.margin = TRUE, label.offset = 0.3, adj = 0.5, use.edge.length = FALSE)

#Find markers between organs
macrophage.integrated <- SetIdent(macrophage.integrated, value = macrophage.integrated@meta.data$orig.ident)
all.markers.organ <- FindAllMarkers(macrophage.integrated)
write.csv(all.markers.organ, "20220303 all genes expression ORGAN - macs only.csv")

#find markers between idents
macrophage.integrated <- SetIdent(macrophage.integrated, value = macrophage.integrated@meta.data$cell.ident)
all.markers.cellident <- FindAllMarkers(macrophage.integrated)
write.csv(all.markers.cellident, "20220303 all genes expression CELL IDENT - macs only.csv")


#highlight each cluster, cell type
{pdf("20220303 aggregate cell identity - macs only - umap.pdf", width = 8.5, height = 6)
  DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Cell Identity")  
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("red", "gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Alveolar")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray","red",  "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Anti-inflammatory")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "red", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Border-Associated")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c( "gray", "gray", "gray", "red","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Cycling")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "red", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Dendritic Cells")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Heart Tissue Resident")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","red", "gray", "gray", "gray","gray", "gray")) + ggtitle("Interstitial")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray","red",  "gray", "gray","gray", "gray")) + ggtitle("Kupffer Cells")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "red", "gray","gray", "gray")) + ggtitle("Microglia")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","red", "gray", "gray")) + ggtitle("Monocyte Derived")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "red", "gray")) + ggtitle("Peritoneal")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "umap", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "red")) + ggtitle("Pro-inflammatory")+ NoLegend()
dev.off()

pdf("20220303 aggregate cell identity - macs only - tsne.pdf", width = 8.5, height = 6)
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Cell Identity")  
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("red", "gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Alveolar")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray","red",  "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Anti-inflammatory")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "red", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Border-Associated")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c( "gray", "gray", "gray", "red","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Cycling")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "red", "gray", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Dendritic Cells")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray")) + ggtitle("Heart Tissue Resident")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","red", "gray", "gray", "gray","gray", "gray")) + ggtitle("Interstitial")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray","red",  "gray", "gray","gray", "gray")) + ggtitle("Kupffer Cells")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "red", "gray","gray", "gray")) + ggtitle("Microglia")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","red", "gray", "gray")) + ggtitle("Monocyte Derived")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "red", "gray")) + ggtitle("Peritoneal")+ NoLegend()
DimPlot(macrophage.integrated, reduction = "tsne", group.by = "cell.ident", cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "red")) + ggtitle("Pro-inflammatory")+ NoLegend()
dev.off()
}

#highlight each cluster, number
{pdf("20220303 aggregate cluster - macs only - umap.pdf", width = 8.5, height = 6)
  DimPlot(macrophage.integrated, reduction = "umap",  label = TRUE, repel = TRUE)
  DimPlot(macrophage.integrated, reduction = "umap", cols = c("red", "gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray","red",  "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "red", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c( "gray", "gray", "gray", "red","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "red", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","red", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray","red",  "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","red", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "red", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray", "red","gray", "gray", "gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray",  "red","gray", "gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray",  "gray","gray", "gray", "red","gray","gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","red", "gray", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "red", "gray")) 
  DimPlot(macrophage.integrated, reduction = "umap",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray",  "gray","gray", "gray", "gray","gray", "gray", "red")) 
dev.off()

pdf("20220303 aggregate cluster - macs only - tsne.pdf", width = 8.5, height = 6)
DimPlot(macrophage.integrated, reduction = "tsne",  label = TRUE, repel = TRUE)
DimPlot(macrophage.integrated, reduction = "tsne", cols = c("red", "gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray","red",  "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "red", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c( "gray", "gray", "gray", "red","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "red", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","red", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray","red",  "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","red", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray"))
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "red", "gray", "gray","gray", "gray", "gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "red", "gray","gray", "gray", "gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray", "red","gray", "gray", "gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray",  "red","gray", "gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray",  "gray","gray", "gray", "red","gray","gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","red", "gray", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "gray", "gray","gray", "red", "gray")) 
DimPlot(macrophage.integrated, reduction = "tsne",  cols = c("gray", "gray", "gray", "gray", "gray", "gray","gray", "gray", "gray","gray", "gray",  "gray","gray", "gray", "gray","gray", "gray", "red")) 
dev.off()
}

#violin plots for each cluster/ident
cluster.features <- c("Btnl9", "Cyp1b1", "Hfe", "H3f3b", "Cda", "Rgs2", "Cobll1", "H2-DMb2", "C5ar1", "Gstm1", "Prss23", "F11r", "Cxcl16", "Tuba1b", "Csf3r", "Lama4", "Pltp", "Mbd1")
cluster.vlnplot <- VlnPlot(macrophage.integrated, cluster.features, stack = TRUE, flip = TRUE, assay = "SCT")
cluster.vlnplot
  

macrophage.integrated <- SetIdent(macrophage.integrated, value = macrophage.integrated@meta.data$cell.ident)
ident.features <- c("Blvrb", "Cxcl3", "Sirpb1a", "Lrrc8c", "Cd3eap", "Tsc22d1", "Cadm1", "Rerg", "Cda", "Plac8", "Retnla", "Csf3r")
ident.vlnplot <- VlnPlot(macrophage.integrated, ident.features, stack = TRUE, flip = TRUE, assay = "SCT")
ident.vlnplot



