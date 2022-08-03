library(Seurat)
library(ggplot2)
macrophage.integrated <- readRDS(file = "~/Selena/20220226 macrophage integrated V4 - other cells removed.rds")

#changed resolution to 0.43, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10
macrophage.integrated <- FindClusters(macrophage.integrated, resolution = 7.5)

#run for each subsequent resolution - 0.43, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10
count.table <- table(macrophage.integrated@meta.data$seurat_clusters, macrophage.integrated@meta.data$orig.ident)
count.table

umap1 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 1.0")
tsne1 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 1.0")
umap0.43 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 0.43")
tsne0.43 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 0.43")
umap1.5 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 1.5")
tsne1.5 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 1.5")
umap2 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 2.0")
tsne2 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 2.0")
umap2.5 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 2.5")
tsne2.5 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 2.5")
umap3 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 3.0")
tsne3 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 3.0")
umap4 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 4.0")
tsne4 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 4.0")
umap5 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 5.0")
tsne5 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 5.0")
umap7.5 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 7.5") + NoLegend()
tsne7.5 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 7.5") + NoLegend()
umap10 <- DimPlot(macrophage.integrated, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Resolution = 10") + NoLegend()
tsne10 <- DimPlot(macrophage.integrated, reduction = "tsne", label = TRUE, repel = TRUE) + ggtitle("Resolution = 10") + NoLegend()


pdf("20220511 integrated macs V4 - changing resolution.pdf", width = 8.5, height = 6)

umap0.43
tsne0.43
umap1
tsne1
umap1.5
tsne1.5
umap2
tsne2
umap2.5
tsne2.5
umap3
tsne3
umap4
tsne4
umap5
tsne5
umap7.5
tsne7.5
umap10
tsne10
dev.off()