library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

macrophage.integrated <- readRDS(file = "~/Selena/20220226 macrophage integrated V4 - other cells removed.rds")

DimPlot(macrophage.integrated, group.by = c("orig.ident", "ident"))

cds <- as.cell_data_set(macrophage.integrated)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = NULL, root_cells = NULL, verbose = FALSE)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p3 <- plot_cells(cds, color_cells_by = "seurat_clusters", graph_label_size = 2, group_label_size = 4, label_leaves = FALSE, label_branch_points = FALSE)
p3
p4 <- plot_cells(cds, color_cells_by = "cell.ident", graph_label_size = 2, group_label_size = 4, label_leaves = FALSE, label_branch_points = FALSE)
p4
wrap_plots(p1, p2)

pdf("20220621 trajectory inf - aggregate - macs only - umap.pdf", width = 8.5, height = 6)
p3
p4
dev.off()

p5 <- plot_cells(cds, color_cells_by = "pseudotime", graph_label_size = 2, group_label_size = 4, show_trajectory_graph = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
p5

#more random stuff??
integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)
FeaturePlot(macrophage.integrated, "monocle3_pseudotime")