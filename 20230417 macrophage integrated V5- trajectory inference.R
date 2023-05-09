if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('satijalab/seurat-data')
install.packages("R.utils")
remotes::install_github('satijalab/seurat-wrappers')

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

macrophage.integrated <- readRDS(file = "20230412 macrophage integrated V5 - other cells removed.rds")

DimPlot(macrophage.integrated, group.by = c("orig.ident", "ident"))

cds <- as.cell_data_set(macrophage.integrated)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = NULL, root_cells = NULL, verbose = FALSE)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p3 <- plot_cells(cds, color_cells_by = "seurat_clusters", graph_label_size = 4, group_label_size = 4, label_leaves = FALSE, label_branch_points = FALSE, label_groups_by_cluster = TRUE, cell_size = 0, cell_stroke = 0.75)
p3
p4 <- plot_cells(cds, color_cells_by = "cell.ident", graph_label_size = 4, group_label_size = 4, label_leaves = FALSE, label_branch_points = FALSE, label_groups_by_cluster = FALSE, trajectory_graph_color = "gray50", cell_size = 0, cell_stroke = 0.75)
p4
wrap_plots(p1, p2)




pdf("20230426 trajectory inf - aggregate - macs only V5.pdf", width = 8.5, height = 6)
p3
p4
p5
dev.off()

ggsave("20230426 trajectory interference - cluster1.png", p3, width = 8.5, height = 6, dpi = 100)

p5 <- plot_cells(cds, color_cells_by = "pseudotime", graph_label_size = 4, group_label_size = 4, show_trajectory_graph = FALSE, label_leaves = FALSE, label_branch_points = FALSE, trajectory_graph_color = "azure3", label_groups_by_cluster = TRUE, cell_size = 0, cell_stroke = 0.75) + theme(legend.position = c(0.95, 0.15)) 
p5

#more random stuff??
integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = TRUE, label_branch_points = TRUE)
FeaturePlot(macrophage.integrated, "monocle3_pseudotime")