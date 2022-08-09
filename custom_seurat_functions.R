plot_integrated_clusters = function (srat) { 
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  
  
  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$cell.ident)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)

  
  colnames(melt_mtx)[2] <- "dataset"

  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value, fill=dataset)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + scale_fill_manual(values = c(brewer.pal(name = "Set3", n = 11), "#FBB4AE", "#B3CDE3", "#E5C494")) +
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
  
  p2
  
}

#  mycolor <- c(brewer.pal(name = "Set3", n = 10), brewer.pal(name = "Pastel1", n = 4))

# ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual("mycolor")
#scale_color_manual(values = "mycolor", labels = c("Other", "Alveolar", "Anti-inflammatory", "Border-Associated", "Dendritic Cells", "Dying", "Heart Tissue Resident", "Interstitial", "Kupffer Cells", "Lymphocytes", "Microglia", "Monocyte-Derived", "Peritoneal", "Pro-inflammatory")) +

