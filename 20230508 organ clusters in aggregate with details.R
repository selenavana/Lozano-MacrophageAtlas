library(Seurat)
library(ggplot2)
library(patchwork)

#Load data
Heart <- readRDS(file = "20210518 R Code SCtransform filtered umap tsne Heart.rds")
Liver <- readRDS(file = "20210517 R Code SCtransform filtered umap tsne Liver.rds")
Brain <- readRDS(file = "20210511 R Code SCtransform filtered umap tsne Brain.rds")
Lung <- readRDS(file = "20210517 R Code SCtransform filtered umap tsne Lung.rds")

#save old cluster idents
Heart[["original.cluster"]] <- Idents(Heart)
Brain[["original.cluster"]] <- Idents(Brain)
Liver[["original.cluster"]] <- Idents(Liver)
Lung[["original.cluster"]] <- Idents(Lung)


#renaming clusters
{Heart <- RenameIdents(Heart, "0" = "0 - Resident Macrophage/Monocyte", "1" = "1 - Proliferating Cells", "2" = "2 - Dendritic Cells", "3" = "3 - Maturing Dendritic Cells", "4" = "4 - Vascular Endothelial Cells", "5" = "5 - Niche Response Mac", 
                       "6" = "6 - Infiltrating Monocytes", "7" = "7 - Monocyte Derived", "8" = "8 - Resident Cardiac Mac",
                       "9" = "9 - Fibroblasts/Myofibroblasts", "10" = "10 - Vascular Smooth Muscle Cells", "11" = "11 - Cardiac Stem Cells", 
                       "12" = "12 - MHCII-hi Monocyte Derived",  "13" = "13 - M2 Monocyte Derived Mac", "14" = "14 - Cardiomyocytes")
}

{Brain <- RenameIdents(Brain, "0" = "0 - Microglia (1)", "1" = "1 - Microglia (2)", "2" = "2 - Microglia (3)", "3" = "3 - Endothelial Cells", "4" = "4 - Microglia (dying)", "5" = "5 - Microglia (4)", 
                       "6" = "6 - Microglia (Cycling)", "7" = "7 - MHCII-hi BAM", "8" = "8 - Microglia (5)", "9" = "9 - MHCII-lo BAM",  "10" = "10 - Microglia (6)", 
                       "11" = "11 - Microglia (7)", "12" = "12 - Smooth Muscle Cells", "13" = "13 - MHCII-lo Monocyte Derived BAM", "14" = "14 - Mixed Types", "15" = "15 - Mitotic Cells")
}

{Liver <- RenameIdents(Liver, "0" = "0 - Kupffer Cells (1)", "1" = "1 - Kupffer Cells (2)", "2" = "2 - Monocyte Derived Ly6c-hi", "3" = "3 - Monocyte Derived Ly6c-lo", "4" = "4 - Sinusoidal Endothelial Cells", 
                       "5" = "5 - Pro-inflammatory Mac", "6" = "6 - Restorative Mac", "7" = "7 - B Cells", "8" = "8 - Neutrophils", "9" = "9 - Dendritic Cells",  "10" = "10 - Cytotoxic T Cells", "11" = "11 - Peritoneal Mac", "12" = "12 - Hepatocytes")
  
}  

{Lung <- RenameIdents(Lung, "0" = "0 - Alveolar Mac (1)", "1" = "1 - Monocyte Derived (1)", "2" = "2 - Lack of Expression (1)", "3" = "3 - Pro-inflammatory M1 (1)", "4" = "4 - Interstitial Mac 1", "5" = "5 - Interstitial Mac 2", 
                      "6" = "6 - Endothelial Cells", "7" = "7 - Monocyte Derived (2)", "8" = "8 - Pro-inflammatory M1 (2)", "9" = "9 - Monocyte Derived (3)",  "10" = "10 - Alveolar Mac (2)", "11" = "11 - Lack of Expression (2)", "12" = "12 - B Cells", "13" = "13 - Alveolar Epithelium", 
                      "14" = "14 - Alveolar Mac (Proliferating)", "15" = "15 - Interstitial Mac 3 (1)", "16" = "16 - Other", "17" = "17 - CD8+ T cell", "18" = "18 - Gland Epithelium", "19" = "19 - Interstitial Mac 3 (2)")
}   

#generating plots

he2 <- DimPlot(Heart, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Heart") 
br2 <- DimPlot(Brain, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Brain")
li2 <- DimPlot(Liver, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Liver")
lu2 <- DimPlot(Lung, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("Lung")

#aggregated analysis - no macrophages
macrophage.integrated <- readRDS(file = "20230412 macrophage integrated V5 - other cells removed.rds")

br1 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", cols = c("red", "gray88", "gray88", "gray88")) + ggtitle("Brain") + theme(legend.position = "none")
he1 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", cols = c("gray88", "red", "gray88", "gray88")) + ggtitle("Heart") + theme(legend.position = "none")
li1 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", cols = c( "gray88", "gray88", "red", "gray88")) + ggtitle("Liver") + theme(legend.position = "none")
lu1 <- DimPlot(macrophage.integrated, reduction = "umap", group.by = "orig.ident", cols = c("gray88", "gray88", "gray88", "red")) + ggtitle("Lung") + theme(legend.position = "none")

pdf("20230508 aggregate organ cells with detail.pdf", width = 17, height = 6)
br1+br2
he1+he2
li1+li2
lu1+lu2
dev.off()