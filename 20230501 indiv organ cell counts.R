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
{Heart <- RenameIdents(Heart, "0" = "Monocyte Derived", "4" = "Other", "6" = "Monocyte Derived", "7" = "Monocyte Derived", 
                      "1" = "Cycling", "9" = "Other", "10" = "Other", "11" = "Other", "14" = "Other", 
                      "2" = "Dendritic Cells", "3" = "Dendritic Cells", 
                      "5" = "Cycling", "12" = "Pro-inflammatory", "8" = "Heart Tissue Resident", "13" = "Anti-inflammatory")
}

{Brain <- RenameIdents(Brain, "0" = "Microglia", "1" = "Microglia", "2" = "Microglia", "3" = "Other", "4" = "Cycling", "5" = "Microglia", 
                       "6" = "Cycling", "7" = "Border-Associated", "8" = "Microglia", "9" = "Border-Associated",  "10" = "Microglia", 
                       "11" = "Microglia", "12" = "Other", "13" = "Border-Associated", "14" = "Other", "15" = "Other")
}

{Liver <- RenameIdents(Liver, "0" = "Kupffer Cells", "1" = "Kupffer Cells", "2" = "Monocyte Derived", "3" = "Monocyte Derived", "4" = "Other", 
                       "5" = "Pro-inflammatory", "6" = "Anti-inflammatory", "7" = "Lymphocytes", "8" = "Other", "9" = "Dendritic Cells",  "10" = "Lymphocytes", "11" = "Peritoneal", "12" = "Other")

}  

{Lung <- RenameIdents(Lung, "0" = "Alveolar", "1" = "Monocyte Derived", "2" = "Cycling", "3" = "Pro-inflammatory", "4" = "Interstitial", "5" = "Interstitial", 
                        "6" = "Other", "7" = "Monocyte Derived", "8" = "Pro-inflammatory", "9" = "Monocyte Derived",  "10" = "Alveolar", "11" = "Cycling", "12" = "Lymphocytes", "13" = "Other", 
                        "14" = "Alveolar", "15" = "Interstitial", "16" = "Other", "17" = "Lymphocytes", "18" = "Other", "19" = "Interstitial")
}   


#Heart table counts
Braincount.table <- table(Brain@meta.data$seurat_clusters, Brain@active.ident)
Braincount.table

Heartcount.table <- table(Heart@meta.data$seurat_clusters, Heart@active.ident)
Heartcount.table

Livercount.table <- table(Liver@meta.data$seurat_clusters, Liver@active.ident)
Livercount.table

Lungcount.table <- table(Lung@meta.data$seurat_clusters, Lung@active.ident)
Lungcount.table
