library(Seurat)
library(patchwork)
library(ggplot2)

M1dotplot <- DotPlot(Heart, assay = "SCT", features = c("Cd38","Ccr7","Cd80","Cd86","Cxcl9","Cybb","Fcgr1","Fpr2","Gpr18","H-2", "Il12a", "Il12b","Il1b", "Il23a","Il4ra","Il6","Marco", "Nos2", "Slamf1","Socs3","Stat1","Tnf")) +  RotatedAxis() + xlab("M1 Macrophage Markers") + ylab("Heart Cluster")
M1dotplot

M2dotplot <- DotPlot(Heart, assay = "SCT", features = c("Arg1", "Chil3", "Egr2", "Il27ra", "Mrc1", "Myc", "Pparg", "Retnla", "Stat6", "CCL17", "Clec7a", "Fn1", "Igf1", "Il10", "Il1r1", "Tgfb1", "Ccl1", "Cd86", "Il12a", "Il12b", "Il1b", "Il6", "Sphk1", "Tnf", "Tnfsf14", "Cd163", "Cxcl13", "Mertk", "Nos2", "Vegfa")) +RotatedAxis() + xlab("M2 Macrophage Markers") + ylab("Heart Cluster")
M2dotplot

Macmonodotplot <- DotPlot(Heart, assay = "SCT", features = c("Apoe","Ccr2","Cd209a","Cd36","Cd74","Cebpa","Cebpd","Csf1r","Cx3cr1","Fos","H2-Aa","H2-Ab1","Itgam","Itgax","Lgals3","Ly6c1","Ly6c2", "Lyz1", "Mef2a","Mmp8","Nr4a1","Pparg", "Akr1b10", "Blvrb", "Camk1", "Glul", "Myo7a", "Nln", "Pcyox1", "Pla2g15", "Pon3", "Slc48a1")) +  RotatedAxis() + xlab("Monocyte / Macrophage Markers") + ylab("Heart Cluster")
Macmonodotplot

DCCytTBCelldotplot <- DotPlot(Heart, assay = "SCT", features = c("Ciita", "Clec9a", "Cst3", "Fcer1a", "Itgae", "Itgam", "Itgax", "Ptprc", "Xcr1", "Il7r", "S100a4", "Cd3d", "Cd3e", "Cd3g", "Cd69", "Cd8a", "Cd8b1", "Itga2", "Cd19", "Cd79a", "Cd79b", "Cr2", "Ebf1", "Ighm", "Igkc", "Iglc1", "Iglc2", "Iglc3", "Ly6d", "Ms4a1","Gnly", "Nkg7",  "Ppbp")) +  RotatedAxis() + xlab("Dendritic Cell / Cytotoxic T Cell / B Cell Markers") + ylab("Heart Cluster")
DCCytTBCelldotplot


pdf("20211228 Heart dot plots 2.pdf", width = 8, height = 5)
M1dotplot
M2dotplot
Macmonodotplot
DCCytTBCelldotplot
dev.off()