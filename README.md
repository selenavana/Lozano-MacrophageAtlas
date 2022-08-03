# Lozano-MacrophageAtlas

The respective raw_feature_bc_matrix files (barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz) for each organ were processed with the appropriately named .R file. The code is the same for each file except for the organ name. 

(organ) filtered sctransform.R --> remove cells with < 200 features, > 5000 features, or >15 mitochondrial reads. Then, data were normalized with SCTransform. PCA was run. FindNeighbors(), RunUMAP(), and RunTSNE() were run with dims 1:40. Louvain clustering was run with resolution 0.43. This scTransformed data was saved as a .rds file under format "yyyymmdd R code SCtransform filtered umap tsne Organ.rds". This is a critical stopping point. These files were used for downstream analyses after confirming that the .rds files gave the same results as running the processing workflow again. Dendrograms, plots, dot plots, and significant markers for each cluster were also generated from this filtered, scTransformed dataset (after the .rds file was saved). 

20220224 code for integration v4.R --> integrate the SCtransformed organs into "macrophage.integrated" using Seurat's SCTIntegration vignette: https://satijalab.org/seurat/archive/v3.1/integration.html 

Prior to integration, each individual cluster was renamed according to our identifications, and idents were saved under "cell.ident" prior to finding new clusters in the integrated data. 

***NEED custom_seurat_functions files !!!





