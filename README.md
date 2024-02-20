# Lozano-MacrophageAtlas


*** Seurat analysis and integration
'(organ) filtered sctransform.R' --> The respective raw_feature_bc_matrix files (barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz) for each organ were processed with the appropriately named .R file. The code is the same for each file except for the organ name. Remove cells with < 200 features, > 5000 features, or >15 mitochondrial reads. Then, data were normalized with SCTransform. PCA was run. FindNeighbors(), RunUMAP(), and RunTSNE() were run with dims 1:40. Louvain clustering was run with resolution 0.43. This scTransformed data was saved as a .rds file under format "yyyymmdd R code SCtransform filtered umap tsne Organ.rds". This is a critical stopping point. These files were used for downstream analyses after confirming that the .rds files gave the same results as running the processing workflow again. Dendrograms, plots, dot plots, and significant markers for each cluster were also generated from this filtered, scTransformed dataset (after the .rds file was saved). 

'20230411 code for integration v5.R' --> integrate the SCtransformed organs into "macrophage.integrated" using Seurat's SCTIntegration vignette: https://satijalab.org/seurat/archive/v3.1/integration.html . Prior to integration, each individual cluster was renamed according to our identifications, and idents were saved under "cell.ident" prior to finding new clusters in the integrated data. "Other" cells were removed from the data, leaving only macrophages and dendritic cells.


***Figure generation
'20230412 aggr macs no other V5 pt 2.R' --> Using only macrophages and dendritic cells, created UMAP and tSNE highlighting all cell types/cluster numbers/organ origin. 

'20230417 macrophage integrated V5- trajectory inference.R' --> using only macrophages and dendritic cells, ran monocle3 pseudotime analysis with the root node set as the most distal point among the pro-inflammatory macrophages. https://cole-trapnell-lab.github.io/monocle3/

'20230417 renaming all organ clusters for final pub.R' --> generated publication quality UMAP and tSNE with all clusters labeled

'20230501 indiv organ cell counts.R' --> output tables of organ cell counts by cluster and identification. 

'20230508 organ clusters in aggregate with details.R' --> generated UMAP with each organ's cell contributions highlights

#need final figure R data 



***Extraneous files: 
'20211228 dot plot code.R' --> dot plots used to help in preliminary identification of macrophages
'20220529 brain subclustering.R' --> created violin plot of BAM markers; attempted to identify subclusters
'20230416 code for integration - all cells - v5.R' --> again using Seurat's SCTIntegration vignette, but left in the "other" cells

