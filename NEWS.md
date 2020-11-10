## version 1.0.0  --- 07.20.2020
 * Add version info
 * filtering add mito&ribo filter
 * config file content to R object
 * function add return code 
 * fisher test ignore 0 cell clusters
 * retain seurat_clusters when use removed_clusters
 * External annotation, replace highest score with mean(x)\*log(length(x)+1)
 * install packages script add R-glue and py-grip(genereate markdown) 
 * viz cluster, fix bug when viz removed or merged clusters
 * viz cluster, no display of 0 cell clusters
 * viz DE, add heatmap
 * viz FeaturePlot add cluster labels
 * viz run shell add option to easily choose which to run
 * viz add external markers violin dotplot and featureplot
 * viz add project name


## version 1.0.1  --- 08.19.2020
 * GO analysis switch in between human and mouse databases
 * Fix source sequencing data cannot load from absolute PATH
 * Fix bug parameter and execution lost after splitting seurat Object during Integration
 * Fix viz when there are no marker genes in some clusters



## version 1.0.2  --- 10.21.2020
 * Fix bug special characters for sheet names failed to produce excel results
 * Fix bug GO analysis lose some results cause by non-GO terms
 * Fix bug mito and ribo genes cannot recognized when analyzing Human organs
 * Change Scale mito&ribo separately
 * Change fishertest for clusters, 1vsothers to 1vs1
 * Change findMarkers to keep all foldchange results for volcano plots
 * Change GO analysis to keep all p adjust value for heatmap
 * Add pathways KEGG,Reatome,Hallmark like GO analysis
 * Add delete mito genes permanently, Dangerous to call!!!
 * Add single sample without integration analysis function
 * Add genesorteR analysis
 * Add DE&GO stage 1 versus the rest
 * Add cluster-wise filter to keep only cells meet the mito or ribo threshold
 * Add Volcano plots and GO heatmap


## version 1.0.3  --- 11.10.2020
 * Fix install_packages.R scMCA from BiocManager to install_github
 * Change cycyle var to regress from CC.Difference to G2M.Score&S.Score
 * Add 3 function with cellcycle regressout
 * Add G1.Score to Metadata G1.Score=1-S.Score-G2M.Score
 * Add violin plot for cell cycle
 * Add tables for quanlity check and clusters

