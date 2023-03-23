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

## version 1.0.4  --- 12.24.2020
 * Fix fishertest failed due to the calculating order
 * Fix more cpu cores used than being set
 * Change organizing of the index file
 * Change Regressout using conf rather than different functions
 * Change Most of viz are based on scrna_phase_comparing.Rds
 * Change execution plan to several groups, make it sample
    * scrna_phase_preprocess -- QC and regressOut mito etc.
    * scrna_phase_clustering -- integration & clustering
    * scrna_phase_comparing -- all other analysis
    * scrna_phase_singleton -- single sample
 * Add static/phase.ini to control each group of analysis functions
    * each section is a group defines which and order to run in this group
    * if set to 0 will omit
 * Add harmony integration
 * Add 2 viz clustering Rmds to exclusively plot harmony and seurat clusters
 * Add G1.Score to clustering viz
 * Add catching exception when comparing has been calculated
 * Add once encountering error, stop for viz
 * Add unimplemented function scrna_HCLannotate

## version 1.0.5  --- 03.23.2023
 * add allinone optioin to decide store calculation result into seurat object or files.
 * explicitly display parameters in run_example.sh
 * compatible with Seurat 4
 * viz add three ways of display features in umap(seurat,nebulosa and schex)
 * viz doublets score
 * viz QCC for a exisiting Seurat object
 * viz improve the acolor scheme
 * viz add references for some plot sections
 * viz add geneset analysis
 * viz tables to show counts and statistics of clusters and samples
 * add PanglaoDB to external markers
 * add ambient detection
 * add doublet detection
 * add scProportion to compare conditions
 * add Rmagic to impute reads  *
 * improve the performance of saving and loading using compression options
 * fix some bugs
