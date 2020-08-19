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
