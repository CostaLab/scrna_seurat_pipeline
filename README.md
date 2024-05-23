# scRNA pipeline

### Documents:
https://costalab.github.io/scrna_seurat_pipeline

#### R sessionInfo

```R
r$> sessionInfo()
R version 4.2.2 (2022-10-31)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Rocky Linux 8.7 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /data/sz753404/miniconda3/envs/R422/lib/libopenblasp-r0.3.21.so

locale:
[1] C

attached base packages:
 [1] parallel  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] genesorteR_0.4.3            SoupX_1.6.2                 celda_1.14.2                SingleCellExperiment_1.20.0
 [5] SummarizedExperiment_1.28.0 GenomicRanges_1.50.0        GenomeInfoDb_1.34.8         MatrixGenerics_1.10.0      
 [9] matrixStats_0.63.0          doParallel_1.0.17           iterators_1.0.14            foreach_1.5.2              
[13] data.table_1.14.8           Matrix_1.5-3                future.apply_1.10.0         future_1.32.0              
[17] progeny_1.20.0              futile.logger_1.4.3         DOSE_3.24.2                 plotly_4.10.1              
[21] digest_0.6.31               WriteXLS_6.4.0              gridExtra_2.3               ggsci_3.0.0                
[25] ggridges_0.5.4              SeuratObject_4.1.3          Seurat_4.3.0                EnhancedVolcano_1.16.0     
[29] ggrepel_0.9.3               ComplexHeatmap_2.14.0       openxlsx_4.2.5.2            stringr_1.5.0              
[33] clustree_0.5.0              ggraph_2.1.0                urltools_1.7.3              reshape2_1.4.4             
[37] cowplot_1.1.1               dplyr_1.1.0                 kableExtra_1.3.4            knitr_1.42                 
[41] Hmisc_5.0-1                 glue_1.6.2                  crayon_1.5.2                optparse_1.7.3             
[45] msigdbr_7.5.1               ReactomePA_1.42.0           org.Hs.eg.db_3.16.0         org.Mm.eg.db_3.16.0        
[49] clusterProfiler_4.7.1.003   GSEABase_1.60.0             graph_1.76.0                annotate_1.76.0            
[53] XML_3.99-0.13               AnnotationDbi_1.60.0        IRanges_2.32.0              S4Vectors_0.36.0           
[57] Biobase_2.58.0              BiocGenerics_0.44.0         scHCL_0.1.1                 scMCA_0.2.0                
[61] ggplot2_3.4.1              

loaded via a namespace (and not attached):
  [1] ica_1.0-3                  svglite_2.1.1              assertive.properties_0.0-5 lmtest_0.9-40             
  [5] MASS_7.3-58.3              nlme_3.1-162               backports_1.4.1            GOSemSim_2.24.0           
  [9] rlang_1.1.0                XVector_0.38.0             HDO.db_0.99.1              ROCR_1.0-11               
 [13] irlba_2.3.5.1              BiocParallel_1.32.5        rjson_0.2.21               bit64_4.0.5               
 [17] pheatmap_1.0.12            sctransform_0.3.5          spatstat.sparse_3.0-1      spatstat.geom_3.1-0       
 [21] tidyselect_1.2.0           fitdistrplus_1.1-8         tidyr_1.3.0                assertive.types_0.0-3     
 [25] zoo_1.8-11                 xtable_1.8-4               magrittr_2.0.3             evaluate_0.20             
 [29] cli_3.6.0                  zlibbioc_1.44.0            rstudioapi_0.14            miniUI_0.1.1.1            
 [33] sp_1.6-0                   rpart_4.1.19               fastmatch_1.1-3            lambda.r_1.2.4            
 [37] treeio_1.22.0              RcppEigen_0.3.3.9.3        shiny_1.7.4                xfun_0.37                 
 [41] clue_0.3-64                gson_0.1.0                 cluster_2.1.4              tidygraph_1.2.3           
 [45] KEGGREST_1.38.0            tibble_3.2.0               ape_5.7-1                  listenv_0.9.0             
 [49] Biostrings_2.66.0          png_0.1-8                  withr_2.5.0                bitops_1.0-7              
 [53] ggforce_0.4.1              plyr_1.8.8                 assertive.base_0.0-9       pillar_1.8.1              
 [57] GlobalOptions_0.1.2        cachem_1.0.7               GetoptLong_1.0.5           graphite_1.44.0           
 [61] vctrs_0.5.2                ellipsis_0.3.2             generics_0.1.3             tools_4.2.2               
 [65] foreign_0.8-84             munsell_0.5.0              tweenr_2.0.2               fgsea_1.24.0              
 [69] DelayedArray_0.24.0        fastmap_1.1.1              compiler_4.2.2             abind_1.4-5               
 [73] httpuv_1.6.9               GenomeInfoDbData_1.2.9     enrichR_3.1                lattice_0.20-45           
 [77] deldir_1.0-6               utf8_1.2.3                 later_1.3.0                jsonlite_1.8.4            
 [81] multipanelfigure_2.1.2     scales_1.2.1               tidytree_0.4.2             pbapply_1.7-0             
 [85] lazyeval_0.2.2             promises_1.2.0.1           goftest_1.2-3              spatstat.utils_3.0-2      
 [89] reticulate_1.25            checkmate_2.1.0            rmarkdown_2.20             webshot_0.5.4             
 [93] Rtsne_0.16                 downloader_0.4             uwot_0.1.14                igraph_1.4.1              
 [97] survival_3.5-5             systemfonts_1.0.4          htmltools_0.5.4            memoise_2.0.1             
[101] graphlayouts_0.8.4         viridisLite_0.4.1          mime_0.12                  rappdirs_0.3.3            
[105] futile.options_1.0.1       RSQLite_2.3.0              yulab.utils_0.0.6          blob_1.2.3                
[109] shinythemes_1.2.0          splines_4.2.2              Formula_1.2-5              RCurl_1.98-1.10           
[113] assertive.numbers_0.0-2    colorspace_2.1-0           base64enc_0.1-3            shape_1.4.6               
[117] assertive.files_0.0-2      aplot_0.1.10               nnet_7.3-18                mclust_6.0.0              
[121] Rcpp_1.0.10                RANN_2.6.1                 circlize_0.4.15            enrichplot_1.18.3         
[125] fansi_1.0.4                parallelly_1.34.0          R6_2.5.1                   lifecycle_1.0.3           
[129] formatR_1.14               zip_2.2.2                  curl_4.3.3                 leiden_0.4.3              
[133] getopt_1.20.3              qvalue_2.30.0              RcppAnnoy_0.0.20           RColorBrewer_1.1-3        
[137] spatstat.explore_3.1-0     htmlwidgets_1.6.1          polyclip_1.10-4            triebeard_0.4.1           
[141] purrr_1.0.1                shadowtext_0.1.2           gridGraphics_0.5-1         reactome.db_1.82.0        
[145] rvest_1.0.3                globals_0.16.2             htmlTable_2.4.1            patchwork_1.1.2           
[149] spatstat.random_3.1-4      progressr_0.13.0           codetools_0.2-19           GO.db_3.16.0              
[153] MCMCprecision_0.4.0        gtable_0.3.1               DBI_1.1.3                  ggfun_0.0.9               
[157] tensor_1.5                 httr_1.4.5                 KernSmooth_2.23-20         stringi_1.7.12            
[161] farver_2.1.1               viridis_0.6.2              magick_2.7.4               ggtree_3.6.2              
[165] DT_0.27                    xml2_1.3.3                 combinat_0.0-8             ggplotify_0.1.0           
[169] scattermore_0.8            bit_4.0.5                  scatterpie_0.1.8           spatstat.data_3.0-1       
[173] pkgconfig_2.0.3            babelgene_22.9            

```
