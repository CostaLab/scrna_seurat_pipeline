#!/bin/bash
RED='\033[0;31m'
NC='\033[0m' # No Color

PROJ=$1
SAVED_DATA="../save${PROJ}"

FUNCS=(
        --QC
        --Singleton
        --DEs
        # QC
        # DEs
        # Clusters
        # Clusters_harmony
        # Clusters_seurat
        # DEGO
        # KEGG
        # hallmark
        # Reactome
        # EXT_MARKERS
        # DEGO_stage
        # hallmark_stage
        # reactome_stage
        # kegg_stage
        # DEGO_1v1
        # hallmark_1v1
        # reactome_1v1
        # kegg_1v1
	--GET_DATA
)



#!!!!!!!!------clusters to choose-------------
# In general, we choose seurat_clusters,
# If you are using removed or merged clusters,
# choose the following:
            # seurat_clusters
            # merged_clusters
            # removed_clusters
            # remove_recluster

#cluster="removed_clusters"
#cluster="remove_recluster"
#cluster="merged_clusters"
#cluster="annotation"

cluster="singleton"

#!!!!!!!!!----------------------------------

echo $PROJ
echo -e "Use cluster slot ${RED} $cluster ${NC}"
mkdir -p report${PROJ}/data
cp -r ../charts${PROJ}/* report${PROJ}/data
python code_generator.py -c $cluster -cf "conf/config${PROJ}.R" -b "../" -p "$PROJ"
grip --export report${PROJ}/index.md

Rscript make_report.R --proj=$PROJ --cluster=$cluster --save_dir=$SAVED_DATA "${FUNCS[@]}"
