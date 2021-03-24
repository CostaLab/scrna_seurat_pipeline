#!/bin/bash



AUTHOR="CostaLab"
DEFAULT_OUTPUT_FOLDER="report_"
MAKE_ELEMENTS=0
MAKE_SINGLE_FILE=0


# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} [-h] [-p PROJECT_NAME] [-a AUTHOR] [-m] [-so] [-s SAVE_PATH] [-e EXT_ANNOT] [-ch CHART_PATH]
Create the scrna report.

    -h, --help                   display this help and exit
    -p, --project PROJECT_NAME   name of your project
    -s, --save_path SAVE_PATH    path where the objects of your analysis were saved, usually "save_[PROJECT_NAME]"
    -c, --config_path CONF_PATH  path to the config file, usually conf/config_[PROJECT_NAME].R
    -ch, --chart_path CHART_PATH path where the chart data of your analysis was saved, usually "chart_[PROJECT_NAME]"
    -e, --ext_annot EXT_ANNOT    path to the external annotation file
    -o, --output OUTPUT_FOLDER   optional, default is "./report_[PROJECT_NAME]", path to the output folder
    -a, --author AUTHOR          optional, default is "CostaLab", define the author of the analysis, this name MUST NOT have spaces
    -m, --make_elements          optional flag, if set, make all report elements (plots/figures) from scratch
    -so, --one_file              optional flag, if set, in addition to the default report, also produce the one file report
EOF
}


while (( $# )); do
  case $1 in
    -h|-\?|--help)
      show_help    # Display a usage synopsis.
      exit 1
      ;;
    -p|--project)
      if [ "$2" ]; then
        PROJECT="$2"
        shift # past argument
        shift # past value
      else
        printf 'ERROR: "-p, --project" requires a non-empty option argument.'
        exit 1
      fi
    ;;
    -a|--author)
      AUTHOR="$2"
      shift # past argument
      shift # past value
    ;;
    -o|--output)
      OUTPUT_FOLDER="$2"
      shift # past argument
      shift # past value
    ;;
    -m|--make_elements)
      MAKE_ELEMENTS=1
      shift # past argument
    ;;
    -so|--one_file)
      MAKE_SINGLE_FILE=1
      shift # past argument
    ;;
    -s|--save_path)
      if [ "$2" ]; then
        SAVE_DATA_PATH="$2"
        shift # past argument
        shift # past value
      else
        printf 'ERROR: "-s, --save_path" requires a non-empty option argument.'
        exit 1
      fi
    ;;
    -c|--config_path)
      if [ "$2" ]; then
        CONFIG_FILE_PATH="$2"
        shift # past argument
        shift # past value
      else
        printf 'ERROR: "-c, --config_path" requires a non-empty option argument.'
        exit 1
      fi
    ;;
    -ch|--chart_path)
      if [ "$2" ]; then
        CHART_DATA_PATH="$2"
        shift # past argument
        shift # past value
      else
        printf 'ERROR: "-ch, --chart_path" requires a non-empty option argument.'
        exit 1
      fi
    ;;
    -e|--ext_annot)
      if [ "$2" ]; then
        EXT_ANNOT_PATH="$2"
        shift # past argument
        shift # past value
      else
        printf 'ERROR: "-e, --ext_annot" requires a non-empty option argument.'
        exit 1
      fi
    ;;
    -?*|-*|--*)
      printf 'ERROR: Unknown option: %s\n' "$1" >&2
      exit 1
    ;;
    *)    # unknown option
      shift # past argument
    ;;
  esac
done


if [ -z "$PROJECT" ]; then
  printf 'ERROR: "-p, --project" need to be defined.'
  exit 1
fi

if [ -z "$CONFIG_FILE_PATH" ]; then
  printf 'ERROR: "-c, --config_path" need to be defined.'
  exit 1
fi

if [ -z "$SAVE_DATA_PATH" ]; then
  printf 'ERROR: "-s, --save_path" need to be defined.'
  exit 1
fi

if [ -z "$CHART_DATA_PATH" ]; then
  printf 'ERROR: "-c, --chart_path" need to be defined.'
  exit 1
fi

if [ -z "$EXT_ANNOT_PATH" ]; then
  printf 'ERROR: "-e, --ext_annot" need to be defined.'
  exit 1
fi

if [[ "$MAKE_SINGLE_FILE" == 1 ]]; then
  MAKE_SINGLE_FILE="TRUE"
else
  MAKE_SINGLE_FILE="FALSE"
fi

if [ -z "$OUTPUT_FOLDER" ]; then
  OUTPUT_FOLDER=($DEFAULT_OUTPUT_FOLDER$PROJECT)
fi


RED='\033[0;31m'
NC='\033[0m' # No Color

FUNCS=(
  #Singleton
  QC
  DEs
  Clusters
  Clusters_harmony
  Clusters_seurat
  EXT_MARKERS
  DEGO
  hallmark
  KEGG
  Reactome
  DEGO_stage
  hallmark_stage
  reactome_stage
  kegg_stage
  DEGO_1v1
  hallmark_1v1
  reactome_1v1
  kegg_1v1
  #intUMAPs
  #GET_DATA
)

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }
join_str=`join_by '", "' ${FUNCS[@]}`
py_exe_list=\[\"${join_str}\"\]


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
#cluster="singleton"
cluster="seurat_clusters"
#!!!!!!!!!----------------------------------

echo $PROJECT
echo -e "Use cluster slot ${RED} $cluster ${NC}"
mkdir -p $OUTPUT_FOLDER/data
cp -r $CHART_DATA_PATH/* $OUTPUT_FOLDER/data

python viz/code_generator.py \
  -c $cluster \
  -cf $CONFIG_FILE_PATH \
  --save_dir "${SAVE_DATA_PATH}" \
  --output_dir "${OUTPUT_FOLDER}" \
  --proj_tag "${PROJECT}" \
  --executing_list "${py_exe_list}"
  # --base_dir "../" \

grip --export $OUTPUT_FOLDER/index.md


# Make report elements
if [[ "$MAKE_ELEMENTS" == 1 ]]; then
  Rscript viz/make_report_elements.R \
    --proj=$PROJECT \
    --cluster=$cluster \
    --save_dir=$SAVE_DATA_PATH \
    --ext_annot=$EXT_ANNOT_PATH \
    --output_dir=$OUTPUT_FOLDER \
    "${FUNCS[@]}"
fi
if [ $? -ne 0 ]; then
  echo "make_report_elements.R has not finished successfully."
  exit 1
fi

# Produce report
Rscript viz/make_report.R \
  --proj=$PROJECT \
  --cluster=$cluster \
  --save_dir=$SAVE_DATA_PATH \
  --make_single_file=$MAKE_SINGLE_FILE \
  --author=$AUTHOR \
  --output_dir=$OUTPUT_FOLDER \
  "${FUNCS[@]}"
if [ $? -ne 0 ]; then
  echo "make_report.R has not finished successfully."
  exit 1
fi
