#!/bin/python3

import datetime
import os
import sys
import argparse
import getopt
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from itertools import combinations
from jinja2 import Template
from requests.utils import requote_uri

def intersect(a, b):
    return set(a).intersection(b)

# GENERATE  R markdown and md according to the config file:
# 1. dego 1v1
# 2. dego stages vs
# 3. index.md


parser = argparse.ArgumentParser(description='Process arguments.')
parser.add_argument('-c', '--cluster_use', default='seurat_clusters', help='clusters to choose')
parser.add_argument('-cf', '--config_file', default='conf/config.R', help='path to config file')
# parser.add_argument('-b', '--base_dir', default='../', help='path to base dir')
parser.add_argument('-o', '--output_dir', default='../', help='path to output dir')
parser.add_argument('-s', '--save_dir', default='../', help='path to save dir')
parser.add_argument('-p', '--proj_tag', default='', help='project name or tag')
parser.add_argument('-l', '--executing_list', default='[]', help='executing list to viz')
args = parser.parse_args()

exe_list = [i.strip() for i in eval(args.executing_list)]


print(
    "==========PLOT LIST===============:\n",
    "\n".join(exe_list),
    "\n==================================\n"
)

viz_dict = {
    "quality": ["QC", "QCC", "AmbientRNA"],
    "batch_clustering": ["Clusters_seurat",
                         "Clusters_harmony",
                         "DEs"],

    "clustering": ["Clusters"],

    "clustersVS": ["EXT_MARKERS",
                   "EXTRA_EXT_MARKERS",
                   "DEGO",
                   "Genesets",
                   "progeny",
                   "hallmark",
                   "KEGG",
                   "Reactome",
                   "intUMAPs"],

    "DEGOstageVS": ["DEGO_stage"],

    "DEGOextrastageVS": ["DEGO_extrastage"],

    "PWstageVS": ["Genesets_stage",
                  "progeny_stage",
                  "hallmark_stage",
                  "reactome_stage",
                  "kegg_stage"],

    "PWextrastageVS": ["hallmark_extrastage",
                       "reactome_extrastage",
                       "kegg_extrastage"],

    "DEGOsampleVS": ["DEGO_1v1"],

    "PWsampleVS": ["Genesets_1v1",
                   "hallmark_1v1",
                   "reactome_1v1",
                   "kegg_1v1"]
}


# DATADIR = args.base_dir

robjects.r['source'](args.config_file)  # load config
names = robjects.r("names(data_src)")
lst_1v1 = list(combinations(names, 2))
stages = robjects.r("stage_lst")
project_name = robjects.r("PROJECT")
genesets_names = robjects.r("MSigDB_Geneset_names")
integration_option = robjects.r("INTEGRATION_OPTION")[0]

seen = set()
u_stages = [x for x in stages if x not in seen and not seen.add(x)]  # remove dup
lst_stages = list(combinations(u_stages, 2))

try:
    robjects.r('if (!exists("extra_stage_cols")) extra_stage_cols <- NULL')
    esc_r = robjects.r("extra_stage_cols")
    if esc_r != robjects.NULL:
        robjects.r('if (!exists("clinical_meta")) clinical_meta <- NULL')
        extra_stage_cols = {}
        for col_name in list(esc_r.names):
            defn = esc_r.rx2(col_name)
            if isinstance(defn, robjects.ListVector):
                extra_stage_cols[col_name] = list(defn.rx2("labels"))
            else:
                cm = robjects.r("clinical_meta")
                if cm != robjects.NULL:
                    vals = list(set(str(v) for v in cm.rx2(str(defn))))
                    extra_stage_cols[col_name] = vals
                else:
                    extra_stage_cols[col_name] = []
    else:
        extra_stage_cols = {}
except Exception:
    extra_stage_cols = {}

lst_extrastages = {}
for col_name, groups in extra_stage_cols.items():
    if len(groups) >= 2:
        lst_extrastages[col_name] = list(combinations(groups, 2))

# FIXME should 'cluster_use' be redefined here?
# cluster_use = "seurat_clusters"
# savedir = os.path.join(DATADIR, "save"+args.proj_tag)#robjects.r("SAVE_DIR")[0]
savedir = args.save_dir

# try:
#     options,args = getopt.getopt(sys.argv[1:],"c:")
# except getopt.GetoptError:
#     print("Error Parameters")
#     sys.exit()
# for name,value in options:
#     if name in "-c":
#         cluster_use = value
#         print("cluster use:", cluster_use)


def generate_md_idx(out):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "index.template"))
    tmplt = tfile.read()
    t = Template(tmplt)
    st = t.render(requote_uri=requote_uri,
                  intersect=intersect,
                  executing_list=exe_list,
                  viz_dict=viz_dict,
                  list_1v1=lst_1v1,
                  list_stages=lst_stages,
                  lst_extrastages=lst_extrastages,
                  integr_option=integration_option,
                  project_name=project_name[0],
                  cluster_use=args.cluster_use)

    fw = open(out, 'w')
    fw.write(st)
# endf generate_md_idx


def generate_report_stagesVS(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "DE-GO-vs.template"))
    tmpl = tfile.read()
    # today = datetime.date.today().strftime("%d%B%Y")

    t = Template(tmpl)
    for x, y in lst_stages:
        r = t.render(
            tX=x,
            tY=y,
            group="Stages"
        )
        fw = open(os.path.join(viz_path, f"4_DE_GO_{x}.vs.{y}_stageVS.Rmd"), "w")
        fw.write("%s\n\n" % r)
        fw.close()
    # end for
# endf generate_groupVS


def generate_report_1v1(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "DE-GO-vs.template"))
    tmpl = tfile.read()
    # today = datetime.date.today().strftime("%d%B%Y")

    t = Template(tmpl)
    for x, y in lst_1v1:
        r = t.render(
            tX=x,
            tY=y,
            group="Samples"
        )
        fw = open(os.path.join(viz_path, f"4_DE_GO_{x}.vs.{y}_1v1.Rmd"), "w")
        fw.write("%s\n\n" % r)
        fw.close()
    # end for
# endf generate_report_1v1


def generate_report_stagesVS_pw(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "pathway_vs.template"))
    tmpl = tfile.read()

    for pw in ["hallmark", "kegg", "reactome"]:
        t = Template(tmpl)
        r = t.render(
            pathway=pw,
            group="Stages",
            lst_group=lst_stages
        )
        fw = open(os.path.join(viz_path, f"4_{pw}_stageVS.Rmd"), "w")
        fw.write("%s\n\n" % r)
        fw.close()
    # endfor
# endf generate_groupVS

def generate_report_1v1_pw(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "pathway_vs.template"))
    tmpl = tfile.read()

    for pw in ["hallmark", "kegg", "reactome"]:
        t = Template(tmpl)
        r = t.render(
            pathway=pw,
            group="1v1",
            lst_group=lst_1v1
        )
        fw = open(os.path.join(viz_path, f"4_{pw}_1v1.Rmd"), "w")
        fw.write("%s\n\n" % r)
        fw.close()
    # endfor
# endf generate_report_1v1_pw



def generate_report_1v1_Genesets(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "Genesets-vs.template"))
    tmpl = tfile.read()

    t = Template(tmpl)
    r = t.render(
            group="1v1",
            lst_group=lst_1v1
    )
    fw = open(os.path.join(viz_path, "4_Genesets_1v1.Rmd"), "w")
    fw.write("%s\n\n" % r)
    fw.close()
# endf generate_1v1


def generate_report_stageVS_Genesets(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "Genesets-vs.template"))
    tmpl = tfile.read()

    t = Template(tmpl)
    r = t.render(
            group="stageVS",
            lst_group=lst_stages
    )
    fw = open(os.path.join(viz_path, "4_Genesets_stageVS.Rmd"), "w")
    fw.write("%s\n\n" % r)
    fw.close()
# endf generate_stageVS



def generate_report_stageVS_progeny(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "progeny-vs.template"))
    tmpl = tfile.read()

    t = Template(tmpl)
    r = t.render(
        group="stageVS",
        lst_group=lst_stages
    )

    fw = open(os.path.join(viz_path, "4_progeny_stageVS.Rmd"), "w")
    fw.write("%s\n\n" % r)
# endf generate_groupVS


def generate_report_extrastagesVS(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "DE-GO-extrastage-vs.template"))
    tmpl = tfile.read()

    for col_name, pairs in lst_extrastages.items():
        prefix = f"extrastage_{col_name}_"
        t = Template(tmpl)
        for x, y in pairs:
            r = t.render(
                tX=x,
                tY=y,
                col_name=col_name,
                prefix=prefix
            )
            fw = open(os.path.join(viz_path, f"4_DE_GO_extrastage_{col_name}_{x}.vs.{y}.Rmd"), "w")
            fw.write("%s\n\n" % r)
            fw.close()
# endf generate_report_extrastagesVS


def generate_report_extrastagesVS_pw(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__), "template", "pathway_extrastage-vs.template"))
    tmpl = tfile.read()

    for pw in ["hallmark", "kegg", "reactome"]:
        for col_name, pairs in lst_extrastages.items():
            prefix = f"extrastage_{col_name}_"
            t = Template(tmpl)
            for x, y in pairs:
                r = t.render(
                    tX=x,
                    tY=y,
                    col_name=col_name,
                    group=col_name,
                    pathway=pw,
                    prefix=prefix
                )
                fw = open(os.path.join(viz_path, f"4_{pw}_extrastage_{col_name}_{x}.vs.{y}.Rmd"), "w")
                fw.write("%s\n\n" % r)
                fw.close()
# endf generate_report_extrastagesVS_pw


def main():

    out_dir = args.output_dir
    viz_dir = os.path.dirname(__file__)

    generate_report_1v1(viz_dir)
    generate_report_stagesVS(viz_dir)
    generate_report_extrastagesVS(viz_dir)
    generate_report_extrastagesVS_pw(viz_dir)

    generate_report_1v1_pw(viz_dir)
    generate_report_stagesVS_pw(viz_dir)

    generate_report_1v1_Genesets(viz_dir)
    generate_report_stageVS_Genesets(viz_dir)

    generate_report_stageVS_progeny(viz_dir)

    generate_md_idx(os.path.join(out_dir, "index.md"))
# endf main


if __name__ == "__main__":
    main()

