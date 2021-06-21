#!/bin/python3

from jinja2 import Template
import datetime
import os
import sys
import argparse
import getopt
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from itertools import combinations
from requests.utils import requote_uri
def intersect(a, b):
    return set(a).intersection(b)


## GENERATE  R markdown and md according to the config file:
#1.  dego 1v1
#2. dego stages vs
#3. index.md

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


print("==========PLOT LIST===============:\n",
        "\n".join(exe_list),
       "\n==================================\n")

viz_dict = { "quality": ["QC", "AmbientRNA"],
             "batch_clustering": ["Clusters_seurat",
                                  "Clusters_harmony",
                                  "DEs"],

             "clustering": ["Clusters"],

             "clustersVS": ["EXT_MARKERS",
                            "DEGO",
                            "Genesets",
                            "progeny",
                            "hallmark",
                            "KEGG",
                            "Reactome",
                            "intUMAPs"],

             "DEGOstageVS": ["DEGO_stage"],

             "PWstageVS": ["Genesets_stage",
                           "progeny_stage",
                           "hallmark_stage",
                           "reactome_stage",
                           "kegg_stage"],

             "DEGOsampleVS":["DEGO_1v1"],

             "PWsampleVS": ["Genesets_1v1",
                            "hallmark_1v1",
                            "reactome_1v1",
                            "kegg_1v1"]
}


# DATADIR = args.base_dir

robjects.r['source'](args.config_file)  ### load config
names = robjects.r("names(data_src)")
lst_1v1 = list(combinations(names, 2))
stages = robjects.r("stage_lst")
project_name = robjects.r("PROJECT")
genesets_names = robjects.r("MSigDB_Geneset_names")
integration_option = robjects.r("INTEGRATION_OPTION")[0]

seen = set()
u_stages = [x for x in stages if x not in seen and not seen.add(x)]  ##remove dup
lst_stages = list(combinations(u_stages, 2))

# FIXME should 'cluster_use' be redefined here?
#cluster_use = "seurat_clusters"
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
    tfile = open(os.path.join(os.path.dirname(__file__),"template/index.template"))
    tmplt = tfile.read()
    t = Template(tmplt)
    st = t.render(requote_uri=requote_uri,
                  intersect=intersect,
                  executing_list=exe_list,
                  viz_dict=viz_dict,
                  list_1v1=lst_1v1,
                  list_stages=lst_stages,
                  integr_option=integration_option,
                  project_name=project_name[0],
                  cluster_use=args.cluster_use)

    fw = open(out, 'w')
    fw.write(st)
#endf generate_md_idx


def generate_report_stagesVS(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/DE-GO-stagesVS.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY=today,
        CC=cc,
        savedir=savedir,
        lst_stages=lst_stages
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS


def generate_report_1v1(viz_path):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/DE-GO-1v1.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")

    t = Template(tmpl)
    for x, y in lst_1v1:
        r = t.render(
            tX = x,
            tY = y
        )
        fw = open(f"{viz_path}/4_DE_GO_{x}.vs.{y}_1v1.Rmd", "w")
        fw.write("%s\n\n" % r)
        fw.close()
    #end for
#endf generate_groupVS


def generate_report_1v1_hallmark(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/hallmark-1v1.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_1v1 = lst_1v1
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS

def generate_report_stageVS_Genesets(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/Genesets-stageVS.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_stage = lst_stages
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS




def generate_report_stageVS_hallmark(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/hallmark-stageVS.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_stage = lst_stages
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS


def generate_report_1v1_reactome(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/reactome-1v1.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_1v1 = lst_1v1
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS

def generate_report_stageVS_reactome(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/reactome-stageVS.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_stage = lst_stages
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS



def generate_report_1v1_Genesets(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/Genesets-1v1.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_1v1 = lst_1v1
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_1v1




def generate_report_1v1_kegg(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/kegg-1v1.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_1v1 = lst_1v1
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_1v1


def generate_report_stageVS_kegg(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/kegg-stageVS.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_stage = lst_stages
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS


def generate_report_stageVS_progeny(out):
    tfile = open(os.path.join(os.path.dirname(__file__),"template/progeny-stageVS.template"))
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY = today,
        CC = cc,
        lst_stage = lst_stages
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS







def main():


    out_dir = args.output_dir
    viz_dir = os.path.dirname(__file__)
    generate_report_1v1(viz_dir)
    generate_report_1v1_Genesets(os.path.join(viz_dir,"Genesets-1v1.Rmd"))
    generate_report_1v1_hallmark(os.path.join(viz_dir,"hallmark-1v1.Rmd"))
    generate_report_1v1_reactome(os.path.join(viz_dir,"reactome-1v1.Rmd"))
    generate_report_1v1_kegg(os.path.join(viz_dir,"kegg-1v1.Rmd"))
    generate_report_stagesVS(os.path.join(viz_dir,"DE-GO-analysis-stagesVS.Rmd"))
    generate_report_stageVS_hallmark(os.path.join(viz_dir,"hallmark-stageVS.Rmd"))
    generate_report_stageVS_reactome(os.path.join(viz_dir,"reactome-stageVS.Rmd"))
    generate_report_stageVS_kegg(os.path.join(viz_dir,"kegg-stageVS.Rmd"))
    generate_report_stageVS_progeny(os.path.join(viz_dir,"progeny-stageVS.Rmd"))
    generate_report_stageVS_Genesets(os.path.join(viz_dir,"Genesets-stageVS.Rmd"))
    generate_md_idx(os.path.join(out_dir, "index.md"))

#endf main


if __name__ == "__main__":
    main()
