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

## GENERATE  R markdown and md according to the config file:
#1.  dego 1v1
#2. dego stages vs
#3. index.md

parser = argparse.ArgumentParser(description='Process arguments.')
parser.add_argument('-c', '--cluster_use', default='seurat_clusters', help='clusters to choose')
parser.add_argument('-cf', '--config_file', default='conf/config.R', help='path to config file')
parser.add_argument('-b', '--base_dir', default='../', help='path to base dir')
parser.add_argument('-p', '--proj_tag', default='', help='project name or tag')
args = parser.parse_args()


DATADIR = args.base_dir

robjects.r['source'](os.path.join(DATADIR, args.config_file))  ### load config
names = robjects.r("names(data_src)")
lst_1v1 = list(combinations(names, 2))
stages = robjects.r("stage_lst")
project_name = robjects.r("PROJECT")


seen = set()
u_stages = [x for x in stages if x not in seen and not seen.add(x)]  ##remove dup
lst_stages = list(combinations(u_stages, 2))

cluster_use = "seurat_clusters"
savedir = os.path.join(DATADIR, "save"+args.proj_tag)#robjects.r("SAVE_DIR")[0]

# try:
#     options,args = getopt.getopt(sys.argv[1:],"c:")
# except getopt.GetoptError:
#     print("Error Parameters")
#     sys.exit()
# for name,value in options:
#     if name in "-c":
#         cluster_use = value
#         print("cluster use:", cluster_use)

def generate_1v1(out):

    fw = open(out, 'w')


    thead = open("template/head.template")
    tmplth = thead.read()
    t = Template(tmplth)
    today = datetime.date.today().strftime("%d%B%Y")
    head = t.render(today = today,
                    svdir = savedir)

    fw.write("%s\n" %head)

    tfile = open("template/DE-GO-1v1.template")
    tmpltr = tfile.read()
    for x,y in lst_1v1:
        vs = "%s vs %s" % (x, y)
        fde =  "%s.vs.%s.de.Rdata" % (x, y)
        fgoup = "%s.vs.%s.goup.Rdata" % (x, y)
        fgodown = "%s.vs.%s.godown.Rdata" % (x, y)
        xlsup = "Go.UP.%s.vs.%s.xlsx" % (x, y)
        xlsdown = "Go.Down.%s.vs.%s.xlsx" % (x, y)

        t = Template(tmpltr)

        one_plots = t.render(Versus= vs,
                             de_file= fde,
                             goup_file= fgoup,
                             goup_xlsx= xlsup,
                             godown_file= fgodown,
                             godown_xlsx= xlsdown,
                             tX = x,
                             tY = y)


        fw.write("%s\n" %one_plots)

#endf generate_1v1



def generate_md_idx(out):
    tfile = open("template/index.template")
    tmplt = tfile.read()
    t = Template(tmplt)
    st = t.render(list_1v1=lst_1v1,
                  list_stages=lst_stages,
                  project_name=project_name[0],
                  cluster_use=args.cluster_use)

    fw = open(out, 'w')
    fw.write(st)
#endf generate_md_idx

def generate_md_idx_pdf(out, folder="report_pdf"):

    if not os.path.exists(folder):
        os.mkdir(folder)

    tfile = open("template/index_pdf.template")
    tmplt = tfile.read()
    t = Template(tmplt)
    st = t.render(list_1v1=lst_1v1,
                  list_stages=lst_stages,
                  project_name=project_name[0],
                  cluster_use = args.cluster_use)

    fw = open(out, 'w')
    fw.write(st)
#endf generate_md_idx




def generate_report_stagesVS(out):
    tfile = open("template/DE-GO-stagesVS.template")
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


def generate_report_1v1(out):
    tfile = open("template/DE-GO-1v1.template")
    tmpl = tfile.read()
    today = datetime.date.today().strftime("%d%B%Y")
    cc = 6

    t = Template(tmpl)
    r = t.render(
        TODAY=today,
        CC=cc,
        savedir=savedir,
        lst_1v1=lst_1v1
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS


def generate_report_1v1_hallmark(out):
    tfile = open("template/hallmark-1v1.template")
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


def generate_report_stageVS_hallmark(out):
    tfile = open("template/hallmark-stageVS.template")
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
    tfile = open("template/reactome-1v1.template")
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
    tfile = open("template/reactome-stageVS.template")
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




def generate_report_1v1_kegg(out):
    tfile = open("template/kegg-1v1.template")
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


def generate_report_stageVS_kegg(out):
    tfile = open("template/kegg-stageVS.template")
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
    out_dir = "report"+args.proj_tag
    generate_report_1v1("DE-GO-analysis-1v1.Rmd")
    generate_report_1v1_hallmark("hallmark-1v1.Rmd")
    generate_report_1v1_reactome("reactome-1v1.Rmd")
    generate_report_1v1_kegg("kegg-1v1.Rmd")
    generate_report_stagesVS("DE-GO-analysis-stagesVS.Rmd")
    generate_report_stageVS_hallmark("hallmark-stageVS.Rmd")
    generate_report_stageVS_reactome("reactome-stageVS.Rmd")
    generate_report_stageVS_kegg("kegg-stageVS.Rmd")

    generate_md_idx(os.path.join(out_dir, "index.md"))

    out_dir = "report_pdf"
    generate_md_idx_pdf(os.path.join(out_dir, "index.md"), out_dir)
#endf main


if __name__ == "__main__":
    main()
