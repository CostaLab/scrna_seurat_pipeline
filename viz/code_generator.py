#!/bin/python3

from jinja2 import Template
import datetime
import os
import sys 
import getopt
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from itertools import combinations 

## GENERATE  R markdown and md according to the config file:
#1.  dego 1v1 
#2. dego stages vs
#3. index.md 

DATADIR = "../"

robjects.r['source'](os.path.join(DATADIR,"conf/config.R")) ### load config
names = robjects.r("names(data_src)")  
lst_1v1 = list(combinations(names, 2))
stages = robjects.r("stage_lst")  
seen = set()
u_stages =  [x for x in stages if x not in seen and not seen.add(x)] ##remove dup
lst_stages = list(combinations(u_stages, 2))

cluster_use = "seurat_clusters"

try:
    options,args = getopt.getopt(sys.argv[1:],"c:")
except getopt.GetoptError:
    print("Erorr Parametes")
    sys.exit()
for name,value in options:
    if name in "-c":
        cluster_use = value
        print("cluster use:", cluster_use) 









savedir = os.path.join(DATADIR, "save/")#robjects.r("SAVE_DIR")[0] 

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
                  cluster_use = cluster_use)
    
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
        TODAY = today,
        CC = cc,
        lst_stages = lst_stages
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
        TODAY = today,
        CC = cc,
        lst_1v1 = lst_1v1
    )
    fw = open(out, "w")
    fw.write("%s\n\n" % r)
#endf generate_groupVS




def main():
    out_dir = "report"
    generate_report_1v1("DE-GO-analysis-1v1.Rmd")
    generate_report_stagesVS("DE-GO-analysis-stagesVS.Rmd")
    generate_md_idx(os.path.join(out_dir, "index.md"))
#endf main


if __name__ == "__main__":
    main()
