---
sort: 1
---

# Getting Started


## System Requirements
OS: `linux`, `MacOS`
R: https://www.r-project.org/
Python: `>3.6`

## Clone repository to local

```shell
git clone https://github.com/CostaLab/scrna_seurat_pipeline.git
```

## Dependent libraries installation

```shell
Rscript packages_install.R
```

## Run Toy example

### data producing
This step will take around 1 hour, depends on the hardware performance.

```shell
cd scrna_seurat_pipeline
cp conf/config.R conf/config_toys.R ## toy configuration file
cp run_example.sh run_toy.sh ## Please edit run_toy.sh to fit your environment
sh run_toy.sh toy
```

### Run visualization
This step will take around 40 minutes, depends on the hardware performance.

```shell
cp run_viz_example.sh run_viz_toy.sh ## Please edit run_viz_toy.sh to fit your environment
sh run_viz_toy.sh toy
```

## Check the report
Go to `toy/report/`, open index.html, you can check the result


