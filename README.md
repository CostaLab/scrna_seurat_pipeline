# scRNA pipeline

### Documents:
https://costalab.github.io/scrna_seurat_pipeline

### A docker image under seurat4
https://costalab.ukaachen.de/open_data/scrna_seurat_pipeline/scrna_pipeline_docker.tar

```shell
docker load < scrna_pipeline_docker.tar
docker run -it docker.io/library/seurat4.0


#directory mapping to link directories outside the container

docker run -v /PATH/conf/:/app/conf -v /PATH/Matrix/:/app/input -it docker.io/library/seurat4.0
```
