### Human gut cell atlas analysis ##   srun -A bch -p bch-largemem -n 1 --qos=largemem --mem=180GB -t 0-06:00 --pty /bin/bash

library(data.table)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(data.table)
library(tidyverse)
library(ggplot2)
library(devtools)
library(harmony)
library(patchwork)
library(viridis)
###Function to read h5ad file, create h5Seurat and recover metadata from dataset
getdata<-function(h5file){
    h5Seurat<-h5file%>%str_replace(".h5ad$",".h5Seurat")
    print(h5Seurat)
    if(file.exists(h5Seurat)==FALSE){
        print("Converting h5file to h5Seurat")
            Convert(h5file, h5Seurat)
    }
    hfile<-Connect(paste0(h5Seurat))
    seurat_object<-LoadH5Seurat(h5Seurat,misc=FALSE)
    seurat_object
}

scdrs_dir="/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/scdrs_analysis"
project="humancellgut_traits"

h5file="/lab-share/IM-Gutierrez-e2/Public/References/humancellatlas/Full_obj_raw_counts_nosoupx_v2.h5ad"
seurat_object<- getdata(h5file)

#LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
#This is then natural-log transformed using log1p.
seurat_object <- NormalizeData(object = seurat_object)
genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", row.names(seurat_object), value = TRUE)
#FindVariable features
seurat_object <- FindVariableFeatures(object = seurat_object )
#Exclude mitocondrial and ribosomal 
myfeatures<-seurat_object@assays$RNA@var.features[!seurat_object@assays$RNA@var.features%in%genes_exclude]
myfeatures%>%length() # 
seurat_object <- ScaleData(object = seurat_object,features = myfeatures)
seurat_object <- RunPCA(object = seurat_object,features = myfeatures)
#UMAP no harmonized
umap_no_harmonized<- RunUMAP(object = seurat_object,dims = 1:20)%>%
  DimPlot(group.by =  "batch",reduction = "umap",pt.size = 1.5,shuffle = TRUE)
umap_no_harmonized<-umap_no_harmonized+ theme(legend.position="none")
#### Harmony #####
# Harmony, tutorial https://portals.broadinstitute.org/harmony/SeuratV3.html ##
#group.by.vars Which variable(s) to remove (character vector).
seurat_object@meta.data$batch<-seurat_object@meta.data$batch%>%droplevels() #droplevels otherwise you get incompatible matrix dimensions for subsets
seurat_object<-harmony::RunHarmony(object = seurat_object,group.by.vars = "batch",plot_convergence=FALSE)
seurat_object <- RunUMAP(object = seurat_object,reduction = "harmony",dims = 1:20)
#LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
seurat_object <- FindNeighbors(object = seurat_object,dims = 1:20)
seurat_object <- FindClusters(object = seurat_object)
umap_harmonized<- DimPlot(object = seurat_object, group.by =  "batch",reduction = "umap",pt.size = 1.5,shuffle = TRUE)+theme(legend.position="none")

umaps<- umap_no_harmonized+umap_harmonized
fs::dir_create(fs::path(scdrs_dir,project,"scplots"))
ggsave(umaps,filename = fs::path(scdrs_dir,project,"scplots/umaps.png"),width = 12,height = 12)

Embeddings_UMAP_1<-Embeddings(object = seurat_object, reduction =  "umap")[,1]
Embeddings_UMAP_2<-Embeddings(object = seurat_object, reduction =  "umap")[,2]
Embeddings_UMAP<-Embeddings(object = seurat_object, reduction =  "umap")[,1:2]
head(Embeddings_UMAP)
fwrite(Embeddings_UMAP%>%as.data.frame()%>%rownames_to_column("cell_id"),file = fs::path(scdrs_dir,project,"UMAP_coordinates.tsv"))
seurat_object<-AddMetaData(object = seurat_object, metadata = Embeddings_UMAP_1,col.name="UMAP_1")
seurat_object<-AddMetaData(object = seurat_object, metadata = Embeddings_UMAP_2,col.name="UMAP_2")
meta.data<-seurat_object@meta.data
head(meta.data)
fwrite(meta.data%>%rownames_to_column("cell_id"),file= fs::path(scdrs_dir,project,"metadata.tsv"),sep="\t")
SaveH5Seurat(seurat_object, filename = "/lab-share/IM-Gutierrez-e2/Public/References/humancellatlas/gutcellatlas_harmonized.h5Seurat")


h5file="/lab-share/IM-Gutierrez-e2/Public/References/humancellatlas/gutcellatlas_harmonized.h5ad"
seurat_object<- getdata(h5file)


