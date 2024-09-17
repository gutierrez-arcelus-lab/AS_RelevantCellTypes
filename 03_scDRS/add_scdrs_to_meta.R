#!/usr/bin/env Rscript
#
# Function: Analize scdrs results on single cell object
# Author: Marcos Chinas
# Last modification: 24/03/2022
#
# Usage:
#   Rscript addres_scdrs.R [scdrs_results] [h5Seurat] [Celltype]
#   
# Input files: 
# -*.all.tsv.gz - summary statistics from eQTL catalogue
# 
# ------------------------------------------------------------------------------
#magma_v1.10/AS_ErmannLab/scdrs_results/default/ /Users/marcosch/Downloads/AS/AS_ErmannLab/seurat/h5ad_files/AS_seurat.h5Seurat SeuratLabel AS_ErmannLab
#Design to read metadata and scdrs results

suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
  make_option(c("-r", "--results"), type ="character" ,default = FALSE,help="Path or file with scdrs results"),
  make_option(c("-i", "--h5seurat"),type ="character" ,default = FALSE,help="File with metadata for h5Seurat"),
  make_option(c("-o", "--outputdir"),type ="character", default=".", help="Dir for outputs"),
  make_option(c("-c", "--celltype"),type ="character" ,default = FALSE, help = "Celltype column on metadata"),
  make_option(c("-s", "--splitby"),type ="character" , default=FALSE,help = "Optional column to split plots and tables"))
args <- parse_args(OptionParser(option_list=option_list))

#stopifnot(length(args) >= 4)
scdrs_path <- args$results
h5Seurat <- args$h5seurat
celltype <- args$celltype
out      <- args$outputdir
split.by <- args$splitby
#print(h5Seurat)
library(Seurat)
library(harmony)
library(tidyverse)
library(SeuratDisk)
library(fs)
library(data.table)
fs::dir_create(out)

#hfile<-Connect(paste0(h5Seurat))
meta.data<-fread(h5Seurat)
head(meta.data)
#ncol(seurat_object)
#head(seurat_object)
#scdrs_path<-"~/magma_v1.10/AS_ErmannLab/scdrs_results/default/"
#out<-"AS_ErmannLab/" 
#celltype<- "SeuratLabel"
#split.by="Status"
#scdrs_path <-"/lab-share/IM-Gutierrez-e2/Public/References/bm_atlas/scdrs_analysis/bm/scdrs_results/"
if(fs::is_dir(scdrs_path)==TRUE){
  files<-dir_ls(scdrs_path,regexp = "\\.score.gz",recurse = TRUE)
  }else{
  files<-scdrs_path
}

add_scdrs<-function(path,metadata,trait="",stats=TRUE){
  library(qvalue)
  if(str_length(trait)>0){
    itrait<-trait
    trait<-paste0("_",trait)
  }
  for(path in path){
  itrait<-basename(path)%>%str_remove("\\..+")
  trait<-paste0("_",basename(path)%>%str_remove("\\..+"))

  scdrs_res<-fread(path)
  colnames(scdrs_res)[1]<-"cell_id"
  #print(trait)
  colnames(scdrs_res)<-paste0(colnames(scdrs_res),trait)
  colnames(scdrs_res)
  pvalname<-paste0("pval",trait)
  fdrname <-paste0("fdr",trait)
  fdr<-qvalue(p=scdrs_res[[pvalname]])
  scdrs_res[[fdrname]]<-fdr$qvalues
  if(stats==TRUE){
    fdr10<-sum(scdrs_res[[fdrname]]<0.1)
    fdr20<-sum(scdrs_res[[fdrname]]<=0.2)
    print(paste("Trait:",itrait,"  cells FDR 10%:",fdr10, "   cells FDR 20%:",fdr20))
  }
  #print(head(scdrs_res))
  metadata<-merge(metadata,scdrs_res,by.x="cell_id",by.y=paste0("cell_id",trait))#%>%column_to_rownames(var = "Row.names")
  }
  metadata
}

meta.data<-add_scdrs(path=files,metadata = meta.data)

#meta.data[[celltype]]%>%aggregate(by=list(celltype),FUN=length)


traits<-paste0(basename(files)%>%str_remove("\\..+"))


p_scdrs<-function(meta.data,split.by=NULL,trait="",threshold=0.1){
  if(str_length(trait)>0){
    itrait<-trait
    trait<-paste0("_",trait)
  }
  meta.data$filtered_norm_score<-NA
  fdrname <-paste0("fdr",trait)
  meta.data$filtered_norm_score[meta.data[[fdrname]]<threshold]<-meta.data[[paste0("norm_score",trait)]][meta.data[[fdrname]]<threshold]
  if("umap"%in%colnames(seurat_object@meta.data)){
  p<-FeaturePlot(seurat_object,features = "filtered_norm_score",label = FALSE,pt.size = 2,split.by = split.by,reduction = "umap",order = TRUE)&
    scale_color_viridis_c(option = "viridis",direction = -1,na.value = "grey")&
    theme(legend.position = c(0.88,0.2))
  print(p)
  ggsave(filename =fs::path(out,paste0("FeaturePlot_","FDR_",threshold,trait,".png")),plot = p)
  }
}

for(trait in traits){
#p_scdrs(seurat_object,split.by = split.by,trait =trait ,threshold = 0.1)
}

mosaic_scdrs_2<-function(meta.data,var1="",var2=""){
  #devtools::install_github("haleyjeppson/ggmosaic")
  if(var2!=FALSE){
  library(ggmosaic)
  library(R.utils)
  meta.data$f1<-meta.data[[var1]]
  meta.data$f2<-meta.data[[var2]]
  p<-ggplot(data = meta.data)+
    geom_mosaic(aes(x = ggmosaic::product(f1,f2),fill= f1),na.rm = TRUE)+theme_mosaic()+xlab(label =var2)+ylab(label = var1)+labs(fill=var1)

    }
}

mosaic_scdrs_2(meta.data = meta.data,var1 = celltype,var2 = split.by)

summary_tab<-function(meta,celltype,trait="",threshold=0.1){
  library(janitor)
  if(str_length(trait)>0){
    itrait<-trait
    trait<-paste0("_",trait)
  }
  meta$FDR_thresh <-if_else(meta[[paste0("fdr",trait)]]<threshold, true = "Significant","No significant")
  meta_ed<-meta%>%dplyr::filter(FDR_thresh%in%"Significant")
  tab<-meta_ed%>%tabyl(all_of(celltype))%>%arrange(desc(n))%>%adorn_totals()
  print(paste("Trait:",itrait))
  print(tab)
  fwrite(x = tab,file = fs::path(out,paste0("tab_", celltype ,"_",trait,"FDR_",threshold,".tsv")),sep = "\t")
  }

summary_tab2vars<-function(meta,var1,var2,trait="",threshold=0.2){
  if(var2!=FALSE){
  library(janitor)
  if(str_length(trait)>0){
    itrait<-trait
    trait<-paste0("_",trait)
  }
  meta$FDR_thresh <-if_else(meta[[paste0("fdr",trait)]]<threshold, true = "Significant","No significant")
  meta_ed<-meta%>%dplyr::filter(FDR_thresh%in%"Significant")
  tab<- as.data.frame.matrix(table(meta_ed[[var1]],meta_ed[[var2]]))%>%
    rownames_to_column("ID")%>%adorn_totals(where = c("row","col"))%>%arrange(desc(Total))
  print(tab)
  fwrite(x = tab,file = fs::path(out,paste0("tab_", celltype,"_" ,var2,trait,"FDR_",threshold,".tsv")),sep = "\t")
  }
}

for(trait in traits){
  summary_tab(meta.data,trait=trait,celltype = celltype)
  summary_tab(meta.data,trait=trait,celltype = celltype,threshold=0.2)
  summary_tab(meta.data,trait=trait,celltype = celltype,threshold=0.4)
  summary_tab(meta.data,trait=trait,celltype = celltype,threshold=0.6)

  summary_tab2vars(meta = meta.data,trait=trait,var1 = celltype,var2=split.by)
}

summary_tabprop<-function(meta,trait,var1,threshold=0.1){
  
    if(str_length(trait)>0){
      itrait<-trait
      trait<-paste0("_",trait)
    }
    meta$FDR_thresh<-if_else(meta[[paste0("fdr",trait)]]<threshold, true = "Significant","No significant")
    tab<-as.data.frame.matrix(prop.table(table(paste0(meta[[var1]]),meta$FDR_thresh),1)*100)%>%
      rownames_to_column("ID")%>%adorn_totals(where = c("col"))%>%arrange(desc(`No significant`))
    print(paste("Trait:",itrait))
    print(tab)
    fwrite(x = tab,file = fs::path(out,paste0("proptab_", celltype ,"_",trait,"FDR_",threshold,".tsv")),sep = "\t")

}

summary_tabprop2<-function(meta,trait,var1,var2,threshold=0.1){
  if(var2!=FALSE){
  if(str_length(trait)>0){
    itrait<-trait
    trait<-paste0("_",trait)
  }
  meta$FDR_thresh<-if_else(meta[[paste0("fdr",trait)]]<threshold, true = "Significant","No significant")
  tab<-as.data.frame.matrix(prop.table(table(paste0(meta[[var1]]),meta$FDR_thresh),1)*100)%>%
    rownames_to_column("ID")%>%adorn_totals(where = c("col"))%>%arrange(desc(Total))
  print(paste("Trait:",itrait))
  print(tab)
  fwrite(x = tab,file = fs::path(out,paste0("proptab_", celltype ,"_",trait,"FDR_",threshold,".tsv")),sep = "\t")
  }
}


for(trait in traits){
summary_tabprop(meta.data,trait = trait,var1 = celltype,threshold=0.1)
summary_tabprop(meta.data,trait = trait,var1 = celltype,threshold=0.2)
summary_tabprop(meta.data,trait = trait,var1 = celltype,threshold=0.3)
summary_tabprop(meta.data,trait = trait,var1 = celltype,threshold=0.4)
summary_tabprop(meta.data,trait = trait,var1 = celltype,threshold=0.6)

}

plot_bar <- function(metadata, category,trait,point_size=0.20,feature, min_cell=15,filter_by = "fdr", threshold = 0.20, outlier_threshold = Inf){
v<-metadata[[paste0(filter_by,"_", feature)]] < threshold

if(sum(v)>=1){
x<-table(metadata[[paste0(category)]],v)
order<-data.frame(CellType=rownames(x),NonSigcell=x[,1],Sigcell=x[,2])%>%arrange(Sigcell)%>%pull(CellType)
keep<- data.frame(CellType=rownames(x),NonSigcell=x[,1],Sigcell=x[,2])%>%filter(Sigcell>=min_cell)%>%pull(CellType)
if(length(keep)>=1){
head(category)

metadata[[paste0(category)]] = factor(metadata[[paste0(category)]], levels=c(order))
#metadata[[paste0("category")]] = factor(metadata[[paste0("category")]], levels=c(legend_order))

#metadata%>%ggplot(aes(x= get(paste0(filter_by,"_" ,feature)),y = Integrated_05))+  geom_bar(stat = "identity", fill = "blue") + theme_classic()
index<-metadata[[category]]%in%keep

#x = get(paste0(filter_by,"_", feature)),
p<-metadata %>% 
  dplyr::filter(get(paste0(filter_by,"_", feature)) < threshold, index)%>%
  ggplot(aes(y = get(celltype),fill=  get(celltype))) +
  geom_bar(stat = "count") +
  labs(x = "Number of disease relevant cells", y = "") +
  theme_classic() + 
  theme(text = element_text(size=20),  plot.title = element_text((hjust = 0)),legend.position ="none" )+ # adjust margin here+
# scale_fill_manual(values = single_cell_palette,name = "Cell types") +
   guides( color=guide_legend(override.aes = list(size=8)))+
  ggtitle(paste0(""))
  #ggtitle(paste0("Significant cells by ", filter_by, " threshold in ", trait))

ggsave(p,filename= fs::path(out , paste0("barplot_",feature,".png")), width=9)
}
}
}

lapply(traits,function(x){
trait<-x#%>%str_replace("ankylosing", "AS")%>%str_replace_all("_"," ")%>%str_replace("immunochip.+","ImmunoChiP")
feature<- x
p<-plot_bar(metadata= meta.data, category = celltype, trait=trait,feature=x, threshold = 0.20)
})

