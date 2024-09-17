#Purpose: Differential selection through Linear Mixed Models 
#Author: Daniela Fernandez Salinas
#Version: 7
#Date: October 17 2022
#Last modification: 
#Usage Examples:
#Dependencies: R:
	#lme4,dyplr,optparse
#Developer Notes:

# 	Packages
require(lme4)
require(dplyr)
require(optparse)
require(tidyverse)

# 	Functions
read_time<-function(name,file){
    t0<-Sys.time()
    temp<-read.table(file, sep = '\t', header = TRUE, row.names=1)
    tf<-Sys.time()-t0
    banner<-paste(file, "...stored in: ",name,". Reading time: ", round(tf,3),"minutes.", sep="")
    print(banner)
    assign(name,temp, envir=parent.frame())
    rm(temp)   
}

selection_model<-function(name,covType,expressed,covF=NULL,covR=NULL){
#Defining regression variables (predictor and covariates)
	#Predictor variable
	x<-as.factor(meta_matrix[,name])
	#Notation for covariates with random effect
	if(covType=="1" || covType=="0"){
	covRandom <- paste(paste("(1|meta_matrix$",covR,")",sep=""),collapse=" + ")
	}
	if(covType=="1"){
	covs <- paste(" ", covRandom, sep="")
	}
	#Notation for covariates with fixed effect
	if(covType=="-1" || covType=="0"){
    covs<- paste(paste("meta_matrix$", covF, sep=""), collapse=" + ")
  	}

  	#Using covariates with multiple effects
  	if(covType=="0"){
  	covs <- paste(covs, covRandom, sep=" + ")
  	}


#Perform test and control regression per peak
	#y<- expressed[1,]
 	results <- apply(expressed,1,function(y){
 		#Differential Model
 	  	form1 <- as.formula(paste("y ~ x +", covs))
 	  	#Control Model
 	  	form0 <- as.formula(paste("y ~", covs))
 	  	#print(form1)
 	  	if(covType=="-1"){
 	  		lm1 <- lm(form1)
 			lm0 <- lm(form0)
 	  	}else{
 	  	lm1 <- lmer(form1, REML=FALSE)
 			lm0 <- lmer(form0, REML=FALSE)
 	  	}
 		#Perform likelihood ratio test using anova function
 		anv <- anova(lm0, lm1)
 		#Output Stats
 		c(summary(lm1)$coefficients[2,], anv$P[2])
 	})           
 	results <- t(results)
 	colnames(results)<-c("x1","Std","x1_t","p_value")

	if(nrow(results)!=nrow(expressed)){
	 		print("WARNING: Input and output lines differ!")
	}

return(results) 	
}


################### 	MAIN 	########################
#Parameter Definition
parser=list(
	make_option(c("--data","-d"),help="Main input containing counts matrix.",metavar="FILE"),
	make_option("--simulation",action="store_true",help="Use if you don't intend to run the entire data set. (Only first 100 lines will be used)",default=FALSE),
	make_option("--oDir",help="Directory where to write the output.",default=".",metavar="/PATH/", type="character"),
	make_option("--oFile",help="Name for the output file.",default="Differential-Output", metavar="FILENAME", type="character"),
	make_option(c("--name","-n"),help="Use to run single analysis mode. Indicate the name of the (column) in your metamatrix for which to select.",type="character",metavar="STRING"),
	make_option(c("--meta","-m"),help="Binary meta_matrix.",type="character",metavar="FILE"),
	make_option(c("--random","-r"),help="Add column names used for covariates with random effect."),
	make_option(c("--fixed","-f"),help="Add column names used for covariates with random effect."),
	make_option(c("--all","-a"),action="store_true",help="Use to test all rows in the counts matrix without filtering."),
	make_option("--number",help="Use to test all instances that have a value greater than 1 in at least N samples."),
	make_option("--sorting_column",help="Provide the name of the column in your metamatrix that matches the column names of the counts matrix for sorting.", default="Library", metavar="STRING"),
	make_option(c("--columns","-c"),type="character",help="Use to indicate columns to keep from your counts matrix.", metavar="FILE",default=NULL)
)
args=parse_args(OptionParser(option_list=parser))
library(fs)
args<-list()
library(data.table)
fwrite(cluster_counts ,file= "cluster_counts_gutimmune.tsv", sep = "\t")
counts<-cluster_counts
args$oDir<-fs::path(pseudo_dir,"pseudo_bulk_LMM")
args$oFile<-"NK_vs_GutImmuneCells"
meta_matrix<-cluster_metadata
meta_matrix<-meta_matrix%>%mutate(contrast=ifelse(group_id=="NK cells",1,0))
head(meta_matrix)
meta_matrix$group_id%>%table()
args$name<-"contrast"
args$r<-"Sample.name" 
args$f<-NA #add sex as fixed effect
args$number<- 1
args$sorting_column<-"Subject"
  
  
#Required Parameters Check
if(is.character(args$data)==FALSE){
	print("ERROR: NO counts matrix provided.")
	stop()
}
if(is.character(args$meta)==FALSE){
	print(args$meta)
	print("ERROR: No meta matrix provided.")
	stop()
}
if(is.character(args$name)==FALSE){
	print("ERROR: No column indicated for selection (name).")
	stop()
}

#Printing Parameters
print("-----------DIFFERENTIAL SELECTION------------")
print(paste0("Counts Matrix:",args$data))
print(paste0("Writting output to:",args$oDir,"/",args$oFile))
print(paste0("File to meta_matrix: ",args$m))
print(paste0("Sorting meta_matrix by: ",args$sorting_column))
print(paste0("Performing selection for: ",args$name))

if(args$simulation){
	print("Running simulation MODE")
}

fixed<-NULL
random<-NULL
covars<-0
if(is.character(args$r)){
	covars<-covars+1
	random<-args$r
}

if(is.character(args$f)){
	covars<-covars-1
	if(grepl(",",args$f)){
		fixed<-str_split(args$f,",",simplify=T)
		#print(paste(paste("meta_matrix$", fixed, sep=""), collapse=" + "))
	}else{
		fixed<-args$f
	}
}


#Reading Input
print("Reading input...")
#meta_matrix<-read.table(args$meta,header=T,sep='\t')
print("Meta matrix successfully read!")
print("Reading counts...")
#read_time("logcounts",args$data)

#Subseting Counts Matrix
if(!is.null(args$columns)){
  print(paste("Subsetting columns from",args$columns,sep=" "))
  selected_columns<-read.table(args$columns)
  logcounts <- logcounts %>% select(selected_columns$V1)
  #meta_matrix<-meta_matrix[which(meta_matrix$sorting_column %in% selected_columns$V1),]
}

#Sorting Metamatrix

if(ncol(logcounts)!=sum(colnames(logcounts) %in% meta_matrix[,args$sorting_column])){
	print("ERROR: Counts column names don't match the meta_matrix.")
	print("Levels in sorting_column:\n")
	print(meta_matrix[,args$sorting_column])
	print("Columns in counts:")
	print(colnames(logcounts))
	stop()
}else{
	row.names(meta_matrix)<-meta_matrix[,args$sorting_column]
	meta_matrix<-meta_matrix[colnames(logcounts),]
	#print(meta_matrix)
}

#Simulation Mode
if(args$simulation){
	logcounts<-logcounts[1:100,]
	print("Subsetting first 100 peaks...")
}

#Defining Output Path
output_path<-args$oDir
if(!dir.exists(output_path)){
 dir.create(output_path)
 print(paste("Created output directory:",output_path,sep=""))
}

if(isTRUE(args$a)){
	log_filtered<-logcounts
	print(paste0("Testing all instances in the counts matrix ",nrow(log_filtered)))
}else{

	minimum<-as.integer(args$number)
	print(paste0("Filtering instances with a value greater than 1 in at least ",minimum," columns..."))
	log_filtered<- logcounts[which(rowSums(logcounts>1)>=minimum),]
	print(paste0("Testing ",nrow(log_filtered)," rows."))
	
}
library(edgeR)
### Normalization ####
#y <- DGEList(cluster_counts, samples=cluster_metadata)
#/ calculate TMM normalization factors:
#y <- calcNormFactors(y)
#/ get the normalized counts:
cluster_metadata$Sample.name%>%unique()%>%length()
####5 counts in at least 20 samples 
metadata$Sample.name%>%unique()%>%length() #41 individual

ncol(cluster_counts)
keep<-rowSums(cluster_counts>=5)>=20
sum(keep) #13667 keep genes
cluster_counts_ed<-cluster_counts[keep,]

log2cpm1 <- cpm(cluster_counts_ed+1, log=TRUE)
head(log2cpm1)
png(fs::path(pseudo_dir,"log2cpm_dist.png"))
p<- plot(density(log2cpm1)) #log2(cpm+1) 
dev.off()

#Running Regression 
head(log2cpm1)
suppressMessages({
model_results<-selection_model(name=args$name,expressed=log2cpm1,covR=random,covF=fixed,covType=as.character(covars))
})
log2cpm1 %>%head()
model_results<-model_results%>%as.data.frame()%>%rownames_to_column("gene")%>%arrange(p_value)
library(qvalue)
pvalname<-paste0("p_value")
fdrname <-paste0("fdr")
fdr<-qvalue(p=model_results[[pvalname]])
model_results$fdr<-fdr$qvalues
#Writing output file
filename_regression<-paste(output_path,"/",args$oFile,"_mixedmodel_filtered_onlyimmune.txt", sep="")
print(filename_regression)
fwrite(log2cpm1%>%as.data.frame()%>%rownames_to_column("ID"),file = fs::path(output_path,"log2cpm1.tsv"), sep = "\t")

print("Done!")
Results_lmm<-model_results
Results_lmm%>%as.data.frame%>%filter(fdr<0.05,x1 > 2.5)
Results_lmm%>%head

qqx<-c()
qqx[1]<-min(Results_lmm$x1)-1 
qqx[2]<-max(Results_lmm$x1)+1 
qqy<- (-log10(Results_lmm$p_value))
qqy<-c(0,  qqy[1]+0.5)
subname<-"NK"
write.table(model_results, fs::path(filename_regression), row.names=F,col.names=T,quote=F, sep = "\t")


library(fs)
Name <- "NK" 
project <- pseudo_dir
get_boxplot_lmm(lognormcounts = log2cpm1%>%as.data.frame,DEA = Results_lmm,contrast = Name, genes = Results_lmm%>%as.data.frame%>%filter(fdr<0.001,x1 > 2.5)%>%pull(gene)%>%head(100),save = TRUE)
log2cpm1%>%as.data.frame%>%head()


get_boxplot_lmm<-function(lognormcounts,genes,DEA,contrast="",subdir="",rows_to_colum=TRUE,save=FALSE){
  dir_create(fs::path(project,"Boxplots"))
  #Assumes that ids match
  #Look for all the genes selected
  if(rows_to_colum){
   lognormcounts <- lognormcounts%>%rownames_to_column(var = "ID")
   DEA<- DEA%>%rownames_to_column(var = "ID")
  }
  filtcounts<-lognormcounts%>%filter(ID%in%genes)
  if(nrow(filtcounts)==0){return(NA)}
  tcounts <- t(filtcounts%>%column_to_rownames(var="ID"))
  
  for(i in filtcounts$ID){
  print(i)
  #Filter DEA by gene 
  gene_DE<-DEA%>%filter(ID%in%i)
  exp<-gene_DE$x1
  P<-gene_DE$fdr
  gene_counts <- tcounts %>% merge(meta_matrix, ., by="row.names")
  gene_counts<-gene_counts%>%mutate(gene=gene_counts[,i])
  gene_counts%>%head()
  gene_counts$condition<- gene_counts$group_id%>%str_replace("_"," ")
  gene_counts$condition<-as.factor(gene_counts$condition)
  gene_counts$condition <- relevel( gene_counts$condition, ref = "NK cell")
  contrast<- contrast%>%str_replace_all("_"," ")
  p<-gene_counts%>%ggplot(aes(x=condition,y=gene,color=condition))+geom_boxplot(outlier.shape = NA)+ geom_jitter(color="black", size=0.4, alpha=0.9)+theme_classic()+xlab(sprintf("P = %s, logFC = %s",signif(as.numeric(P), 2),
        signif(as.numeric(exp), 2)))+ylab("log2 +1 DESeq2 norm counts")+theme(text = element_text(size = 22),plot.title = element_text(hjust=0.5,face = "italic"),axis.text.x = element_text(angle = 70, hjust =1))+labs(title=contrast,subtitle = i)
  # Get counts 
    #plot(p)
    if(save){ 
    dir_create(fs::path(project,"Boxplots",subdir))
    dir_create(fs::path(project,"Boxplots",subdir,Name))
    ggsave(filename = fs::path(project,"Boxplots",subdir,Name, paste0(i,".png")),plot = p,width = 10,height = 8)  
    } 
      
  }
  
}
