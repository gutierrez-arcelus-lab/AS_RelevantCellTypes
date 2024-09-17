#Action: tell the reviewer that it was using both restins and stimulating, and repeat the differential accessibility analysis correcting for stimulation status, then re-run LDSC-seg. And add these results to the MS as supplementary material. 

#Purpose: Differential Selection of Open Chromatin Peaks through Lineal Mixed Model for all cell-states in Log2Normalized_RPKM_atac_pks_v2.txt
#Author: Daniela Fernandez Salinas
#Version: 3
#Date: October 25 2021
#Last modification: November 11, 2021
#Usage examples: 
	#Rscript Peak_Selection_full_dataset3.R --simulation --name "NKs" --compare "Monocytes" -o "Comparison"
	#Rscript Peak_Selection_full_dataset3.R --peaks 2nd_Selection/Selected_cellstates1.csv -o /home/ch226756/mylab/ATAC-seq/2nd_Selection/Cellstates1
#Dependencies: R:
	#lme4,ggplot2,gt
#Developer Notes: 	
#	(Modify) Three running modes: ALL (Generate entire meta matrix and run a model for each column.), SINGLE (Contrast only a single column against the rest of the dataset, like the first run with B Cells.), COMPARE (B Cells vs NKs mode)		
#	Indicate the name or the default output.
#Required input:
#	* Peaks count matrix
#	* List of groups to divide by.
#	* Meta matrix with predictor and covariates.
#	* Cutoff value(s)
##################################################
# 	Packages
require(lme4)
require(ggplot2)
require(gt)
require(dplyr)
require(optparse)
library(future.apply)
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

get_indexes<-function(linker,df){
	indexes<-lapply(linker, function(g){
	  wrapper<-sapply(g,function(n){
	    n2<-paste("\\.",n,sep = "")
	    grep(n2,colnames(df))
	  })
	  if(is.list(wrapper)){
	    wrapper<-unlist(wrapper)
	  }else{
	    wrapper<-as.vector(wrapper)
	  }
	  return(wrapper)
	})
	return(indexes)
}
#name<-"NKs"

selection_model<-function(name,write_pks=FALSE,log_atac_pks,atlas_indexes, peaks_name=NULL){

	#Selecting peaks that are greater than the mean for at least half the samples
	expressed_i<-which(apply(log_atac_pks[,atlas_indexes[[name]]],1,function(r){sum(r>atlas_summary[name,"Mean"])})>(length(atlas_indexes[[name]])/2)-1)
	expressed<-log_atac_pks[expressed_i,]
	print(paste("Subset peaks:",length(expressed_i),sep=" "))
	print(head(expressed))
	if (write_pks==TRUE){
		if(is.null(peaks_name)){
			filename_peaks<-paste(output_path,"/",name,"_peaks.txt",sep="")
		}else{
			filename_peaks<-paste(output_path,"/",peaks_name,sep="")
		}
		print("Writting down expressed peaks...")
		write.table(expressed, filename_peaks, sep='\t', quote=FALSE, row.names=TRUE,col.names=TRUE)
	}
   
#Defining regression variables (predictor and covariates)
	x<-as.factor(meta_matrix[,name])
    meta_m<-meta_matrix
    vx<-x
    #length(x)
	cov.rand<-c("Donor")
	covRandom <- paste(paste("(1|meta_matrix$",cov.rand,")",sep=""),collapse=" + ")
	
	covRandom <- paste(paste("(1|",cov.rand,")",sep=""),collapse=" + ")
	covs <- paste(covRandom, sep=" ")
    mycov<- covs
#Perform test and control regression per peak
    #y= expressed[,1]
    #length(x)) {y}, 
 	results <- future_apply(expressed,1,function(y, x= vx,meta_matrix = meta_m,covs = mycov){
 		#Differential Model
        form1 <- as.formula(paste("y ~ x + Stim +", covs))
 		lm1 <- lmer(form1, REML=FALSE,data= meta_matrix)
 		#Control Model
		form0 <- as.formula(paste("y ~ Stim +", covs))
 		lm0 <- lmer(form0, REML=FALSE,data= meta_matrix)
 		#Perform likelihood ratio test using anova function
 		anv <- anova(lm0, lm1)
 		#Output Stats
 		c(summary(lm1)$coefficients[2,], anv$P[2])
 	}) 

	#results1<-results
	#identical(results,results1)
 	results <- t(results)
 	colnames(results)<-c("x1","Std","x1_t","p_value")
 	filename_regression<-paste(output_path,"/",name,"_mixedmodel.txt", sep="")
 	write.table( results, filename_regression, row.names=T,col.names=T,quote=F, sep = "\t")
 
 	if(nrow(results)!=nrow(expressed)){
 		print("WARNING: Input and output lines differ!")
 	}
 	results<-as.data.frame(results)
 
 	###### 		Graphical Summary of Results 	######
 	halloween<-c("orange","black","purple") 	#Color palette
 	v_title<-paste("Differential Peaks in",name, sep=" ")
 	h_title<-paste("Histogram of",name,"regression p_values",sep=" ")
 
 	h_name<-paste(output_path, "/",name,"_hist.png",sep = "")
 	hist1<-ggplot(results, aes(x=p_value))+
   		geom_histogram(color="#e9ecef", alpha=0.6, position="identity", fill= "#69b3a2")+
   		theme_minimal()+
   		ggtitle(h_title)
   	ggsave(h_name,plot =hist1, dpi = "screen",height = 7,width = 10)
 
 	cutoff<-(0.05/nrow(results))
 	print(paste("Corrected cutoff for p-values:",cutoff, sep=" "))
 	results$Significant<-"No"
 	results$Significant[results$p_value < 0.05]<- "Yes"
 	results$Significant[results$p_value < cutoff]<-"Corrected"
 
 	v1_name<-paste(output_path, "/",name,"_v1.png",sep = "")
 	volcano1<-ggplot(results, aes(x=x1, y=-log10(p_value), col= Significant))+ 
 	  	geom_point()+
 	  	scale_color_manual(values = c("orange","black","darkred"))+
 	  	theme_minimal()+
 	  	ggtitle(v_title)+
 	  	xlab(expression(beta))+
 	  	geom_hline(yintercept = -log10(cutoff), color="red")+
 	  	geom_hline(yintercept= -log10(0.05), color= "blue")
 	ggsave(v1_name,plot =volcano1, dpi = "screen",height = 7,width = 10)
 
 	results$Differential<-"NO"
 	results$Differential[(results$p_value<cutoff & results$x1 > 2)]<-"UP"
 	results$Differential[(results$p_value<cutoff & results$x1 < -2)]<-"DOWN"
 
 	v2_name<-paste(output_path, "/",name,"_v2.png",sep = "")
 	volcano2<-ggplot(results, aes(x=x1, y=-log10(p_value), col= Differential))+ 
 	  geom_point()+
 	  scale_color_manual(values = halloween)+
 	  theme_minimal()+
 	  ggtitle(v_title)+
 	  xlab(expression(beta))+
 	  geom_hline(yintercept = -log10(cutoff), color="red")+
 	  geom_hline(yintercept= -log10(0.05), color= "blue")
 	ggsave(v2_name,plot =volcano2, dpi = "screen",height = 7,width = 10)
 	
 	#Results in numbers
 	total<-nrow(log_atac_pks)
 	peak_stats<-as.matrix(
 		c(Expressed=nrow(results),
 		Significant=sum(results$Significant!="No"),
 		Corrected=sum(results$Significant=="Corrected"),
 		Differential=sum(results$Differential!="NO"),
 		Up=sum(results$Differential=="UP"),
 		Down=sum(results$Differential=="DOWN"),
 		No=sum(results$Differential=="NO")))

 		peak_stats<-cbind(peak_stats,round((peak_stats[,1]/total)*100,1))
 		colnames(peak_stats)<-c("Peaks","Percentage")
 	stats_filename<-paste(output_path,"/",name,"_stats.txt",sep="")
 	write.table(peak_stats,stats_filename, sep='\t', quote=FALSE, row.names=TRUE,col.names=TRUE)
 	
 	#Nice table (PNG)
 	#peak_stats$Percentage<-paste(peak_stats$Percentage,"%",sep="")
	#peak_stats<-rbind(Total=c(829942,"100%"),stats)
	#gt_stats <- gt(stats, rownames_to_stub = T)
	#gt_name<-paste(output_path,"/",name,"_stats.png",sep="")
	#print(gt_name)
	#gtsave(gt_stats,gt_name)
 }

################### 	MAIN 	########################
#Parameter Definition
parser=list(
	make_option(c("--peaks","-p"),action="store",type="character",help="Use to indicate columns (cell-states) to use and generate a subset of the Log2Normalized_RPKM_atac_pks_v2.txt peaks.", metavar="FILE",default=NULL),
	make_option("--simulation",action="store_true",help="Use if you don't intend to run the entire data set. (Only first 100 lines will be used)",default=FALSE),
	make_option(c("--output","-o"),action="store",help="Directory where to write the output.",default="Differential-Output", metavar="NAME", type="character"),
	make_option("--name",action="store",help="Name of the celltype for which to run the analysis.",type="character",metavar="CELLTYPE"),
	make_option(c("--compare", "-c"),help="Use if you want to compare two groups in particular instead to controlling against the entire dataset."))
args=parse_args(OptionParser(option_list=parser))
plan(multisession, workers = 100) ## Run in parallel on local computer

#print(length(args))
#output_path<-"/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg/balanced_groups"
output_path=args$output
read_time("log_atac_pks","/lab-share/IM-Gutierrez-e2/Public/Dany/Calderon_ATAC/Counts-Matrix/Log2Normalized_RPKM_atac_pks_v2.txt")
read_time("atlas_summary", "/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg/Balanced_Log2Normalized_summary.txt")

if(args$simulation){
	log_atac_pks<-log_atac_pks[1:100,]
}

#args$peaks<- "/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/LMM_LDSC_seg/balanced_groups.csv" 
if(!is.null(args$peaks)){
  selected_cellstates<-read.csv(args$peaks)
  subset_indexes<-unlist(get_indexes(selected_cellstates,log_atac_pks))
  log_atac_pks<-log_atac_pks[subset_indexes]
}

#colnames(log_atac_pks[subset_indexes])
#Atlas definition (Dictionary)

atlas<-list(
  NKs=c("Mature_NK.S", "Mature_NK.U", "Memory_NK.U", "Immature_NK.U"),
  TCells=c("CD8pos_T.S","CD8pos_T.U","Naive_Teffs.U","Memory_Teffs.U"),
  BCells=c("Bulk_B.U","Mem_B.U","Naive_B.U","Bulk_B.S"),
  Monocytes=c("Monocytes.S","Monocytes.U"),
  DC1=("Myeloid_DCs.U"),
  Plasma=("Plasmablasts.U"),
  DC2=("pDCs.U"))

#Retrieving Indexes
atlas_indexes<-get_indexes(atlas,log_atac_pks)

#sapply(atlas_indexes, function(x){
#log_atac_pks[,x]%>%apply(1,mean)%>%summary
#})%>%t()%>%write.table(file = "LMM_LDSC_seg/Balanced_Log2Normalized_summary.txt",sep="\t",quote=FALSE)

#Creating output directory
if(!dir.exists(output_path)){
 dir.create(output_path)
}
#Subsetting peaks for comparative mode (One on one comparison, E.g: Monocytes vs NK)
args_name<-args$name
args_compare<-args$compare

if (is.character(args_compare)){
 	if (is.character(args_name)){
 		subset_i<-c(atlas_indexes[[args_name]],atlas_indexes[[args_compare]])
 		log_subset<-log_atac_pks[,subset_i]
 		new_indexes<-get_indexes(list(atlas[[args_name]],atlas[[args_compare]]),log_subset)
 		names(new_indexes)<-c(args_name,args_compare)
 		# Writting Meta-Matrix
		meta_matrix<-NULL
	  	binary<-rep(0,ncol(log_subset))
	  	binary[unlist(new_indexes[[args_name]])]<-1
	  	meta_matrix<-as.matrix(cbind(meta_matrix,binary))
	  	colnames(meta_matrix)<-args_name
		meta_matrix<-cbind.data.frame(meta_matrix,Donor=as.array(matrix(unlist(strsplit(colnames(log_subset),"\\.")), ncol = 3, byrow = T)[,1]))
		meta_matrix<-cbind.data.frame(meta_matrix,Library=colnames(log_subset))
		print(meta_matrix)
		selection_model(args_name,FALSE,log_subset,new_indexes)
 		}else{
	 		print("ERROR: NO main name provided.")
	 		quit()
 		}
}else{
	# Writting Meta-Matrix
	meta_matrix<-NULL
	for (c in 1:length(atlas_indexes)){
	  binary<-rep(0,ncol(log_atac_pks))
	  binary[unlist(atlas_indexes[[c]])]<-1
	  meta_matrix<-as.matrix(cbind(meta_matrix,binary))
	  colnames(meta_matrix)[c]<-names(atlas_indexes)[c]
	}
	meta_matrix<-cbind.data.frame(meta_matrix,Donor=as.array(matrix(unlist(strsplit(colnames(log_atac_pks),"\\.")), ncol = 3, byrow = T)[,1]))
   meta_matrix<-cbind.data.frame(meta_matrix,Stim=as.array(matrix(unlist(strsplit(colnames(log_atac_pks),"\\.")), ncol = 3, byrow = T)[,3]))

	meta_matrix<-cbind.data.frame(meta_matrix,Library=colnames(log_atac_pks))
	meta_name<-paste(output_path,"/Full_meta_matrix.txt",sep="")
	write.table(meta_matrix,meta_name, quote=F, col.names=T, row.names=F, sep='\t')


	#Peak Selection
	lapply(names(atlas),function(a){
		print(a)
		selection_model(name = a,TRUE,log_atac_pks,atlas_indexes)
	})
}	

