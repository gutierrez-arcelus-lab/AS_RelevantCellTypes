library(readxl)
library(biomaRt)
library(tidyverse)
library(data.table)
# Specify the file path
file_path <- "/lab-share/IM-Gutierrez-e2/Public/Marcos/Projects/AS/AS_Li2017extendedBiomart_gwas_forChinasMS.xlsx"

# Read the first sheet from the Excel file
data <- read_excel(file_path)

annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens", 
              "gencode.v41lift37.annotation.gtf") %>% 
    read_tsv(comment = "#", col_types = "c-cii-c-c",
             col_names = c("chr", "feature", "start", "end", "strand", "info"))

bed <- annotations %>%
    filter(feature == "gene") %>%
    mutate(tss = case_when(strand == "+" ~ start, 
			   strand == "-" ~ end, 
			   TRUE ~ NA_integer_)) %>%
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = sub("^(ENSG\\d+)\\.\\d+((_PAR_Y)?)$", "\\1\\2", gene_id)) %>%
    dplyr::select(chr, start, end, gene_id, gene_name)
extract_genes_in_window <- function(chromosome, position, window = 250000, bed) {
  # Define the start and end positions with the window
  start_pos <- position - window
  end_pos <- position + window
  
  # Filter the genes in the specified region
  genes_in_region <- bed %>%
    filter(chr == chromosome & start <= end_pos & end >= start_pos)
  
  return(genes_in_region)
}

data_gtf<-lapply(1:nrow(data), function(x){
chromosome <-  data[x,] %>%pull("CHR") 
snp_position <- data[x,] %>%pull("POS")
SNP_ID <- data[x,] %>%pull ("SNP_ID")
genes_in_window <- extract_genes_in_window(chromosome, snp_position, bed = bed)
})

SNPseagenes<- data_gtf %>% rbindlist()
fwrite(SNPseagenes , file = "02_SNPsea_GSE124731/AS_Li2017_genesnearby.txt")

SNPseagenes<- fread("02_SNPsea_GSE124731/AS_Li2017_genesnearby.txt")

dir<-"02_SNPsea_GSE124731/"
args<-list()

args$oDir<-fs::path(dir, "pseudo_bulk_LMM")
args$oFile<-"NK_vsNK_others_all"

# LMM NK versus T cells 
filename<-paste0(args$oDir,"/",args$oFile,"_mixedmodel_filtered.txt")
NKvsTcells_bulkRNAseq<- fread(filename)
"TAGAP"%in%genes
"IL23R"%in%genes
NKvsTcells_bulkRNAseq %>%filter(GENE_NAME %in%"IL23R")
model_results %>%filter(GENE_NAME %in%"TAGAP")
model_results %>%filter(GENE_NAME %in%"TAGAP")

NKvsTcells_bulkRNAseq %>%filter(GENE_NAME %in%genes)%>%filter(fdr< 0.05) %>%arrange(fdr)
fwrite( NKvsTcells_bulkRNAseq %>%filter(GENE_NAME %in%genes)%>%filter(fdr<=0.05)%>%arrange(fdr)%>%dplyr::select(-gene, -GENE_TYPE)%>%dplyr::select(GENE_ID,GENE_NAME,everything()),paste0("AS_Figures/supplementary_NKvsTcells.tsv"),sep="\t")

NKvsTcells_bulkRNAseq_ASriskloci <-NKvsTcells_bulkRNAseq%>%filter(GENE_NAME %in% SNPseagenes$gene_name)
sig_genes1 <- NKvsTcells_bulkRNAseq_ASriskloci %>%filter(abs(x1) > 1 , fdr < 0.05) 
sig_genes1 <- NKvsTcells_bulkRNAseq_ASriskloci %>%filter(abs(x1) > 1 , fdr < 0.05) %>%arrange(fdr)

magma_genes<-read.table(file = "scdrs_analysis/gs_file_gutcellatlas/AS_ImmunoChip.gs",header = TRUE)
genes<-magma_genes$GENESET%>%str_split(",")%>%unlist()%>%str_remove(":.+")%>%head(1000)
score<-magma_genes$GENESET%>%str_split(",")%>%unlist()%>%str_remove(".+:")
magma_df<-data.frame(GENE= genes,SCORE = score)%>%rownames_to_column("rank") 

NKvsImmuneGutCells_scRNAseq<-fread("pseudobulk/pseudo_bulk_LMM/NK_vs_Others_mixedmodel_filtered.txt")
NKvsImmuneGutCells_scRNAseq_ASriskloci <- NKvsImmuneGutCells_scRNAseq %>% filter(gene %in% genes)
sig_genes2<-NKvsImmuneGutCells_scRNAseq_ASriskloci %>% filter(gene %in% genes) %>%filter(fdr < 0.05 , x1 >1)


library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(Seurat)
library(SeuratDisk)
library(viridis)
library(data.table)

# x1 to
merged_df <- merge(NKvsTcells_bulkRNAseq_ASriskloci, NKvsAllGutCells_scRNAseq_ASriskloci, by.x = "GENE_NAME", by.y = "gene", suffixes = c("_NKvsTcells_bulkRNA-seq", "_NKvsImmuneGutCells_scRNA-seq"), all = FALSE)
merged_df_sig <- merged_df %>%filter(`x1_NKvsTcells_bulkRNA-seq` > 1, `x1_NKvsImmuneGutCells_scRNA-seq` > 1,
                   `fdr_NKvsTcells_bulkRNA-seq` < 0.05,  `fdr_NKvsImmuneGutCells_scRNA-seq` < 0.05)
colnames(merged_df_sig)
merged_df_sig <- merged_df_sig %>%dplyr::select(Gene = GENE_NAME , `beta_NKvsTcells_bulkRNA-seq` = `x1_NKvsTcells_bulkRNA-seq` ,`Std_NKvsTcells_bulkRNA-seq`, `p_value_NKvsTcells_bulkRNA-seq`,   
 `fdr_NKvsTcells_bulkRNA-seq`, `beta_NKvsImmuneGutCells_scRNA-seq` = `x1_NKvsImmuneGutCells_scRNA-seq` , `Std_NKvsImmuneGutCells_scRNA-seq` ,`p_value_NKvsImmuneGutCells_scRNA-seq` ,`fdr_NKvsImmuneGutCells_scRNA-seq`   
 )
fwrite(merged_df_sig , "AS_Figures/Table_S2_overlaploci.tsv", sep = "\t")
fwrite(merged_df_sig %>%arrange( desc(`beta_NKvsTcells_bulkRNA-seq`)) , "AS_Figures/supp/Table_S2_overlaploci.csv")
#merged_df_sig %>%arrange( desc(`beta_NKvsTcells_bulkRNA-seq`))%>%pull(Gene)  
# "FCGR3A" "SLAMF7" "TBX21"  "APOBR"  "NOTCH1" "RUNX3"  "IL18R1" "GPR65" 
# Replace NA values with 0 (assuming that missing means not significant in the other dataset)
merged_df_sig[is.na(merged_df_sig)] <- 0

# Prepare data for the heatmap
heatmap_data <- as.matrix(merged_df_sig %>% dplyr::select( `beta_NKvsTcells_bulkRNA-seq`, `beta_NKvsImmuneGutCells_scRNA-seq` ))
rownames(heatmap_data) <- merged_df_sig$Gene
library(pheatmap)
library(viridis)
# Create heatmap
png(file = "pheatmap genes.png")
# Set values lower than 1 to NA (so they will be white in the heatmap)
breaks <- c(-Inf, 1, 2, Inf)  # Define the boundaries for color blocks
colors <- c("white", "#FDBE85", "#FD8D3C", "#E6550D")  # White, Light Orange, Orange, Dark Orange

# Create heatmap with discrete color blocks
pheatmap(heatmap_data, cluster_rows = FALSE, cluster_cols = FALSE, # Disable clustering
         display_numbers = TRUE, color = colors, breaks = breaks, na_col = "white",
         main = "Differential Expression Comparison (Discrete Color Blocks)",
         labels_col = c("beta_NKvsTcells_bulkRNA-seq", ""))
dev.off()



GO<-enrichGO(gene = genes_intersect,keyType = "SYMBOL",ont = "BP",OrgDb = org.Hs.eg.db)
entrezid_data<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=as.character(genes_intersect), columns=c("ENTREZID","SYMBOL", "GENENAME"), keytype="SYMBOL")
entrezid<-AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=as.character(genes_intersect), columns=c("ENTREZID","SYMBOL", "GENENAME"), keytype="SYMBOL")%>%select(ENTREZID)%>%na.omit()%>%pull(ENTREZID)
as.data.frame(entrezid) %>% fwrite(file = "entrezid_NK.txt")

kegg<- enrichKEGG(gene = as.numeric(entrezid),organism = "hsa",pvalueCutoff = 1)

kegg<- enrichKEGG(gene = as.numeric(entrezid),organism = "hsa",keyType  = "ncbi-geneid",pvalueCutoff = 1)
library("pathview")
getIDfrompathway<-function(enrichR,IDs){
    enrichR@result%>%filter(ID%in%IDs)%>%pull(geneID)%>%str_split("/")
}
kegg@result%>%filter(p.adjust<0.05)


library(ggvenn)

# Example gene sets
gene_list <- list(
  Set1 = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE"),
  Set2 = c("GeneB", "GeneC", "GeneF", "GeneG", "GeneH"),
  Set3 = c("GeneC", "GeneI", "GeneJ", "GeneK", "GeneL")
)
# Create the Venn diagram
p<- ggvenn(gene_list, 
       fill_color = c("#FF9999", "#99CC99", "#9999FF"),
       stroke_size = 0.5, 
       set_name_size = 4, 
       text_size = 4 ,   show_elements = TRUE)
ggsave(p , file = "gene_intersect.png")
install.packages("ggVennDiagram")

library(ggVennDiagram)
# Create the Venn diagram showing gene names
# Create the Venn diagram
venn <- ggVennDiagram(gene_list, label = "none") +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# Extract the region labels to display genes
region_labels <- get_venn_region_data(venn)

# Create a custom Venn plot with gene names as labels
venn_plot <- ggVennDiagram(gene_list, label = "none") +
  geom_text(data = region_labels, aes(x = x, y = y, label = name), size = 3) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

# Display the Venn diagram
ggsave(venn_plot , file = "gene_intersect.png")
