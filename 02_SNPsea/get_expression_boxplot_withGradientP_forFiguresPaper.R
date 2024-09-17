library(mhg)
library(liger)
library(org.Hs.eg.db)
library(RUVSeq)
library(RColorBrewer)
library(qvalue)
library(gplots)
library(ggplot2)
#wd <- "/Users/mgutierr/Documents/work/ITC/rnaseq/kallisto_lincRNAs/validations/"
#wd <- "/Users/mgutierr/Documents/work/ITC/rnaseq/kallisto_lincRNAs/validations/"

load("./ITC_rnaseq_qc_norm_coreCellT.rda")
ls()

match <- read.table(file = "./gene_ID_name_match_gencodeV24.txt", header = T, stringsAsFactors = F )
ids <- sapply(match$GENE_ID, function(gene){
  strsplit(gene, "\\.")[[1]][1]
})
rownames(match) <- ids

PLOT_GENE <- function(gene){
  reorder <- c(which(m$cell_type %in% "CD4"), which(m$cell_type %in% "CD8"), which(m$cell_type %in% "MAIT"), which(m$cell_type %in% "NKT"), which(m$cell_type %in% "Vd1"), which(m$cell_type %in% "Vd2"), which(m$cell_type %in% "NK"))
  m <- m[reorder,]
  tpm <- tpm[,reorder]
  log2tpm <- log2tpm[,reorder]
  cbPalette <- c("#D7A4EA", "#D7A4EA", "#D7A4EA", "#D7A4EA","#D7A4EA", "#D7A4EA",  "#FF8020")
  geneID <- rownames(match)[which(match$GENE_NAME %in% gene)]
   x <- as.numeric(log2tpm[geneID,])
  #pval <- format.pval(grad$Pvalue[which(grad$GENE_NAME %in% gene)])
  #beta <- round(grad$Beta[which(grad$GENE_NAME %in% gene)], digits = 2)
  png(filename = paste0(gene,".png"),res = 600,units = "in",height = 6,width = 8 )
  boxplot(list(CD4=x[which(m$cell_type %in% "CD4")], CD8=x[which(m$cell_type %in% "CD8")], MAIT=x[which(m$cell_type %in% "MAIT")], NKT=x[which(m$cell_type %in% "NKT")], Vd1=x[which(m$cell_type %in% "Vd1")], Vd2=x[which(m$cell_type %in% "Vd2")], NK=x[which(m$cell_type %in% "NK")]), col = cbPalette, ylab = "log2(tpm+1)", cex.names = 1.5, cex.lab = 1.5, cex.axis = 1.5, main = paste(gene, sep = ""), cex.main = 1.5)
  dev.off()
  boxplot(list(CD4=x[which(m$cell_type %in% "CD4")], CD8=x[which(m$cell_type %in% "CD8")], MAIT=x[which(m$cell_type %in% "MAIT")], NKT=x[which(m$cell_type %in% "NKT")], Vd1=x[which(m$cell_type %in% "Vd1")], Vd2=x[which(m$cell_type %in% "Vd2")], NK=x[which(m$cell_type %in% "NK")]), col = cbPalette, ylab = "log2(tpm+1)", cex.names = 1.5, cex.lab = 1.5, cex.axis = 1.5, main = paste(gene, sep = ""), cex.main = 1.5)
  
}
PLOT_GENE("RUNX3")
PLOT_GENE("NPEPPS")
PLOT_GENE("TNFRSF1A")
PLOT_GENE("BTN3A2")
PLOT_GENE("FCGR2A")


grad <- read.table("../gradient/geneNames_results_gradient.txt", header = T, stringsAsFactors = F)

pdf("boxplot_RPS26_RPL21.pdf", width = 7, height = 4)
PLOT_GENE("RPS26")
PLOT_GENE("RPL21")
dev.off()


PLOT_GENE("NKG7")
PLOT_GENE("CD4")
PLOT_GENE("CD8A")
PLOT_GENE("CD8B")

PLOT_GENE("CCR7")

PLOT_GENE("IL7")
PLOT_GENE("IL15RA")
PLOT_GENE("IL2")
PLOT_GENE("IL2R")


PLOT_GENE("PRDM1")
PLOT_GENE("RORC")
PLOT_GENE("MAF")
PLOT_GENE("BCL6")
PLOT_GENE("GZMK")
PLOT_GENE("GZMA")
PLOT_GENE("GZMM")
PLOT_GENE("GZMB")
PLOT_GENE("PRF1")
PLOT_GENE("ITGAM")
PLOT_GENE("ITGAM")
PLOT_GENE("TBXAS1")
PLOT_GENE("LGALS1")
PLOT_GENE("CD52")
PLOT_GENE("CLU")
PLOT_GENE("JAKMIP1")
PLOT_GENE("HLA-DRB1")
PLOT_GENE("CIITA")
PLOT_GENE("CIITA")


pdf("boxplot_lincRNA_TRgammaDeltaGenes.pdf", width = 7, height = 4)
PLOT_GENE("LINC00299")
PLOT_GENE("ID2")
PLOT_GENE("ZBTB16")
PLOT_GENE("TRGC2")
PLOT_GENE("TRGC1")
PLOT_GENE("TRGV2")
PLOT_GENE("TRGV3")
PLOT_GENE("TRDC")
PLOT_GENE("TRDV1")
PLOT_GENE("TRDV2")
PLOT_GENE("TRDC")
dev.off()


match[grep("TRD", match$GENE_NAME),]

PLOT_GENE("BCL6")
PLOT_GENE("ID2")
PLOT_GENE("TCF3")
PLOT_GENE("TCF7")
PLOT_GENE("ID3")
PLOT_GENE("CXCR5")
PLOT_GENE("IL6R")
match[grep("IL6", match$GENE_NAME),]

PLOT_GENE("STAT5B")
PLOT_GENE("STAT5A")

pdf("boxplot_LINC00299_ID2_PLZF.pdf", width = 7, height = 4)
PLOT_GENE("LINC00299")
PLOT_GENE("ID2")
PLOT_GENE("ZBTB16")
dev.off()
PLOT_GENE("TBX21")
PLOT_GENE("TCF")
PLOT_GENE("TBX21")
PLOT_GENE("IL6ST")
PLOT_GENE("ZEB2")
PLOT_GENE("EOMES")

head(match)
match[grep("LINC00299", match$GENE_NAME),]
PLOT_GENE("ORC4")
PLOT_GENE("ORC2")
PLOT_GENE("ORC1")
PLOT_GENE("CDC6")

PLOT_GENE("RORC")
PLOT_GENE("MTOR")

dev.off()
PLOT_GENE("TRDC")
PLOT_GENE("TRDJ1")
PLOT_GENE("TRDJ1")
PLOT_GENE("TRDV2")
PLOT_GENE("TRGJP")
PLOT_GENE("TRGC1")
PLOT_GENE("TRGC2")
PLOT_GENE("TRGV9")
PLOT_GENE("TRGV9")

PLOT_GENE("SLAMF6")
PLOT_GENE("CIITA")
PLOT_GENE("HLA-DQB1")
PLOT_GENE("HLA-DRB1")
PLOT_GENE("HLA-DQA1")

x11(width = 12.5, height = 8)
pdf("boxplots_TCF7_LEF1.pdf", width = 12.5, height = 8)
par(mar = c(3,5,5,2), mfrow = c(2,2) )
PLOT_GENE("TCF7")
PLOT_GENE("LEF1")
dev.off()

x11(width = 12.5, height = 4)
pdf("boxplots_ID2_BHLHE40.pdf", width = 7.5, height = 8)
par(mar = c(3,5,5,2), mfrow = c(2,1) )
PLOT_GENE("ID2")
PLOT_GENE("BHLHE40")
dev.off()


PLOT_GENE("IL7R")

PLOT_GENE("RP11-340E6.1")


x11(width = 12.5, height = 8)
pdf("boxplots_PLZFposITCupregGenesInPaper.pdf", width = 12.5, height = 8)
par(mar = c(3,5,5,2), mfrow = c(2,2) )
PLOT_GENE("ID2")
PLOT_GENE("BHLHE40")
PLOT_GENE("LTK")
PLOT_GENE("ME1")
dev.off()
#
PLOT_GENE("MSC")
PLOT_GENE("FOXP3")

match[grep("RP11-340E", match$GENE_NAME),]

PLOT_GENE("TIGIT")
PLOT_GENE("CD226")
PLOT_GENE("PVRL2")
PLOT_GENE("PVR")
PLOT_GENE("HAVCR2")
PLOT_GENE("BTLA")
PLOT_GENE("CTLA4")
PLOT_GENE("CD8A")
PLOT_GENE("CD8B")
PLOT_GENE("PDCD1")
PLOT_GENE("CD86")
PLOT_GENE("TRAIL")

PLOT_GENE("EOMES")
PLOT_GENE("TBX21")
PLOT_GENE("IFNG")
PLOT_GENE("GZMB")
PLOT_GENE("PRF1")

PLOT_GENE("LEF1")
PLOT_GENE("TCF7")

match[grep("PVR", match$GENE_NAME),]

PLOT_GENE("PRDM1")
PLOT_GENE("ZNF683")
PLOT_GENE("BHLHE40")

PLOT_GENE("LTK")

PLOT_GENE("IRF4")
PLOT_GENE("BACH2")


PLOT_GENE("CCR1")
PLOT_GENE("CCR2")
PLOT_GENE("CCR5")
PLOT_GENE("CXCR4")

PLOT_GENE("CXCR1")
PLOT_GENE("CXCR2")

PLOT_GENE("CCR3")

PLOT_GENE("EOMES")
PLOT_GENE("TBX21")
PLOT_GENE("ZNF683")
PLOT_GENE("PRDM1")
PLOT_GENE("IRF4")
PLOT_GENE("ZNF365")


PLOT_GENE("TRGC1")
PLOT_GENE("TRGC2")
PLOT_GENE("IFNG")
dev.print("boxplot_TRGC2_ITCs.pdf", dev = pdf)
match[grep("PRG4", match$GENE_NAME),]
PLOT_GENE("PRG4")
PLOT_GENE("COL1A2")
dev.print("boxplot_COL1A2_ITCs.pdf", dev = pdf)

PLOT_GENE("RPL22")

PLOT_GENE("CCL4")
PLOT_GENE("XCL2")
PLOT_GENE("CCL3")
PLOT_GENE("CCL5")
PLOT_GENE("XCL1")
PLOT_GENE("TNFSF14")
PLOT_GENE("PRF1")
PLOT_GENE("GZMB")
PLOT_GENE("TYROBP")
PLOT_GENE("KLRK1")

x11(width = 7.5, height = 5)
pdf("boxplots_TCRgenes_examples_adaptiveness.pdf", width = 7.5, height = 5)
par(mar = c(3,5,5,2))
PLOT_GENE("TRAV25")
PLOT_GENE("TRBV23-1")
dev.off()


x11(width = 7.5, height = 5)
pdf("boxplots_topTFs_innateness.pdf", width = 7.5, height = 5)
par(mar = c(3,5,5,2))
PLOT_GENE("TBX21")
PLOT_GENE("ZEB2")
PLOT_GENE("HOPX")
dev.off()

pdf("boxplots_topLincRNAs_innateness.pdf", width = 7.5, height = 5)
par(mar = c(3,5,5,2))
PLOT_GENE("RP11-121A8.1")
PLOT_GENE("RP11-305L7.3")
dev.off()


pdf("boxplots_topLincRNAs_innateness_possibleTargets.pdf", width = 7.5, height = 5)
par(mar = c(3,5,5,2))
PLOT_GENE("TRGC2")
PLOT_GENE("TRGV8")
PLOT_GENE("TRGV9")
PLOT_GENE("NFIL3")
dev.off()
PLOT_GENE("IFNG")
gos.signif

###

x11(width = 7.5, height = 5)
par(mar = c(3,5,5,2))
PLOT_GENE("G6PD")
dev.print("boxplot_G6PD.pdf", dev = pdf)
  
## HLA genes

grad.hla <- grad[grep("HLA-", grad$GENE_NAME),]
grad.hla.filt <- grad.hla[which(grad.hla$GENE_TYPE %in% "protein_coding"),]
dim(grad.hla.filt)
bonf <- 0.05/nrow(grad)
length(which(grad.hla.filt$Pvalue < bonf))# 11, of which 6 are class I and 5 are class II

grad.hla.filt[which(grad.hla.filt$Pvalue < bonf),]

pdf("boxplots_HLAgenes.pdf", width = 17, height = 8.5)
par(mfrow = c(2,3))
for(gene in grad.hla.filt$GENE_NAME){
  PLOT_GENE(gene)
}
dev.off()

PLOT_GENE("HLA-DRB1")
#PLOT_GENE("RPL22")
PLOT_GENE("HLA-DRB5")
PLOT_GENE("HLA-DRB6")
PLOT_GENE("HLA-DQA2")
PLOT_GENE("HLA-DMA")

PLOT_GENE("HLA-DQB1")

PLOT_GENE("KLRB1")
PLOT_GENE("IL12RB1")
PLOT_GENE("IL12RB2")

match[grep("")]

#PLOT_GENE("CD80")
PLOT_GENE("CD86")
#adhesion receptors (help MHC? for antigen presentation?)
PLOT_GENE("ITGAL")#CD11a
PLOT_GENE("ITGB2")#CD18
PLOT_GENE("ICAM1")#CD54

PLOT_GENE("HLA-A")
PLOT_GENE("HLA-B")
PLOT_GENE("HLA-C")
PLOT_GENE("HLA-E")

PLOT_GENE("HLA-DRA")
PLOT_GENE("HLA-DOA")
PLOT_GENE("HLA-DPA1")
PLOT_GENE("HLA-DPB1")
PLOT_GENE("HLA-DPB2")
PLOT_GENE("HLA-DMB")
PLOT_GENE("HLA-DMA")
PLOT_GENE("HLA-DOB")

PLOT_GENE("TYROBP")
PLOT_GENE("FGR")
PLOT_GENE("CLIC3")
PLOT_GENE("CCL4")
PLOT_GENE("TBX21")
PLOT_GENE("SH2D1B")

PLOT_GENE("LEF1")
PLOT_GENE("MAL")
PLOT_GENE("TIGIT")
PLOT_GENE("KLRB1")#CD161
PLOT_GENE("BTN3A1")


match[grep("CD1C", match$GENE_NAME),]
PLOT_GENE("CD1D")
PLOT_GENE("CD1A")
PLOT_GENE("CD1B")
PLOT_GENE("CD1C")
PLOT_GENE("MR1")
PLOT_GENE("CD4")
PLOT_GENE("CD8A")
PLOT_GENE("CD8B")

x11()
PLOT_GENE("LINC00299")
PLOT_GENE("ID2")
PLOT_GENE("MBOAT2")
PLOT_GENE("KIDINS220")
PLOT_GENE("RSAD2")
PLOT_GENE("ASAP2")
PLOT_GENE("CMPK2")

pdf("boxplots_TCRgenes_negBetaPassBonfGradient.pdf", width = 17, height = 8.5)
par(mfrow = c(2,3))
for(gene in neg.tr$GENE_NAME){
    PLOT_GENE(gene)
}
dev.off()
wd

PLOT_GENE("TRDJ1")
PLOT_GENE("TRDC")
PLOT_GENE("KLRB1")
PLOT_GENE("IL18RAP")

PLOT_GENE("CCR7")
PLOT_GENE("TRBV7-9")
PLOT_GENE("TRAV12-3")

PLOT_GENE("ZEB2-AS1")

PLOT_GENE("TRBV23-1")
PLOT_GENE("TRBV30")

PLOT_GENE("TRGV8")

PLOT_GENE("TRAV2")
PLOT_GENE("TRBV2")
PLOT_GENE("TRBV28")
PLOT_GENE("TRAV12-1")
PLOT_GENE("TRAV41")
PLOT_GENE("TRGC2")
PLOT_GENE("TRGC1")

PLOT_GENE("MAL")
PLOT_GENE("RP11-61O1.2")
PLOT_GENE("LRRN3")
PLOT_GENE("CCR7")
PLOT_GENE("TRAJ16")
PLOT_GENE("AQP3")
PLOT_GENE("FAM102A")
PLOT_GENE("CAMK4")

#lincrNAs innateness
PLOT_GENE("RP11-121A8.1")# in the middle of TCRgama locus!
PLOT_GENE("TRGV9")
PLOT_GENE("TRGV8")
PLOT_GENE("TRGC2")

PLOT_GENE("TRGV10")
PLOT_GENE("TRGV11")
PLOT_GENE("TRGV5")
PLOT_GENE("TRGV4")


PLOT_GENE("LINC00299")# more PLZF+ pattern

PLOT_GENE("RP11-305L7.3")
PLOT_GENE("AUH")
PLOT_GENE("SYK")
PLOT_GENE("NFIL3")

PLOT_GENE("ROR2")
PLOT_GENE("DIRAS2")

PLOT_GENE("BCL2")
PLOT_GENE("BAX")
PLOT_GENE("APAF1")
PLOT_GENE("CASP9")

PLOT_GENE("CYBA")
PLOT_GENE("RAB27A")
PLOT_GENE("ICAM1")
PLOT_GENE("GADD45A")
PLOT_GENE("DDAH2")
PLOT_GENE("GADD45A")
PLOT_GENE("GSTP1")
PLOT_GENE("ATP2B4")
PLOT_GENE("VIMP")
PLOT_GENE("GLA")
PLOT_GENE("IFNG")

#other clock genes
PLOT_GENE("PER1")
PLOT_GENE("NFIL3")

pdf("boxplots_ClockGenes_perIndivProfile.pdf", width = 17, height = 8.5)
par(mfrow = c(2,3))
PLOT_GENE("ARNTL")
PLOT_GENE("PER1")
PLOT_GENE("RORA")
PLOT_GENE("NR1D1")
PLOT_GENE("CRY1")
dev.off()

PLOT_GENE("CD27")
PLOT_GENE("CD38")


PLOT_GENE("IL18RAP")
PLOT_GENE("TBX21")

pdf("boxplots_MAITmarkers_reviewGapin2014.pdf", width = 14, height = 11)
par(mfrow = c(2,2))
#MAIT genes according to review
PLOT_GENE("TRAV1-2")
PLOT_GENE("TRAJ33")
PLOT_GENE("KLRB1")#CD161
PLOT_GENE("ABCD1")
PLOT_GENE("IL12RB1")
PLOT_GENE("IL12RB2")
PLOT_GENE("IL18R1")
PLOT_GENE("IL18RAP")
PLOT_GENE("IL23R")
PLOT_GENE("CCR6")#these are involved in traficking to peripheral tissues, particularly the instestine and the liver, so also NKT and Vd2 go there?
PLOT_GENE("CXCR6")
PLOT_GENE("ZBTB16")#PLZF TF
PLOT_GENE("RORC")#TF
#upone TCR stim, MAITS produce IFNG, TNF alpha and IL2, as well as IL17 after PMAI stim. 
#they are cytotoxic, killin bacteria-infected epithelial cells
#it is unclear if mouse MAIT cells are functionally and phenotypically equivalent to MAIT cells in humans
#in mice not so studied bc they are not so frequent, unless big infection in lungs? or other...
#later shown that only bacteria with riboflavin pathway stim MAIT (produce vitamin B2)
#MAIT recognized APCs infected with some bacteria and yeast but not virus.
PLOT_GENE("TRAJ12")# 5-15% of human MAIT cells were found to use two other TRAJ segments (TRAJ12 and TRAJ20) (ref 7)
PLOT_GENE("TRAJ20")
PLOT_GENE("TRBV20-1")#The TCRalpha chain is found paired with a limitd number of TCRBeta chains (TRBV20, TRBV6, TRBV2)
PLOT_GENE("TRBV6-1")#i am plotting all trbv6 genes bc unclear which one is the one referred in the review, but Ildiko's review says TRBV6-1
PLOT_GENE("TRBV6-4")
PLOT_GENE("TRBV6-5")
PLOT_GENE("TRBV6-6")
PLOT_GENE("TRBV6-7")
PLOT_GENE("TRBV6-8")
PLOT_GENE("CCL20")#MAIT rapidly secrete IFNG, IL17 and CCL20 upon Ag recognition.
PLOT_GENE("IFNG")
PLOT_GENE("TNF")
#PLOT_GENE("IL17")#is not considered expressed
#based on cytokine secretion profile, MAIT cells might be good targets for new vaccines under development (how?)
#many diseases mentioned in review (IBD, HIV infection, MS, etc)
dev.off()

PLOT_GENE("LL22NC03-75H12.2")
dev.print("boxplots_MAITspecificLincRNA_LL22NC03-75H12.2.pdf", dev = pdf)

PLOT_GENE("IRF4")
match[grep("TNFA", match$GENE_NAME),]

pdf("boxplots_log2tpm_ROSgenesMaketal2017.pdf", width = 14, height = 11)
par(mfrow = c(2,2))
PLOT_GENE("SLC3A2")#CD98
PLOT_GENE("GCLC")
PLOT_GENE("GCLM")
PLOT_GENE("IL2RA")
PLOT_GENE("IL2RB")
PLOT_GENE("IL2RG")
PLOT_GENE("MYC")
PLOT_GENE("ME1")
PLOT_GENE("G6PD")
dev.off()

PLOT_GENE("GLYCTK")
PLOT_GENE("PGM1")
PLOT_GENE("CYBA")
PLOT_GENE("G6PD")
dev.print("boxplots_CYBA_GSPD.pdf", dev = pdf)

PLOT_GENE("CYBA")

pdf("boxplots_lincRNAs.pdf", width = 14, height = 11)
par(mfcol = c(2,2))
gn <- "LINC00402"
PLOT_GENE(gn)
gn <- "LINC00299"
PLOT_GENE(gn)
dev.off()

pdf("boxplots_log2tpm_TFs_passBonfFCgradient.pdf", width = 14, height = 11)
par(mfrow = c(2,2))
PLOT_GENE("HOPX")
PLOT_GENE("ZEB2")
PLOT_GENE("MYC")
PLOT_GENE("LEF1")
dev.off()

pdf("boxplots_log2tpm_genesAssocZEB2_CTBP1_SMADs_TGFB.pdf",  width = 14, height = 11)
par(mfrow = c(2,2))
PLOT_GENE("TGFBR3")
PLOT_GENE("TGFBR1")
PLOT_GENE("TGFB1")
PLOT_GENE("TGFBRAP1")
PLOT_GENE("TGFB3")
PLOT_GENE("TGFBI")
PLOT_GENE("TGFBR2")
PLOT_GENE("CTBP1")
PLOT_GENE("SMAD5")
PLOT_GENE("SMAD7")
PLOT_GENE("ID2")
PLOT_GENE("ID3")
PLOT_GENE("MAPK1")
PLOT_GENE("ZFYVE16")
PLOT_GENE("CREBBP")
PLOT_GENE("RHOA")
PLOT_GENE("NOG")
dev.off()

PLOT_GENE("HSP90AA1")
dev.print("boxplot_HSP90AA1.pdf", dev = pdf)

pdf("boxplots_log2tpm_genesAssocLEF_TLE1_UBC_NLK.pdf",  width = 8, height = 7)
par(mfrow = c(2,2))
PLOT_GENE("TLE1")
PLOT_GENE("UBC")
PLOT_GENE("NLK")
PLOT_GENE("TLE4")
dev.off()

pdf("boxplots_log2tpm_genesAssocLEF_top15TRAgenesNegBeta.pdf",  width = 8, height = 7)
par(mfrow = c(2,2))
PLOT_GENE("TRAV8-3")
PLOT_GENE("TRAV12-1")
PLOT_GENE("TRAV13-1")
PLOT_GENE("TRAV21")
PLOT_GENE("TRAV12-3")
PLOT_GENE("TRAV9-2")
PLOT_GENE("TRAV4")
PLOT_GENE("TRAV16")
PLOT_GENE("TRAV1-2")
PLOT_GENE("TRAV8-6")
PLOT_GENE("TRAJ20")
PLOT_GENE("TRAV2")
PLOT_GENE("TRAV17")
PLOT_GENE("TRAV8-2")
PLOT_GENE("TRAJ12")
dev.off()

tr <- grad[grep("TRA", grad$GENE_NAME),]
tr[order(tr$Beta),]

pdf("boxplots_Tbet_helios.pdf", width = 6, height = 4.5)
PLOT_GENE("TBX21")
PLOT_GENE("TBX21")
PLOT_GENE("IKZF2")
dev.off()

PLOT_GENE("TNFAIP8")#negative mediator of apoptosis
PLOT_GENE("HDAC9")

pdf("boxplots_TbetInteractingGenes_IL12_JAK_STAT.pdf", width = 12, height = 9)
par(mfrow = c(2,2))
PLOT_GENE("IL12RB1")
PLOT_GENE("IL12RB2")
PLOT_GENE("IL12A")
PLOT_GENE("STAT1")
PLOT_GENE("JAK2")
PLOT_GENE("UBC")
PLOT_GENE("TYK2")
PLOT_GENE("GATA3")
PLOT_GENE("SMARCA2")
PLOT_GENE("SMARCA4")
dev.off()

PLOT_GENE("NKG7")
PLOT_GENE("MAL")

pdf("boxplots_log2tpm_MAIT_TCRgenes_IldikoRev.pdf", width = 5, height = 2.8)
PLOT_GENE("TRAV1-2")
PLOT_GENE("TRAJ33")
PLOT_GENE("TRBV20-1")
PLOT_GENE("TRBV6-1")
dev.off()

pdf("boxplots_log2tpm_NKT_TCRgenes_IldikoRev.pdf", width = 7, height = 5.5)
PLOT_GENE("TRAV10")
PLOT_GENE("TRAJ18")
PLOT_GENE("TRBV25-1")
dev.off()

match[grep("TRG", match$GENE_NAME),]
match[grep("TRD", match$GENE_NAME),]
pdf("boxplots_log2tpm_gammaDelta_TCRgenes.pdf", width = 7, height = 5.5)
PLOT_GENE("TRDV1")
PLOT_GENE("TRDV2")
PLOT_GENE("TRDJ1")
PLOT_GENE("TRDJ2")
PLOT_GENE("TRDJ3")
PLOT_GENE("TRDJ4")
PLOT_GENE("TRDC")
PLOT_GENE("TRDV3")

PLOT_GENE("TRGC2")
PLOT_GENE("TRGJ2")
PLOT_GENE("TRGJP2")
PLOT_GENE("TRGC1")
PLOT_GENE("TRGJP")
PLOT_GENE("TRGJP1")
PLOT_GENE("TRGV11")
PLOT_GENE("TRGV10")
PLOT_GENE("TRGV9")
PLOT_GENE("TRGV8")
PLOT_GENE("TRGV5")
PLOT_GENE("TRGV4")
PLOT_GENE("TRGV3")
PLOT_GENE("TRGV2")
PLOT_GENE("TRGV1")

dev.off()




pdf("boxplots_MYC_log2tpm.pdf", width = 7, height = 5.5)
PLOT_GENE("MYC")
dev.off()

pdf("boxplots_IFNG_log2tpm.pdf", width = 7, height = 5.5)
PLOT_GENE("IFNG")
dev.off()

PLOT_GENE("CCL5")
PLOT_GENE("CCL4")
PLOT_GENE("CCL3")
PLOT_GENE("CX3CR1")
PLOT_GENE("SPON2")
PLOT_GENE("MAPK1")
PLOT_GENE("CCR5")
PLOT_GENE("XCL2")
PLOT_GENE("XCL1")
PLOT_GENE("CCL3L3")
PLOT_GENE("VIMP")
PLOT_GENE("IRF8")
PLOT_GENE("IRAK1")
PLOT_GENE("FGR")
PLOT_GENE("CD58")
PLOT_GENE("CASP1")
PLOT_GENE("STX11")
PLOT_GENE("CD300A")
####
PLOT_GENE("PKC")


pdf("boxplots_riboProt_initFact_log2tpm.pdf", width = 7, height = 5.5)
PLOT_GENE("RPL36")
PLOT_GENE("RPL22")
PLOT_GENE("RPL36A")
PLOT_GENE("EIF3E")
PLOT_GENE("EIF3H")
PLOT_GENE("EIF3L")
dev.off()

PLOT_GENE("TLR3")
PLOT_GENE("TLR1")
PLOT_GENE("TRAF6")
PLOT_GENE("ECSIT")
PLOT_GENE("SDHA")
PLOT_GENE("SDHB")

PLOT_GENE("NSDHL")
PLOT_GENE("NLRP7")
PLOT_GENE("NLRC3")
PLOT_GENE("NLRP2")
PLOT_GENE("NLRX1")

pdf("boxplot_caspases.pdf", width = 9, height = 9)
par(mfrow = c(2,1))
PLOT_GENE("CASP6")
PLOT_GENE("CASP8")
PLOT_GENE("CASP4")
PLOT_GENE("CASP9")
PLOT_GENE("CASP3")
PLOT_GENE("CASP8AP2")
PLOT_GENE("CASP7")
PLOT_GENE("CASP1")
dev.off()

PLOT_GENE("TLR4")
PLOT_GENE("TLR10")
PLOT_GENE("TLR7")
PLOT_GENE("TLR5")

#signif?
PLOT_GENE("TRAF3IP3")
PLOT_GENE("TRAF1")

##
PLOT_GENE("PVRIG")#P=4.2e-10, polio vrius receptory, inhibits T cell proliferation?


x11(width = 9, height = 9)
par(mfrow = c(2,1))
PLOT_GENE("SIRT2")
dev.print("barplot_SIRT2.pdf")

PLOT_GENE("NOX1")

PLOT_GENE("EGR2")
tpm["ENSG00000122877",]
##### Pol I (transcribes rRNA)

gn <- "EGR"
gn <- "POL"
gn <- "POLR1"
gn <- "TAF1"
gn <- "RRN3"
gn <- "SIRT2"
gn <- "NOX"
match[grep(gn, match$GENE_NAME),]

gn <- 
gn <- "POLR1B"
x11(width = 7, height = 4)
PLOT_GENE("POLR1A")
PLOT_GENE("POLR1B")
PLOT_GENE("POLR1C")
PLOT_GENE("POLR1D")
PLOT_GENE("POLR1E")
PLOT_GENE("TAF1C")
PLOT_GENE("TAF1D")
PLOT_GENE("TAF1A")
PLOT_GENE("TAF1B")
PLOT_GENE("TAF1C")
PLOT_GENE("RRN3")
PLOT_GENE("MYC")
dev.print("barplot_MYC.pdf", dev = pdf)
PLOT_GENE("IFNG")
dev.print("barplot_IFNG.pdf", dev = pdf)
PLOT_GENE("TNF")
PLOT_GENE("TNFRSF1B")
PLOT_GENE("TNFSF14")


PLOT_GENE("IL2RG")
PLOT_GENE("IL2RB")
PLOT_GENE("IL2RA")
PLOT_GENE("IL7R")
PLOT_GENE("IL15RA")
PLOT_GENE("HLA-DRB1")

PLOT_GENE("PPARA")
PLOT_GENE("TP53")

pdf("barplots_tpm_houseKeepingGenesForRTpcr.pdf", width = 9, height = 8)
par(mfrow = c(2,1))
PLOT_GENE("GAPDH")
PLOT_GENE("RPL13A")
PLOT_GENE("ACTB")
PLOT_GENE("TUBA1A")#okish
PLOT_GENE("ALB")# the best but very lowly expressed
PLOT_GENE("TBP")#perfect, beta = -0.01, P = 0.62
PLOT_GENE("HPRT1")#ok but 3 outlier samples... beta = -0.005, P = 0.84
PLOT_GENE("POLR2H")#okish
PLOT_GENE("POLR2K")#innate signif
PLOT_GENE("POLR2G")#innate signif
dev.off()

pdf("barplots_tpm_STATgenes.pdf", width = 9, height = 9)
par(mfrow = c(2,1))
PLOT_GENE("STAT1")
PLOT_GENE("STAT4")
PLOT_GENE("STAT5B")
PLOT_GENE("STAT5A")
PLOT_GENE("STAT3")
PLOT_GENE("STAT2")
PLOT_GENE("STAT6")
dev.off()

x11(width = 9, height = 9)
par(mfrow = c(2,1))
PLOT_GENE("FOS")
PLOT_GENE("RAF1")
PLOT_GENE("CDKN1B")
PLOT_GENE("CDKN1B")
PLOT_GENE("IGFBP3")

pdf( "barplot_assocNegInnatenessGradient_pcg.pdf", width = 9, height = 9)
par(mfrow = c(2,1))
PLOT_GENE("TBC1D4")#induces glucose upake, within top hits neg assoc innateness grad, and also huge neg beta PLZF+ITCs vs adap
PLOT_GENE("LRRN3")#unkown function?
PLOT_GENE("FAM102A")# possibly estrogen receptor?
PLOT_GENE("MYC")
PLOT_GENE("NGFRAP1")#May be a signaling adapter molecule involved in p75NTR-mediated apoptosis induced by NGF
PLOT_GENE("RCAN3")#Inhibits calcineurin-dependent transcriptional responses by binding to the catalytic domain of calcineurin A. Calcineurin activates nuclear factor of activated T cell, cytoplasmic (NFATc), a transcription factor, by dephosphorylating it. The activated NFATc is then translocated into the nucleus, where it upregulates the expression of interleukin 2 (IL-2), which, in turn, stimulates the growth and differentiation of T cell response. Calcineurin is the target of a class of drugs called calcineurin inhibitors, which includes cyclosporin, voclosporin, pimecrolimus and tacrolimus.
PLOT_GENE("IL6ST")#Binding of IL6 to IL6R induces IL6ST homodimerization and formation of a high-affinity receptor complex, which activates Janus kinases (PubMed:2261637). That causes phosphorylation of IL6ST tyrosine residues which in turn activates STAT3
PLOT_GENE("CD5")#The encoded protein contains three SRCR domains and may act as a receptor to regulate T-cell proliferation
dev.off()

pdf( "barplot_MAPKs.pdf", width = 9, height = 9)
par(mfrow = c(2,1))
PLOT_GENE("MAPK14")#p38 subunit
PLOT_GENE("MAPKAPK5-AS1")
PLOT_GENE("MAPKAPK5")
PLOT_GENE("MAPK13")#p38 subunit! regulates p53, gradient beta = 0.106, P = 1e-06
PLOT_GENE("MAPK1")#ERK, upstream of MYC (activates it?), beta=0.3, P = 4.8e-20
dev.off()

PLOT_GENE("AQP3")#aquaporin, for water transport...? doesn't pass FC threshold
PLOT_GENE("ACTN1")#binds to actin, same story as TBC1D4

gn <- "PKR"
##### cytokines

gn <- "IL12A"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "IL12B"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "IL18"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "IL4"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "IFNG"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "TNF"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "IKZF2"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "TBX21"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)


gn <- "KIR2DL3"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)


### lincRNA

gn <- "LINC00402"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

gn <- "LINC00299"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

x11()
gn <- "ZBTB16"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)


x11(width = 9, height = 8)
par(mfrow = c(2,1))
PLOT_GENE(gn)
dev.print("barplot_LINC00402.pdf", dev = pdf)
dev.off()


gn <- "AC009495.2"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)

x11(width = 9, height = 8)
par(mfrow = c(2,1))
PLOT_GENE(gn)

gn <- "LINC00299"
match[grep(gn, match$GENE_NAME),]
PLOT_GENE(gn)
dev.print("barplot_2lincRNAsCorrelNegPC2.pdf", dev = pdf)
dev.off()

#CARKL
gn <- "SHPK"
match[grep(gn, match$GENE_NAME),]
x11(width = 9, height = 8)
par(mfrow = c(2,1))
PLOT_GENE(gn)
dev.print("barplot_SHPK_CARKLkinaseForPentPhP.pdf", dev = pdf)
dev.off()


#CARKL
gn <- "HUWE1"
match[grep(gn, match$GENE_NAME),]
x11(width = 9, height = 8)
par(mfrow = c(2,1))
PLOT_GENE(gn)
dev.print("barplot_HUWE1_ubuiquitinaseForExtraRibosomes.pdf", dev = pdf)
dev.off()

gn <- "POLR3"
match[grep(gn, match$GENE_NAME),]

# translation initiation factor

gn <- "EIF4G1"
match[grep(gn, match$GENE_NAME),]
x11(width = 9, height = 8)
par(mfrow = c(2,1))
PLOT_GENE(gn)
###################  OLD ############
## TFs from motif analysis
match[grep("NR5A2", match$GENE_NAME),]
match[grep("CTCFL", match$GENE_NAME),]
match[grep("NR5A2", match$GENE_NAME),]
pdf("barplot_expr_TFs.pdf", width = 12, height = 5)
#x11(width = 12, height = 5)
par(mar = c(4,5,5,1))
PLOT_GENE("NR5A2")

PLOT_GENE("BORIS")
PLOT_GENE("CTCFL")
PLOT_GENE("FOXL2")
  dev.off()



### genes for Yang
match[grep("CD1A", match$GENE_NAME),]
match[grep("CD1B", match$GENE_NAME),]
match[grep("CD1C", match$GENE_NAME),]
match[grep("CD1D", match$GENE_NAME),]

pdf("barplot_expr_CD1ABC_forYang.pdf", width = 12, height = 5)
#x11(width = 12, height = 5)
par(mar = c(4,5,5,1))
PLOT_GENE("CD1A")
PLOT_GENE("CD1B")
PLOT_GENE("CD1C")
dev.off()

pdf("barplot_expr_CD1D_ITC_forYang.pdf", width = 12, height = 5)
#x11(width = 12, height = 5)
par(mar = c(4,5,5,1))
PLOT_GENE("CD1D")
dev.off()

match[grep("RASA2", match$GENE_NAME),]
match[grep("RNF7", match$GENE_NAME),]
pdf("barplot_expr_RASA2_RNF7_forYang.pdf", width = 12, height = 5)
#x11(width = 12, height = 5)
par(mar = c(4,5,5,1))
PLOT_GENE("RASA2")
PLOT_GENE("RNF7")
dev.off()



head(tpm)
summary(as.numeric(tpm))
class(tpm)
summary(as.numeric(as.matrix(tpm)))

load("../mixedModels3/ITC_rnaseq_qc_norm3.rda")
summary(as.numeric(as.matrix(tpm)))
summary(as.numeric(as.matrix(log2tpm)))


match[grep("IL7", match$GENE_NAME),]

pdf("barplot_expr_STAT5_AB.pdf", width = 12, height = 5)
#x11(width = 12, height = 5)
par(mar = c(4,5,5,1))
PLOT_GENE("STAT5A")
PLOT_GENE("STAT5B")
dev.off()

match[grep("STARD3NL", match$GENE_NAME),]

pdf("barplot_expr_genes_TCRgammaLocus.pdf", width = 12, height = 5)
#x11(width = 12, height = 5)
par(mar = c(4,5,5,1))
PLOT_GENE("TRGC2")
PLOT_GENE("TRGJ2")
PLOT_GENE("TRGC1")
PLOT_GENE("TRGJP2")
PLOT_GENE("TRGJP")
PLOT_GENE("TRGJP")
PLOT_GENE("TRGV5")
PLOT_GENE("TRGV4")
PLOT_GENE("TRGV3")
PLOT_GENE("TRGV2")
PLOT_GENE("TRGV1")
PLOT_GENE("TRGV11")
PLOT_GENE("TRGV10")
PLOT_GENE("TRGV9")
PLOT_GENE("TRGV9")
PLOT_GENE("AMPH")
PLOT_GENE("STARD3NL")
PLOT_GENE("VPS41")
PLOT_GENE("EPDR1")
PLOT_GENE("NME8")
PLOT_GENE("SFRP4")
PLOT_GENE("GPR141")
PLOT_GENE("POU6F2")
PLOT_GENE("ELMO1")
PLOT_GENE("YAE1D1")
dev.off()

pdf("barplot_expr_genes_variousIL7.pdf", width = 12, height = 10)
#x11(width = 12, height = 5)
par(mar = c(4,5,5,1), mfrow = c(2,1))
PLOT_GENE("IFNG")
PLOT_GENE("TRGC2")
PLOT_GENE("IL7R")
PLOT_GENE("IL7")
PLOT_GENE("BACH2")
PLOT_GENE("ZBTB16")

dev.off()



##############

match[grep("TRG", match$GENE_NAME),]

write.table(cbind(rownames(match[grep("TRG", match$GENE_NAME),])), file = "TCRgamma_genes_IDs.txt", quote = F, row.names = F, col.names = F)
wd
