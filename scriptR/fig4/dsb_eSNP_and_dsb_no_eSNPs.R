#Sébastien AUBER
#13/10/2025

#################################################################################
# Analyze overlaps between DSB-associated SNPs and eSNPs:
# Annotate dsb eSNPs and dsb non eSNPs  
# Generate pie charts of genomic annotations and Venn diagrams of overlaps
################################################################################


DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)

# LOAD LIBRARIES
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(ChIPpeakAnno)
library(ChIPseeker)
library(readr)
library(tidyr)
library(ggVennDiagram)
library(dplyr)




# LOAD INHOUSE LIBRARIES
source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")


# SETUP HG19 GENOM

load(file="data/eQTL_GTEx.RData")

load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
SNP_DSB.GR

#FIND OVERLAP WITH eSNPs
Hits<- findOverlaps(SNP_DSB.GR, eQTL.GR)
dsb_eSNP<-unique(SNP_DSB.GR[queryHits(Hits)])
dsb_no_eSNP<-SNP_DSB.GR[!SNP_DSB.GR$ID %in% dsb_eSNP$ID]



#export
write_bed(dsb_eSNP,"results/allelic_imbalance/SNP_DSB/dsb_eSNP.bed")
write_bed(dsb_no_eSNP,"results/allelic_imbalance/SNP_DSB/dsb_no_eSNP.bed")




###########################Annotate#############################################
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
annot_dsb_eSNP=annotatePeak(dsb_eSNP,tssRegion=c(-3000, 3000), TxDb=txdb)
annot_dsb_no_eSNP=annotatePeak(dsb_no_eSNP,tssRegion=c(-3000, 3000), TxDb=txdb)

#plot
pdf("results/allelic_imbalance/SNP_DSB/pie_annot_dsb_eSNP.pdf")
plotAnnoPie(annot_dsb_eSNP)
dev.off()


pdf("results/allelic_imbalance/SNP_DSB/pie_annot_dsb_no_eSNP.pdf")

plotAnnoPie(annot_dsb_no_eSNP)
dev.off()


############################Venn Diagram########################################


#genrate common names based on position
dsb<-SNP_DSB.GR%>%as_tibble()%>%mutate(n=paste0(seqnames,"_",start))

eqtl<-eQTL.GR%>%as_tibble()%>%mutate(n=paste0(seqnames,"_",start))

#generate list
x <- list(
  dsb = dsb$n,
  eqtl = eqtl$n
)

#plot

pdf("results/allelic_imbalance/overlapps/dsbSNP_eSNP.pdf",3,3)
ggVennDiagram(x, label_alpha = 0)+  ggplot2::scale_fill_gradient(low="white",high = "red")+theme(legend.position = "none")

dev.off()


