# Raphael Mourad
# 15/06/2023


# This R script is used to do enrichment analyses for dsbSNPs at non B DNA motifs and DRIP peaks (R loops)

# SETUP DIRECTORY
setwd("/media/sauber/Elements/DSBsnp_project_seb")


# LOAD LIBRARIES
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(ggplot2)
library(gwascat)
library(rtracklayer)
library(data.table)
library(ChIPseeker)
library(ChIPpeakAnno)
#library(bnlearn)
library(Rgraphviz)
library(liftOver)
library(GWASTools)
library(sigminer)
#library(seqplots)
library(LDlinkR)
library(readxl)
library(clusterProfiler)


# LOAD INHOUSE LIBRARIES
source("scriptR/comp_NonBDNA_SNP.R")
source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")


# SETUP HG19 GENOME
Chr.V=c(paste0("chr",1:22),"chrX")
seqinfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)[Chr.V]


# PARAMETER
mod="all" # "all" "RAD51" "BRCA2" "53BP1"
 
###
# LOAD SNPS

load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_GR.RData")
SNP_DSB.GR

load(file="results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR


###
# GENE SET ENRICHMENT ANALYSIS (NO ENRICHMENT FOUND)

dataGenes=read.table("data/genes/ucsc_genes_ensembl_hg19.csv",header=F,sep='\t')
genes.GR=GRanges(dataGenes[,2],IRanges(dataGenes[,4],dataGenes[,5]),ENSG=dataGenes[,9])
SNP_DSB.GR2=SNP_DSB.GR
start(SNP_DSB.GR2)=start(SNP_DSB.GR2)-5e3
end(SNP_DSB.GR2)=end(SNP_DSB.GR2)+5e3
olGenes=findOverlaps(SNP_DSB.GR2,genes.GR)
genes_DSB.GR=unique(genes.GR[subjectHits(olGenes)])

ids_genes_DSB <- bitr(genes_DSB.GR$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego <- enrichGO(gene=ids_genes_DSB$ENTREZID, OrgDb='org.Hs.eg.db', ont="BP")
ego
barplot(ego)



###
# Non B DNA

kCtrlSNPs=length(SNP_ctrl.GR)/length(SNP_DSB.GR) #

# All non B DNA
bed_nonBDNA=list.files("data/nonBDNA/")
matEnrichNonBDNA=matrix(NA,length(bed_nonBDNA),4)
for(i in 1:length(bed_nonBDNA)){
 nonBDNAi.GR=import.bed(paste0("data/nonBDNA/",bed_nonBDNA[i]))
 start(nonBDNAi.GR)=start(nonBDNAi.GR)-1
 nonBDNAi.GR=GRanges(seqnames(nonBDNAi.GR),ranges(nonBDNAi.GR))
 olnonBDNAi=findOverlaps(SNP_DSB.GR,nonBDNAi.GR)
 SNP_DSB_nonBDNAi.GR=SNP_DSB.GR[queryHits(olnonBDNAi)]
 SNP_DSB_nonBDNAi.GR=unique(SNP_DSB_nonBDNAi.GR)

 olnonBDNAi_ctrl=findOverlaps(SNP_ctrl.GR,nonBDNAi.GR)
 SNP_ctrl_nonBDNAi.GR=SNP_ctrl.GR[queryHits(olnonBDNAi_ctrl)]
 SNP_ctrl_nonBDNAi.GR=unique(SNP_ctrl_nonBDNAi.GR)
 
 tabi=rbind(c(length(SNP_DSB_nonBDNAi.GR),length(SNP_DSB.GR)),c(length(SNP_ctrl_nonBDNAi.GR)/kCtrlSNPs,length(SNP_DSB.GR)))
 fti=fisher.test(tabi)
 
 matEnrichNonBDNA[i,1]=length(SNP_DSB_nonBDNAi.GR)
 matEnrichNonBDNA[i,2]=length(SNP_ctrl_nonBDNAi.GR)/kCtrlSNPs
 matEnrichNonBDNA[i,3]=length(SNP_DSB_nonBDNAi.GR)/matEnrichNonBDNA[i,2]
 matEnrichNonBDNA[i,4]=fti$p.value
  
 print(i)
}
rownames(matEnrichNonBDNA)=c("AR","DR","G4","IR","MR","STR","Z-DNA")
colnames(matEnrichNonBDNA)=c("DSB_SNP","CTRL_SNP","FC","p")
matEnrichNonBDNA=data.frame(matEnrichNonBDNA)
write.table(matEnrichNonBDNA,file="results/allelic_imbalance/nonBDNA/matEnrich_nonBDNA.csv",row.names=T,col.names=T,sep='\t',quote=F)

matEnrichNonBDNAPlot=data.table(NonBDNA=rownames(matEnrichNonBDNA),NbSNP=c(matEnrichNonBDNA[,1],matEnrichNonBDNA[,2]),type=c(rep("DSB_SNP",7),rep("Ctrl_SNP",7)))
matEnrichNonBDNAPlot$type=factor(matEnrichNonBDNAPlot$type,levels=c("DSB_SNP","Ctrl_SNP"))
######################## get rloops from bW drip
sequ_list=c("DRIP")
w_list=c(100)

for(sequ_name in sequ_list){
  print(sequ_name)
  f<-list.files(paste0("dsbSNP_Av_pr/",sequ_name),pattern ="narrowPeak")
  
  
  sequ<-read_narrowpeaks(paste0("dsbSNP_Av_pr/",sequ_name,"/",f))
  
  
  
  for (w in w_list)  {
    dsb<-SNP_DSB.GR%>%resize(fix = "center",width = w+1)

    
    ctrl<-SNP_ctrl.GR%>%resize(fix = "center",width = w+1)
    
    dsb_ov<-dsb[queryHits(findOverlaps(dsb,sequ))]
    
    ctrl_ov<-ctrl[queryHits(findOverlaps(ctrl,sequ))]
    
    #ratio_eSNP<-mean(dsb_ov$score)/mean(ctrl_ov$score)
    
    
    # dsb_no_eSNP_ov<-sequ[subjectHits(findOverlaps(dsb_no_eSNP,sequ))]
    # 
    # ctrl_no_eSNP_ov<-sequ[subjectHits(findOverlaps(ctrl_no_eSNP,sequ))]
    # 
    # ratio_no_eSNP<-mean(dsb_no_eSNP_ov$score)/mean(ctrl_no_eSNP_ov$score)
    
    
    df<-data.frame(sequ_name,w)
    #df$ratio_no_eSNP<-ratio_no_eSNP
    df$Ctrl_SNP<-length(ctrl_ov)/10
    
    df$DSB_SNP<-length(dsb_ov)
    #f$dsb_no_eSNP_ov<-mean(dsb_no_eSNP_ov$score)
    df$w=w
    if(sequ_name==sequ_list[1]&w==w_list[1]){
      f.df=df
      
    }else{
      f.df=rbind(f.df,df)
    }
  }
}
f.df<-f.df%>%select(-w)%>%rename(NonBDNA=sequ_name)%>%mutate(NonBDNA="Rloop")%>% pivot_longer(c(Ctrl_SNP,DSB_SNP),values_to ="NbSNP",names_to = "type" )
matEnrichNonBDNAPlot<-rbind(matEnrichNonBDNAPlot,f.df)
##################### print plots

file_enrich_nonBDNA="results/allelic_imbalance/nonBDNA/barplot_enrich_nonBDNA.pdf"
pdf(file_enrich_nonBDNA,6,3)
ggplot(data=matEnrichNonBDNAPlot, aes(x=NonBDNA, y=NbSNP, fill=type)) +
geom_bar(stat="identity", position=position_dodge())
dev.off()

file_enrich_FC_nonBDNA="results/allelic_imbalance/nonBDNA/barplot_enrich_FC_nonBDNA.pdf"
pdf(file_enrich_FC_nonBDNA,6,4)
barplot(matEnrichNonBDNA$FC,names.arg=rownames(matEnrichNonBDNA),
	ylab="Fold-change (dsbSNP/ctrlSNP)",ylim=c(0,4))
dev.off()
