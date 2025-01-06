###############################################################################
# This R script is used to test the enrichment of dsbSNPs at eSNPs and other SNPs.
#############################################################################

# SETUP DIRECTORY
#setwd("/media/mourad/diskSave/MCF_Toulouse/recherche/DSBPsnp")
DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)

# LOAD LIBRARIES
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(stringr)


# LOAD INHOUSE LIBRARIES
source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")


# SETUP HG19 GENOME
Chr.V=c(paste0("chr",1:22),"chrX")
seqinfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)[Chr.V]



###
# LOAD SNPS

load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
SNP_DSB.GR
#write_bed(SNP_DSB.GR,'results/allelic_imbalance/SNP_DSB/SNP_DSB_deeptools.bed')

load(file="results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR
#write_bed(SNP_ctrl.GR,'results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb.bed')





###
# Compute the set of eQTL SNPs and save it as a GRanges object
# GTEx (tissues) / PANCANQTL (cancer)

eQTLproject="GTEx" #"pancanQTL"


if(F){ # RUN ONLY ONCE (ALREADY COMPUTED)
lf=list.files("data/GTEx_Analysis_v7_eQTL",pattern="pairs")
eQTL.GR=NULL
for(i in 1:length(lf)){
 dataEQTLi=as.data.frame(fread(paste0("data/GTEx_Analysis_v7_eQTL/",lf[i]),header=T,sep='\t'))
 phenoi=strsplit(lf[i],'[.]')[[1]][1]
 vid=t(sapply(as.character(dataEQTLi[,1]),function(x){as.character(strsplit(x,'_')[[1]][1:4])}))
 rownames(vid)=rep("",nrow(vid))
 eQTL.GRi=GRanges(paste0("chr",vid[,1]),IRanges(as.numeric(vid[,2]),as.numeric(vid[,2])),ref=vid[,3],alt=vid[,4],
	slope=dataEQTLi[,8],pheno=rep(phenoi,nrow(dataEQTLi)),tss_distance=dataEQTLi$tss_distance,
	gene=dataEQTLi$gene_id)
 if(i==1){
  eQTL.GR=eQTL.GRi
 }else{
  eQTL.GR=c(eQTL.GR,eQTL.GRi)
 }
 print(lf[i])
}
save(eQTL.GR,file="data/eQTL_GTEx.RData")
}
if(F){ # RUN ONLY ONCE (ALREADY COMPUTED)
dataPancan=as.data.frame(fread("data/PancanQTL/cis_eQTLs_all_re.gz",sep='\t',header=T))
vid=t(sapply(as.character(dataPancan[,5]),function(x){as.character(strsplit(x,'/')[[1]])}))
gid1=t(sapply(as.character(dataPancan[,7]),function(x){as.character(strsplit(x,':')[[1]][1:3])}))
dataPancan=dataPancan[dataPancan[,3]==gid1[,1] & nchar(as.character(gid1[,3]))==1,]
vid=t(sapply(as.character(dataPancan[,5]),function(x){as.character(strsplit(x,'/')[[1]])}))
gid1=t(sapply(as.character(dataPancan[,7]),function(x){as.character(strsplit(x,':')[[1]][1:3])}))
gid2=t(sapply(gid1[,2],function(x){as.character(strsplit(x,'-')[[1]])}))
strandGene=sapply(gid1[,3],function(x){as.character(strsplit(x,';')[[1]][1])})
gene.GR=GRanges(gid1[,1],IRanges(as.numeric(gid2[,1]),as.numeric(gid2[,2])),strand=strandGene)
eQTL.GR=GRanges(dataPancan[,3],IRanges(dataPancan[,4],dataPancan[,4]),rs=dataPancan[,2],
	ref=vid[,1],alt=vid[,2],slope=dataPancan$beta,gene=dataPancan$gene,pheno=factor(dataPancan[,1]))
dist=rep(NA,length(gene.GR))
dist[as.vector(strand(gene.GR)=='+')]=start(gene.GR[as.vector(strand(gene.GR)=='+')])-start(eQTL.GR[as.vector(strand(gene.GR)=='+')])
dist[as.vector(strand(gene.GR)=='-')]=start(eQTL.GR[as.vector(strand(gene.GR)=='-')])-end(gene.GR[as.vector(strand(gene.GR)=='-')])
eQTL.GR$tss_distance=dist
save(eQTL.GR,file="data/eQTL_pancanQTL.RData")
}


###
# Load eQTL and test enrichment

if(eQTLproject=="GTEx"){
 load(file="data/eQTL_GTEx.RData")
}else if(eQTLproject=="pancanQTL"){
 load(file="data/eQTL_pancanQTL.RData")
}

TSS=T
# Filter at TSS
if (TSS==T){
eQTL.GR=eQTL.GR[eQTL.GR$tss_distance>-2000 & eQTL.GR$tss_distance<0]
#eQTL.GR=eQTL.GR[abs(eQTL.GR$slope)>0.5]
}
# Overlap with eQTLs
olGTEx=findOverlaps(SNP_DSB.GR,eQTL.GR)
SNP_DSB_eQTL.GR=SNP_DSB.GR[queryHits(olGTEx)]
SNP_eQTL_DSB.GR=eQTL.GR[subjectHits(olGTEx)]
SNP_eQTL_DSB_uniq.GR=unique(SNP_eQTL_DSB.GR)

# Overlap of ctrlSNPs with eQTLs
olGTEx_ctrl=findOverlaps(SNP_ctrl.GR,eQTL.GR)
SNP_eQTL_ctrl.GR=eQTL.GR[subjectHits(olGTEx_ctrl)]
SNP_eQTL_ctrl_uniq.GR=unique(SNP_eQTL_ctrl.GR)

# Enrichment test
kCtrlSNPs=length(SNP_ctrl.GR)/length(SNP_DSB.GR)
tabSNPeQTL=rbind(c(length(SNP_eQTL_DSB_uniq.GR),length(SNP_DSB.GR)),c(length(SNP_eQTL_ctrl_uniq.GR)/kCtrlSNPs,length(SNP_DSB.GR)))
tabSNPeQTL

f=fisher.test(tabSNPeQTL)

FC=tabSNPeQTL[1,1]/tabSNPeQTL[2,1]
if (TSS==T){
file_barplot_enrichment_eSNPs=paste0("results/allelic_imbalance/eQTL/barplot_enrichment_eSNPs_Promoter_TSS_2kb",eQTLproject,".pdf")
}else{
  file_barplot_enrichment_eSNPs=paste0("results/allelic_imbalance/eQTL/barplot_enrichment_eSNPs_Promoter_All",eQTLproject,".pdf")
  
}


pdf(file_barplot_enrichment_eSNPs,3,4)
barplot(round(tabSNPeQTL[,1]),names.arg=c("dsbSNPs","ctrl SNPs"),ylab="Number of SNPs matching with eSNPs",
	col=c("red","blue"))
title(paste0("p=",format(f$p.value,digits=3)," FC=",  format(FC,digits=3)))

dev.off()


###
# Test enrichment by tissue

# Tissues
matSNP_eQTL=data.frame(pheno=names(table(eQTL.GR$pheno)),eQTl_DSB=as.numeric(table(SNP_eQTL_DSB.GR$pheno)),
	eQTL=as.numeric(table(eQTL.GR$pheno)))
matSNP_eQTL$perc=matSNP_eQTL$eQTl_DSB/matSNP_eQTL$eQTL
matSNP_eQTL=matSNP_eQTL[order(matSNP_eQTL$perc,decreasing=T),]
matSNP_eQTL=matSNP_eQTL[matSNP_eQTL$eQTl_DSB>=10,]

file_barplot_eSNP=paste0("results/allelic_imbalance/eQTL/barplot_eSNP_pheno_",eQTLproject,".pdf")
pdf(file_barplot_eSNP,9,6)
par(mar = c(bottom=5, left=15, top=2, right=15))
# barplot(matSNP_eQTL$perc[20:1]*100,names=matSNP_eQTL$pheno[20:1],las=2,horiz = TRUE,xlab="Overlap of dsbSNPs with eSNPs located at promoters (%)")
barplot(matSNP_eQTL$perc[10:1]*100,names=matSNP_eQTL$pheno[10:1],las=2,horiz = TRUE,xlab="Overlap of dsbSNPs with eSNPs located at promoters (%)")

dev.off()

print(matSNP_eQTL)


brain<-c()
for (i in(1:length(matSNP_eQTL$pheno))){
if(("Brain")%in%(str_split(matSNP_eQTL$pheno,"_")[[i]])|("brain")%in%(str_split(matSNP_eQTL$pheno,"_")[[i]])){
  brain[i]="brain"
}else{
  brain[i]="no_brain"
}
}
matSNP_eQTL$brain<-brain
matSNP_eQTL_brain<-matSNP_eQTL[matSNP_eQTL$brain=="brain",]
nrow(matSNP_eQTL_brain)

file_barplot_eSNP=paste0("results/allelic_imbalance/eQTL/barplot_eSNP_pheno_Brain_",eQTLproject,".pdf")
pdf(file_barplot_eSNP,9,6)
par(mar = c(bottom=5, left=15, top=2, right=15))
barplot(matSNP_eQTL_brain$perc[13:1]*100,names=matSNP_eQTL_brain$pheno[13:1],las=2,horiz = TRUE,xlab="Overlap of dsbSNPs with eSNPs located at promoters (%)")
dev.off()

matSNP_eQTL_nobrain<-matSNP_eQTL[matSNP_eQTL$brain=="no_brain",]
pheno="Brain_tissues_mean_percentage"
mean_brain<-as.data.frame(pheno)
mean_brain$eQTL_DSB<-"eqtl_DSB"

mean_brain$eQTL<-"eqtl"
mean_brain$perc<-mean(matSNP_eQTL_brain$perc)
mean_brain$brain<-"brain"
colnames(mean_brain)=colnames(matSNP_eQTL_nobrain)


matSNP_eQTL_nobrain<-rbind(mean_brain,matSNP_eQTL_nobrain)
matSNP_eQTL_nobrain<-matSNP_eQTL_nobrain[order(matSNP_eQTL_nobrain$perc,decreasing = T),]

file_barplot_eSNP=paste0("results/allelic_imbalance/eQTL/barplot_eSNP_pheno_noBrain_",eQTLproject,".pdf")
pdf(file_barplot_eSNP,9,6)
par(mar = c(bottom=5, left=15, top=2, right=15))

barplot(matSNP_eQTL_nobrain$perc[20:1]*100,names=matSNP_eQTL_nobrain$pheno[20:1],las=2,horiz = TRUE,xlab="Overlap of dsbSNPs with eSNPs located at promoters (%)")
dev.off()


