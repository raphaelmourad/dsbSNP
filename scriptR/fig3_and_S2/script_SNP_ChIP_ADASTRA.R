# Raphael Mourad
# 15/06/2023


# This R script is used to extract the dsbSNPs mapped from ADASTRA database.
# ADASTRA SNPs were originally mapped to hg38, and then we liftovered them to hg19.


# SET DIRECTORY


DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)
# LOAD R PACKAGES
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(rtracklayer)
library(data.table)
library(ChIPseeker)
library(ChIPpeakAnno)
library(liftOver)
library(fastDummies)

# LOAD INHOUSE R SCRIPTS
source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")




###
# PROCESS SNP FROM ADASTRA DATABASE

# Mapped on hg38
PvalT=0.1

protNames=c("ATM","BRCA1","BRCA2","NBN","RAD51") # dsbSNPs
for(i in 1:length(protNames)){
 file_ChIPi=paste0("data/ADASTRA_allelic_imbalance_database/adastra.bill_cipher/release_dump/TF/",protNames[i],"_HUMAN.tsv")
 dataSNPi=read.csv(file_ChIPi,sep='\t',header=T)
 dataSNPi=dataSNPi[dataSNPi$fdrp_bh_ref<PvalT | dataSNPi$fdrp_bh_alt<PvalT,]
 
 # If there are significant SNPs in the file
 if(nrow(dataSNPi)>1){
  SNP_DSB_Adastra_hg38.GRi=GRanges(dataSNPi[,1],IRanges(dataSNPi[,2],dataSNPi[,2]),ID=dataSNPi[,3],ref=dataSNPi$ref,alt=dataSNPi$alt,prot=protNames[i])
  ESi=as.vector(apply(dataSNPi,1,function(x){wm=which.min(x[c("fdrp_bh_ref","fdrp_bh_alt")]);y=x[c("es_mean_ref","es_mean_alt")];as.numeric(y[wm])}))
  prefAli=as.vector(apply(dataSNPi,1,function(x){wm=which.min(x[c("fdrp_bh_ref","fdrp_bh_alt")]);y=x[c("ref","alt")];y[wm]}))  
  SNP_DSB_Adastra_hg38.GRi$ES=ESi
  SNP_DSB_Adastra_hg38.GRi$prefAl=prefAli  
  SNP_DSB_Adastra_hg38.GRi$ES[prefAli==SNP_DSB_Adastra_hg38.GRi$ref]=-SNP_DSB_Adastra_hg38.GRi$ES[prefAli==SNP_DSB_Adastra_hg38.GRi$ref]
  
  if(i==1){
   SNP_DSB_Adastra_hg38.GR=SNP_DSB_Adastra_hg38.GRi
  }else{
   SNP_DSB_Adastra_hg38.GR=c(SNP_DSB_Adastra_hg38.GR,SNP_DSB_Adastra_hg38.GRi)
  }
 }
}


# Liftover to hg19
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch38to19 = import.chain(path)
SNP_DSB_Adastra.GR = liftOver(SNP_DSB_Adastra_hg38.GR, ch38to19)
SNP_DSB_Adastra.GR=unlist(SNP_DSB_Adastra.GR)

# Make occurrence table
matEXPS=dummy_cols(SNP_DSB_Adastra.GR$prot)[,-1]
colnames(matEXPS)=sapply(colnames(matEXPS),function(x){strsplit(x,'_')[[1]][2]})
values(SNP_DSB_Adastra.GR)=cbind(values(SNP_DSB_Adastra.GR)[-c(4,6)],matEXPS)

# Number of unique dsbSNPs from ADASTRA
length(unique(SNP_DSB_Adastra.GR))

# Number of dsbSNPs by protein
sort(colSums(matEXPS),decreasing=T)



###
# REMOVE SUBTELOMERIC + CENTROMERIC REGIONS

SNP_DSB_Adastra.GR=filterSNPCentroTelo(SNP_DSB_Adastra.GR)
save(SNP_DSB_Adastra.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Adastra_GR.RData")
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Adastra_GR.RData")
length(unique(SNP_DSB_Adastra.GR[SNP_DSB_Adastra.GR$ATM>0]))
write_bed(SNP_DSB_Adastra.GR,"results/allelic_imbalance/SNP_DSB/SNP_DSB_Adastra_GR.bed")
### 
# PIE PROTEINS

pie(colSums(matEXPS))

file_pie_SNPs="results/allelic_imbalance/ADASTRA_SNP/pie_SNP.pdf"
pdf(file_pie_SNPs,4,4)
pie(colSums(matEXPS))
dev.off()


### 
# ANNOTATION

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
annot_SNP_DSB=annotatePeak(SNP_DSB_Adastra.GR,tssRegion=c(-3000, 3000), TxDb=txdb)
file_annot_SNPs="results/allelic_imbalance/ADASTRA_SNP/annot_SNP.pdf"
pdf(file_annot_SNPs,6,6)
plotAnnoPie(annot_SNP_DSB)
dev.off()

######### GET CONTROL SNPs

kCtrlSNPs=10
SNP_DSB.GR=SNP_DSB_Adastra.GR
# Control SNPs as SNPs in LD with dsbSNPs (500kb)
dataSNP1KG=fread("data/SNP/1KG_SNPs_filt.bed",header=F,sep='\t')
dataSNP1KG=dataSNP1KG[dataSNP1KG$V7>0.1,] # SNPs whose minor allele frequency > 10%
dataSNP1KG=data.frame(dataSNP1KG)
SNP_ctrl.GR=GRanges(dataSNP1KG[,1],IRanges(dataSNP1KG[,3],dataSNP1KG[,3]),rs=dataSNP1KG[,4])

distLD=sample(c(-1,1),length(SNP_DSB.GR),replace=T)*250e3
SNP_DSB_s.GR=GenomicRanges::shift(SNP_DSB.GR,shift=distLD)
SNP_DSB_s.GR=resize(SNP_DSB_s.GR,width=500e3,fix="center")
olLD_ctrl=findOverlaps(SNP_ctrl.GR,SNP_DSB_s.GR)
SNP_ctrl.GR=SNP_ctrl.GR[queryHits(olLD_ctrl)]

SNP_ctrl.GR=filterSNPCentroTelo(SNP_ctrl.GR)

olDSB_ctrl=findOverlaps(SNP_ctrl.GR,SNP_DSB.GR)
if(length(olDSB_ctrl)>0){
  SNP_ctrl.GR=SNP_ctrl.GR[-queryHits(olDSB_ctrl)]
}

SNP_ctrl.GR=sort(sample(SNP_ctrl.GR,length(SNP_DSB.GR)*kCtrlSNPs))
save(SNP_ctrl.GR,file="results/allelic_imbalance/SNP_DSB/ADASTRA_SNP_ctrl_dist500kb_GR.RData")


write_bed(SNP_ctrl.GR,file="results/allelic_imbalance/SNP_DSB/ADASTRA_SNP_ctrl_dist500kb_GR.bed")




