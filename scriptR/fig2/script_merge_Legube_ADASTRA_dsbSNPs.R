# Raphael Mourad
# 15/06/2023


# This R script is used to merge dsbSNPs from Legube and ADASTRA database into one set of dsbSNPs.


# SETUP DIRECTORY
setwd("/media/sauber/Elements/DSBsnp_project_seb")


# LOAD LIBRARIES
library(GenomicRanges)





###
# LOAD SNPS

load("results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_ES_corrected_GR.RData")
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Adastra_GR.RData")

# mv(from = "SNP_DSB.GR",to =paste0("test","SNP_DSB.GR"))
# testSNP_DSB.GR
###
# OVERLAP BETWEEN LEGUBE DSB SNPS AND ADASTRA DSB SNPS

olSNPS=findOverlaps(unique(SNP_DSB_Legube_DIVA.GR),unique(SNP_DSB_Adastra.GR))
length(olSNPS)


###
# COMBINE OUR SNPS WITH ADASTRA SNPS

SNP_DSB.GR=c(SNP_DSB_Legube_DIVA.GR,SNP_DSB_Adastra.GR)
SNP_DSB.GR=sort(unique(SNP_DSB.GR))
valSNP=as.data.frame(values(SNP_DSB.GR))
valSNP[is.na(valSNP)]=0
values(SNP_DSB.GR)=valSNP
save(SNP_DSB.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")


write.table(as.data.frame(SNP_DSB.GR)[,1:3],file="results/allelic_imbalance/SNP_DSB/SNP_DSB_GR_corrected.bed",row.names=F,col.names=F,quote=F,sep='\t')



