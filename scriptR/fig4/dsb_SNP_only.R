

#This script compute dsb SNP that are only dsbSNPs

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(ChIPpeakAnno)
library(ChIPseeker)
library(readr)
library(UpSetR)
library("ggupset")

#load the SNPs 
load(file="results/allelic_imbalance/SNP_DSB/SNP_ATAC_Legube_DIVA_GR.RData")
load("/media/sauber/Elements/DSBsnp_project_seb/results/allelic_imbalance/SNP_DSB/SNP_CTCF_Legube_DIVA_GR.RData")
load("/media/sauber/Elements/DSBsnp_project_seb/results/allelic_imbalance/SNP_DSB/SNP_ATAC_Legube_DIVA_GR.RData")
load("/media/sauber/Elements/DSBsnp_project_seb/results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
load("/media/sauber/Elements/DSBsnp_project_seb/results/allelic_imbalance/SNP_DSB/SNP_DRIP_Legube_DIVA_GR.RData")

load("data/eQTL_GTEx.RData")
ueQTL.GR<-unique(eQTL.GR)
SNP_DSB.GR$name=paste0(seqnames(SNP_DSB.GR),ranges(SNP_DSB.GR))
ueQTL.GR$name=paste0(seqnames(ueQTL.GR),ranges(ueQTL.GR))
SNP_CTCF_Legube_DIVA.GR$name=paste0(seqnames(SNP_CTCF_Legube_DIVA.GR),ranges(SNP_CTCF_Legube_DIVA.GR))
SNP_ATAC_Legube_DIVA.GR$name=paste0(seqnames(SNP_ATAC_Legube_DIVA.GR),ranges(SNP_ATAC_Legube_DIVA.GR))
SNP_DRIP_Legube_DIVA.GR$name=paste0(seqnames(SNP_DRIP_Legube_DIVA.GR),ranges(SNP_DRIP_Legube_DIVA.GR))
load("data/eQTL_pancanQTL.RData")
u_cancer_eQTL.GR<-unique(eQTL.GR)
u_cancer_eQTL.GR$name=paste0(seqnames(u_cancer_eQTL.GR),ranges(u_cancer_eQTL.GR))

x <- list(
  dsbSNP=SNP_DSB.GR$name,
  DRIP=SNP_DRIP_Legube_DIVA.GR$name,
  
  CTCF_SNP=SNP_CTCF_Legube_DIVA.GR$name,
  atacSNP=SNP_ATAC_Legube_DIVA.GR$name,
  eSNP=ueQTL.GR$name,
  cancer_eQTL=u_cancer_eQTL.GR$name
)
ggVennDiagram(x, label_alpha = 0)+  ggplot2::scale_fill_gradient(low="white",high = "red")+theme(legend.position = "none")
#install.packages('UpSetR')


ueQTL.GR$type="eSNP"
SNP_DSB.GR$type="dsbSNP"
SNP_CTCF_Legube_DIVA.GR$type="CTCF"
SNP_ATAC_Legube_DIVA.GR$type="atacSNP"
u_cancer_eQTL.GR$type="cancer_eQTL"
SNP_DRIP_Legube_DIVA.GR$type="DRIP"

ordered<-rbind(
  SNP_DRIP_Legube_DIVA.GR  %>%as.data.frame()%>%dplyr::select(name,type),
  u_cancer_eQTL.GR%>%as.data.frame()%>%dplyr::select(name,type),
  ueQTL.GR%>%as.data.frame()%>%dplyr::select(name,type),
  SNP_DSB.GR%>%as.data.frame()%>%dplyr::select(name,type),
  SNP_CTCF_Legube_DIVA.GR %>%as.data.frame()%>%dplyr::select(name,type),
  SNP_ATAC_Legube_DIVA.GR  %>%as.data.frame()%>%dplyr::select(name,type)
  
  
)








ordered %>%
  group_by(name) %>%
  summarize(type = list(type)) %>%
  distinct(name, .keep_all=TRUE) %>%
  ggplot(aes(x=type)) +
  geom_bar() +
  scale_x_upset()

library(tidyr)
library(dplyr)

types=unique(ordered$type)
ordered$value=1
unique(ordered$name)
or<-unique(ordered)%>%pivot_wider(names_from = type,values_from = value,values_fill=0)

library('ComplexUpset')
or[types] = or[types] == 1
t(head(or[types], 3))
ComplexUpset::upset(or, types, name='type', width_ratio=0.1)
upset(or, types,mode="exclusive_intersection", name='type', width_ratio=0.1, max_size=4158)
?upset
#dsb snp only
or2<-or
or2[or2$dsbSNP==FALSE,]=NA
or2<-na.omit(or2)
ComplexUpset::upset(or2, types,mode="exclusive_intersection", name='type', width_ratio=0.1
                    
                    
ComplexUpset::upset(or2, types,mode="exclusive_intersection", name='type', width_ratio=0.1, max_degree=2)



or3<-unique(ordered)%>%pivot_wider(names_from = type,values_from = value,values_fill=0)
or3[or3$dsbSNP==0,]=NA
or3<-na.omit(or3)
or3<-or3[rowSums(or3[types])==1,]

dsb_only<-SNP_DSB.GR[SNP_DSB.GR$name %in%or3$name]


############"
library(ChIPseeker)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
annot_SNP_DSB=annotatePeak(dsb_only,tssRegion=c(-3000, 3000), TxDb=txdb)
plotAnnoPie(annot_SNP_DSB)




matEXP<-values(dsb_only)[repairProtExpes]%>%as.matrix()
matEXP[matEXP>1]=1
pie(colSums(matEXP))





