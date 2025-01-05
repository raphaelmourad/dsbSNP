# Raphael Mourad
# 15/06/2023


# This R script is used to extract the SNPs mapped from Legube data.
# NB: our SNPs are mapped in hg19 (original bam files were mapped to hg19).



# SET DIRECTORY
setwd("/media/sauber/Elements/DSBsnp_project_seb")

# LOAD R PACKAGES
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(ggplot2)
library(gplots)
library(grDevices)

library(proxy)

library(rtracklayer)
library(plyranges)

library(data.table)
library(ChIPseeker)
library(ChIPpeakAnno)
library(Rgraphviz)
library(liftOver)

#library(sigminer)
#library(seqplots)
library(ggcorrplot)
library(tidyverse)
library(dplyr)
library(qqman)
library(hexbin)





# LOAD INHOUSE R SCRIPTS
source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")


###
# PROCESS SNP FROM LEGUBE DATABASE

# Pick experience
expe_type="ChIP-seq"
date_expe="6000SNPs_Nov2023"
# expe_type="CUT-Tag"
 #expe_type="ATAC-seq"
expe_type="DRIP-seq"
setwd("/media/sauber/Elements/DSBsnp_project_seb/")

# Threshold
PvalT=0.1
MinCoverage=10

# Read SNPs from ADASTRA pipeline
# file_SNP_expe_type=paste0("data/Legube_allelic_imbalance/",date_expe,"/mixalime_pvalues_18_05_23_renamed/",expe_type,"_noOHT_pvs.tsv")
# file_SNP_expe_type=paste0("data/Legube_allelic_imbalance/",date_expe,"/mixalime_pvalues_corrected_renamed/",expe_type,"_noOHT_pvs.tsv")
file_SNP_expe_type=paste0("data/Legube_allelic_imbalance/",date_expe,"/mixalime_pvalues_renamed/",expe_type,"_noOHT_pvs.tsv")

dataExpe=read.table(file_SNP_expe_type,sep='\t',header=T)
dataExpe=dataExpe[dataExpe$fdr_comb_pval<PvalT & dataExpe$max_cover>MinCoverage,]

# Get all experiment names
all_expes=NULL
for(i in 1:nrow(dataExpe)){
 stri=strsplit(dataExpe[i,"scorefiles"],',')[[1]]
 expei=as.vector(sapply(stri,function(x){strsplit(x,"@")[[1]][1]}))
 all_expes=c(all_expes,expei)
}
all_expes=unique(all_expes)
print(all_expes)

# Get presence/absence of SNP for each experiment
tab_list=list()
mat_SNP=NULL
for(i in 1:nrow(dataExpe)){
 stri=strsplit(dataExpe[i,"scorefiles"],',')[[1]]
 expei=as.vector(sapply(stri,function(x){strsplit(x,"@")[[1]][1]}))
 tabe=rep(0,length(setdiff(all_expes,expei)))
 names(tabe)=setdiff(all_expes,expei)
 tabi=table(expei)
 tabie=c(tabe,tabi)
 tabie=tabie[order(names(tabie))]
 mat_SNP=rbind(mat_SNP,tabie)
}
all_expes=colnames(mat_SNP)
mat_SNP=mat_SNP[,!colnames(mat_SNP)%in%c("input","tot","Tot")]
mat_SNP[mat_SNP>1]=1
colnames(mat_SNP)=gsub("-","_",colnames(mat_SNP))
if(expe_type=="ChIP-seq"){
 mat_SNP[,"pATM"]=mat_SNP[,"pATM"]+mat_SNP[,"pATM_siCtrl"]
 mat_SNP=mat_SNP[,colnames(mat_SNP)!="pATM_siCtrl"]
}
if(expe_type=="ChIP-seq"){
  mat_SNP[,"53BP1"]=mat_SNP[,"53BP1"]+mat_SNP[,"53BP1_S"]
  mat_SNP=mat_SNP[,colnames(mat_SNP)!="53BP1_S"]
}
colSums(mat_SNP)


# Number of SNP by experience
sort(colSums(mat_SNP),decreasing=T)



#correction ES
dataExpe.GR=GRanges(dataExpe[,"chr"],IRanges(dataExpe[,"end"],dataExpe[,"end"]),
                    ID=dataExpe[,"id"],ref=dataExpe[,"ref"],alt=dataExpe[,"alt"],ES=dataExpe[,"comb_es"],pref=dataExpe[,"pref_allele"])
dataExpe.GR$ES[dataExpe.GR$pref=="ref"]=-dataExpe.GR$ES[dataExpe.GR$pref=="ref"]
dataExpe.GR<-dataExpe.GR%>% as_tibble()
dataExpe.GR<-dataExpe.GR[,-10]%>%as_granges()
values(dataExpe.GR)=cbind(values(dataExpe.GR),mat_SNP)

###
# REMOVE SUBTELOMERIC + CENTROMERIC REGIONS

dataExpe.GR=filterSNPCentroTelo(dataExpe.GR)
#save(dataExpe.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Adastra_GR.RData")


###
# SAVE SNPs
if(expe_type=="ChIP-seq"){

# dsbSNPs

repairProtExpes=c("53BP1","Lig4","MDC1","NBN","pATM","RPA","SETX")

mat_dataExpe=values(dataExpe.GR)
mat_dataExpe_colsel=as.data.frame(mat_dataExpe[,colnames(mat_dataExpe)%in%repairProtExpes])
mat_dataExpe_DSB=mat_dataExpe_colsel[rowSums(mat_dataExpe_colsel)>0,]

SNP_DSB_Legube_DIVA.GR=dataExpe.GR[rowSums(mat_dataExpe_colsel)>0]

# dsbSNPs from perturbed experiments (KD, drugs)
repairProtPertubExpes=c("pATM_DRB","pATM_siSCC1","pATM_curaxin","pATM_DRB")
mat_dataExpe_colselpertub=as.data.frame(mat_dataExpe[,colnames(mat_dataExpe)%in%repairProtPertubExpes])
SNP_DSB_pertub_Legube_DIVA.GR=dataExpe.GR[rowSums(mat_dataExpe_colselpertub)>0]
#save(SNP_DSB_pertub_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_pertub_Legube_DIVA_GR.RData")

# non dsbSNPs
SNP_nonDSB_Legube_DIVA.GR=dataExpe.GR[rowSums(mat_dataExpe_colsel)==0]
#save(SNP_nonDSB_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_nonDSB_Legube_DIVA_GR.RData")

}else if(expe_type=="ATAC-seq"){

SNP_ATAC_Legube_DIVA.GR=dataExpe.GR
#save(SNP_ATAC_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_ATAC_Legube_DIVA_GR.RData")

}else if(expe_type=="DRIP-seq"){

SNP_DRIP_Legube_DIVA.GR=dataExpe.GR
#save(SNP_DRIP_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DRIP_Legube_DIVA_GR.RData")

}else if(expe_type=="CUT-Tag"){

repairProtExpes=c("pATM","RPA")
mat_dataExpe=values(dataExpe.GR)
mat_dataExpe_colsel=as.data.frame(mat_dataExpe[,colnames(mat_dataExpe)%in%repairProtExpes])

SNP_CUTnTag_Legube_DIVA.GR=dataExpe.GR[rowSums(mat_dataExpe_colsel)>0]
#save(SNP_CUTnTag_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_CUTnTag_Legube_DIVA_GR.RData")

}
#masked AsiSi

SNP_DSB_Legube_DIVA.GR
sites= read.table("data/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv",header=T)

sites<-sites[1:400,]
sites.GR<- GRanges(sites[,1],IRanges(sites[,2],sites[,3]))


prot_list=list(
  "pATM"=2000,
  "53BP1"=1000000,
  "Lig4"=20000,
  "MDC1"=1000000
)


for (i in 1: length(prot_list)){
  sites.GR<-sites.GR%>% resize(width = as.numeric(prot_list[[i]]),fix="center")
  
  test<-SNP_DSB_Legube_DIVA.GR%>% filter(get(names(prot_list)[i])==1)
  
  
  
  finalGR <- test[!test %over% sites.GR,]
  if (i==1){
    ultimeGR <-as_tibble(finalGR)
    
  }else{
    
    ultimeGR<- rbind(ultimeGR,as_tibble(finalGR))
    
  }
}

ultimeGR<-unique(ultimeGR)%>%as_granges()
#ultimeGR <- SNP_DSB_Legube_DIVA.GR[!SNP_DSB_Legube_DIVA.GR %over% sites.GR,]

length(unique(ultimeGR))


SNP_DSB_Legube_DIVA.GR<-ultimeGR
SNP_DSB_Legube_DIVA.GR
save(SNP_DSB_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_ES_corrected_GR.RData")



load("results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_ES_corrected_GR.RData")


mat_dataExpe_all=values(SNP_DSB_Legube_DIVA.GR)
mat_dataExpe_DSB=as.data.frame(mat_dataExpe_all[,!colnames(mat_dataExpe_all)%in%c("ID","ref","alt","ES","pATM_DRB","pATM_siSCC1","pATM_curaxin")])
mat_dataExpe_DSB[mat_dataExpe_DSB>1]=1
repairProtExpes=c("X53BP1","Lig4","NBN","MDC1","pATM","RPA","SETX")
repairProtExpes=c("X53BP1","Lig4","NBN","pATM","RPA","SETX")

mat_dataExpe_DSB<-mat_dataExpe_DSB%>% dplyr::select(any_of(repairProtExpes))
counts<-table(rowSums(mat_dataExpe_DSB))
counts/ length(SNP_DSB_Legube_DIVA.GR)*100
p1<-counts%>% as.data.frame()%>% ggplot()+aes(y=Freq,x=Var1)+geom_bar(stat="identity")+xlab("Number of DSB repair proteins \nsupporting a SNP")+ylab("Number of SNPs")#+
 # theme(axis.title.x = element_text(colour = "black"),   axis.title.y = element_text(colour = "black"))+
  # theme(text = element_text(size = 8,family="Arial",colour = "black"),
  #       axis.text = element_text(element_text(size = 8,family="Arial",colour = "black")))
file_coocc_SNPs="results/allelic_imbalance/Legube_SNP/plot_cooccur_gg_SNP.pdf"
pdf(file_coocc_SNPs,4,2)

print(p1)
dev.off()
### 
# PIE PROTEINS

file_pie_SNPs="results/allelic_imbalance/Legube_SNP/pie_SNP_bigersize.pdf"
pdf(file_pie_SNPs,8,8)
pie(colSums(mat_dataExpe_DSB))
dev.off()

colSums(mat_dataExpe_DSB)


### 
# NUMBER CHIPS SUPPORTING A SNP
file_coocc_SNPs="results/allelic_imbalance/Legube_SNP/plot_cooccur_SNP.pdf"
pdf(file_coocc_SNPs,4,4)
barplot(table(rowSums(mat_dataExpe_DSB[rowSums(mat_dataExpe_DSB)>0,])),ylab="Number of ChIPs\nsupporting a SNP")+
dev.off()

###########
library(ggupset)

mat_dataExpe_DSB%>%as.data.frame()%>%mutate(row=rownames(mat_dataExpe_DSB))%>%
  pivot_longer(cols=(colnames(mat_dataExpe_DSB)))%>%filter(value==1)%>%
  group_by(row) %>%
  summarize(name = list(name)) %>%
  ggplot(aes(x = name)) +
  geom_bar() +
  scale_x_upset()
### 
# ANNOTATION
library(ggupset)
library(ggimage)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
annot_SNP_DSB=annotatePeak(SNP_DSB_Legube_DIVA.GR,tssRegion=c(-3000, 3000), TxDb=txdb)
file_annot_SNPs="results/allelic_imbalance/Legube_SNP/annot_SNP.pdf"
pdf(file_annot_SNPs,6,6)
plotAnnoPie(annot_SNP_DSB)
upsetplot(annot_SNP_DSB,vennpie=T)
vennpie(annot_SNP_DSB)
dev.off()


                
load("results/allelic_imbalance/SNP_DSB/Legube_SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR
Annot_SNP_ctrl=annotatePeak(SNP_ctrl.GR,tssRegion=c(-3000, 3000), TxDb=txdb)
file_annot_SNPs="results/allelic_imbalance/Legube_SNP/annot_ctrlSNP.pdf"
pdf(file_annot_SNPs,6,6)
plotAnnoPie(Annot_SNP_ctrl)
dev.off()
f="Promoter (<=1kb)"
for (f in unique(annot_SNP_DSB@annoStat$Feature)){

p_dsb=annot_SNP_DSB@annoStat$Frequency [annot_SNP_DSB@annoStat$Feature==f]*length(SNP_DSB_Legube_DIVA.GR)/100
ptot_dsb=length(SNP_DSB_Legube_DIVA.GR)
p_ctrl=Annot_SNP_ctrl@annoStat$Frequency [Annot_SNP_ctrl@annoStat$Feature==f]*length(SNP_ctrl.GR)/100
ptot_ctrl=length(SNP_ctrl.GR)
tabSNP=rbind(c(p_dsb,ptot_dsb),c(p_ctrl,ptot_ctrl))
t<-fisher.test(tabSNP)
df<-data_frame(p_dsb,ptot_dsb)
df$p_ctrl=p_ctrl
df$ptot_ctrl=ptot_ctrl
df$pval=format(t$p.value,digits=2)
df$OR=t$estimate
df$feature=f
if(f==unique(annot_SNP_DSB@annoStat$Feature)[1]){
dff<-df
}else{
  
  dff<-rbind(dff,df)
}
}

#all intron

p_dsb=annot_SNP_DSB@annoStat$Frequency [annot_SNP_DSB@annoStat$Feature=="1st Intron"]*length(SNP_DSB_Legube_DIVA.GR)/100+annot_SNP_DSB@annoStat$Frequency [annot_SNP_DSB@annoStat$Feature=="Other Intron"]*length(SNP_DSB_Legube_DIVA.GR)/100
ptot_dsb=length(SNP_DSB_Legube_DIVA.GR)

p_ctrl=Annot_SNP_ctrl@annoStat$Frequency [Annot_SNP_ctrl@annoStat$Feature=="1st Intron"]*length(SNP_ctrl.GR)/100+Annot_SNP_ctrl@annoStat$Frequency [Annot_SNP_ctrl@annoStat$Feature=="Other Intron"]*length(SNP_ctrl.GR)/100
ptot_ctrl=length(SNP_ctrl.GR)
tabSNP=rbind(c(p_dsb,ptot_dsb),c(p_ctrl,ptot_ctrl))
t<-fisher.test(tabSNP)
df<-data_frame(p_dsb,ptot_dsb)
df$p_ctrl=p_ctrl
df$ptot_ctrl=ptot_ctrl
df$pval=format(t$p.value,digits=2)
df$OR=t$estimate
df$feature=f
df

##################Venn diagrams
library(nVennR)

SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$X53BP1==1] 
SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$pATM==1]


myV<-plotVenn(list(
  
  
  pATM=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$pATM==1],

  SETX=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$SETX==1],
  
  X53BP1=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$X53BP1==1] ,
   # RPA=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$RPA==1],
    MDC1=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$MDC1==1],
    Lig4=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$Lig4==1]
  
#

#CTCF=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$Lig4==1]

),outFile="test")
listVennRegions(myV)


myV<-plotVenn(list(
  
  
  pATM=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$pATM==1],
  
 #pATM_siSCC1=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$pATM_siSCC1==1]
  
  #pATM_DRB=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$pATM_DRB==1] 
  pATM_curaxin=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$pATM_curaxin==1] 

  
  #
  
  #CTCF=SNP_DSB_Legube_DIVA.GR$ID[SNP_DSB_Legube_DIVA.GR$Lig4==1]
  
),outFile="test")


SNP_DSB_Legube_DIVA.GR[SNP_DSB_Legube_DIVA.GR$CTCF==1]

library(nVennR)


load("/media/sauber/Elements/DSBsnp_project_seb/results/allelic_imbalance/SNP_DSB/SNP_CTCF_Legube_DIVA_GR.RData")

SNP_DSB_Legube_DIVA.GR$ID 


myV<-plotVenn(list(
  
  
  dsbSNP=SNP_DSB_Legube_DIVA.GR$ID,
  
CTCF_SNP=SNP_CTCF_Legube_DIVA.GR$ID
  
),outFile="testCTCF_DSB")



######export BED


load("results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_GR.RData")


write_bed(SNP_DSB_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_GR.bed")


# genNullSeqs(inputBedFN="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_GR.bed",nMaxTrials=10,xfold=1.5,genomeVersion="hg19",
#             outputPosFastaFN="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA.fasta",outputBedFN="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_GR.bed",outputNegFastaFN="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_CTRL.fasta",length_match_tol=0)

df<-read.table("results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_GR.bed")
df$V3<- df$V3+1

write_tsv(df,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_GR_deeptools.bed",col_names = F)



###
# CONTROL SNPS FROM 1000k GENOME
kCtrlSNPs=10
SNP_DSB.GR=SNP_DSB_Legube_DIVA.GR
# Control SNPs as SNPs in LD with dsbSNPs (500kb away)
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
save(SNP_ctrl.GR,file="results/allelic_imbalance/SNP_DSB/Legube_SNP_ctrl_dist500kb_GR.RData")


write_bed(SNP_ctrl.GR,file="results/allelic_imbalance/SNP_DSB/Legube_SNP_ctrl_dist500kb_GR.bed")

df<-read.table("results/allelic_imbalance/SNP_DSB/Legube_SNP_ctrl_dist500kb_GR.bed")
df$V3<- df$V3+1

write_tsv(df,file="results/allelic_imbalance/SNP_DSB/Legube_SNP_ctrl_dist500kb_GR_deeptools.bed",col_names = F)
