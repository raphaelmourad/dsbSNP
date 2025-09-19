#Sébastien AUBER
#06/08/2025
#This script Allows to make pair pannels

library(dplyr)
library(tidyr)
library(matsbyname)
library(psych)

expe_type="ChIP-seq"
date_expe="6000SNPs_Nov2023"


setwd("/media/sauber/Elements/DSBsnp_project_seb/")

source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")

suppressPackageStartupMessages(library(plyranges))
# Threshold
PvalT=0.1
MinCoverage=1

DSB=T #stay at dsbSNP

# set aggregation of ES on keeping the hifher absolute value
multi_ES_process="biger"#"mean"#

# Read SNPs from ADASTRA pipeline
file_SNP_expe_type=paste0("data/Legube_allelic_imbalance/",date_expe,"/mixalime_pvalues_renamed/",expe_type,"_noOHT_pvs.tsv")

dataExpe=read.table(file_SNP_expe_type,sep='\t',header=T)
dataExpe=dataExpe[dataExpe$fdr_comb_pval<PvalT & dataExpe$max_cover>MinCoverage,]


all_expes=NULL
for(i in 1:nrow(dataExpe)){
  stri=strsplit(dataExpe[i,"scorefiles"],',')[[1]]
  expei=as.vector(sapply(stri,function(x){strsplit(x,"@")[[1]][1]}))
  all_expes=c(all_expes,expei)
}
all_expes=unique(all_expes)
print(all_expes)
tab_list=list()
mat_SNP=NULL
#i=55
for(i in 1:nrow(dataExpe)){
  print((i/nrow(dataExpe))*100)
  stri=strsplit(dataExpe[i,"scorefiles"],',')[[1]]
  expei=as.vector(sapply(stri,function(x){strsplit(x,"@")[[1]][1]}))
  tabe=rep(0,length(setdiff(all_expes,expei)))
  names(tabe)=setdiff(all_expes,expei)
  if(dataExpe[i,"pref_allele"]=="ref"){
    ESi=-as.numeric(strsplit(dataExpe[i,"ref_es"],',')[[1]])
    #  ESi=as.numeric(strsplit(dataExpe[i,"ref_es"],',')[[1]])
    
  }else{
    ESi=as.numeric(strsplit(dataExpe[i,"alt_es"],',')[[1]])
    
  }
  dup= which(duplicated(expei))
  if (length(dup)>0){
    if (multi_ES_process=="biger"){
      for(j in length(dup)) 
        if(abs(ESi[dup[[j]]])>abs(ESi[dup[[j]]-1])){
          ESi=ESi[-(dup[[j]]-1)]
          expei=expei[-(dup[[j]]-1)]
        }else{
          ESi=ESi[-(dup[[j]])]
          expei=expei[-(dup[[j]])]
        }
    }
    
    if (multi_ES_process=="mean"){
      for(j in length(dup)) {
        ESi[j]=mean(ESi[ which(expei==expei[dup[[j]]])])
        expei=expei[-(dup[[j]])]
        
      }
    }
  }
  ESi= as.matrix(t(as.data.frame(ESi)))
  names(ESi)=expei
  tabi=ESi
  tabie=c(tabe,tabi)
  tabie=tabie[order(names(tabie))]
  mat_SNP=rbind(mat_SNP,tabie)
  
}

all_expes=colnames(mat_SNP)
mat_SNP=mat_SNP[,!colnames(mat_SNP)%in%c("input","tot","Tot")]
#mat_SNP[mat_SNP>1]=1
colnames(mat_SNP)=gsub("-","_",colnames(mat_SNP))
colSums(mat_SNP)

# Number of SNP by experience
sort(colMeans(mat_SNP),decreasing=T)
dataExpe.GR=GRanges(dataExpe[,"chr"],IRanges(dataExpe[,"end"],dataExpe[,"end"]),
                    ID=dataExpe[,"id"],ref=dataExpe[,"ref"],alt=dataExpe[,"alt"],ES=dataExpe[,"comb_es"])
values(dataExpe.GR)=cbind(values(dataExpe.GR),mat_SNP)
###
# REMOVE SUBTELOMERIC + CENTROMERIC REGIONS

dataExpe.GR=filterSNPCentroTelo(dataExpe.GR)

dataExpe.GR$pATM[dataExpe.GR$pATM_siCtrl!=0]=(dataExpe.GR$pATM[dataExpe.GR$pATM_siCtrl!=0]+dataExpe.GR$pATM_siCtrl[dataExpe.GR$pATM_siCtrl!=0])/2
dataExpe.GR<-dataExpe.GR%>%select(-pATM_siCtrl)
if(DSB==T){
  if(expe_type=="ChIP-seq"){
    
    # dsbSNPs
    repairProtExpes=c("53BP1","Lig4","pATM","RPA","SETX")
    mat_dataExpe=values(dataExpe.GR)
    mat_dataExpe_colsel=as.data.frame(mat_dataExpe[,colnames(mat_dataExpe)%in%repairProtExpes])
    mat_dataExpe_DSB=mat_dataExpe_colsel[rowSums(mat_dataExpe_colsel)!=0,]
    
    dataExpe.GR=dataExpe.GR[rowSums(mat_dataExpe_colsel)!=0]

  }
}
# #masked AsiSi


sites= read.table("data/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv",header=T)

sites<-sites[1:400,]
sites.GR<- GRanges(sites[,1],IRanges(sites[,2],sites[,3]))


prot_list=list(
  "pATM"=2000,
  "53BP1"=1000000
)


for (i in 1: length(prot_list)){
  sites.GR<-sites.GR%>% resize(width = as.numeric(prot_list[[i]]),fix="center")
  
  test<-dataExpe.GR#%>% filter(get(names(prot_list)[i])!=0)
  
  
  
  finalGR <- test[!test %over% sites.GR,]
  if (i==1){
    ultimeGR <-as_tibble(finalGR)
    
  }else{
    
    ultimeGR<- rbind(ultimeGR,as_tibble(finalGR))
    
  }
}

ultimeGR<-unique(ultimeGR)%>%as_granges()
length(unique(ultimeGR))

dataExpe.GR<-ultimeGR

if(DSB==T){ #select only dsbSNP
  load("results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_ES_corrected_GR.RData")
  hits<-findOverlaps(dataExpe.GR,SNP_DSB_Legube_DIVA.GR)
  dataExpe.GR<-dataExpe.GR[queryHits(hits)]
}

## create ES matrix
matExpe<-dataExpe.GR%>%as_tibble()%>%dplyr::select(-c(seqnames,start,end,ID,ref,ES,width,strand,alt))

matExpe<-matExpe%>%select(c("X53BP1","pATM","Lig4","RPA","SETX"))
matExpe[matExpe==0]<-NA

matExpe<-matExpe[,c( "pATM","X53BP1","Lig4","RPA","SETX")] 

### compute and save pair pannels (scatter and correlation)
dir.create("results/allelic_imbalance/correlation")
pdf("results/allelic_imbalance/correlation/heatmap_scaterplot_dsbSNP.pdf")
pairs.panels(matExpe, lm=TRUE,ellipses = F,density=F,smoother=F,stars=T,hist.border="white",hist.col = "white")
dev.off()# 
