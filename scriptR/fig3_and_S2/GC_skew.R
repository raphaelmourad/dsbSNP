#Sébastien AUBER
#07/08/205

#This script compute delta GCskew for both allele of a SNP within a given window
# at promoters (-500,+1500)

library(stringr)
library(rstatix)
library(plyranges)
library(ChIPseeker)
library(tidyr)
library(dplyr)
library(ggplot2)

DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)

#load inhouse functions
source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")


getAllelesSeq=function(region.GR,window=201,genome=BSgenome.Hsapiens.UCSC.hg19){
  #this fonction compute sequence for both allele
  region.GRr=resize(region.GR,fix="center",width=window)
  region.seq=as.character(getSeq(genome, names=seqnames(region.GRr), 
                                 start=start(region.GRr), end=end(region.GRr)))
  
  SNPposrel=ceiling((window+1)/2)
  region.seqRef=region.seq
  substring(region.seqRef,SNPposrel,SNPposrel)=as.character(region.GRr$ref)
  region.seqAlt=region.seq
  substring(region.seqAlt,SNPposrel,SNPposrel)=as.character(region.GRr$alt)
  
  return(list(as.character(region.seqRef),as.character(region.seqAlt)))
}


#get promoters 
prom<- getPromoters(upstream = 500, downstream = 1500)
#prom<- getPromoters()

getAllelesSeq=function(region.GR,window=201,genome=BSgenome.Hsapiens.UCSC.hg19){
  
  region.GRr=resize(region.GR,fix="center",width=window)
  region.seq=as.character(getSeq(genome, names=seqnames(region.GRr), 
                                 start=start(region.GRr), end=end(region.GRr)))
  
  SNPposrel=ceiling((window+1)/2)
  region.seqRef=region.seq
  substring(region.seqRef,SNPposrel,SNPposrel)=as.character(region.GRr$ref)
  region.seqAlt=region.seq
  substring(region.seqAlt,SNPposrel,SNPposrel)=as.character(region.GRr$alt)
  
  return(list(as.character(region.seqRef),as.character(region.seqAlt)))
}

#load SNPs
load(file="results/allelic_imbalance/SNP_DSB/SNP_pATM_ES_corrected_GR.RData")
load("results/allelic_imbalance/SNP_DSB/SNP_DRIP_Legube_DIVA_GR.RData") 


#define SNPs to use
list_mode=c("pATM","DRIP")

for (mode in list_mode){
  
  DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
  setwd(DSB_path)
#Define SNP of interest accordinf to mode  
  if(mode=="DRIP"){
    
SNP_interest.GR=SNP_DRIP_Legube_DIVA.GR[SNP_DRIP_Legube_DIVA.GR$DRIP>0|SNP_DRIP_Legube_DIVA.GR$DRIP_C1>0]#pATM
  }
  
  if(mode=="pATM"){
    pATM<-pATM[!pATM$ES==0]
    SNP_interest.GR=pATM#pATM
  }
  



  
  w_list=c(1501)#define a list of window on xhivh compute GCskew
  for(window in w_list){
    #get DNA sequences
  seqSNPs.seq=getAllelesSeq(SNP_interest.GR,window=window,genome=BSgenome.Hsapiens.UCSC.hg19)
  seqSNPs.seqRef=DNAStringSet(seqSNPs.seq[[1]])
  seqSNPs.seqAlt=DNAStringSet(seqSNPs.seq[[2]])
  
  seqSNPs.seqAlt<-seqSNPs.seqAlt%>%as.data.frame()
  seqSNPs.seqAlt$rsID<-SNP_interest.GR$ID
  
  seqSNPs.seqRef<-seqSNPs.seqRef%>%as.data.frame()
  seqSNPs.seqRef$rsID<-SNP_interest.GR$ID
  
for(i in 1:nrow(seqSNPs.seqRef))  {
  
  #compute GCskew for ref allele
 ID=seqSNPs.seqRef$rsID[i] 
 G_ref=length(which(str_split(seqSNPs.seqRef$x[i],"")[[1]]=="G"))
 C_ref=length(which(str_split(seqSNPs.seqRef$x[i],"")[[1]]=="C"))
 skew_ref=(G_ref-C_ref)/(G_ref+C_ref)

  #compute GCskew for alt allele
  G_alt=length(which(str_split(seqSNPs.seqAlt$x[i],"")[[1]]=="G"))
  C_alt=length(which(str_split(seqSNPs.seqAlt$x[i],"")[[1]]=="C"))
  skew_alt=(G_alt-C_alt)/(G_alt+C_alt)

  #stock in a dataframe
dat<-data.frame(ID,G_ref)
dat$C_ref<-C_ref
dat$skew_ref<-skew_ref
dat$G_alt<-G_alt
dat$C_alt<-C_alt

dat$skew_alt<-skew_alt

if(i==1){
  df=dat
}else{
  df=rbind(df,dat)
}
}
  
  
  for (t in c(0)){
  
  #merge with SNP of interest data
  merged<-merge.data.frame(SNP_interest.GR%>%as_data_frame()%>%dplyr::select(seqnames,start,end,ID,ES,alt,ref),df  )
  ###compute delta GCskew
  merged$delta=merged$skew_alt-merged$skew_ref
  
  merged$ratio=merged$skew_alt/merged$skew_ref
  #class according ES
  merged<-merged[!(merged$ES<t & merged$ES > -t),]
  merged<-merged[merged$ES!=t,]
  merged$class="NA"
  merged$class[merged$ES>t]=paste0("ES > ",t)
  merged$class[merged$ES< -t]=paste0("ES < ",-t)

#  
  nofilter_p<-merged%>%wilcox_test(delta ~ class)
  nofilter_p$p
 no_filter<- merged%>%
    ggplot()+aes(x=class,y=delta,fill=reorder(class,delta))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2)+#+facet_wrap(~highES)
  ggtitle(paste0(mode,"SNPs\nGC Skew","\np=",  nofilter_p$p
))
  pdf(paste0("results/allelic_imbalance/GC_skew/delta_nofilter_",t,"_window",window,"_",mode,".pdf"),4,6 )
 print(no_filter)
 dev.off()
  
 
 ######Compute at promoters only
 
 merged_GR<- merged%>%as_granges()
 #overlap with promoter
 merged_prom<-merged_GR[queryHits(findOverlaps(merged_GR,prom))]
 #class according ES
 merged_prom$class[merged_prom$ES> t]=paste0(merged_prom$class[merged_prom$ES> t]," n=",length(merged_prom$class[merged_prom$ES> t]))
 merged_prom$class[merged_prom$ES< -t]=paste0(merged_prom$class[merged_prom$ES< -t]," n=",length(merged_prom$class[merged_prom$ES< -t]))
 
 #compute wilcoxon test
 wilc_prom<-merged_prom%>%as_tibble()%>%
   wilcox_test(delta ~ class)
 
 wilc_prom

   #plot
   p<- merged_prom%>%as_tibble()%>%ggplot()+aes(x=class,y=delta,fill = class)+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2)+#+facet_wrap(~highES)
 ggtitle(paste0(mode," SNPs in promoters (-500;+1500)\n","window=",window,"\np=",wilc_prom$p,"\n"))
 print(p)
 
 #plot at several scale and export pdf
 pdf(paste0("results/allelic_imbalance/GC_skew/delta_nofilter_",t,"window_",window,"_promoters_",mode,".pdf"),4,6 )
 print(p)
 print(p+coord_cartesian(ylim = c(-0.03,0.03)))
 print(p+coord_cartesian(ylim = c(-0.01,0.01)))
 
 print(p+coord_cartesian(ylim = c(-0.000001*window,0.000001*window)))
 
 dev.off()
 
 
 }
}
}

