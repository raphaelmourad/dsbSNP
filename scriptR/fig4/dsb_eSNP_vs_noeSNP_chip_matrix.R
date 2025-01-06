#########################################################################
##############This script allows to compute average signal enrichment for 
#########given bigwig against the same control SNPs at dsb eSNP and dsb noeSNPs
#########################################################################
library(plyranges)
getwd()
setwd("/media/sauber/Elements/DSBsnp_project_seb/")

#load SNPs
dsb_eSNP<-read_bed("results/allelic_imbalance/SNP_DSB/dsb_eSNP.bed")
dsb_no_eSNP<-read_bed("results/allelic_imbalance/SNP_DSB/dsb_no_eSNP.bed")






###AD a widow same control for everyone


w=2000


#load control SNPs

load("results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR=sort(sample(SNP_ctrl.GR,7643))
ctrl_eSNP<-SNP_ctrl.GR
ctrl_no_eSNP<-SNP_ctrl.GR

#define list of BW

sequ_list=c("repliseq","ATAC","Bless","pATM_seb_all","H3K27ac","H3K4me3","TTseq","RAD51","lig4","Ser2P","Ser5P","DRIP","H3K9me3","H3K4me1")#H4K20me3


w_list=c(100)# can compute for several windows sizes
for(sequ_name in sequ_list){
  print(sequ_name)
  f<-list.files(paste0("dsbSNP_Av_pr/",sequ_name,"/bw/"),pattern ="bw")
  
  
  sequ<-read_bigwig(paste0("dsbSNP_Av_pr/",sequ_name,"/bw/",f))
  
  
  
  for (w in w_list)  {
    dsb_eSNP<-dsb_eSNP%>%resize(fix = "center",width = w+1)
    dsb_no_eSNP<-dsb_no_eSNP%>%resize(fix = "center",width = w+1)
    
    
    ctrl_eSNP<-ctrl_eSNP%>%resize(fix = "center",width = w+1)
    ctrl_no_eSNP<-ctrl_no_eSNP%>%resize(fix = "center",width = w+1)
    dsb_eSNP_ov<-sequ[subjectHits(findOverlaps(dsb_eSNP,sequ))]
    
    ctrl_eSNP_ov<-sequ[subjectHits(findOverlaps(ctrl_eSNP,sequ))]
    
    ratio_eSNP<-mean(dsb_eSNP_ov$score)/mean(ctrl_eSNP_ov$score)
    
    
    dsb_no_eSNP_ov<-sequ[subjectHits(findOverlaps(dsb_no_eSNP,sequ))]
    
    ctrl_no_eSNP_ov<-sequ[subjectHits(findOverlaps(ctrl_no_eSNP,sequ))]
    
    ratio_no_eSNP<-mean(dsb_no_eSNP_ov$score)/mean(ctrl_no_eSNP_ov$score)
    
    
    df<-data.frame(sequ_name,ratio_eSNP)
    df$ratio_no_eSNP<-ratio_no_eSNP
    
    df$dsb_eSNP_ov<-mean(dsb_eSNP_ov$score)
    df$dsb_no_eSNP_ov<-mean(dsb_no_eSNP_ov$score)
    df$w=w
    if(sequ_name==sequ_list[1]&w==w_list[1]){
      f.df=df
      
    }else{
      f.df=rbind(f.df,df)
    }
  }
}
#f.df2<-f.df
#f.df<-rbind(f.df2,f.df)
library(ggplot2)
library(tidyr)
library(dplyr)
f.df%>% pivot_longer(c(ratio_eSNP,ratio_no_eSNP))%>% ggplot(aes(x=name,y=sequ_name,fill=log2(value)))+geom_tile()+
  scale_fill_gradient2(low="blue4",mid="white", high = "red4",midpoint = 0)+facet_grid(~w)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

f.df%>% filter(w==100)%>%pivot_longer(c(ratio_eSNP,ratio_no_eSNP))%>% ggplot(aes(x=name,y=sequ_name,fill=log2(value)))+geom_tile()+
  scale_fill_gradient2(low="blue",mid="white", high = "red4",midpoint = 0)+facet_grid(~w)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+theme(text = element_text(size = 20),
                                                                                                                                                                  axis.text.y = element_text(angle = 0, hjust = 1)) 
pdf("results/allelic_imbalance/Heatmap_dsb_eSNP_vs_dsb_noeSNP.pdf",7,10)
f.df%>% filter(w==100)%>%pivot_longer(c(ratio_eSNP,ratio_no_eSNP))%>% ggplot(aes(x=name,y=sequ_name,fill=log2(value)))+geom_tile()+
  scale_fill_gradient2(low="blue",mid="white", high = "red4",midpoint = 0)+facet_grid(~w)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+theme(text = element_text(size = 35),
                                                                                                                                                                  axis.text.y = element_text(angle = 0, hjust = 1)) 

dev.off()

f.df%>% pivot_longer(c(dsb_eSNP_ov,dsb_no_eSNP_ov))%>% ggplot(aes(x=name,y=sequ_name,fill=(value)))+geom_tile()+
  scale_fill_gradient2(low="white",mid="red",high = "red4",midpoint = max(c(f.df$dsb_eSNP_ov,f.df$dsb_no_eSNP_ov))/2)+
  facet_grid(~w)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

f.df%>% pivot_longer(c(dsb_eSNP_ov,dsb_no_eSNP_ov))%>% ggplot(aes(x=name,y=sequ_name,fill=log2(value)))+geom_tile()+#scale_fill_gradient(low="white",high = "red")
  scale_fill_gradient2(low="white",mid="red",high = "red4",midpoint = max(c(log2(f.df$dsb_eSNP_ov),log2(f.df$dsb_no_eSNP_ov)))/2)+
  facet_grid(~w)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


f.df%>% filter(w==10000)%>%pivot_longer(c(dsb_eSNP_ov,dsb_no_eSNP_ov))%>% ggplot(aes(x=name,y=sequ_name,fill=(value)))+geom_tile()+#scale_fill_gradient(low="white",high = "red")
  scale_fill_gradient2(low="white",mid="red",high = "red4",midpoint = max(c((f.df$dsb_eSNP_ov),(f.df$dsb_no_eSNP_ov)))/2)+
  facet_grid(~w)


f.df%>% filter(w==10000)%>%pivot_longer(c(dsb_eSNP_ov,dsb_no_eSNP_ov))%>% ggplot(aes(x=name,y=sequ_name,fill=log2(value)))+geom_tile()+#scale_fill_gradient(low="white",high = "red")
  scale_fill_gradient2(low="white",mid="red",high = "red4",midpoint = max(c(log2(f.df$dsb_eSNP_ov),log2(f.df$dsb_no_eSNP_ov)))/2)+
  facet_grid(~w)
######################Parralell version

library(parallel)
dsb_eSNP<-read_bed("results/allelic_imbalance/SNP_DSB/dsb_eSNP.bed")
dsb_no_eSNP<-read_bed("results/allelic_imbalance/SNP_DSB/dsb_no_eSNP.bed")


ctrl_eSNP<-read_bed("results/allelic_imbalance/SNP_DSB/ctrl_eSNP.bed")
ctrl_no_eSNP<-read_bed("results/allelic_imbalance/SNP_DSB/ctrl_no_eSNP.bed")


load("results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR=sort(sample(SNP_ctrl.GR,7643))
ctrl_eSNP<-SNP_ctrl.GR
ctrl_no_eSNP<-SNP_ctrl.GR
sequ_name="repliseq"
sequ_list=c("repliseq","ATAC","Bless","pATM_seb_all","H3K27ac","H3K4me3","TTseq","RAD51","lig4","Ser2P","Ser5P","DRIP")
w_list=c(1,2000,10000)

f.df<-mclapply(sequ_list,mc.cores = 3,function(sequ_name){
  print(sequ_name)
  f<-list.files(paste0("dsbSNP_Av_pr/",sequ_name,"/bw/"),pattern ="bw")
  
  
  sequ<-read_bigwig(paste0("dsbSNP_Av_pr/",sequ_name,"/bw/",f))
  
  
  
  lapply(w_list,function(w) {
    dsb_eSNP<-dsb_eSNP%>%resize(fix = "center",width = w+1)
    dsb_no_eSNP<-dsb_no_eSNP%>%resize(fix = "center",width = w+1)
    
    
    ctrl_eSNP<-ctrl_eSNP%>%resize(fix = "center",width = w+1)
    ctrl_no_eSNP<-ctrl_no_eSNP%>%resize(fix = "center",width = w+1)
    dsb_eSNP_ov<-sequ[subjectHits(findOverlaps(dsb_eSNP,sequ))]
    
    ctrl_eSNP_ov<-sequ[subjectHits(findOverlaps(ctrl_eSNP,sequ))]
    
    ratio_eSNP<-mean(dsb_eSNP_ov$score)/mean(ctrl_eSNP_ov$score)
    
    
    dsb_no_eSNP_ov<-sequ[subjectHits(findOverlaps(dsb_no_eSNP,sequ))]
    
    ctrl_no_eSNP_ov<-sequ[subjectHits(findOverlaps(ctrl_no_eSNP,sequ))]
    
    ratio_no_eSNP<-mean(dsb_no_eSNP_ov$score)/mean(ctrl_no_eSNP_ov$score)
    
    
    df<-data.frame(sequ_name,ratio_eSNP)
    df$ratio_no_eSNP<-ratio_no_eSNP
    
    df$dsb_eSNP_ov<-mean(dsb_eSNP_ov$score)
    df$dsb_no_eSNP_ov<-mean(dsb_no_eSNP_ov$score)
    df$w=w
 df
  })
})
library(tidyr)
library(dplyr)

f.df[[1]]#%>%bind_rows()
f.df<-f.df%>%bind_rows()
f.df%>% pivot_longer(c(ratio_eSNP,ratio_no_eSNP))%>% ggplot(aes(x=name,y=sequ_name,fill=log2(value)))+geom_tile()+scale_fill_gradient2(low="blue4",mid="white", high = "red",midpoint = 0)+facet_grid(~w)

