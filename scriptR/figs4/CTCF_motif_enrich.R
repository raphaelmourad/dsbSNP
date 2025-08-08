# Sébastien AUBER
# 07/05/2025
#this script compute the enrichement for a given motif 'CTCF) between 
#control and dsbSNPs

DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)

library("Signac")
library("BSgenome.Hsapiens.UCSC.hg19")
library("JASPAR2024")
library("TFBSTools")
library("motifmatchr")
library('dplyr')
library(tidyr)
library('BRGenomics')
library(gplots)
library(proxy)
library(ggplot2)
library(plyranges)
library(parallel)
library(dplyr)
library(plyranges)

#get jaspar 2024 motifs
JASPAR2024 <- JASPAR2024()

pfm<-getMatrixSet(x=db(JASPAR2024),
                  opts=list(species=9606))
load("results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")

load("results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist200_500kb_GR.RData")


########################################################################"



getSequences=function(region.GR,window=201,genome=BSgenome.Hsapiens.UCSC.hg19){
  #get sequence for at given coordinates
  region.GRr=resize(region.GR,fix="center",width=window)
  region.seq=as.character(getSeq(genome, names=seqnames(region.GRr), 
                                 start=start(region.GRr), end=end(region.GRr)))
  
  SNPposrel=ceiling((window+1)/2)
  region.seqRef=region.seq
 
  return(list(as.character(region.seqRef)))
}

win=25
get_df_all_motif<-function(SNP_DSB.GR,mode="no_mode", win=25,pm=  5e-03){
  ###get all motif for a GR
  #resize
  SNP_DSB.GR_W<-SNP_DSB.GR%>% anchor_center()%>%mutate(width = win) 
  #get sequance
  seqSNPs.seq=getSequences(SNP_DSB.GR,window=win,genome=BSgenome.Hsapiens.UCSC.hg19)
  seqSNPs.seqRef=DNAStringSet(seqSNPs.seq[[1]])
  names(seqSNPs.seqRef)<-SNP_DSB.GR$ID
  #get motifs
  Annotated_DSB_SNP.seqRef_pos<-matchMotifs(pfm,seqSNPs.seqRef,out="position",ranges=SNP_DSB.GR_W,p.cutoff = pm ) #ajouter le range diminue beaucoup le temps de calcul
  #unlist and add ID
  U_Annotated_DSB_SNP.seqRef_pos<-unlist(Annotated_DSB_SNP.seqRef_pos)
  U_Annotated_DSB_SNP.seqRef_pos$motifs<-names(U_Annotated_DSB_SNP.seqRef_pos)
  
  motif_df_ref<-join_overlap_left(SNP_DSB.GR,U_Annotated_DSB_SNP.seqRef_pos)
  motif_df_ref<-motif_df_ref%>%as_tibble()
  
  motif_df_ref$score[is.na(motif_df_ref$score)]=0
  motif_df_ref<-motif_df_ref%>%dplyr::rename(motif_score_ref=score)
  
  
  
  motif_df_final<-(motif_df_ref)#%>%mutate(motif_delta_scores_alt_minus_ref=motif_score_alt-motif_score_ref)
  motif_df_final$motif_score_ref[is.na(motif_df_final$motif_score_ref)]=0
  
  motif_df_final$motifs[is.na(motif_df_final$motifs)]="noMotif"
  

  #get TF names
  nam<-unique(motif_df_final$motifs[motif_df_final$motifs!="noMotif"])

  namlist<-mclapply(nam,mc.cores=8,
                    function(m){ getMatrixByID (db(JASPAR2024),ID=m )})
  

  motifs<-c()
  motif_tf<-c()
  for (i in (1:length(namlist))){
    print(i)
    motif_tf[i]<-TFBSTools:: name(namlist[[i]])
    motifs[i]<-ID(namlist[[i]])
  }
  
  id_motifs<-data.frame(motifs,motif_tf)
  #join name and motif df
  motif_df_final<-left_join(motif_df_final,id_motifs)
  
  motif_df_final$motif_tf[is.na(motif_df_final$motif_tf)]="noMotif"
  motif_df_final$motif_score_ref=round(motif_df_final$motif_score_ref,3)
  motif_df_final$prot=mode
  colnames(motif_df_final)
  
  motif_df_final<-motif_df_final%>%select("seqnames",   "start",     "end", "width", "strand", "ID",   "motifs","motif_score_ref","prot","motif_tf")
    df_final_motifs<-motif_df_final
  
  #return dataframe with motif, motif score,and delta
  return(df_final_motifs)
}


#SNP_ctrl.GR<-SNP_ctrl.GR%>%as_tibble()%>%dplyr::rename(ID=rs)%>%as_granges()
#get dsb motifs
df_dsb<-get_df_all_motif(SNP_DSB.GR,mode="dsb",win=25,pm=5e-3)

#get ctrl motifs
df_ctrl<-get_df_all_motif(SNP_ctrl.GR,mode="ctrl",win=25,5e-3)

fact=length(SNP_ctrl.GR)/length(SNP_DSB.GR) 
df<-rbind(df_dsb,df_ctrl)
df_summ<-df%>%group_by(prot,motif_tf)%>%mutate(val=1)%>%summarise(n=sum(val))
df_summ$n[df_summ$prot=="ctrl"]=df_summ$n[df_summ$prot=="ctrl"]/fact
#select one motif (CTCF)
motif="CTCF"
datas<-df%>%group_by(prot,motif_tf)%>%filter(motif_tf==motif)

#compute n
uniq_sum<-df%>%select(ID,prot,motif_tf)%>%filter(motif_tf==motif)%>%unique()%>%mutate(val=1)%>%group_by(prot,motif_tf)%>%summarise(n=sum(val))
uniq_sum$n[uniq_sum$prot=="ctrl"]=uniq_sum$n[uniq_sum$prot=="ctrl"]/fact #correct sizes

#compute prop test
prop.test(c(df_summ$n[df_summ$prot=="dsb"&df_summ$motif_tf==motif],df_summ$n[df_summ$prot=="ctrl"&df_summ$motif_tf==motif]),c(length(SNP_DSB.GR),length(SNP_ctrl.GR)/fact))
pt<-prop.test(c(uniq_sum$n[uniq_sum$prot=="dsb"&uniq_sum$motif_tf==motif],uniq_sum$n[uniq_sum$prot=="ctrl"&uniq_sum$motif_tf==motif]),c(length(SNP_DSB.GR),length(SNP_ctrl.GR)/fact),alternative = "greater")

#comptute Fold change
FC=(uniq_sum$n[uniq_sum$prot=="dsb"&uniq_sum$motif_tf==motif]/length(SNP_DSB.GR))/ (uniq_sum$n[uniq_sum$prot=="ctrl"&uniq_sum$motif_tf==motif]/(length(SNP_ctrl.GR)/fact))
FC
#plot
p<-uniq_sum%>%filter(motif_tf==motif)%>%ggplot()+aes(x=prot,y=n,fill=reorder(prot,-n))+geom_bar(stat = "identity")+ggtitle(paste0(motif," p=" ,pt$p.value))
p
length(SNP_DSB.GR)-uniq_sum$n[uniq_sum$prot=="dsb"&uniq_sum$motif_tf==motif]

#export pdf
pdf("results/allelic_imbalance/TF_motif/CTCF_motifs_n.pdf")
p
dev.off()

