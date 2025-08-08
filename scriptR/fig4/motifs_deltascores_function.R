##################################################################
#Sébastien AUBER
#10/07/2025
#This script is made to compute dsbSNP effect on TF motif score
###################################################################

library("Seurat")
#library("GenomeInfoDb")
library("Signac")
library("BSgenome.Hsapiens.UCSC.hg19")
library("JASPAR2024")
library("TFBSTools")
library("motifmatchr")
library('dplyr')
library(tidyr)
#BiocManager::install("motifmatchr")
library('BRGenomics')
library(gplots)
library(proxy)
library(parallel)
library(ggplot2)
library(readr)

library(plyranges)

#load JASPAR 2024 Library


JASPAR2024 <- JASPAR2024()
db(JASPAR2024)

#install.packages(("/home/sauber/Downloads/JASPAR2022_0.99.7.tar.gz"), repos=NULL)

DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)
# load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
# load(file="results/allelic_imbalance/SNP_DSB/SNP_pATM_ES_corrected_GR.RData")
# load(file="results/allelic_imbalance/SNP_DSB/SNP_RAD51_Adastra_GR.RData")
# load("results/allelic_imbalance/SNP_DSB/SNP_X53BP1_ES_corrected_GR.RData")
# mode="pATM"
# mode="RAD51"
# mode="53BP1"
list_mode=c("pATM","RAD51","53BP1")
list_mode=c("pATM")

for (mode in list_mode){
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
load(file="results/allelic_imbalance/SNP_DSB/SNP_pATM_ES_corrected_GR.RData")
load(file="results/allelic_imbalance/SNP_DSB/SNP_RAD51_Adastra_GR.RData")
load("results/allelic_imbalance/SNP_DSB/SNP_X53BP1_ES_corrected_GR.RData")
#load the corresponding GR
if(mode=="pATM"){
  pATM<-pATM[!pATM$ES==0]
  SNP_DSB.GR=pATM#pATM
}
if(mode=="53BP1"){
  X53BP1<-X53BP1[!X53BP1$ES==0]
  SNP_DSB.GR=X53BP1#pATM
}

if(mode=="RAD51"){
  SNP_DSB.GR=SNP_RAD51_Adastra.GR
}
#SNP_DSB.GR=pATM#pATM

#SNP_DSB.GR_w<-SNP_DSB.GR%>% anchor_center()%>%mutate(width = 100) 


#list of motif
pfm<-getMatrixSet(x=db(JASPAR2024),
                  opts=list(species=9606))


library(dplyr)
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
# SNP1.GR=GRanges("chr17",IRanges(31095150),REF="C",ALT="A",)
# 
# SNP2.GR=GRanges("chr7",IRanges(117559592,117559594),REF="CTT",ALT="-",rsid="rs113993960")
# 
# SNP3.GR=GRanges("chr11",IRanges(31664397),REF="C",ALT="A")

# SNPs.GR=c(SNP1.GR,SNP2.GR,SNP3.GR)
DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)

window=25
SNP_DSB.GR_W<-SNP_DSB.GR%>% anchor_center()%>%mutate(width = window) 

#get Allele specific sequence on the defined widow
seqSNPs.seq=getAllelesSeq(SNP_DSB.GR,window=window,genome=BSgenome.Hsapiens.UCSC.hg19)
seqSNPs.seqRef=DNAStringSet(seqSNPs.seq[[1]])
seqSNPs.seqAlt=DNAStringSet(seqSNPs.seq[[2]])
names(seqSNPs.seqRef)<-SNP_DSB.GR$ID
names(seqSNPs.seqAlt)<-SNP_DSB.GR$ID

############################# 

#pm=  5e-04
pm=  5e-03 # p value cutoff for motifMatchr

# get motif scores for each allele
Annotated_DSB_SNP.seqRef_pos<-matchMotifs(pfm,seqSNPs.seqRef,out="position",ranges=SNP_DSB.GR_W,p.cutoff = pm ) #ajouter le range diminue beaucoup le temps de calcul
Annotated_DSB_SNP.seqAlt_pos<-matchMotifs(pfm,seqSNPs.seqAlt,out="position",ranges=SNP_DSB.GR_W,p.cutoff = pm )
#?matchMotifs


U_Annotated_DSB_SNP.seqRef_pos<-unlist(Annotated_DSB_SNP.seqRef_pos)
U_Annotated_DSB_SNP.seqRef_pos$motifs<-names(U_Annotated_DSB_SNP.seqRef_pos)

score_table_ref<-join_overlap_left(SNP_DSB.GR,U_Annotated_DSB_SNP.seqRef_pos)
score_table_ref<-score_table_ref%>%as_tibble()

score_table_ref$score[is.na(score_table_ref$score)]=0
score_table_ref<-score_table_ref%>%dplyr::rename(motif_score_ref=score)



U_Annotated_DSB_SNP.seqAlt_pos<-unlist(Annotated_DSB_SNP.seqAlt_pos)
U_Annotated_DSB_SNP.seqAlt_pos$motifs<-names(U_Annotated_DSB_SNP.seqAlt_pos)

score_table_alt<-join_overlap_left(SNP_DSB.GR,U_Annotated_DSB_SNP.seqAlt_pos)# select the motifs that only overlaps dsbSNP in side the 25bk window arround the SNP
score_table_alt<-score_table_alt%>%as_tibble()

score_table_alt$score[is.na(score_table_alt$score)]=0
score_table_alt<-score_table_alt%>%dplyr::rename(motif_score_alt=score)

score_table_final<-left_join(score_table_ref,score_table_alt)#%>%mutate(motif_delta_scores_alt_minus_ref=motif_score_alt-motif_score_ref)
score_table_final$motif_score_alt[is.na(score_table_final$motif_score_alt)]=0
score_table_final$motif_score_ref[is.na(score_table_final$motif_score_ref)]=0

score_table_final$motifs[is.na(score_table_final$motifs)]="noMotif"

score_table_final<-score_table_final%>%mutate(motif_delta_scores_alt_minus_ref=motif_score_alt-motif_score_ref)


#get protein name corresponding to motif

nam<-unique(score_table_final$motifs[score_table_final$motifs!="noMotif"])

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



score_table_final<-left_join(score_table_final,id_motifs)

score_table_final$motif_tf[is.na(score_table_final$motif_tf)]="noMotif"
score_table_final$motif_score_ref=round(score_table_final$motif_score_ref,3)
score_table_final$motif_score_alt=round(score_table_final$motif_score_alt,3)
score_table_final$motif_delta_scores_alt_minus_ref=round(score_table_final$motif_delta_scores_alt_minus_ref,3)
score_table_final$prot=mode
colnames(score_table_final)
score_table_final<-score_table_final%>%select("seqnames",   "start",     "end", "width", "strand", "ID",      "ref",   "alt"   ,   "ES","motifs","motif_score_ref","motif_score_alt","motif_delta_scores_alt_minus_ref","prot","motif_tf")
if (mode==list_mode[1]){
df_final_motifs<-score_table_final

}else{
  
df_final_motifs<-rbind(df_final_motifs,score_table_final)
  
}


}
#score_table_final_motifs=score_table_final_motifs[score_table_final_motifs$motif_tf==motif,]    
score_table_final_motifs=df_final_motifs
mean(score_table_final_motifs$motif_delta_scores_alt_minus_ref)
mean(score_table_final_motifs$ES  )




motif="all" #compute for all motifs
t=0
score_table_final_motifs<-score_table_final_motifs[!(score_table_final_motifs$ES<t & score_table_final_motifs$ES > -t),]
score_table_final_motifs<-score_table_final_motifs[score_table_final_motifs$ES!=t,]
score_table_final_motifs$class="NA"
score_table_final_motifs$class[score_table_final_motifs$ES>t]=paste0("ES > ",t)
score_table_final_motifs$class[score_table_final_motifs$ES< -t]=paste0("ES < ",-t)
#wt=wilcox.test(score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES>t] , score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES< -t])
wt=""
for (prot in unique(score_table_final_motifs$prot)){ #compute wilcoxon pval for every protein
  
w=wilcox.test(score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES>t&score_table_final_motifs$prot==prot] , score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES< -t&score_table_final_motifs$prot==prot])
line=paste0(prot," p=",format(w$p.value,,digits=3),"\n")
wt<-paste0(wt,line)
}
 

p0<-score_table_final_motifs %>% ggplot()+aes(x= class,y=motif_delta_scores_alt_minus_ref,fill=reorder(class,-motif_delta_scores_alt_minus_ref))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2) +
  ggtitle(paste0(motif," \np=",(wt)))+#+geom_jitter()
  theme(legend.position = "none")+xlab("")+ylab("Motifmatchr delta Score (Alt-Ref)")+coord_cartesian(ylim=c(-10,10))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(~prot)
print(p0)


p01<-score_table_final_motifs %>% ggplot()+aes(x= class,y=motif_delta_scores_alt_minus_ref,fill=reorder(class,-motif_delta_scores_alt_minus_ref))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2) +
  ggtitle(paste0(motif," \np=",(wt)))+#+geom_jitter()
  theme(legend.position = "none")+xlab("")+ylab("Motifmatchr delta Score (Alt-Ref)")+coord_cartesian(ylim=c(-5,5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(~prot)
print(p01)

p02<-score_table_final_motifs %>% ggplot()+aes(x= class,y=motif_delta_scores_alt_minus_ref,fill=reorder(class,-motif_delta_scores_alt_minus_ref))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2) +
  ggtitle(paste0(motif," \np=",(wt)))+#+geom_jitter()
  theme(legend.position = "none")+xlab("")+ylab("Motifmatchr delta Score (Alt-Ref)")+coord_cartesian(ylim=c(-1.5,1.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(~prot)
print(p02)

pdf(paste0("results/allelic_imbalance/TF_motif/3_prots_ES_",motif,"_motimatchrp_",pm,".pdf"),1.72,6.34)
print(p0)
print(p01)
print(p02)
dev.off()
pdf(paste0("results/allelic_imbalance/TF_motif/3_prots_ES_",motif,"_motimatchrp_5to5_",pm,".pdf"),1.72,6.34)
print(p01)

dev.off()


motif_list=c("USF2","BHLHE40","TFEC")
motif_list=c("SP1")
motif_list=c("SRF", "MEF2A")
motif_list=c("ORC2") 
motif_list=c("CTCF") 
score_table_final_motifs=df_final_motifs

dup<-score_table_final_motifs[which(duplicated(score_table_final_motifs$ID[score_table_final_motifs$prot=="pATM"])),]
counted<-dup%>%group_by(ID)%>%mutate(val=1)%>%summarise(sum(val))
max(counted$`sum(val)`)
counted$ID[counted$`sum(val)`==max(counted$`sum(val)`)]



#compute for specific motifs

for (motif in motif_list){
  score_table_final_motifs=df_final_motifs
  
  score_table_final_motifs=score_table_final_motifs[score_table_final_motifs$motif_tf==motif,]              
  mean(score_table_final_motifs$motif_delta_scores_alt_minus_ref)
  mean(score_table_final_motifs$ES  )
  
  
#  wt=wilcox.test(score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES>t] , score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES< -t])
  
  
  
  t=0
  score_table_final_motifs<-score_table_final_motifs[!(score_table_final_motifs$ES<t & score_table_final_motifs$ES > -t),]
  score_table_final_motifs<-score_table_final_motifs[score_table_final_motifs$ES!=t,]
  score_table_final_motifs$class="NA"
  score_table_final_motifs$class[score_table_final_motifs$ES>t]=paste0("ES > ",t)
  score_table_final_motifs$class[score_table_final_motifs$ES< -t]=paste0("ES < ",-t)
  
  
  score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES > ",t)] =paste(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES > ",t)]," \n n=",length(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES > ",t)] ) )
  score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES < ",-t)] =paste(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES < ",-t)]," \n n=",length(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES < ",-t)] ) )
  #compute wilcoxon pval for every protein
   wt=""
  for (prot in unique(score_table_final_motifs$prot)){
    
    w=wilcox.test(score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES>t&score_table_final_motifs$prot==prot] , score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES< -t&score_table_final_motifs$prot==prot])
    line=paste0(prot," p=",format(w$p.value,,digits=3),"; ",
                "n ES >",t,"=",
                length(score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES>t&score_table_final_motifs$prot==prot]), "; ",
                "n ES <",-t,"=",
                length(score_table_final_motifs$motif_delta_scores_alt_minus_ref[score_table_final_motifs$ES< -t&score_table_final_motifs$prot==prot]),"\n"
                
    )
    wt<-paste0(wt,line)
  }
  
  score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES > ",t)] =paste(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES > ",t)]," \n n=",length(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES > ",t)] ) )
  score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES < ",-t)] =paste(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES < ",-t)]," \n n=",length(score_table_final_motifs$class[score_table_final_motifs$class==paste0("ES < ",-t)] ) )
  
  #plot at several scales
  p1<-score_table_final_motifs %>% ggplot()+aes(x= class,y=motif_delta_scores_alt_minus_ref,fill=reorder(class,-motif_delta_scores_alt_minus_ref))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2) +
    ggtitle(paste0(motif," \n",(wt)))+#+geom_jitter()
    theme(legend.position = "none")+xlab("")+ylab("Motifmatchr delta Score (Alt-Ref)")+coord_cartesian(ylim=c(-10,10))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(~prot)
  print(p1)
  
  p2<-score_table_final_motifs %>% ggplot()+aes(x= class,y=motif_delta_scores_alt_minus_ref,fill=reorder(class,-motif_delta_scores_alt_minus_ref))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2) +
    ggtitle(paste0(motif," \n",(wt)))+#+geom_jitter()
    theme(legend.position = "none")+xlab("")+ylab("Motifmatchr delta Score (Alt-Ref)")+coord_cartesian(ylim=c(-2,2))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(~prot)
  print(p2)
  p3<-score_table_final_motifs %>% ggplot()+aes(x= class,y=motif_delta_scores_alt_minus_ref,fill=reorder(class,-motif_delta_scores_alt_minus_ref))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2) +
    ggtitle(paste0(motif," \n",(wt)))+#+geom_jitter()
    theme(legend.position = "none")+xlab("")+ylab("Motifmatchr delta Score (Alt-Ref)")+coord_cartesian(ylim=c(-2.5,0.5))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+facet_wrap(~prot)
  print(p3)
  pdf(paste0("results/allelic_imbalance/TF_motif/",paste0(list_mode,collapse = "_"),"_ES_",motif,"delta_motimatchrp_",pm,".pdf"),1.72*length(list_mode),6.34)
  print(p1)
  print(p2)
  print(p3)
  
  dev.off()
  

}





