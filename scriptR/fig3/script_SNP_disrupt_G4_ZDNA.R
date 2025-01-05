# Raphael Mourad
# 15/06/2023


# This R script is used to do test if dsbSNPs affect G4s or Z-DNAs.


# SETUP DIRECTORY
setwd("/media/sauber/Elements/DSBsnp_project_seb/")
#install.packages("keras")
library(keras)
#reticulate::use_condaenv("DeepG4")
reticulate::use_condaenv("DeepG4")
#devtools::install_github("morphos30/DeepG4")

library(reticulate)
#use_virtualenv("r-tensorflow")
library(tensorflow)
library(keras)
#reticulate::use_condaenv("DeepG4")
#import("scipy")
# LOAD LIBRARIES
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(data.table)
library(vioplot)
library(rtracklayer)
library(keras)

library(DeepG4)

library(tidyverse)
library(ggplot2)
# LOAD INHOUSE LIBRARIES
source("scriptR/comp_NonBDNA_SNP.R")
source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")


# SETUP HG19 GENOME
Chr.V=c(paste0("chr",1:22),"chrX")
seqinfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)[Chr.V]


# PARAMETER
mod="all" # "all" "RAD51" "BRCA2" "53BP1"

###
# LOAD SNPS

load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
SNP_DSB.GR
SNP_DSB.GR[SNP_DSB.GR$ES<0]
load(file="results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR

# blacklist<-read.table("/media/ElissarDisk/SEB/blacklist_hg19.bed",sep='\t')
# blacklist<-GRanges(paste0(blacklist$V1),
#                    IRanges(blacklist$V2 ,blacklist$V3))
# SNP_DSB_blacklist<-findOverlaps(SNP_DSB.GR,blacklist)
# SNP_DSB_blacklist<-SNP_DSB.GR[SNP_DSB.GR %outside% blacklist,]
################deja triés par Raph
###
# Look if dsbSNPs disrupt G4 motifs
# 

# The main reason is that there is only 24 dsbSNPs that locate inside a G4
G4.GR=import.bed("data/nonBDDNA/g-quadruplex_forming_repeats_hg19.bed")
G4.GR=GRanges(seqnames(G4.GR),ranges(G4.GR))
start(G4.GR)=start(G4.GR)-1
olG4=findOverlaps(SNP_DSB.GR,G4.GR)
SNP_DSB_G4.GR=SNP_DSB.GR[queryHits(olG4)]
SNP_DSB_G4.GR=unique(SNP_DSB_G4.GR)

# Compute delta G4Hunter score (Not significant)
compg4=comp_G4Hunter_SNP(SNP_DSB.GR,window=50)
deltaScoreG4=compg4$G4scoreAlt-compg4$G4scoreRef
t.test(deltaScoreG4[SNP_DSB.GR$ES>0],deltaScoreG4[SNP_DSB.GR$ES<0])
length(deltaScoreG4)
SNP_DSB.GR$dG4h_ref<-compg4$G4scoreRef
SNP_DSB.GR$dG4h_Alt<-compg4$G4scoreAlt

SNP_DSB.GR$deltaScoreG4hunter=deltaScoreG4
SNP_DSB_G4h.GR<-SNP_DSB.GR

SNP_DSB_G4h.GR$G4h_strand[SNP_DSB.GR$dG4h_ref<0|SNP_DSB.GR$dG4h_Alt<0]<-"-"
SNP_DSB_G4h.GR$G4h_strand[SNP_DSB.GR$dG4h_ref>0|SNP_DSB.GR$dG4h_Alt>0]<-"+"
SNP_DSB_G4h.GR$dG4h_ref=abs(SNP_DSB_G4h.GR$dG4h_ref)
SNP_DSB_G4h.GR$dG4h_Alt=abs(SNP_DSB_G4h.GR$dG4h_Alt)
SNP_DSB_G4h.GR$deltaScoreG4hunter=round(SNP_DSB_G4h.GR$dG4h_Alt-SNP_DSB_G4h.GR$dG4h_ref,3)

SNP_DSB_G4h.GR$dG4h_ref=round(SNP_DSB_G4h.GR$dG4h_ref,3)
SNP_DSB_G4h.GR$dG4h_Alt=round(SNP_DSB_G4h.GR$dG4h_Alt,3)



path_save=("results/allelic_imbalance/G4hunter/")
dir.create(path_save)
save(SNP_DSB_G4h.GR,file=paste0(path_save,"SNP_DSB_G4h.GR"))


########################################
ESt=0
SNP_DSB_G4h.GR$class="NA"
SNP_DSB_G4h.GR$class[SNP_DSB_G4h.GR$ES>ESt]<-paste0("ES > ",ESt,"\n (n=",length(SNP_DSB_G4h.GR$ES[SNP_DSB_G4h.GR$ES>ESt]),")")
SNP_DSB_G4h.GR$class[SNP_DSB_G4h.GR$ES< -ESt]<-paste0("ES < ",-ESt,"\n (n=",length(SNP_DSB_G4h.GR$ES[SNP_DSB_G4h.GR$ES< -ESt]),")")

wil<-wilcox.test(SNP_DSB_G4h.GR$deltaScoreG4hunter[SNP_DSB_G4h.GR$class==paste0("ES > ",ESt,"\n (n=",length(SNP_DSB_G4h.GR$ES[SNP_DSB_G4h.GR$ES>ESt]),")")],
                 SNP_DSB_G4h.GR$deltaScoreG4hunter[SNP_DSB_G4h.GR$class==paste0("ES < ",-ESt,"\n (n=",length(SNP_DSB_G4h.GR$ES[SNP_DSB_G4h.GR$ES< -ESt]),")")])

p<-SNP_DSB_G4h.GR%>%as.tibble()%>%filter(class!="NA")%>%ggplot()+aes(x = reorder(class,ES),y=deltaScoreG4hunter,fill= reorder(class,-ES))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2)+coord_cartesian(ylim=c(-0.03, 0.03))+theme(legend.position="none")+
  coord_cartesian(ylim=c(-0.025, 0.025))+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  ggtitle(paste0("Delta (Alt-Ref) of\nG4hunter score\n(Effect Size above ",ESt,"\nand under -",ESt,")\n","wilcox p=",format(wil$p.value,digits=3)))+
  ylab("G4hunter score (Alt-Ref)")
p

SNP_DSB_G4h.GR%>%as.tibble()%>%filter(class!="NA")%>%ggplot()+aes(x = reorder(class,ES),y=abs(deltaScoreG4hunter),fill= reorder(class,-ES))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2)+theme(legend.position="none")+
 # coord_cartesian(ylim=c(-0.025, 0.025))+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  ggtitle(paste0("Delta (Alt-Ref) of\nG4hunter score\n(Effect Size above ",ESt,"\nand under -",ESt,")\n","wilcox p=",format(wil$p.value,digits=3)))+
  ylab("G4hunter score (Alt-Ref)")
p
pdf("results/allelic_imbalance/DeepG4/SNP_DSB_G4h.pdf",2,4)
print(p)
dev.off()


SNP_DSB_G4h.GR%>%as.tibble()%>%ggplot()+aes(x = abs(ES),y=abs(deltaScoreG4hunter))+geom_point()

###########################################
ESt=2.2
t.test(deltaScoreG4[SNP_DSB.GR$ES> ESt],deltaScoreG4[SNP_DSB.GR$ES< -ESt])
wilcox.test(deltaScoreG4[SNP_DSB.GR$ES> ESt],deltaScoreG4[SNP_DSB.GR$ES< -ESt])
mean(deltaScoreG4[SNP_DSB.GR$ES> ESt])
mean(deltaScoreG4[SNP_DSB.GR$ES< -ESt])

vioplot(deltaScoreG4[SNP_DSB.GR$ES> ESt],deltaScoreG4[SNP_DSB.GR$ES< -ESt])




###############################################################
# Do the same with DeepG4
# TODO Sebastien
#######fasta for alt, fasta for ref
window=201
region.GRr=resize(SNP_DSB.GR,window,fix="center")
region.seq=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, names=seqnames(region.GRr), 
                               start=start(region.GRr), end=end(region.GRr)))

SNPposrel=ceiling((window+1)/2)
region.seqRef=region.seq
substring(region.seqRef,SNPposrel,SNPposrel)=as.character(region.GRr$ref)
region.seqAlt=region.seq
substring(region.seqAlt,SNPposrel,SNPposrel)=as.character(region.GRr$alt)

resRef <- DeepG4(region.seqRef)
resAlt<-  DeepG4(region.seqAlt)

names(region.seqRef)<- paste0((region.GRr$ID),"_",(seqnames(region.GRr)))
names(region.seqAlt)<- paste0((region.GRr$ID),"_",(seqnames(region.GRr)))
write(region.seqRef,file="fasta_DSB_SNP_seq_w201_ref.fasta")
write(region.seqAlt,file="fasta_DSB_SNP_seq_w201_Alt.fasta")


resRef
resAlt

SNP_DSB.GR$resRef<- resRef
SNP_DSB.GR$resAlt<- resAlt
SNP_DSB.GR$DeepG4Comp<- round(SNP_DSB.GR$resAlt-SNP_DSB.GR$resRef,3)
SNP_DSB.GR$resRef<- round(SNP_DSB.GR$resRef,3)
SNP_DSB.GR$resAlt<- round(SNP_DSB.GR$resAlt,3)
SNP_DSB_deepG4.GR<-SNP_DSB.GR%>%as_tibble()

path_save="results/allelic_imbalance/DeepG4/"
dir.create(path_save)
write_tsv(SNP_DSB_deepG4.GR,paste0(path_save,"SNP_DSB_DeepG4.tsv"))

t.test(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES>1],SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES<1])





ESt=2
t.test(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES> ESt],SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES< -ESt])
wil<-wilcox.test(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES> ESt],SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES< -ESt])
wil$p.value
mean(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES> ESt])
mean(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES< -ESt])

vioplot(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES> ESt],SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES< -ESt])

df1<-tibble()
#nrow(df1)=length(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES> ESt])
df2<-data.frame()


df1= as_tibble(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES> ESt])
df1$name=paste0("ES > ",ESt,"\n(n= ",(length(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES> ESt])),")")
names(df1)=c("Delta DeepG4","names")

df2= as_tibble(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES< -ESt])
df2$name<- paste0("ES < -",ESt,"\n(n= ",length(SNP_DSB.GR$DeepG4Comp[SNP_DSB.GR$ES< -ESt]),")")
names(df2)=c("Delta DeepG4","names")

df3<-rbind(df1,df2)
p<-df3%>% ggplot()+aes(y =`Delta DeepG4`,x=reorder(`names`,`Delta DeepG4`),fill=reorder(`names`,-`Delta DeepG4`))+geom_violin()+geom_boxplot(width=0.2, color="black", alpha=0.2) +
  ggtitle(paste0("Delta (Alt-Ref) of\nDeepG4 score\n(Effect Size above ",ESt,"\nand under -",ESt,")\n","wilcox p=",format(wil$p.value,digits=3)))+xlab("")+theme(legend.position="none")+
  coord_cartesian(ylim=c(-0.010, 0.010))+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ ylab("DeepG4 (Alt-Ref)")
print(p)

pdf("results/allelic_imbalance/DeepG4/SNP_DSB_DeepG4.pdf",2,4)
print(p)
dev.off()



SNP_DSB.GR%>%as_tibble()%>% #filter(ES)%>% 
  ggplot()+aes(x =log(ES),y= log(DeepG4Comp))+geom_point()+geom_smooth()

SNP_DSB.GR%>%as_tibble()%>%#filter(ES)%>% 
  ggplot()+aes(x =ES,y= resAlt)+geom_point()+geom_smooth()


SNP_DSB.GR%>%as_tibble()%>%
  ggplot()+aes(x =ES,y=DeepG4Comp )+geom_point()+geom_smooth()





