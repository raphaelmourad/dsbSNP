
##################################################################

#This script is made to compute dsbSNp effect on TF motif score
###################################################################
library("Seurat")
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

DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
SNP_DSB.GR<-SNP_DSB.GR%>% anchor_center()%>%mutate(width = 100) 

JASPAR2024 <- JASPAR2024()
JASPAR2024=db(JASPAR2024)
#list of motif
pfm<-getMatrixSet(x=JASPAR2024,
                  opts=list(species=9606)
)

 
 heads(motif)


Annotated_DSB_SNP<-matchMotifs(pfm,SNP_DSB.GR,genome='hg19',out = "position",bg = "genome")
head((motif.matrix@data))
Annotated_DSB_SNP<-unlist(Annotated_DSB_SNP)
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_GR.RData")
Annotated_DSB_SNP.ti<-Annotated_DSB_SNP%>%anchor_center()%>%mutate(width=1)%>%as_tibble()%>%
  mutate(key=paste0(seqnames,"_",start))%>%mutate(motif=names(Annotated_DSB_SNP))%>%select(c(-score,-strand))%>%unique()%>%
  mutate(value=1)%>% pivot_wider(names_from = motif,values_from = value,values_fn = length,values_fill = 0)

SNP_DSB.ti<-SNP_DSB.GR%>%anchor_center()%>%mutate(width=1)%>%as_tibble()%>%mutate(key=paste0(seqnames,"_",start))

SNP_DSB.ti$key %in%Annotated_DSB_SNP.ti$key



Annotated_DSB_SNP.ti$key %in%SNP_DSB.ti$key
final_Annotated_DSB_SNP.ti<-na.omit(left_join(SNP_DSB.ti,Annotated_DSB_SNP.ti))
final_Annotated_DSB_SNP.GR<-final_Annotated_DSB_SNP.ti%>%as_granges()



mat_dataExpe_DSB_all<-values(final_Annotated_DSB_SNP.GR)
row.names(mat_dataExpe_DSB_all)<-mat_dataExpe_DSB_all$ID
mat_dataExpe_DSB_all=as.data.frame(mat_dataExpe_DSB_all[,!colnames(mat_dataExpe_DSB_all)%in%c("seqnames","start","end","key","width")])
colnames(mat_dataExpe_DSB_all)
mat_dataExpe_DSB_all<-mat_dataExpe_DSB_all[,-c(1:44)]#on retire les dsb prot
colnames(mat_dataExpe_DSB_all)

#mat_dataExpe_DSB_all<-mat_dataExpe_DSB_all[colSums(mat_dataExpe_DSB_all)>40]
mat_dataExpe_DSB_all<-mat_dataExpe_DSB_all[colSums(mat_dataExpe_DSB_all)!=0]

cormat_DSB_all=as.matrix(proxy::simil(mat_dataExpe_DSB_all, by_rows = F, method = "Jaccard",diag = T, upper = T))
#?p
col<- colorRampPalette(c("white", "red"))(256)
heatmap.2(cormat_DSB_all,scale = "none",trace = "none", density.info = "none",col=col)

max(colSums(mat_dataExpe_DSB_all))
colnames(mat_dataExpe_DSB_all[colSums(mat_dataExpe_DSB_all)==0])
sum(mat_dataExpe_DSB_all$MA0071.1)

heatmap.2(as.matrix(mat_dataExpe_DSB_all),scale = "row",trace = "none", density.info = "none",Rowv = F,Colv = F, col=col)


##########################################################################################

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(universalmotif))
#BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm3")
library(tidyr)
library("memes")
#data("example_chip_summits", package = "memes")
library(readr)


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
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")
#resize dsb SNP
window=25
SNP_DSB.GR_W<-SNP_DSB.GR%>% anchor_center()%>%mutate(width = window) 

seqSNPs.seq=getAllelesSeq(SNP_DSB.GR,window=window,genome=BSgenome.Hsapiens.UCSC.hg19)
seqSNPs.seqRef=DNAStringSet(seqSNPs.seq[[1]])
seqSNPs.seqAlt=DNAStringSet(seqSNPs.seq[[2]])
names(seqSNPs.seqRef)<-SNP_DSB.GR$ID
names(seqSNPs.seqAlt)<-SNP_DSB.GR$ID
pfm<-getMatrixSet(x=JASPAR2024,
                  
                  opts=list(species=9606)
)


Annotated_DSB_SNP.seqRef<-matchMotifs(pfm,seqSNPs.seqRef,ranges=SNP_DSB.GR_W)
Annotated_DSB_SNP.seqAlt<-matchMotifs(pfm,seqSNPs.seqAlt,ranges=SNP_DSB.GR_W)
head(motifMatches(Annotated_DSB_SNP.seqRef))

Annotated_DSB_SNP.seqRef_score<-matchMotifs(pfm,seqSNPs.seqRef,out="score",ranges=SNP_DSB.GR_W)
Annotated_DSB_SNP.seqAlt_score<-matchMotifs(pfm,seqSNPs.seqAlt,out="score",ranges=SNP_DSB.GR_W)
test <-(motifScores(Annotated_DSB_SNP.seqRef_score))%>%as.data.frame()
head(motifCounts(Annotated_DSB_SNP.seqRef_score))
head(motifMatches(Annotated_DSB_SNP.seqRef_score))


#get motif score at dsbSNP position
Annotated_DSB_SNP.seqRef_pos<-matchMotifs(pfm,seqSNPs.seqRef,out="position",ranges=SNP_DSB.GR_W)
Annotated_DSB_SNP.seqAlt_pos<-matchMotifs(pfm,seqSNPs.seqAlt,out="position",ranges=SNP_DSB.GR_W)

U_Annotated_DSB_SNP.seqRef_pos<-unlist(Annotated_DSB_SNP.seqRef_pos)
U_Annotated_DSB_SNP.seqRef_pos$motifs<-names(U_Annotated_DSB_SNP.seqRef_pos)

data_f_ref<-join_overlap_left(SNP_DSB.GR,U_Annotated_DSB_SNP.seqRef_pos)
data_f_ref<-data_f_ref%>%as_tibble()

data_f_ref$score[is.na(data_f_ref$score)]=0
data_f_ref<-data_f_ref%>%rename(motif_score_ref=score)



U_Annotated_DSB_SNP.seqAlt_pos<-unlist(Annotated_DSB_SNP.seqAlt_pos)
U_Annotated_DSB_SNP.seqAlt_pos$motifs<-names(U_Annotated_DSB_SNP.seqAlt_pos)



data_f_alt<-join_overlap_left(SNP_DSB.GR,U_Annotated_DSB_SNP.seqAlt_pos)
data_f_alt<-data_f_alt%>%as_tibble()

data_f_alt$score[is.na(data_f_alt$score)]=0
data_f_alt<-data_f_alt%>%rename(motif_score_alt=score)

data_f_final<-left_join(data_f_ref,data_f_alt)#%>%mutate(motif_delta_scores_alt_minus_ref=motif_score_alt-motif_score_ref)
data_f_final$motif_score_alt[is.na(data_f_final$motif_score_alt)]=0
data_f_final$motif_score_ref[is.na(data_f_final$motif_score_ref)]=0

data_f_final$motifs[is.na(data_f_final$motifs)]="noMotif"

data_f_final<-data_f_final%>%mutate(motif_delta_scores_alt_minus_ref=motif_score_alt-motif_score_ref)


#get motif names
nam<-unique(data_f_final$motifs[data_f_final$motifs!="noMotif"])
#nam<-nam[1:5]
namlist<-mclapply(nam,mc.cores=8,
                  function(m){ getMatrixByID (JASPAR2024,ID=m )})

motifs<-c()
motif_tf<-c()
for (i in (1:length(namlist))){
  print(i)
motif_tf[i]<-name(namlist[[i]])
motifs[i]<-ID(namlist[[i]])
}

#data_f_final$m
data_f_final<-left_join(data_f_final,id_motifs)

data_f_final$motif_tf[is.na(data_f_final$motif_tf)]="noMotif"
data_f_final$motif_score_ref=round(data_f_final$motif_score_ref,3)
data_f_final$motif_score_alt=round(data_f_final$motif_score_alt,3)
data_f_final$motif_delta_scores_alt_minus_ref=round(data_f_final$motif_delta_scores_alt_minus_ref,3)

data_f_final_expand<-data_f_final


data_f_final_collapsed<-data_f_final%>%group_by(seqnames,start,end,ID,motifs)%>%summarise(summary=paste(motifs,"|",motif_tf,"|",motif_score_ref,"|",motif_score_alt,"|",motif_delta_scores_alt_minus_ref,collapse = ""))
data_f_final_collapsed<-left_join(data_f_final_collapsed,SNP_DSB.GR%>%as_tibble())
data_f_final_ultra_collapsed<-data_f_final_collapsed%>%group_by(ID)%>%summarise(paste(summary_motif=summary,collapse=" _ "))
data_f_final_ultra_collapsed<-left_join(data_f_final_ultra_collapsed,SNP_DSB.GR%>%as_tibble())


data_f_final_collapsedV2<-data_f_final%>%group_by(seqnames,start,end,ID,motifs,motif_score_ref,motif_score_alt,motif_delta_scores_alt_minus_ref)%>%summarise(summary=paste(motifs,"_",motif_tf))
data_f_final_ultra_collapsedV2.1<-data_f_final_collapsedV2%>%group_by(ID)%>%summarise(paste(summary_motif=summary,collapse=";"))
data_f_final_ultra_collapsedV2.2<-data_f_final_collapsedV2%>%group_by(ID)%>%summarise(motif_score_ref=paste(motif_score_ref,collapse=";"))
data_f_final_ultra_collapsedV2.3<-data_f_final_collapsedV2%>%group_by(ID)%>%summarise(motif_score_alt=paste(motif_score_alt,collapse=";"))
data_f_final_ultra_collapsedV2.4<-data_f_final_collapsedV2%>%group_by(ID)%>%summarise(motif_delta_scores_alt_minus_ref=paste(motif_delta_scores_alt_minus_ref,collapse=";"))

data_f_final_ultra_collapsedV2<-left_join(data_f_final_ultra_collapsedV2.1,data_f_final_ultra_collapsedV2.2)%>%left_join(data_f_final_ultra_collapsedV2.3)%>%left_join(data_f_final_ultra_collapsedV2.4)

data_f_final_ultra_collapsedV2<-left_join(data_f_final_ultra_collapsedV2,SNP_DSB.GR%>%as_tibble())

#data_f_final_collapsedV2<-data_f_final%>%group_by(seqnames,start,end,ID,motifs,motif_score_ref,motif_score_alt,motif_delta_scores_alt_minus_ref)%>%summarise(summary=paste(motifs,"_",motif_tf))




path_save="results/allelic_imbalance/DSB_SNP_motifs/"
dir.create(path_save)

write_tsv(data_f_final_expand,paste0(path_save,"DSB_SNP_motif_score_expanded.tsv"))
write_tsv(data_f_final_collapsed,paste0(path_save,"DSB_SNP_motif_score_collapsed.tsv"))
write_tsv(data_f_final_ultra_collapsed,paste0(path_save,"DSB_SNP_motif_score_ultra_collapsed.tsv"))
write_tsv(data_f_final_ultra_collapsedV2,paste0(path_save,"DSB_SNP_motif_score_ultra_collapsedV2.tsv"))



data_f_final_expand$ES
data_f_final_expand$motif_delta_scores_alt_minus_ref
library(ggplot2)
data_f_final_expand_plot<-data_f_final_expand



data_f_final_expand_plot%>% filter(motif_score_ref>0)%>% filter(motif_score_alt<0)%>%
  ggplot(aes(x=(as.numeric(ES)),y=(as.numeric(motif_score_ref))))+geom_point()



data_f_final_expand_plot_miror=data_f_final_expand_plot
data_f_final_expand_plot_miror$motif_score_alt[data_f_final_expand_plot_miror$ES<0]=-data_f_final_expand_plot_miror$motif_score_alt[data_f_final_expand_plot_miror$ES<0]
data_f_final_expand_plot_miror$motif_score_ref[data_f_final_expand_plot_miror$ES<0]=-data_f_final_expand_plot_miror$motif_score_ref[data_f_final_expand_plot_miror$ES<0]
data_f_final_expand_plot_miror$motif_delta_scores_alt_minus_ref=data_f_final_expand_plot_miror$motif_score_ref-data_f_final_expand_plot_miror$motif_score_alt
data_f_final_expand_plot_miror$ES[data_f_final_expand_plot_miror$ES<0]=-data_f_final_expand_plot_miror$ES[data_f_final_expand_plot_miror$ES<0]

#data_f_final_expand_plot_miror$motif_delta_scores_alt_minus_ref[data_f_final_expand_plot_miror$ES<0]=-data_f_final_expand_plot_miror$motif_delta_scores_alt_minus_ref[data_f_final_expand_plot_miror$ES<0]




data_f_final_expand_plot_miror%>% #filter(motif_delta_scores_alt_minus_ref!=0)%>% 
  ggplot(aes(x=(as.numeric(ES)),y=(as.numeric(motif_delta_scores_alt_minus_ref))))+geom_point()

data_f_final_expand_plot_miror%>% filter(motif_score_alt!=0)%>% 
  ggplot(aes(x=(as.numeric(ES)),y=(as.numeric(motif_score_alt))))+geom_point()

data_f_final_expand_plot_miror%>% filter(motif_score_ref!=0)%>% 
  ggplot(aes(x=(as.numeric(ES)),y=(as.numeric(motif_score_ref))))+geom_point()




data_f_final_expand_plot<-data_f_final_expand



thr=0
data_f_final_expand_plot$over<-"oui"
data_f_final_expand_plot$over[data_f_final_expand_plot$ES>=thr]<-paste0(" dsbSNP\neffect >  ",thr)
data_f_final_expand_plot$over[data_f_final_expand_plot$ES<=-thr]<-paste0("dsbSNP\neffect < ",thr)

w<-wilcox.test(abs(data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0(" dsbSNP\neffect >  ",thr)]),abs(data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0("dsbSNP\neffect < ",thr)]))
pabs<-w$p.value
w<-wilcox.test((data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0(" dsbSNP\neffect >  ",thr)]),(data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0("dsbSNP\neffect < ",thr)]))
p<-w$p.value

data_f_final_expand_plot%>%filter(over!="oui")%>%
#  filter( motif_delta_scores_alt_minus_ref!=0)%>%
  ggplot( aes(x=over, y=abs(motif_delta_scores_alt_minus_ref), fill=over)) +
  geom_violin()  +geom_boxplot() + #coord_cartesian(ylim=c(0,0.01))+ 
  ggtitle(paste0("Aboslute value of Delta (Alt-Ref) of \nJASPAR2020 score (Effect Size above ",thr,"\nand under -",thr,")\nWilcoxon pval=",pabs))+ 
  ylab("ABS(Alt-Ref)")

p<-data_f_final_expand_plot%>%filter(over!="oui")%>%
  #filter( motif_delta_scores_alt_minus_ref!=0)%>%
  ggplot( aes(x=reorder(over,motif_delta_scores_alt_minus_ref), y=(motif_delta_scores_alt_minus_ref), fill=over)) +
  geom_violin()  +geom_boxplot(width=0.2, color="black", alpha=0.2) + 
  theme(legend.position = "none")+coord_cartesian(ylim=c(-10,10))+ 
  ggtitle(paste0("All motifs Delta score \np=",format(p,digits=3)))+ 
 # ggtitle(paste0(" Delta (Alt-Ref) of \nJASPAR2020 score (Effect Size above ",thr,"\nand under -",thr,")\nWilcoxon pval=",p))+ 
  ylab("Delta score(Alt-Ref)")+xlab("dsbSNP\neffect")
p
path_save="results/allelic_imbalance/DSB_SNP_motifs/"

pdf(paste0(path_save,"boxplot_score_ES_under_and_over_0.pdf"),2,4)
p
dev.off()

##############################################only!=0
data_f_final_expand_plot<-data_f_final_expand
thr=0
data_f_final_expand_plot$over<-"oui"
data_f_final_expand_plot$over[data_f_final_expand_plot$ES>=thr]<-paste0("over ",thr)
data_f_final_expand_plot$over[data_f_final_expand_plot$ES<=-thr]<-paste0("under -",thr)
data_f_final_expand_plot<-data_f_final_expand_plot%>%filter( motif_delta_scores_alt_minus_ref!=0)
w<-wilcox.test(abs(data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0("over ",thr)]),abs(data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0("under -",thr)]))
pabs<-w$p.value
w<-wilcox.test((data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0("over ",thr)]),(data_f_final_expand_plot$motif_delta_scores_alt_minus_ref [data_f_final_expand_plot$over==paste0("under -",thr)]))
p<-w$p.value

data_f_final_expand_plot%>%filter(over!="oui")%>%
  #  filter( motif_delta_scores_alt_minus_ref!=0)%>%
  ggplot( aes(x=over, y=abs(motif_delta_scores_alt_minus_ref), fill=over)) +
  geom_violin()  +geom_boxplot() + #coord_cartesian(ylim=c(0,1))+ 
  ggtitle(paste0("Aboslute value of Delta (Alt-Ref!=0) of \nJASPAR2024 score (Effect Size above ",thr,"\nand under -",thr,")\nWilcoxon pval=",pabs))+ 
  ylab("ABS(Alt-Ref)")

data_f_final_expand_plot%>%filter(over!="oui")%>%
  #  filter( motif_delta_scores_alt_minus_ref!=0)%>%
  ggplot( aes(x=over, y=(motif_delta_scores_alt_minus_ref), fill=over)) +
  geom_violin()  +geom_boxplot() +# coord_cartesian(ylim=c(-1,1))+ 
  ggtitle(paste0(" Delta (Alt-Ref!=0) of \nJASPAR2024 score (Effect Size above ",thr,"\nand under -",thr,")\nWilcoxon pval=",p))+ 
  ylab("(Alt-Ref)")
#################################################
################################################################################
data_f_final_expand_protein<-data_f_final_expand
colnames(data_f_final_expand_protein)

data_f_final_expand_protein<-data_f_final_expand_protein%>%pivot_longer(cols=colnames(data_f_final_expand_protein)[10:59], names_to="protein")

data_f_final_expand_protein=data_f_final_expand_protein[data_f_final_expand_protein$value==1,]

unique(data_f_final_expand_protein)%>%
  ggplot(aes(x=motif_tf,y=protein,fill=motif_delta_scores_alt_minus_ref))+geom_tile()

unique(data_f_final_expand_protein$motif_tf)

score_t=0

data_f_final_expand_protein_plot<-data_f_final_expand_protein[13:15]
data_f_final_expand_protein_plot$motif_delta_scores_alt_minus_ref=as.numeric(data_f_final_expand_protein_plot$motif_delta_scores_alt_minus_ref)
f_data_f_final_expand_protein_plot<-unique(data_f_final_expand_protein_plot)%>%group_by(protein,motif_tf)%>%summarise(mean_delta_score=mean(motif_delta_scores_alt_minus_ref))#
#f_data_f_final_expand_protein_plot<-f_data_f_final_expand_protein_plot[f_data_f_final_expand_protein_plot$mean_delta_score>=score_t|f_data_f_final_expand_protein_plot$mean_delta_score<= -score_t,]
f_data_f_final_expand_protein_plot<-f_data_f_final_expand_protein_plot%>%pivot_wider(names_from=motif_tf,values_from = (mean_delta_score),values_fill =as.double(0))

 m_data_f_final_expand_protein_plot<-f_data_f_final_expand_protein_plot[-1]%>%as.matrix()
 rownames(m_data_f_final_expand_protein_plot)<-f_data_f_final_expand_protein_plot$protein
 
 col<- colorRampPalette(c("blue","white", "red"))(256)
 
heatmap.2(as.matrix(m_data_f_final_expand_protein_plot))
heatmap.2(m_data_f_final_expand_protein_plot,scale = "none",trace = "none", density.info = "none",col=col)



mean(unique(data_f_final_expand_protein_plot)$motif_delta_scores_alt_minus_ref[unique(data_f_final_expand_protein_plot)$motif_tf=="CTCF"])

####prot DSB only
score_t=0
DSB_prot=c("X53BP1","X53BP1_S","MDC1","pATM","BRCA2","BRCA1","NBN","Lig4","ATM","RAD51","RPA","SETX")

data_f_final_expand_protein_plot<-data_f_final_expand_protein[13:15]
data_f_final_expand_protein_plot<-data_f_final_expand_protein_plot[data_f_final_expand_protein_plot$protein %in% DSB_prot,]
data_f_final_expand_protein_plot$motif_delta_scores_alt_minus_ref=as.numeric(data_f_final_expand_protein_plot$motif_delta_scores_alt_minus_ref)
f_data_f_final_expand_protein_plot<-unique(data_f_final_expand_protein_plot)%>%group_by(protein,motif_tf)%>%summarise(mean_delta_score=mean(motif_delta_scores_alt_minus_ref))#
#f_data_f_final_expand_protein_plot<-f_data_f_final_expand_protein_plot[f_data_f_final_expand_protein_plot$mean_delta_score>=score_t|f_data_f_final_expand_protein_plot$mean_delta_score<= -score_t,]
f_data_f_final_expand_protein_plot<-f_data_f_final_expand_protein_plot%>%pivot_wider(names_from=motif_tf,values_from = (mean_delta_score),values_fill =as.double(0))

m_data_f_final_expand_protein_plot<-f_data_f_final_expand_protein_plot[-1]%>%as.matrix()

rownames(m_data_f_final_expand_protein_plot)<-f_data_f_final_expand_protein_plot$protein

col<- colorRampPalette(c("blue","white", "red"))(256)

#heatmap.2(as.matrix(m_data_f_final_expand_protein_plot))
heatmap.2(m_data_f_final_expand_protein_plot,scale = "none",trace = "none", density.info = "none",col=col)
path_save="results/allelic_imbalance/DSB_SNP_motifs/"

pdf(paste0(path_save,"heatmap_motif_DSB_prot_6000DSB_SNP.pdf"),60,10)
heatmap.2(m_data_f_final_expand_protein_plot,scale = "none",trace = "none", density.info = "none",col=col)
dev.off()
pdf(paste0(path_save,"heatmap_motif_DSB_prot_6000_DSB_SNP_2.pdf"),10,60)
heatmap.2(t(m_data_f_final_expand_protein_plot),scale = "none",trace = "none", density.info = "none",col=col)
dev.off()

##########thres of mean
score_t=2
m_data_f_final_expand_protein_plot<-f_data_f_final_expand_protein_plot[-1]#%>%as.matrix()

max(colMeans(m_data_f_final_expand_protein_plot))
min(colMeans(m_data_f_final_expand_protein_plot))
colMeans(m_data_f_final_expand_protein_plot)
m_data_f_final_expand_protein_plot<-m_data_f_final_expand_protein_plot[colMeans(m_data_f_final_expand_protein_plot)>score_t|colMeans(m_data_f_final_expand_protein_plot)< -score_t]
m_data_f_final_expand_protein_plot<-as.matrix(m_data_f_final_expand_protein_plot)


rownames(m_data_f_final_expand_protein_plot)<-f_data_f_final_expand_protein_plot$protein

col<- colorRampPalette(c("blue","white", "red"))(256)

#heatmap.2(as.matrix(m_data_f_final_expand_protein_plot))
heatmap.2(m_data_f_final_expand_protein_plot,scale = "none",trace = "none", density.info = "none",col=col)


###############################################################################


data_f_final_ultra_collapsedV3<-data_f_final_ultra_collapsedV2%>%pivot_longer(cols=colnames(data_f_final_expand)[10:49], names_to="protein")
data_f_final_ultra_collapsedV3<-data_f_final_ultra_collapsedV3%>%rename( summary_motifs=`paste(summary_motif = summary, collapse = ";")`)
data_f_final_ultra_collapsedV3<-data_f_final_ultra_collapsedV3[data_f_final_ultra_collapsedV3$value==1,]

data_f_final_ultra_collapsedV3<-data_f_final_ultra_collapsedV3%>%group_by(ID,summary_motifs,motif_score_ref,motif_score_alt,motif_delta_scores_alt_minus_ref,seqnames,start,end,width,strand,ref,alt)%>%
  summarise(protein_summary=paste(protein,collapse = ";"))
write_tsv(data_f_final_ultra_collapsedV3,paste0(path_save,"DSB_SNP_motif_score_ultra_collapsedV3.tsv"))






