###########
#This script is made to compute  TF motif enrichment at dsbSNPs

############



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

library(plyranges)
JASPAR2024<-JASPAR2024()
JASPAR2024<-db(JASPAR2024)
DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_GR.RData")
SNP_DSB.GR<-SNP_DSB.GR%>% anchor_center()%>%mutate(width = 100) 
pfm<-getMatrixSet(x=JASPAR2024,
                  opts=list(species=9606)
)



Annotated_DSB_SNP<-matchMotifs(pfm,SNP_DSB.GR,genome='hg19',out = "position",bg = "genome")
#head((motif.matrix@data))
Annotated_DSB_SNP<-unlist(Annotated_DSB_SNP)


Annotated_DSB_SNP$motif=names(Annotated_DSB_SNP)


load(file="results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR<-sort(sample(SNP_ctrl.GR,length(SNP_DSB.GR)))
SNP_ctrl.GR<-SNP_ctrl.GR%>% anchor_center()%>%mutate(width = 100) 

Annotated_ctrl_SNP<-matchMotifs(pfm,SNP_ctrl.GR,genome='hg19',out = "position",bg = "genome")
#head((motif.matrix@data))
Annotated_ctrl_SNP<-unlist(Annotated_ctrl_SNP)
Annotated_ctrl_SNP$motif=names(Annotated_ctrl_SNP)

load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")



SNPID<-SNP_DSB.GR[,(-2: -length(colnames(mcols(SNP_DSB.GR))))]
# hits<-queryHits(findOverlaps(Annotated_DSB_SNP,SNP_DSB.GR))
# Annotated_DSB_SNP[hits]
final_DSB<-join_overlap_inner(SNPID,Annotated_DSB_SNP)



SNPID<-SNP_ctrl.GR%>% anchor_center()%>%mutate(width = 1) 
# hits<-queryHits(findOverlaps(Annotated_ctrl_SNP,SNP_ctrl.GR))
# Annotated_ctrl_SNP[hits]
final_ctrl<-join_overlap_inner(SNPID,Annotated_ctrl_SNP)

length(final_ctrl)

tablectrl<-table(final_ctrl$rs,final_ctrl$motif)
tableSNP<-table(final_DSB$ID,final_DSB$motif)


SNP_ctrl<-colSums(tablectrl[,which(colnames(tablectrl)%in%colnames(tableSNP))])
SNP_ctrl<-(as.data.frame(SNP_ctrl))
SNP_ctrl$motif<-rownames(SNP_ctrl)

SNP_DSB<-colSums(tableSNP[,which(colnames(tableSNP)%in%colnames(tablectrl))])

SNP_DSB<-(as.data.frame(SNP_DSB))
SNP_DSB$motif<-rownames(SNP_DSB)

matEnrichTF<-merge.data.frame(SNP_DSB,SNP_ctrl,by="motif")
matEnrichTF$FC=matEnrichTF$SNP_DSB/matEnrichTF$SNP_ctrl


library(pbmcapply)
nam<-unique(matEnrichTF$motif)
#nam<-nam[1:5]
namlist<-mclapply(nam,mc.cores=8,
                  function(m){ getMatrixByID (JASPAR2024,ID=m )})
motif<-c()
motif_tf<-c()
for (i in (1:length(namlist))){
  print(i)
  motif_tf[i]<-name(namlist[[i]])
  motif[i]<-ID(namlist[[i]])
}

id_motifs<-data.frame(motif,motif_tf)
#alors_final$m
matEnrichTF<-left_join(matEnrichTF,id_motifs)

matEnrichTF$len<-length(SNP_ctrl.GR)

for (i in 1:nrow(matEnrichTF)){
dff<-data.frame(matEnrichTF[i,2],matEnrichTF[i,3])
dff[2,]<-c(length(SNP_ctrl.GR),length(SNP_ctrl.GR))

library(data.table) 
setDT(dff)

p<-fisher.test(dff)


matEnrichTF[i,7]<-p$p.value
}
matEnrichTF<-matEnrichTF%>%rename(p=V7)
matEnrichTF$fdr<- p.adjust(matEnrichTF$p)

write.table(matEnrichTF,"results/allelic_imbalance/TF_motif/matriceTF.txt")

file_enrich_TF_motif="results/allelic_imbalance/TF_motif/plot_enrich_TF_motif.pdf"
pdf(file_enrich_TF_motif,7,5)
plot(matEnrichTF$SNP_DSB,matEnrichTF$FC,xlab="Number of dsbSNPs in TF motifs",ylab="Enrichment (Fold-change)",
     log="y",pch=3)
abline(0,0,col="red")
labelTF=(matEnrichTF$motif_tf)
labelTF[matEnrichTF$FC<10&matEnrichTF$SNP_DSB<500]=""
text(matEnrichTF$SNP_DSB,matEnrichTF$FC,labels=labelTF)
dev.off()


