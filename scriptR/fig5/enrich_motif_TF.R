###########
#Sébastien AUBER
#10/07/2025
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
library(pbmcapply)

library(universalmotif)
library(TFBSTools)
library(plyranges)
DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)
load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_GR.RData")

########################## Create non redondant human pfm##################################



jaspar_dir <- "/media/sauber/Elements/DSBsnp_project_seb/data/JASPAR_non_redond/"
jaspar_files <- list.files(jaspar_dir, pattern = "\\.jaspar$", full.names = TRUE)

# Conversion universalmotif -> PFMatrix
umotif_to_PFMatrix <- function(um) {
  PFMatrix(
    ID = um@name,
    name = um@name,
    profileMatrix = um@motif,   # PFM brute
    strand = "+",
    bg = c(A=0.25, C=0.25, G=0.25, T=0.25)
  )
}

#read motifs
pfms <- list()
for (f in jaspar_files) {
  ums <- read_jaspar(f)
  
  if (inherits(ums, "universalmotif")) ums <- list(ums)
  
  for (um in ums) {
    pfms[[length(pfms)+1]] <- umotif_to_PFMatrix(um)
  }
}

# Create PFMatrixList
pfm_list <- do.call(PFMatrixList, pfms)

#add names
names(pfm_list) <- sapply(pfms, function(x) x@name)

# check
pfm_list
names(pfm_list)


# load Jaspar db
jaspar_db <- JASPAR2024()
jaspar_db <- db(jaspar_db)
JASPAR2024 <- JASPAR2024()
db(JASPAR2024)#for name calling


ids <- sapply(pfm_list, function(x) x@ID)

# EXTRACT ONLY HUMAN IDs
human_pfms <- getMatrixSet(jaspar_db, opts = list(species = 9606, matrix_id = ids))
human_ids <- names(human_pfms)  # ou sapply(human_pfms, function(x) x@ID)
pfm_list_human <- pfm_list[names(pfm_list) %in% human_ids]




DSB_path="/media/sauber/Elements/DSBsnp_project_seb"
setwd(DSB_path)

list_mode=c("pATM","RAD51","53BP1")
#list_mode=c("pATM")



pfm_list<-pfm_list_human

#############################compute enrich###############################################

SNP_DSB.GR<-SNP_DSB.GR%>% anchor_center()%>%mutate(width = 100) #take 100 bp arround a SNP



#######find motif at SNP regions
Annotated_DSB_SNP<-matchMotifs(pfm_list,SNP_DSB.GR,genome='hg19',out = "position",bg = "genome")
#head((motif.matrix@data))
Annotated_DSB_SNP<-unlist(Annotated_DSB_SNP)


Annotated_DSB_SNP$motif=names(Annotated_DSB_SNP)

########do the same for control SNPs
load(file="results/allelic_imbalance/SNP_DSB/SNP_ctrl_dist500kb_GR.RData")
SNP_ctrl.GR<-sort(sample(SNP_ctrl.GR,length(SNP_DSB.GR)))
SNP_ctrl.GR<-SNP_ctrl.GR%>% anchor_center()%>%mutate(width = 100) 

Annotated_ctrl_SNP<-matchMotifs(pfm_list,SNP_ctrl.GR,genome='hg19',out = "position",bg = "genome")
#head((motif.matrix@data))
Annotated_ctrl_SNP<-unlist(Annotated_ctrl_SNP)
Annotated_ctrl_SNP$motif=names(Annotated_ctrl_SNP)

load(file="results/allelic_imbalance/SNP_DSB/SNP_DSB_ES_corrected_GR.RData")



SNPID<-SNP_DSB.GR[,(-2: -length(colnames(mcols(SNP_DSB.GR))))]
final_DSB<-join_overlap_inner(SNPID,Annotated_DSB_SNP) #ADD the SNP IDs



SNPID<-SNP_ctrl.GR%>% anchor_center()%>%mutate(width = 1) 
final_ctrl<-join_overlap_inner(SNPID,Annotated_ctrl_SNP)#ADD the SNP IDs

length(final_ctrl)

tablectrl<-table(final_ctrl$rs,final_ctrl$motif) #get motif table per SNP
tableSNP<-table(final_DSB$ID,final_DSB$motif)  #get motif table per SNP

#############take common motifs
SNP_ctrl<-colSums(tablectrl[,which(colnames(tablectrl)%in%colnames(tableSNP))])
SNP_ctrl<-(as.data.frame(SNP_ctrl))
SNP_ctrl$motif<-rownames(SNP_ctrl)

SNP_DSB<-colSums(tableSNP[,which(colnames(tableSNP)%in%colnames(tablectrl))])

SNP_DSB<-(as.data.frame(SNP_DSB))
SNP_DSB$motif<-rownames(SNP_DSB)

############compute fold change
matEnrichTF<-merge.data.frame(SNP_DSB,SNP_ctrl,by="motif")
matEnrichTF$FC=matEnrichTF$SNP_DSB/matEnrichTF$SNP_ctrl


#########get protein name corresponding to motif
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
matEnrichTF<-left_join(matEnrichTF,id_motifs)

matEnrichTF$len<-length(SNP_ctrl.GR)

#compute fisher exact test p.values
for (i in 1:nrow(matEnrichTF)){
dff<-data.frame(matEnrichTF[i,2],matEnrichTF[i,3])
dff[2,]<-c(length(SNP_ctrl.GR),length(SNP_ctrl.GR))

library(data.table) 
setDT(dff)

p<-fisher.test(dff)


matEnrichTF[i,7]<-p$p.value
}

matEnrichTF<-matEnrichTF%>%rename(p=V7)
matEnrichTF$fdr<- p.adjust(matEnrichTF$p)#get adjust pvalues

write.table(matEnrichTF,"results/allelic_imbalance/TF_motif/matriceTF.txt")

###plot
file_enrich_TF_motif="results/allelic_imbalance/TF_motif/plot_enrich_TF_motif.pdf"
pdf(file_enrich_TF_motif,7,5)
plot(matEnrichTF$SNP_DSB,matEnrichTF$FC,xlab="Number of dsbSNPs in TF motifs",ylab="Enrichment (Fold-change)",
     log="y",pch=3)
abline(0,0,col="red")
labelTF=(matEnrichTF$motif_tf)
labelTF[matEnrichTF$FC<10&matEnrichTF$SNP_DSB<500]=""
text(matEnrichTF$SNP_DSB,matEnrichTF$FC,labels=labelTF)
dev.off()


