library(rtracklayer)
library(tibble)
library(gkmSVM)
library(BSgenome.Hsapiens.UCSC.hg19)
library("BRGenomics")

args = commandArgs(trailingOnly=TRUE)


dir.create("hg19_blacklisted")
dir.create("hg19_blacklisted_400bestAsiSI_masked")
dfname<-args[1]
datai<- read.table(dfname)
GR<- GRanges(datai[,1],IRanges(datai[,2],datai[,3]))


values(GR)<- datai[,c(-1,-2,-3)]  

hg19_blacklist<-read.csv("../blacklist_hg19.bed",sep="\t")
hg19_blacklist<- GRanges(hg19_blacklist[,1],IRanges(hg19_blacklist[,2],hg19_blacklist[,3]))
names(GR)<-datai[,6]

finalGR <- GR[!GR %over% hg19_blacklist,]

write.table(as_tibble(finalGR),file =  paste0("hg19_blacklisted/",args[2],"_hg19_blacklisted.bed"),sep="\t", row.names = F, col.names=F)
hg19=BSgenome.Hsapiens.UCSC.hg19
myseq=getSeq(hg19,finalGR )
#names(myseq)<- (paste0(names(finalGR),".",seqnames(finalGR)))


writeXStringSet(myseq,file= paste0("hg19_blacklisted/",args[2],"_hg19_blacklisted.fasta"))



##########mask 400 best clived asiSI sites

sites= read.table("../DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv",header=T)

sites<-sites[1:400,]
sites.GR<- GRanges(sites[,1],IRanges(sites[,2],sites[,3]))
sites.GR<-sites.GR%>% resize(width = as.numeric(args[3]),fix="center")


ultimeGR <- finalGR[!finalGR %over% sites.GR,]


write.table(as_tibble(ultimeGR),file =  paste0("hg19_blacklisted_400bestAsiSI_masked/",args[2],"hg19_blacklisted_400bestAsiSI_masked_",args[3],".bed"),sep="\t", row.names = F, col.names=F)
myseq=getSeq(hg19,ultimeGR )
#names(myseq)<- (paste0(names(finalGR),".",seqnames(finalGR)))


writeXStringSet(myseq,file= paste0("hg19_blacklisted_400bestAsiSI_masked/",args[2],"_hg19_blacklisted_400bestAsiSI_masked_",args[3],".fasta"))




