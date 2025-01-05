library(dplyr)
library(tidyr)
expe_type="ChIP-seq"
date_expe="6000SNPs_Nov2023"
# expe_type="CUT-Tag"
# expe_type="ATAC-seq"
# expe_type="DRIP-seq"
setwd("/media/sauber/Elements/DSBsnp_project_seb/")

source("scriptR/helper_functions.R")
source("scriptR/misc_functions.R")

suppressPackageStartupMessages(library(plyranges))
# Threshold
PvalT=1
MinCoverage=1

DSB=T

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
for(i in 1:nrow(dataExpe)){
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
    for(j in length(dup)) 
      if(abs(ESi[dup[[j]]])>abs(ESi[dup[[j]]-1])){
        ESi=ESi[-(dup[[j]]-1)]
        expei=expei[-(dup[[j]]-1)]
      }else{
        ESi=ESi[-(dup[[j]])]
        expei=expei[-(dup[[j]])]
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
if(expe_type=="ChIP-seq"){
  mat_SNP[,"pATM"]=mat_SNP[,"pATM"]+mat_SNP[,"pATM_siCtrl"]
  mat_SNP=mat_SNP[,colnames(mat_SNP)!="pATM_siCtrl"]
  mat_SNP[,"Ser2P"]=mat_SNP[,"Ser2P"]+mat_SNP[,"Pol2S2P"]+mat_SNP[,"S2P"]
  mat_SNP=mat_SNP[,colnames(mat_SNP)!="Pol2S2P"]
  mat_SNP=mat_SNP[,colnames(mat_SNP)!="S2P"]
}
colSums(mat_SNP)

# Number of SNP by experience
sort(colMeans(mat_SNP),decreasing=T)
dataExpe.GR=GRanges(dataExpe[,"chr"],IRanges(dataExpe[,"end"],dataExpe[,"end"]),
                    ID=dataExpe[,"id"],ref=dataExpe[,"ref"],alt=dataExpe[,"alt"],ES=dataExpe[,"comb_es"])
values(dataExpe.GR)=cbind(values(dataExpe.GR),mat_SNP)

###
# REMOVE SUBTELOMERIC + CENTROMERIC REGIONS

dataExpe.GR=filterSNPCentroTelo(dataExpe.GR)



if(DSB==T){ 
  if(expe_type=="ChIP-seq"){

    # dsbSNPs
    repairProtExpes=c("53BP1","Lig4","NBN","pATM","RPA","SETX")
    mat_dataExpe=values(dataExpe.GR)
    mat_dataExpe_colsel=as.data.frame(mat_dataExpe[,colnames(mat_dataExpe)%in%repairProtExpes])
    mat_dataExpe_DSB=mat_dataExpe_colsel[rowSums(mat_dataExpe_colsel)!=0,]

    dataExpe.GR=dataExpe.GR[rowSums(mat_dataExpe_colsel)!=0]
    #save(SNP_DSB_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_GR.RData")

    # dsbSNPs from perturbed experiments (KD, drugs)
    #repairProtPertubExpes=c("pATM_DRB","pATM_siSCC1","pATM_curaxin","pATM_DRB")
    #mat_dataExpe_colselpertub=as.data.frame(mat_dataExpe[,colnames(mat_dataExpe)%in%repairProtPertubExpes])
    #SNP_DSB_pertub_Legube_DIVA.GR=dataExpe.GR[rowSums(mat_dataExpe_colselpertub)>0]
    #save(SNP_DSB_pertub_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DSB_pertub_Legube_DIVA_GR.RData")

    # non dsbSNPs
    #SNP_nonDSB_Legube_DIVA.GR=dataExpe.GR[rowSums(mat_dataExpe_colsel)==0]
    #save(SNP_nonDSB_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_nonDSB_Legube_DIVA_GR.RData")

  }else if(expe_type=="ATAC-seq"){

    SNP_ATAC_Legube_DIVA.GR=dataExpe.GR
    #save(SNP_ATAC_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_ATAC_Legube_DIVA_GR.RData")

  }else if(expe_type=="DRIP-seq"){

    SNP_DRIP_Legube_DIVA.GR=dataExpe.GR
    #save(SNP_DRIP_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_DRIP_Legube_DIVA_GR.RData")

  }else if(expe_type=="CUT-Tag"){

    repairProtExpes=c("pATM","RPA")
    mat_dataExpe=values(dataExpe.GR)
    mat_dataExpe_colsel=as.data.frame(mat_dataExpe[,colnames(mat_dataExpe)%in%repairProtExpes])

    SNP_CUTnTag_Legube_DIVA.GR=dataExpe.GR[rowSums(mat_dataExpe_colsel)>0]
    #save(SNP_CUTnTag_Legube_DIVA.GR,file="results/allelic_imbalance/SNP_DSB/SNP_CUTnTag_Legube_DIVA_GR.RData")

  }
}
# #masked AsiSi


sites= read.table("data/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv",header=T)

sites<-sites[1:400,]
sites.GR<- GRanges(sites[,1],IRanges(sites[,2],sites[,3]))


prot_list=list(
  "pATM"=2000,
  "53BP1"=1000000,
  "Lig4"=20000,
  "MDC1"=1000000
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

if(DSB==T){
  load("results/allelic_imbalance/SNP_DSB/SNP_DSB_Legube_DIVA_ES_corrected_GR.RData")
  hits<-findOverlaps(dataExpe.GR,SNP_DSB_Legube_DIVA.GR)
  dataExpe.GR<-dataExpe.GR[queryHits(hits)]
  }
matExpe<-dataExpe.GR%>%as_tibble()%>%dplyr::select(-c(seqnames,start,end,ID,ref,ES,width,strand,alt))

library(corrplot)
matExpe<-matExpe%>%dplyr::select(-c(pATM_siSCC1,pATM_DRB,pATM_curaxin))

matExpe[matExpe==0]<-NA

cor.info <-cor(matExpe,use="pairwise.complete.obs")


matExpe_count=matExpe
matExpe_count[matExpe_count!=0]=1
# on eneleve la ou on a pas assez de points
one_liste<-c()
two_liste<-c()
for(one in colnames(matExpe_count)){
  print(one)
  for (two in colnames(matExpe_count)){
    print(two)
    tabl<-matExpe_count%>%dplyr::select(c(any_of(one),any_of(two)))%>%table()
    if(one!=two){
    if(tabl["1","1"]<20){
      one_liste=c(one_liste,one)
      two_liste=c(two_liste,two)
    }
    }
  }
}
df_under_th<-data.frame(one_liste,two_liste)
for (i in 1:nrow(df_under_th)){
  
cor.info[df_under_th[i,1],df_under_th[i,2]]<-NA
}

cor.info_no_NA<-cor.info
cor.info_no_NA[is.na(cor.info_no_NA)]<-0
clust<-hclust(as.dist(1 -cor.info_no_NA),method = c("ward.D"))
ord<-clust$order
plot((clust))
#?dist
# ord<-order.dendrogram(as.dendrogram(clust))
# plot(as.dendrogram(clust))

# 
# cor.info[clust]
# cor.info[order(row.names(cor.info)[ord]),order(colnames(cor.info)[ord])]

data=cor.info
#on reordonne:
library(ggdendro)
# data <- factor(x = rownames(data) ,
#                     levels = rownames(data)[ord], 
#                     ordered = TRUE)
# 
# data$names2<-as.factor(data$names2)
# data$names2 <- factor(x = data$names2 ,
#                       levels = data$name[clustering_order], 
#                       ordered = TRUE)
# 
# 
# cor.info[(ord)]
library("matsbyname")
ordered_corinfo<-sort_rows_cols( cor.info,roworder = row.names(cor.info)[clust$order], colorder = colnames(cor.info)[clust$order])
#ordered_corinfo+g

# library("textshape")
# ordered_corinfo<-cluster_matrix(cor.info,method = "complete",dim="both")

# ordered_corinfo[ordered_corinfo>0.5]=0.75
# ordered_corinfo[ordered_corinfo< -0.5]=-0.75
# min(na.omit(ordered_corinfo))
corrplot(ordered_corinfo,  order="original", tl.cex=0.7, cl.cex=0.5, tl.col="black", addrect=6,is.corr = F,#col = COL2('RdBu', 8),
         na.label = "square", na.label.col = "grey")#,col.lim = c(-0.5,0.5))

name_l<-c()
for (i in 1:length(colnames(cor.info))){
# if(length(cor.info[,i][is.na(cor.info[,i])])>=length(colnames(cor.info))-20){
  if(length(cor.info[,i][is.na(cor.info[,i])])>=30){
    
name_l<-c(name_l,colnames(cor.info)[i])  
}
}
#length(cor.info[,15][is.na(cor.info[,15])])
cor.info2<-cor.info%>%as_tibble%>%
  dplyr::select(-c(any_of(name_l)))
rownames(cor.info2)=rownames(cor.info)
cor.info2<-cor.info2%>%as_tibble%>%
  filter((!row.names(cor.info2)%in%(name_l)))
row.names(cor.info2)=row.names(cor.info)[!row.names(cor.info)%in%(name_l)]

cor.info3<- as.matrix(cor.info2)
row.names(cor.info3)<-row.names(cor.info2)



cor.info3_no_NA<-cor.info3
cor.info3_no_NA[is.na(cor.info3_no_NA)]<-0
clust<-hclust(as.dist(1 -cor.info3_no_NA),method = c("ward.D"))
ord<-clust$order
plot((clust))

#?dist
library(ggplot2)
ordered_corinfo3<-sort_rows_cols( cor.info3,roworder = row.names(cor.info3)[clust$order], colorder = colnames(cor.info3)[clust$order])
ordered_corinfo4<-sort_rows_cols( cor.info3,roworder = row.names(cor.info3)[clust$order], colorder = colnames(cor.info3)[clust$order])
ordered_corinfo4<-ordered_corinfo3%>%as_tibble%>%dplyr::select(pATM,SETX,RPA,X53BP1,Lig4)%>%as.matrix()
row.names(ordered_corinfo4)<-row.names(ordered_corinfo3)
if(DSB==T){
  ordered_corinfo3[ordered_corinfo3>0.5]=0.5
   ordered_corinfo3[ordered_corinfo3< -0.5]=-0.5
   max(na.omit(ordered_corinfo3))
  corrplot(ordered_corinfo3,  order="original", tl.cex=0.7, cl.cex=0.5, tl.col="black", addrect=6,is.corr = F,col = COL2('RdBu', 100),
           na.label = "square", na.label.col = "grey",col.lim = c(-0.5,0.5))
  
}


p<-corrplot(cor.info3,  order="original", tl.cex=1.3, cl.cex=1, tl.col="black", addrect=6,is.corr = F,col = COL2('RdBu', 8),
            na.label = "square", na.label.col = "grey")


#ordered_corinfo3<-sort_rows_cols( cor.info3,roworder = row.names(cor.info3)[clust$order], colorder = colnames(cor.info3)[clust$order])
p<-corrplot(t(ordered_corinfo4),  order="original", tl.cex=1.3, cl.cex=1, tl.col="black", addrect=6,is.corr = T,col = COL2('RdBu', 8),
         na.label = "square", na.label.col = "grey")#+#,col.lim = c(-0.5,0.5))


pdf("All_warD_pearson.pdf",11,8)

corrplot(ordered_corinfo3,  order="original", tl.cex=1.3, cl.cex=1, tl.col="black", addrect=6,is.corr = F,col = COL2('RdBu', 8),
            na.label = "square", na.label.col = "grey")#+#,col.lim = c(-0.5,0.5))

dev.off()


matExpe_count%>%pivot_longer(any_of(colnames(matExpe)))

cor.info2=cor.info%>%as.matrix()
cor.info2[cor.info2>0.5]=0.5
cor.info2[cor.info2<-0.5]=-0.5
cor.info2[is.na(cor.info2)]=0
cor.info[is.na(cor.info)]=0
corrplot((cor.info2),  order="hclust", tl.cex=0.7, cl.cex=0.5, tl.col="black", addrect=2,is.corr = F,col = COL2('RdBu', 10), hclust.method = ("complete"))
corrplot(cor.info%>%as.matrix(),  order="hclust", tl.cex=0.7, cl.cex=0.5, tl.col="black", addrect=6,is.corr = F,col = COL2('RdBu', 10), hclust.method = "complete")
corrplot(cor.info_no_NA%>%as.matrix(),  order="hclust", tl.cex=0.7, cl.cex=0.5, tl.col="black", addrect=6,is.corr = F,col = COL2('RdBu', 10), hclust.method = "complete")


plot(clust)

library(gplots)
heatmap.2(
  as.matrix(cor.info), Rowv=as.dendrogram(hclust(as.dist(1-cor.info_no_NA))),
  Colv=as.dendrogram(hclust(as.dist(1-cor.info_no_NA))))
table(cor.info)

cor.info <- factor( cor.info,
                     levels = row.names(cor.info)[clust$order] )
matExpe%>%mutate(r_n=rownames(matExpe))%>%pivot_longer(any_of(colnames(matExpe)))



dataExpe.ti<-dataExpe.GR%>% as_tibble()


library(ggpmisc)
library(ggpubr)


######################""
##hexbin
########################

dataExpe.ti<-dataExpe.GR%>% as_tibble()


inspect_after_stat <- function(p, i = 1L) {
  ._env <- environment()
  .out <- NULL
  suppressMessages({
    trace(
      what = "ggplot_build.ggplot",
      tracer = substitute(assign(".out", data[[i]], envir = ._env), ._env),
      at = 19L,
      print = FALSE,
      where = asNamespace("ggplot2")
    )
  })
  ggplot_build(p)
  suppressMessages({
    untrace("ggplot_build.ggplot", where = asNamespace("ggplot2"))
  })
  .out
}

?scale_fill_gradient2()
?geom_hex
library(RColorBrewer)

my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
get_hex<-function(request){
  
  plot.ti<-dataExpe.ti%>%dplyr::select(c(any_of(request)))%>%filter(get(request[1])!=0)%>%filter(get(request[2])!=0)#%>%pivot_longer(c(any_of(request)))
  print(plot.ti)
  cort=cor.test(plot.ti%>%dplyr::select(any_of(request[1]))%>%as.vector()%>%unlist,plot.ti%>%dplyr::select(any_of(request[2]))%>%as.vector()%>%unlist,method = "pearson")
  # cort$estimate*cort$estimate
  print((nrow(plot.ti)/30)/2) #
  
  
  p<-plot.ti%>%ggplot(aes(x=get(request[1]),y=get(request[2])))+geom_hex(bins=30)+
    scale_fill_gradient2(low="white",mid = "red",high="red4",midpoint = 15)+
    
   # scale_fill_gradient(low="navyblue",high="red")+
   # scale_color_brewer(palette='Set1')+
    geom_smooth(method=lm)+
    ylab(request[2])+xlab(request[1])+
    theme(
      axis.title.x = element_text(size = 18*2),
      axis.text.x = element_text(size = 16*2),
      axis.text.y = element_text(size = 16*2),
      
      axis.title.y = element_text(size = 18*2),
      title = element_text(size = 30)
    )+ggtitle(paste0("R=",format(cort$estimate,digits=2), ",\np=",format(cort$p.value,digits=2)))+
  #  ggtitle(after_stat(count))+
    coord_cartesian(ylim = c(-2,2),xlim=c(-2,2))+xlab(paste0(request[1], " Effect Size"))+ylab(paste0(request[2], " Effect Size"))
  # stat_poly_line() +
  # stat_poly_eq() + 
  # stat_cor(size=13)#+ 
  #stat_regline_equation() 
  
  
  mid=max(inspect_after_stat(p)$count)/2
  #mid=median(inspect_after_stat(p1)$count)
  print(mid)
  p<-p+ scale_fill_gradient2(low="navy",mid = "cyan4",high="yellow",midpoint = mid)+ theme(legend.text=element_text(size=20))
  
  return(p)
}


p1=get_hex(c("pATM","Ser5P"))
p2=get_hex(c("pATM","Ser7P"))
p3=get_hex(c("pATM","CTCF"))
p4=get_hex(c("pATM","H2AZ"))
p5=get_hex(c("pATM","SUN1"))
p6=get_hex(c("pATM","SUN2"))
p7=get_hex(c("RPA","CTCF"))
p8=get_hex(c("pATM","T1P"))
p9=get_hex(c("pATM","H3K27Ac"))
p10=get_hex(c("pATM","H3K4me3"))
p11=get_hex(c("pATM","H3K4me1"))
p12=get_hex(c("RPA","CTCF"))
p13=get_hex(c("RPA","Ser2P"))
p14=get_hex(c("Lig4","CTCF"))
p15=get_hex(c("Lig4","SUN1"))

if (DSB==T){
  pdf("results/allelic_imbalance/hexbin_es_prot_dsb_only_square.pdf",5.5,5.5)
  }else{
pdf("results/allelic_imbalance/hexbin_es_prot_square.pdf",5.5,5.5)
}
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
print(p11)
print(p12)
print(p13)
print(p14)
print(p15)


dev.off()

