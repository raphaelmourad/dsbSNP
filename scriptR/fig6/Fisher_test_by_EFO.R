#09/10/2025
#Sébastien AUBER

#this script compares enrichement of dsbSNP vs ctrlSNP regions using Fisher’s test on harmonized p-values."

require(tidyr)
require(dplyr)
require("parallel")
library(tidyverse)
library(Matrix)
library(MASS)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(liftOver)
library(rtracklayer)
library(rlist)
library(poolr)
library(pbmcapply)
args = commandArgs(trailingOnly=TRUE)
thresP=1e-08
setwd ("/mnt/DATA/MCD/eq_legube/SEB/GWAS/")

#define SNP type
SNP_type_list=strsplit(read_lines(args[1])," ")

#define number of pvalues under10-6
n_pval=args[2]

#remoove study with issues and define EFO list
listStudy<- read_tsv(paste0("data/list_study_more_than_",n_pval,"_pvalue_under10Eminus6.tsv"))
EFO_list=unique(listStudy$EFO[listStudy$EFO!="EFO_0007937.h.tsv.gz"&listStudy$EFO!="EFO_0005105.h.tsv.gz"&listStudy$EFO!= "EFO_0007056.h.tsv.gz"&listStudy$EFO!= "EFO_0004526.h.tsv.gz"])


print(EFO_list)
print(paste0("n EFO"=length(EFO_list)))

#define list of control
ctrl_list= c("RDM","CTRL","PROM","ENH","SE")

for (ctrl_type in ctrl_list){
  for (SNP_type in SNP_type_list) {
    df<-pbmclapply(mc.cores=99,EFO_list ,function(EFO){
      
      print(EFO)  
      
      diseases= listStudy$Disease[listStudy$EFO==EFO]
      diseases <- gsub('/','_',diseases)
      
      
      
      print(SNP_type)
      
      
      setwd (paste0("/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_noTAG_multi_",n_pval,"pval_under_10minus6/",SNP_type,"/"))
     
      index=c()  
      listPval=list()
      listES=list()
      listBeta=list()
      path =(paste0("/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_noTAG_multi_",n_pval,"pval_under_10minus6/",SNP_type,"/"))
      
      for ( i in (1:length(diseases))) {
        data_list_DSB<- data.frame()
        data_list_CTRL<- data.frame()
        print(diseases[i])
        print(i)
        if(file.exists(paste0(path,"harmonised_tab_lists_All_rsID2/listPval_sup1000_harmonised_All_rsID_",diseases[i],"_SNP"))){
          load(file=paste0("harmonised_tab_lists_All_rsID2/listPval_sup1000_harmonised_All_rsID_",diseases[i],"_SNP"))
          
          print(length( data_list_DSB$p_SNP))
          
          load(file=paste0("harmonised_tab_lists_All_rsID2/listPval_sup1000_harmonised_All_rsID_",diseases[i],"_",ctrl_type))
          
          
          
          if(ctrl_type=="CTRL"){
            (listPval[[i]]=list(p_dsb=data_list_DSB$p_SNP,p_ctrl=data_list_CTRL$p_CTRL))
          }
          if(ctrl_type=="PROM"){
            (listPval[[i]]=list(p_dsb=data_list_DSB$p_SNP,p_ctrl=data_list_prom$p_PROM))
          }
          if(ctrl_type=="ENH"){
            (listPval[[i]]=list(p_dsb=data_list_DSB$p_SNP,p_ctrl=data_list_ENH$p_ENH))
          }
          if(ctrl_type=="SE"){
            (listPval[[i]]=list(p_dsb=data_list_DSB$p_SNP,p_ctrl=data_list_SE$p_SE))
          }
          if(ctrl_type=="RDM"){
            (listPval[[i]]=list(p_dsb=data_list_DSB$p_SNP,p_ctrl=data_list_RDM$p_RDM))
          }
        }
      }
      
      print(paste0(ctrl_type," list ok"))  
     
    #compute contagency table across replicate for one EFO  
      l_ftab<-lapply(c(1:length(diseases)),function(i){
        listPvali=listPval[[i]]
        
        listPvali$p_ctrl<- listPvali$p_ctrl%>%as.numeric()
        listPvali$p_dsb<- listPvali$p_dsb%>%as.numeric()
        
        f_tab<-rbind(c( length(listPvali$p_dsb[listPvali$p_dsb<thresP]), length(listPvali$p_ctrl[listPvali$p_ctrl<thresP])/10),#10 time more control at the origin 
                     c( length(listPvali$p_dsb[listPvali$p_dsb>=thresP]), length(listPvali$p_ctrl[listPvali$p_ctrl>=thresP])/10))
        
        f_tab        
      })
      
      tab_f_test<- Reduce("+",l_ftab)  
      f_test=fisher.test(tab_f_test+1, alternative = "greater") #compute fisher test with a pseudo count
      f_pval=f_test$p.value
      f_OR = f_test$estimate
      c(f_pval,f_OR)
      
      
      name=diseases[1]%>%str_split("_GCST")%>%dplyr::first()%>%dplyr::first()
      n_p=tab_f_test[1,1]
      n_d<-length(diseases)
      mean_n<-n_p/n_d
      
      tot_p<-tab_f_test[1,1]+tab_f_test[2,1]
      mean_tot_p<-  tot_p/n_d
      
      #ctrl
      n_c=tab_f_test[1,2]
      n_d<-length(diseases)
      mean_n_c<-n_c/n_d 
      tot_c<-tab_f_test[1,2]+tab_f_test[2,2]
      mean_c_tot<-  tot_c/n_d
      
      
      list(p=f_pval,OR=as.numeric(f_OR),name=name,ctrl_type=ctrl_type,EFO=EFO,n_mean=mean_n,n_ctrl_mean=mean_n_c,mean_tot=mean_tot_p,mean_tot_control=mean_c_tot)
    })#%>%bind_rows()
    print(df)
    df<-df%>%bind_rows()
    
    matMedPval=df
    setwd (paste0("/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_noTAG_multi_",n_pval,"pval_under_10minus6/",SNP_type,"/"))
    #create directories
    dir.create("Analysis_pval_fisher_greater_fused_dup_pseudo_count_at_1_with_tot_RDM")
    dir.create(paste0("Analysis_pval_fisher_greater_fused_dup_pseudo_count_at_1_with_tot_RDM/",ctrl_type))
    dir.create(paste0("Analysis_pval_fisher_greater_fused_dup_pseudo_count_at_1_with_tot_RDM/",ctrl_type,"/thresP"))
    write.csv(matMedPval,paste0("Analysis_pval_fisher_greater_fused_dup_pseudo_count_at_1_with_tot_RDM/",ctrl_type,"/thresP","/matMedPval_sup1000_harmonised_6000SNP_noTAG_",ctrl_type,"_thresP_",thresP,".csv"),row.names = F)
  
    #create corrected matrix
    matCorr<- matMedPval[matMedPval$p!= "NA in allele or location",]
    nrow(matCorr)
    
    matCorr <- matCorr[matCorr$p!= "n DSB under threshold limit",]
    matCorr <- matCorr[matCorr$p!= "n CTRL under threshold limit",]
    matCorr <- matCorr[matCorr$p!= "check not ok",]
    
    
   
    matCorr$p<- as.numeric(as.character(matCorr$p))
    matCorr$fdr_with_all_1314_trait<- p.adjust(as.numeric(matCorr$p), method = "fdr",n=nrow(matMedPval))
    matCorr$fdr_filterwith_nrows_trait<-p.adjust(as.numeric(matCorr$p), method = "fdr",n=nrow(matCorr))
    #write table
    write.table(matCorr,paste0("Analysis_pval_fisher_greater_fused_dup_pseudo_count_at_1_with_tot_RDM/",ctrl_type,"/thresP","/matMedPvalcorrected_sup1000_harmonised_filtered_with_fdr_6000SNP_",ctrl_type,"_thresP_",thresP,"_noTAG.tsv"),sep="\t",row.names = F)
    
  }
}


