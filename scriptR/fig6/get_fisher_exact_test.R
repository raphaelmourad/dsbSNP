############################################################################
#this script allows to compute fisher exact test between interest SNP types and
#their multiple controls

#it requires the output of get_tab_lists_TAG.R
############################################################################


require(tidyr)
require(dplyr)
require(parallel)
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
args = commandArgs(trailingOnly = TRUE)



setwd ("/mnt/DATA/MCD/eq_legube/SEB/GWAS/")


#get the Type of SNPs to request (eg. dsb_no_eSNP)

SNP_type_list = strsplit(read_lines(args[1]), " ")

n_pval = args[2]#minimum number of pvals under 10e-6 for a GWAS to be included


###############################################################################
#get interest GWAS trait and diseases

listStudy <- read_tsv(paste0(
  "data/list_study_more_than_",
  n_pval,
  "_pvalue_under10Eminus6.tsv"
))
diseases = listStudy$Disease
diseases <- gsub('/', '_', diseases)

for (SNP_type in SNP_type_list) {
  print(SNP_type)
  
  
  ctrl_list = c("CTRL", "PROM", "ENH", "SE")#define the type of controls
  
  setwd (
    paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_TAG_multi_",
      n_pval,
      "pval_under_10minus6/",
      SNP_type,
      "/"
    )
  )
  
  
  
  #parse interest SNPs
  for (ctrl_type in ctrl_list) {
    listPval = list()
    
    path = (
      paste0(
        "/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_TAG_multi_",
        n_pval,
        "pval_under_10minus6/",
        SNP_type,
        "/"
      )
    )
    
    for (i in (1:length(diseases))) {
      data_list_interest_SNPs <- data.frame()
      data_list_CTRL <- data.frame()
      print(diseases[i])
      print(i)
      if (file.exists(
        paste0(
          path,
          "harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
          diseases[i],
          "_SNP"
        )
      )) {
        load(
          file = paste0(
            "harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
            diseases[i],
            "_SNP"
          )
        )
        
        print(length(data_list_interest_SNPs$p_SNP))
        
        load(
          file = paste0(
            "harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
            diseases[i],
            "_",
            ctrl_type
          )
        )
        
        
        if (ctrl_type == "CTRL") {
          (
            listPval[[i]] = list(
              p_interest = data_list_interest_SNPs$p_SNP,
              p_ctrl = data_list_CTRL$p_CTRL
            )
          )
        }
        if (ctrl_type == "PROM") {
          (
            listPval[[i]] = list(
              p_interest = data_list_interest_SNPs$p_SNP,
              p_ctrl = data_list_prom$p_PROM
            )
          )
        }
        if (ctrl_type == "ENH") {
          (
            listPval[[i]] = list(
              p_interest = data_list_interest_SNPs$p_SNP,
              p_ctrl = data_list_ENH$p_ENH
            )
          )
        }
        if (ctrl_type == "SE") {
          (
            listPval[[i]] = list(
              p_interest = data_list_interest_SNPs$p_SNP,
              p_ctrl = data_list_SE$p_SE
            )
          )
        }
        
      }
    }
    
    print(paste0(ctrl_type, " list ok"))
    ###############################################################################
    
    #################################################################################
    
    #perform fisher exact test using several threshold of GWAS association significance
    
    
    for (thresP in c(
      0.00000001,
      0.0000001,
      0.0000001,
      0.000001,
      0.00001,
      0.0001,
      0.001,
      0.01,
      0.05,
      0.1,
      1
    )) {
      matMedPval = matrix(NA, length(diseases), 4)
      #matMedPval[,c(3,5)]=1
      matP = NULL
      
      matMedPval <- pbmclapply(c(1:length(diseases)), mc.cores = 99, function(i) {
        print(ctrl_type)
        print(i)
        print((i / length(diseases)) * 100)
        print(((i / length(diseases)) / which(ctrl_list == ctrl_type)) *
                100)
        print(diseases[i])
        if (!("NA in allele or location" %in% listPval[[i]]$p_interest)) {
          if (!("check not ok" %in% listPval[[i]]$p_interest)) {
            if (!("harmonised file has a 0 size" %in% listPval[[i]]$p_interest)) {
              listPvali = listPval[[i]]
              print("pvali")
              listPvali$p_interest = as.numeric(listPvali$p_interest)
              print("asnum dsb")
              listPvali$p_ctrl = as.numeric(listPvali$p_ctrl)
              print("as.num ctrl")
              #thresP=0.001
              if (length(listPvali$p_interest < thresP) > 10) {
                if (length(listPvali$p_ctrl < thresP) > 10) {
                  print("under thresP ok")
                  matMedPval[i, 1] = poolr::fisher(listPvali$p_interest)$p
                  
                  print(length(listPvali$p_ctrl))
                  print(length(listPvali$p_interest))
                  
                  if (length(listPvali$p_ctrl) > length(listPvali$p_interest)) {
                    matMedPval[i, 2] = poolr::fisher(listPvali$p_ctrl[sample(1:length(listPvali$p_ctrl),
                                                                             length(listPvali$p_interest))])$p
                  } else{
                    matMedPval[i, 2] = poolr::fisher(listPvali$p_ctrl)$p
                    
                  }
                  
                  f_tab <- rbind(c(
                    length(listPvali$p_interest[listPvali$p_interest < thresP]),
                    length(listPvali$p_interest[listPvali$p_ctrl < thresP])
                  ), c(
                    length(listPvali$p_interest[listPvali$p_interest >= thresP]),
                    length(listPvali$p_interest[listPvali$p_ctrl >= thresP])
                  ))
                  f_test = fisher.test(f_tab, alternative = "greater")
                  
                  matMedPval[i, 3] = f_test$p.value
                  matMedPval[i, 4] = f_test$estimate
                  
                  
                } else{
                  print(paste0("under thresPctrl not  ok"))
                  
                  matMedPval[i, 1] = "n CTRL under threshold limit"
                  matMedPval[i, 2] = "n CTRL under threshold limit"
                  matMedPval[i, 3] = "n CTRL under threshold limit"
                  matMedPval[i, 4] = "n CTRL under threshold limit"
                  
                  
                  
                }
              } else{
                print(paste0("under thresP DSB not ok "))
                
                matMedPval[i, 1] = "n DSB under threshold limit"
                matMedPval[i, 2] = "n DSB under threshold limit"
                matMedPval[i, 3] = "n DSB under threshold limit"
                matMedPval[i, 4] = "n DSB under threshold limit"
              }
              
            } else{
              matMedPval[i, 1] = "harmonised file has a 0 size"
              matMedPval[i, 2] = "harmonised file has a 0 size"
              matMedPval[i, 3] = "harmonised file has a 0 size"
              matMedPval[i, 4] = "harmonised file has a 0 size"
              
              
            }
          } else{
            matMedPval[i, 1] = "check not ok"
            matMedPval[i, 2] = "check not ok"
            matMedPval[i, 3] = "check not ok"
            matMedPval[i, 4] = "check not ok"
            
            
            
          }
        } else{
          matMedPval[i, 1] = "NA in allele or location"
          matMedPval[i, 2] = "NA in allele or location"
          matMedPval[i, 3] = "NA in allele or location"
          
        }
        return(na.omit(matMedPval))
      })
      
      print(nrow(matMedPval))
      
      
      print(head(matMedPval))
      
      df <- data.frame()
      for (j in c(1:length(diseases))) {
        df[j, 1] <-  matMedPval[[j]][1, 1]
        
        df[j, 2] <-  matMedPval[[j]][1, 2]
        df[j, 3] <-  matMedPval[[j]][1, 3]
        df[j, 4] <-  matMedPval[[j]][1, 4]
        
      }
      
      
      matMedPval <- df
      colnames(matMedPval) = c(
        "Fisher_aggregartion_p_interestSNP",
        "FisherP_aggregartion_p_ctrlSNP",
        "p",
        "OR"
      )
      ###############################################################
      
      dir.create("Analysis_pval_fisher_greater")
      dir.create(paste0("Analysis_pval_fisher_greater/", ctrl_type))
      dir.create(paste0("Analysis_pval_fisher_greater/", ctrl_type, "/thresP"))
      write.csv(
        matMedPval,
        paste0(
          "Analysis_pval_fisher_greater/",
          ctrl_type,
          "/thresP",
          "/matMedPval_sup1000_harmonised_6000SNP_TAG_",
          ctrl_type,
          "_thresP_",
          thresP,
          ".csv"
        ),
        row.names = F
      )
      
      
      
      #correct fisher p_values (using fdr)
matCorrected<- matMedPval[matMedPval$p!= "NA in allele or location",]
nrow(matCorrected)

matCorrected <- matCorrected[matCorrected$p!= "n DSB under threshold limit",]
matCorrected <- matCorrected[matCorrected$p!= "n CTRL under threshold limit",]
matCorrected <- matCorrected[matCorrected$p!= "check not ok",]



matCorrected$p<- as.numeric(as.character(matCorrected$p))
matCorrected$fdr_with_all_1314_trait<- p.adjust(as.numeric(matCorrected$p), method = "fdr",n=nrow(matMedPval))
matCorrected$fdr_filterwith_nrows_trait<-p.adjust(as.numeric(matCorrected$p), method = "fdr",n=nrow(matCorrected))
write.table(matCorrected,paste0("Analysis_pval_fisher_greater/",ctrl_type,"/thresP","/matMedPvalcorrected_sup1000_harmonised_filtered_with_fdr_6000SNP_",ctrl_type,"_thresP_",thresP,"_TAG.tsv"),sep="\t",row.names = F)

    }
  }#)
}

