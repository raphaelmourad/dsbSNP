#09/10/2025
#Sébastien AUBER

# Tis script computes genomic inflation (lambda GC) 
# for DSB and control SNPs across GWAS traits

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
library(readr)
args = commandArgs(trailingOnly=TRUE)

#set working directory
setwd("/mnt/DATA/MCD/eq_legube/SEB/GWAS/")

# Load SNP types from input file
SNP_type_list = strsplit(read_lines(args[1]), " ")

# Load number of p-values under threshold
n_pval = args[2]

# Load study list and clean disease names
listStudy <- read_tsv(paste0("data/list_study_more_than_", n_pval, "_pvalue_under10Eminus6.tsv"))
diseases = listStudy$Disease
diseases <- gsub('/', '_', diseases) # replace '/' with '_'
listStudy$diseases = diseases
# diseases <- diseases[1:10]  # optional subset for testing

# Iterate over SNP types
for (SNP_type in SNP_type_list) {
  print(SNP_type)
  
  # Define list of control types
  ctrl_list = c("RDM","CTRL","PROM","ENH","SE")
  
  # Set working directory for results of current SNP type
  setwd(paste0("/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_noTAG_multi_", n_pval, "pval_under_10minus6/", SNP_type,"/"))
  
  # Base path for harmonized results
  path = paste0("/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_noTAG_multi_", n_pval, "pval_under_10minus6/", SNP_type,"/")
  
  # Iterate over diseases
  for (i in (1:length(diseases))) {
    # Initialize empty data frames for DSB and CTRL SNPs
    data_list_DSB <- data.frame()
    data_list_CTRL <- data.frame()
    
    print(diseases[i])
    print(i)
    
    # Load harmonized p-value data if file exists
    if(file.exists(paste0(path, "harmonised_tab_lists_All_rsID2/listPval_sup1000_harmonised_All_rsID_", diseases[i], "_SNP"))) {
      load(file=paste0("harmonised_tab_lists_All_rsID2/listPval_sup1000_harmonised_All_rsID_", diseases[i], "_SNP"))
      print(length(data_list_DSB$p_SNP))
      
      # Load control SNP p-values for all control types
      for (ctrl_type in ctrl_list){        
        load(file=paste0("harmonised_tab_lists_All_rsID2/listPval_sup1000_harmonised_All_rsID_", diseases[i], "_", ctrl_type))
      }
      
      # Extract numeric p-values and remove NA
      pvals <- as.numeric(data_list_DSB$p_SNP) %>% na.omit()
      pvals_CTRL <- as.numeric(data_list_CTRL$p_CTRL) %>% na.omit()
      pvals_RDM <- as.numeric(data_list_RDM$p_RDM) %>% na.omit()
      pvals_PROM <- as.numeric(data_list_prom$p_PROM) %>% na.omit()
      pvals_ENH <- as.numeric(data_list_ENH$p_ENH) %>% na.omit()
      pvals_SE <- as.numeric(data_list_SE$p_SE) %>% na.omit()
      
      # Convert p-values to chi-squared statistics
      chi2_stats <- qchisq(1 - pvals, df = 1)
      chi2_stats_RDM <- qchisq(1 - pvals_RDM, df = 1)
      chi2_stats_CTRL <- qchisq(1 - pvals_CTRL, df = 1)
      chi2_stats_PROM <- qchisq(1 - pvals_PROM, df = 1)
      chi2_stats_ENH <- qchisq(1 - pvals_ENH, df = 1)
      chi2_stats_SE <- qchisq(1 - pvals_SE, df = 1)
      
      # Compute median chi-squared values
      median_chi2 <- median(chi2_stats, na.rm = TRUE)
      median_chi2_RDM <- median(chi2_stats_RDM, na.rm = TRUE)
      median_chi2_CTRL <- median(chi2_stats_CTRL, na.rm = TRUE)
      median_chi2_PROM <- median(chi2_stats_PROM, na.rm = TRUE)
      median_chi2_ENH <- median(chi2_stats_ENH, na.rm = TRUE)
      median_chi2_SE <- median(chi2_stats_SE, na.rm = TRUE)
      
      # Compute lambda GC relative to each control
      lambda_gc_over_RDM <- median_chi2 / median_chi2_RDM
      lambda_gc_over_CTRL <- median_chi2 / median_chi2_CTRL
      lambda_gc_over_PROM <- median_chi2 / median_chi2_PROM
      lambda_gc_over_ENH <- median_chi2 / median_chi2_ENH
      lambda_gc_over_SE <- median_chi2 / median_chi2_SE
      
      print(paste0(diseases[i], ":lambda_gc=", lambda_gc_over_RDM))
      
      # Compute absolute lambda GC (using expected median of chi-squared with df=1)
      lambda_gc_SNP <- median_chi2 / qchisq(0.5, df = 1)
      lambda_gc_RDM <- median_chi2_RDM / qchisq(0.5, df = 1)
      lambda_gc_CTRL <- median_chi2_CTRL / qchisq(0.5, df = 1)
      lambda_gc_PROM <- median_chi2_PROM / qchisq(0.5, df = 1)
      lambda_gc_ENH <- median_chi2_ENH / qchisq(0.5, df = 1)
      lambda_gc_SE <- median_chi2_SE / qchisq(0.5, df = 1)
      
      # Get EFO code for current disease
      EFO = listStudy$EFO[listStudy$diseases == diseases[i]]
      
      # Create data frame with lambda GC and median chi2 values
      df = data.frame(lambda_gc_over_RDM, lambda_gc_over_CTRL, lambda_gc_over_PROM,
                      lambda_gc_over_ENH, lambda_gc_over_SE,
                      lambda_gc_SNP, lambda_gc_RDM, lambda_gc_CTRL, lambda_gc_PROM, 
                      lambda_gc_ENH, lambda_gc_SE,
                      median_chi2, median_chi2_RDM, median_chi2_CTRL, median_chi2_PROM, 
                      median_chi2_ENH, median_chi2_SE,
                      EFO)
      df$diseases = diseases[i]
      
      # Combine results across diseases
      if(i == 1){
        dff = df
      } else {
        dff = rbind(dff, df)
      } 
      
    }
  }
}

# Final output
dff

# Create output directory and write results
dir.create(paste0("/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_noTAG_multi_", n_pval, "pval_under_10minus6/", SNP_type,"/lambda"))
write_tsv(dff, paste0("/mnt/DATA/MCD/eq_legube/SEB/GWAS/results_noTAG_multi_", n_pval, "pval_under_10minus6/", SNP_type,"/lambda/lamba_score_RDM_All_col.tsv"))

