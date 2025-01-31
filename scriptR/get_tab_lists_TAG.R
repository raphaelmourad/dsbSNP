######################################
#THE MAIN GOAL OF THIS SCRIPT IS TO OBTAIN THE OVERLAPPING dsbSNPs or ctrlSNPs
#AND THEIR P-VAL IN GWAS FOR EVERY INTEREST GWAS

######################################

require(gwascat)
require(tidyr)
require(dplyr)
library(tidyverse)

library(Matrix)
library(MASS)
require(parallel)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(liftOver)
library(rtracklayer)
library(plyranges)
library(rlist)
library(poolr)
args = commandArgs(trailingOnly = TRUE)


#this function test that a correct download of the gwas
md5check <- function(path = path, id_) {
  gzfile <- list.files(paste0(path, id_), pattern = "h.tsv.gz")
  
  md5 <- list.files(paste0(path, id_), pattern = "md5sum.txt")
  if (identical(md5, character(0))) {
    return("not_OK")
    write(
      paste0(id_, " md5sum.txt not available"),
      file = paste0(path, "fail_download.txt"),
      append = TRUE
    )
  } else{
    if (identical(gzfile, character(0))) {
      return("not_OK")
      write(
        paste0(id_, " not downloaded"),
        file = paste0(path, "fail_download.txt"),
        append = TRUE
      )
    } else{
      system(command = paste0(
        "md5sum ",
        path,
        id_,
        gzfile,
        " >",
        path,
        id_,
        "nchecksum.txt"
      ))
      ncs <- read.csv(paste0(path, id_, "nchecksum.txt"),
                      sep = " ",
                      header = F)
      cs <- read.csv(paste0 (path, id_, md5),
                     sep = " ",
                     header = F)
      ncs$V1[ncs$V3 == paste0(path, id_, gzfile)]
      cs$V1[cs$V2 == paste0(gzfile)]
      if (ncs$V1[ncs$V3 == paste0(path, id_, gzfile)] != cs$V1[cs$V2 == paste0(gzfile)]) {
        print (paste0 ("download incorrect", gzfile))
        return("not_OK")
        write(
          paste0(id_, " no checksum"),
          file = paste0(path, "fail_download.txt"),
          append = TRUE
        )
        
      } else{
        print(paste0(gzfile, " md5check ok"))
        write(id_,
              file = paste0(path, "succeded_download.txt"),
              append = TRUE)
        
        return("OK")
      }
    }
  }
}


#fix the r2 threshold of TAG SNPs
thres = 0.99



#get the Type of SNPs to request (eg. dsb_no_eSNP)
SNP_type_list = strsplit(read_lines(args[1]), " ")
n_pval = args[2]#minimum number of pvals under 10e-6 for a GWAS to be included


#load interest SNPs and convert to hg38
for (SNP_type in SNP_type_list) {
  load(
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_GR.RDATA"
    )
  )
  
  
  SNP.GR = SNP_TAG.GR
  
  
  load(
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_CTRL_dist200_500kb.RDATA"
    )
  )
  load(
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_PROM_dist200_500kb.RDATA"
    )
  )
  load(
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_ENH_dist200_500kb.RDATA"
    )
  )
  load(
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_SE_dist200_500kb.RDATA"
    )
  )
  
  SNP_ctrl.GR = SNP_TAG_CTRL.GR
  SNP_prom.GR = SNP_TAG_PROM.GR
  SNP_ENH.GR = SNP_TAG_ENH.GR
  SNP_SE.GR = SNP_TAG_SE.GR
  
  length(SNP.GR)
  length(SNP_ctrl.GR)
  length(SNP_prom.GR)
  length(SNP_ENH.GR)
  length(SNP_SE.GR)
  
  
  lpath = "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/hg19ToHg38.over.chain"
  ch19to38 = import.chain(lpath) #on va tout passer en hg 38
  
  
  #  length(SNP_TAG_sup_0.99_sample56k.GR)
  #length(SNP_ctrl.GR )
  #length(SNP_prom.GR)
  #length(SNP_ENH.GR )
  
  #length(SNP_SE.GR )
  
  SNP.GR = liftOver(SNP.GR, ch19to38)
  SNP.GR = unlist(SNP.GR)
  
  
  SNP_ctrl.GR = liftOver(SNP_ctrl.GR, ch19to38)
  SNP_ctrl.GR = unlist(SNP_ctrl.GR)
  
  SNP_prom.GR = liftOver(SNP_prom.GR, ch19to38)
  SNP_prom.GR = unlist(SNP_prom.GR)
  
  SNP_ENH.GR = liftOver(SNP_ENH.GR, ch19to38)
  SNP_ENH.GR = unlist(SNP_ENH.GR)
  
  
  SNP_SE.GR = liftOver(SNP_SE.GR, ch19to38)
  SNP_SE.GR = unlist(SNP_SE.GR)
  
  length(SNP.GR)
  length(SNP_ctrl.GR)
  length(SNP_prom.GR)
  length(SNP_ENH.GR)
  length(SNP_SE.GR)
  
  save(
    SNP.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_dist200_500kb_hg38.RDATA"
    )
  )
  save(
    SNP_ctrl.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_CTRL_dist200_500kb_hg38.RDATA"
    )
  )
  save(
    SNP_prom.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_PROM_dist200_500kb_hg38.RDATA"
    )
  )
  save(
    SNP_ENH.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_ENH_dist200_500kb_hg38.RDATA"
    )
  )
  save(
    SNP_SE.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_SE_dist200_500kb_hg38.RDATA"
    )
  )
  
  
  write_bed(
    SNP.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_dist200_500kb_hg38.bed"
    )
  )
  write_bed(
    SNP_ctrl.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_CTRL_dist200_500kb_hg38.bed"
    )
  )
  write_bed(
    SNP_prom.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_PROM_dist200_500kb_hg38.bed"
    )
  )
  write_bed(
    SNP_ENH.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_ENH_dist200_500kb_hg38.bed"
    )
  )
  write_bed(
    SNP_SE.GR,
    file = paste0(
      "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
      SNP_type,
      "/",
      "TAG/",
      thres,
      "/SNP_TAG_sup_0.99_SE_dist200_500kb_hg38.bed"
    )
  )
  
  
}






path = "/mnt/DATA/MCD/eq_legube/SEB/GWAS/"

setwd(path)

####
#get interest GWAS trait and diseases

listStudy <- read_tsv(paste0(
  "./data/list_study_more_than_",
  n_pval,
  "_pvalue_under10Eminus6.tsv"
))
diseases = listStudy$Disease
diseases <- gsub('/', '_', diseases)
print(diseases)
df <- tibble()
df.f <- tibble()
listPval = list()
listES = list()
listBeta = list()
numPvalLow = rep(NA, length(diseases))
print(SNP_type_list)


for (SNP_type in SNP_type_list) {
  dir.create(paste0("results_TAG_multi_", n_pval, "pval_under_10minus6/"))
  dir.create(paste0(
    "results_TAG_multi_",
    n_pval,
    "pval_under_10minus6/",
    SNP_type
  ))
  dir.create(
    paste0(
      "results_TAG_multi_",
      n_pval,
      "pval_under_10minus6/",
      SNP_type,
      "/harmonised_tab_lists_TAG_rsID2"
    )
  )
  dir.create(
    paste0(
      "results_TAG_multi_",
      n_pval,
      "pval_under_10minus6/",
      SNP_type,
      "/qqPlots_TAG_rsID2"
    )
  )
}


#get the tables of overlapping pvalues and qqplots

for (i in (1:length(diseases))) {
  listPval = list()
  data_list_DSB <- data.frame()
  data_list_CTRL <- data.frame()
  rsID <- data.frame()
  disease = diseases[i]
  print(i)
  print(paste0(i / length(diseases) * 100))
  id = listStudy$study_id[i]
  
  #load GWAS summary statistics
  
  #  if (try(md5check(path="/mnt/DATA/MCD/eq_legube/SEB/GWAS/SumStat/",id_=paste0(path,"SumStat/",id,"/")))=="OK"){
  gzfile <- list.files(paste0(path, "SumStat/", id, "/"), pattern = "h.tsv.gz")
  dataSumStat = fread(
    file = paste0(path, "SumStat/", id, "/", gzfile),
    header = T,
    sep = "\t"
  )
  
  
  for (SNP_type in SNP_type_list) {
    #load interest SNPs
    load(
      paste0(
        "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
        SNP_type,
        "/",
        "TAG/",
        thres,
        "/SNP_TAG_sup_0.99_dist200_500kb_hg38.RDATA"
      )
    )
    load(
      paste0(
        "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
        SNP_type,
        "/",
        "TAG/",
        thres,
        "/SNP_TAG_sup_0.99_CTRL_dist200_500kb_hg38.RDATA"
      )
    )
    load(
      paste0(
        "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
        SNP_type,
        "/",
        "TAG/",
        thres,
        "/SNP_TAG_sup_0.99_PROM_dist200_500kb_hg38.RDATA"
      )
    )
    load(
      paste0(
        "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
        SNP_type,
        "/",
        "TAG/",
        thres,
        "/SNP_TAG_sup_0.99_ENH_dist200_500kb_hg38.RDATA"
      )
    )
    load(
      paste0(
        "/mnt/DATA/MCD/eq_legube/SEB/GWAS/data/",
        SNP_type,
        "/",
        "TAG/",
        thres,
        "/SNP_TAG_sup_0.99_SE_dist200_500kb_hg38.RDATA"
      )
    )
    
    
    print(length(SNP.GR))
    print(length(SNP_ctrl.GR))
    print(length(SNP_prom.GR))
    print(length(SNP_ENH.GR))
    print(length(SNP_SE.GR))
    
    
    
    #get corresonding rsID
    if (length(dataSumStat) != 0) {
      # dataSumStat=fread(file=fileSumStats,header=T,sep="\t")
      
      if (sum(colnames(dataSumStat) %in% "beta") == 0) {
        dataSumStat$beta = log(dataSumStat$odds_ratio)
      }
      
      if (length(which(is.na(dataSumStat$base_pair_location))) != nrow(dataSumStat) &
          length(which(is.na(dataSumStat$base_pair_location))) == 0 &
          length(which(is.na(dataSumStat$other_allele))) == 0) {
        if (!is.na(dataSumStat$hm_rsid[1]) &
            length(which(!is.na(dataSumStat$hm_rsid))) == nrow(dataSumStat)) {
          print("hm rsid")
          sumStat.GR = GRanges(
            paste0("chr", dataSumStat$chromosome),
            IRanges(
              dataSumStat$base_pair_location,
              dataSumStat$base_pair_location
            ),
            ref = dataSumStat$other_allele,
            alt = dataSumStat$effect_allele,
            p = dataSumStat$p_value,
            beta = dataSumStat$beta,
            rsID = dataSumStat$hm_rsid
          )
          
          
        } else{
          if (!is.na(dataSumStat$variant_id[1]) &
              length(which(!is.na(dataSumStat$variant_id))) == nrow(dataSumStat)) {
            print("variant rsid")
            sumStat.GR = GRanges(
              paste0("chr", dataSumStat$chromosome),
              IRanges(
                dataSumStat$base_pair_location,
                dataSumStat$base_pair_location
              ),
              ref = dataSumStat$other_allele,
              alt = dataSumStat$effect_allele,
              p = dataSumStat$p_value,
              beta = dataSumStat$beta,
              rsID = dataSumStat$variant_id
            )
            
          } else{
            print("no rsid")
            sumStat.GR = GRanges(
              paste0("chr", dataSumStat$chromosome),
              IRanges(
                dataSumStat$base_pair_location,
                dataSumStat$base_pair_location
              ),
              ref = dataSumStat$other_allele,
              alt = dataSumStat$effect_allele,
              p = dataSumStat$p_value,
              beta = dataSumStat$beta,
              rsID = "noID"
            )
            
          }
        }
        
        
        
        numPvalLow[i] = sum(sumStat.GR$p < 1e-6)
        
        
        #get Overlapps
        olSumStat = findOverlaps(SNP.GR, sumStat.GR)
        SNp_SNP_sumStats.GR = SNP.GR[queryHits(olSumStat)]
        sumStat_DSB.GR = sumStat.GR[subjectHits(olSumStat)]
        sumStat_DSB.GR$fdr = p.adjust(sumStat_DSB.GR$p, method = "fdr")
        sumStat_DSB.GR$pbonf = p.adjust(sumStat_DSB.GR$p, method = "bonferroni")
        
        
        olSumStat_ctrl = findOverlaps(SNP_ctrl.GR, sumStat.GR)
        SNP_ctrl_sumStats.GR = SNP_ctrl.GR[queryHits(olSumStat_ctrl)]
        sumStat_ctrl.GR = sumStat.GR[subjectHits(olSumStat_ctrl)]
        sumStat_ctrl.GR$fdr = p.adjust(sumStat_ctrl.GR$p, method = "fdr")
        
        olSumStat_prom = findOverlaps(SNP_prom.GR, sumStat.GR)
        SNP_prom_sumStats.GR = SNP_prom.GR[queryHits(olSumStat_prom)]
        sumStat_prom.GR = sumStat.GR[subjectHits(olSumStat_prom)]
        sumStat_prom.GR$fdr = p.adjust(sumStat_prom.GR$p, method = "fdr")
        
        
        olSumStat_ENH = findOverlaps(SNP_ENH.GR, sumStat.GR)
        SNP_ENH_sumStats.GR = SNP_ENH.GR[queryHits(olSumStat_ENH)]
        sumStat_ENH.GR = sumStat.GR[subjectHits(olSumStat_ENH)]
        sumStat_ENH.GR$fdr = p.adjust(sumStat_ENH.GR$p, method = "fdr")
        
        olSumStat_SE = findOverlaps(SNP_SE.GR, sumStat.GR)
        SNP_SE_sumStats.GR = SNP_SE.GR[queryHits(olSumStat_SE)]
        sumStat_SE.GR = sumStat.GR[subjectHits(olSumStat_SE)]
        sumStat_SE.GR$fdr = p.adjust(sumStat_SE.GR$p, method = "fdr")
        
        #export qqplots
        file_enrich_GWAS = paste0(
          path,
          "results_TAG_multi_",
          n_pval,
          "pval_under_10minus6/",
          SNP_type,
          "/qqPlots_TAG_rsID2/CTRL_SNP_",
          disease,
          "_GWAS.pdf"
        )
        try(pdf(file_enrich_GWAS, 4, 4))
        try(qqplot(
          -log10(as.numeric(sumStat_ctrl.GR$p)),
          -log10(as.numeric(sumStat_DSB.GR$p)),
          xlim = c(0, 15),
          ylim = c(0, 15),
          xlab = "Ctrl SNP -log10(p)",
          ylab = "dsbSNP -log10(p)",
          pch = 18
        ))
        try(abline(0, 1, col = "red"))
        try(dev.off())
        
        
        file_enrich_GWAS = paste0(
          path,
          "results_TAG_multi_",
          n_pval,
          "pval_under_10minus6/",
          SNP_type,
          "/qqPlots_TAG_rsID2/Promoter_SNP_",
          disease,
          "_GWAS.pdf"
        )
        try(pdf(file_enrich_GWAS, 4, 4))
        try(qqplot(
          -log10(as.numeric(sumStat_prom.GR$p)),
          -log10(as.numeric(sumStat_DSB.GR$p)),
          xlim = c(0, 15),
          ylim = c(0, 15),
          xlab = "Ctrl promoter SNP -log10(p)",
          ylab = "dsbSNP -log10(p)",
          pch = 18
        ))
        try(abline(0, 1, col = "red"))
        try(dev.off())
        
        
        
        
        
        
        file_enrich_GWAS = paste0(
          path,
          "results_TAG_multi_",
          n_pval,
          "pval_under_10minus6/",
          SNP_type,
          "/qqPlots_TAG_rsID2/ENH_SNP_",
          disease,
          "_GWAS.pdf"
        )
        try(pdf(file_enrich_GWAS, 4, 4))
        try(qqplot(
          -log10(as.numeric(sumStat_ENH.GR$p)),
          -log10(as.numeric(sumStat_DSB.GR$p)),
          xlim = c(0, 15),
          ylim = c(0, 15),
          xlab = "Ctrl_enhancer SNP -log10(p)",
          ylab = "dsbSNP -log10(p)",
          pch = 18
        ))
        try(abline(0, 1, col = "red"))
        try(dev.off())
        
        
        
        file_enrich_GWAS = paste0(
          path,
          "results_TAG_multi_",
          n_pval,
          "pval_under_10minus6/",
          SNP_type,
          "/qqPlots_TAG_rsID2/SE_SNP_",
          disease,
          "_GWAS.pdf"
        )
        try(pdf(file_enrich_GWAS, 4, 4))
        try(qqplot(
          -log10(as.numeric(sumStat_SE.GR$p)),
          -log10(as.numeric(sumStat_DSB.GR$p)),
          xlim = c(0, 15),
          ylim = c(0, 15),
          xlab = "Ctrl Super_enhancer SNP -log10(p)",
          ylab = "dsbSNP -log10(p)",
          pch = 18
        ))
        try(abline(0, 1, col = "red"))
        try(dev.off())
        
        
        listPval[[i]] = list(
          p_SNP = sumStat_DSB.GR$p,
          p_ctrl = sumStat_ctrl.GR$p,
          p_prom = sumStat_prom.GR$p,
          p_ENH = sumStat_ENH.GR$p,
          p_SE = sumStat_SE.GR$p
        )
        rsID <- as.data.frame(sumStat_DSB.GR$rsID)
        
        # print(round(i/length(diseases),1))
        
        
        
        
        
        
      } else{
        listPval[[i]] = "NA in allele or location"
        
        listPval[[i]]$p_SNP = "NA in allele or location"
        listPval[[i]]$p_ctrl = "NA in allele or location"
        listPval[[i]]$p_prom = "NA in allele or location"
        listPval[[i]]$p_ENH = "NA in allele or location"
        listPval[[i]]$p_SE = "NA in allele or location"
        listES[[i]] = "NA in allele or location"
        listBeta[[i]] = "NA in allele or location"
        rsID = "NA in allele or location"
        sumStat_DSB.GR = data_frame()
        sumStat_DSB.GR$rsID <- NULL
      }
      
      
    } else{
      print("harmonised file has a 0 size")
      listPval[[i]] = "harmonised file has a 0 size"
      
      listPval[[i]]$p_SNP = "harmonised file has a 0 size"
      listPval[[i]]$p_ctrl = "harmonised file has a 0 size"
      listPval[[i]]$p_prom = "harmonised file has a 0 size"
      listPval[[i]]$p_SE = "harmonised file has a 0 size"
      listPval[[i]]$p_ENH = "harmonised file has a 0 size"
      
      listES[[i]] = "harmonised file has a 0 size"
      listBeta[[i]] = "harmonised file has a 0 size"
      
      rsID = "harmonised file has a 0 size"
      
      sumStat_DSB.GR = data_frame()
      sumStat_DSB.GR$rsID <- NULL
      #   try(system(command= paste0("rm -r ", path,id)))
      
    }
    
    data_list_DSB <- as.data.frame(listPval[[i]]$p_SNP)
    
    if (is.null(sumStat_DSB.GR$rsID)) {
      data_list_DSB$listID <- "no_rsID"
      print("no rsID")
    } else{
      data_list_DSB$listID <- rsID
    }
    colnames(data_list_DSB) <- c("p_SNP", "rsID")#,"listES","listBeta")
    
    data_list_CTRL <- as.data.frame(listPval[[i]]$p_ctrl)
    colnames(data_list_CTRL) <- c("p_CTRL")
    
    data_list_prom <- as.data.frame(listPval[[i]]$p_prom)
    colnames(data_list_prom) <- c("p_PROM")
    
    data_list_ENH <- as.data.frame(listPval[[i]]$p_ENH)
    colnames(data_list_ENH) <- c("p_ENH")
    
    data_list_SE <- as.data.frame(listPval[[i]]$p_SE)
    colnames(data_list_SE) <- c("p_SE")
    print(nrow(data_list_DSB))
    print(nrow(data_list_CTRL))
    
    #save list of pval (tab_list)
    save(
      data_list_DSB,
      file = paste0(
        "results_TAG_multi_",
        n_pval,
        "pval_under_10minus6/",
        SNP_type,
        "/harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
        disease,
        "_SNP"
      )
    )
    save(
      data_list_CTRL,
      file = paste0(
        "results_TAG_multi_",
        n_pval,
        "pval_under_10minus6/",
        SNP_type,
        "/harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
        disease,
        "_CTRL"
      )
    )
    save(
      data_list_prom,
      file = paste0(
        "results_TAG_multi_",
        n_pval,
        "pval_under_10minus6/",
        SNP_type,
        "/harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
        disease,
        "_PROM"
      )
    )
    
    save(
      data_list_ENH,
      file = paste0(
        "results_TAG_multi_",
        n_pval,
        "pval_under_10minus6/",
        SNP_type,
        "/harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
        disease,
        "_ENH"
      )
    )
    save(
      data_list_SE,
      file = paste0(
        "results_TAG_multi_",
        n_pval,
        "pval_under_10minus6/",
        SNP_type,
        "/harmonised_tab_lists_TAG_rsID2/listPval_sup1000_harmonised_TAG_rsID_",
        disease,
        "_SE"
      )
    )
    
    i = i + 1
  }
}

