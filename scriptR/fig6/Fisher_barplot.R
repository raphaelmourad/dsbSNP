#################################################################################
# 09/10/2025
# Sébastien AUBER
# Assemble a comprehensive table of all diseases with SNP counts, percentages, 
# OR, FDR, and lambda statistics, adding disease labels from EFO codes, 
# generate bar plots of percentages and counts ,
# and export final TSV/XLSX files for downstream analysis.
#################################################################################



library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggpubr)
library(stringr)


DSB_path <- "/media/sauber/Elements/DSBsnp_project_seb/results/allelic_imbalance/slurm_outputs/"
setwd(DSB_path)

SNP_type_list <- c("SNP_DSB")
ctrl_list <- c("CTRL","RDM")
thresP <- 1e-08
TAG <- "noTAG"

# read files
dff <- NULL
if(TAG == "noTAG"){
  for(ctrl_type in ctrl_list){
    for(SNP_type in SNP_type_list){
      df <- read_tsv(paste0(
        "results_noTAG_multi_500pval_under_10minus6/", SNP_type, 
        "/Analysis_pval_fisher_greater_fused_dup_pseudo_count_at_1_with_tot_RDM/",
        ctrl_type, "/thresP/matMedPvalcorrected_sup1000_harmonised_filtered_with_fdr_6000SNP_",
        ctrl_type,"_thresP_1e-08_noTAG.tsv"
      ))
      df$SNP_type <- SNP_type
      df$CTRL <- ctrl_type
      dff <- if(is.null(dff)) df else rbind(dff, df)
    }
  }
}

# Read eSNP
eSNP <- read_tsv(paste0(
  "results_noTAG_multi_500pval_under_10minus6/7000_eSNP_no_dsb/",
  "Analysis_pval_fisher_greater_fused_dup_pseudo_count_at_1_with_tot_RDM/", ctrl_type,
  "/thresP/matMedPvalcorrected_sup1000_harmonised_filtered_with_fdr_6000SNP_", ctrl_type, "_thresP_1e-08_noTAG.tsv"
))

#  annotation
setwd("/media/sauber/Elements/DSBsnp_project_seb/")
Annotation <- read_tsv("All_harmonised_GWAS_pvalue_under1e-6_above500_gwas_with_traits_annotated_by_class_all_diseaese_currated.tsv")
Annotation$EFO <- Annotation$EFO %>% str_split(":") %>% sapply(paste, collapse="_") %>% paste0(".h.tsv.gz")
annotated_dff <- merge(Annotation, dff, by="EFO") %>% na.omit()


###keep only GWAS in disease classes

DiseaseList <- c("cancer.or.benign.tumor",
                 "cardiovascular.disease",
                 "chronic.disease",
                 "connective.tissue.disease",
                 "digestive.system.disease",
                 "disorder.of.visual.system",
                 "endocrine.system.disease",
                 "genetic.disorder",
                 "head.and.neck.disorder",
                 "immune.system.disease",
                 "infectious.disease",
                 "inflammatory.disease",
                 "integumentary.system.disease",
                 "intestinal.disease",
                 "metabolic.disease",
                 "musculoskeletal.system.disease",
                 "nervous.system.disease",
                 "otorhinolaryngologic.disease",
                 "psychiatric.disorder",
                 "reproductive.system.disease",
                 "respiratory.system.disease",
                 "urinary.system.disease")

ouais<-data.frame()
final_annontation_dff<-data.frame()
for (clas in DiseaseList){
  ouais<-(annotated_dff[annotated_dff$name.x==clas,])
  ouais$class<- rep(clas,nrow(ouais))
  final_annontation_dff<-(rbind(final_annontation_dff,ouais))
  
}
dff<-dff[dff$EFO %in% final_annontation_dff$EFO,]

#  FDR
annotated_dff <- annotated_dff %>%
  group_by(ctrl_type) %>%
  mutate(fdr = p.adjust(p, method = "fdr")) %>%
  ungroup()


#  eSNP
eS <- eSNP %>%
  select(name, n_mean, mean_tot, EFO) %>%
  rename(name.y = name) %>%
  mutate(perc_eSNP = (n_mean / mean_tot) * 100) %>%
  select(-EFO, -n_mean, -mean_tot)

# percentage computation
annotated_dff$perc_n <- (annotated_dff$n_mean / annotated_dff$mean_tot) * 100
annotated_dff$perc_c <- (annotated_dff$n_ctrl_mean / annotated_dff$mean_tot_control) * 100

# Sélection et fusion pour graphique
plot_df <- annotated_dff %>%
  filter(ctrl_type == "RDM") %>%
  select(OR, fdr, name.x, name.y, perc_n, perc_c) %>%
  filter(fdr < 0.05) %>%
  rename(perc_n_no_dist = perc_n, perc_c_no_dist = perc_c) %>%
  unique() %>%
  left_join(
    annotated_dff %>% filter(ctrl_type == "CTRL") %>% select(name.x, name.y, perc_n, perc_c)
  )


plot_df20<- plot_df
#merge with eSNP
unique(plot_df20$name.y)
eS<-eSNP%>%dplyr::select(name,n_mean,mean_tot,EFO)%>%
  dplyr:: rename(name.y=name)%>%
  #  mutate(ctrl_type="eSNP")%>%
  mutate(perc_eSNP=(n_mean/mean_tot)*100)%>%dplyr::select(-EFO ,n_mean ,mean_tot)

plot_df20<-left_join(plot_df20,eS)


plot_df20<-plot_df20%>%#mutate(perc_c=perc_c/10+1)%>%mutate(perc_n=perc_n+1)%>%
  dplyr:: rename(dsbSNP=perc_n)%>%dplyr:: rename(CTRL=perc_c)%>%dplyr::rename(CTRL_no_dist=perc_c_no_dist)%>%dplyr::rename(eSNP=perc_eSNP)%>%
  pivot_longer(cols =c(CTRL,dsbSNP,CTRL_no_dist,eSNP), values_to= "mean"
  )
p2<-plot_df20%>% ggplot(aes(x=forcats:: fct_reorder(name.y,-OR)  ,y= (mean),fill=(name)))+geom_bar(position="dodge", stat="identity")+ 
  scale_fill_manual(values = c( "blue","lightblue","red","violet"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=7,colour = "black"))+
  theme(axis.text.y = element_text(size=7,colour = "black"))+
  ylab("Percentage of significant SNPs")+xlab('')+ggtitle(paste0(SNP_type," no pseudocompte"))
p2

#order BARs
plot_df20$name <- factor(plot_df20$name, 
                         levels = c( "dsbSNP", "CTRL_no_dist","CTRL", "eSNP"))
#plot
p2 <- plot_df20 %>%
  ggplot(aes(x = forcats::fct_reorder(name.y, -OR), 
             y = mean, 
             fill = name)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("CTRL" = "blue",
                               "CTRL_no_dist" = "lightblue",
                               "dsbSNP" = "red",
                               "eSNP" = "violet")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7, colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black")) +
  ylab("Percentage of significant SNPs") +
  xlab('') +
  ggtitle(paste0(SNP_type, " no pseudocompte"))
p2

#save
pdf(paste0("results/allelic_imbalance/top",top,"_fisher_fused_tables_fdr_1e-08_order_by_OR_percentage_RDM_eSNP.pdf"),9,5)
p2
dev.off()

################order BAR

plot_df20$name <- factor(plot_df20$name, 
                         levels = c("dsbSNP", "CTRL_no_dist","CTRL", "eSNP"))

p2 <- plot_df20 %>%
  ggplot(aes(x = forcats::fct_reorder(name.y, -OR), 
             y = mean, 
             fill = name)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("CTRL" = "blue",
                               "CTRL_no_dist" = "lightblue",
                               "dsbSNP" = "red",
                               "eSNP" = "violet")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7, colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black")) +
  ylab("Percentage of significant SNPs") +
  xlab('') +
  ggtitle(paste0(SNP_type, " no pseudocompte"))
p2

pdf(paste0("results/allelic_imbalance/top",top,"_fisher_fused_tables_fdr_1e-08_order_by_OR_percentage_RDM_eSNP.pdf"),9,5)
p2
dev.off()

################# Count table


plot_df_count<-annotated_dff[annotated_dff$name.y %in% plot_df20$name.y,]

plot_df_count<-plot_df_count%>%dplyr::select(name.y,n_mean,n_ctrl_mean,ctrl_type)%>%unique()%>%dplyr::rename(dsbSNP=n_mean)%>%
  left_join(.,eSNP%>%dplyr::select(name,n_mean)%>%unique()%>%dplyr::rename(eSNP=n_mean)%>%dplyr::rename(name.y=name))%>%
  pivot_wider(names_from = ctrl_type,values_from = n_ctrl_mean)%>%
  pivot_longer(c(dsbSNP, eSNP ,   RDM , CTRL),values_to = "Count")


ordre_x <- plot_df20 %>%
  dplyr::arrange(desc(OR)) %>%
  dplyr::pull(name.y) %>%
  unique()

plot_df_count$name.y <- factor(plot_df_count$name.y, levels = ordre_x)

pc <- plot_df_count %>%
  ggplot(aes(x = name.y, y = Count, fill = name)) +
  geom_bar(position = "dodge", stat = "identity") + 
  scale_fill_manual(values = c("blue","lightblue","violet","red")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7, colour = "black")) +
  theme(axis.text.y = element_text(size = 7, colour = "black")) +
  ylab("Count of significant SNPs") +
  xlab('') +
  ggtitle(paste0(SNP_type, " no pseudocompte"))

pc
p2

#### Merge count and percentage

plot_df_count <- plot_df_count %>%
  mutate(name = ifelse(name == "RDM", "CTRL_no_dist", name))

df_merge <- plot_df20 %>%
  left_join(plot_df_count, by = c("name.y", "name")) %>%
  mutate(Count = round(Count, 0))

ordre_x <- df_merge %>%
  arrange(desc(OR)) %>%
  pull(name.y) %>%
  unique()

df_merge$name.y <- factor(df_merge$name.y, levels = ordre_x)

p2_count <- ggplot(df_merge, aes(x = name.y, y = mean, fill = name)) +
  geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c("blue", "lightblue", "violet", "red")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7, colour = "black")) +
  theme(axis.text.y = element_text(size = 7, colour = "black")) +
  ylab("Percentage of significant SNPs") +
  xlab("") +
  ggtitle(paste0(SNP_type, " no pseudocompte"))

p2_count

#### FDR ordering

df_merge2 <- df_merge %>%
  arrange(fdr)

ordre_x <- df_merge2 %>%
  pull(name.y) %>%
  unique()

df_merge2$name.y <- factor(df_merge2$name.y, levels = ordre_x)

p2_fdr <- ggplot(df_merge2, aes(x = name.y, y = mean, fill = name)) +
  geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c("blue", "lightblue", "violet", "red")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7, colour = "black")) +
  theme(axis.text.y = element_text(size = 7, colour = "black")) +
  ylab("Percentage of significant SNPs") +
  xlab("") +
  ggtitle(paste0(SNP_type, " no pseudocompte"))

csv <- df_merge2 %>%
  dplyr::select(name.y, mean, name) %>%
  unique() %>%
  rename(percentage = mean) %>%
  pivot_wider(names_from = name, values_from = percentage)

write_tsv(csv, "top20_lambda.tsv")

library(ggpubr)
ggarrange(p2, pc, ncol = 1, nrow = 2, align = "v")

pdf(paste0("results/allelic_imbalance/top",top,"_fisher_fused_tables_count_1_1e-08_percentage_RDM_eSNP.pdf"),9,5)
p2_count
dev.off()

top20_list <- unique(df_merge$name.y)
save(top20_list, file = paste0("results/allelic_imbalance/top",top,"diseases_",thresP,".RDATA"))

print(p2_fdr)


top20_list_fdr <- unique(df_merge2$name.y)
save(top20_list_fdr, file = paste0("results/allelic_imbalance/top",top,"diseases_fdr_",thresP,".RDATA"))

top20_list_fdr_EFO <- unique(dff$EFO[dff$name %in% top20_list_fdr])
order_table <- dff[dff$name %in% top20_list_fdr,]

order_table <- left_join(
  order_table %>% dplyr::select(name, EFO) %>% dplyr::rename(name.y = name),
  df_merge2 %>% dplyr::select(name.y, fdr, OR)
)

DSB_path <- "/media/sauber/Elements/DSBsnp_project_seb/"
setwd(DSB_path)

save(top20_list_fdr_EFO, file = paste0("results/allelic_imbalance/top",top,"_EFO_diseases_fdr_",thresP,".RDATA"))
save(order_table, file = paste0("results/allelic_imbalance/table_fdr_top",top,"_EFO_diseases_fdr_",thresP,".RDATA"))
############################
# Build beautiful TSV

# Only significant


plot_df<- annotated_dff%>%filter(ctrl_type=="RDM")%>%dplyr::select(OR,fdr,name.x,name.y,perc_n,perc_c)%>%filter(fdr<0.05)%>%
  
  dplyr:: rename(perc_n_no_dist=perc_n)%>%dplyr::rename(perc_c_no_dist=perc_c)  %>%
  #filter(OR>2)%>%
  unique()%>%#pivot_wider(names_from = name.x,values_from = value)%>%
  left_join(. , annotated_dff%>%filter(ctrl_type=="CTRL")%>%dplyr::select(name.x,name.y,perc_n,perc_c))
#top=20
perc_tsv<-left_join(plot_df,eS)%>%#mutate(perc_c=perc_c/10+1)%>%mutate(perc_n=perc_n+1)%>%
  dplyr::rename(dsbSNP=perc_n)%>%  dplyr::rename(CTRL=perc_c)%>%dplyr::rename(RDM=perc_c_no_dist)%>%dplyr::rename(eSNP=perc_eSNP)%>%
  pivot_longer(cols =c(CTRL,dsbSNP,RDM,eSNP), values_to= "Percentage"
  )%>%dplyr::select(-c(perc_n_no_dist, n_mean,mean_tot))%>%dplyr::select(-name.x)%>%unique()

tsv_count<-annotated_dff[annotated_dff$name.y %in% perc_tsv$name.y,]

tsv_count<-tsv_count%>%dplyr::select(name.y,n_mean,n_ctrl_mean,ctrl_type)%>%unique()%>%dplyr::rename(dsbSNP=n_mean)%>%
  left_join(.,eSNP%>%dplyr::select(name,n_mean)%>%unique()%>%dplyr::rename(eSNP=n_mean)%>%dplyr::rename(name.y=name))%>%
  pivot_wider(names_from = ctrl_type,values_from = n_ctrl_mean)%>%#dplyr::select(-c(ENH,SE,PROM))%>%
  pivot_longer(c(dsbSNP, eSNP ,   RDM , CTRL),values_to = "Count")

tsv_tot<-annotated_dff[annotated_dff$name.y %in% perc_tsv$name.y,]


tsv_tot<-tsv_tot%>%dplyr::select(name.y,mean_tot,mean_tot_control,ctrl_type)%>%unique()%>%dplyr::rename(dsbSNP=mean_tot)%>%
  left_join(.,eSNP%>%dplyr::select(name,mean_tot)%>%unique()%>%dplyr::rename(eSNP=mean_tot)%>%dplyr::rename(name.y=name))%>%
  pivot_wider(names_from = ctrl_type,values_from = mean_tot_control)%>%#dplyr::select(-c(ENH,SE,PROM))%>%
  pivot_longer(c(dsbSNP, eSNP ,   RDM , CTRL),values_to = "tot")


left_join(tsv_count%>%unique(),tsv_tot%>%unique())%>%left_join(.,perc_tsv%>%unique())%>%
  dplyr:: rename(fdr_over_CTRL=fdr)%>%dplyr::rename(OR_over_CTRL=OR)%>% #rename(class=name.x)%>% 
  dplyr::rename(diseases=name.y)%>%pivot_wider(values_from = c(Count,Percentage),names_from = name)



tsv_count <- tsv_count %>%  dplyr::rename(diseases = name.y)
tsv_tot   <- tsv_tot   %>%  dplyr::rename(diseases = name.y)
perc_tsv  <- perc_tsv  %>%  dplyr::rename(diseases = name.y)
left_join(tsv_count, tsv_tot, by = c("diseases", "name")) %>%
  left_join(perc_tsv, by = c("diseases", "name"))


tsv<-left_join(tsv_count%>%unique(),tsv_tot%>%unique())%>%left_join(.,perc_tsv%>%unique())%>%
  dplyr::rename(fdr_over_CTRL=fdr)%>%dplyr::rename(OR_over_CTRL=OR)%>% #rename(class=name.x)%>% 
  pivot_wider(values_from = c(Count,Percentage,tot),names_from = name)









write_tsv(tsv,"results/allelic_imbalance/only_signif_beautifulltsv_fisher.tsv")

###all even non signinf
# #no filter on fdr

plot_df<- annotated_dff%>%filter(ctrl_type=="RDM")%>%dplyr::select(OR,p,fdr,name.x,name.y,perc_n,perc_c)%>%#filter(fdr<0.05)%>%
  
  dplyr:: rename(perc_n_no_dist=perc_n)%>%dplyr::rename(perc_c_no_dist=perc_c)  %>%
  #filter(OR>2)%>%
  unique()%>%#pivot_wider(names_from = name.x,values_from = value)%>%
  left_join(. , annotated_dff%>%filter(ctrl_type=="CTRL")%>%dplyr::select(name.x,name.y,perc_n,perc_c))
#top=20
perc_tsv<-left_join(plot_df,eS)%>%#mutate(perc_c=perc_c/10+1)%>%mutate(perc_n=perc_n+1)%>%
  dplyr::rename(dsbSNP=perc_n)%>% dplyr::rename(CTRL=perc_c)%>%dplyr::rename(RDM=perc_c_no_dist)%>%dplyr::rename(eSNP=perc_eSNP)%>%
  pivot_longer(cols =c(CTRL,dsbSNP,RDM,eSNP), values_to= "Percentage"
  )%>%dplyr::select(-c(perc_n_no_dist, n_mean,mean_tot))%>%dplyr::select(-name.x)%>%unique()

tsv_count<-annotated_dff[annotated_dff$name.y %in% perc_tsv$name.y,]

tsv_count<-tsv_count%>%dplyr::select(name.y,n_mean,n_ctrl_mean,ctrl_type)%>%unique()%>%dplyr::rename(dsbSNP=n_mean)%>%
  left_join(.,eSNP%>%dplyr::select(name,n_mean)%>%unique()%>%dplyr::rename(eSNP=n_mean)%>%dplyr::rename(name.y=name))%>%
  pivot_wider(names_from = ctrl_type,values_from = n_ctrl_mean)%>%#dplyr::select(-c(ENH,SE,PROM))%>%
  pivot_longer(c(dsbSNP, eSNP ,   RDM , CTRL),values_to = "Count")

tsv_tot<-annotated_dff[annotated_dff$name.y %in% perc_tsv$name.y,]


tsv_tot<-tsv_tot%>%dplyr::select(name.y,mean_tot,mean_tot_control,ctrl_type)%>%unique()%>%dplyr::rename(dsbSNP=mean_tot)%>%
  left_join(.,eSNP%>%dplyr::select(name,mean_tot)%>%unique()%>%dplyr::rename(eSNP=mean_tot)%>%dplyr::rename(name.y=name))%>%
  pivot_wider(names_from = ctrl_type,values_from = mean_tot_control)%>%#dplyr::select(-c(ENH,SE,PROM))%>%
  pivot_longer(c(dsbSNP, eSNP ,   RDM , CTRL),values_to = "tot")


left_join(tsv_count%>%unique(),tsv_tot%>%unique())%>%left_join(.,perc_tsv%>%unique())%>%
  dplyr::rename(fdr_over_CTRL=fdr)%>%dplyr::rename(OR_over_CTRL=OR)%>% #rename(class=name.x)%>%
  dplyr:: rename(diseases=name.y)%>%pivot_wider(values_from = c(Count,Percentage),names_from = name)


tsv<-tsv_count%>%left_join(.,perc_tsv%>%unique())%>%
  dplyr::rename(fdr_over_CTRL=fdr)%>%dplyr::rename(OR_over_CTRL=OR)%>% #rename(class=name.x)%>%
  pivot_wider(values_from = c(Count,Percentage),names_from = name)

tsv_final<-tsv%>% left_join(.,tsv_tot%>%mutate(name=paste0("tot_",name))%>%
                              pivot_wider(names_from = name,values_from = tot))

ncol(tsv_final)
tsv_final<-tsv_final[,c(1:4,5:8,13:16,9:12)]

write_tsv(tsv_final,"results/allelic_imbalance/All_diseases_beautifulltsv_fisher.tsv")

tsv_final_EFO<-left_join(tsv_final, annotated_dff%>%select(name.y,EFO)%>%unique())

write_tsv(tsv_final_EFO,"results/allelic_imbalance/All_diseases_beautifulltsv_fisher_EFO.tsv")




#################################################################################
#assemble in a nice tsv 
#this part must be run after the lambda computation
#############################################################################


fisher<-read_tsv("results/allelic_imbalance/All_diseases_beautifulltsv_fisher_EFO.tsv")
lambda<-read_tsv("results/allelic_imbalance/slurm_outputs/lambda/dsbSNP_eSNP_CTRL_RDM_all_diseases.tsv")
#fisher<-left_join(fisher,annotated_dff%>%dplyr::select(name.y,EFO)%>%unique()%>%dplyr::rename(diseases=name.y))



mean_lambda<-lambda%>%na.omit()%>%dplyr::select(-diseases)%>%pivot_longer(cols = starts_with("lambda_gc"),)%>%dplyr::group_by(EFO,name)%>%
  summarise(mean_lambda=mean(value))%>%
  mutate(name=paste0("Mean_",name))%>%
  pivot_wider(names_from = name,values_from = mean_lambda)
ultime_tsv<-left_join(fisher,mean_lambda)

library(httr)
library(jsonlite) 
get_disease_label <- function(code) {
  ontology <- "efo"
  iri <- paste0("http://www.ebi.ac.uk/efo/", code)
  
  url <- paste0("https://www.ebi.ac.uk/ols/api/ontologies/", ontology,
                "/terms?iri=", URLencode(iri, reserved = TRUE))
  
  res <- GET(url)
  if (status_code(res) == 200) {
    data <- fromJSON(content(res, "text", encoding = "UTF-8"))
    if (length(data[["_embedded"]]$terms) > 0) {
      return(data[["_embedded"]]$terms$label[1])
    }
  }
  return(NA)
}

#get EFO and corresponding diseases label
ultime_tsv <- ultime_tsv %>%
  mutate(code = str_extract(EFO, "EFO_[0-9]+"))

unique_codes <- unique(ultime_tsv$code)
lookup <- tibble(
  code = unique_codes,
  disease_label = sapply(unique_codes, get_disease_label)
)

ultime_tsv <- ultime_tsv %>%
  left_join(lookup, by = "code")

ultime_tsv %>% select(EFO, disease_label)

ultime_tsv$disease_label[is.na(ultime_tsv$disease_label)]=ultime_tsv$name.y[is.na(ultime_tsv$disease_label)]


#ultime_tsv<-left_join(fisher,lambda%>%select(diseases,EFO)%>%unique())%>%dplyr::select(-name.y)%>%unique()
df_unique <- ultime_tsv %>%
  distinct()

##name columns

df_unique <- df_unique %>%
  dplyr:: rename(
    P=p,
    FDR_P=fdr_over_CTRL,
    OR=OR_over_CTRL,
    n_dsbSNP = Count_dsbSNP,
    n_eSNP = Count_eSNP,
    `n_CTRL(no_dist)` = Count_RDM,
    `n_ctrlSNP(200-500kb)` = Count_CTRL,
    n_tot_dsbSNP = tot_dsbSNP,
    n_tot_eSNP = tot_eSNP,
    `n_tot_CTRL(no_dist)` = tot_RDM,
    `n_tot_ctrlSNP(200-500kb)` = tot_CTRL,
    `Percentage_CTRL(no_dist)` = Percentage_RDM,
    `Percentage_ctrlSNP(200-500kb)` = Percentage_CTRL,
    `Mean_lambda_gc_dsbSNP` = Mean_lambda_gc_dsbSNP,
    `Mean_lambda_gc_eSNP` = Mean_lambda_gc_eSNP,
    `Mean_lambda_gc_CTRL(no_dist)` = Mean_lambda_gc_RDM,
    `Mean_lambda_gc_ctrlSNP(200-500kb)` = Mean_lambda_gc_CTRL
  )

df_final <- df_unique %>%
  select(
    disease_label,          # label de la maladie
    EFO,
    OR,           # OR
    P,
    FDR_P,          # FDR
    
    # Comptages dans l'ordre dsbSNP, CTRL(no_dist), CTRL(200-500kb), eSNP
    n_dsbSNP, `n_CTRL(no_dist)`, `n_ctrlSNP(200-500kb)`, n_eSNP,   
    # Totaux
    n_tot_dsbSNP, `n_tot_CTRL(no_dist)`, `n_tot_ctrlSNP(200-500kb)`, n_tot_eSNP,  
    # Pourcentages
    Percentage_dsbSNP, `Percentage_CTRL(no_dist)`, `Percentage_ctrlSNP(200-500kb)`, Percentage_eSNP,  
    # Lambda
    Mean_lambda_gc_dsbSNP, `Mean_lambda_gc_CTRL(no_dist)`, `Mean_lambda_gc_ctrlSNP(200-500kb)`, Mean_lambda_gc_eSNP#,  
    # Identifiants
    #code,         
  )
library(dplyr)

df_final <- df_final %>%
  dplyr::rename(
    n_CTRL = `n_CTRL(no_dist)`,
    n_tot_CTRL = `n_tot_CTRL(no_dist)`,
    Percentage_CTRL = `Percentage_CTRL(no_dist)`,
    Mean_lambda_gc_CTRL = `Mean_lambda_gc_CTRL(no_dist)`
  )


df_final <- df_final %>%
  dplyr::rename(
    Disease=disease_label,
    # FDR_P = fdr_over_CTRL,
    # OR = OR_over_CTRL,
    
    # dsbSNP
    NumSignif_dsbSNP = n_dsbSNP,
    NumTot_dsbSNP = n_tot_dsbSNP,
    PercSignif_dsbSNP = Percentage_dsbSNP,
    Lambda_dsbSNP = Mean_lambda_gc_dsbSNP,
    
    # CTRL
    NumSignif_ctrlSNP = n_CTRL,
    NumTot_ctrlSNP = n_tot_CTRL,
    PercSignif_ctrlSNP = Percentage_CTRL,
    Lambda_ctrlSNP = Mean_lambda_gc_CTRL,
    
    # CTRL(200-500kb)
    `NumSignif_ctrlSNP(200-500kb)` = `n_ctrlSNP(200-500kb)`,
    `NumTot_ctrlSNP(200-500kb)` = `n_tot_ctrlSNP(200-500kb)`,
    `PercSignif_ctrlSNP(200-500kb)` = `Percentage_ctrlSNP(200-500kb)`,
    `Lambda_ctrlSNP(200-500kb)` = `Mean_lambda_gc_ctrlSNP(200-500kb)`,
    
    # eSNP
    NumSignif_eSNP = n_eSNP,
    NumTot_eSNP = n_tot_eSNP,
    PercSignif_eSNP = Percentage_eSNP,
    Lambda_eSNP = Mean_lambda_gc_eSNP
  ) %>%
  mutate(Disease = stringr::str_to_title(Disease))

df_final

df_final_no_EFO<-df_final%>%select(-EFO)
DSB_path="/media/sauber/Elements/DSBsnp_project_seb/"
setwd(DSB_path)

#1 NA occured during process, fix it with
SNP<-read_tsv("results/allelic_imbalance/slurm_outputs/lambda/lamba_score_RDM_All_col.tsv")
SNP<-SNP%>%dplyr::select(lambda_gc_SNP,lambda_gc_RDM,lambda_gc_CTRL,diseases,EFO)%>%dplyr::rename(lambda_gc_dsbSNP=lambda_gc_SNP)
eSNP<-read_tsv("results/allelic_imbalance/slurm_outputs/lambda/lambda_eSNP_RDM.tsv")
eSNP<-eSNP%>%select(lambda_gc,diseases)%>%dplyr::rename(lambda_gc_eSNP=lambda_gc)
lamb<-left_join(eSNP,SNP)
lamb<-lamb %>%
  pivot_longer(
    cols = starts_with("lambda_gc"),
    names_to = "snp_type",
    values_to = "lambda"
  ) %>%
  mutate(
    snp_type = gsub("lambda_gc_", "", snp_type) 
  )

Na<-df_final$EFO[is.na(df_final$Lambda_dsbSNP)]

df_final[df_final$EFO==Na,]
Na_df<-lamb[lamb$EFO==Na,]%>%group_by(diseases,EFO,snp_type) %>%summarise(lambda=mean(lambda))%>%pivot_wider(names_from =snp_type,values_from = c(lambda))

# Renommer les colonnes de Na_df pour matcher df_final
Na_df_renamed <- Na_df %>%
  dplyr::rename(
    Lambda_dsbSNP_new = dsbSNP,
    Lambda_ctrlSNP_new = CTRL,
    `Lambda_ctrlSNP(200-500kb)_new` = RDM,   # si RDM correspond à CTRL(no_dist)
    Lambda_eSNP_new = eSNP
  )

# Faire le remplacement dans df_final
df_final <- df_final %>%
  left_join(Na_df_renamed %>% select(-diseases), by = "EFO") %>%
  mutate(
    Lambda_dsbSNP = coalesce(Lambda_dsbSNP, Lambda_dsbSNP_new),
    Lambda_ctrlSNP = coalesce(Lambda_ctrlSNP, Lambda_ctrlSNP_new),
    `Lambda_ctrlSNP(200-500kb)` = coalesce(`Lambda_ctrlSNP(200-500kb)`, `Lambda_ctrlSNP(200-500kb)_new`),
    Lambda_eSNP = coalesce(Lambda_eSNP, Lambda_eSNP_new)
  ) %>%
  select(-ends_with("_new"))%>%select(-diseases)


df_final <- df_final %>%
  mutate(
    EFO = str_replace(EFO, "_", ":"),        # remplacer le premier _ par :
    EFO = str_remove(EFO, "\\.h\\.tsv\\.gz$")  # enlever la terminaison .h.tsv.gz
  )

####remove class and duplicates
df_final<-df_final%>%
  filter(Disease!= "Connective Tissue Disease")%>%#class
  filter(Disease!= "Metabolic Disease")%>%#class
  filter(Disease!= "Type 1 Diabetes (Spa Correction)")%>%#duplicate
  filter(Disease!= "Endocrine System Disease")#duplicate
#filter(Disease!= "Type 1 Diabetes (Spa Correction)")#duplicate


write_tsv(df_final,"results/allelic_imbalance/All_diseases_final_table.tsv")
write_tsv(ultime_tsv_noEFO,"results/allelic_imbalance/All_diseases_final_table_no_EFO.tsv")


library(writexl)

# Export with EFo
write_xlsx(df_final, "results/allelic_imbalance/All_diseases_final_table.xlsx")

# Export without EFO
write_xlsx(df_final_no_EFO, "results/allelic_imbalance/All_diseases_final_table_no_EFO.xlsx")


