# Load required libraries
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(viridis)
library("pheatmap")
library("ggplot2")
library("ggdendro")
library("ggalign")
library("tidyr")
library("grid")
library(httr)
library(jsonlite)
library(stringr)
library(dplyr)

# Set path and working directory
DSB_path="/media/sauber/Elements/DSBsnp_project_seb/results/allelic_imbalance/"
setwd(DSB_path)

# Number of top associated EFOs to use before trimming
top=31 

# Read lambda GC scores for SNPs
SNP<-read_tsv("slurm_outputs/lambda/lamba_score_RDM_All_col.tsv")
SNP<-SNP %>%
  dplyr::select(lambda_gc_SNP,lambda_gc_RDM,lambda_gc_CTRL,diseases,EFO) %>%
  rename(lambda_gc_dsbSNP=lambda_gc_SNP)

# Read lambda GC scores for eSNPs
eSNP<-read_tsv("slurm_outputs/lambda/lambda_eSNP_RDM.tsv")
eSNP<-eSNP %>%
  select(lambda_gc,diseases) %>%
  rename(lambda_gc_eSNP=lambda_gc)

# Merge SNP and eSNP data
df<-left_join(eSNP,SNP)

# Pivot the data to long format for plotting
df_f<-df %>%
  pivot_longer(
    cols = starts_with("lambda_gc"),
    names_to = "snp_type",
    values_to = "lambda"
  ) %>%
  mutate(
    snp_type = gsub("lambda_gc_", "", snp_type) 
  )

# Load top 31 EFO diseases
load("top31_EFO_diseases_fdr_1e-08.RDATA")
top20_list_fdr_EFO

# Filter top 20 diseases by FDR and reorder
df_top20 <- df_f %>%
  filter(EFO %in% top20_list_fdr_EFO) %>%
  mutate(
    EFO = factor(EFO, levels = top20_list_fdr_EFO),
    diseases = forcats::fct_reorder(diseases, as.numeric(EFO))
  )

# Calculate mean lambda GC per EFO and SNP type
df_mean_lambda <- df_top20 %>%
  mutate(disease_name=diseases) %>%
  group_by(EFO, snp_type) %>%
  summarise(mean_lambda = mean(lambda, na.rm = TRUE), .groups = "drop")

# Function to get disease labels from EFO codes
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

# Map EFO codes to disease labels
df_top20_mean <- df_mean_lambda %>%
  mutate(code = str_extract(EFO, "EFO_[0-9]+"))

unique_codes <- unique(df_top20_mean$code)
lookup <- tibble(
  code = unique_codes,
  disease_label = sapply(unique_codes, get_disease_label)
)

df_top20_mean <- df_top20_mean %>%
  left_join(lookup, by = "code")

df_top20_mean %>% select(EFO, disease_label, snp_type, mean_lambda)

# Fill missing disease labels
df_top20_mean$disease_name=df_top20_mean$disease_label
EFO_na_list<-df_top20_mean$EFO[is.na(df_top20_mean$disease_name)]

df_top20_mean <- df_top20_mean %>%
  left_join(
    df_f %>%
      filter(EFO %in% EFO_na_list) %>%
      select(diseases, EFO) %>%
      distinct() %>%
      group_by(EFO) %>%
      slice_head(n = 1),
    by = "EFO"
  ) %>%
  mutate(diseases = sub("_.*", "", diseases)) %>%
  mutate(
    disease_name  = if_else(is.na(disease_name), diseases, disease_name),
    disease_label = if_else(is.na(disease_label), diseases, disease_label)
  )

# Load top EFO FDR table for plotting OR
load(file=paste0("table_fdr_top31_EFO_diseases_fdr_1e-08.RDATA"))
df_fdr<-order_table

# Merge lambda GC with OR and FDR
plot_OR<-left_join(df_top20_mean%>%select(-diseases)%>%unique(), df_fdr%>%select(EFO,OR,fdr)) %>% unique()

# Plot mean lambda GC ordered by OR
p__lamnda_OR<-plot_OR %>%
  mutate(snp_type = factor(snp_type, levels = c("dsbSNP", "RDM","CTRL", "eSNP"))) %>%
  ggplot(aes(x = reorder(disease_label,-OR), y = mean_lambda, fill = snp_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("CTRL" = "blue", "RDM" = "lightblue", "dsbSNP" = "red", "eSNP" = "violet")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  coord_cartesian(ylim = c(1, max(df_top20_mean$mean_lambda) * 1.05)) +
  ylab("Mean Lambda GC") +
  xlab("Diseases") +
  ggtitle("Top 20 Associated Diseases by SNP Type (ordered by OR)")

# Plot mean lambda GC ordered by FDR
p_lambda_fdr<-plot_OR %>%
  mutate(snp_type = factor(snp_type, levels = c("dsbSNP", "RDM","CTRL",  "eSNP"))) %>%
  ggplot(aes(x = reorder(disease_label,fdr), y = mean_lambda, fill = snp_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("CTRL" = "blue", "RDM" = "lightblue", "dsbSNP" = "red", "eSNP" = "violet")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  coord_cartesian(ylim = c(1, max(df_top20_mean$mean_lambda) * 1.05)) +
  ylab("Mean Lambda GC") +
  xlab("Diseases") +
  ggtitle("Top 20 Associated Diseases by SNP Type (ordered by FDR)")

# Save plots to PDF
pdf(paste0("slurm_outputs/lambda/top_",top,"_lambda_fused_tables_count_1_1e-08_dsbSNP_CTRL_RDM_eSNP.pdf"),9,5)
print(p__lamnda_OR)
print(p_lambda_fdr)
dev.off()

# Save CSV for all diseases
dff<-df
DSB_path="/media/sauber/Elements/DSBsnp_project_seb/"
setwd(DSB_path)

Annotation<- read_tsv("All_harmonised_GWAS_pvalue_under1e-6_above500_gwas_with_traits_annotated_by_class_all_diseaese_currated.tsv")
Annotation$EFO<-Annotation$EFO %>% str_split(":") %>% sapply(paste, collapse ="_") %>% paste0(".h.tsv.gz")
annotated_dff<-merge.data.frame(Annotation,dff,by="EFO")
annotated_dff<-na.omit(annotated_dff)

# Define disease classes
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

# Annotate diseases with class
ouais<-data.frame()
final_annontation_dff<-data.frame()
for (clas in DiseaseList){
  ouais<-(annotated_dff[annotated_dff$name==clas,])
  ouais$class<- rep(clas,nrow(ouais))
  final_annontation_dff<-(rbind(final_annontation_dff,ouais))
}
dff<-dff[dff$EFO %in% final_annontation_dff$EFO,]

print(dff)
dff_tsv<-dff[,c(2,3,1,4,5,6)]
write_tsv(dff_tsv,"results/allelic_imbalance/slurm_outputs/lambda/dsbSNP_eSNP_CTRL_RDM_all_diseases.tsv")

###############################################################################"

#######################plot boxplot

##############################################################################

############################################################################

# Copy the original data
dff_class <- dff_tsv

# Join annotations and reshape data
dff_class <- left_join(dff_class, annotated_dff %>% select(EFO, name) %>% unique()) %>%
  pivot_longer(
    cols = starts_with("lambda_gc"), 
    names_to = "snp_type", 
    values_to = "lambda_gc"
  ) %>%
  mutate(
    snp_type = sub("lambda_gc_", "", snp_type),   # keep only dsbSNP, eSNP, RDM, CTRL
    snp_type = factor(snp_type, levels = c("dsbSNP", "CTRL", "RDM", "eSNP"))  # desired order
  )

# Count number of diseases per gene
class_count <- dff_class %>% select(diseases, name) %>% unique() %>%
  group_by(name) %>%
  summarise(n_maladies = n_distinct(diseases)) %>%
  arrange(desc(n_maladies))

mode = "ALL"
mode = 20
#mode = 30

# Filter by mode
if (mode == "ALL"){ 
  dff_class = dff_class
} else {
  dff_class <- dff_class[dff_class$name %in% class_count$name[class_count$n_maladies >= mode], ]
  dff_class <- dff_class[!dff_class$name %in% c("genetic.disorder", "integumentary.system.disease"), ]
}

unique(dff_class$name)

# Perform Wilcoxon test per group
wilcox_results <- dff_class %>%
  group_by(name, snp_type) %>%
  summarise(
    test = list(wilcox.test(lambda_gc, mu = 1)),
    .groups = "drop"
  )

# Extract test results
wilcox_results <- wilcox_results %>%
  mutate(
    p_value = sapply(test, function(x) x$p.value),
    W       = sapply(test, function(x) x$statistic),
    p_adj   = p.adjust(p_value, method = "fdr")
  ) %>%
  select(name, snp_type, W, p_value, p_adj)

wilcox_results$p_value %>% as.numeric()

# Quick plots per SNP type
left_join(dff_class, wilcox_results) %>% filter(snp_type == "dsbSNP") %>%
  ggplot() + aes(x = reorder(name, p_adj), y = lambda_gc) + geom_violin() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0.9, 1.8))

left_join(dff_class, wilcox_results) %>% filter(snp_type == "eSNP") %>%
  ggplot() + aes(x = reorder(name, p_adj), y = lambda_gc) + geom_violin() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0.9, 1.8))

left_join(dff_class, wilcox_results) %>% filter(snp_type == "RDM") %>%
  ggplot() + aes(x = reorder(name, p_adj), y = lambda_gc) + geom_violin() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0.9, 1.8))

left_join(dff_class, wilcox_results) %>% filter(snp_type == "CTRL") %>%
  ggplot() + aes(x = reorder(name, p_adj), y = lambda_gc) + geom_violin() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0.9, 1.8))

plot_wilcox <- left_join(dff_class, wilcox_results)

# Add dsbSNP p_adj for ordering
plot_wilcox <- left_join(
  plot_wilcox,
  plot_wilcox %>% filter(snp_type == "dsbSNP") %>% mutate(dsb_p_adj = p_adj) %>% select(diseases, name, dsb_p_adj)
)

# Basic violin plot
plot_wilcox %>%
  ggplot() + aes(x = reorder(name, dsb_p_adj), y = lambda_gc, fill = snp_type) + geom_violin() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0.9, 1.8)) +
  scale_fill_manual(values = c("CTRL" = "blue", "RDM" = "lightblue", "dsbSNP" = "red", "eSNP" = "violet"))

# Filter significant dsbSNP
plot_wilcox %>%
  filter(dsb_p_adj < 0.05) %>%
  ggplot() + aes(x = reorder(name, dsb_p_adj), y = lambda_gc, fill = snp_type) + geom_violin() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0.9, 1.8)) +
  scale_fill_manual(values = c("CTRL" = "blue", "RDM" = "lightblue", "dsbSNP" = "red", "eSNP" = "violet"))

# With facets
plot_wilcox %>%
  ggplot() + aes(x = reorder(name, dsb_p_adj), y = lambda_gc, fill = snp_type) + geom_violin() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  coord_cartesian(ylim = c(0.9, 1.8)) +
  facet_wrap(~snp_type, nrow = 2) +
  scale_fill_manual(values = c("CTRL" = "blue", "RDM" = "lightblue", "dsbSNP" = "red", "eSNP" = "violet")) +
  ggtitle(paste0("Diseases class association (n>", mode, ")"))

# Prepare labels
plot_wilcox_labels <- plot_wilcox %>%
  group_by(name, snp_type) %>%
  summarise(
    lambda_max = max(lambda_gc, na.rm = TRUE),
    padj_label = signif(first(p_adj), digits = 3),
    .groups = "drop"
  )





# Height of the band (bottom of the band)
band_low <- 1.30

# Set ymax = Inf so the band reaches the top of the panel
label_y   <- band_low + 0.03   # vertical position of labels (just above the bottom of the band)

plot_class_wilcox_stars <- plot_wilcox_clean %>%
  mutate(snp_type = factor(snp_type, levels = c("dsbSNP","RDM","CTRL","eSNP"))) %>%
  ggplot(aes(x = name_clean, y = lambda_gc, fill = snp_type)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.7) +
  
  # White band that goes up to the top of the panel
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = band_low, ymax = Inf,
           fill = "white", color = NA) +
  
  # Horizontal labels, all at the same y (inside the band)
  geom_text(
    data = plot_wilcox_labels,
    aes(x = name_clean, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 3, angle = 0, hjust = 0.5, vjust = 0.5
  ) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        clip = "off",          # allow annotations outside the panel if needed
        strip.clip = "off",
        plot.margin = margin(t = 8, r = 6, b = 6, l = 6)) +
  facet_wrap(~snp_type, nrow = 2) +
  scale_fill_manual(values = c(
    "dsbSNP" = "red",
    "RDM"    = "lightblue",
    "CTRL"   = "blue",
    "eSNP"   = "violet"
  )) +
  coord_cartesian(ylim = c(0.95, 1.4)) +   # adjust if needed
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ggtitle(paste0("Diseases class association (n>", mode, ")"))

plot_class_wilcox_stars

# Save plot to PDF
pdf(paste0("results/allelic_imbalance/slurm_outputs/lambda/n_over_",mode,"_class_lambda_fused_tables_count_1_1e-08_dsbSNP_CTRL_RDM_eSNP_ordered_by_median_with_white.pdf"),10,8)
print(plot_class_wilcox_stars)
dev.off()


median_dsb <- plot_wilcox %>%
  filter(snp_type == "dsbSNP") %>%
  group_by(name) %>%
  summarise(median_dsb = median(lambda_gc, na.rm = TRUE))

library(dplyr)
library(stringr)

# clean labels
plot_wilcox_clean <- plot_wilcox %>%
  mutate(
    name_clean = gsub("\\.", " ", name),
    name_clean = str_to_sentence(name_clean)
  )

# median
order_dsb <- plot_wilcox_clean %>%
  filter(snp_type == "dsbSNP") %>%
  group_by(name_clean) %>%
  summarise(med_dsb = median(lambda_gc, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_dsb)

# define order
plot_wilcox_clean <- plot_wilcox_clean %>%
  mutate(name_clean = factor(name_clean, levels = order_dsb$name_clean))

plot_wilcox_labels <- plot_wilcox_labels %>%
  mutate(
    name_clean = gsub("\\.", " ", name),
    name_clean = str_to_sentence(name_clean),
    name_clean = factor(name_clean, levels = order_dsb$name_clean)
  )

# Plot
plot_class_wilcox <- plot_wilcox_clean %>%  mutate(
  snp_type = factor(snp_type, levels = c("dsbSNP", "RDM","CTRL", "eSNP"))
) %>%
  ggplot(aes(x = name_clean, y = lambda_gc, fill = snp_type)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.7) +
  geom_text(
    data = plot_wilcox_labels,
    aes(x = name_clean, y = 1.2, label = label),
    inherit.aes = FALSE,
    size = 3,
    angle = 45,
    hjust = 0
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        clip = "off") +
  facet_wrap(~snp_type, nrow = 2) +
  scale_fill_manual(values = c("CTRL" = "blue",
                               "RDM" = "lightblue",
                               "dsbSNP" = "red",
                               "eSNP" = "violet")) +
  coord_cartesian(ylim = c(0.95, 1.4)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  expand_limits(y = max(plot_wilcox_labels$lambda_max) + 0.05) +
  ggtitle(paste0("Diseases class association (n>", mode, ")"))

plot_wilcox_clean
plot_class_wilcox



# Same plot but with 1 row of facets
plot_class_wilcox_stars <- plot_wilcox_clean %>%
  mutate(snp_type = factor(snp_type, levels = c("dsbSNP","RDM","CTRL","eSNP"))) %>%
  ggplot(aes(x = name_clean, y = lambda_gc, fill = snp_type)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.7) +
  
  # White band that goes to the top of the panel
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = band_low, ymax = Inf,
           fill = "white", color = NA) +
  
  # Horizontal labels
  geom_text(
    data = plot_wilcox_labels,
    aes(x = name_clean, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 3, angle = 0, hjust = 0.5, vjust = 0.5
  ) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        clip = "off",          # allow annotations outside the panel
        strip.clip = "off",
        plot.margin = margin(t = 8, r = 6, b = 6, l = 6)) +
  facet_wrap(~snp_type, nrow = 1) +
  scale_fill_manual(values = c(
    "dsbSNP" = "red",
    "RDM"    = "lightblue",
    "CTRL"   = "blue",
    "eSNP"   = "violet"
  )) +
  coord_cartesian(ylim = c(0.95, 1.4)) +   # adjust if needed
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ggtitle(paste0("Diseases class association (n>", mode, ")"))

plot_class_wilcox_stars

# Save larger PDF version
pdf(paste0("results/allelic_imbalance/slurm_outputs/lambda/n_over_",mode,"_class_lambda_fused_tables_count_1_1e-08_dsbSNP_CTRL_RDM_eSNP_ordered_by_median_with_white_large.pdf"),10*1.5,5)
print(plot_class_wilcox_stars)
dev.off()


# Update band for a higher range
band_low <- 2.05
# ymax = Inf so that band goes to top
label_y <- band_low + 0.03   # vertical position of labels (just above the bottom of the band)


plot_class_wilcox_stars2 <- plot_wilcox_clean %>%
  mutate(snp_type = factor(snp_type, levels = c("dsbSNP","RDM","CTRL","eSNP"))) %>%
  ggplot(aes(x = name_clean, y = lambda_gc, fill = snp_type)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.7) +
  
  # White band to top of panel
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = band_low, ymax = Inf,
           fill = "white", color = NA) +
  
  # Horizontal labels inside the band
  geom_text(
    data = plot_wilcox_labels,
    aes(x = name_clean, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 3, angle = 0, hjust = 0.5, vjust = 0.5
  ) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        clip = "off",          # allow annotations outside the panel
        strip.clip = "off",
        plot.margin = margin(t = 8, r = 6, b = 6, l = 6)) +
  facet_wrap(~snp_type, nrow = 1) +
  scale_fill_manual(values = c(
    "dsbSNP" = "red",
    "RDM"    = "lightblue",
    "CTRL"   = "blue",
    "eSNP"   = "violet"
  )) +
  # coord_cartesian(ylim = c(0.95, 1.4)) +   # optional zoom
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ggtitle(paste0("Diseases class association (n>", mode, ")"))

plot_class_wilcox_stars2

# Save unzoomed large version
pdf(paste0("results/allelic_imbalance/slurm_outputs/lambda/n_over_",mode,"_class_lambda_fused_tables_count_1_1e-08_dsbSNP_CTRL_RDM_eSNP_ordered_by_median_with_white_large_unzoomed.pdf"),10*1.5,5)
print(plot_class_wilcox_stars2)
dev.off()

