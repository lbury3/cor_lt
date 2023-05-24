library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(readr)
library(msigdbr)
library(stringr)
library(ggplot2)
# choose your msigdb collection of interest
hs_gsea_c5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5",
                      subcategory = "BP"
) %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
hs_gsea_c5<- hs_gsea_c5 %>% mutate(gs_name = str_replace_all(gs_name, 'GOBP_', ''))

#for  Cycling
y_n = read_tsv("ranks_EN_deep_w10_correction_yellow_vs_nocolor_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
y_n.gsea <- y_n$rank
names(y_n.gsea) <- as.character(y_n$GeneName)
y_n.gsea <- sort(y_n.gsea, decreasing = T)

#for  IPC
g_n = read_tsv("ranks_EN_deep_w10_correction_green_vs_nocolor_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
g_n.gsea <- g_n$rank
names(g_n.gsea) <- as.character(g_n$GeneName)
g_n.gsea <- sort(g_n.gsea, decreasing = T)

#for  EN_upper_w10
r_n = read_tsv("ranks_EN_deep_correction_red_vs_nocolor_w10_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
r_n.gsea <- r_n$rank
names(r_n.gsea) <- as.character(r_n$GeneName)
r_n.gsea <- sort(r_n.gsea, decreasing = T)

#for  EN_upper_w10
y_r = read_tsv("ranks_EN_deep_correction_yellow_vs_red_w10_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
y_r.gsea <- y_r$rank
names(y_r.gsea) <- as.character(y_r$GeneName)
y_r.gsea <- sort(y_r.gsea, decreasing = T)

#for  EN_upper_w10
g_r = read_tsv("ranks_EN_deep_w10_correction_green_vs_red_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
g_r.gsea <- g_r$rank
names(g_r.gsea) <- as.character(g_r$GeneName)
g_r.gsea <- sort(g_r.gsea, decreasing = T)

#for  EN_upper_w10
y_g = read_tsv("ranks_EN_deep_w10_correction_yellow_vs_green_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
y_g.gsea <- y_g$rank
names(y_g.gsea) <- as.character(y_g$GeneName)
y_g.gsea <- sort(y_g.gsea, decreasing = T)




inputList <- list(
  "Yellow vs no color" = y_n.gsea,
  "Green vs no color" = g_n.gsea,
  "Red vs no color" = r_n.gsea,
  "Yellow vs Red" = y_r.gsea,
  "Green vs Red" = g_r.gsea,
  "Yellow vs Green" = y_g.gsea
  
)


set.seed(42)
test.out<- compareCluster(geneCluster=inputList, fun="GSEA", pvalueCutoff=0.05,
                          pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                          minGSSize = 10, maxGSSize = 500,
                          seed = T)
saveRDS(test.out, "w10_correction_EN_deep_6pairwise_correction_comparecluster_10_500_0.05_3.10.23.rds")
test.out <- readRDS("w10_arch_EN_deep_6pairwise_arch_comparecluster_10_500_0.05_3.10.23.rds")
write.csv(test.out, "w10_arch_EN_deep_6pairwise_arch_comparecluster_10_500_0.05_3.10.23.csv")

selected_GO <- c("SYNAPTIC_SIGNALING", 
               "RIBOSOME_BIOGENESIS",
               "ACTION_POTENTIAL",
               "SYNAPTIC_TRANSMISSION_GLUTAMATERGIC",
               "REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
               "CYTOPLASMIC_TRANSLATION",
               "NEURONAL_ACTION_POTENTIAL",
               "SYNAPTIC_TRANSMISSION_GLUTAMATERGIC",
               "ATP_METABOLIC_PROCESS",
               "MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN",
               "ELECTRON_TRANSPORT_CHAIN",
               "MITOCHONDRIAL_TRANSLATION",
               "MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
               "ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
               "REGULATION_OF_SYNAPTIC_PLASTICITY"
               
)


dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60, includeAll = T) + ggplot2::facet_grid(.sign~.,scales="free") + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) 



