library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(readr)
library(msigdbr)
library(stringr)
library(ggplot2)
renv::restore(project = "~/SCP_env")

# choose your msigdb collection of interest
hs_gsea_c5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5",
                      subcategory = "BP"
) %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
hs_gsea_c5<- hs_gsea_c5 %>% mutate(gs_name = str_replace_all(gs_name, 'GOBP_', ''))

#for  Cycling
y_n = read_tsv("ranks_EN_upper_w20_correction_yellow_vs_nocolor_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
y_n.gsea <- y_n$rank
names(y_n.gsea) <- as.character(y_n$GeneName)
y_n.gsea <- sort(y_n.gsea, decreasing = T)

#for  IPC
g_n = read_tsv("ranks_EN_upper_w20_correction_green_vs_nocolor_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
g_n.gsea <- g_n$rank
names(g_n.gsea) <- as.character(g_n$GeneName)
g_n.gsea <- sort(g_n.gsea, decreasing = T)

#for  EN_upper_w10
r_n = read_tsv("ranks_EN_upper_w20_correction_red_vs_nocolor_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
r_n.gsea <- r_n$rank
names(r_n.gsea) <- as.character(r_n$GeneName)
r_n.gsea <- sort(r_n.gsea, decreasing = T)

#for  EN_upper_w10
y_r = read_tsv("ranks_EN_upper_w20_correction_yellow_vs_red_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
y_r.gsea <- y_r$rank
names(y_r.gsea) <- as.character(y_r$GeneName)
y_r.gsea <- sort(y_r.gsea, decreasing = T)

#for  EN_upper_w10
g_r = read_tsv("ranks_EN_upper_w20_correction_green_vs_red_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
g_r.gsea <- g_r$rank
names(g_r.gsea) <- as.character(g_r$GeneName)
g_r.gsea <- sort(g_r.gsea, decreasing = T)

#for  EN_upper_w10
y_g = read_tsv("ranks_EN_upper_w20_correction_yellow_vs_green_recorrect_umi_F_sct_cca_wilcox_lfc.rnk")
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
saveRDS(test.out, "w20_correction_EN_upper_6pairwise_comparecluster_10_500_0.05_3.10.23.rds")
test.out <- readRDS("w20_arch_EN_upper_6pairwise_comparecluster_10_500_0.05_3.10.23.rds")

write.csv(test.out, "w20_arch_EN_upper_6pairwise_comparecluster_10_500_0.05_3.10.23.csv")
selected_GO <- c("SYNAPTIC_SIGNALING", 
                 "G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",
               "REGULATION_OF_MEMBRANE_POTENTIAL",
               "REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",
               "ACTION_POTENTIAL",
               "SYNAPTIC_MEMBRANE_ADHESION",
               "REGULATION_OF_CELL_ADHESION",
               "REGULATION_OF_JNK_CASCADE",
               "JNK_CASCADE",
               "SEMAPHORIN_PLEXIN_SIGNALING_PATHWAY",
               "REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
               "CYTOPLASMIC_TRANSLATION",
               "NEURONAL_ACTION_POTENTIAL",
               "ATP_METABOLIC_PROCESS",
          
               "REGULATION_OF_SYNAPTIC_PLASTICITY"
               
)

dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60, includeAll = T) + ggplot2::facet_grid(.sign~.,scales="free") + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) 



