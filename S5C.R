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

#for  yellow vs no clor
y_n = read_tsv("ranks_EN_upper_w16_clay_yellow_vs_nocolor_recorrect_umi_F_sct_cca_wilcox_lfc_new.rnk")
y_n.gsea <- y_n$rank
names(y_n.gsea) <- as.character(y_n$GeneName)
y_n.gsea <- sort(y_n.gsea, decreasing = T)



#for  red vs no color
r_n = read_tsv("ranks_EN_upper_w16_clay_red_vs_nocolor_recorrect_umi_F_sct_cca_wilcox_lfc_new.rnk")
r_n.gsea <- r_n$rank
names(r_n.gsea) <- as.character(r_n$GeneName)
r_n.gsea <- sort(r_n.gsea, decreasing = T)

#for  yellow vs red
y_r = read_tsv("ranks_EN_upper_w16_clay_yellow_vs_red_recorrect_umi_F_sct_cca_wilcox_lfc_new.rnk")
y_r.gsea <- y_r$rank
names(y_r.gsea) <- as.character(y_r$GeneName)
y_r.gsea <- sort(y_r.gsea, decreasing = T)








inputList <- list(
  "Yellow vs no color" = y_n.gsea,
  "Red vs no color" = r_n.gsea,
  "Yellow vs Red" = y_r.gsea

)


set.seed(42)
test.out<- compareCluster(geneCluster=inputList, fun="GSEA", pvalueCutoff=0.05,
                          pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                          minGSSize = 10, maxGSSize = 500,
                          seed = T)
saveRDS(test.out, "w16_Clay_EN_upper_3pairwise_comparecluster_10_500_0.05_4.5.23.rds")
test.out <- readRDS("w16_Clay_EN_upper_3pairwise_comparecluster_10_500_0.05_4.5.23.rds")

#week 16_clay
selected_GO <- c("SYNAPTIC_SIGNALING", 
                 "REGULATION_OF_CHROMATIN_ORGANIZATION",
               "REGULATION_OF_SYNAPTIC_PLASTICITY",
               "G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",
               "REGULATION_OF_MEMBRANE_POTENTIAL",
               "NEURON_DEATH",
               "DNA_CONFORMATION_CHANGE",
               "CHROMATIN_ORGANIZATION",
               "ATP_METABOLIC_PROCESS",
               "ACTION_POTENTIAL",
               "REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
               "CYTOPLASMIC_TRANSLATION",
               "NEURONAL_ACTION_POTENTIAL",
               "ATP_METABOLIC_PROCESS",
           
               "REGULATION_OF_SYNAPTIC_PLASTICITY"
               
)

dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60, includeAll = T) + ggplot2::facet_grid(.sign~.,scales="free") + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) 



