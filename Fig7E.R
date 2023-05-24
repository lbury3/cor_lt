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
# NO color
nocolor = read_tsv("ranks_EN_deep_nocolor_arch_vs_correction_w10_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
nocolor.gsea <- nocolor$rank
names(nocolor.gsea) <- as.character(nocolor$GeneName)
nocolor.gsea <- sort(nocolor.gsea, decreasing = T)

#for  green
green = read_tsv("ranks_EN_deep_green_arch_vs_correction_w10_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
green.gsea <- green$rank
names(green.gsea) <- as.character(green$GeneName)
green.gsea <- sort(green.gsea, decreasing = T)

#for  red
red = read_tsv("ranks_EN_deep_red_arch_vs_correction_w10_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
red.gsea <- red$rank
names(red.gsea) <- as.character(red$GeneName)
red.gsea <- sort(red.gsea, decreasing = T)

#for  yellow
yellow = read_tsv("ranks_EN_deep_yellow_arch_vs_correction_w10_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
yellow.gsea <- yellow$rank
names(yellow.gsea) <- as.character(yellow$GeneName)
yellow.gsea <- sort(yellow.gsea, decreasing = T)


inputList <- list(
  "no color" = nocolor.gsea,
  "green" = green.gsea,
  "red" = red.gsea,
  "yellow" = yellow.gsea
)


set.seed(42)
test.out<- compareCluster(geneCluster=inputList, fun="GSEA", pvalueCutoff=0.05,
                          pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                          minGSSize = 10, maxGSSize = 500,
                          seed = T)
saveRDS(test.out, "w10_EN_deep_arch_corr_noBRAIN_color_cellbender_comparecluster_3.5.2023.rds")
test.out <- readRDS("EN_deep_w10_ALL_arch_vs_correction_recorrect_umi_T.rds")



selected_GO <- c("OXIDATIVE_PHOSPHORYLATION",
                 "CHROMATIN_REMODELING",
                 "RNA_SPLICING",
                 "SYNAPTIC_SIGNALING",
                 "REGULATION_OF_MEMBRANE_POTENTIAL",
                 "HISTONE_H3_K36_METHYLATION",
                 "CHROMOSOME_SEGREGATION",
                 "RETINOIC_ACID_RECEPTOR_SIGNALING_PATHWAY",
                 "REGULATION_OF_TELOMERASE_ACTIVITY",
                 "SYNAPSE_ASSEMBLY",
                 "POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION",
                 "G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY",
                 "GAMMA_AMINOBUTYRIC_ACID_SIGNALING_PATHWAY",
                 "CYTOPLASMIC_TRANSLATION",
                 "INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
                 "REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS")





dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60, includeAll = T) + ggplot2::facet_grid(.sign~.,scales="free") + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) 



