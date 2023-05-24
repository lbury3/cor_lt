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

#for  EN_upper_nocolor
EN_upper_no = read_tsv("ranks_EN_upper_nocolor_arch_vs_correction_w20_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
EN_upper_no.gsea <- EN_upper_no$rank
names(EN_upper_no.gsea) <- as.character(EN_upper_no$GeneName)
EN_upper_no.gsea <- sort(EN_upper_no.gsea, decreasing = T)

#for  EN_upper_red
EN_upper_red = read_tsv("ranks_EN_upper_red_arch_vs_correction_w20_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
EN_upper_red.gsea <- EN_upper_red$rank
names(EN_upper_red.gsea) <- as.character(EN_upper_red$GeneName)
EN_upper_red.gsea <- sort(EN_upper_red.gsea, decreasing = T)

#for  EN_upper_green
EN_upper_green = read_tsv("ranks_EN_upper_green_arch_vs_correction_w20_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
EN_upper_green.gsea <- EN_upper_green$rank
names(EN_upper_green.gsea) <- as.character(EN_upper_green$GeneName)
EN_upper_green.gsea <- sort(EN_upper_green.gsea, decreasing = T)

#for  EN_deep_yellow
EN_upper_yellow = read_tsv("ranks_EN_upper_yellow_arch_vs_correction_w20_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
EN_upper_yellow.gsea <- EN_upper_yellow$rank
names(EN_upper_yellow.gsea) <- as.character(EN_upper_yellow$GeneName)
EN_upper_yellow.gsea <- sort(EN_upper_yellow.gsea, decreasing = T)


inputList <- list(
 
  "no color" = EN_upper_no.gsea,
  "green" = EN_upper_green.gsea,
 "red" = EN_upper_red.gsea,
 "yellow" = EN_upper_yellow.gsea
)


set.seed(42)
library(clusterProfiler)
test.out<- compareCluster(geneCluster=inputList, fun="GSEA", pvalueCutoff=0.05,
                          pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                          minGSSize = 10, maxGSSize = 500,
                          seed = T)
saveRDS(test.out, "w20_EN_upper_LTS_color_arch_vs_correction_cellbender_comparecluster_10_500_0.05_3.10.23.rds")

test.out <- readRDS("w20_EN_upper_LTS_color_arch_vs_correction_cellbender_comparecluster_10_500_0.05_3.10.23.rds")
dotplot(test.out, showCategory=4, split=".sign", label_format = 60, includeAll = T) + facet_grid(.sign~.,scales="free") + theme(axis.text.x=element_text(angle=45, hjust=1)) 
write.csv(test.out, "w20_EN_upper_LTS_color_arch_vs_correction_cellbender_comparecluster_10_500_0.05_3.10.23.csv")


selected_GO <- c("CYTOPLASMIC_TRANSLATION",
                 "NEURON_PROJECTION_GUIDANCE",
                 "NEUROPEPTIDE_SIGNALING_PATHWAY",
                 "ATP_METABOLIC_PROCESS",
                 "SYNAPTIC_SIGNALING",
                 "REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
                 "NEURON_MIGRATION",
                 "REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY",
                 "REGULATION_OF_RNA_SPLICING",
                 "RNA_SPLICING",
                 "ADRENERGIC_RECEPTOR_SIGNALING_PATHWAY",
                 "ADENYLATE_CYCLASE_ACTIVATING_ADRENERGIC_RECEPTOR_SIGNALING_PATHWAY",
                 "ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY"
)



dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60, includeAll = T) + facet_grid(.sign~.,scales="free") + theme(axis.text.x=element_text(angle=45, hjust=1)) 



