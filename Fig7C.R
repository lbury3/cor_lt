EN_deep_w10 = read_tsv("ranks_EN_deep_ALL_arch_vs_correction_w10_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
EN_deep_w10.gsea <- EN_deep_w10$rank
names(EN_deep_w10.gsea) <- as.character(EN_deep_w10$GeneName)
EN_deep_w10.gsea <- sort(EN_deep_w10.gsea, decreasing = T)

library(msigdbr)
library(stringr)
# choose your msigdb collection of interest
hs_gsea_c5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5",
                      subcategory = "BP"
) %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
hs_gsea_c5<- hs_gsea_c5 %>% mutate(gs_name = str_replace_all(gs_name, 'GOBP_', ''))
set.seed(42)
#library(clusterProfiler)
test.out<- GSEA( EN_deep_w10.gsea, pvalueCutoff=0.05,
                 pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                 minGSSize = 10, maxGSSize = 500,
                 seed = T)

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




dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60) + ggplot2::facet_grid(.sign~.,scales="free") + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1))
