EN_upper_w20 = read_tsv("ranks_EN_upper_all_arch_vs_correction_w20_recorrect_umi_T_sct_cca_wilcox_lfc.rnk")
EN_upper_w20.gsea <- EN_upper_w20$rank
names(EN_upper_w20.gsea) <- as.character(EN_upper_w20$GeneName)
EN_upper_w20.gsea <- sort(EN_upper_w20.gsea, decreasing = T)

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
test.out<- GSEA( EN_upper_w20.gsea, pvalueCutoff=0.05,
                 pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                 minGSSize = 10, maxGSSize = 500,
                 seed = T)


dotplot(test.out, showCategory=5, split=".sign", label_format = 60) + ggplot2::facet_grid(.sign~.,scales="free") + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1))
