
library(Matrix)
library(Seurat)
library(irlba)
library(matrixStats)
###
#install.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

arch_20wk_ali.data <- Read10X_h5("arch_20wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")
arch_20wk_ali <- CreateSeuratObject(counts = arch_20wk_ali.data, project = "arch_w20_WT/Q76*")

arch_20wk_ali[["percent.mt"]] <- PercentageFeatureSet(arch_20wk_ali, pattern = "^MT-")
arch_20wk_ali <- subset(arch_20wk_ali, subset = nFeature_RNA > 1250 & nFeature_RNA < 7500 & percent.mt < 15) 

arch_20wk_ali <- SCTransform(arch_20wk_ali, vst.flavor = "v2", vars.to.regress = "percent.mt")


correction_20wk_ali.data <- Read10X_h5("correction_20wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")
correction_20wk_ali <- CreateSeuratObject(counts = correction_20wk_ali.data, project = "arch_w20_WT/WT")
correction_20wk_ali[["percent.mt"]] <- PercentageFeatureSet(correction_20wk_ali, pattern = "^MT-")
correction_20wk_ali <- subset(correction_20wk_ali, subset = nFeature_RNA > 1250 & nFeature_RNA < 7500 & percent.mt < 15) 

correction_20wk_ali <- SCTransform(correction_20wk_ali, vst.flavor = "v2", vars.to.regress = "percent.mt")

merge.list <- list(correction_20wk_ali, arch_20wk_ali)
features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 3004)
features <- features[-grep(pattern = "ZsGreen|TdTomato|Cre|Dre", x = features)]
merge.list <- PrepSCTIntegration(object.list = merge.list, anchor.features = features)
merge.anchors <- FindIntegrationAnchors(object.list = merge.list, normalization.method = "SCT", anchor.features = features)
merge.combined <- IntegrateData(anchorset = merge.anchors, normalization.method = "SCT")
merge.combined <- RunPCA(merge.combined, verbose = FALSE)
ElbowPlot(merge.combined, 50)

merge.combined <- FindNeighbors(merge.combined, dims = 1:40)
merge.combined <- FindClusters(merge.combined, resolution = 1.5)
merge.combined <- RunUMAP(merge.combined, dims = 1:40)

p1 <- DimPlot(merge.combined, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 1)



saveRDS(merge.combined, file = "arch_arch_corr_w20_ali_integrated_4.18.24.rds")

DefaultAssay(object = merge.combined) <- "SCT"
merge.combined<- PrepSCTFindMarkers(merge.combined)

markers <- FindAllMarkers(merge.combined, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(markers, "markers_w20_prepSCT_findallmarkers_arch_archcorr_4.18.24.rds")
write.csv(markers, file = "markers_w20_prepSCT_findallmarkers_arch_archcorr_4.18.24.csv")

# annotating the cell types for each cluster based on marker gene expression. 
library(Seurat)
merge.combined<- readRDS("w20_arch_archcorr_annotated_4.18.24.rds")

Idents(merge.combined) <- "seurat_clusters"

new_ident_1 <- setNames(c("RG/Ast", #0  
                          "RG/Ast",  #1  
                          "DL_neuron",#2 
                          "IPC", #3 RG/Ast_PGK1
                          "UL_neuron", #4
                          "UL_neuron",  #5  
                          "UL_neuron", #6 
                          "DL_neuron", #7
                          "Other", #8
                          "UL_neuron", #9  
                          "IPC", #10 , 
                          "DL_neuron", #11, 
                          "Cycling_NPC", #12 
                          "Other", #13
                          "Cycling_NPC",#14 
                          "IPC",  #15 
                          "Other", #16  
                          "IPC", #17
                          "RG/Ast", #18
                          "Other",#19
                          "DL_neuron",#20
                          "UL_neuron",#21
                          "RG/Ast",#22 Lhx1 and Lhx5 are expressed in CR cells 
                          "Other",#23
                          "RG/Ast"#24
                         
),
levels(merge.combined))
organoid.combined.cca.sct_1_rename <- RenameIdents(merge.combined, new_ident_1)
DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)
Idents(organoid.combined.cca.sct_1_rename) <- "CellType"

organoid.combined.cca.sct_1_rename$CellType <- Idents(organoid.combined.cca.sct_1_rename)
saveRDS(organoid.combined.cca.sct_1_rename, "w20_arch_archcorr_reannotated_4.23.24.rds")
organoid.combined.cca.sct_1_rename<- readRDS("w20_arch_archcorr_reannotated_4.23.24.rds")
DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)
# For lineage proportion calculation deep layer neuron.
# Among lineage DEG generation for Arch WT/WT week 20 upper layer neuron. 
EN_upper <- subset(organoid.combined.cca.sct_1_rename, idents = "UL_neuron")
DefaultAssay(EN_upper) <- "RNA"
Idents(EN_upper) <- "orig.ident"
EN_upper_correction <- subset(EN_upper, idents = "arch_w20_WT/WT")

EN_upper_correction$lineage <- "no_color"  #18
EN_upper_correction$lineage[WhichCells(EN_upper_correction, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow" #335
EN_upper_correction$lineage[WhichCells(EN_upper_correction, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"  #814
EN_upper_correction$lineage[WhichCells(EN_upper_correction, expression = TdTomato == 0 & ZsGreen >0)] <- "green" #1
# total EN_upper_correction  1202
# no color   14
# yellow     473
# red        714
# green      1
# UL_neuron correction  yellow vs red 
Idents(EN_upper_correction) <- "lineage"
DefaultAssay(EN_upper_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_upper_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "red",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_upper_w20_correction_yellow_vs_red_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_correction_yellow_vs_red_recorrect_umi_F.csv")
# UL_neuron correction yellow vs no color
myTopHits.df <- FindMarkers(EN_upper_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_upper_w20_correction_yellow_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_correction_yellow_vs_nocolor_recorrect_umi_F.csv")
# UL_neuron correction red vs no color
myTopHits.df <- FindMarkers(EN_upper_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "red", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_upper_w20_correction_red_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_correction_red_vs_nocolor_recorrect_umi_F.csv")
#cannot do yellow vs green NA   
#cannot do green vs nocolor
#cannot do green vs red



# Among lineage DEG generation for Arch WT/Q76* week 20 upper layer neuron.
organoid.combined.cca.sct_1_rename<- readRDS("w20_arch_archcorr_reannotated_4.23.24.rds")
EN_upper <- subset(organoid.combined.cca.sct_1_rename, idents = "UL_neuron")
DefaultAssay(EN_upper) <- "RNA"
Idents(EN_upper) <- "orig.ident"
EN_upper_arch <- subset(EN_upper, idents = "arch_w20_WT/Q76*")

EN_upper_arch$lineage <- "no_color"
EN_upper_arch$lineage[WhichCells(EN_upper_arch, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_upper_arch$lineage[WhichCells(EN_upper_arch, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_upper_arch$lineage[WhichCells(EN_upper_arch, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
#  EN_upper_arch   1568
# no color          3
# yellow          1163
# red             401
# green           1
### UL_neuron correction  yellow vs red 
Idents(EN_upper_arch) <- "lineage"
DefaultAssay(EN_upper_arch) <- "SCT"
myTopHits.df <- FindMarkers(EN_upper_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "red",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_upper_w20_arch_yellow_vs_red_recorrect_umi_F.rds")
library(tidyverse)
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_arch_yellow_vs_red_recorrect_umi_F.csv")
# UL_neuron correction yellow vs no color
myTopHits.df <- FindMarkers(EN_upper_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_upper_w20_arch_yellow_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_arch_yellow_vs_nocolor_recorrect_umi_F.csv")
# UL_neuron correction red vs no color
myTopHits.df <- FindMarkers(EN_upper_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "red", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_upper_w20_arch_red_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_arch_red_vs_nocolor_recorrect_umi_F.csv")
#these are for Fig S17B.
##### Among lineage DE analysis for week 10 deep layer neurons in Arch WT/WT
organoid.combined.cca.sct_1_rename<- readRDS("w20_arch_archcorr_reannotated_4.23.24.rds")
EN_deep <- subset(organoid.combined.cca.sct_1_rename, idents = "DL_neuron")
DefaultAssay(EN_deep) <- "RNA"
Idents(EN_deep) <- "orig.ident"
EN_deep_correction <- subset(EN_deep, idents = "arch_w20_WT/WT")

EN_deep_correction$lineage <- "no_color"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
  #  EN_deep_correction   1039
  # no color                 22
  # yellow                  115
  # red                     900
  # green                    2

##### Among lineage DE analysis for week 10 deep layer neurons in Arch WT/Q76*

EN_deep_arch <- subset(EN_deep, idents = "arch_w20_WT/Q76*")

EN_deep_arch$lineage <- "no_color"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
  #  EN_deep_arch        1031
# no color                 22
# yellow                  494
# red                     510
# green                    5

##Calculation of 





organoid.combined.cca.sct_1_rename<- readRDS("w20_arch_archcorr_reannotated_4.23.24.rds")
srt_UL_neuron <- subset(organoid.combined.cca.sct_1_rename, idents = "UL_neuron")
DefaultAssay(srt_UL_neuron) <- "RNA"

srt_UL_neuron_Green_yes_red_yes <- subset(srt_UL_neuron,  ZsGreen > 0 &  TdTomato >0) 
saveRDS(srt_UL_neuron_Green_yes_red_yes, "srt_w20_UL_neuron_yellow.rds")



Idents(srt_UL_neuron_Green_yes_red_yes) <- "orig.ident"
DefaultAssay(srt_UL_neuron_Green_yes_red_yes) <- "SCT"
srt_UL_neuron_Green_yes_red_yes<- PrepSCTFindMarkers(srt_UL_neuron_Green_yes_red_yes)
myTopHits.df <- FindMarkers(srt_UL_neuron_Green_yes_red_yes, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "arch_w20_WT/Q76*", ident.2 = "arch_w20_WT/WT",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "UL_neuron_w20_yellow_arch_vs_correction_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "UL_neuron_w20_yellow_arch_vs_correction_recorrect_umi_F.csv")



srt_UL_neuron_Green_no_red_yes <- subset(srt_UL_neuron,  ZsGreen == 0 &  TdTomato >0) 
saveRDS(srt_UL_neuron_Green_no_red_yes, "srt_w20_UL_neuron_red.rds")



Idents(srt_UL_neuron_Green_no_red_yes) <- "orig.ident"
DefaultAssay(srt_UL_neuron_Green_no_red_yes) <- "SCT"
srt_UL_neuron_Green_no_red_yes<- PrepSCTFindMarkers(srt_UL_neuron_Green_no_red_yes)
myTopHits.df <- FindMarkers(srt_UL_neuron_Green_no_red_yes, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "arch_w20_WT/Q76*", ident.2 = "arch_w20_WT/WT",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "UL_neuron_w20_red_arch_vs_correction_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "UL_neuron_w20_red_arch_vs_correction_recorrect_umi_F.csv")

#rnk file generation
#rnk file generation
Chap_p_vs_chap <- readRDS("UL_neuron_w20_yellow_arch_vs_correction_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- Chap_p_vs_chap %>%
  as_tibble(rownames = "geneID")

#create ranks file
ranks_chap_p.chap = Chap_p_vs_chap.df$avg_log2FC
ranks_chap_p.chap <- cbind(Chap_p_vs_chap.df$geneID, ranks_chap_p.chap) 
colnames(ranks_chap_p.chap) <- c("GeneName","rank")
head(ranks_chap_p.chap)
#sort ranks in decreasing order
ranks_chap_p.chap <- ranks_chap_p.chap[order(as.numeric(ranks_chap_p.chap[,2]),decreasing = TRUE),]
write.table(ranks_chap_p.chap,file="UL_neuron_w20_yellow_arch_vs_correction.rnk",col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)
