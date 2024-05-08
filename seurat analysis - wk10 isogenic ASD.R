# Install the remotes package 
install.packages('remotes')
# Replace '4.1.0' with your desired version
#remotes::install_version(package = 'Seurat', version = package_version('4.1.0'))
#remotes::install_version("Matrix", version = "1.6-3")
library(Matrix)
library(Seurat)
library(irlba)
library(matrixStats)
library(tidyverse)
###
#install.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

arch_10wk_ali.data <- Read10X_h5("arch_10wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")
arch_10wk_ali <- CreateSeuratObject(counts = arch_10wk_ali.data, project = "arch_w10_WT/Q76*")

arch_10wk_ali[["percent.mt"]] <- PercentageFeatureSet(arch_10wk_ali, pattern = "^MT-")
arch_10wk_ali <- subset(arch_10wk_ali, subset = nFeature_RNA > 1250 & nFeature_RNA < 7500 & percent.mt < 15) 

arch_10wk_ali <- SCTransform(arch_10wk_ali, vst.flavor = "v2", vars.to.regress = "percent.mt")


correction_10wk_ali.data <- Read10X_h5("correction_10wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")
correction_10wk_ali <- CreateSeuratObject(counts = correction_10wk_ali.data, project = "arch_w10_WT/WT")
correction_10wk_ali[["percent.mt"]] <- PercentageFeatureSet(correction_10wk_ali, pattern = "^MT-")
correction_10wk_ali <- subset(correction_10wk_ali, subset = nFeature_RNA > 1250 & nFeature_RNA < 7500 & percent.mt < 15) 

correction_10wk_ali <- SCTransform(correction_10wk_ali, vst.flavor = "v2", vars.to.regress = "percent.mt")

merge.list <- list(correction_10wk_ali, arch_10wk_ali)
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



saveRDS(merge.combined, file = "arch_arch_corr_w10_ali_integrated_4.18.24.rds")
DefaultAssay(object = merge.combined) <- "SCT"
merge.combined<- PrepSCTFindMarkers(merge.combined)

markers <- FindAllMarkers(merge.combined, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(markers, "markers_w10_prepSCT_findallmarkers_arch_archcorr_4.18.24.rds")
write.csv(markers, file = "markers_w10_prepSCT_findallmarkers_arch_archcorr_4.18.24.csv")
library(Seurat)
merge.combined <- readRDS("w10_arch_archcorr_annotated_4.18.24.rds")

# annotating the cell types for each cluster based on marker gene expression. 
Idents(merge.combined) <- "seurat_clusters"

#no brain 
new_ident_1 <- setNames(c("RG/Ast", #0  
                          "RG/Ast",  #1  
                          "DL_neuron",#2 
                          "Other", #3 RG/Ast_PGK1
                          "DL_neuron", #4
                          "IPC",  #5  
                          "RG/Ast", #6 
                          "DL_neuron", #7
                          "Other", #8
                          "DL_neuron", #9  
                          "Cycling_NPC", #10 , 
                          "DL_neuron", #11, 
                          "DL_neuron", #12 
                          "Cycling_NPC", #13
                          "Other",#14 
                          "DL_neuron",  #15 
                          "Other", #16  
                          "IPC", #17
                          "Other", #18
                          "Cycling_NPC",#19
                          "Cycling_NPC",#20
                          "RG/Ast",#21
                          "Other",#22 Lhx1 and Lhx5 are expressed in CR cells 
                          "Other",#23
                          "Other",#24
                          "Cycling_NPC",#25
                          "Cycling_NPC",#26
                          "IPC",#27
                          "Cycling_NPC",#28
                          "Other",#29
                          "Other"#30

                          
),
levels(merge.combined))
organoid.combined.cca.sct_1_rename <- RenameIdents(merge.combined, new_ident_1)
DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)
organoid.combined.cca.sct_1_rename$CellType <- Idents(organoid.combined.cca.sct_1_rename)
saveRDS(organoid.combined.cca.sct_1_rename, "w10_arch_archcorr_annotated_4.22.24.rds")
library(Seurat)
organoid.combined.cca.sct_1_rename <- readRDS("w10_arch_archcorr_annotated_4.22.24.rds")
DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)
EN_deep <- subset(organoid.combined.cca.sct_1_rename, idents = "DL_neuron")
DefaultAssay(EN_deep) <- "RNA"
Idents(EN_deep) <- "orig.ident"
EN_deep_correction <- subset(EN_deep, idents = "arch_w10_WT/WT")

EN_deep_correction$lineage <- "no_color"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
# For lineage proportion for DL_neuron week 10 Arch WT/WT
#  total                   4067
# no color                 121
# yellow                   326
# red                      3611
# green                    9
##For week 10 Arch WT/WT Yellow vs red DE analysis.
Idents(EN_deep_correction) <- "lineage"
DefaultAssay(EN_deep_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "red",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_correction_yellow_vs_red_recorrect_umi_F.rds")
library(tidyverse)
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_correction_yellow_vs_red_recorrect_umi_F.csv")
# DL_neuron_w10_correction_yellow_vs_red_recorrect_umi_F.csv for Supplemeental Figure 17C, 
# For lineage proportion for DL_neuron week 10 Arch WT/WT

EN_deep_arch <- subset(EN_deep, idents = "arch_w10_WT/Q76*")
    
EN_deep_arch$lineage <- "no_color"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
  # For lineage proportion for DL_neuron week 10 Arch WT/Q76*

#  TOTAL                2917
# no color              86
# yellow                968
# red                   1832
# green                 31
##For week 10 Arch WT/Q76* Yellow vs red DE analysis.
Idents(EN_deep_arch) <- "lineage"
DefaultAssay(EN_deep_arch) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "red",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_arch_yellow_vs_red_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_arch_yellow_vs_red_recorrect_umi_F.csv")
# "DL_neuron_w10_correction_yellow_vs_red_recorrect_umi_F.csv" for Supplemeental Figure 17D, 

# for Supple  Figure 17B.
# DL_neuron correction. yellow vs no color.
Idents(EN_deep_correction) <- "lineage"
DefaultAssay(EN_deep_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_correction_yellow_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_correction_yellow_vs_nocolor_recorrect_umi_F.csv")
# DL_neuron correction. yellow vs green
Idents(EN_deep_correction) <- "lineage"
DefaultAssay(EN_deep_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "green",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_correction_yellow_vs_green_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_correction_yellow_vs_green_recorrect_umi_F.csv")
# DL_neuron correction. green vs nocolor
Idents(EN_deep_correction) <- "lineage"
DefaultAssay(EN_deep_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "green", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_correction_green_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_correction_green_vs_nocolor_recorrect_umi_F.csv")
# DL_neuron correction. red vs nocolor
Idents(EN_deep_correction) <- "lineage"
DefaultAssay(EN_deep_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "red", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_correction_red_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_correction_red_vs_nocolor_recorrect_umi_F.csv")
# DL_neuron correction. green vs red
Idents(EN_deep_correction) <- "lineage"
DefaultAssay(EN_deep_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "green", ident.2 = "red",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_correction_green_vs_red_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_correction_green_vs_red_recorrect_umi_F.csv")
# DL_neuron Arch.  yellow vs no color.
Idents(EN_deep_arch) <- "lineage"
DefaultAssay(EN_deep_arch) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_arch_yellow_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_arch_yellow_vs_nocolor_recorrect_umi_F.csv")
# DL_neuron arch . yellow vs green
Idents(EN_deep_arch) <- "lineage"
DefaultAssay(EN_deep_arch) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "green",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_arch_yellow_vs_green_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_arch_yellow_vs_green_recorrect_umi_F.csv")
# DL_neuron correction. green vs nocolor
Idents(EN_deep_arch) <- "lineage"
DefaultAssay(EN_deep_arch) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "green", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_arch_green_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_arch_green_vs_nocolor_recorrect_umi_F.csv")
# DL_neuron arch red vs nocolor
Idents(EN_deep_arch) <- "lineage"
DefaultAssay(EN_deep_arch) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "red", ident.2 = "no_color",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_arch_red_vs_nocolor_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_arch_red_vs_nocolor_recorrect_umi_F.csv")
# DL_neuron correction. green vs red
Idents(EN_deep_arch) <- "lineage"
DefaultAssay(EN_deep_arch) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_arch, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "green", ident.2 = "red",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_arch_green_vs_red_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_arch_green_vs_red_recorrect_umi_F.csv")



#rnk file generation
Chap_p_vs_chap <- readRDS("DL_neuron_w10_correction_yellow_vs_red_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- Chap_p_vs_chap %>%
  as_tibble(rownames = "geneID")

#create ranks file
ranks_chap_p.chap = Chap_p_vs_chap.df$avg_log2FC
ranks_chap_p.chap <- cbind(Chap_p_vs_chap.df$geneID, ranks_chap_p.chap) 
colnames(ranks_chap_p.chap) <- c("GeneName","rank")
head(ranks_chap_p.chap)
#sort ranks in decreasing order
ranks_chap_p.chap <- ranks_chap_p.chap[order(as.numeric(ranks_chap_p.chap[,2]),decreasing = TRUE),]
write.table(ranks_chap_p.chap,file="DL_neuron_w10_correction_yellow_vs_red.rnk",col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)

####yellow lineage arch vs correction
organoid.combined.cca.sct_1_rename <- readRDS("w10_arch_archcorr_annotated_4.18.24.rds")
DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)
DL_neuron <- subset(organoid.combined.cca.sct_1_rename, idents = "DL_neuron")
DefaultAssay(DL_neuron) <- "RNA"


#225 WT/WT, #706 WT/Q76*
srt_DL_neuron_Green_yes_red_yes <- subset(DL_neuron,  ZsGreen > 0 &  TdTomato >0) 
saveRDS(srt_DL_neuron_Green_yes_red_yes, "srt_DL_neuron_Green_yes_red_yes.rds")
Idents(srt_DL_neuron_Green_yes_red_yes) <- "orig.ident"
DefaultAssay(srt_DL_neuron_Green_yes_red_yes) <- "SCT"
srt_DL_neuron_Green_yes_red_yes<- PrepSCTFindMarkers(srt_DL_neuron_Green_yes_red_yes)
myTopHits.df <- FindMarkers(srt_DL_neuron_Green_yes_red_yes, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "arch_w10_WT/Q76*", ident.2 = "arch_w10_WT/WT",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_yellow_arch_vs_correction_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_yellow_arch_vs_correction_recorrect_umi_F.csv")

#red lineage week10
srt_DL_neuron_Green_no_red_yes <- subset(DL_neuron,  ZsGreen == 0 &  TdTomato >0) 
saveRDS(srt_DL_neuron_Green_no_red_yes, "srt_DL_neuron_Green_no_red_yes_w10.rds")

srt_DL_neuron_Green_no_red_yes<- readRDS("srt_DL_neuron_Green_no_red_yes_w10.rds")
##### 4.21.24 stop here. 
library(tidyverse)
Idents(srt_DL_neuron_Green_no_red_yes) <- "orig.ident"
DefaultAssay(srt_DL_neuron_Green_no_red_yes) <- "SCT"
srt_DL_neuron_Green_no_red_yes<- PrepSCTFindMarkers(srt_DL_neuron_Green_no_red_yes)
myTopHits.df <- FindMarkers(srt_DL_neuron_Green_no_red_yes, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "arch_w10_WT/Q76*", ident.2 = "arch_w10_WT/WT",
                            test.use = "wilcox",
                            recorrect_umi = F
)
saveRDS(myTopHits.df, "DL_neuron_w10_red_arch_vs_correction_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "DL_neuron_w10_red_arch_vs_correction_recorrect_umi_F.csv")

#for wk10 comparison of RG/Ast cells between CTNNB1 WT/WT and CTNNB1 WT/Q76* (Supplemental Table 6)

DefaultAssay(organoid.combined.cca.sct_1_rename) <- "RNA"
Idents(object = organoid.combined.cca.sct_1_rename) <- "seurat_clusters"
rg <- subset(x = organoid.combined.cca.sct_1_rename, idents = c(0,1,6,21))
Idents(object = rg) <- "orig.ident"
DefaultAssay(rg) <- "SCT"
zsgreen_ULneuron.markers <- FindMarkers(object = rg, logfc.threshold = 0, assay = "SCT", ident.1 = "arch_w10_WT/WT", ident.2 = "arch_w10_WT/Q76*", recorrect_umi = FALSE)
write.csv(zsgreen_ULneuron.markers, file = "arch_archcorrection_10wk_rg_comparison_markers.csv")


