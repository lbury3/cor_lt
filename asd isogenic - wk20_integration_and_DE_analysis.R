library(SCP)
library(hdf5r)
library(Seurat)

#loading the h5 file for Arch correction week 20 dataset. 
correction.w20.data <- Read10X_h5("correction_20wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")


correction.w20 <- CreateSeuratObject(counts = correction.w20.data, assay = 'RNA', project = 'Arch WT/WT.w20')
correction.w20 <- RunCellQC(srt = correction.w20)


correction.w20[["percent.mt"]] <- PercentageFeatureSet(correction.w20, pattern = "^MT-")

##### CUTOFF undecided. 
subset_correction.w20 <- subset(correction.w20, db.scDblFinder_class == "singlet"  & percent.mito < 20)
VlnPlot(subset_correction.w20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#### loading arch CTNNB1 WT/Q76* week 20 dataset
arch.w20.data <- Read10X_h5("arch_20wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")

arch.w20 <- CreateSeuratObject(counts = arch.w20.data, assay = 'RNA', project = 'Arch WT/Q76*.w20')
arch.w20[["percent.mt"]] <- PercentageFeatureSet(arch.w20, pattern = "^MT-")

arch.w20 <- RunCellQC(srt = arch.w20)
# cutoff using 15% leaves only 5500 cells.  decided to use 20%
subset_arch.w20 <- subset(arch.w20, db.scDblFinder_class == "singlet" & percent.mito < 20)
VlnPlot(subset_arch.w20, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# perform normalization with SCTransform v2. 
subset_correction.w20 <- SCTransform(subset_correction.w20, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)
#saveRDS(subset_correction.w20, "subset_correction.w20_SCT_v2.rds")

subset_arch.w20 <- SCTransform(subset_arch.w20, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)
#saveRDS(subset_arch.w20, "subset_arch.w20_SCT_v2.rds")

s_list = list( subset_correction.w20, subset_arch.w20)

features <- SelectIntegrationFeatures(object.list = s_list, nfeatures = 3000)
s_list <- PrepSCTIntegration(object.list = s_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = s_list, normalization.method = "SCT", 
                                  anchor.features = features)

organoid.combined.cca.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

organoid.combined.cca.sct <- RunPCA(organoid.combined.cca.sct, verbose = T)

# Plot PCA
PCAPlot(organoid.combined.cca.sct,
        split.by = "orig.ident")
organoid.combined.cca.sct <- RunUMAP(organoid.combined.cca.sct, reduction = "pca", dims = 1:30)
#organoid.combined.cca.sct <- RunTSNE(organoid.combined.cca.sct, reduction = "pca", dims = 1:30)

organoid.combined.cca.sct <- FindNeighbors(organoid.combined.cca.sct, dims = 1:30)
organoid.combined.cca.sct_1 <- FindClusters(organoid.combined.cca.sct, resolution = 1.5)
DimPlot(organoid.combined.cca.sct, reduction = 'umap', label = F, repel = T, split.by = "orig.ident", ncol = 3) 
DimPlot(organoid.combined.cca.sct_1, reduction = 'umap', label = T, repel = T) 


DefaultAssay(object = organoid.combined.cca.sct_1) <- "SCT"
organoid.combined.cca.sct_1<- PrepSCTFindMarkers(organoid.combined.cca.sct_1)
# marker genes for each cluster. 
markers<- FindAllMarkers(organoid.combined.cca.sct_1, assay = "SCT", 
                         only.pos = TRUE)


library(dplyr)
top200_ <- 
  markers %>%
  dplyr::group_by(cluster) %>%
  dplyr:: top_n(n = 500, wt = avg_log2FC)  #
write.csv(top200_,"top500_umap_31_markers_w20_prepSCT_findallmarkers_arch_archcorr_nobrain_3.6.23.csv", row.names = T)
# annotating the cell types for each cluster based on marker gene expression. 
Idents(organoid.combined.cca.sct_1) <- "seurat_clusters"


#no brain 
new_ident_1 <- setNames(c("oRG", #0  
                          "tRG",  #1  
                          "Pericytes",#2 
                          "EN_upper", #3
                          "IPC", #4
                          "EN_deep",  #5  
                          "EN_deep", #6 
                          "EN_upper", #7
                          "EN_upper", #8
                          "oRG/Astrocyte", #9  
                          "IPC", #10 , 
                          "U1", #11, 
                          "EN_upper", #12 
                          "oRG/Astrocyte", #13
                          "EN_upper",#14 
                          "Cycling",  #15  
                          "IPC", #16  
                          "Cycling", #17 
                          "EN_upper", #18
                          "Glyc",#19
                          "Choroid",#20
                          "U2_RPS",#21 
                          "EN_deep",#22
                          "U3_RPS",#23
                          "IPC",#24
                          "U4_RPS",#25
                          "oRG/Astrocyte",#26 
                          "U5", #27
                          "U6_TTR", #28
                          "Choroid", #29
                          "Pericytes" #30
),
levels(organoid.combined.cca.sct_1))
organoid.combined.cca.sct_1_rename <- RenameIdents(organoid.combined.cca.sct_1, new_ident_1)
DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)
organoid.combined.cca.sct_1_rename$CellType <- Idents(organoid.combined.cca.sct_1_rename)
saveRDS(organoid.combined.cca.sct_1_rename, "w20_arch_archcorr_nobrain.sct.cca_db_rm_mt_15_3.6.23.rds")
organoid.combined.cca.sct_1_rename<- readRDS("w20_arch_archcorr_nobrain.sct.cca_db_rm_mt_15_3.6.23.rds")


# Among lineage DEG generation for Arch WT/WT week 20 upper layer neuron. 
EN_upper <- subset(organoid.combined.cca.sct_1_rename, idents = "EN_upper")
DefaultAssay(EN_upper) <- "RNA"
Idents(EN_upper) <- "orig.ident"
EN_upper_correction <- subset(EN_upper, idents = "Arch WT/WT.w20")

EN_upper_correction$lineage <- "no_color"
EN_upper_correction$lineage[WhichCells(EN_upper_correction, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_upper_correction$lineage[WhichCells(EN_upper_correction, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_upper_correction$lineage[WhichCells(EN_upper_correction, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
      
    Idents(EN_upper_correction) <- "lineage"
    DefaultAssay(EN_upper_correction) <- "SCT"
    myTopHits.df <- FindMarkers(EN_upper_correction, 
                                logfc.threshold =0 ,  # equivalent to using 0
                                assay = "SCT",
                                min.pct = 0.1,
                                ident.1 = "yellow", ident.2 = "green",
                                test.use = "wilcox",
                                recorrect_umi = F
                                # max.cells.per.ident = 36
    )
saveRDS(myTopHits.df, "EN_upper_w20_correction_yellow_vs_green_recorrect_umi_F.rds")
    Chap_p_vs_chap.df <- myTopHits.df %>%
      as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_correction_yellow_vs_green_recorrect_umi_F.csv")

# Among lineage DEG generation for Arch WT/Q76* week 20 upper layer neuron. 
EN_upper_arch <- subset(EN_upper, idents = "Arch WT/Q76*.w20")

EN_upper_arch$lineage <- "no_color"
EN_upper_arch$lineage[WhichCells(EN_upper_arch, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_upper_arch$lineage[WhichCells(EN_upper_arch, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_upper_arch$lineage[WhichCells(EN_upper_arch, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
      
    Idents(EN_upper_arch) <- "lineage"
    DefaultAssay(EN_upper_arch) <- "SCT"
    #EN_deep_correction<- PrepSCTFindMarkers(EN_deep_correction)
    myTopHits.df <- FindMarkers(EN_upper_arch, 
                                logfc.threshold =0 ,  # equivalent to using 0
                                assay = "SCT",
                                min.pct = 0.1,
                                ident.1 = "yellow", ident.2 = "green",
                                test.use = "wilcox",
                                recorrect_umi = F
                                # max.cells.per.ident = 36
    )
saveRDS(myTopHits.df, "EN_upper_w20_arch_yellow_vs_green_recorrect_umi_F.rds")
myTopHits.df <- readRDS("EN_upper_w20_yellow_arch_vs_correction_recorrect_umi_T.rds")
    Chap_p_vs_chap.df <- myTopHits.df %>%
      as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_upper_w20_yellow_arch_vs_correction_recorrect_umi_T.csv")

###### creating rnk files for GSEA analysis. 
myTopHits.df <- readRDS("EN_upper_w20_correction_yellow_vs_red_recorrect_umi_F.rds")


Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
#create ranks file for GSEA 
ranks_chap_p.chap = Chap_p_vs_chap.df$avg_log2FC
ranks_chap_p.chap <- cbind(Chap_p_vs_chap.df$geneID, ranks_chap_p.chap) 
colnames(ranks_chap_p.chap) <- c("GeneName","rank")
head(ranks_chap_p.chap)
#sort ranks in decreasing order
ranks_chap_p.chap <- ranks_chap_p.chap[order(as.numeric(ranks_chap_p.chap[,2]),decreasing = TRUE),]
write.table(ranks_chap_p.chap,file="ranks_EN_upper_w20_correction_yellow_vs_red_recorrect_umi_F_sct_cca_wilcox_lfc.rnk",col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)





#### for within lineage CTNNB1 WT/Q76* vs WT/WT comparison for upper layer neurons. 


srt_EN_upper <- subset(organoid.combined.cca.sct_1_rename, idents = "EN_upper")
DefaultAssay(srt_EN_upper) <- "RNA"
#41 WT/WT, #25 WT/Q76
srt_EN_upper_Green_no_red_no <- subset(srt_EN_upper, ZsGreen == 0 & TdTomato ==0) 
saveRDS(srt_EN_upper_Green_no_red_no, "srt_w20_EN_upper_Green_no_red_no.rds")
#7 WT/WT, #17 WT/Q76* no enrichment
srt_EN_upper_Green_yes_red_no <- subset(srt_EN_upper, ZsGreen > 0 & TdTomato ==0) 
saveRDS(srt_EN_upper_Green_yes_red_no, "srt_w20_EN_upper_Green_yes_red_no.rds")

#944 WT/WT, 487 WT/Q76*
srt_EN_upper_Green_no_red_yes <- subset(srt_EN_upper,  ZsGreen == 0 &  TdTomato >0) 
saveRDS(srt_EN_upper_Green_no_red_yes, "srt_w20_EN_upper_Green_no_red_yes.rds")

#371 WT/WT, #987 WT/Q76*
srt_EN_upper_Green_yes_red_yes <- subset(srt_EN_upper,  ZsGreen > 0 &  TdTomato >0) 
saveRDS(srt_EN_upper_Green_yes_red_yes, "srt_w20_EN_upper_Green_yes_red_yes.rds")



Idents(srt_EN_upper) <- "orig.ident"
DefaultAssay(srt_EN_upper) <- "SCT"
srt_EN_upper<- PrepSCTFindMarkers(srt_EN_upper)
myTopHits.df <- FindMarkers(srt_EN_upper, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "Arch WT/Q76*.w20", ident.2 = "Arch WT/WT.w20",
                            test.use = "wilcox",
                            recorrect_umi = T
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_upper_w20_all_arch_vs_correction_recorrect_umi_T.rds")
library(tidyverse)
myTopHits.df <- readRDS("EN_deep_w20_red_arch_vs_correction_recorrect_umi_T.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_deep_w20_red_arch_vs_correction_recorrect_umi_T.csv")

#create ranks file for GSEA analysis. 
ranks_chap_p.chap = Chap_p_vs_chap.df$avg_log2FC
ranks_chap_p.chap <- cbind(Chap_p_vs_chap.df$geneID, ranks_chap_p.chap) 
colnames(ranks_chap_p.chap) <- c("GeneName","rank")
head(ranks_chap_p.chap)
#sort ranks in decreasing order
ranks_chap_p.chap <- ranks_chap_p.chap[order(as.numeric(ranks_chap_p.chap[,2]),decreasing = TRUE),]
write.table(ranks_chap_p.chap,file="ranks_EN_upper_all_arch_vs_correction_w20_recorrect_umi_T_sct_cca_wilcox_lfc.rnk",col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)


