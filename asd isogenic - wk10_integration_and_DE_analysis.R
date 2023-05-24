library(SCP)
library(hdf5r)


correction.w10.data <- Read10X_h5("correction_10wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")


correction.w10 <- CreateSeuratObject(counts = correction.w10.data, assay = 'RNA', project = 'Arch WT/WT.w10')
correction.w10 <- RunCellQC(srt = correction.w10)
correction.w10[["percent.mt"]] <- PercentageFeatureSet(correction.w10, pattern = "^MT-")
VlnPlot(subset_correction.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
VlnPlot(subset_correction.w10, features = c("TBR1", "HOPX", "EOMES", "SATB2", "ZsGreen", "TdTomato"), ncol = 3, pt.size = 0.1)


subset_correction.w10 <- subset(correction.w10, db.scDblFinder_class == "singlet"  & percent.mito < 15)
#######
arch.w10.data <- Read10X_h5("arch_10wk_ali_dob20220603_raw_feature_bc_matrix_cellbender_filtered.h5")

arch.w10 <- CreateSeuratObject(counts = arch.w10.data, assay = 'RNA', project = 'Arch WT/Q76*.w10')

VlnPlot(subset_arch.w20, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

arch.w10 <- RunCellQC(srt = arch.w10)

subset_arch.w10 <- subset(arch.w10, db.scDblFinder_class == "singlet" & percent.mito < 15)
VlnPlot(subset_arch.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)




subset_correction.w10 <- SCTransform(subset_correction.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)
saveRDS(subset_correction.w10, "subset_correction.w10_SCT_v2.rds")

subset_arch.w10 <- SCTransform(subset_arch.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)
saveRDS(subset_arch.w10, "subset_arch.w10_SCT_v2.rds")

######### without braindataset. 
s_list = list( subset_correction.w10, subset_arch.w10)

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
saveRDS(organoid.combined.cca.sct, "w10_arch_archcorr_LTS_nobrain_3.5.23.rds")
organoid.combined.cca.sct <- readRDS("w10_arch_archcorr_LTS_nobrain_3.5.23.rds")

organoid.combined.cca.sct <- FindNeighbors(organoid.combined.cca.sct, dims = 1:30)
organoid.combined.cca.sct_1 <- FindClusters(organoid.combined.cca.sct, resolution = 1)
DimPlot(organoid.combined.cca.sct, reduction = 'umap', label = F, repel = T, split.by = "orig.ident", ncol = 3) 
DimPlot(organoid.combined.cca.sct_1, reduction = 'umap', label = T, repel = T) 


DefaultAssay(object = organoid.combined.cca.sct_1) <- "SCT"
organoid.combined.cca.sct_1<- PrepSCTFindMarkers(organoid.combined.cca.sct_1)

markers<- FindAllMarkers(organoid.combined.cca.sct_1, assay = "SCT", 
                         only.pos = TRUE)

saveRDS(markers, "markers_w10_prepSCT_findallmarkers_arch_archcorr_nobrain_3.5.23.rds")

library(dplyr)
top200_ <- 
  markers %>%
  dplyr::group_by(cluster) %>%
  dplyr:: top_n(n = 1000, wt = avg_log2FC)  #
write.csv(top200_,"top1000_markers_w10_prepSCT_findallmarkers_arch_archcorr_nobrain_3.5.23.csv", row.names = T)
# annotating the cell types for each cluster based on marker gene expression. 
Idents(organoid.combined.cca.sct_1) <- "seurat_clusters"

#no brain 
new_ident_1 <- setNames(c("vRG", #0  RPS
                          "EN_deep",  #1  
                          "oRG",#2 RPS
                          "oRG", #3
                          "Glyc", #4
                          "EN_deep",  #5  
                          "oRG", #6 also cycling
                          "EN_deep", #7
                          "IPC", #8
                          "EN_deep", #9  
                          "Microglia", #10 , 
                          "oRG_Cycling", #11, 
                          "EN_deep", #12 
                          "IPC", #13
                          "IPC",#14 
                          "U2_RPS",  #15 
                          "Cycling", #16  
                          "IN", #17
                          "EN_deep", #18
                          "Cycling",#19
                          "Cycling",#20
                          "U3_RPS",#21
                          "Cycling",#22
                          "Cycling",#23
                          "Pericytes",#24
                          "OPC",#25
                          "oRG_RPS"#26
                          
),
levels(organoid.combined.cca.sct_1))
organoid.combined.cca.sct_1_rename <- RenameIdents(organoid.combined.cca.sct_1, new_ident_1)
DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)
organoid.combined.cca.sct_1_rename$CellType <- Idents(organoid.combined.cca.sct_1_rename)
saveRDS(organoid.combined.cca.sct_1_rename, "w10_arch_archcorr_nobrain.sct.cca_db_rm_mt_15_3.5.23.rds")
organoid.combined.cca.sct_1_rename<- readRDS("w10_arch_archcorr_nobrain.sct.cca_db_rm_mt_15_3.5.23.rds")
##### Among lineage DE analysis for week 10 deep layer neurons in Arch WT/WT

EN_deep <- subset(organoid.combined.cca.sct_1_rename, idents = "EN_deep")
DefaultAssay(EN_deep) <- "RNA"
Idents(EN_deep) <- "orig.ident"
EN_deep_correction <- subset(EN_deep, idents = "Arch WT/WT.w10")

EN_deep_correction$lineage <- "no_color"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_deep_correction$lineage[WhichCells(EN_deep_correction, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
    
Idents(EN_deep_correction) <- "lineage"
DefaultAssay(EN_deep_correction) <- "SCT"
myTopHits.df <- FindMarkers(EN_deep_correction, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "yellow", ident.2 = "green",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_deep_w10_correction_yellow_vs_green_recorrect_umi_F.rds")
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
write.csv(Chap_p_vs_chap.df, "EN_deep_w10_correction_yellow_vs_green_recorrect_umi_F.csv")

##### Among lineage DE analysis for week 10 deep layer neurons in Arch WT/Q76*

EN_deep_arch <- subset(EN_deep, idents = "Arch WT/Q76*.w10")

EN_deep_arch$lineage <- "no_color"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato > 0 & ZsGreen >0)] <- "yellow"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato > 0 & ZsGreen ==0)] <- "red"
EN_deep_arch$lineage[WhichCells(EN_deep_arch, expression = TdTomato == 0 & ZsGreen >0)] <- "green"
      
    Idents(EN_deep_arch) <- "lineage"
    DefaultAssay(EN_deep_arch) <- "SCT"
    myTopHits.df <- FindMarkers(EN_deep_arch, 
                                logfc.threshold =0 ,  # equivalent to using 0
                                assay = "SCT",
                                min.pct = 0.1,
                                ident.1 = "green", ident.2 = "red",
                                test.use = "wilcox",
                                recorrect_umi = F
                                # max.cells.per.ident = 36
    )
    saveRDS(myTopHits.df, "EN_deep_w10_arch_green_vs_red_recorrect_umi_F.rds")

    myTopHits.df <- readRDS("EN_deep_w10_correction_yellow_vs_green_recorrect_umi_F.rds")


Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")
#create ranks file
ranks_chap_p.chap = Chap_p_vs_chap.df$avg_log2FC
ranks_chap_p.chap <- cbind(Chap_p_vs_chap.df$geneID, ranks_chap_p.chap) 
colnames(ranks_chap_p.chap) <- c("GeneName","rank")
head(ranks_chap_p.chap)
#sort ranks in decreasing order
ranks_chap_p.chap <- ranks_chap_p.chap[order(as.numeric(ranks_chap_p.chap[,2]),decreasing = TRUE),]
write.table(ranks_chap_p.chap,file="ranks_EN_deep_w10_correction_yellow_vs_green_recorrect_umi_F_sct_cca_wilcox_lfc.rnk",col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)







#### For within lineage DE analysis
DefaultAssay(srt_EN_deep) <- "RNA"
#166 WT/WT, #110 WT/Q76*
srt_EN_deep_Green_no_red_no <- subset(srt_EN_deep, ZsGreen == 0 & TdTomato ==0) 
srt_EN_deep_Green_no_red_no <- readRDS("srt_w10_EN_deep_Green_no_red_no.rds")
saveRDS(srt_EN_deep_Green_no_red_no, "srt_w10_EN_deep_Green_no_red_no.rds")

#12 WT/WT, 37 WT/Q76*
srt_EN_deep_Green_yes_red_no <- subset(srt_EN_deep, ZsGreen > 0 & TdTomato ==0) 
saveRDS(srt_EN_deep_Green_yes_red_no, "srt_w10_EN_deep_Green_yes_red_no.rds")

#3000 WT/WT, #1512 WT/Q76*
srt_EN_deep_Green_no_red_yes <- subset(srt_EN_deep,  ZsGreen == 0 &  TdTomato >0) 
saveRDS(srt_EN_deep_Green_no_red_yes, "srt_w10_EN_deep_Green_no_red_yes.rds")

#225 WT/WT, #706 WT/Q76*
srt_EN_deep_Green_yes_red_yes <- subset(srt_EN_deep,  ZsGreen > 0 &  TdTomato >0) 
saveRDS(srt_EN_deep_Green_yes_red_yes, "srt_w10_EN_deep_Green_yes_red_yes.rds")


Idents(srt_EN_deep) <- "orig.ident"
DefaultAssay(srt_EN_deep) <- "SCT"
srt_EN_deep<- PrepSCTFindMarkers(srt_EN_deep)
myTopHits.df <- FindMarkers(srt_EN_deep, 
                            logfc.threshold =0 ,  # equivalent to using 0
                            assay = "SCT",
                            min.pct = 0.1,
                            ident.1 = "Arch WT/Q76*.w10", ident.2 = "Arch WT/WT.w10",
                            test.use = "wilcox",
                            recorrect_umi = F
                            # max.cells.per.ident = 36
)
saveRDS(myTopHits.df, "EN_deep_w10_ALL_arch_vs_correction_recorrect_umi_F.rds")
library(tidyverse)
Chap_p_vs_chap.df <- myTopHits.df %>%
  as_tibble(rownames = "geneID")

#create ranks file
ranks_chap_p.chap = Chap_p_vs_chap.df$avg_log2FC
ranks_chap_p.chap <- cbind(Chap_p_vs_chap.df$geneID, ranks_chap_p.chap) 
colnames(ranks_chap_p.chap) <- c("GeneName","rank")
head(ranks_chap_p.chap)
#sort ranks in decreasing order
ranks_chap_p.chap <- ranks_chap_p.chap[order(as.numeric(ranks_chap_p.chap[,2]),decreasing = TRUE),]
write.table(ranks_chap_p.chap,file="ranks_EN_deep_correction_red_vs_nocolor_w10_recorrect_umi_F_sct_cca_wilcox_lfc.rnk",col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)



