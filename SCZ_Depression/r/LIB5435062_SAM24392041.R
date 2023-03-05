#!/usr/bin/env Rscript
##################################################################################
######  NGS3422 Brd9HSC by Ximo Pechuan i Jorge 14/02/2023 
#################################################################################

######## 1. Set-up 
# Libraries
library(scRNAseq)
library(Seurat)
library(SingleR)
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)
library(scater)
library(EmbolcallRNAseq)
source("/Users/pechuanj/Documents/Workstation/stamp/STAMP_ngs_analyses/src/singlecell/STAMP_scRNAseq.R")

setwd("~/Downloads/NGS3422_Brd9/")
# read data
cell_bender_mat = Read_CellBender_h5_Mat("~/Downloads/LIB5435062_SAM24392041_filtered.h5")
cell_bender_mat = cell_bender_mat[,!(colSums(cell_bender_mat) == 0)]
cell_bender_mat = as(cell_bender_mat,"dgCMatrix")
seu = CreateSeuratObject(counts = cell_bender_mat)

# Add ribo and mito
seu[["percent.mt"]] =  PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.ribo"]] =  PercentageFeatureSet(seu, pattern = "^Rp[sl]")

# Add Metadata
seu@meta.data$Genotype = "Brd9N216Y"
seu@meta.data$HITS = "LIB5435062_SAM24392041"
seu@meta.data$Tissue = "BM"
seu@meta.data$SAMID = "SAM24392041"
seu@meta.data$PaperID = "Brd9N216Y_Bm_2"
name = "LIB5435062_SAM24392041"

########################### 2. QC metrics
# Limits
nFeathLow = 2500
nCountHigh = 30000
riboT = 30
mitoT = 5
# Manually decide on the qc threshold nFeature_RNA , nCount_RNA
VlnPlot(seu,"nFeature_RNA", pt.size = 0, log = T) & 
  geom_hline(yintercept = nFeathLow)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = F) & 
  geom_hline(yintercept = nCountHigh)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = T) & 
  geom_hline(yintercept = nFeathLow)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = T) & 
  geom_hline(yintercept = nFeathLow)
VlnPlot(seu,"percent.ribo", pt.size = 0, log = F) & 
  geom_hline(yintercept = riboT)
VlnPlot(seu,"percent.mt", pt.size = 0, log = F) +NoLegend()+
  geom_hline(yintercept = mitoT)

seu  = subset(seu, subset = nCount_RNA < nCountHigh &  percent.mt < mitoT & nFeature_RNA > nFeathLow & percent.ribo < riboT)
sce = as.SingleCellExperiment(seu)
plotHighestExprs(sce, exprs_values = "counts",n = 100)

########################### 2. Final Clean-up and clustering
# Quick SCT
seu = NormalizeData(seu)
seu  = CellCycleScoring(seu, s.features = firstup(tolower(cc.genes$s.genes)), g2m.features = firstup(tolower(cc.genes$g2m.genes)), set.ident = F)
seu = SCTransform(seu, verbose = T, conserve.memory = T, vars.to.regress = c("S.Score","G2M.Score","percent.ribo","percent.mt"))
seu = RunPCA(seu,npcs = 100)
ndims = 40
TalPlot(seu,ndims)
seu = FindNeighbors(seu, reduction = "pca", dims = 1:ndims)
seu = RunUMAP(seu, reduction = "pca", dims = 1:ndims)
DimPlot(seu, group.by = "Phase")
FeaturePlot(seu,"percent.mt")
FeaturePlot(seu,"percent.ribo")
FeaturePlot(seu,"nCount_RNA")
# Resolutions
resolutions = c(0.1,0.2,0.3,0.4,0.5,0.6)
for (i in 1:length(resolutions)) {
  res = resolutions[i]
  seu = FindClusters(seu, verbose = FALSE,resolution=res)
}
# Assess clustering
seu = BuildClusterTree(seu)
clustree(seu,prefix="SCT_snn_res.", node_colour_aggr = "median")
Idents(seu) = seu@meta.data$SCT_snn_res.0.6
DimPlot(seu)

########################################### 3. Immgen from SingleR
# load
immgen = ImmGenData(ensembl=FALSE)
sce = as.SingleCellExperiment(seu,assay = "RNA")

# Run Single R
predicted = SingleR(test = sce , ref = immgen ,
                    labels = colData(immgen)$label.fine,
                    assay.type.ref = "logcounts"
                    
)

# Add metadata to Seurat
immgenLabels = predicted$labels
names(immgenLabels) = rownames(predicted)
seu = AddMetaData(seu, immgenLabels, col.name = "immgen")

Markers = FindAllMarkers(seu,only.pos = T)
Markers = Markers %>% filter(p_val_adj<0.05)
Markers %>% filter(cluster == 0) %>% top_n(avg_log2FC,n=10)
DimPlot(seu, label = T, label.box = T)

# For each relevant phenotype
pdf("ImmGenSingleR_umap.pdf",width = 12, height = 12)
DimPlot(seu,group.by = "immgen",label = T,
        pt.size = 1,repel = T,label.size = 2, label.box = T) + theme_pubr(base_size = 10) +
  ggtitle("ImmGen Clusters")+NoLegend()
dev.off()
table(seu@meta.data$immgen)
saveRDS(seu,"Seurat_LIB5435062_SAM24392041.rds")
