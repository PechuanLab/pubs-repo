#!/usr/bin/env Rscript
##################################################################################
######  NGS2180 Brain by Ximo Pechuan i Jorge 21/11/2022
#################################################################################

####################################################################### 1. Set-up 
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

setwd("~/Downloads/")
# read data
cell_bender_mat = Read_CellBender_h5_Mat("~/Downloads/LIB5428352_HITS5447650_filtered.h5",)
cell_bender_mat = cell_bender_mat[,!(colSums(cell_bender_mat) == 0)]
cell_bender_mat = as(cell_bender_mat,"dgCMatrix")
seu = CreateSeuratObject(counts = cell_bender_mat)
seu[["percent.mt"]] =  PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] =  PercentageFeatureSet(seu, pattern = "^RP[SL]")
# Add Metadata
seu@meta.data$PatientID = "DEP3"
seu@meta.data$HITS = "LIB5428352_HITS5447650"
seu@meta.data$Tissue = "rSMG"
seu@meta.data$Diagnosis = "DEP"
seu@meta.data$Operator = "JC"
seu@meta.data$SampleName = "S00-0176"
seu@meta.data$Triad = "T219"
seu@meta.data$DetailedID = "DEP_rSMG_3"
name = "LIB5428352_HITS5447650"

########################### 2. QC metrics
nFeathLow = 500
nCountHigh = 20000
riboT = 1
mitoT = 1
# Manually decide on the qc threshold nFeature_RNA , nCount_RNA
VlnPlot(seu,"nFeature_RNA", pt.size = 0, log = T) & 
  geom_hline(yintercept = nFeathLow)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = F) & 
  geom_hline(yintercept = nCountHigh)
VlnPlot(seu,"percent.ribo", pt.size = 0, log = F) & 
  geom_hline(yintercept = riboT)
VlnPlot(seu,"percent.mt", pt.size = 0, log = F) +NoLegend()+
  geom_hline(yintercept = mitoT)

seu  = subset(seu, subset = nCount_RNA < nCountHigh &  percent.mt < mitoT & nFeature_RNA > nFeathLow & percent.ribo < riboT)
sce = as.SingleCellExperiment(seu)
plotHighestExprs(sce, exprs_values = "counts")

########################### 2. Final Clean-up and clustering
# Quick SCT
seu = NormalizeData(seu, scale.factor = median(seu@meta.data$nCount_RNA))
median(seu@meta.data$nCount_RNA)
seu  = CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
seu = SCTransform(seu, verbose = T, conserve.memory = T, vars.to.regress = c("S.Score","G2M.Score","percent.ribo","percent.mt"))
seu = RunPCA(seu,npcs = 100)
ndims = 20
TalPlot(seu,ndims)
seu = FindNeighbors(seu, reduction = "pca", dims = 1:ndims)
seu = RunUMAP(seu, reduction = "pca", dims = 1:ndims)
DimPlot(seu)

# Resolutions
resolutions = c(0.1,0.2,0.3,0.4,0.5,0.6)
for (i in 1:length(resolutions)) {
  res = resolutions[i]
  seu = FindClusters(seu, verbose = FALSE,resolution=res)
}
# Assess clustering
seu = BuildClusterTree(seu)
clustree(seu,prefix="SCT_snn_res.", node_colour_aggr = "median")

Idents(seu) =  seu@meta.data$SCT_snn_res.0.1
DimPlot(seu,label=T)
DimPlot(seu)
DimPlot(seu, group.by = "Phase")
FeaturePlot(seu, features  = "nFeature_RNA")
FeaturePlot(seu,"percent.mt")
FeaturePlot(seu,"percent.ribo")
FeaturePlot(seu, features  = "PLP1")
FeaturePlot(seu, features  = "DRD2")
FeaturePlot(seu, features  = "DRD1")
FeaturePlot(seu,"ROBO1")
FeaturePlot(seu,"CADPS2")
FeaturePlot(seu,"TMEM106B")

########################################### 3. Grubman from SingleR
# load
# Libraries
library(DataSetDB)
library(gp.sa.core)
library(gp.sa.diff)
library(gp.sa)
library(tidyverse)
library(gp.sa.core)
library(gp.sa.solo)
library(scRNAseq)
library(ggpubr)
# Get dataset from datasetDB
DS = "DS000010394"
se = getDatasetAsSE(DS, experiment = 1)
se = scuttle::logNormCounts(se)
plotUMAP(se, colour_by = "clusterName")

# Change to Symbol
library(scater)
rownames(se) = uniquifyFeatureNames(rowData(se)$ID, rowData(se)$symbol)
head(rownames(se))

sce = as.SingleCellExperiment(seu,assay = "RNA")

# Run Single R
predicted = SingleR(test = sce , ref = se ,
                    labels = colData(se)$clusterName,
                    assay.type.ref = "logcounts"
                    
)

# Add metadata to Seurat
GrubmanLabels = predicted$labels
names(GrubmanLabels) = rownames(predicted)
seu = AddMetaData(seu, GrubmanLabels, col.name = "GrubmanLabels")
table(seu@meta.data$GrubmanLabels)

# For each relevant phenotype
pdf(paste0(name,"IGrubmanLabels_umap.pdf"),width = 12, height = 12)
DimPlot(seu,group.by = "GrubmanLabels",label = T,
        pt.size = 1,repel = T,label.size = 8, label.box = T) + theme_pubr(base_size = 10) +
  ggtitle(paste0(name,":",unique(seu@meta.data$DetailedID),":",unique(seu@meta.data$Operator)))+NoLegend()
dev.off()

# Qucik FindMarkers
# Markers from the cluster
Idents(seu) =  seu@meta.data$SCT_snn_res.0.1
pdf(paste0(name,"umap.pdf"),width = 12, height = 12)
DimPlot(seu,label=T)
dev.off()

Markers = FindAllMarkers(seu,only.pos = T,logfc.threshold = 0.5)
Markers = Markers %>% filter(p_val_adj<0.05)
write.csv(Markers,paste0(name,"markers.csv"))
topMarkers =  Markers %>% group_by(cluster) %>% top_n(avg_log2FC,n=10)

pdf(paste(name,"topMarkers.pdf",sep=""))
DoHeatmap(seu,unique(topMarkers$gene), group.colors  = c25,size = 2.5, draw.lines = T,raster=T) +
  theme(axis.text.y = element_text(size = 3)) +  scale_fill_gradientn(colors =  c("steelblue1", "white", "tomato"))
dev.off()


saveRDS(seu,"Seurat_LIB5428352_HITS5447650.rds")
