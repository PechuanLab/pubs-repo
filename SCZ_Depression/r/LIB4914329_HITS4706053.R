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
source("/gstore/scratch/u/pechuanj/Ishizuka/scRNAseq_AnalysisFunctions.R")

setwd("/gstore/scratch/u/pechuanj/NGS2180_Brain/BrainJasmine/LIB4914329_HITS4706053/outs")
# read data
cell_bender_mat = Read_CellBender_h5_Mat("HITS4706053_filtered.h5",)
cell_bender_mat = cell_bender_mat[,!(colSums(cell_bender_mat) == 0)]
cell_bender_mat = as(cell_bender_mat,"dgCMatrix")
seu = CreateSeuratObject(counts = cell_bender_mat)
seu[["percent.mt"]] =  PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] =  PercentageFeatureSet(seu, pattern = "^RP[SL]")
# Add Metadata
seu@meta.data$PatientID = "SCZ3"
seu@meta.data$HITS = "LIB4914329_HITS4706053"
seu@meta.data$Tissue = "rCH"
seu@meta.data$Diagnosis = "Schizophrenia"
seu@meta.data$Operator = "JC"
seu@meta.data$SampleName = "S10-0257"
seu@meta.data$Triad = "T219"
name = "LIB4914329_HITS4706053"

########################### 2. QC metrics
# Limits
nFeathLow = 500
nCountHigh = 40000
riboT = 5
mitoT = 2
# Manually decide on the qc threshold nFeature_RNA , nCount_RNA
VlnPlot(seu,"nFeature_RNA", pt.size = 0, log = T) & 
  geom_hline(yintercept = nFeathLow)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = F) & 
  geom_hline(yintercept = nCountHigh)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = T) & 
  geom_hline(yintercept = nFeathLow)
VlnPlot(seu,"percent.ribo", pt.size = 0, log = F) & 
  geom_hline(yintercept = riboT)
VlnPlot(seu,"percent.mt", pt.size = 0, log = F) +NoLegend()+
  geom_hline(yintercept = mitoT)

seu  = subset(seu, subset = nCount_RNA < nCountHigh &  percent.mt < mitoT & nFeature_RNA > nFeathLow & percent.ribo < riboT)
sce = as.SingleCellExperiment(seu)
plotHighestExprs(sce, exprs_values = "counts")

########################### 2. Final Clean-up and clustering
# Quick SCT
seu = SCTransform(seu, verbose = T, conserve.memory = F)
seu = RunPCA(seu,npcs = 100)
ndims = 27
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

Idents(seu) =  seu@meta.data$SCT_snn_res.0.2
DimPlot(seu,label=T)

seu  = CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)
DimPlot(seu, group.by = "Phase")
FeaturePlot(seu, features  = "nFeature_RNA")
FeaturePlot(seu, features  = "nCount_RNA")
FeaturePlot(seu, features  = "PLP1")
FeaturePlot(seu, features  = c("DRD1","DRD2"))
FeaturePlot(seu, features  = "percent.mt")
FeaturePlot(seu, features  = "CADPS2")

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


# For each relevant phenotype
pdf(paste0(name,"IGrubmanLabels_umap.pdf"),width = 12, height = 12)
DimPlot(seu,group.by = "GrubmanLabels",label = T,
        pt.size = 1,repel = T,label.size = 2, label.box = T) + theme_pubr(base_size = 10) +
  ggtitle("GrubmanLabels Clusters")+NoLegend()
dev.off()

table(seu@meta.data$GrubmanLabels)

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
DoHeatmap(seu,unique(topMarkers$gene), group.colors  = c25,size = 5, draw.lines = T,raster=T) +
  theme(axis.text.y = element_text(size = 5)) +  scale_fill_gradientn(colors =  c("steelblue1", "white", "tomato"))
dev.off()

FeaturePlot(seu,"CADPS2")
VlnPlot(seu,"CADPS2")
saveRDS(seu,"Seurat_HITS4706053.rds")
