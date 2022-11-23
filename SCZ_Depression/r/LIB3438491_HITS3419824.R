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

setwd("/gstore/scratch/u/pechuanj/NGS2180_Brain/BrainJasmine/LIB3438491_HITS3419824/outs")
# read data
cell_bender_mat = Read_CellBender_h5_Mat("HITS3419824.h5",)
cell_bender_mat = cell_bender_mat[,!(colSums(cell_bender_mat) == 0)]
cell_bender_mat = as(cell_bender_mat,"dgCMatrix")
seu = CreateSeuratObject(counts = cell_bender_mat)
seu[["percent.mt"]] =  PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] =  PercentageFeatureSet(seu, pattern = "^RP[SL]")
# Add Metadata
seu@meta.data$PatientID = "SCZ1"
seu@meta.data$HITS = "LIB3438491_HITS3419824"
seu@meta.data$Tissue = "rSMG"
seu@meta.data$Diagnosis = "Schizophrenia"
seu@meta.data$Operator = "HI"
seu@meta.data$SampleName = "S10-0290"
seu@meta.data$Triad = "T216"
name = "LIB3438491_HITS3419824"
########################### 2. QC metrics

nFeathLow = 500
nCountHigh = 10000
riboT = 2
mitoT = 2.5
# Manually decide on the qc threshold nFeature_RNA , nCount_RNA
VlnPlot(seu,"nFeature_RNA", pt.size = 0, log = T) & 
  geom_hline(yintercept = nFeathLow)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = F) & 
  geom_hline(yintercept = nCountHigh)
VlnPlot(seu,"nCount_RNA", pt.size = 0, log = F) & 
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
ndims = 21
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

seu  = CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
DimPlot(seu)
DimPlot(seu, group.by = "Phase")
FeaturePlot(seu, features  = "nFeature_RNA")
FeaturePlot(seu, features  = "PLP1")

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

FeaturePlot(seu,"ROBO1")
saveRDS(seu,"Seurat_HITS3419824.rds")
