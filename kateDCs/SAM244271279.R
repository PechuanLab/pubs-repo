# Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(SingleR)

setwd("~/Downloads/NGS4936_KateDCs/data/")

########################################### 1. Read and prepare data
# Antibody info
abs = read.csv("../NGS4936_BarcodeFeature.csv")
htos = abs$name[1:7]
cites =  abs$name[8:16]
# Read data
data = Read10X(data.dir = "SAM24427129_CellRanger/")
data_hto = data$`Antibody Capture`[htos,]
data_adt = data$`Antibody Capture`[cites,]
# Change ADTs Rownames
rownames(data_hto) = c("dLN1_1","dLN1_2","dLN2_2","dLN3_2","dLN3_1","dLN2_1")
rownames(data_adt) = abs$id[8:16]

seu = CreateSeuratObject(counts = data$`Gene Expression`)
seu[["HTO"]] = CreateAssayObject(counts =data_hto)
seu[["CiteSeq"]] = CreateAssayObject(counts =data_adt)

# Add experiment info
seu@meta.data$Experiment = "dLN"

########################### 2. HTO Demux
# Demultiplex
seu = NormalizeData(seu, assay = "HTO", normalization.method = "CLR",verbose = T)
seu = HTODemux(seu, assay = "HTO", positive.quantile = 0.99)
table(seu$HTO_classification.global)

Idents(seu) = "HTO_classification.global"
VlnPlot(seu, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
HTOHeatmap(seu, assay = "HTO", ncells = 5000)

# Keep only singlets
seu = subset(seu, idents = "Singlet")

########################### 2. QC metrics
# QC
seu[["percent.mt"]] =  PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.ribo"]] =  PercentageFeatureSet(seu, pattern = "^Rp[sl]")

nFeathLow = 500
nCountHigh = 15000
riboT = 40
mitoT = 5
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

########################### 3. Quick SCT
# Quick SCT
seu = NormalizeData(seu, scale.factor = median(seu@meta.data$nCount_RNA))
median(seu@meta.data$nCount_RNA)
s.genes = cc.genes$s.genes %>% tolower() %>% firstup()
g2m.genes = cc.genes$g2m.genes %>% tolower() %>% firstup()
seu  = CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
seu = SCTransform(seu, verbose = T, conserve.memory = T, vars.to.regress = c("S.Score","G2M.Score","percent.ribo","percent.mt"))
seu = RunPCA(seu,npcs = 100)
ndims = 25

# umaps
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

Idents(seu) =  seu@meta.data$SCT_snn_res.0.6
DimPlot(seu,label=T)
DimPlot(seu)
DimPlot(seu, group.by = "Phase")
FeaturePlot(seu, features  = "nFeature_RNA")
FeaturePlot(seu,"percent.mt")
FeaturePlot(seu,"percent.ribo")
FeaturePlot(seu, features  = "Foxp3")
FeaturePlot(seu, features  = "Cd8a")
FeaturePlot(seu, features  = "Cd4")
FeaturePlot(seu,"Foxp3")
FeaturePlot(seu,"Krt19")
FeaturePlot(seu,"S100a8")

# Map to STAMP atlas
stamp = readRDS("~/Documents/Workstation/stamp/NGS2751_STAMP_Atlas/data/STAMP_atlas_August2022.Rdata")
DimPlot(stamp, group.by = "ClusterDetail")


se = as.SingleCellExperiment(stamp,assay = "RNA")
sce = as.SingleCellExperiment(seu,assay = "RNA")

# Run Single R
predicted = SingleR(test = sce , ref = se ,
                    labels = colData(se)$ClusterDetail,
                    assay.type.ref = "logcounts"
                    
)



# Add metadata to Seurat
GrubmanLabels = predicted$labels
names(GrubmanLabels) = rownames(predicted)
seu = AddMetaData(seu, GrubmanLabels, col.name = "STAMPLabels")
table(seu@meta.data$STAMPLabels)

# For each relevant phenotype
pdf(paste0(name,"_umap.pdf"),width = 12, height = 12)
DimPlot(seu,group.by = "STAMPLabels",label = T,
        pt.size = 1,repel = T,label.size = 3, label.box = T)+NoLegend()
dev.off()

# load
immgen = ImmGenData(ensembl=FALSE)

# Run Single R
predicted = SingleR(test = sce , ref = immgen ,
                    labels = colData(immgen)$label.fine,
                    assay.type.ref = "logcounts"
                    
)


# Add metadata to Seurat
GrubmanLabels = predicted$labels
names(GrubmanLabels) = rownames(predicted)
seu = AddMetaData(seu, GrubmanLabels, col.name = "ImmGenLabels")
table(seu@meta.data$ImmGenLabels)

# For each relevant phenotype
pdf(paste0(name,"_umap.pdf"),width = 12, height = 12)
DimPlot(seu,group.by = "ImmGenLabels",label = T,
        pt.size = 1,repel = T,label.size = 3, label.box = T) + theme_pubr(base_size = 10) +
  ggtitle("ImmGen Labels")+NoLegend()
dev.off()


# 


seu = NormalizeData(seu, assay = "CiteSeq", normalization.method = "CLR",verbose = T)
abs
FeaturePlot(seu,"aPDL1")

# Save
library(Seurat)
library(SeuratData)
library(SeuratDisk)

SaveH5Seurat(seu, filename = "SAM24427129.h5Seurat")
Convert("SAM24427129.h5Seurat", dest = "h5ad")
saveRDS(seu,"SAM24427129.rds")

