####################################################################################################
######  Lung MERFISH  Progeny Pathway Activity Analysis for SPEC Ximo Pechuan i Jorge 8/03/2023
###################################################################################################

########################## 1. Set-up
# Libraries and Packages
library(Seurat)
library(progeny)
library(tidyverse)
library(ComplexHeatmap)
library(EnhancedVolcano)

# Wrapper function

ProgenyWrapper <- function(seu,organism="Mouse",identities){
  
  CellsClusters = data.frame(Cell = names(identities), 
                             CellType = as.character(identities),
                             stringsAsFactors = FALSE)
  ## We compute the Progeny activity scores and add them to our Seurat object
  ## as a new assay called Progeny. 
  seuLite = seu
  seuLite = progeny(seuLite, scale=FALSE, organism=organism, top=500, perm=1, 
                    return_assay = TRUE)
  
  ## We can now directly apply Seurat functions in our Progeny scores. 
  ## For instance, we scale the pathway activity scores. 
  seuLite = Seurat::ScaleData(seuLite, assay = "progeny") 
  
  ## We transform Progeny scores into a data frame to better handling the results
  progeny_scores_df =
    as.data.frame(t(GetAssayData(seuLite, slot = "scale.data", 
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 
  
  ## We match Progeny scores with the cell clusters.
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
  
  ## We summarize the Progeny scores by cellpopulation
  summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  ## We prepare the data for the plot
  summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
  paletteLength = 100
  
  myColor = colorRampPalette(c("steelblue1", "white", "tomato"))(paletteLength)
  
  progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(summarized_progeny_scores_df)/paletteLength, 
                        max(summarized_progeny_scores_df), 
                        length.out=floor(paletteLength/2)))
  
  progenymat = summarized_progeny_scores_df
  return(progenymat)
}

# Colore Palette
c25  = c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "coral3", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# Prefix for saving plots
name = "SPEXPathway"

setwd("~/Downloads/")
# Read the Seurat object
seu = readRDS("seurat_set.RDS")
DimPlot(seu)

# Read the output of the SPEX Communities Module
spex = read.csv("SPEXCommunityAnalysisOutput.csv",
                colClasses = 'character')

# Add it to the Seurat object as metadata
rownames(spex) = spex$cell_id
spex$SpexCluster = paste0("SpexCluster_",spex$phenograph_label)
seu = AddMetaData(seu,spex)
# Make it a factor for plotting
seu@meta.data$Cluster = factor(seu@meta.data$Cluster) 
seu@meta.data$SpexCluster = factor(seu@meta.data$SpexCluster) 


########################## 2. Fundamental exploratory analysis plots

# At this point, it is useful to have a visualizations umaps where you can change the label set. The umap is comming from
# the expression features and should be computed outside

DimPlot(seu, group.by = "SpexCluster",label=T, label.box = T, repel = T)+NoLegend()
DimPlot(seu, group.by = "Cluster",label=T, label.box = T, repel = T)+NoLegend()

#  same UMAPs but prettier
pdf(paste(name,"cluster_umap.pdf",sep=""),width = 12)
DimPlot(seu,group.by = "Cluster",label = T,cols = c25,
        pt.size = 1.2,repel = T,label.size = 4, label.box = T) + 
        scale_color_manual(labels = paste(levels(seu@meta.data$Cluster),table(seu@meta.data$Cluster),sep=": "),
                           values=c25)+
        ggtitle(paste0(name," Clusters"))
dev.off()

pdf(paste(name,"Spexcluster_umap.pdf",sep=""),width = 12)
DimPlot(seu,group.by = "SpexCluster",label = T,cols = c25,
        pt.size = 1.2,repel = T,label.size = 4, label.box = T) + 
        scale_color_manual(labels = paste(levels(seu@meta.data$SpexCluster),table(seu@meta.data$SpexCluster),sep=": "),
                           values=c25)+
        ggtitle(paste0(name," SpexCluster"))
dev.off()

# Slide properties, it is important to have visualizations for the slide coordinates that show gene expression, spex clusters 
# and clusters themselves
df = seu@meta.data

# SPEX clusters
pdf("SpexC_on_slide.pdf", width = 10)
ggplot(df) + aes(x, y, col = SpexCluster) + geom_point(size=0.1)+ theme_bw()+
        theme(axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.text.y = element_blank())+ 
        scale_color_manual(values=c(c25,c25))+
        guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

# Clusters
pdf("Cluster_on_slide.pdf", width = 8)
ggplot(df) + aes(x, y, col = Cluster) + geom_point(size=0.1)+ theme_bw()+
        theme(axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),axis.text.y = element_blank())+ 
        scale_color_manual(values=c25)+
        guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

# SPEX cluster composition  heatmap

mat = table(seu@meta.data$SpexCluster, seu@meta.data$Cluster) %>% scale() %>% t()

pdf(paste0(name,"_SpexComposition.pdf"),width = 12)
Heatmap(mat)
dev.off()



########################## 3. Progeny Pathway analyses
# Choose the level, here spex clusters
Idents(seu) = seu@meta.data$SpexCluster
identities = Idents(seu) 
progenymat = ProgenyWrapper(seu, identities = identities,organism = "Human")

pdf(paste0(name,"_Progeny_pathway_by_Cluster.pdf"),width = 30)
Heatmap(t(progenymat),
        cluster_columns = F,
        cluster_rows = T,
        row_split = 5, 
        column_names_gp = gpar(fontsize = 10),
        border_gp = gpar(col = "gray", lty = 2),
        row_dend_gp = gpar(col = "gray"), column_dend_gp = gpar(col = "gray"))
dev.off()

# Heatmap of the pathway analysis of each cell type in the context of each SPEX community

Idents(seu) = paste0(seu@meta.data$Cluster,"_",seu@meta.data$SpexCluster)
identities = Idents(seu) 
progenymat = ProgenyWrapper(seu, identities = identities,organism = "Human")


celltypes = lapply(rownames(progenymat), 
                   function(x) as.character(unlist(str_split(x,pattern= "_Spex"))[1])) %>%
        as.character()
labs = lapply(rownames(progenymat), 
                   function(x) as.character(unlist(str_split(x,pattern= "_Spex"))[2])) %>%
        as.character()      

pdf(paste0(name,"_Progeny_pathway_by_Cluster.pdf"),width = 30)
Heatmap(t(progenymat),
        cluster_columns = F,
        cluster_rows = T,
        row_split = 5, 
        column_names_gp = gpar(fontsize = 10),
        column_split = celltypes,
        column_labels = labs,
        border_gp = gpar(col = "gray", lty = 2),
        row_dend_gp = gpar(col = "gray"), column_dend_gp = gpar(col = "gray"))
dev.off()

########################## 4. Marker Differential Expression Analysis

# This code is the same for the CLQ but the labels are different, avoidant/attractive
name1 = "Fibroblasts_SpexCluster_2.0"
name2 = "Fibroblasts_SpexCluster_1.0"
# Markers from the cluster will work if  Idents(seu) = "...."  is set to the right fields
Markers = FindMarkers(seu,ident.1 = name1,
                      ident.2 = name2,
                      logfc.threshold = 0,min.pct = 0)

EnhancedVolcano(Markers,
                lab = rownames(Markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.2,
                pCutoff=0.05,
                title = paste0(name1,"vs",name2),
                subtitle = "",
                titleLabSize = 10)
