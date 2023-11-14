########################################################################################
######  AMO1 Ire1 Xbp1 time depency data analysis by Ximo Pechuan i Jorge 5/09/2023
#######################################################################################

#################### 1. Set-up
# Libraries
library(EmbolcallRNAseq)
library(edgeR)
library(PCAtools)
library(Biobase)
library(tidyverse)
library(plyr)
library(sva)

# Working Directory
path.to.data = "~/Downloads/AshkenaziLab/Iratxe/NGS3940_XBP1vsIRE1_TimeSeries/"
# Gene Sets
species = "Homo sapiens"
ListSupport(c(wp2gene,msigDB,dbGSEA,keggdb),LoadGSEADB(species))
name = "ExploratoryAnalysis_Sept2023"
dir.new = paste0(path.to.data,name)
dir.create(dir.new)
setwd(dir.new)

################### 2. Loading and prepare data
# Read ESET 
es.all = readRDS("../data/NGS3940_ESET.Rdata")
pheno = pData(es.all) 
# DGEList 
y.all = DGEList(
  exprs(es.all),genes=fData(es.all))
y.all$samples = cbind(y.all$samples, pData(es.all))
y.all = calcNormFactors(y.all)

################################################ 3. Filter lowly expressed genes
cpms = 15/(median(y.all$samples$lib.size)/10^6)
keep.exprs = rowSums(edgeR::cpm(y.all) >= cpms) > 3
# We keep
sum(keep.exprs)/length(keep.exprs)
yf =  y.all[keep.exprs,, keep.lib.sizes=FALSE]
yf = calcNormFactors(yf,method="TMM")


############################ 4. Obtain expression matrices and annotations
#pheno
annotations = pheno[,c("shRNA","clone","DayTime")]
annotations$DayTime = annotations$DayTime %>% as.character()

# Aesthetics
ClonePalette = c("#ED574EFF","#F37062FF","#BED8EBFF","#2E5A87FF")
shRNAPalette = c("firebrick","#4674A0FF")
DayPalette = c("gray90","gray50","gray28","black")


# Top bar
names(shRNAPalette) = levels(factor(annotations$shRNA))
names(DayPalette) = levels(factor(annotations$DayTime))
names(ClonePalette) = levels(factor(annotations$clone, levels = 
                               c("clone1","clone3","clone1_2",
                                 "clone18")))

ann_colors = list(shRNA = shRNAPalette,
                  clone = ClonePalette,DayTime = DayPalette)
top_bar=HeatmapAnnotation(df=annotations,col=ann_colors)

# Expression matrices
ExprMat = cpm(yf,log=T) # After calling edger::calcNormFactors() is TMM normal
write.csv(ExprMat,"../data/TMMLogCPM.csv")
write.csv(pheno,"../data/Pheno.csv")

############################ 5. PCA
# How many PCAs in this dataset
scaled_mat = t(scale(t(ExprMat)))

p = pca(scaled_mat, metadata = pheno, removeVar = 0)
biplot(p, colby="Complete")

# Determine the rellevan PCs by talus plot
EmbolcallRNAseq::TalusPlot(p,10)
# Calculate Projection Score
ps.df = EmbolcallRNAseq::ThetaProjectionScore(DataMatrix = ExprMat,
                                              NPCs =  10, 
                                              nboot = 100,
                                              thetas =c(0.005,0.01,0.02,0.03,
                                                        0.05,0.06,0.07,0.075,0.08,0.085,0.09,
                                                        0.1,0.12, 0.15,0.17,
                                                        0.18,
                                                        0.2,0.23,0.25,0.28,
                                                        0.3,0.33,0.34,0.35,0.37,0.38,
                                                        0.4,0.45,0.48,
                                                        0.5,0.55,
                                                        0.65,0.75,0.8,1))
EmbolcallRNAseq::ProjectionScorePlot(ps.df)

# Filter the top fraction of the genes
nvar = floor(0.06*nrow(ExprMat))
nvar
topdat = EmbolcallRNAseq::NthMostVariableFeatures(ExprMat,nvar)
# scale
scaled_mat = t(scale(t(topdat))) # PCA
p = pca(scaled_mat, metadata = pheno, removeVar = 0)
biplot(p,colby = "Complete", showLoadings = T )

# Look at the top nvar genes Heatmap
scaled_mat = t(scale(t(topdat)))

# Annotate by Go Terms WIKI
g = GOHeatMap(scaled_mat,Gomet = "Wiki",nk_Go = 16,
              organismRef = "Homo sapiens",nHL = 5,Csplit=5, top_annotation = top_bar)
pdf("Wiki_TopVariable_genes.pdf", width = 15)
g
dev.off()

# Annotate by Go Terms Kegg
g = GOHeatMap(scaled_mat,Gomet = "Kegg",nk_Go = 16,
              organismRef = "Homo sapiens",nHL = 5,Csplit=5,top_annotation = top_bar)
pdf("Kegg_TopVariable_genes.pdf", width = 15)
g
dev.off()

############################ 6. PCA unsupervised analysis
# Quick exploration
p = pca(scaled_mat, metadata = pheno, removeVar = 0)
# How many PCAs?
horn = PCAtools::parallelPCA(scaled_mat)
elbow = findElbowPoint(p$variance)

pdf(paste0(name,"_ScreePlot.pdf"))
screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))
dev.off()

pdf(paste0(name,nvar,"_PC1_PC2.pdf"))
biplot(p, showLoadings = TRUE,
       labSize = 3, pointSize = 2, 
       sizeLoadingsNames = 3,
       colLoadingsNames = 'red4',
       fillBoxedLoadings = alpha("white", 2/4),
       colby="clone",
       colkey = ClonePalette,
       ntopLoadings=5,
       lab = NULL,
       alphaLoadingsArrow=0.8,
       colLegendTitle = 'Clone: ',
       # Ellipse
       ylim = c(-30,30),
       xlim = c(-50,50),
       ellipse =  T,
       hline = 0, vline = c(0),
       title = paste0("PCA of Top ",nvar," Variable Genes"),
       legendPosition = 'top', 
       legendLabSize = 10, 
       legendIconSize = 8.0)
dev.off()

pdf(paste0(name,nvar,"_PC1_PC3.pdf"))
biplot(p, 
       x= "PC1",
       y ="PC3",
       showLoadings = TRUE,
       labSize = 3, pointSize = 2, 
       sizeLoadingsNames = 3,
       colLoadingsNames = 'red4',
       fillBoxedLoadings = alpha("white", 2/4),
       colby="Complete",
       #colkey = c25,
       ntopLoadings=5,
       lab = NULL,
       alphaLoadingsArrow=0.8,
       colLegendTitle = 'Clone: ',
       # encircle config
       encircle = TRUE,
       encircleFill = TRUE,
       hline = 0, vline = c(0),
       title = paste0("PCA of Top ",nvar," Variable Genes"),
       legendPosition = 'top', 
       legendLabSize = 5, 
       legendIconSize = 5)
dev.off()

pdf(paste0(name,nvar,"RainBowCarnival_PC1_PC2.pdf"))
biplot(p, showLoadings = TRUE,
       labSize = 3, pointSize = 2, 
       sizeLoadingsNames = 3,
       colLoadingsNames = 'red4',
       fillBoxedLoadings = alpha("white", 2/4),
       colby="Complete",
       #colkey = c25,
       ntopLoadings=5,
       lab = NULL,
       alphaLoadingsArrow=0.8,
       colLegendTitle = 'Clone: ',
       # encircle config
       encircle = TRUE,
       encircleFill = TRUE,
       hline = 0, vline = c(0),
       title = paste0("PCA of Top ",nvar," Variable Genes"),
       legendPosition = 'top', 
       legendLabSize = 5, 
       legendIconSize = 5)
dev.off()


pdf("PCA_CorrelationPlot.pdf")
eigencorplot(p,
             metavars = c('DayTime','shRNA',"clone"))
dev.off()

# Sanity Check Rest,Ptprn,Gabrq
gensmbl="IRF4"
pdf(paste0(gensmbl,"CPM.pdf"))
print(stripchart(ExprMat[gensmbl,]~pheno$shRNA*pheno$DayTime,vertical=TRUE,
                 las=2,cex.axis=0.65,pch=16,cex=1.1,
                # col=c25,
                 method="jitter",ylab="TMM Log CPM",
                 main=gensmbl)
)
dev.off()


#################################### 6.  Pathway Analysis
library(progeny)
pathways = progeny(ExprMat,scale = T, organism = "Human", top = 100,
                   verbose = FALSE)

pdf("Progeny_Analysis.pdf", width = 12)
Heatmap(t(pathways),
        show_column_names = F,
        cluster_columns = F,
        show_column_dend = T,
        column_names_gp = gpar(fontsize = 10),
        column_split = pheno$shRNA,
        #column_order = sam_order,
        # column_order = col_order_id,
        border_gp = gpar(col = "gray", lty = 2),
        column_title ="PROGENy  Pathways",
        column_title_side = "top",
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        show_row_dend = TRUE,
        row_dend_side = "left",
        row_split = 5,
        #cluster_rows =cluster_within_group(t(scaled_mat), group$cluster),
        cluster_row_slices = TRUE,
        cluster_column_slices = T,
        row_title_rot = 0,
        #row_title = cluster_name,
        name="Scaled Pathway Score",
        top_annotation = top_bar
)
dev.off()


# Heatmap for S phase of cell cycle
library(Seurat)
s.genes = cc.genes$s.genes

Ifn_jeremy = ExprMat[(rownames(ExprMat) %in% s.genes),]


pdf("S_Cell_Cycle.pdf", width = 12)
Heatmap(t(scale(t(Ifn_jeremy))),
        show_column_names = F,
        cluster_columns = F,
        show_column_dend = T,
        column_names_gp = gpar(fontsize = 10),
        column_split =  annotations$shRNA,
        #column_order = sam_order,
        # column_order = col_order_id,
        border_gp = gpar(col = "gray", lty = 2),
        column_title ="S Phase Score",
        column_title_side = "top",
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        show_row_dend = TRUE,
        row_dend_side = "left",
        #row_split = 5,
        #cluster_rows =cluster_within_group(t(scaled_mat), group$cluster),
        cluster_row_slices = TRUE,
        cluster_column_slices = T,
        row_title_rot = 0,
        #row_title = cluster_name,
        name="Scaled TMM Log(cpm)",
        top_annotation = top_bar
)
dev.off()

# G-Score
scoredf = yf$samples
scoredf$pcSigJeremie = NULL
pcSig = gsScore(Ifn_jeremy)
scoredf$SPhaseScore = pcSig

# Plot
#my_comparisons = list(c("Tumor_Day1","Skin") ,
                  #    c("Tumor_Day3","Skin"))

pdf("SPhase_Boxplot_Score.pdf")
ggboxplot(scores, x = "Complete", y = "pcSigJeremie",
          color = "Complete", palette = c25 ,add="jitter")+
 # stat_compare_means(comparisons = my_comparisons,method ="t.test",
                    # method.args = list(alternative = "two.sided"),
                    # p.adjust.method = "fdr",var.equal=F,
                    # label="p.signif")+
  ylab("S Phase Gene Score")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+
  font("x.text", size = 11)+font("y.text", size = 11)
dev.off()

# Heatmap for G2M phase of cell cycle
g2m.genes = cc.genes$g2m.genes

Ifn_jeremy = ExprMat[(rownames(ExprMat) %in% g2m.genes),]


pdf("G2M_Cell_Cycle.pdf", width = 12)
Heatmap(t(scale(t(Ifn_jeremy))),
        show_column_names = F,
        cluster_columns = F,
        show_column_dend = T,
        column_names_gp = gpar(fontsize = 10),
        column_split =  pheno$shRNA,
        #column_order = sam_order,
        # column_order = col_order_id,
        border_gp = gpar(col = "gray", lty = 2),
        column_title ="G2M Phase Score",
        column_title_side = "top",
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        show_row_dend = TRUE,
        row_dend_side = "left",
        #row_split = 3,
        #cluster_rows =cluster_within_group(t(scaled_mat), group$cluster),
        cluster_row_slices = TRUE,
        cluster_column_slices = T,
        row_title_rot = 0,
        #row_title = cluster_name,
        name="Scaled TMM Log(cpm)",
        top_annotation = top_bar
)
dev.off()

# G-Score
pcSig = gsScore(Ifn_jeremy)
scoredf$G2MPhaseScore = pcSig
write.csv(scoredf,"CellCycleScores.csv")
# Plot
#my_comparisons = list(c("Tumor_Day1","Skin") ,
#    c("Tumor_Day3","Skin"))

pdf("G2M_Phase_Boxplot_InfScore.pdf")
ggboxplot(scores, x = "Complete", y = "pcSigJeremie",
          color = "Complete", palette = c25 ,add="jitter")+
  # stat_compare_means(comparisons = my_comparisons,method ="t.test",
  # method.args = list(alternative = "two.sided"),
  # p.adjust.method = "fdr",var.equal=F,
  # label="p.signif")+
  ylab("G2M Phase Gene Score")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+
  font("x.text", size = 11)+font("y.text", size = 11)
dev.off()



# HeatmapProgeny
progenymodel =  get("model_human_full", envir = .GlobalEnv)
subPro = progenymodel %>% filter(pathway == "p53") %>% top_n(weight,n = 100)
p53 = subPro$gene

Ifn_jeremy = ExprMat[(rownames(ExprMat) %in% p53),]


pdf("TP53.pdf", width = 12, height = 20)
Heatmap(t(scale(t(Ifn_jeremy))),
        show_column_names = F,
        cluster_columns = F,
        show_column_dend = T,
        column_names_gp = gpar(fontsize = 10),
        column_split =  pheno$shRNA,
        #column_order = sam_order,
        # column_order = col_order_id,
        border_gp = gpar(col = "gray", lty = 2),
        column_title ="Top 100 TP53 Pathway Genes",
        column_title_side = "top",
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        show_row_dend = TRUE,
        row_dend_side = "left",
        row_split = 6,
        #cluster_rows =cluster_within_group(t(scaled_mat), group$cluster),
        cluster_row_slices = TRUE,
        cluster_column_slices = T,
        row_title_rot = 0,
        #row_title = cluster_name,
        name="Scaled TMM Log(cpm)",
        top_annotation = top_bar
)
dev.off()

# G-Score
pcSig = gsScore(Ifn_jeremy)
yf$samples$pcSigJeremie = pcSig
scores = yf$samples
scores$Complete = factor(scores$Complete)

# Plot
#my_comparisons = list(c("Tumor_Day1","Skin") ,
#    c("Tumor_Day3","Skin"))

pdf("G2M_Phase_Boxplot_InfScore.pdf")
ggboxplot(scores, x = "Complete", y = "pcSigJeremie",
          color = "Complete", palette = c25 ,add="jitter")+
  # stat_compare_means(comparisons = my_comparisons,method ="t.test",
  # method.args = list(alternative = "two.sided"),
  # p.adjust.method = "fdr",var.equal=F,
  # label="p.signif")+
  ylab("G2M Phase Gene Score")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+
  font("x.text", size = 11)+font("y.text", size = 11)
dev.off()

