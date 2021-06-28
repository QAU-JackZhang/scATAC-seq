############################################################
###--------------- Analysis for GSM4643732 --------------###
###---------------- by jack_zhang_qau -------------------###
library(Seurat)
library(ggplot2)
library(ggpubr)
library(stringr)
set.seed(111)

setwd("~/ATAC-data/Result-scRNA/GSM484632/")


gr11 = Read10X("~/ATAC-data/GSM4643732/GR11/")
gr11 = CreateSeuratObject(gr11, project = "gr11")
gr11[["percent.mt"]] <- PercentageFeatureSet(gr11, pattern = "^mt-")

gr12 = Read10X("~/ATAC-data/GSM4643732/GR12/")
gr12 = CreateSeuratObject(gr12, project = "gr12")
gr12[["percent.mt"]] <- PercentageFeatureSet(gr12, pattern = "^mt-")




aggr = merge(x= gr11, y= gr12)

VlnPlot(aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)


fivenum(aggr)

aggr <- subset(x = aggr, 
                   subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000 &
                   percent.mt < 10)


### PCA + UMAP
aggr <- NormalizeData(object = aggr)
aggr <- FindVariableFeatures(object = aggr, nfeatures = 2000)
aggr <- ScaleData(object = aggr)
aggr <- RunPCA(object = aggr, npcs = 50, verbose = FALSE)
ElbowPlot(aggr,ndims = 50)

aggr <- FindNeighbors(aggr, dims = 1:30)
aggr <- FindClusters(aggr, resolution = 0.6, algorithm = 3)
aggr <- RunUMAP(object = aggr, dims = 1:30, n.neighbors = 35)

DimPlot(object = aggr, label = TRUE, pt.size = 0.5)

DimPlot(object = aggr, label = TRUE, pt.size = 0.5, 
        group.by = "orig.ident")


VlnPlot(aggr, features = c("Ddx4","Dazl","Dppa3"),pt.size = 0, ncol = 1, log = T)

# 查看每个cluster 数量
barplot(table(aggr$seurat_clusters))

# 
FeaturePlot(object = aggr,
            features = c( 'Ddx4','Dazl','Dppa3'),
            pt.size = 0.1,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

# Filter/Extract
gc <- subset(aggr, idents = c(11))

# UMAP
gc <- NormalizeData(object = gc)
gc <- FindVariableFeatures(object = gc, nfeatures = 2000)
gc <- ScaleData(object = gc)
gc <- RunPCA(object = gc, npcs = 50, verbose = FALSE)
ElbowPlot(gc,ndims = 50)

gc <- FindNeighbors(gc, dims = 1:15)
gc <- FindClusters(gc, resolution = 0.6, algorithm = 3)
gc <- RunUMAP(object = gc, dims = 1:15, n.neighbors = 35)

DimPlot(object = gc, label = TRUE, pt.size = 2)
DimPlot(object = gc, label = TRUE, pt.size = 2, 
        group.by = "orig.ident")


# UMAP
gc <- NormalizeData(object = gc)
gc <- FindVariableFeatures(object = gc, nfeatures = 2000)
gc <- ScaleData(object = gc)
gc <- RunPCA(object = gc, npcs = 50, verbose = FALSE)
ElbowPlot(gc,ndims = 50)

gc <- FindNeighbors(gc, dims = 1:15)
gc <- FindClusters(gc, resolution = 0.6, algorithm = 3)
gc <- RunUMAP(object = gc, dims = 1:15, n.neighbors = 35)

DimPlot(object = gc, label = TRUE, pt.size = 2)
DimPlot(object = gc, label = TRUE, pt.size = 2, 
        group.by = "orig.ident")


## 计算基因表达
cluster.averages <- AverageExpression(gc)

write.csv(cluster.averages,"mean_Exp.csv")




##### 清除变量
rm(list=ls())

sessionInfo()