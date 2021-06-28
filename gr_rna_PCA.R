############################################################
#----------------- Analysis for Aggr_rna --------------###
#---------------- by jack_zhang_qau -------------------###
library(Seurat)
library(ggplot2)
library(stringr)
set.seed(111)

counts <- Read10X_h5(filename = "~/ATAC-data/scRNA-seq/cellranger/RNA_aggr/outs/filtered_feature_bc_matrix.h5")
gr_rna <- CreateSeuratObject(counts , project = 'aggr', min.cells = 10)

gr_rna$state = as.factor(ifelse(str_detect(names(Idents(gr_rna)), "[1-2]$") == "TRUE", ifelse(str_detect(names(Idents(gr_rna)),"2$") == "TRUE",'GR12','GR11'),'GR13'))

gr_rna[["percent.mt"]] <- PercentageFeatureSet(gr_rna, pattern = "^mt-")

VlnPlot(gr_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = 'blue',pt.size = 0)

FeatureScatter(object = gr_rna, 
               feature1 = "nFeature_RNA", 
               feature2 = "percent.mt",
               cols = 'deepskyblue',pt.size = 0.1)+ 
               geom_vline(xintercept=c(1000,6000), lty=4,col="black",lwd=0.6)+ 
               geom_hline(yintercept = c(0,8),lty=4,col="black",lwd=0.6)

FeatureScatter(object = gr_rna, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               cols = 'deepskyblue',pt.size = 0.1)+ 
               geom_vline(xintercept=c(1500),lty=4,col="black",lwd=0.6)+ 
               geom_hline(yintercept = c(1000,6000), lty=4,col="black",lwd=0.6)

### ------ Filter--------#####
gr_rna <- subset(x = gr_rna, 
                 subset = nFeature_RNA > 1000 & 
                   nFeature_RNA < 6000 &
                   percent.mt < 8)


### PCA + UMAP
gr_rna <- NormalizeData(object = gr_rna)
gr_rna <- FindVariableFeatures(object = gr_rna, nfeatures = 2000)
gr_rna <- ScaleData(object = gr_rna)

gr_rna <- RunPCA(object = gr_rna, npcs = 50, verbose = FALSE)
## 根据各成分解释的方差百分比对主要成分进行排序
ElbowPlot(gr_rna,ndims = 50)

gr_rna <- FindNeighbors(gr_rna, dims = 1:15)
gr_rna <- FindClusters(gr_rna, resolution = 0.6, algorithm = 3)
gr_rna <- RunUMAP(object = gr_rna, dims = 1:15, n.neighbors = 35)

DimPlot(object = gr_rna, label = TRUE, pt.size = 0.5)
DimPlot(object = gr_rna, label = TRUE, pt.size = 0.5, 
        group.by = "state", cols = c("#D94DFF","#E68AB8","#fdd51e"))


#########  剔除某组
gr_rna <- subset(gr_rna, idents = 12, invert = TRUE)

#### 查看每个cluster 数量
table(gr_rna$seurat_clusters)
barplot(table(gr_rna$seurat_clusters))


### 系统发育分析
gr_rna <- BuildClusterTree(gr_rna)
Tool(object = gr_rna, slot = 'BuildClusterTree')
PlotClusterTree(gr_rna)

## Markers
VlnPlot(gr_rna, features = c("Ddx4","Dazl","Dppa3","Wnt4","Wnt6"),pt.size = 0, ncol = 1, log = T)

VlnPlot(gr_rna, features = c("Upk3b","Krt19","Col1a2","Bgn","Pecam1"),pt.size = 0, ncol = 1, log = T)

VlnPlot(gr_rna, features = c("Kdr","Alas2","Cx3cr1","Car2","Cd52"),pt.size = 0, ncol = 1, log = T)


### Germ cells
FeaturePlot(object = gr_rna,
            features = c( 'Dppa3','Ddx4','Dazl'),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

### preGranulosa cells
FeaturePlot(object = gr_rna,
            features = c( 'Wnt4','Wnt6' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

### Epithelial cells (上皮细胞)
FeaturePlot(object = gr_rna,
            features = c('Upk3b','Krt19' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

### Interstitial cells (间质细胞)
FeaturePlot(object = gr_rna,
            features = c( 'Col1a2','Bgn' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

### Endothelial cells (内皮细胞)
FeaturePlot(object = gr_rna,
            features = c( 'Pecam1','Kdr' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

### Immune cells
FeaturePlot(object = gr_rna,
            features = c( 'Cd52','Car2' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

### Erythroid cells (红细胞)
FeaturePlot(object = gr_rna,
            features = c( 'Alas2','Cx3cr1' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c( "#7e5ab8","#d79cd2", "#fdd51e"))

###------Cluster marker_gene
cluster9.markers <- FindMarkers(gr_rna, ident.1 = 9, min.pct = 0.25)

Markers <- FindAllMarkers(gr_rna, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Markers, file = "~/ATAC-data/Result-scRNA/Aggr/markers.csv")

### gr_rna_aggr 
new.cluster.ids <- c(
  'Pre_Granulosa',
  'Pre_Granulosa',
  'Epithelial',
  'Pre_Granulosa',
  'Germ_cells',
  'Interstitial',
  'Interstitial',
  'Interstitial',
  'Epithelial',
  'Pre_Granulosa',
  'Epithelial',
  'Endothelial',
  'Blood_related',
  'Immune_related'
  )


names(x = new.cluster.ids) <- levels(x = gr_rna)

### 重命名
gr_rna <- RenameIdents(gr_rna, new.cluster.ids)

gr_rna$celltype <- Idents(gr_rna)

DimPlot(object = gr_rna, label = TRUE, pt.size = 0.5) 


# 提取特定cluster，继续后续分析。
ident_df <- data.frame(cell=names(Idents(gr_rna)), cluster=Idents(gr_rna))
gc_rna <- subset(gr_rna, cells=as.vector(ident_df[ident_df$cluster=="Germ_cells",1]))

saveRDS(gc_rna,file = "~/ATAC-data/Result-scRNA/gc_rna_PCA/gc_rna_pca.rds")



##########################################################
###------ Building trajectories with Monocle 3 --------###
##########################################################

library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)

Germline.cds <- as.cell_data_set(gr_rna)

Germline.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(gc_rna[["RNA"]])
Germline.cds <- estimate_size_factors(Germline.cds)

Germline.cds <- cluster_cells(cds = Germline.cds, reduction_method = "UMAP")
Germline.cds <- learn_graph(Germline.cds, use_partition = TRUE)
Germline.cds <- order_cells(Germline.cds, reduction_method = "UMAP")

plot_cells(
  cds = Germline.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, 
  cell_size = 1, 
  graph_label_size = 5)

gc_rna <- AddMetaData(
  object = gc_rna,
  metadata = Germline.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Germline")

######  plot gene trajectory ########
plot_cells(
  Germline.cds, genes=c("Stra8"),
  show_trajectory_graph=T,
  cell_size = 1,
  graph_label_size = 6,
  label_leaves=T)

AFD_genes <- c("Stra8")
AFD_lineage_cds <- Germline.cds[rowData(Germline.cds)$gene_short_name %in% AFD_genes]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5)


##### 清除变量
rm(list=ls())

sessionInfo()