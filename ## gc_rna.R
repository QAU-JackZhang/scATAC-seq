###-------------- Analysis for gc_rna -------------------###
###---------------- by jack_zhang_qau --------------------###

#-------Seurat Object ------------------------
rm(list=ls())
gc_rna <- readRDS("~/ATAC-data/Result-scRNA/RNA_rds/gc_rna.rds")

library(Seurat)
set.seed(111)

# PCA + UMAP  
gc_rna <- NormalizeData(object = gc_rna)
gc_rna <- FindVariableFeatures(object = gc_rna, nfeatures = 2000)
gc_rna <- ScaleData(object = gc_rna)
gc_rna <- RunPCA(object = gc_rna, npcs = 50, verbose = FALSE)
ElbowPlot(gc_rna,ndims = 50)

gc_rna <- FindNeighbors(gc_rna, dims = 1:20)
gc_rna <- FindClusters(gc_rna, resolution = 1, algorithm = 3)
gc_rna <- RunUMAP(object = gc_rna, dims = 1:20, n.neighbors = 25)

DimPlot(object = gc_rna, label = TRUE, pt.size = 1)
VlnPlot(gc_rna, features = c("nFeature_RNA", "nCount_RNA"),pt.size = 0)

# Tree 
gc_rna <- BuildClusterTree(gc_rna)
Tool(object = gc_rna, slot = 'BuildClusterTree')
PlotClusterTree(gc_rna)

# Filter_low_quality_cluster 
gc_rna <- subset(gc_rna, idents = c(6,8), invert = TRUE)

# Re-UMAP 
gc_rna <- NormalizeData(object = gc_rna)
gc_rna <- FindVariableFeatures(object = gc_rna, nfeatures = 2000)
gc_rna <- ScaleData(object = gc_rna)
gc_rna <- RunPCA(object = gc_rna, npcs = 50, verbose = FALSE)
ElbowPlot(gc_rna,ndims = 50)

gc_rna <- FindNeighbors(gc_rna, dims = 1:20)
gc_rna <- FindClusters(gc_rna, resolution = 1, algorithm = 3)
gc_rna <- RunUMAP(object = gc_rna, dims = 1:20, n.neighbors = 25)

DimPlot(object = gc_rna, label = TRUE, pt.size = 1)
DimPlot(object = gc_rna, label = TRUE, pt.size = 1, 
        group.by = "state", cols = c("purple","deepskyblue","gold"))
DimPlot(object = gc_rna, label = TRUE, pt.size = 1, 
        split.by = "state", cols = c("purple","deepskyblue","gold"))


### 查看每个cluster 数量
table(gc_rna$seurat_clusters)
barplot(table(gc_rna$seurat_clusters))

FindMarkers(gc_rna, ident.1 = 8, min.pct = 0.25)

###  查看不同细胞类型百分比 ## ggpubr包内ggpie函数
ggpie(as.data.frame(prop.table(table(Idents(gc_rna)))),'Freq',fill = 'Var1')
a= as.matrix(prop.table(table(Idents(gc_rna), gc_rna$state), margin = 2))

# Stage Markers
VlnPlot(gc_rna, features = c("Pou5f1","Utf1", "Sycp3","Syce3"), pt.size = 0,log = T, ncol = 1)
VlnPlot(gc_rna, features = c("Stra8","Dmc1","Ugt8a","Meioc"), pt.size = 0,log = T, ncol = 1)


VlnPlot(gc_rna, features = c("Pou5f1","Utf1","Sycp3","Syce3","Stra8"),
        pt.size = 0,log = T, ncol = 1,
        group.by = "state", cols = c("purple","deepskyblue","gold"))

VlnPlot(gc_rna, features = c("Dmc1","Ugt8a","Meioc","Cdk2","Tcfl5"),
        pt.size = 0,log = T, ncol = 1,
        group.by = "state", cols = c("purple","deepskyblue","gold"))

#  命名
new.cluster.ids <- c(
  'Mitotic',
  'Transition_phase',
  'Transition_phase',
  'Leptotene',
  'Transition_phase',
  'preLeptotene',
  'preLeptotene',
  'Zygotene',
  'Transition_phase'
  )

names(x = new.cluster.ids) <- levels(x = gc_rna)

# 重命名
gc_rna <- RenameIdents(gc_rna, new.cluster.ids)
gc_rna$celltype <- Idents(gc_rna)

DimPlot(object = gc_rna, label = TRUE, pt.size = 1) 
barplot(table(gc_rna$celltype))

# 鉴定 marker
markers <- FindAllMarkers(gc_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = "~/ATAC-data/Result-scRNA/gc_rna/markers.csv")

# average gene expression
aver = AverageExpression(gc_rna)
write.csv(aver, file = "~/ATAC-data/Result-scRNA/gc_rna/aver_genes.csv")

#------ Building trajectories with Monocle 3 --------
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(111)

gc_rna = readRDS("~/ATAC-data/Result-scRNA/RNA_rds/gc_name.rds")
Germline.cds <- as.cell_data_set(gc_rna)

Germline.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(gc_rna[["RNA"]])
Germline.cds <- estimate_size_factors(Germline.cds)
Germline.cds <- preprocess_cds(Germline.cds)

Germline.cds <- cluster_cells(cds = Germline.cds, reduction_method = "UMAP")
Germline.cds <- learn_graph(Germline.cds, use_partition = TRUE)
Germline.cds <- order_cells(Germline.cds, reduction_method = "UMAP")

plot_cells(cds = Germline.cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE, 
           cell_size = 1, 
           graph_label_size = 6,
           trajectory_graph_color = "black",
           trajectory_graph_segment_size = 1) & scale_color_gradient(low = '#FF8C00', high = 'blue')


gc_rna <- AddMetaData(
  object = gc_rna,
  metadata = Germline.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Germline")

######  plot gene trajectory 
plot_cells(
  Germline.cds, genes=c("Stra8","Rec8","Tcfl5"),
  show_trajectory_graph=T,
  cell_size = 0.5,
  graph_label_size = 6,
  label_leaves=T)

genes <- c("Stra8","Tcfl5")
lineage_cds <- Germline.cds[rowData(Germline.cds)$gene_short_name %in% genes]

plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5,
                         cell_size = 1.5) & scale_color_gradient(low = '#FF8C00', high = 'blue')
# 3d 
cds_3d <- reduce_dimension(Germline.cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d)

plot_cells_3d(cds_3d, color_cells_by="pseudotime")

# sub_branch 
cds_sub <- choose_graph_segments(cds)
cds_sub <- estimate_size_factors(cds_sub)
cds_sub = preprocess_cds(cds = cds_sub, method = "PCA")
cds_sub = reduce_dimension(cds = cds_sub, reduction_method = "UMAP")
cds_sub = cluster_cells(cds = cds_sub, reduction_method = "UMAP")
cds_sub = learn_graph(cds = cds_sub, use_partition = T)
cds_sub = order_cells(cds = cds_sub, reduction_method = "UMAP")

plot_cells(
cds = cds_sub,
color_cells_by = "pseudotime",
show_trajectory_graph = TRUE, cell_size = 1, graph_label_size = 5)


# 画基因分化热图
ciliated_cds_pr_test_res <- graph_test(Germline.cds, 
                                       neighbor_graph="principal_graph", 
                                       cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(Germline.cds[pr_deg_ids,], resolution=0.6)

cell_group_df <- tibble::tibble(cell=row.names(colData(Germline.cds)),
                                cell_group=colData(Germline.cds)$celltype)

agg_mat <- aggregate_gene_expression(Germline.cds, gene_module_df, cell_group_df)

row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))

pheatmap::pheatmap(agg_mat,scale="column", 
                   clustering_method="ward.D2",
                   colorRampPalette(c('#FF8C00','blue'))(255),
                   border_color =F)

# 挑选分化关键基因
a =as.data.frame(gene_module_df)
a = dplyr::filter(a, module == '53'|module =='18'|module =='44'|module =='1'|module =='51'|module =='25'|module =='8'|module =='31'|module =='52'|module =='9'|module =='82'|module =='64'|module =='73'|module =='48'|module =='83'|module =='62'|module =='22'|module =='77')

# 功能注释
library(clusterProfiler)
library(org.Mm.eg.db)
ego<-enrichGO(a$id,
              OrgDb = org.Mm.eg.db,
              keyType = 'SYMBOL',
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.05)

barplot(ego,showCategory = 10)
write.csv(ego,"go_diff_core_gene.csv")


#-------Monocle2 -----------------------
detach("package:monocle3", unload = TRUE)
library(monocle)

a <- as.CellDataSet(gc_rna)

marker_genes <- row.names(subset(fData(a), gene_short_name %in% c("Stra8", "Rec8")))

# 归一化
a <- estimateSizeFactors(a)
a <- estimateDispersions(a)
a <- detectGenes(a, min_expr = 3 )
print(head(fData(a)))
expressed_genes <- row.names(subset(fData(a), num_cells_expressed >= 10))
print(head(pData(a)))


#  The ordering workflow 
diff_test_res <- differentialGeneTest(a[expressed_genes,], fullModelFormulaStr = "~num_genes_expressed")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

a <- setOrderingFilter(a, ordering_genes)
plot_ordering_genes(a)

a <- reduceDimension(a, max_components = 2, method = 'DDRTree')

a <- orderCells(a)

plot_cell_trajectory(a, color_by = "seurat_clusters")
plot_cell_trajectory(a, color_by = "celltype")
plot_cell_trajectory(a, color_by = "state")
plot_cell_trajectory(a, color_by = "Pseudotime")


###  “状态”只是单片树的术语，表示这段树。
###  下面的函数便于识别包含时间为零的大多数细胞的状态。
###  然后我们可以把它传递给orderCells
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

a1 <- orderCells(a, root_state = GM_state(a))
plot_cell_trajectory(a1, color_by = "Pseudotime")

### 按照不同分类画图
plot_cell_trajectory(a, color_by = "state") + facet_wrap(~state, nrow = 1)
plot_cell_trajectory(a, color_by = "celltype") + facet_wrap(~celltype, nrow = 1)

### 拟时间分化热图
diff_test_res <- differentialGeneTest(a, fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))

p = plot_pseudotime_heatmap(a[ordering_genes,],
                          num_clusters = 3,
                          cores = 1,return_heatmap=T,
                          show_rownames = T)

## 提取DEGs 

clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.table(clustering, 
            file = "~/ATAC-data/Result-scRNA/gc_rna_PCA/Filter-7-8/monocle2/mono_heatmap.txt", 
            sep = "\t",
            row.names = T, 
            quote = F)

####------Multi-Factorial Differential Expression Analysis

to_be_tested <-
  row.names(subset(fData(a),
                   gene_short_name %in% c("Pou5f1","Stra8", "Sycp3","Dmc1","Ugt8a")))

cds_subset <- a[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~celltype + state",
                                      reducedModelFormulaStr = "~state")
diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_jitter(cds_subset,
                  grouping = "state", color_by = "celltype", plot_trend = TRUE) +
  facet_wrap( ~ feature_label, scales= "free_y")




FeaturePlot(
  object = gc_rna,
  features = c("Tcfl5","Zbtb14","Zbtb7a","Hinfp"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.5,
  cols = c('#fee0d2','#fc9272','#de2d26')
)




rm(list=ls())
sessionInfo()