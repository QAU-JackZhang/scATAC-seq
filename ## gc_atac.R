###-------------- Analysis for gc_atac -------------------###
###---------------- by jack_zhang_qau --------------------###
###------------------- Seurat Object ------------------------ #####
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(stringr)
set.seed(111)

gc_atac <- readRDS("~/ATAC-data/Results-Signac/atac_rds/gc_atac.rds")
DefaultAssay(gc_atac) <- 'peaks'

VlnPlot(
  object = gc_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'nCount_peaks', 'nFeature_peaks'), pt.size = 0)


# Normalization and linear dimensional reduction
gc_atac <- RunTFIDF(gc_atac)
gc_atac <- FindTopFeatures(gc_atac, min.cutoff = 'q5')
gc_atac <- RunSVD(gc_atac)

DepthCor(gc_atac, n = 50) 
ElbowPlot(gc_atac,ndims = 50,reduction = 'lsi')

gc_atac <- FindNeighbors(gc_atac, reduction = 'lsi', dims = 2:30)
gc_atac <- FindClusters(gc_atac, resolution = 1, algorithm = 3)
gc_atac <- RunUMAP(gc_atac, reduction = 'lsi', dims = 2:30, n.neighbors = 25)

# Filter low_quality
gc_atac <- subset(gc_atac, idents = 6, invert = TRUE)


DimPlot(gc_atac, label = TRUE, pt.size = 1) 
DimPlot(gc_atac, label = TRUE, pt.size = 1, 
        group.by = "state",cols = c("#E6005C","#22C32E","#007FFF"))

#### 查看每个cluster 数量
table(gc_atac$seurat_clusters)
barplot(table(gc_atac$seurat_clusters))


# Integration with scRNA-seq 
gc_rna <- readRDS("~/ATAC-data/Result-scRNA/RNA_rds/gc_name.rds")
DefaultAssay(gc_atac) <- 'RNA'

gc_rna <- FindVariableFeatures(gc_rna,nfeatures = 5000)

transfer.anchors <- FindTransferAnchors(
  reference = gc_rna,  
  query = gc_atac, 
  reduction = 'cca')

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = gc_rna$celltype,
  weight.reduction = gc_atac[['lsi']])

gc_atac <- AddMetaData(object = gc_atac, metadata = predicted.labels)

plot1 <- DimPlot(
  object = gc_rna,pt.size = 1,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE)  + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = gc_atac,pt.size = 1,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE)  + ggtitle('scATAC-seq')

plot1 + plot2


#  命名
new.cluster.ids <- c(
  'preLeptotene',
  'Zygotene',
  'Meiotic_interphase',
  'Mitotic',
  'preLeptotene',
  'Meiotic_interphase',
  'Meiotic_interphase',
  'Leptotene',
  'Meiotic_interphase'
  )

names(x = new.cluster.ids) <- levels(x = gc_atac)

# 重命名
gc_atac <- RenameIdents(gc_atac, new.cluster.ids)
gc_atac$celltype <- Idents(gc_atac)

DimPlot(object = gc_atac, label = TRUE, pt.size = 1) 
barplot(table(gc_atac$celltype))

# 展示不同群之间的峰 
DefaultAssay(gc_atac) <- 'peaks'

## Dazl
CoveragePlot(gc_atac, 
             region = "chr17-50279393-50296000", 
             features = "Dazl",
             annotation = TRUE,
             peaks = TRUE,
             tile = TRUE,
             links = TRUE)
VlnPlot(gc_rna, features = c("Dazl"),pt.size = 0)

## Stra8
CoveragePlot(gc_atac, 
             region = "chr6-34917000-34940000", 
             features = "Stra8",
             annotation = TRUE,
             peaks = TRUE,
             tile = TRUE,
             links = TRUE)
VlnPlot(gc_rna, features = c("Stra8"),pt.size = 0)

## Tcfl5
CoveragePlot(gc_atac, 
             region = "chr2-180628000-180645000", 
             features = "Tcfl5",
             annotation = TRUE,
             peaks = F,
             tile = F,
             links = F)
VlnPlot(gc_rna, features = c("Stra8"),pt.size = 0)





## Rec8
CoveragePlot(gc_atac, 
             region = "chr14-55617000-55624000", 
             features = "Rec8",
             annotation = TRUE,
             peaks = TRUE,
             tile = TRUE,
             links = TRUE)
VlnPlot(gc_rna, features = c("Rec8"),pt.size = 0)


## Ugt8a
CoveragePlot(gc_atac, 
             region = "chr3-125880000-125950000", 
             features = "Ugt8a",
             annotation = TRUE,
             peaks = TRUE,
             tile = TRUE,
             links = TRUE)
VlnPlot(gc_rna, features = c("Ugt8a"),pt.size = 0)

# Identify peak markers
DefaultAssay(gc_atac) <- 'peaks'

peaks <- FindAllMarkers(gc_atac, min.pct = 0.25,  
                        only.pos = TRUE, 
                        logfc.threshold = 0.25,
                        test.use = 'LR', 
                        latent.vars = 'nCount_peaks')
write.csv(peaks, file = "~/ATAC-data/Results-Signac/gc-atac/peaks/Marker-gc.csv")

list <- c('Mitotic','Meiotic_interphase','preLeptotene','Leptotene','Zygotene')
for (i in list) {
  VlnPlot(
  object = gc_atac,
  features = rownames(subset(peaks, peaks$cluster==i))[1],
  pt.size = 0)
}

for (i in list) {
  FeaturePlot(
  object = gc_atac,
  features = rownames(subset(peaks, peaks$cluster==i)[1],
  pt.size = 1,
  max.cutoff = 'q95',
  cols = c('#D8BFD8','#B0E0E6','#0033FF')))
}


open_peaks <- rownames(peaks[peaks$avg_logFC > 0.25, ])
write.table(open_peaks,
            file = "~/ATAC-data/Results-Signac/gc-atac/peaks/open-peak.bed",
            row.names = F, 
            sep = " ",
            col.names = F,
            quote = F)

closest_peaks <- ClosestFeature(gc_atac, rownames(peaks))
write.csv(closest_peaks,file = "~/ATAC-data/Results-Signac/gc-atac/peaks/open-peak-gene.csv")

###------------------ Motif Analysis ------------------------ ###########################
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
set.seed(111)

gc_atac <- readRDS("~/ATAC-data/Results-Signac/atac_rds/gc_name_atac.rds")
DefaultAssay(gc_atac) <- 'peaks'

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(gc_atac),
  pwm = pfm,
  genome = 'mm10',
  use.counts = FALSE
)

# Create a new Motif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
gc_atac <- SetAssayData(
  object = gc_atac,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)

gc_atac[["peaks"]]

gc_atac <- RegionStats(object = gc_atac, genome = BSgenome.Mmusculus.UCSC.mm10)


-------------------------  1 -------  #默认assay 为 peaks
# Finding over-represented motifs 
gc_peaks <- FindMarkers(
  object = gc_atac,
  ident.1 = 'Zygotene',
  ident.2 = 'Mitotic', 
  only.pos = TRUE, 
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
  )


# get top differentially accessible peaks
top.gc.peak <- rownames(gc_peaks[gc_peaks$p_val < 0.005, ])

# find peaks open in Zyg
open.peaks <- AccessiblePeaks(gc_atac, idents = c("Zygotene", "Mitotic"))

# match the overall GC content in the peak set
peaks.matched <- MatchRegionStats(
  meta.feature = GetAssayData(gc_atac, assay = "peaks", slot = "meta.features")[open.peaks, ],
  query.feature = GetAssayData(gc_atac, assay = "peaks", slot = "meta.features")[top.gc.peak, ],
  n = 50000
  )

# test enrichment
enriched.motifs <- FindMotifs(
  object = gc_atac,
  features = top.gc.peak,
  background = peaks.matched
  )

head(enriched.motifs)
enriched.motifs <- filter(enriched.motifs, fold.enrichment > 1 & pvalue < 0.05) 
enriched.motifs <- enriched.motifs[order(enriched.motifs$fold.enrichment, decreasing= T), ]
write.csv(enriched.motifs,file = "~/ATAC-data/Result-scRNA/DEG-analysis/motif/Lep-preL.csv")

MotifPlot(
  object = gc_atac,
  motifs = rownames(enriched.motifs)[1:20]
)

# Computing motif activities 
gc_atac <- RunChromVAR(
  object = gc_atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)



###------  Building trajectories with Monocle 3 ------------- -------------
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)

set.seed(111)

gc_atac <- readRDS("~/ATAC-data/Results-Signac/atac_rds/gc_motif_atac.rds")

DefaultAssay(gc_atac) <- "peaks"

Germline.cds <- as.cell_data_set(gc_atac)

Germline.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(gc_atac[["RNA"]])
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
           trajectory_graph_segment_size = 1) & scale_color_gradient(low = '#FFBF00', high = 'purple')

gc_atac <- AddMetaData(
  object = gc_atac,
  metadata = Germline.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Germline")

FeaturePlot(gc_atac, "Germline", pt.size = 1) & scale_color_viridis_c()

# 使用monocle.3估计转录因子分化路径
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

a = pheatmap::pheatmap(agg_mat,scale="column", 
                   clustering_method="ward.D2",
                   colorRampPalette(c('#FFBF00','purple'))(255),
                   border_color = F)
# 挑选分化关键基因
clusters <- cutree(a$tree_row, k = 2)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_module"
table(clustering)

up_module = dplyr::filter(clustering,clustering$Gene_Clusters=='2')
up_module$module = rownames(up_module)

# 将gene_module_df 写出后在module列前添加Module以便使用dplyr操作
mod_peak = read.csv("~/mod.csv")

# 合并两数据
all = dplyr::left_join(mod_peak,up_module, by = 'module')

# 挑选指定module对应peak
up_peak_monocle =  dplyr::filter(all, all$Gene_Clusters =='2')

# 使用signac中ClosestFeature函数将peak转化为gene
peak_to_gene = ClosestFeature(gc_atac, up_peak_monocle$id)

# 使用clusterProfiler注释基因功能
ego<-clusterProfiler::enrichGO(peak_to_gene$gene_id,
              OrgDb      = org.Mm.eg.db,
              keyType    = 'ENSEMBL',
              ont        = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.05)

write.csv(ego,"~/ATAC-data/Results-Signac/gc-atac/monocle/peak_module_go.csv")

## peak pesudotime 
Germline.cds_lin <- Germline.cds[,is.finite(pseudotime(Germline.cds))]
plot_accessibility_in_pseudotime(Germline.cds_lin[c("chr6-34920019-34921077" )]) 



###-------------  Footprinting Analysis --------------------- ---------------------
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(111)

gc_atac = readRDS("~/ATAC-data/Results-Signac/atac_rds/gc_atac_motif.rds")
DefaultAssay(gc_atac) = 'peaks'

pwm <- getMatrixSet(x = JASPAR2020,
                    opts = list(species = 9606, all_versions = FALSE))

motif.positions <- matchMotifs(
                   pwms = pwm,
                   subject = granges(gc_atac),
                   out = 'positions',
                   genome = 'mm10')
  
motif <- CreateMotifObject(
         positions = motif.positions,
         pwm = pwm)
  
gc_atac <- SetAssayData(
           object = gc_atac,
           slot = 'motifs',
           new.data = motif) 
  
gc_atac <- Footprint(
       object = gc_atac,
       motif.name = c("MA0632.2","MA1650.1","MA0750.2","MA0131.2"),
       genome = BSgenome.Mmusculus.UCSC.mm10)

PlotFootprint(gc_atac, features = c("MA0632.2","MA1650.1","MA0750.2","MA0131.2")) 

###------------------- CCANs Cicero ------------------------- ------------------------
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
set.seed(111)

# convert to CellDataSet format and make the cicero object
gc.cds <- as.cell_data_set(x = gc_atac)
gc.cicero <- make_cicero_cds(gc.cds, reduced_coordinates = reducedDims(gc.cds)$UMAP)

genome <- seqlengths(gc_atac)
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# 共表达得分计算
conns <- run_cicero(gc.cicero, genomic_coords = genome.df, sample_num = 100)
write.table(conns, file="~/ATAC-data/Results-Signac/gc-atac/cicero/conns-score.txt", 
            sep = " ", quote = F, row.names = F, col.names = T)

# Building CCANs
ccans <- generate_ccans(conns)
write.table(ccans, file="~/ATAC-data/Results-Signac/gc-atac/cicero/ccans.txt", 
            sep = " ", quote = F, row.names = F, col.names = T)

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(gc_atac) <- links

# Annotation of peaks to genes 
closest_peaks <- ClosestFeature(gc_atac, rownames(ccans))
write.csv(closest_peaks, file = "~/ATAC-data/Results-Signac/gc-atac/cicero/ccans-genes.csv")


CoveragePlot(gc_atac, 
             region = "Tcfl5", 
             features = "Tcfl5",
             annotation = TRUE,
             peaks = TRUE,
             tile = TRUE,
             links = TRUE)


rm(list=ls())
sessionInfo()