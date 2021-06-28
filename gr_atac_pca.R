############################################################
###--------------  Analysis for gr_atac_atac   -------------
###----------------    by jack_zhang_qau  ------------------
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(stringr)
set.seed(111)

counts <- Read10X_h5(filename = "~/ATAC-data/CellRanger-ATAC/Aggr/outs/filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "~/ATAC-data/CellRanger-ATAC/Aggr/outs/singlecell.csv",
  header = TRUE,
  row.names = 1)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = '~/ATAC-data/CellRanger-ATAC/Aggr/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
  )

gr_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
  )

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(gr_atac) <- annotations

#----------- Computing Q.C Metrics -----------
# compute nucleosome signal score per cell
gr_atac <- NucleosomeSignal(object = gr_atac)

# compute TSS enrichment score per cell
gr_atac <- TSSEnrichment(object = gr_atac, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
gr_atac$pct_reads_in_peaks <- gr_atac$peak_region_fragments / gr_atac$passed_filters * 100
gr_atac$blacklist_ratio <- gr_atac$blacklist_region_fragments / gr_atac$peak_region_fragments

gr_atac$high.tss <- ifelse(gr_atac$TSS.enrichment > 2.5, 'High', 'Low')
TSSPlot(gr_atac, group.by = 'high.tss') 

gr_atac$nucleosome_group <- ifelse(gr_atac$nucleosome_signal > 4, 'nucleosome_signal > 4', 'nucleosome_signal < 4')
FragmentHistogram(object = gr_atac, group.by = 'nucleosome_group',region = 'chr1-1-10000000')

VlnPlot(
  object = gr_atac,
  features = c('pct_reads_in_peaks', 
               'peak_region_fragments',
               'TSS.enrichment', 
               'blacklist_ratio', 
               'nucleosome_signal'),
  pt.size = 0,
  cols = '#7B68EE',
  ncol = 5
)

FeatureScatter(object = gr_atac,
              feature1 = "pct_reads_in_peaks",
              feature2 = "peak_region_fragments",
              cols = 'hotpink',pt.size = 0.5)+
              geom_vline(xintercept=c(30), lty=4,col="black",lwd=0.6)+
              geom_hline(yintercept = c(3000,40000),lty=4,col="black",lwd=0.6)

FeatureScatter(object = gr_atac,
               feature1 = "blacklist_ratio",
               feature2 = "peak_region_fragments",
               cols = 'hotpink',pt.size = 0.5)+
               geom_vline(xintercept=c(0,0.01), lty=4,col="black",lwd=0.6)+
               geom_hline(yintercept = c(3000,40000),lty=4,col="black",lwd=0.6)

FeatureScatter(object = gr_atac,
               feature1 = "nucleosome_signal",
               feature2 = "peak_region_fragments",
               cols = 'hotpink',pt.size = 0.5)+
               geom_vline(xintercept=c(0,4), lty=4,col="black",lwd=0.6)+
               geom_hline(yintercept = c(3000,40000),lty=4,col="black",lwd=0.6)

FeatureScatter(object = gr_atac,
               feature1 = "pct_reads_in_peaks",
               feature2 = "TSS.enrichment",
               cols = 'hotpink',pt.size = 0.5)+
               geom_vline(xintercept=30, lty=4,col="black",lwd=0.6)+
               geom_hline(yintercept = 2.5,lty=4,col="black",lwd=0.6)

FeatureScatter(object = gr_atac,
               feature1 = "nCount_peaks",
               feature2 = "nFeature_peaks",
               cols = 'hotpink',pt.size = 0.5)

FeatureScatter(object = gr_atac,
               feature1 = "DNase_sensitive_region_fragments",
               feature2 = "peak_region_fragments",
               cols = 'hotpink',pt.size = 0.5)            

FeatureScatter(object = gr_atac,
               feature1 = "enhancer_region_fragments",
               feature2 = "peak_region_fragments",
               cols = 'hotpink',pt.size = 0.5)

FeatureScatter(object = gr_atac,
               feature1 = "promoter_region_fragments",
               feature2 = "peak_region_fragments",
               cols = 'hotpink',pt.size = 0.5)

#-------------- View Low quality data ----------------------
FeaturePlot(gr_atac, features = "blacklist_region_fragments")

#-------------- High quality data --------------------------
gr_atac <- subset(
    x = gr_atac,
    subset = peak_region_fragments > 3000 &
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 30 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2.5
)

#--------------------- Add state info. ----------------------------
gr_atac$state = as.factor(ifelse(str_detect(names(Idents(gr_atac)), "[1-2]$") == "TRUE", ifelse(str_detect(names(Idents(gr_atac)),"2$") == "TRUE",'GR12','GR11'),'GR13'))

#--------- Normalization and linear dimensional reduction ---------
gr_atac <- RunTFIDF(gr_atac)
gr_atac <- FindTopFeatures(gr_atac, min.cutoff = 'q5')
gr_atac <- RunSVD(gr_atac)

DepthCor(gr_atac, n = 50) + geom_point(size = 3) 
ElbowPlot(gr_atac,ndims = 50,reduction = 'lsi')

gr_atac <- FindNeighbors(gr_atac, reduction = 'lsi', dims = 2:30)
gr_atac <- FindClusters(gr_atac, resolution = 0.6, algorithm = 3)
gr_atac <- RunUMAP(gr_atac, reduction = 'lsi', dims = 2:30, n.neighbors = 35)

DimPlot(gr_atac, label = TRUE, pt.size = 0.5) 
DimPlot(gr_atac, label = TRUE, pt.size = 0.5, 
        group.by = "state",
        cols = c("hotpink","#FFD700","deepskyblue"))

#----------------- View No. of cluster ----------------------
table(gr_atac$seurat_clusters)
barplot(table(gr_atac$seurat_clusters))

#------------ Output UMAP info. ---------------------------- 
umap = cbind("Barcode" = rownames(Embeddings(object = gr_atac, reduction = "umap")), Embeddings(object = gr_atac, reduction = "umap"))
write.table(umap, file="~/ATAC-data/Results-Signac/Aggr/Analysis/umap_atac.csv", sep = ",", quote = F, row.names = F, col.names = T)

#------------- Create a gene activity matrix  ----------------
gene.activities <- GeneActivity(gr_atac)

#----- add the gene activity matrix to the Seurat object------------ 
#-------------- as a new assay and normalize it --------------------
gr_atac[['RNA']] <- CreateAssayObject(counts = gene.activities)

gr_atac <- NormalizeData(
  object = gr_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(gr_atac$nCount_RNA)
  )

DefaultAssay(gr_atac) <- 'RNA'

#### Save ###
saveRDS(gr_atac,file = "~/ATAC-data/Results-Signac/atac_rds/gr_atac_RNA.rds")

### Germ cells
FeaturePlot(object = gr_atac,
            features = c( 'Dppa3','Ddx4','Dazl' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c('#D8BFD8','#B0E0E6','#FF00FF'))
### PG
FeaturePlot(object = gr_atac,
            features = c( 'Wnt4','Wnt6' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c('#D8BFD8','#B0E0E6','#0033FF'))

### Epithelial cells (间皮细胞)
FeaturePlot(object = gr_atac,
            features = c( 'Upk3b','Krt19' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c('#D8BFD8','#B0E0E6','#0033FF'))

### Interstitial cells (间质细胞)
FeaturePlot(object = gr_atac,
            features = c( 'Col1a2','Bgn' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c('#D8BFD8','#B0E0E6','#0033FF'))

### Endothelial cells (内皮细胞)
FeaturePlot(object = gr_atac,
            features = c( 'Pecam1','Kdr' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c('#D8BFD8','#B0E0E6','#0033FF'))

### immune cells
FeaturePlot(object = gr_atac,
            features = c( 'Cd52','Car2' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c('#D8BFD8','#B0E0E6','#0033FF'))

### Erythroid cells (红细胞)
FeaturePlot(object = gr_atac,
            features = c( 'Alas2','Cldn5','Cx3cr1' ),
            pt.size = 0.01,
            max.cutoff = 'q90',
            cols = c('#D8BFD8','#B0E0E6','#0033FF'))

#-------------- Integration with scRNA-seq ----------------
gr_rna <- readRDS("~/ATAC-data/Result-scRNA/RNA_rds/aggr_pca.rds")

#----------使用之前的variableFeatures   ------------
# gr_atac_rna <- FindVariableFeatures(object = gr_rna, nfeatures = 5000)

transfer.anchors <- FindTransferAnchors(
  reference = gr_rna,  
  query = gr_atac, 
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = gr_rna$celltype,
  weight.reduction = gr_atac[['lsi']],
  )

gr_atac <- AddMetaData(object = gr_atac, metadata = predicted.labels)

plot1 <- DimPlot(
  object = gr_rna,pt.size = 0.1,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE)  + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = gr_atac,pt.size = 0.1,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE)  + ggtitle('scATAC-seq')

plot1 + plot2

### 给ATAC数据命名
for(i in levels(gr_atac)) {
  cells_to_reid <- WhichCells(gr_atac, idents = i)
  newid <- names(sort(table(gr_atac$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(gr_atac, cells = cells_to_reid) <- newid
}

gr_atac$celltype <- Idents(gr_atac)

DimPlot(object = gr_atac, label = TRUE, pt.size = 0.5) 
barplot(table(gr_atac$celltype))

##### 展示不同群之间的峰 #######
DefaultAssay(gr_atac) <- 'peaks'

## Dazl
CoveragePlot(gr_atac, region = "chr17-50279393-50300000", features = "Dazl")
VlnPlot(gr_rna, features = c("Dazl"),pt.size = 0)

## Wnt4
CoveragePlot(gr_atac, region = "chr4-137270000-137300000", features = "Wnt4")
VlnPlot(gr_rna, features = c("Wnt4"),pt.size = 0)

## Wnt6
CoveragePlot(gr_atac, region = "chr1-74770000-74785000", features = "Wnt6")
VlnPlot(gr_rna, features = c("Wnt6"),pt.size = 0)

## Upk3b
CoveragePlot(gr_atac, region = "chr5-136037000-136044000", features = "Upk3b")
VlnPlot(gr_rna, features = c("Upk3b"),pt.size = 0)

## Bgn
CoveragePlot(gr_atac, region = "chrX-73481000-73496000", features = "Bgn")
VlnPlot(gr_rna, features = c("Bgn"),pt.size = 0)

## Kdr
CoveragePlot(gr_atac, region = "chr5-75940000-75985000", features = "Kdr")
VlnPlot(gr_rna, features = c("Kdr"),pt.size = 0)

## Cd52
CoveragePlot(gr_atac, region = "Cd52", features = "Cd52")
VlnPlot(gr_rna, features = c("Cd52"),pt.size = 0)


#----------------- Identify peak markers --------------------
DefaultAssay(gr_atac) <- 'peaks'

peaks <- FindAllMarkers(gr_atac, min.pct = 0.25,  
                        only.pos = TRUE, 
                        logfc.threshold = 0.25,
                        test.use = 'LR', 
                        latent.vars = 'nCount_peaks')

write.csv(peaks,file = "~/ATAC-data/Results-Signac/Aggr/Analysis/Peak_markers.csv")
barplot(table(peaks$cluster))

VlnPlot(
  object = gr_atac,
  features = rownames(subset(peaks, peaks$cluster=="Germ_cells"))[1],
  pt.size = 0)

FeaturePlot(
  object = gr_atac,
  features = rownames(subset(peaks, peaks$cluster=="Germ_cells"))[1],
  pt.size = 0.1,
  max.cutoff = 'q95',
  cols = c('#D8BFD8','#B0E0E6','#0033FF')
)

closest_peaks <- ClosestFeature(gr_atac, rownames(peaks))

write.table(open_peaks,
            file = "~/ATAC-data/Results-Signac/Aggr/Peaks/open-peak.bed",
            row.names = F, 
            sep = " ",
            col.names = F,
            quote = F)

write.csv(closest_peaks,file = "~/ATAC-data/Results-Signac/Aggr/Analysis/peak-marker-gene.csv")

# 提取特定cluster，继续后续分析。
gc_atac <- subset(gr_atac, idents ="Germ_cells")
saveRDS(gc_atac, "~/ATAC-data/Results-Signac/atac_rds/gc_atac.rds")



#----------------------- Motif 分析 -----------------------------
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)

gr_atac = readRDS("~/ATAC-data/Results-Signac/atac_rds/gr_atac_RNA.rds")
DimPlot(gr_atac, label = TRUE, pt.size = 1) 

#-----------------------  The Motif class------------------------------------
#--- Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# ---------Scan the DNA sequence of each peak for the presence of each motif--
motif.matrix <- CreateMotifMatrix(
  features = granges(gr_atac),
  pwm = pfm,
  genome = 'mm10',
  use.counts = FALSE
)

# ----------------Create a new Motif object to store the results--------------
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
gr_atac <- SetAssayData(
  object = gr_atac,
  assay = 'peaks',
  slot = 'motifs',
  new.data = motif
)
gr_atac[["peaks"]]


gr_atac <- RegionStats(object = gr_atac, genome = BSgenome.Mmusculus.UCSC.mm10)

#------------------- Finding over-represented motifs ----------------
gc_peaks <- FindMarkers(
  object = gr_atac,
  ident.1 = 'Germ_cells',
  only.pos = TRUE, 
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)


# --------------get top differentially accessible peaks-----------
top.gc.peak <- rownames(gc_peaks[gc_peaks$p_val < 0.005, ])

# ---------------find peaks open in GCs---------------------------
open.peaks <- AccessiblePeaks(gr_atac, idents = "Germ_cells")

# ---------------match the overall GC content in the peak set-----
peaks.matched <- MatchRegionStats(
  meta.feature = GetAssayData(gr_atac, assay = "peaks", slot = "meta.features")[open.peaks, ],
  regions = top.gc.peak,
  n = 50000
)

# ----------------- test enrichment ------------------------------
enriched.motifs <- FindMotifs(
  object = gr_atac,
  features = top.gc.peak,
  background = peaks.matched
)

knitr::kable(head(enriched.motifs))
write.csv(enriched.motifs,file = "~/ATAC-data/Results-Signac/Aggr_PCA_atac/motifs/gc_motifs.csv")

MotifPlot(
  object = gr_atac,
  motifs = head(rownames(enriched.motifs))
)

#------------------- Computing motif activities -------------------
gr_atac <- RunChromVAR(
  object = gr_atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(gr_atac) <- 'chromvar'

#-----------------------cluster motif -----------------------------
motifs <- FindAllMarkers(gr_atac, min.pct = 0.3, 
                        test.use = 'LR', 
                        latent.vars = 'peak_region_fragments')

write.csv(motifs, file = "~/ATAC-data/Results-Signac/Aggr_PCA_atac/motifs/all-marker_motifs.csv")

# look at the activity of GCs's motif
FeaturePlot(
  object = gr_atac,
  features = "MA0003.4",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

#--------------- Building trajectories with Monocle 3 ----------------
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(EnsDb.Mmusculus.v79)

DefaultAssay(gr_atac) <- "peaks"

Germline.cds <- as.cell_data_set(gr_atac)
Germline.cds <- cluster_cells(cds = Germline.cds, reduction_method = "UMAP")
Germline.cds <- learn_graph(Germline.cds, use_partition = TRUE)
Germline.cds <- order_cells(Germline.cds, reduction_method = "UMAP")

plot_cells(
  cds = Germline.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, cell_size = 1, graph_label_size = 5)

gc_atac <- AddMetaData(
  object = gc_atac,
  metadata = Germline.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Germline")

FeaturePlot(gc_atac, "Germline", pt.size = 1) & scale_color_viridis_c()





rm(list=ls())

sessionInfo()