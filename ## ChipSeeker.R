#---------- Peak visualization --------------
#---------- by Jack_Zhang_QAU ---------------
library(ChIPseeker)
library(ggplot2)
library(ggupset)
library(ggimage)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
set.seed(111)

#----------- input bed data ---------------------
peak <- readPeakFile("~/ATAC-data/CellRanger-ATAC/Aggr/outs/peaks.bed")


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

# 
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="#E680FF")

plotAvgProf(tagMatrix,  xlim=c(-3000, 3000),  
            xlab="Genomic Region (5'->3')",  
            ylab = "Read Count Frequency")

plotAvgProf2(peak, TxDb=txdb, 
             upstream = 3000, downstream = 3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)

#--------------- peak Annotation -----------------
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
plotAnnoPie(peakAnno)
upsetplot(peakAnno,vennpie=TRUE)

#---------------- nearest gene annotation --------------------------
#指peak最近的基因：
#不管peak落在内含子、基因间区还是其他位置，
#按照peak相对于转录起始位点的距离，都能找到一个离它最近的基因

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

#----------------- Peak to gene annotation & KEGG ------------------
genes <- seq2gene(peak, c(-1000, 3000), 3000, TxDb=txdb)

cc = compareCluster(geneClusters = genes, 
                    fun="enrichKEGG", organism="mmu")

dotplot(cc, showCategory=10)

# -------------- visualization of peak length < 1000 ----------------
# peaksLength=abs(peakAnno_df$end-peakAnno_df$start)

peaksinfo <- as.data.frame(peak)
peaksLength <- peaksinfo[peaksinfo$width < 1000,]
hist(peaksLength, breaks =50, col = "lightblue", xlim=c(0,1000), xlab = "peak length", main="Histogram of peak length")






## ------  Cluster specific peak Annotation ----
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
set.seed(111)

germ = readPeakFile("~/ATAC-data/Results-Signac/ChipSeeker/peak_bed/Germ_cells.bed")
preg = readPeakFile("~/ATAC-data/Results-Signac/ChipSeeker/peak_bed/PreG.bed")
immune = readPeakFile("~/ATAC-data/Results-Signac/ChipSeeker/peak_bed/Immune.bed")
blood = readPeakFile("~/ATAC-data/Results-Signac/ChipSeeker/peak_bed/Bllod.bed")
inter = readPeakFile("~/ATAC-data/Results-Signac/ChipSeeker/peak_bed/Interstitial.bed")
endo = readPeakFile("~/ATAC-data/Results-Signac/ChipSeeker/peak_bed/Endo.bed")
epi = readPeakFile("~/ATAC-data/Results-Signac/ChipSeeker/peak_bed/Epithelial.bed")

peaks = list(germ,preg,immune,blood,inter,endo,epi)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# heatmap
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000),color = rainbow(length(peaks)))

# Annotation to genes
peakAnnoList <- lapply(peaks, annotatePeak, 
                       TxDb=txdb,
                       tssRegion=c(-3000, 3000), 
                       flankDistance = 5000, 
                       verbose=FALSE,
                       addFlankGeneInfo=TRUE,
                       annoDb="org.Mm.eg.db")

plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList,title="Peaks relative to TSS")

for (i in 1:7) {
  plotAvgProf(tagMatrixList[i], xlim=c(-3000, 3000), 
            conf = 0.95, resample = 1000)
}



# Create a list with genes from each sample
gene = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# Run GO enrichment analysis 
ego <- enrichGO(gene = gene[[1]], 
                keyType = "ENTREZID", 
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)


###### TCFL5 activable gene ##########
TCFL5-peak = readPeakFile("~/TCFL5.bed")

peakAnno <- annotatePeak(TCFL5-peak, tssRegion=c(-3000,3000), 
                         TxDb=txdb, annoDb="org.Mm.eg.db")

peakAnno_df <- as.data.frame(peakAnno)

TCFL5_activatable_gene = dplyr::filter(peakAnno_df, peakAnno_df$distanceToTSS < 3000) 




