## 差异表达分析
## jack zhang

library(Seurat)
library(Signac)
library(clusterProfiler)
library(org.Mm.eg.db)
set.seed(111)
setwd("~/ATAC-data/Result-scRNA/DEG-analysis/")

gc_rna <- readRDS("~/ATAC-data/Result-scRNA/RNA_rds/gc_name.rds")
gc_atac <- readRDS("~/ATAC-data/Results-Signac/atac_rds/gc_motif_atac.rds")

###------------ 关键转录因子挑选   RNA & peak-motif & motif score ---------------
#----1---- Differentially expressed genes
DEG = FindMarkers(gc_rna, 
                  ident.1 = "Leptotene", 
                  ident.2 = "preLeptotene", 
                  min.pct = 0.25,
                  logfc.threshold = 0)

# 挑选DEG
up_rna <- dplyr::filter(DEG,DEG$avg_logFC > 0 & DEG$p_val < 0.05)
down_rna <- dplyr::filter(DEG,DEG$avg_logFC < 0 & DEG$p_val < 0.05)
# 结果导出
write.csv(DEG, file = "DEG-Lep-preL.csv")


#----2---- Differentially expressed peaks
DefaultAssay(gc_atac) = 'peaks'

DEP = FindMarkers(gc_atac, 
                  ident.1 = "Leptotene", 
                  ident.2 = "preLeptotene", 
                  min.pct = 0.25,
                  logfc.threshold = 0,
                  test.use = 'LR', 
                  latent.vars = 'nCount_peaks')

# 挑选差异peak
up_atac = dplyr::filter(DEP,DEP$avg_logFC > 0 & DEP$p_val < 0.05)
down_atac = dplyr::filter(DEP,DEP$avg_logFC < 0 & DEP$p_val < 0.05)
# 结果导出
write.csv(DEP, file = "DEP-Lep-preL.csv")

# find peaks open in Lep
open.peaks <- AccessiblePeaks(gc_atac, idents = c("Leptotene", "preLeptotene"))
up.peak <- rownames(up_atac)
down.peak <- rownames(down_atac)

# match the overall GC content in the peak set
meta.feature <- GetAssayData(gc_atac, assay = "peaks", slot = "meta.features")

peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[up.peak, ],
  n = 50000
)

# test enrichment
enriched.motifs <- FindMotifs(
  object = gc_atac,
  features = up.peak,
  background = peaks.matched
)

# 写出结果
write.csv(enriched.motifs,file = "Motif_from_peak_Lep-preL.csv")

# 转化大小写
up_TFs = Hmisc::capitalize(tolower(enriched.motifs$motif.name))


#----3---- Differentially motif
DefaultAssay(gc_atac) = 'chromvar'

DEM = FindMarkers(gc_atac, 
                  ident.1 = "Leptotene", 
                  ident.2 = "preLeptotene", 
                  min.pct = 0.25,
                  logfc.threshold = 0,
                  test.use = 'LR', 
                  latent.vars = 'nCount_peaks')

# 挑选差异motif
up_motif = dplyr::filter(DEM,DEM$avg_logFC > 0 & DEM$p_val < 0.05)
down_motif = dplyr::filter(DEM,DEM$avg_logFC < 0 & DEM$p_val < 0.05)
# 写出结果
write.csv(DEM,file = "DEM_Lep-preL.csv")

# 转化大小写
up_motif$motif = rownames(up_motif)
a = dplyr::left_join(up_motif,enriched.motifs,by="motif")
up3 = Hmisc::capitalize(tolower(a$motif.name))


# RNA & ATAC & Motif 交集
b = intersect(rownames(up_rna), up_TFs)
c = intersect(b, up3)

write.csv(c,file = "Co_up_TFs.csv")

# RNA & ATAC 交集
peak_to_gene <- ClosestFeature(gc_atac, rownames(up_atac))
d = intersect(rownames(up_rna), peak_to_gene$gene_name)

# ATAC & Motif 交集
e = intersect(up3, peak_to_gene$gene_name)

# 韦恩图展示
library(eulerr)
vd <- euler(c(RNA = 708, ATAC = 70, Motif=0,
              "RNA&ATAC"=1,"RNA&Motif"=0,"ATAC&Motif"=549,
              "RNA&ATAC&Motif"=14))

plot(vd,
     fills = list(fill = c("#fbb4ae", "#b3cde3","red"), alpha = 0.8),
     labels = list(col = "black", font = 5),
     quantities = TRUE)

# --------- DEM 火山图 
# 设置阈值
DEM$threshold = as.factor(ifelse(DEM$p_val < 0.05 & abs(DEM$avg_logFC) >= 0, 
                                 ifelse(DEM$avg_logFC > 0 ,'Up','Down'),'NoSignifi'))
# 挑选前10位基因
DEM$lable=""
DEM=DEM[order(DEM$p_val),]
up=head(rownames(DEM)[which(DEM$threshold == "Up")],10)
down=head(rownames(DEM)[which(DEM$threshold == "Down")],10)
top10=c(as.character(up),as.character(down))
DEM$lable[match(top10,rownames(DEM))]=top10

# 火山图展示
ggplot(data = DEM, aes(x = avg_logFC, y = -log10(p_val), colour = threshold))+ 
  geom_point(alpha=0.7, size=2) +
  scale_color_manual(values=c("deepskyblue", "grey","#d94dff")) +
  geom_text(aes(label = lable), size =3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1))



###------------ 差异基因功能注释------

---------------- DEP---------------------
# 挑选差异peak
up_atac = dplyr::filter(DEP,DEP$avg_logFC > 0 & DEP$p_val < 0.05)
down_atac = dplyr::filter(DEP,DEP$avg_logFC < 0 & DEP$p_val < 0.05)

# DEP注释至临近Gene
peak_to_gene <- ClosestFeature(gc_atac, rownames(up_atac))

# 统计上调与下调个数
table(DEP$avg_logFC > 0 & DEP$p_val < 0.05)
table(DEP$avg_logFC < 0 & DEP$p_val < 0.05)

# peak_to_gene GO注释
ego <- enrichGO(peak_to_gene$gene_name,
                OrgDb      = org.Mm.eg.db,
                keyType    = 'SYMBOL',
                ont        = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
barplot(ego)
write.csv(ego,file = "GO_ATAC-zyg-lep.csv")


-------------- DEG -------------------
# DEG 火山图
DEG = FindMarkers(gc_rna, 
                  ident.1 = "Leptotene", 
                  ident.2 = "preLeptotene", 
                  min.pct = 0.25,
                  logfc.threshold = 0)

# 设置阈值
DEG$threshold = as.factor(ifelse(DEG$p_val < 0.05 & abs(DEG$avg_logFC) >= 0, 
                                 ifelse(DEG$avg_logFC > 0 ,'Up','Down'),'NoSignifi'))
# 挑选前10位基因
DEG$lable=""
DEG=DEG[order(DEG$p_val),]
up=head(rownames(DEG)[which(DEG$threshold == "Up")],10)
down=head(rownames(DEG)[which(DEG$threshold == "Down")],10)
top10=c(as.character(up),as.character(down))
DEG$lable[match(top10,rownames(DEG))]=top10

# 火山图展示
ggplot(data = DEG, aes(x = avg_logFC, y = -log10(p_val), colour = threshold))+ 
  geom_point(alpha=0.7, size=2) +
  scale_color_manual(values=c("deepskyblue", "grey","#d94dff")) +
  geom_text(aes(label = lable), size =3)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1))

write.csv(DEG, file = "DEG-Lep-preL.csv")

# 挑选DEG
up_rna <- dplyr::filter(DEG,DEG$avg_logFC > 0 & DEG$p_val < 0.05)
down_rna <- dplyr::filter(DEG,DEG$avg_logFC < 0 & DEG$p_val < 0.05)

# 统计上调与下调个数
table(DEG$avg_logFC > 0 & DEG$p_val < 0.05)
table(DEG$avg_logFC < 0 & DEG$p_val < 0.05)

# GO注释
ego <- enrichGO(rownames(up_rna),
                OrgDb      = org.Mm.eg.db,
                keyType    = 'SYMBOL',
                ont        = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

barplot(ego,showCategory = 10)
write.csv(ego, file = "GO_up_RNA-Lep-preL.csv")

# 统计peak_to_gene与DEG交集
gene = intersect(rownames(up_rna), peak_to_gene$gene_name)


# DEP注释至临近Gene
peak_to_gene <- ClosestFeature(gc_atac, rownames(up_atac))

# 统计上调与下调个数
table(DEP$avg_logFC > 0 & DEP$p_val < 0.05)
table(DEP$avg_logFC < 0 & DEP$p_val < 0.05)

# peak_to_gene GO注释
ego <- enrichGO(peak_to_gene$gene_name,
                OrgDb      = org.Mm.eg.db,
                keyType    = 'SYMBOL',
                ont        = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
barplot(ego)
write.csv(ego,file = "GO_ATAC-zyg-lep.csv")

# rna & atac 交集
common_gene <- intersect(rownames(up_rna), peak_to_gene$gene_name)

library(eulerr)

# 韦恩图展示
vd <- euler(c(RNA = (722-45), ATAC = (872-45), "RNA&ATAC" = 45))
plot(vd,
     fills = list(fill = c("#d94dff", "#d94dff"), alpha = 0.6),
     labels = list(col = "white", font = 4), 
     quantities = TRUE)

###------------ 细胞周期减分关键基因染色质开放展示----
a = readRDS("~/ATAC-data/Results-Signac/atac_rds/gc_ccan_atac.rds")

------ Stra8 -------------
CoveragePlot(a, 
             region = "chr6-34917000-34940000", 
             features = "Stra8",
             annotation = TRUE,
             peaks = TRUE,
             tile = TRUE,
             links = TRUE)
VlnPlot(gc_rna, features = c("Stra8"),pt.size = 0)

------ E2f1 -------------
CoveragePlot(a, 
             region = "chr2-154560000-154575000", 
             features = "E2f1",
             annotation = TRUE,
             peaks = TRUE,
             tile = TRUE,
             links = TRUE)
VlnPlot(gc_rna, features = c("E2f1"),pt.size = 0)


 
  

###------------ 候选core TFs 展示 ------
#---------- top10 motif ferom peak 
MotifPlot(
  object = gc_atac,
  motifs = rownames(enriched.motifs)[1:10]
)

#------- Candidate core TFs expession 
DotPlot(gc_rna, features = core.tfs,
        cols = c('#FF8C00','blue'),
        dot.scale = 10)

core.tfs = c('Stra8',
             'Tcfl5','Vezf1','Foxn3','Tfdp1','Batf3','Fos','Cebpg',
             'Jund','Gbx2','Sp1',
             'E2f1','E2f6',
             'E2f2','E2f8')


#### ---------- sessionInfo --------
rm(list=ls())
sessionInfo()