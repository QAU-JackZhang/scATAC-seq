############################################################
###---- Analysis for CellChat   -------------
###---------------- by jack_zhang_qau  ---------------------
library(iTALK)
library(Hmisc)
library(stringr)
library(Seurat)
set.seed(111)

gr_rna = readRDS("~/ATAC-data/Result-scRNA/RNA_rds/aggr_pca.rds")

cell.meta <- subset(gr_rna@meta.data, select=c("state","celltype"))
names(cell.meta) <- c("compare_group", "cell_type")
cell.expr <- data.frame(t(as.matrix(gr_rna@assays$RNA@counts)), check.names=F)
data.italk <- merge(cell.expr, cell.meta, by=0)
rownames(data.italk) <- data.italk$Row.names
data.italk <- data.italk[,-1]

##设置绘图颜色和其他变量
mycolor <- c("#FFA500","#FF2400","#22C32E","#007FFF","#F400A1","#0000FF","#8CE600")
cell_type <- unique(data.italk$cell_type)
cell_col <- structure(mycolor[1:length(cell_type)], names=cell_type)
#通讯类型变量
comm_list <- c('growth factor','other','cytokine')
#圈图展示配体-受体对的数量
PN=100

path <- "~/italk/"


MouseSymbol <- colnames(data.italk)
HumanSymbol <- capitalize(toupper(MouseSymbol))
HumanSymbol[18402:18403]=c("compare_group", "cell_type")
colnames(data.italk) = HumanSymbol


#------------------
## 查看某个样品中
data1 <- subset(data.italk, subset=data.italk$compare_group=="GR11")
highly_exprs_genes <- rawParse(data1, top_genes=100, stats="mean")

## 查看总
highly_exprs_genes <- rawParse(data.italk, top_genes=100, stats="mean")

res<-NULL

for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  res_cat <- res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs, decreasing=T),]
  write.csv(res_cat, paste0(path,'LRpairs_Overview_',comm_type,'.csv'))
  #绘制细胞通讯关系网络图
  pdf(paste0(path,'LRpairs_Overview_',comm_type,'_.pdf'), width=10, height=10)
  NetView(res_cat, 
          col=cell_col, 
          vertex.label.cex=1.2, 
          edge.label.cex=1, 
          vertex.size=20, 
          arrow.width=1, 
          edge.max.width=2, 
          margin = 0.2)
  dev.off()
  #绘制topPN的配体-受体圈图
  pdf(paste0(path,'LRpairs_Overview_',comm_type,'_circ.pdf'), width=10, height=10)
  LRPlot(res_cat[1:PN,],
         datatype='mean count',
         link.arr.lwd=res_cat$cell_from_mean_exprs[1:PN],
         link.arr.width=res_cat$cell_to_mean_exprs[1:PN], 
         link.arr.col = 'skyblue',
         text.vjust = "0.35cm")
  dev.off()
  res<-rbind(res,res_cat)
}


res <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs, decreasing=T),]
write.csv(res, paste0(path,'LRpairs_Overview.csv'))
res.top <- res[1:PN,]
write.csv(res.top, paste0(path,'res_top.csv'))

NetView(res, 
        col=cell_col, 
        vertex.label.cex=1.2, 
        edge.label.cex=0.9, 
        vertex.size=20, 
        arrow.width=1, 
        edge.max.width=3, 
        margin=0.2)


LRPlot(res.top, 
       datatype='mean count', 
       link.arr.lwd=res.top$cell_from_mean_exprs,
       link.arr.width=res.top$cell_to_mean_exprs,
       link.arr.col = '#FFB3E6'
       )

#-------------- 差异分析 -----------
data2 <- subset(data.italk, subset=data.italk$compare_group=="GR12"|data.italk$compare_group=="GR13")
deg_all = NULL
for (i in unique(data2$cell_type)){
      a = dplyr::filter(data2, data2$cell_type==i)
      deg = DEG(a,method='DESeq2',contrast = c('GR12','GR13'))
      deg_all = rbind(deg_all,deg)
}


res<-NULL
for(comm_type in comm_list){
  #多个细胞类型之间显著表达的配体-受体，结果会过滤同一细胞类型内的配体-受体关系
  res_cat <- FindLR(deg_all, datatype='DEG', comm_type=comm_type)
  res_cat <- res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC, decreasing=T),]
  res<-rbind(res,res_cat)
}

##所有配体-受体分类一起作图
res<-res[order(res$cell_from_logFC*res$cell_to_logFC, decreasing=T),]


NetView(res, 
        col=cell_col, 
        vertex.label.cex=1.2, 
        edge.label.cex=0.9, 
        vertex.size=20, 
        arrow.width=1, 
        edge.max.width=5, 
        margin = 0.2)

## 此函数有bug  
LRPlot(res, 
       datatype='DEG', 
       link.arr.lwd=res$cell_from_logFC,
       link.arr.width=res$cell_to_logFC)

###### 使用以下代码

#When you've got your object "res", just run my codes and you'll get the plot######
library(tidyr)
library(dplyr)
library(circlize)

###first, set the parameter
data <- res[1:100,]
datatype<-'DEG'
gene_col = NULL 
transparency = 0.5 
link.arr.lwd = res[1:100,]$cell_from_logFC
link.arr.lty = NULL 
link.arr.col = NULL
link.arr.width = res[1:100,]$cell_to_logFC 
link.arr.type = NULL 
facing = "clockwise" 
cell_col = cell_col
print.cell = TRUE
track.height_1 = uh(2,"mm")
track.height_2 = uh(12, "mm")
annotation.height_1 = 0.01
annotation.height_2 = 0.01 
text.vjust = "0.4cm" 

###Second, get the correct order
###Here is the bug. The object "genes" only contains genes but not the "cell"
cell_group <- unique(c(data$cell_from, data$cell_to))
genes <- c(structure(data$ligand, names = data$cell_from), 
           structure(data$receptor, names = data$cell_to))
genes <- genes[!duplicated(paste(names(genes), genes))]
genes <- genes[order(names(genes))]
##So We creat a new matrix  
cell_group <- unique(c(data$cell_from, data$cell_to))
genes_2 <- c(structure(cbind(paste(data$cell_from, data$ligand)), names =data$cell_from ), 
             structure(cbind(paste(data$cell_to, data$receptor)), names = data$cell_to ))
genes_2 <- genes_2[!duplicated(paste(names(genes_2), genes_2))]
genes_2 <- genes_2[order(names(genes_2))]
genes_matrix <- data.frame(names=names(genes_2),genes_2=genes_2)
genes_matrix <- separate(genes_matrix,genes_2,into = c("cell","gene"),sep = " ")
genes_matrix$genes_2 <- genes_2




##Third, set the parameters of the lines in the plot. Just run the code
if (is.null(link.arr.lty)) {
  if (datatype == "mean count") {
    link.arr.lty = "solid"
  }else if (datatype == "DEG") {
    link.arr.lty = structure(ifelse(data$cell_from_logFC == 
                                      1e-04, "dashed", "solid"), names = paste(data$cell_from, 
                                                                               data$receptor))
  }else {
    print("invalid datatype")
  }
}

if (is.null(link.arr.col)) {
  if (datatype == "mean count") {
    data <- data %>% mutate(link_col = "black")
  }else if (datatype == "DEG") {
    data <- data %>% mutate(link_col = ifelse(cell_from_logFC == 
                                                1e-04, ifelse(cell_to_logFC > 0, "#d73027", 
                                                              "#00ccff"), ifelse(cell_to_logFC == 1e-04, 
                                                                                 ifelse(cell_from_logFC > 0, "#d73027", 
                                                                                        "#00ccff"), ifelse(cell_from_logFC > 
                                                                                                             0, ifelse(cell_to_logFC > 0, "#d73027", 
                                                                                                                       "#dfc27d"), ifelse(cell_to_logFC > 0, 
                                                                                                                                          "#9933ff", "#00ccff")))))
  }else {
    print("invalid datatype")
  }
}else {
  data$link_col = link.arr.col
}

if (is.null(link.arr.type)) {
  if (datatype == "mean count") {
    link.arr.type = "triangle"
  }else if (datatype == "DEG") {
    link.arr.type = structure(ifelse(data$cell_to_logFC == 
                                       1e-04, "ellipse", "triangle"), names = paste(data$cell_from, 
                                                                                    data$receptor))
  }else {
    print("invalid datatype")
  }
}

if (is.null(gene_col)) {
  comm_col <- structure(c("#99ff99", "#99ccff", 
                          "#ff9999", "#ffcc99"), names = c("other", 
                                                           "cytokine", "checkpoint", "growth factor"))
  gene_col <- structure(c(comm_col[data$comm_type], rep("#073c53", 
                                                        length(data$receptor))), names = c(data$ligand, data$receptor))
}
if (is.null(cell_col)) {
  cell_col <- structure(randomColor(count = length(unique(names(genes))), 
                                    luminosity = "dark"), names = unique(names(genes)))
}
if (is.null(link.arr.lwd)) {
  data <- data %>% mutate(arr_width = 1)
}else if (max(abs(link.arr.lwd)) - min(abs(link.arr.lwd)) == 
          0 && all(link.arr.lwd != 1e-04)) {
  data <- data %>% mutate(arr_width = ifelse(abs(link.arr.lwd < 
                                                   5), abs(link.arr.lwd), 5))
}else {
  data <- data %>% mutate(arr_width = ifelse(link.arr.lwd == 
                                               1e-04, 2, 1 + 5/(max(abs(link.arr.lwd)) - min(abs(link.arr.lwd))) * 
                                               (abs(link.arr.lwd) - min(abs(link.arr.lwd)))))
}

if (length(cell_group) != 1) {
  gap.degree <- do.call("c", lapply(table(names(genes)), 
                                    function(i) c(rep(1, i - 1), 8)))
}else {
  gap.degree <- do.call("c", lapply(table(names(genes)), 
                                    function(i) c(rep(1, i))))
}
circos.par(gap.degree = gap.degree)
if (length(gene_col) == 1) {
  grid.col = gene_col
}else {
  grid.col = gene_col[genes]
  names(grid.col) <- paste(names(genes), genes)
}


if (is.null(link.arr.width)) {
  data <- data %>% mutate(link.arr.width = data$arr_width/10)
}else if (max(abs(link.arr.width)) - min(abs(link.arr.width)) == 
          0 && all(link.arr.width != 1e-04)) {
  data <- data %>% mutate(link.arr.width = ifelse(abs(link.arr.width) < 
                                                    0.5, abs(link.arr.width), 0.5))
}else {
  data <- data %>% mutate(link.arr.width = ifelse(link.arr.width == 
                                                    1e-04, 0.2, (1 + 5/(max(abs(link.arr.width)) - min(abs(link.arr.width))) * 
                                                                   (abs(link.arr.width) - min(abs(link.arr.width))))/10))
}



###Forth,chordDiagram the plot. 
chordDiagram(as.data.frame(cbind(paste(data$cell_from, data$ligand), 
                                 paste(data$cell_to, data$receptor))), 
             order = genes_matrix$genes_2, ####We replace "paste(names(genes),genes)" with "genes_matrix$genes_2"
             grid.col = grid.col, 
             transparency = transparency, 
             directional = 1, 
             direction.type = "arrows", 
             link.arr.lwd = data$arr_width, 
             link.arr.lty = link.arr.lty, 
             link.arr.type = link.arr.type, 
             link.arr.width = data$link.arr.width, 
             link.arr.col = data$link_col, 
             col = "#00000000", 
             annotationTrack = c("grid"), 
             preAllocateTracks = list(list(track.height = track.height_1),list(track.height = track.height_2)), 
             annotationTrackHeight = c(annotation.height_1,annotation.height_2))

###Add the label of genes
circos.trackPlotRegion(track.index = 2, panel.fun = function(x,y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = genes[get.cell.meta.data("sector.numeric.index")]
  circos.text(mean(xlim), mean(ylim), sector.index, col = "black", 
              cex = 0.7, facing = facing, niceFacing = TRUE)
}, bg.border = 0)



###Add the label of cells
for (c in unique(genes_matrix$names)) {
  genes_matrix_2 = genes_matrix[which(genes_matrix$names==c),]
  highlight.sector(sector.index = genes_matrix_2$genes_2, track.index = 1, 
                   col = ifelse(length(cell_col) == 1, cell_col, 
                                cell_col[c]), text = genes_matrix_2$cell, text.vjust = text.vjust, 
                   niceFacing = TRUE, lwd = 1)
}  

rm(list=ls())