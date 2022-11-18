library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(scran)
library(viridis)
library(ggforce)
library(gghalves)
library(scDblFinder)
library(ggpubr)
library(SingleCellExperiment)
library(Matrix)
library(scRecover)
library(BiocParallel)
library(SAVER)
library(monocle3)
library(enrichR)
library(clusterProfiler)
library(ggsci)

custom_colors <- list()
custom_colors$c15 <- c('#ED4C67', '#ff6348', '#F79F1F','#f9ca24',  '#A3CB38', 
                       '#20bf6b', '#4b7bec', '#2d98da', '#12CBC4', '#1289A7', 
                       '#778ca3', '#82589F', '#D980FA', '#f78fb3', '#f8a5c2')
custom_colors$c4 <- c('#1289A7', '#ffda79', '#a4b0be', '#ba4b7f')


#data pre-processing and imputation ####
file <- c("/disk1/snRNAseq/mouse/mouse_result/Mouse_hypo_C_01_02_Nuc/outs/filtered_feature_bc_matrix", 
          "/disk1/snRNAseq/mouse/mouse_result/Mouse_hypo_C_04_Nuc/outs/filtered_feature_bc_matrix", 
          "/disk1/snRNAseq/mouse/mouse_result/Mouse_hypo_T_01_03_Nuc/outs/filtered_feature_bc_matrix",
          "/disk1/snRNAseq/mouse/mouse_result/Mouse_hypo_T_04_Nuc/outs/filtered_feature_bc_matrix")


names(file) <- c("C1", "C4", "T1", "T4")
set.seed(73)

param <- MulticoreParam(workers = 6, progressbar = TRUE)
register(param)

data.list <- list()
data.list.scr <- list()
qcplot.list <- list()
for (i in 1:length(file)){
  data <- Read10X(data.dir = file[i])
  
  cells <- tibble(
    cell = colnames(data),
    nCount = colSums(data),
    nFeature = colSums(data != 0))
  
  cells$percent_ribo <- Matrix::colSums(data[c("Gm42418", "AY036118"),]) / Matrix::colSums(data)
  cells$percent_mt <- Matrix::colSums(data[grep("^mt-", rownames(data)),]) / Matrix::colSums(data)
  
  #doublet removal
  sce <- SingleCellExperiment(assays = list(counts = data))
  sce <- scDblFinder(sce)
  cells$multiplet_class <- colData(sce)$scDblFinder.class
  table(cells$multiplet_class)
  
  
  #prefiltering
  median_nCount <- median(cells$nCount)
  mad_nCount <- mad(cells$nCount)
  median_nFeature <- median(cells$nFeature)
  mad_nFeature <- mad(cells$nFeature)
  median_percent_ribo <- median(cells$percent_ribo)
  mad_percent_ribo <- mad(cells$percent_ribo)
  
  thresholds_nCount <- c(0, 7000)
  thresholds_nFeature <- c(0, 3000)
  thresholds_percent_ribo <- c(0, 0.1)
  thresholds_percent_mt <- c(0, 0.1)
  
  p1 <- ggplot(cells, aes(nCount, nFeature, color = percent_ribo)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = thresholds_nFeature, color = 'red') +
    geom_vline(xintercept = thresholds_nCount, color = 'red') +
    scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
    scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
    theme_bw() +
    scale_color_viridis(
      name = 'Percent Ribo\ntranscripts',
      limits = c(0,0.1),
      labels = scales::percent,
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
    )
  
  
  p2 <- ggplot(cells, aes(nCount, nFeature, color = percent_mt)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = thresholds_nFeature, color = 'red') +
    geom_vline(xintercept = thresholds_nCount, color = 'red') +
    scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
    scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
    theme_bw() +
    scale_color_viridis(
      name = 'Percent Mt\ntranscripts',
      limits = c(0,0.1),
      labels = scales::percent,
      guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
    )
  
  
  cells_filtered <- cells %>%
    dplyr::filter(
      nCount >= thresholds_nCount[1],
      nCount <= thresholds_nCount[2],
      nFeature >= thresholds_nFeature[1],
      nFeature <= thresholds_nFeature[2],
      percent_ribo >= thresholds_percent_ribo[1],
      percent_ribo <= thresholds_percent_ribo[2]
    )
  
  cells_to_keep <- cells_filtered$cell
  length(cells_to_keep)
  
  temp.object <- CreateSeuratObject(counts = data[,cells_to_keep], project = names(file)[i], min.cells = 3, min.features = 10)
  temp.object[["percent.ribo"]] <- PercentageFeatureSet(temp.object, features = c("Gm42418", "AY036118"))
  temp.object[["percent.mito"]] <- PercentageFeatureSet(temp.object, pattern = "^mt-")
  data.list[[paste0("data.",names(file)[i])]] <- temp.object
  qcplot.list[[paste0("data.",names(file)[i])]] <- p1+p2
  
  
  temp.object <- scRecover(counts= as.matrix(data[,cells_to_keep]), Kcluster= 32, SAVER= T, UMI = F, parallel = T, BPPARAM=param, outputDir = names(file)[i])
  temp.object <- CreateSeuratObject(counts = read.csv(file= paste0(names(file)[i],"scRecover+SAVER.csv"), row.names = 1), project = names(file)[i], min.cells = 3, min.features = 10)
  temp.object[["percent.ribo"]] <- PercentageFeatureSet(temp.object, features = c("Gm42418", "AY036118"))
  temp.object[["percent.mito"]] <- PercentageFeatureSet(temp.object, pattern = "^mt-")
  data.list.scr[[paste0("data.",names(file)[i])]] <- temp.object
}

#https://www.biorxiv.org/content/10.1101/657726v1.full <- justification for 18S rRNA


ggarrange(qcplot.list[[1]], qcplot.list[[2]], qcplot.list[[3]], qcplot.list[[4]], ncol = 1, nrow=4)
ggsave("./images/qcplot.pdf", width = 8, height = 12)

data.list[[1]]$group <- rep("Con", ncol(data.list[[1]]))
data.list[[2]]$group <- rep("Con", ncol(data.list[[2]]))
data.list[[3]]$group <- rep("Tx", ncol(data.list[[3]]))
data.list[[4]]$group <- rep("Tx", ncol(data.list[[4]]))

data.list[[1]]$batch <- rep("run1", ncol(data.list[[1]]))
data.list[[2]]$batch <- rep("run2", ncol(data.list[[2]]))
data.list[[3]]$batch <- rep("run1", ncol(data.list[[3]]))
data.list[[4]]$batch <- rep("run2", ncol(data.list[[4]]))

data.list.scr[[1]]$group <- rep("Con", ncol(data.list.scr[[1]]))
data.list.scr[[2]]$group <- rep("Con", ncol(data.list.scr[[2]]))
data.list.scr[[3]]$group <- rep("Tx", ncol(data.list.scr[[3]]))
data.list.scr[[4]]$group <- rep("Tx", ncol(data.list.scr[[4]]))

data.list.scr[[1]]$batch <- rep("run1", ncol(data.list.scr[[1]]))
data.list.scr[[2]]$batch <- rep("run2", ncol(data.list.scr[[2]]))
data.list.scr[[3]]$batch <- rep("run1", ncol(data.list.scr[[3]]))
data.list.scr[[4]]$batch <- rep("run2", ncol(data.list.scr[[4]]))


#making seurat object non inputed ####
data.list <- lapply(data.list, function(x) {SCTransform(x, vst.flavor="v2", vars.to.regress =)})
features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features, normalization.method = "SCT")
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT") 
combined <- RunPCA(combined, npcs = 50, verbose = FALSE) %>% FindNeighbors(combined, reduction = "pca", dims = 1:50)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:50, n.neighbors = 50)
combined <- FindClusters(combined, resolution = 0.5, algorithm = 4)
head(Idents(combined), 5)

umap1 <- DimPlot(combined, reduction = "umap", label = T, cols = custom_colors$c15 ) & NoAxes() & NoLegend()
umap2 <- DimPlot(combined, reduction = "umap", label = T, cols = custom_colors$c15, split.by = "orig.ident", ncol=2) & NoAxes() & NoLegend()
umap3 <- DimPlot(combined, reduction = "umap", label = F, cols = custom_colors$c4, group.by ="orig.ident" ) & NoAxes() 
umap4 <- DimPlot(combined, reduction = "umap", label = F, cols = custom_colors$c4, group.by = "group" ) & NoAxes()

umap.all <- ggarrange(umap1, umap2, umap3, umap4, ncol = 2, nrow=2)
umap.all 
ggsave("./images/umap_qc_no_imput.pdf", width = 3, height = 3, scale = 2)

#batch effect QC
a <- data.frame(file = combined@meta.data$orig.ident, batch = combined@meta.data$batch, ribo = combined@meta.data$percent.ribo, mito = combined@meta.data$percent.mito)
p1 <- ggplot(data = a, aes(x= batch, y= ribo, col = file)) + geom_boxplot() + scale_color_npg() + theme_classic()
p2 <- ggplot(data = a, aes(x= batch, y= mito, col = file)) + geom_boxplot() + scale_color_npg() + theme_classic()
ggarrange(p1,p2)
ggsave("./images/batch_effect.pdf", width =6, height = 3)


combined@active.assay <- "RNA"
combined <- NormalizeData(combined)
combined <- SCTransform(combined, vst.flavor="v2", vars.to.regress = "batch") %>% PrepSCTFindMarkers()

FeaturePlot(combined, c("P2ry12", "C1qa", "Cx3cr1", "Tmem119", "Snap25", "Gad1", "Slc17a6", "Aqp4",  "Slc1a3", "Pdgfra", "Olig2", "Mbp", "Plp1", "Fgf10", "Rax", "Col23a1")) & NoLegend() & NoAxes()
ggsave("./images/featureplot_no_imput.pdf", scale = 2)

saveRDS(combined, "combined.rds")

View( all.markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC))



#seurat object scr ####
data.list.scr <- lapply(data.list.scr, function(x) {SCTransform(x, vst.flavor="v2", vars.to.regress = )})
features <- SelectIntegrationFeatures(object.list = data.list.scr)
data.list.scr <- PrepSCTIntegration(object.list = data.list.scr, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = data.list.scr, anchor.features = features, normalization.method = "SCT")
combined.scr <- IntegrateData(anchorset = anchors, normalization.method = "SCT") 
combined.scr <- RunPCA(combined.scr, npcs = 50, verbose = FALSE) %>% FindNeighbors(combined, reduction = "pca", dims = 1:50)
combined.scr <- RunUMAP(combined.scr, reduction = "pca", dims = 1:50, n.neighbors = 50, min.dist = )
combined.scr <- FindClusters(combined.scr, resolution = 0.5, algorithm = 4)
head(Idents(combined.scr), 5)

labels <- forcats::fct_collapse(combined.scr$seurat_clusters, "GABA" = as.character(c(1,5,11,12,18,19,21,23,27,29)), "Glut" = as.character(c(2,6,7,8,10,13,16,17,22,24,26,28,30)),
                                "Oligo" = "3", "Ast" = "4", "OPC" = "9", "Tan" = c("14","20"), "Mg" = "15", "Peri" = "25")
combined.scr <- AddMetaData(combined.scr, labels, col.name = "celltype")
Idents(combined.scr) <- "celltype"

a <- data.frame(celltype = combined.scr$celltype, group = combined.scr$group)
a %>% group_by(group, celltype) %>% tally()

ggplot(a, aes(x= group, fill = celltype)) +geom_bar() + scale_fill_manual(values = colorspace::darken(custom_colors$c14, 0.05)) +coord_flip() + theme_classic()
ggsave("./images/cellnumber.pdf", scale = 2)
round(prop.table(table(Idents(combined.scr), combined.scr$group), margin = 2) * 100,2)

umap1 <- DimPlot(combined.scr, reduction = "umap", label = T, cols = custom_colors$c14 ) & NoAxes() & NoLegend()
umap2 <- DimPlot(combined.scr, reduction = "umap", label = T, cols = custom_colors$c14, split.by = "orig.ident", ncol=2) & NoAxes() & NoLegend()
umap3 <- DimPlot(combined.scr, reduction = "umap", label = F, cols = custom_colors$c4, group.by ="orig.ident" ) & NoAxes() 
umap4 <- DimPlot(combined.scr, reduction = "umap", label = F, cols = custom_colors$c14, group.by = "group" ) & NoAxes()

umap.all.before.drop <- ggarrange(umap1, umap2, umap3, umap4, ncol = 2, nrow=2)
umap.all.before.drop
ggsave("./images/umap_qc_rc.scr.pdf", scale = 2)

DimPlot(combined.scr, reduction = "umap", label = F, cols = custom_colors$c14 ) & NoAxes()
ggsave("./images/dimplot.pdf", width = 6, height = 4.5, scale = 1)

combined.scr@active.assay <- "RNA"
combined.scr <- NormalizeData(combined.scr)
combined.scr <- FindVariableFeatures(combined.scr, selection.method = "vst", nfeatures = 5000)
combined.scr <- ScaleData(combined.scr, features = , vars.to.regress = "batch")
combined.scr <- SCTransform(combined.scr, vst.flavor="v2", vars.to.regress = "batch") %>% PrepSCTFindMarkers(W)

combined.scr@active.assay <- "RNA"
FeaturePlot(combined.scr, c("Snap25", "Gad1", "Slc17a6", "Slc1a3", "Agt", "Mbp", "Plp1", "Olig2", "Pdgfra" ,"P2ry12", "Cx3cr1", "Rax", "Col23a1", "Flt1", "Slc47a1"), 
            min.cutoff = "q10" , max.cutoff = "q99", cols = c('#dfe4ea','#6F1E51') ) & NoLegend() & NoAxes()
ggsave("./images/featureplot_rc.pdf", scale = 2)


all.markers <- FindAllMarkers(combined.scr, assay = "RNA", slot= "scale.data", logfc.threshold = 0.25, min.pct = 0.5, only.pos = T)
top50 <- (all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_diff))

DoHeatmap(subset(combined.scr, downsample = 20), features= as.character(top50$gene), group.colors = custom_colors$c14, size = 3, disp.min = -1, disp.max = 2) + 
  scale_fill_gradient2(low = "lightsteelblue1", mid= "white", high= "brown" , midpoint =) + theme_classic()
ggsave("./images/DEG.pdf", scale=1, height = 5)

saveRDS(combined.scr, "combined.scr.RDS")

#DEG
Idents(combined.scr) <- "celltype"
a <-combined.scr@meta.data$group
b <- combined.scr@meta.data$celltype
c <- paste(b,a,sep = "_")
combined.scr@meta.data$deg.compare <- as.factor(c)
table <- data.frame("con" = paste(levels(combined.scr), "Con", sep = "_"), "tx" = paste(levels(combined.scr), "Tx", sep="_") )
Idents(combined.scr) <- "deg.compare"

# DE ---------------------------------------------------------------------


deg.per.cluster.wilcox <- list()
for(i in 1:nrow(table)) {
  deg.per.cluster.wilcox[[i]] <-  FindMarkers(combined.scr, ident.2 = table[i,1], ident.1 = table[i,2], test.use = , 
                                              pseudocount.use = 1, slot="scale.data", assay = "RNA")}
View(deg.per.cluster.wilcox[[7]])

deg.per.cluster.wilcox.sig.up <- lapply(deg.per.cluster.wilcox, function(x) {
  return(x %>% filter(p_val_adj < 0.01) %>% filter(avg_diff > 0.5))})
deg.per.cluster.wilcox.sig.down <- lapply(deg.per.cluster.wilcox, function(x) {
  return(x %>% filter(p_val_adj < 0.01) %>% filter(avg_diff < -0.5))})

a <- data.frame(up = unlist(lapply(deg.per.cluster.wilcox.sig.up, nrow)), 
                down = -1 * unlist(lapply(deg.per.cluster.wilcox.sig.down, nrow)),
                celltype = as.character(levels(combined.scr$celltype)))

a <- unlist(lapply(deg.per.cluster.wilcox.sig.up, nrow))
b <- -1*unlist(lapply(deg.per.cluster.wilcox.sig.down, nrow))
c <- data.frame(count = c(a,b),
                cell = rep(as.character(levels(combined.scr$celltype)),2),
                dir = rep(c('up','down'), each = 8))


ggplot(c, aes(x=factor(cell, levels = cell[1:8]), y=count, fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + scale_x_discrete(name = "Cluster") + scale_fill_npg()
ggsave("./images/deg_number.pdf", width = 5.3)



VlnPlot(subset(combined.scr, seurat_clusters == c(3,15)), "Nfkb1", assay = "RNA", split.by = "group", cols = custom_colors$c4)


temp <- subset(combined.scr, seurat_clusters == c(4,14,15,20,25))
temp2 <- as.data.frame(temp@assays$RNA@data)
ggplot(data = t(as.data.frame(subset(combined.scr, seurat_clusters == c(4,14,15,20,25))@assays$RNA@data)), aes(x = temp@meta.data$deg.compare, y = Nfkb1 )) + geom_boxplot()


# external database -------------------------------------------------------

lopes.up <- read.csv("Lopes 2022.csv", row.names = ) %>% filter(logFC > 0.01) %>% filter(adj.P.Val <= 0.01)
lopes.down <- read.csv("Lopes 2022.csv", row.names = ) %>% filter(logFC < -0.01) %>% filter(adj.P.Val <= 0.01)
galatro.up <- read.csv("Galatro 2017.csv", row.names = )[,-3] %>% filter(logFC > 0.01) %>% filter(P.Value <= 0.00001)
galatro.down <- read.csv("Galatro 2017.csv", row.names = )[,-3] %>% filter(logFC < -0.01) %>% filter(P.Value <= 0.00001)
chen.mg.up <- read.csv('Chen 2022 MG.csv') %>% filter(avg_log2FC > 0.5) %>% filter(p_val_adj <= 0.001)
chen.mg.down <- read.csv('Chen 2022 MG.csv') %>% filter(avg_log2FC < -0.5) %>% filter(p_val_adj <= 0.001)
ximerakis.up <- read.csv("Ximerakis 2019.csv", row.names = ) %>% filter(logFC_Young_to_Old > 0.01) %>% filter(padj <= 0.001)
ximerakis.down <- read.csv("Ximerakis 2019.csv", row.names = ) %>% filter(logFC_Young_to_Old < -0.01) %>% filter(padj <= 0.001)

chen.tc.up <- read.csv('Chen 2022 TC.csv') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj <= 0.05)
chen.tc.down <- read.csv('Chen 2022 TC.csv') %>% filter(avg_log2FC < 0) %>% filter(p_val_adj <= 0.05)
chen.od.up <- read.csv('Chen 2022 OD.csv') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj <= 0.05)
chen.od.down <- read.csv('Chen 2022 OD.csv') %>% filter(avg_log2FC < 0) %>% filter(p_val_adj <= 0.05)
chen.opc.up <- read.csv('Chen 2022 OPC.csv') %>% filter(avg_log2FC > 0) %>% filter(p_val_adj <= 0.05)
chen.opc.down <- read.csv('Chen 2022 OPC.csv') %>% filter(avg_log2FC < 0) %>% filter(p_val_adj <= 0.05)


lopes.up2 <- intersect(rownames(combined.scr), stringr::str_to_title(lopes.up$symbol))
lopes.down2 <- intersect(rownames(combined.scr), stringr::str_to_title(lopes.down$symbol))
galatro.up2 <- intersect(rownames(combined.scr), stringr::str_to_title(galatro.up$external_gene_name))
galatro.down2 <- intersect(rownames(combined.scr), stringr::str_to_title(galatro.down$external_gene_name))
chen.mg.up2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.mg.up$gene))
chen.mg.down2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.mg.down$gene))
ximerakis.up2 <- intersect(rownames(combined.scr), stringr::str_to_title(ximerakis.up$Gene))
ximerakis.down2 <- intersect(rownames(combined.scr), stringr::str_to_title(ximerakis.down$Gene))

chen.tc.up2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.tc.up$gene))
chen.tc.down2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.tc.down$gene))
chen.od.up2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.od.up$gene))
chen.od.down2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.od.down$gene))
chen.opc.up2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.opc.up$gene))
chen.opc.down2 <- intersect(rownames(combined.scr), stringr::str_to_title(chen.opc.down$gene))


test <- c("lopes.up2", "lopes.down2", "galatro.up2", "galatro.down2", "chen.mg.up2", "chen.mg.down2", "ximerakis.up2", "ximerakis.down2")

for(i in 1:length(test)) {
  combined.scr <- AddModuleScore(combined.scr,
                                 features = list(get(test[i])),
                                 name= test[i])
}

VlnPlot(combined.scr, paste0(test,"1"), group.by = "group", idents = "Mg", pt.size = 0, cols = custom_colors$c4[3:4], ncol = 4) & stat_compare_means(label = "p.signif")
ggsave("./images/comparison.pdf", width = 8, height = 4)

module.stats <- data.frame(deg.compare = combined.scr$deg.compare,
                           celltype = combined.scr$celltype,
                           lopes.up = combined.scr$lopes.up21,
                           lopes.down = combined.scr$lopes.down21,
                           galatro.up = combined.scr$galatro.up21,
                           galatro.down= combined.scr$galatro.down21,
                           chen.mg.up= combined.scr$chen.mg.up21,
                           chen.mg.down= combined.scr$chen.mg.down21,
                           ximerakis.up= combined.scr$ximerakis.up21,
                           ximerakis.down= combined.scr$ximerakis.down21)

module.stats %<>% group_by(deg.compare)

wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(lopes.up)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(lopes.up)) %>% .[,2])
wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(lopes.down)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(lopes.down)) %>% .[,2])

wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(galatro.up)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(galatro.up)) %>% .[,2])
wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(galatro.down)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(galatro.down)) %>% .[,2])

wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(chen.mg.up)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(chen.mg.up)) %>% .[,2])
wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(chen.mg.down)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(chen.mg.down)) %>% .[,2])

wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(ximerakis.up)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(ximerakis.up)) %>% .[,2])
wilcox.test(as.data.frame(module.stats %>% filter(deg.compare == "Mg_Con") %>% select(ximerakis.down)) %>% .[,2], 
            as.data.frame(module.stats %>% filter(deg.compare == "Mg_Tx") %>% select(ximerakis.down)) %>% .[,2])



# neuron subset -----------------------------------------------------------

##gaba
gaba <- subset(combined.scr, subset = celltype == "GABA")
gaba <- NormalizeData(gaba)
gaba <- ScaleData(gaba, features = rownames(gaba), vars.to.regress = "batch")
gaba <- SCTransform(gaba, vst.flavor="v2", vars.to.regress = "batch") %>% PrepSCTFindMarkers()
gaba <- RunPCA(gaba, npcs = 50, verbose = FALSE) %>% FindNeighbors(combined, reduction = "pca", dims = 1:50)
gaba <- RunUMAP(gaba, reduction = "pca", dims = 1:50, n.neighbors = 10)
labels <- forcats::fct_recode(gaba$seurat_clusters, "GABA-1" = as.character(1), "GABA-2" = as.character(5), "GABA-3" = as.character(11),
                              "GABA-4" = as.character(12),"GABA-5" = as.character(18),"GABA-6" = as.character(19),"GABA-7" = as.character(21), 
                              "GABA-8" = as.character(23),"GABA-9" = as.character(27),"GABA-10" = as.character(29))
labels <- forcats::fct_drop(labels)
gaba <- AddMetaData(gaba, labels, col.name = "celltype")
Idents(gaba) <- "celltype"
head(Idents(gaba), 5)
DimPlot(gaba, cols = custom_colors$c15[c(2,4:13)]) + NoAxes()
ggsave("./images/gaba_dimplot.pdf")

gaba@active.assay <- "RNA"
gaba.all.marker <- FindAllMarkers(gaba, assay = "RNA", slot = "scale.data", only.pos = T)
View(gaba.all.marker %>% group_by(cluster) %>% top_n(20, wt= avg_diff)) 

features <- gaba.all.marker %>% group_by(cluster) %>% top_n(4, wt= avg_diff) %>% select(gene)
DotPlot(gaba, features = rev(features$gene) , cols = c('#dfe4ea', custom_colors$c14[1]), 
        col.min = 0.5, col.max = 2, dot.min = 0.25, dot.scale = 2, assay = "RNA" ) + coord_flip() + theme_bw() + scale_y_discrete(guide = guide_axis(angle= -90))
ggsave("./images/potential suppl/gaba marker.pdf", width = 5, height = 7)

a <- gaba@meta.data$group
b <- gaba@meta.data$celltype
c <- paste(b,a,sep = "_")
gaba@meta.data$deg.compare <- as.factor(c)
Idents(gaba) <- "celltype"
table <- data.frame("con" = paste(levels(gaba), "Con", sep= "_"), "tx" = paste(levels(gaba), "Tx", sep="_") )
Idents(gaba) <- "deg.compare"


deg.gaba <- list()
for(i in 1:nrow(table)) {
  deg.gaba[[i]] <-  FindMarkers(gaba, ident.2 = table[i,1], ident.1 = table[i,2], test.use = , 
                                pseudocount.use = 1, slot="scale.data", assay = "RNA")}
View(deg.gaba[[10]])

deg.gaba.sig.up <- lapply(deg.gaba, function(x) {
  return(x %>% filter(p_val_adj < 0.01) %>% filter(avg_diff > 0.5))})
deg.gaba.sig.down <- lapply(deg.gaba, function(x) {
  return(x %>% filter(p_val_adj < 0.01) %>% filter(avg_diff < -0.5))})

a <- data.frame(up = unlist(lapply(deg.gaba.sig.up, nrow)), 
                down = -1 * unlist(lapply(deg.gaba.sig.down, nrow)),
                celltype = as.character(levels(gaba$seurat_clusters)))

a <- unlist(lapply(deg.gaba.sig.up, nrow))
b <- -1*unlist(lapply(deg.gaba.sig.down, nrow))
c <- data.frame(count = c(a,b),
                cell = rep(as.character(levels(gaba$celltype)),2),
                dir = rep(c('up','down'), each = 10))


ggplot(c, aes(x=factor(cell, levels = cell[1:10]), y=count, fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + scale_x_discrete(name = "Cluster") + scale_fill_npg()
ggsave("./images/deg_number_gaba.pdf", width = 9)


ego.gaba.up <- list()
ego.gaba.down <- list()
for(i in 6){
  x <- enrichGO(gene          = rownames(deg.gaba.sig.up[[i]]),
                universe      = rownames(gaba),
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = F,
                minGSSize = 10)
  
  x %<>%
    mutate(GO_levels = sapply(ID, function(ID) {
      which(
        unlist(lapply(ont_list, function(s) ID %in% s))
      )
    }),
    
    first_GO_level = unlist(lapply(GO_levels, min))
    )
  
  ego.gaba.up[[i]] <- clusterProfiler::simplify(x, cutoff = 0.5, by = "first_GO_level")
  
  x <- enrichGO(gene          = rownames(deg.gaba.sig.down[[i]]),
                universe      = rownames(gaba),
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = F,
                minGSSize = 10)
  
  x %<>%
    mutate(GO_levels = sapply(ID, function(ID) {
      which(
        unlist(lapply(ont_list, function(s) ID %in% s))
      )
    }),
    
    first_GO_level = unlist(lapply(GO_levels, min))
    )
  
  ego.gaba.down[[i]] <- clusterProfiler::simplify(x, cutoff = 0.5, by = "first_GO_level") 
}


##glut
glut <- subset(combined.scr, subset = celltype == "Glut")
glut <- NormalizeData(glut)
glut <- ScaleData(glut, features = rownames(glut), vars.to.regress = "batch")
glut <- SCTransform(glut, vst.flavor="v2", vars.to.regress = "batch") %>% PrepSCTFindMarkers()
glut <- RunPCA(glut, npcs = 50, verbose = FALSE) %>% FindNeighbors(combined, reduction = "pca", dims = 1:50)
glut <- RunUMAP(glut, reduction = "pca", dims = 1:50, n.neighbors = 10)
labels <- forcats::fct_recode(glut$seurat_clusters, "Glut-1" = as.character(2), "Glut-2" = as.character(6), "Glut-3" = as.character(7),
                              "Glut-4" = as.character(8),"Glut-5" = as.character(10),"Glut-6" = as.character(13),"Glut-7" = as.character(16), 
                              "Glut-8" = as.character(17),"Glut-9" = as.character(22),"Glut-10" = as.character(24), "Glut-11" = as.character(26), 
                              "Glut-12" = as.character(28), "Glut-13" = as.character(30))
labels <- forcats::fct_drop(labels)
glut <- AddMetaData(glut, labels, col.name = "celltype")
Idents(glut) <- "celltype"
head(Idents(glut), 5)
DimPlot(glut, cols = custom_colors$c15[c(15,13:1)]) + NoAxes()
ggsave("./images/glut_dimplot.pdf")
glut@active.assay <- "RNA"
glut.all.marker <- FindAllMarkers(glut, assay = "RNA", slot = "scale.data", only.pos = T)
View(glut.all.marker %>% group_by(cluster) %>% top_n(20, wt= avg_diff)) 



features <- glut.all.marker %>% group_by(cluster) %>% top_n(4, wt= avg_diff) %>% select(gene)
DotPlot(glut, features = rev(unique(features$gene)) , cols = c('#dfe4ea', custom_colors$c14[2]), 
        col.min = 0.5, col.max = 2, dot.min = 0.25, dot.scale = 2, assay = "RNA" ) + coord_flip() + theme_bw() + scale_y_discrete(guide = guide_axis(angle= -90))
ggsave("./images/potential suppl/glut marker.pdf", width = 5, height = 7)



a <- glut@meta.data$group
b <- glut@meta.data$celltype
c <- paste(b,a,sep = "_")
glut@meta.data$deg.compare <- as.factor(c)
Idents(glut) <- "celltype"
table <- data.frame("con" = paste(levels(glut), "Con", sep= "_"), "tx" = paste(levels(glut), "Tx", sep="_") )
Idents(glut) <- "deg.compare"


deg.glut <- list()
for(i in 1:nrow(table)) {
  deg.glut[[i]] <-  FindMarkers(glut, ident.2 = table[i,1], ident.1 = table[i,2], test.use = , 
                                pseudocount.use = 1, slot="scale.data", assay = "RNA")}
View(deg.glut[[13]])

deg.glut.sig.up <- lapply(deg.glut, function(x) {
  return(x %>% filter(p_val_adj < 0.01) %>% filter(avg_diff > 0.5))})
deg.glut.sig.down <- lapply(deg.glut, function(x) {
  return(x %>% filter(p_val_adj < 0.01) %>% filter(avg_diff < -0.5))})

a <- data.frame(up = unlist(lapply(deg.glut.sig.up, nrow)), 
                down = -1 * unlist(lapply(deg.glut.sig.down, nrow)),
                celltype = as.character(levels(glut$celltype)))

a <- unlist(lapply(deg.glut.sig.up, nrow))
b <- -1*unlist(lapply(deg.glut.sig.down, nrow))
c <- data.frame(count = c(a,b),
                cell = rep(as.character(levels(glut$celltype)),2),
                dir = rep(c('up','down'), each = 13))


ggplot(c, aes(x=factor(cell, levels = cell[1:13]), y=count, fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + scale_x_discrete(name = "Cluster") + scale_fill_npg()
ggsave("./images/deg_number_glut.pdf", width = 9)



ego.glut.up <- list()
ego.glut.down <- list()
for(i in 1:11){
  x <- enrichGO(gene          = rownames(deg.glut.sig.up[[i]]),
                universe      = rownames(glut),
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = F,
                minGSSize = 10,
  )
  
  x %<>%
    mutate(GO_levels = sapply(ID, function(ID) {
      which(
        unlist(lapply(ont_list, function(s) ID %in% s))
      )
    }),
    
    first_GO_level = unlist(lapply(GO_levels, min))
    )
  
  ego.glut.up[[i]] <- clusterProfiler::simplify(x, cutoff = 0.5, by = "first_GO_level")
  
  x <- enrichGO(gene          = rownames(deg.glut.sig.down[[i]]),
                universe      = rownames(glut),
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = F,
                minGSSize = 10)
  
  x %<>%
    mutate(GO_levels = sapply(ID, function(ID) {
      which(
        unlist(lapply(ont_list, function(s) ID %in% s))
      )
    }),
    
    first_GO_level = unlist(lapply(GO_levels, min))
    )
  
  ego.glut.down[[i]] <- clusterProfiler::simplify(x, cutoff = 0.5, by = "first_GO_level") 
}

ggarrange(plotlist = list(VlnPlot(gaba, "Lgr4", idents = "GABA-4", split.by = "group", pt.size = 0) + stat_compare_means(label = "p.signif"),
                          VlnPlot(gaba, "Cd63", idents = "GABA-4", split.by = "group", pt.size = 0) + stat_compare_means(label = "p.signif"),
                          VlnPlot(gaba, "Nr3c1", idents = "GABA-6", split.by = "group", pt.size = 0) + stat_compare_means(label = "p.signif"),
                          VlnPlot(glut, "Nkx2-1", idents = "Glut-5", split.by = "group", pt.size = 0) + stat_compare_means(label = "p.signif"),
                          VlnPlot(glut, "Vgf", idents = "Glut-9", split.by = "group", pt.size = 0) + stat_compare_means(label = "p.signif"),
                          VlnPlot(glut, "Pmch", idents = "Glut-4", split.by = "group",pt.size = 0) + stat_compare_means(label = "p.signif")), ncol = 3, nrow = 2)
ggsave("./images/potential suppl/DEG neuron results.pdf", height =4, width = 8)


# CellChat ----------------------------------------------------------------

#CellChat
library(CellChat)
Idents(combined.scr) <- "seurat_clusters"
Idents(combined.scr)



cc.1 <- GetAssayData(subset(combined.scr, subset = group == 'Con'), assay = "RNA", slot = "data")
labels <- subset(combined.scr, subset = group == 'Con')$celltype
meta <- data.frame(cluster = labels, row.names = names(labels))
cc.1 <- createCellChat(object = cc.1, meta = meta, group.by = "cluster")


cc.2 <- GetAssayData(subset(combined.scr, subset = group == 'Tx'), assay = "RNA", slot = "data")
labels <- subset(combined.scr, subset = group == 'Tx')$celltype
meta <- data.frame(cluster = labels, row.names = names(labels))
cc.2 <- createCellChat(object = cc.2, meta = meta, group.by = "cluster")

CellChatDB.use <- CellChatDB.mouse
cc.1@DB <- CellChatDB.use
cc.2@DB <- CellChatDB.use

cc.1 <- subsetData(cc.1)
cc.2 <- subsetData(cc.2)

future::plan("multisession", workers = 11)
cc.1 <- identifyOverExpressedGenes(cc.1)
cc.1 <- identifyOverExpressedInteractions(cc.1)
cc.2 <- identifyOverExpressedGenes(cc.2)
cc.2 <- identifyOverExpressedInteractions(cc.2)
cc.1 <- projectData(cc.1, PPI.mouse)
cc.2 <- projectData(cc.2, PPI.mouse)
cc.1 <- computeCommunProb(cc.1, population.size = F)
cc.2 <- computeCommunProb(cc.2, population.size = F)
cc.1 <- computeCommunProbPathway(cc.1)
cc.2 <- computeCommunProbPathway(cc.2)
cc.1 <- aggregateNet(cc.1)
cc.2 <- aggregateNet(cc.2)
cc.1 <- netAnalysis_computeCentrality(cc.1, slot.name = "netP") 
cc.2 <- netAnalysis_computeCentrality(cc.2, slot.name = "netP")

cc.1@netP

pathways.show="VEGF"
netVisual_heatmap(cc.1, signaling = pathways.show, color.heatmap = "Reds")
netVisual_heatmap(cc.2, signaling = pathways.show, color.heatmap = "Reds")

par(mfrow=c(1,2))
netVisual_aggregate(cc.1, signaling = pathways.show, layout = "chord")
netVisual_aggregate(cc.2, signaling = pathways.show, layout = "chord")


vertex.receiver = seq(1,8) # a numeric vector. 
netVisual_aggregate(cc.1, signaling = pathways.show,  vertex.receiver = vertex.receiver)
netVisual_aggregate(cc.2, signaling = pathways.show,  vertex.receiver = vertex.receiver)


groupSize <- as.numeric(table(cc.2@idents))
netVisual_circle(cc.2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cc.2@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}




#merge cellchat

object.list <- list(Con = cc.1, Tx = cc.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

pathways.show <- gg1$data %>% filter(pvalues < .05) %>% 
  filter(contribution.relative.1 > 1.05 | contribution.relative.1 < 0.95) %>% rownames() %>% .[1:29]


color <- custom_colors$c14[1:8]
names(color) <- levels(cc.1@idents)

pathway.union <- c("VEGF","PMCH", "ncWNT", "GDF")
i=1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], 
                                        width = 4, height = 2, color.use = color, color.heatmap = "YlGnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], 
                                        width = 4, height = 2, color.use = color, color.heatmap = "YlGnBu")
ht1+ht2

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
i=1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], 
                                        width = 4, height = 16, color.use = color, color.heatmap = "YlGnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], 
                                        width = 4, height = 16, color.use = color, color.heatmap = "YlGnBu")
ht1+ht2


netVisual_bubble(cellchat, sources.use = 1:8, targets.use = 8,  comparison = c(1, 2), angle.x = 45)

pathways.show <- gg1$data %>% filter(pvalues < .001) %>% 
  filter(contribution.relative.1 > 1.25 | contribution.relative.1 < 0.75) %>% rownames() %>% .[1:29]

for(i in 3:length(pathways.show)){
  netVisual_aggregate(cc.1, signaling = pathways.show[i],  vertex.receiver = vertex.receiver)
  netVisual_aggregate(cc.2, signaling = pathways.show[i],  vertex.receiver = vertex.receiver)
}

pathways.show

par(mfrow=c(1,2))
a <- 'ncWNT'
netVisual_aggregate(cc.1, signaling = a,  vertex.receiver = vertex.receiver, layout = "chord", color.use = custom_colors$c14, vertex.label.cex = 1.2) 
netVisual_aggregate(cc.2, signaling = a,  vertex.receiver = vertex.receiver, layout = "chord", color.use = custom_colors$c14, vertex.label.cex = 1.2)

