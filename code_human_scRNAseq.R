library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(scran)
library(viridis)
library(ggforce)
library(gghalves)
library(ggridges)
library(scDblFinder)
library(ggpubr)
library(SingleCellExperiment)
library("leidenAlg")


custom_colors$c15 <- c('#ED4C67', '#ff6348', '#F79F1F','#f9ca24',  '#A3CB38', 
                       '#20bf6b', '#4b7bec', '#2d98da', '#12CBC4', '#1289A7', 
                       '#778ca3', '#82589F', '#D980FA', '#f78fb3', '#f8a5c2')
custom_colors$c4 <- c('#1289A7', '#ffda79', '#a4b0be', '#ba4b7f')

# pre-processing ----------------------------------------------------------



data <- Read10X(data.dir = "/home/alan/Hypo_scRNAseq/filtered_feature_bc_matrix")
at_least_one <- apply(data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

cells <- tibble(
  cell = colnames(data),
  nCount = colSums(data),
  nFeature = colSums(data != 0)
)

cells$percent_MT <- Matrix::colSums(data[grep("^MT-", rownames(data)),]) / Matrix::colSums(data)


#doublet removal

sce <- SingleCellExperiment(assays = list(counts = data))
set.seed(100)
sce <- scDblFinder(sce)
cells$multiplet_class <- colData(sce)$scDblFinder.class
table(cells$multiplet_class)

p1 <- ggplot(cells, aes(x = multiplet_class, y = nCount, fill = multiplet_class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$c4) +
  scale_x_discrete(limits = (levels(cells$multiplet_class))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Transcripts', subtitle = 'Counts') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )

p2 <- ggplot(cells, aes(x = multiplet_class, y = nFeature, fill = multiplet_class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$c4) +
  scale_x_discrete(limits = (levels(cells$multiplet_class))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Genes', subtitle = 'Counts') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) 


ggarrange(p1 , p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

ggsave("doubletQC.jpg")


cells.unfiltered <- cells

cells <- filter(cells, multiplet_class == 'singlet')

#prefiltering
median_nCount <- median(cells$nCount)
mad_nCount <- mad(cells$nCount)
median_nFeature <- median(cells$nFeature)
mad_nFeature <- mad(cells$nFeature)
median_percent_MT <- median(cells$percent_MT)
mad_percent_MT <- mad(cells$percent_MT)

thresholds_nCount <- c(2000, 40000)
thresholds_nFeature <- c(200, 10000)
thresholds_percent_MT <- c(0, 0.05)

p3 <- ggplot(cells, aes(nCount, nFeature, color = percent_MT)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  geom_vline(xintercept = thresholds_nCount, color = 'red') +
  scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
  theme_bw() +
  scale_color_viridis(
    name = 'Percent MT\ntranscripts',
    limits = c(0,0.25),
    labels = scales::percent,
    guide = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')
  )

p3
ggsave("cellexclude.jpg")

cells_filtered <- cells %>%
  dplyr::filter(
    nCount >= thresholds_nCount[1],
    nCount <= thresholds_nCount[2],
    nFeature >= thresholds_nFeature[1],
    nFeature <= thresholds_nFeature[2],
    percent_MT >= thresholds_percent_MT[1],
    percent_MT <= thresholds_percent_MT[2]
  )

cells_to_keep <- cells_filtered$cell
length(cells_to_keep)

hypo <- CreateSeuratObject(counts = data[,cells_to_keep], project = "hypo10k", min.cells = 10, min.features = 0)
plot <- p1+p2+p3
plot
ggsave("QCplot.pdf", width = 12)

hypo <- SCTransform(hypo)
hypo <- FindVariableFeatures(hypo)
hypo <- RunPCA(hypo, npcs = 50, verbose = FALSE)
hypo <- FindNeighbors(hypo, reduction = "pca", dims = 1:50)
hypo <- RunUMAP(hypo, reduction = "pca", dims = 1:50, n.neighbors = 30, min.dist = 0.3)
hypo <- FindClusters(hypo, resolution = 0.55, algorithm = 4)
head(Idents(hypo), 5)

umap.before.drop <- DimPlot(hypo, reduction = "umap", label = T, cols = custom_colors$c15) & NoAxes() & NoLegend()
umap.before.drop
ggsave("umap.before.drop.pdf")

temp_labels <- hypo@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()

p1 <- ggplot() +
  geom_boxplot(
    data = hypo@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
    show.legend = FALSE,
  ) +
  geom_text(
    data = temp_labels,
    aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',',)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_color_manual(values = custom_colors$c15) +
  scale_fill_manual(values = custom_colors$c15) +
  scale_y_continuous(name = 'Number of transcripts', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank()
  )

p2 <- ggplot() +
  geom_boxplot(
    data = hypo@meta.data, aes(seurat_clusters, nFeature_RNA, fill = seurat_clusters),
    show.legend = FALSE
  ) +
  geom_text(
    data = temp_labels,
    aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',')), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_color_manual(values = custom_colors$c15) +
  scale_fill_manual(values = custom_colors$c15) +
  scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma, expand = c(0.08, 0)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank()
  )

cluster.drop <- ggarrange(p1,p2, nrow=2)
cluster.drop
ggsave("cluster.drop.pdf", width = 7)

#subset the drop-out cells
hypo <- subset(hypo, idents = 12, invert = T)
table(hypo@active.ident)
cellname <- c("Int1","Neuron1","Ast3","Ast1","Ast2","Int2","Neuron2","Tan","IPC","RG","APC","NSC","PGG")
names(cellname) <- levels(hypo)
hypo <- RenameIdents(hypo, cellname)

DimPlot(hypo, reduction = "umap", label = T, cols = custom_colors$c15) & NoAxes() & NoLegend()
ggsave("dimplot.pdf")


lookup <- data.frame(number = as.factor(levels(hypo$seurat_clusters))[-12], cell = as.factor(levels(hypo@active.ident)))
lookup
hypo[["cell_type"]] <- Idents(hypo)


# cluster DEG -------------------------------------------------------------

all.markers <- FindAllMarkers(hypo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(print(all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)))
hypo.markers <- (print(all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)))
write.csv(hypo.markers, file = "./hypo.markers.csv")



VlnPlot(hypo, features = c("SNAP25","SYP", 'AQP4', 'GFAP', 'SLC1A3', 'GPC3', 'RAX', 'AGT', 'ZIC1', 'TOP2A', 'COL1A1', 'BGN'), pt.size = 0, ncol = 15, log = F, cols = custom_colors$c15)  & 
  NoAxes() & coord_flip() &   
  theme(title= element_text(angle = 90, face="plain", size=8), panel.spacing = unit(0,'lines'), axis.text.x = element_text(size=7))

ggsave("vlnplot.pdf")


FeaturePlot(hypo, features = c("SNAP25","SYP",'STMN2','STMN4', 'SLC17A6','SLC32A1', 'GFAP', 'SLC1A3', 'CD44', 'CRYM', 'RAX', 'TOP2A', 'MKI67', 'COL1A1', 'BGN'), pt.size = 0, ncol = 3, max.cutoff = "q90", cols = c('#dfe4ea','#6F1E51')  ) & NoAxes() & NoLegend()
ggsave("featureplot.pdf", height = 7)


# monocle3 ----------------------------------------------------------------


library(remotes)
library(monocle3)
library(SeuratWrappers)
library(circlize)
hypo@active.assay <- "RNA"
hypo.cds <- as.cell_data_set(hypo)
hypo.cds <- estimate_size_factors(hypo.cds) #there is a bug in SeuratWrapper, this step is to prevent further error
hypo.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(hypo[["RNA"]]) #there is a bug in SeuratWrapper, this step is to prevent further error
hypo.cds <- cluster_cells(cds = hypo.cds, reduction_method = "UMAP", cluster_method = "leiden", k=10)
hypo.cds <- learn_graph(hypo.cds, use_partition = F)
hypo.cds <- order_cells(hypo.cds, reduction_method = "UMAP")

plot_cells(
  cds = hypo.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = F
)

plot_cells(
  cds = hypo.cds,
  genes = "SNAP25", show_trajectory_graph = F)

plot_cells(
  cds = hypo.cds, color_cells_by = "ident" , show_trajectory_graph = F)

#add pseudotime metadata back to seurat
hypo<- AddMetaData(
  object = hypo,
  metadata = hypo.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

FeaturePlot(hypo, c("pseudotime"), pt.size = 0.01) & scale_color_viridis_c() & NoAxes()
ggsave("pseudotype_hypo.pdf")

neuron.cluster.pathway <- as.factor(c('NSC','Int1','Int2','Neuron1','Neuron2'))
tan.cluster.pathway <- as.factor(c('NSC','Int1', 'APC', 'Ast2', "Tan"))


#monocle deg followed by pseudotimes regualted gene
hypo.cds.res <- graph_test(hypo.cds, neighbor_graph="principal_graph", cores=10)
pr_deg_ids <- hypo.cds.res[order(hypo.cds.res$q_value),]
plot_cells(hypo.cds, genes = rownames(pr_deg_ids[1:5,]), show_trajectory_graph = F)

deg.pt = rownames(pr_deg_ids[1:5,])
deg.pt = c("NKX2-1","GFAP","PMCH","GAD1","DCX")

#removing pseudotime outlier
z <- hypo.cds[rowData(hypo.cds)$gene_short_name %in% c("SNAP25", "SYP", "ASCL1", "TUBB3", "CDK1", "MAPT", "HES5", "RBFOX3", "NEUROD1"),
              colData(hypo.cds)$ident %in% neuron.cluster.pathway & z@principal_graph_aux@listData$UMAP$pseudotime <30]

plot_genes_in_pseudotime(z, min_expr=0.5, color_cells_by="ident", ncol = 2) + scale_color_manual(values = custom_colors$c15[c(1,2,6,7,12)])

ggsave("genes_in_pt_neuron.pdf")

z <- hypo.cds[rowData(hypo.cds)$gene_short_name %in% c("GFAP","APOE","GPC3","NFIB","SPARCL1", "RAX","FGF10","COL25A1"),
              colData(hypo.cds)$ident %in% tan.cluster.pathway & z@principal_graph_aux@listData$UMAP$pseudotime <25]
plot_genes_in_pseudotime(z, min_expr=0.5, color_cells_by="ident", ncol = 2) + scale_color_manual(values = custom_colors$c15[c(1,5,8,11,12)]) 


#monocle heatmap


pt.matrix <- as.matrix(z@assays@data$counts[,order(pseudotime(z))])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rowData(z)[,1];
ht <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 12),
  row_km = 3,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
)
print(ht)

# voxhunt -----------------------------------------------------------------


devtools::install_github('quadbiolab/voxhunt')
library(voxhunt)
load_aba_data('/disk1/voxhunt_data')
genes_use <- variable_genes(stage='P14', nfeatures=150)$gene
vox_map <- voxel_map(hypo, genes_use=genes_use, stage = 'P14')
p1 <-plot_map(vox_map) & no_legend()
p1 + voxhunt::plot_annotation('P14', show_legend=F) 
ggsave("voxhunt_p14.pdf")

plot_map(vox_map, view='slice', slices=seq(from=20, to= 54, by=4)) & no_legend()
ggsave("voxhunt_p14_slice.pdf")

plot_map(vox_map, view='slice', slices=seq(from=20, to= 54, by=4)) 



# reclustering neuron -----------------------------------------------------


neuron <- subset(hypo, idents = c("Neuron1","Neuron2"))
neuron@active.assay <- "SCT"
neuron <- RunPCA(neuron, npcs = 30)
neuron <- FindNeighbors(neuron, dims = 1:30)
neuron <- FindClusters(neuron, resolution = 2, algorithm = 4)
head(Idents(neuron), 5)
neuron <- RunUMAP(neuron, dims = 1:30, min.dist = 0.3, n.neighbors = 10)
Idents(neuron) <- "seurat_clusters"
DimPlot(neuron, reduction = "umap", label = F, cols = rep(custom_colors$c15,2)) + NoAxes()
ggsave("./neuron2/neuron.dimplot.pdf", width = 6, height = 4)

neuron@active.assay <- "RNA"
neuron <- NormalizeData(neuron)
neuron <- ScaleData(neuron)
neuron <- FindVariableFeatures(neuron)
all.markers.neuron <- FindAllMarkers(neuron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.9, min.diff.pct = 0.25, test.use = )
View(all.markers.neuron %>% group_by(cluster)  %>% arrange(desc(avg_log2FC), .by_group = T))

FeaturePlot(neuron, features = c("SCGN", "PMCH", "PNOC", "TRH", "NR3C1", "FEZF2", "PITX2", "ONECUT2", "VGF"), max.cutoff = "q90", min.cutoff = "q1", cols = c('#dfe4ea','#6F1E51')  ) & NoAxes() & NoLegend()
ggsave("./neuron2/neuron.featureplot.pdf", width = 4, height = 4)

neuron <- BuildClusterTree(neuron)
PlotClusterTree(neuron)
DimPlot(neuron, label = T, cols = custom_colors$c15)

data.tree <- Tool(object = neuron, slot = "BuildClusterTree")
ape::plot.phylo(x = data.tree, direction = "downwards", edge.width = 3, use.edge.length = F, show.tip.label = F, edge.color = "#6D214F")
nodes <- unique(neuron@tools$BuildClusterTree$edge[,1])


# pyscenic ----------------------------------------------------------------

#subset protein coding genes for pyscenic
pc <- read.delim(file = "protein_coding.txt")
pcgenes <- intersect(pc[,2], rownames(hypo))

subset.matrix <- hypo@assays$RNA@data[pcgenes, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
hypo.pc <- CreateSeuratObject(subset.matrix, project = 'hypo_pc') # Create a new Seurat object with just the genes of interest
hypo.pc@meta.data$cell_type <- hypo@active.ident
Idents(hypo.pc) <- 'cell_type'
hypo.pc@active.ident


subset.matrix <- neuron@assays$RNA@data[pcgenes, ] 
neuron.pc <- CreateSeuratObject(subset.matrix, project = 'neuron_pc')
neuron.pc@meta.data$cell_type <- neuron@active.ident
Idents(neuron.pc) <- 'cell_type'
neuron.pc@active.ident

#exporting to pyscenic
library(remotes)
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
SaveLoom(hypo.pc)
SaveLoom(neuron.pc)

# importing from pyscenic
hypo.pyscenic <- hypo
regulonAUC <- read.csv(file = "/home/alan/data/pyscenic/auc_mtx.csv", row.names = 1)
colnames(regulonAUC) <- gsub("\\.\\.\\.", "\\(\\+\\)", colnames(regulonAUC))

table(apply(regulonAUC, 2, function(x) {sum(x>0)} >= 20)) #check whether the regulon has at least X cells with active regulons, too few cells cannot make treshold
a <- names(which(apply(regulonAUC, 2, function(x) {sum(x>0)}) >= 20))
regulonAUC <- regulonAUC[,a]
hypo.pyscenic <- AddMetaData(hypo, regulonAUC.hypo)

FeaturePlot(hypo.pyscenic, features=c("ISL1", "ISL1.reg"), min.cutoff = "q20", max.cutoff = "q99")
FeaturePlot(hypo.pyscenic, features=colnames(regulonAUC)[141:160]) & NoLegend() & NoAxes()
FeaturePlot(hypo.pyscenic, features=c("EGR3", "EGR3.reg"))
FeaturePlot(hypo.pyscenic, features=c("STAT1.reg", "BCL3.reg"), min.cutoff = "q20", max.cutoff = "q95") & NoLegend()


a <- hist(regulonAUC$ISL1..., breaks = 200, plot = F)
hist(regulonAUC$SF1..., breaks = 200, freq= T, plot = T, xlab = "T")
lines(density(regulonAUC$ISL1), lty=3, col="Red")

regulonAUC.hypo <- regulonAUC
regulonAUC <- regulonAUC.hypo

colnames(regulonAUC.hypo) <- gsub("\\.reg", "\\(\\+\\)", colnames(regulonAUC.hypo))

table(apply(regulonAUC.hypo, 2, function(x) {sum(x>0)} >= 10)) #check whether the regulon has at least X cells with active regulons, too few cells cannot make treshold

hypo.regulon <- CreateSeuratObject(t(regulonAUC.hypo))
hypo.regulon <- ScaleData(hypo.regulon)
table(colnames(hypo) == colnames(hypo.regulon))
hypo.regulon <- FindVariableFeatures(hypo.regulon)
hypo.regulon <- RunPCA(hypo.regulon)
Idents(hypo.regulon) <- Idents(hypo)
hypo.regulon <- RunUMAP(hypo.regulon, dims = 1:20, n.neighbors = 10, min.dist =0.01)
p1 <- DimPlot(hypo.regulon, reduction = "umap", label = T) + NoLegend()
p2 <- DimPlot(hypo, label = T, cols = custom_colors$c15)
p1+p2




all.markers.hypo.regulon <- FindAllMarkers(hypo.regulon, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.01, min.diff.pct = 0.05)
tip.marker.regulon <- print(all.markers.hypo.regulon %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% arrange(match(cluster, tip)))
tip.marker.regulon$gene


View(print(all.markers.hypo.regulon %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)))




#neuron pyscenic

neuron.pyscenic <- neuron
regulonAUC <- read.csv(file = "~/data/pyscenic2/neuron2/auc_mtx.csv", row.names = 1)
colnames(regulonAUC) <- gsub("\\.\\.\\.", "\\(\\+\\)", colnames(regulonAUC))

table(apply(regulonAUC, 2, function(x) {sum(x>0)} >= 10)) #check whether the regulon has at least X cells with active regulons, too few cells cannot make treshold
a <- names(which(apply(regulonAUC, 2, function(x) {sum(x>0)}) >= 20))
regulonAUC <- regulonAUC[,a]
neuron.pyscenic <- AddMetaData(neuron, regulonAUC)
regulonAUC.neuron <- regulonAUC

neuron.regulon <- CreateSeuratObject(t(regulonAUC))
neuron.regulon <- ScaleData(neuron.regulon)
table(colnames(neuron) == colnames(neuron.regulon))
neuron.regulon <- FindVariableFeatures(neuron.regulon)
neuron.regulon <- RunPCA(neuron.regulon)
Idents(neuron.regulon) <- Idents(neuron)
neuron.regulon <- RunUMAP(neuron.regulon, dims = 1:20, n.neighbors = 10, min.dist =0.01)
p1 <- DimPlot(neuron.regulon, reduction = "umap", label = T) + NoLegend()
p2 <- DimPlot(neuron, label = F, cols = custom_colors$c15)
p1+p2


neuron.regulon <- BuildClusterTree(neuron.regulon)
PlotClusterTree(neuron.regulon)
PlotClusterTree(neuron)
DimPlot(neuron, label = T)
tip <- c(8,12,5,13,6,11,2,4,10,3,1,7,9)


nodes <- unique(neuron.regulon@tools$BuildClusterTree$edge[,1])
tree_markers <- list()
for(i in 1:length(nodes)) {
  tree_markers[[i]] <- FindMarkers(
    hypo, ident.1 = "clustertree", ident.2 = nodes[i], only.pos = T, slot = "scale.data")
}

mako.col <- mako(10)
rocket.col <- rocket(10)


all.markers.neuron.regulon <- FindAllMarkers(neuron.regulon, only.pos = TRUE, min.pct = 0.1, slot = "scale.data")
View(print(all.markers.neuron.regulon %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff) %>% arrange(match(cluster, tip))))
tip.marker.regulon <- print(all.markers.neuron.regulon %>% group_by(cluster) %>% top_n(n = 2, wt = avg_diff) %>% arrange(match(cluster, tip)))
tip.marker.regulon$gene

neuron.regulon@active.ident <- factor(neuron.regulon@active.ident, levels = tip)
DotPlot(neuron.regulon, features = rev(unique(tip.marker.regulon$gene)), cols = c('#dfe4ea','#1e3799'), col.min = 0.2, col.max = 2, dot.min = 0, dot.scale =3 ) +
  scale_color_gradientn(colours = rev(rocket.col)[1:9]) & coord_flip() & NoLegend()



neuron <- BuildClusterTree(neuron)
PlotClusterTree(neuron)

p1 <- ColorDimSplit(neuron, node=15)
p2 <- DimPlot(neuron, label=T)
p1+p2

FeaturePlot(neuron, features = genes_use[1:16], pt.size = 0, ncol = 3, max.cutoff = "q90", cols = c('#dfe4ea','#6F1E51')  ) & NoAxes() & NoLegend()
all.markers.neuron %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


tip.marker <- print(all.markers.neuron %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) %>% arrange(match(cluster, tip)))
tip.marker$gene


Idents(neuron) <- factor(neuron@active.ident, levels = tip)
Idents(neuron.regulon) <- factor(neuron.regulon@active.ident, levels = tip)


FeaturePlot(neuron.regulon, features = unique(tip.marker.regulon$gene), pt.size = 0, ncol = 3, max.cutoff = "q90", cols = c('#dfe4ea','#6F1E51')  ) & NoAxes() & NoLegend()

VlnPlot(neuron, features = c('SLC17A6','SLC32A1','TH'), pt.size = 0, ncol = 1, log =T, cols = custom_colors$c15[tip])
ggsave("./neuron2/vln.pdf", width = 7, height = 4.3)

ape::plot.phylo(x = data.tree, direction = "downwards", edge.width = 3, use.edge.length = F, show.tip.label = F, edge.color = "#6D214F")

DotPlot(neuron, features = rev(unique(tip.marker$gene)), cols = c('#dfe4ea','#6F1E51'), col.min = 0, col.max = , dot.min = 0, dot.scale = 2) +
  scale_color_gradientn(colours = rev(mako.col)[1:8]) &coord_flip() & NoLegend()
ggsave("./neuron2/feature.pdf", width = 7, height = 4)

DotPlot(neuron.regulon, features = rev(unique(tip.marker.regulon$gene)), cols = c('#dfe4ea','#1e3799'), col.min = 0, col.max = 2, dot.min = 0, dot.scale = 2, scale = T ) +
  scale_color_gradientn(colours = rev(rocket.col)[1:9]) & coord_flip() & NoLegend()
ggsave("./neuron2/regulon.pdf", width = 7, height = 4)



# ast reclustering --------------------------------------------------------

ast <- subset(hypo, idents = c("Ast1","Ast2", "Ast3", "Tan", "APC"))
ast <- FindVariableFeatures(ast)
ast <- ScaleData(ast)
ast <- RunPCA(ast)
ast <- RunTSNE(ast, seed.use = 2)
ast <- RunUMAP(ast, dims = 1:30, min.dist = 0.3, n.neighbors = 50, n.components = 2)
DimPlot(ast, reduction = "umap", cols = custom_colors$c15[c(3,4,5,8,11)], label = ) & NoAxes()
FeaturePlot(ast, features = "pseudotime") & scale_color_viridis_c() & NoAxes()
ast <- BuildClusterTree(ast)
PlotClusterTree(ast)

astro.markers <- c("GFAP", "SOX9", "ID3", "BHLHE41", "MEIS2", "APOE", "AQP4", "SPARC", "SPARCL1", "SLC1A3", "SLCO1C1", "CRYM", "SCG2", "MANF", "TRIL", "RSAD2", "ISG15", "IFIT1", "STAT1", "CLDN5", "ALCAM", "COL1A2", "RAX", "FGF10", "COL25A1", "GPC3", "DIO2", "LEPR")

VlnPlot(ast, features = astro.markers, pt.size = 0, col = custom_colors$c15[c(3,4,5,8,11)], ncol = 1 ) & NoAxes() & NoLegend()
ggsave("./neuron2/manf.pdf", width = 3, height = 25)

all.markers.ast <- FindAllMarkers(ast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 

View(print(all.markers.ast %>% group_by(cluster) %>% filter(cluster == c(3,4,5,8,11))%>%top_n(n = 10, wt = avg_log2FC)))
astro.deg <- print(all.markers.ast %>% group_by(cluster) %>% filter(cluster == 3 | 4|5|8|11)%>%top_n(n = 10, wt = avg_log2FC))

VlnPlot(ast, c("LEPR","SLCO1C1","CRYM", "SCG2", "ISG15", "IFIT1", "STAT1", "COL25A1", "RAX", "FGF10", "SCN7A"), pt.size = 0)


write_delim(as.data.frame(GetAssayData(hypo, slot = "counts")), delim = "\t", "count_matrix.txt")

