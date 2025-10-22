{ library(dplyr)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(tidyverse)
  library(scran)
  library(viridis)
  library(ggforce)
  library(gghalves)
  library(ggpubr)
  library(scDblFinder)
  library(ggsci)
  library(SingleCellExperiment)
  library(Matrix)
  library(clusterProfiler)
  library(scCustomize)
  library(paletteer)
  library(GOSemSim)
  library(rrvgo)
  library(purrr)
  library(furrr)
}

# plan(multisession, workers = parallel::detectCores() - 4)

custom_colors <- list()
custom_colors$c15 <- c('#bf2c41', '#F79F1F', '#2d98da', '#82589F', '#778ca3', 
                       '#20bf6b', '#1289A7', '#f9ca24', '#A3CB38', '#f8a5c2',
                       '#325A9B', '#822E1C', '#90AD1C', '#B00068', 
                       '#48115e', '#4b7bec')

custom_colors$c4 <- c('#1289A7', '#F79F1F', '#778ca3', '#B00068')
custom_colors$c8 <- c('#4b7bec', '#F79F1F','#82589F', '#20bf6b', '#B00068', '#bf2c41', '#f8a5c2',"#822E1C" )
custom_colors$c3 <- c('#2d98da', '#778ca3', '#B00068')

pal <- viridis(n = 10, option = "D", direction = -1)


#data pre-processing ####
file <- c("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/cellbender_output/11-week__GEX-output_filtered.h5", 
          "/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/cellbender_output/18-month__GEX-output_filtered.h5")

names(file) <- c("young", "old")
set.seed(73)

# param <- MulticoreParam(workers = 6, progressbar = TRUE)
# register(param)

data.list <- list()
qcplot.list <- list()
for (i in 1:length(file)){
  data <- Read_CellBender_h5_Mat(file_name = file[i])
  
  cells <- tibble(
    cell = colnames(data),
    nCount = colSums(data),
    nFeature = colSums(data != 0),
    complexity = log10(colSums(data != 0)/colSums(data)))
  
  cells$percent_ribo <- Matrix::colSums(data[grep("^Rp[sl]", rownames(data)),]) / Matrix::colSums(data)
  cells$percent_mt <- Matrix::colSums(data[grep("^mt-", rownames(data)),]) / Matrix::colSums(data)
  
  #doublet removal
  sce <- SingleCellExperiment(assays = list(counts = data))
  sce <- scDblFinder(sce)
  cells$multiplet_class <- colData(sce)$scDblFinder.class
  
  #prefiltering
  
  thresholds_nCount <- c(3000, 30000)
  thresholds_nFeature <- c(2000, 8000)
  thresholds_percent_mt <- 0.05
  tresholds_percent_ribo <- 0.02
  
  p1 <- ggplot(cells, aes(nCount, nFeature, color = percent_ribo)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = thresholds_nFeature, color = 'red') +
    geom_vline(xintercept = thresholds_nCount, color = 'red') +
    scale_x_continuous(name = 'Number of transcripts', labels = scales::comma) +
    scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma) +
    theme_bw() +
    scale_color_viridis(
      name = 'Percent Ribo\ntranscripts',
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
      percent_mt <= thresholds_percent_mt,
      percent_ribo >= tresholds_percent_ribo,
      multiplet_class == "singlet"
    )
  
  cells_to_keep <- cells_filtered$cell
  
  temp.object <- CreateSeuratObject(counts = data[,cells_to_keep], project = names(file)[i], min.cells = 10, min.features = 5)
  temp.object[["percent.ribo"]] <- PercentageFeatureSet(temp.object, pattern = "^Rp[sl]")
  temp.object[["percent.mito"]] <- PercentageFeatureSet(temp.object, pattern = "^mt-")
  data.list[[paste0("data.",names(file)[i])]] <- temp.object
  qcplot.list[[paste0("data.",names(file)[i])]] <- p1+p2
}


ggarrange(plotlist = qcplot.list, ncol = 1, nrow = 3)
ggsave("./images/supp/qcplot.pdf", width = 8, height = 8)

data.list[[1]]$group <- rep("young", ncol(data.list[[1]])) ; data.list[[2]]$group <- rep("old", ncol(data.list[[2]]))

#making seurat object
merge <- merge(data.list[[1]] , data.list[2] , add.cell.ids = c("young", "old"))
merge@active.assay <- "RNA"
merge <- Add_Cell_QC_Metrics(object = merge, species = "mouse")
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge)
merge <- ScaleData(merge)
merge <- RunPCA(merge)

options(future.globals.maxSize= 5 *1024 ^3)
combined <- IntegrateLayers(object = merge, method = "HarmonyIntegration", new.reduction = "integrated", orig.reduction = "pca",
                            verbose = FALSE)

combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
combined$group <- factor(combined$group, levels = c("young", "old"))
combined <- FindNeighbors(combined, reduction = "integrated", dims = 1:35)
combined <- RunUMAP(combined, reduction = "integrated", graph = "RNA_snn", min.dist = , umap.method = "umap-learn", 
                    n.epochs = 200) 
combined <- FindClusters(combined, resolution = 0.28, algorithm = 2, graph.name = "RNA_snn", random.seed = 73)

{umap1 <- DimPlot_scCustom(combined, colors_use = custom_colors$c15, aspect_ratio = 1, repel = F) & NoAxes() & NoLegend()
  umap2 <- DimPlot_scCustom(combined, colors_use = custom_colors$c15, aspect_ratio = 1, label = F) & NoAxes() 
  umap3 <- DimPlot_scCustom(combined, colors_use = custom_colors$c15, aspect_ratio = 1, label = F, split.by = "orig.ident", ncol=2) & NoAxes() & NoLegend()
  umap4 <- DimPlot_scCustom(combined, colors_use = custom_colors$c4, aspect_ratio = 1, label = F, group.by = "group", shuffle = T, alpha = 0.6 ) & NoAxes()
  umap.all <- ggpubr::ggarrange(umap1, umap2, umap3, umap4, ncol = 2, nrow=2)}
umap.all 

ggsave("./images/main/umap_qc.pdf", width = 9, height = 9)

# saveRDS(combined, "combined.rds")


Stacked_VlnPlot(seurat_object = combined, features = c("Pecam1", "Cx3cr1", "Gpr50", "Ccdc153", "Kcnj8", "Acta2", "Agt", "Tcf21",
                                                       "Pf4", "Snap25", "Plvap", "Cd8a"),
                x_lab_rotate = TRUE,
                colors_use = custom_colors$c15)

ggsave("./images/supp/global_vp.pdf", scale = 1, height = 6, width = 6)

FeaturePlot_scCustom(combined, c("Pecam1", "Cx3cr1", "Gpr50", "Ccdc153", "Kcnj8", "Acta2", "Agt", "Tcf21",
                                 "Pf4", "Snap25", "Plvap", "Cd8a"), colors_use = viridis_magma_dark_high, na_cutoff = 0.4, alpha_exp = 1, aspect_ratio = 1) & NoAxes()
ggsave("./images/supp/global_fp.pdf", scale = 1, height = 7, width = 8.5)

QC_Plots_Combined_Vln(combined, plot_boxplot = T, plot_median = F, log = F, pt.size = 0, colors_use = custom_colors$c4, group.by = "group")
ggsave("./images/supp/qc1.jpg", width = 6, height = 8)
QC_Plots_Combined_Vln(combined, plot_boxplot = T, plot_median = F, log = F, pt.size = 0, colors_use = custom_colors$c15)
ggsave("./images/supp/qc2.jpg", width = 15, height = 8)


chemokines <- sort(rownames(combined)[grepl("Ccl|Cxcl", rownames(combined))])
Stacked_VlnPlot(combined, chemokines, split.by = "group", colors_use = custom_colors$c4)
ggsave("images/supp/chemokines.pdf", width = 5.8, height = 11)


# sc-type -----------------------------------------------------------------

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = 'ScTypeDB_HypoMap.xlsx'
gs_list = gene_sets_prepare(db_, "Brain") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
es.max = sctype_score(scRNAseqData = combined@assays$RNA$scale.data, scaled = T, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

cL_resutls = do.call("rbind", lapply(unique(Idents(combined)), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(combined@meta.data[Idents(combined)==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Idents(combined)==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/2] = "Unknown"
print(sctype_scores %>% arrange(cluster))

lapply(c("ggraph","igraph","tidyverse"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); 
edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
{nodes_lvl1 = sctype_scores[,c("cluster", "ncells")] 
  nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
  ccolss= custom_colors$c15
  for (i in 1:length(unique(cL_resutls$cluster))){
    dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; 
    nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), 
                                              ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
  }
  nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1; nodes$ncells = nodes$ncells
  files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); 
  nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
  nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; 
  nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
  mygraph <- graph_from_data_frame(edges, vertices=nodes)
}
# Make the graph

cell_anno <- ggraph(mygraph, layout = 'circlepack', weight=I(sqrt(ncells))) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#ffffff"), colour=I("#dddfeb")), alpha=1) + 
  geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I(Colour)), alpha=0.6) +
  theme_void() + 
  geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#000000"), size = I(sapply(log(ncells,10) , function(x) ifelse(x > 2, x, 0)))))+ 
  geom_node_label(aes(filter=ord==1,  label=realname, colour=I("#000000"), size = I(5), fill=), fontface = "bold", repel = T, force = 80, 
                  max.overlaps = 3, alpha = 0.8)

cell_anno
ggsave("./images/supp/sc_type.pdf", height = 6, width = 6)

#EC subtype = https://www.sciencedirect.com/science/article/pii/S0092867420300623
#artery most specific marker = "Gkn3", vein = "Slc38a5". Smooth muscle only in artery but not in vein.
#plvap = fenestrated endo (https://pmc.ncbi.nlm.nih.gov/articles/PMC10172233); https://www.biorxiv.org/content/10.1101/2021.04.26.441465v1.full#F6

labels <- forcats::fct_collapse(combined$seurat_clusters, "Ven" = 0, "Art" = 1, "Mg" = 2,  "Tan" = 3, "Edy" = 4,  "Peri" = 5, "SMC" = 6,  "Ast" = 7, 
                                "Fib" = 8, "BAM" = 9, "Neu" = 10, "FEC" = 11, "CD8" = 12)
combined <- AddMetaData(combined, labels, col.name = "celltype")
Idents(combined) <- "celltype"

Cluster_Stats_All_Samples(seurat_object = combined, group_by_var = "group")
Proportion_Plot(seurat_object = combined, plot_type = "pie", split.by = "group", colors_use = custom_colors$c15)
Proportion_Plot(seurat_object = combined, plot_type = "bar", split.by = "group", colors_use = custom_colors$c15)
ggsave("./images/supp/cellnumber_bar.pdf", scale = 1.2, width = 3, height = 4)


a <- data.frame(celltype = Idents(combined), group = combined$group)
a %>% group_by(group, celltype) %>% tally()
ggplot(a, aes(x= group, fill = celltype)) +geom_bar(position = "fill") + 
  scale_fill_manual(values = colorspace::darken(custom_colors$c15, 0.05)) + coord_flip() + theme_classic()

ggplot(a, aes(x= factor(celltype, levels = levels(celltype) [1:17]), fill = group)) + geom_bar(position = "fill") + 
  scale_fill_manual(values = colorspace::darken(custom_colors$c4, 0.05)) + theme_classic() + 
  geom_hline(yintercept = 0.5, lty = 2, alpha = 0.9) + ylab("Proportion") + xlab("Cluster") + guides(fill = guide_legend(reverse = F)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
ggsave("./images/main/cellnumber.pdf", scale = 1, width = 5.5, height = 3.5)

a <- as.data.frame(round(prop.table(table(Idents(combined), combined$group), margin = 2) * 100,2)) %>% group_by(Var2) %>%
  mutate(ypos = 100 - (cumsum(Freq)- 0.5*Freq ))
a <- a %>% mutate(Freq2 = ifelse(Freq > 5, Freq, ""))

ggplot(a %>% filter(Var2 == "young"), aes(x= "", y = Freq, fill = Var1)) +geom_bar(stat = "identity") + 
  scale_fill_manual(values = custom_colors$c15) + coord_polar("y", start=0, direction = -1)+ theme_void() +
  geom_text(aes(y = ypos, label = Freq2), color = "white", size=4) +
  ggplot(a %>% filter(Var2 == "old"), aes(x= "", y = Freq, fill = Var1)) +geom_bar(stat = "identity") + 
  scale_fill_manual(values = custom_colors$c15) + coord_polar("y", start=0, direction = -1)+ theme_void() +
  geom_text(aes(y = ypos, label = Freq2), color = "white", size=4)

ggsave("./images/supp/cellnumber_pie.pdf", scale = 1, width = 7, height = 4)


DotPlot(combined, list(Pan = c("Pecam1", "Cdh5", "Cd34","Eng"),
                       Vein = c("Slc38a5", "Tmsb10", "Ctla2a", "Ctsc"),
                       Artery = c("Gkn3", "Mgp", "Sema3g" ,"Stmn2"), 
                       FEC = c("Plvap", "Col15a1", "Exoc3l2", "Esm1")),
        dot.min = 0.1, dot.scale = 5, scale = T, idents = c("Ven", "Art", "FEC")) + 
  scale_color_viridis(option = "plasma", direction = 1, limits = c(0,1.5), oob = scales::squish) + scale_x_discrete(name = NULL, guide = guide_axis(angle = 90))
ggsave("./images/supp/artery.pdf", scale = 1, width = 7.2, height = 3.5)

# https://www.nature.com/articles/s42255-022-00674-x#code-availability


#ifn_go
ifnb_go_genes <- unlist(strsplit(read.delim("REACTOME_INTERFERON_GAMMA_SIGNALING.v2025.1.Mm.tsv")[17,2], ","))
ifnb_go_genes <-  list(intersect(rownames(combined),ifnb_go_genes))
combined <- AddModuleScore(combined, features = ifnb_go_genes)
VlnPlot(combined, "Cluster1", split.by = "group") + stat_compare_means()

df <- data.frame(celltype = combined$celltype, group = combined$group, ifn = combined$Cluster1)
df %>%
  group_by(celltype, group) %>%
  summarise(
    n = n(),
    mean = mean(ifn),
    median = median(ifn),
    sd = sd(ifn),
    min = min(ifn),
    max = max(ifn),
    .groups = "drop"
  )






#cluster level DEG

all.markers <- FindAllMarkers(combined, assay = "RNA", slot = "data", min.pct = 0.25, only.pos = T)
top50 <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
openxlsx::write.xlsx(top50, "./supp/cluster_marker_top50.xlsx")


DoHeatmap(subset(combined, downsample = 50), features= as.character(top50$gene), group.colors = custom_colors$c15, slot = "scale.data", 
          size = 3, disp.min = -4, disp.max = 3, raster = F, label = F) + scale_fill_viridis(option = "D") + theme(axis.text.y = element_blank())
ggsave("./images/supp/cluster_marker_heatmap.pdf", scale= 1, height = 8, width = 6)


VlnPlot(combined, c("Ccl3", "Ccl4"), idents = "Mg", pt.size = 0, cols = custom_colors$c3, split.by = "group")
ggsave("./images/main/sc_mg_ccl34.pdf", width = 4, height = 3)


# DE ---------------------------------------------------------------------
Idents(combined) <- "celltype"
a <-combined@meta.data$group
b <- Idents(combined)
c <- paste(b,a,sep = "_")
combined@meta.data$deg.compare <- as.factor(c)
table <- data.frame("young" = paste(levels(combined), "young", sep = "_"), "old" = paste(levels(combined), "old", sep="_") )
Idents(combined) <- "deg.compare"


deg.per.cluster <- list()
for(i in 1:nrow(table)) {
  deg.per.cluster[[i]] <-  FindMarkers(combined, ident.2 = table[i,1], ident.1 = table[i,2], test.use = ,
                                       logfc.threshold = , only.pos = F, min.pct = 0.1)}

deg.per.cluster.sig.up <- lapply(deg.per.cluster, function(x) {
  return(x %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.585))})
deg.per.cluster.sig.down <- lapply(deg.per.cluster, function(x) {
  return(x %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < -0.585))})

a <- unlist(lapply(deg.per.cluster.sig.up, nrow))
b <- -1*unlist(lapply(deg.per.cluster.sig.down, nrow))
c <- data.frame(count = c(a,b),
                cell = rep(as.character(levels(combined$celltype))[1:13],2),
                dir = rep(c('up','down'), each = 13))


ggplot(c, aes(x=factor(cell, levels = cell[1:13]), y=count, fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + 
  scale_x_discrete(name = "Cluster", guide = guide_axis(angle = 90)) + scale_fill_jama() + geom_hline(yintercept = 0)
ggsave("./images/main/deg_count_mast.pdf", width = 3.3, height = 3)

ggplot(c, aes(x=factor(cell, levels = cell[1:13]), y=abs(count), fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + 
  scale_x_discrete(name = "Cluster") + scale_fill_jama() + geom_hline(yintercept = 0)
ggsave("./images/main/deg_count_abs_mast.pdf", width = 5.3)

names(deg.per.cluster.sig.up) <- levels(combined$celltype)
names(deg.per.cluster.sig.down) <- levels(combined$celltype)
openxlsx::write.xlsx(deg.per.cluster.sig.up, "./supp/combined_DEG_up.xlsx", rowNames = T)
openxlsx::write.xlsx(deg.per.cluster.sig.down, "./supp/combined_DEG_down.xlsx", rowNames = T)


# GO Identity --------------------------------------------------------------
Idents(combined) <- "celltype"


#function for calculating simplified GO
run_enrich_rrvgo <- safely(function(deg_obj, universe, orgdb = "org.Mm.eg.db", threshold = 0.5) {
  if (is.null(deg_obj) || nrow(deg_obj) == 0) {
    return(list(ego = NULL, reduced = tibble()))
  }
  
  gene_list <- rownames(deg_obj)
  
  x <- enrichGO(
    gene          = gene_list,
    universe      = universe,
    OrgDb         = "org.Mm.eg.db",
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = FALSE
  )
  
  df <- as.data.frame(x)
  
  if (!is.null(df) && nrow(df) > 2) {
    
    simMatrix <- calculateSimMatrix(
      df$ID,
      orgdb = orgdb,
      ont = "BP",
      method = "Rel"
    )
    
    scores <- setNames(-log10(df$p.adjust), df$ID)
    
    reducedTerms <- reduceSimMatrix(
      simMatrix,
      scores = scores,
      threshold = threshold,
      orgdb = orgdb
    )
    
    # Add back p.adjust and Count / geneID
    reducedTerms <- reducedTerms %>%
      mutate(p.adjust = 10^(-score)) %>%
      left_join(df %>% dplyr::select(ID, Count, geneID), by = c("go" = "ID"))
    
    return(list(ego = x, reduced = reducedTerms))
    
  } else {
    return(list(ego = x, reduced = tibble()))
  }
})

clusters <- 1:length(levels(combined))
universe_genes <- rownames(combined)

#function for extracting reduced GO
extract_reduced <- possibly(
  function(x) x$result$reduced,
  otherwise = NULL
)


all.markers.list <- all.markers %>% dplyr::group_by(cluster, .add = T) %>% 
  top_n(200, wt = avg_log2FC) %>%   group_split()
all.markers.list <- lapply(all.markers.list, function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x$gene
  return(x)
})

all.go.list <- future_map(clusters, ~ run_enrich_rrvgo(
  deg_obj = all.markers.list[[.x]],
  universe = universe_genes
), .progress = TRUE)


all.go <- future_map(all.go.list, possibly(function(x)
{x <- as.data.frame(x$result$ego)
return(x)}))
names(all.go) <- levels(combined)


all.go2 <- imap_dfr(all.go, ~ mutate(.x, celltype = .y)) %>% mutate(celltype = factor(celltype, levels = levels(combined))) %>%
  group_by(celltype) %>% slice_min(n = 4, p.adjust, with_ties = F)

data <- data.frame(
  term = all.go2$Description,
  significance = -log10(all.go2$p.adjust),
  count = all.go2$Count,
  col = factor(all.go2$celltype, levels = levels(combined$celltype)),
  term2 = paste(all.go2$Description, all.go2$celltype, sep = "_")
)
data <- data[nrow(data):1,]


legend_col <- custom_colors$c15[1:13]
names(legend_col) <- levels(combined$celltype)


ggplot(data, aes(x= fct_inorder(term2), y=significance)) +
  geom_segment(aes(xend=term2, y=0, yend= significance), color = rep(custom_colors$c15[13:1], each = 4), lwd = 1) +
  geom_point(aes(size= count, color = col)) + scale_color_manual(values = legend_col) +
  coord_flip() + theme_bw() + scale_x_discrete(labels = data$term, name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2))
ggsave("./images/supp/go_all_style1.pdf", width = 12, height = 8)

ggplot(data, aes(x=factor(term2, levels = term2), y=significance, color = term)) +
  geom_segment(aes(xend=term2, y=0, yend=significance), color = rep(custom_colors$c15[13:1], each = 4), lwd = 1.5) +
  geom_point(aes(size= count), color = rep(custom_colors$c15[13:1], each = 4)) +
  coord_flip() + theme_classic() +  facet_grid(col~., scales = "free", space = "free") +
  scale_x_discrete(labels=setNames(data$term, data$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") +
  labs(y="-log10(p-value)", x="GO Term") + 
  theme(legend.position = "bottom", strip.text = element_text(color = "black", face = "bold", size = 7), 
        strip.background = element_rect(fill = "gray90"))
ggsave("./images/supp/go_all_style2.pdf", width = 10, height = 11)




# GO DEG ------------------------------------------------------------------

ego.up.list <- future_map(clusters, ~ run_enrich_rrvgo(
  deg_obj = deg.per.cluster.sig.up[[.x]],
  universe = universe_genes
), .progress = TRUE)

ego.down.list <- future_map(clusters, ~ run_enrich_rrvgo(
  deg_obj = deg.per.cluster.sig.down[[.x]],
  universe = universe_genes
), .progress = TRUE)

ego.up.reduced <- future_map(ego.up.list, extract_reduced)
ego.down.reduced <- future_map(ego.down.list, extract_reduced)

ego.up.reduced <- map(ego.up.reduced, ~ {
  if (nrow(.x) == 0) return(.x)
  
  .x %>%
    filter(termDispensability < 0.5) %>%
    distinct(parentTerm, .keep_all = TRUE)
})
ego.down.reduced <- map(ego.down.reduced, ~ {
  if (nrow(.x) == 0) return(.x)
  
  .x %>%
    filter(termDispensability < 0.5) %>%
    distinct(parentTerm, .keep_all = TRUE)
})

names(ego.up.reduced) <- levels(combined)[seq_along(1:length(ego.up.reduced))]
names(ego.down.reduced) <- levels(combined)[seq_along(1:length(ego.down.reduced))]

openxlsx::write.xlsx(ego.up.reduced, "./supp/combined_GO_up.xlsx")
openxlsx::write.xlsx(ego.down.reduced, "./supp/combined_GO_down.xlsx")

color_lookup <- data.frame(cluster = levels(combined$celltype), col = custom_colors$c15[1:13]) 

dataup <- imap_dfr(ego.up.reduced, ~ mutate(.x, celltype = .y)) %>%  mutate(term2 = paste(parentTerm, celltype, sep = "_")) 
datadown <- imap_dfr(ego.down.reduced, ~ mutate(.x, celltype = .y)) %>%  mutate(term2 = paste(parentTerm, celltype, sep = "_")) 

data1 <- dataup %>%
  mutate(
    celltype = factor(celltype, levels = levels(combined$celltype)),
    term2 = stringr::str_wrap(term2, width = 70),
    term2 = forcats::fct_inorder(term2)
  ) %>%
  group_by(celltype) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(term2 = forcats::fct_rev(term2))

data2 <- datadown %>%
  mutate(
    celltype = factor(celltype, levels = levels(combined$celltype)),
    term2 = stringr::str_wrap(term2, width = 70),
    term2 = forcats::fct_inorder(term2)
  ) %>%
  group_by(celltype) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(term2 = forcats::fct_rev(term2))


legend_col <- custom_colors$c15
names(legend_col) <- levels(combined$celltype)

ggplot(data1, aes(x = fct_reorder(term2, -p.adjust), y = -log10(p.adjust))) +
  geom_segment(aes(xend=term2, y=0, yend=-log10(p.adjust), color = celltype), lwd = 1) +
  geom_point(aes(size= Count, color = celltype)) + scale_color_manual(values = legend_col) +
  scale_x_discrete(labels = setNames(sub("_(.*)", "", data1$term2), data1$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") + coord_flip() +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2)) +
  facet_grid(celltype~., scales = "free", space = "free") + theme_classic()
ggsave("./images/supp/go_deg_up_style1.pdf", width = 12, height = 15)

ggplot(data2, aes(x = fct_reorder(term2, -p.adjust), y = -log10(p.adjust))) +
  geom_segment(aes(xend=term2, y=0, yend=-log10(p.adjust), color = celltype), lwd = 1) +
  geom_point(aes(size= Count, color = celltype)) + scale_color_manual(values = legend_col) +
  scale_x_discrete(labels = setNames(sub("_(.*)", "", data2$term2), data2$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") + coord_flip() +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2)) +
  facet_grid(celltype~., scales = "free", space = "free") + theme_classic()
ggsave("./images/supp/go_deg_down_style1.pdf", width = 12, height = 12)

data3 <- rbind(data1 %>% mutate(dir = "up"), data2 %>% mutate(dir = "down"))
data3$dir <- factor(data3$dir, levels = unique(data3$dir))

ggplot(data3, aes(x=term2, y=-log10(p.adjust))) +
  geom_segment(aes(xend=term2, y=0, yend=-log10(p.adjust), color = celltype), lwd = 1) +
  geom_point(aes(size= Count, color = celltype)) + scale_color_manual(values = legend_col) +
  coord_flip() + theme_bw() + scale_x_discrete(labels = setNames(sub("_(.*)", "", data3$term2), data3$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2)) + 
  facet_grid(dir~., scales = "free", space = "free") 
ggsave("./images/supp/go_deg_combined_style1.pdf", width = 12, height = 30)


data4 <- data1 %>% filter(celltype == "Mg" | celltype == "Tan") %>% mutate(across(everything(), ~ .[rev(seq_along(.))]))
data4$term2 <- fct_inorder(data4$term2)

ggplot(data4, aes(x=term2, y=-log10(p.adjust))) +
  geom_segment(aes(xend=term2, y=0, yend=-log10(p.adjust), color = celltype), lwd = 1.5) +
  geom_point(aes(size= Count, color = celltype)) +  scale_color_manual(values = legend_col) +
  coord_flip() + theme_classic() + scale_x_discrete(labels = setNames(sub("_(.*)", "", data4$term2), data4$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + facet_grid(celltype~., scales = "free", space = "free") 

ggsave("./images/main/go_deg_tanmg_style1.pdf", width = 7.5, height = 5.5)


# Global GO dotplot -------------------------------------------------------

all_ego_up <- map_dfr(seq_along(ego.up.list), function(i) {
  res <- ego.up.list[[i]]
  if (is.null(res$result$ego)) return(tibble())
  df <- as.data.frame(res$result$ego)
  if (nrow(df) == 0) return(tibble())
  df$cluster_id <- i
  df$direction <- "UP"
  return(df)
})
all_ego_down <- map_dfr(seq_along(ego.down.list), function(i) {
  res <- ego.down.list[[i]]
  if (is.null(res$result$ego)) return(tibble())
  df <- as.data.frame(res$result$ego)
  if (nrow(df) == 0) return(tibble())
  df$cluster_id <- i
  df$direction <- "DOWN"
  return(df)
})
all_ego_df <- bind_rows(all_ego_up, all_ego_down)
all_go_terms <- unique(all_ego_df$ID)

simMatrix_global <- calculateSimMatrix(
  all_go_terms,
  orgdb = "org.Mm.eg.db",
  ont = "BP",
  method = "Wang"
)

global_scores <- all_ego_df %>%
  group_by(ID) %>%
  summarise(score = -log10(min(p.adjust))) %>%
  deframe()

reducedTerms_global <- reduceSimMatrix(
  simMatrix_global,
  scores = global_scores,
  threshold = 0.92,
  orgdb = "org.Mm.eg.db"
) %>%
  mutate(p.adjust = 10^(-score))



dist_mat <- as.dist(1 - simMatrix_global[reducedTerms_global$go,reducedTerms_global$go])
hc <- hclust(dist_mat, method = "complete")
go_order <- hc$labels[hc$order]

id_to_parent <- reducedTerms_global %>%
  dplyr::select(go, parentTerm) %>%
  distinct()

parent_order <- id_to_parent %>%
  filter(go %in% go_order) %>%
  arrange(match(go, go_order)) %>%
  pull(parentTerm) %>%
  make.unique(sep = " ")

plot_data <- reducedTerms_global %>%
  left_join(
    all_ego_df %>%
      dplyr::select(ID, cluster_id, direction, Count) %>%
      distinct(ID, cluster_id, direction, Count),
    by = c("go" = "ID")
  ) %>%
  mutate(
    parentTerm = factor(parentTerm, levels = parent_order),
    cluster = factor(cluster_id, 
                     levels = seq_along(levels(combined$celltype)), 
                     labels = levels(combined$celltype)),
    direction = factor(direction, levels = c("UP", "DOWN"))
  ) %>%
  group_by(parentTerm, cluster, direction) %>%
  summarise(
    p.adjust = min(p.adjust, na.rm = TRUE),
    Count    = max(Count, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(logFDR = pmin(-log10(p.adjust), 6))

ggplot(plot_data, aes(x = cluster, y = fct_inorder(stringr::str_wrap(parentTerm, width = 70)))) +
  geom_point(aes(size = Count^1.5, color = logFDR)) +
  scale_color_viridis_c(option = "viridis", direction = 1, name = "-log10(FDR)", limits = c(1, 6)) +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Cluster",
    y = "GO term (RRVGO parent, HC ordered)",
    title = "Global GO terms - Reduced by RRVGO - HC ordered Y"
  ) +
  facet_wrap(~ direction) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  )
ggsave("./images/supp/global_go.pdf", height = 10, width = 12)


ggplot(plot_data, aes(x = cluster, y = fct_inorder(stringr::str_wrap(parentTerm, width = 70)))) +
  geom_point(aes(size = Count^1.5, color = logFDR)) +
  scale_color_viridis_c(option = "viridis", direction = 1, name = "-log10(FDR)", limits = c(1, 6)) +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Cluster",
    y = "GO term (RRVGO parent, HC ordered)",
    title = "Global GO terms - Reduced by RRVGO - HC ordered Y"
  ) +
  facet_wrap(~ direction) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  )
ggsave("./images/supp/global_go.pdf", height = 10, width = 12)


selected_terms <- c("defense response to virus", "hypersensitivity", "response to interferon-beta", "cellular response to interferon-beta", 
                    "positive regulation of type II hypersensitivity", "peptide antigen assembly with MHC class I protein complex",
                    "antigen processing and presentation of endogenous peptide antigen via MHC class I", "negative regulation of viral process",
                    "positive regulation of natural killer cell cytokine production")

ggplot(plot_data %>% filter(parentTerm %in% selected_terms), aes(x = cluster, y = fct_inorder(stringr::str_wrap(parentTerm, width = 70)))) +
  geom_point(aes(size = Count^1.5, color = logFDR)) +
  scale_color_viridis_c(option = "viridis", direction = 1, name = "-log10(FDR)", limits = c(1, 6)) +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Cluster",
    y = "GO term (RRVGO parent, HC ordered)",
    title = "Global GO terms - Reduced by RRVGO - HC ordered Y"
  ) +
  facet_wrap(~ direction) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  )
ggsave("./images/main/global_go_selected.pdf", height = 5, width = 12)

VlnPlot(combined, idents = c("Mg", "Tan"), c("B2m", "H2-K1", "H2-D1", "Tap2", "Bst2", "Ifitm3", "Ifi1", "Irf1"), split.by = "group", 
        cols = custom_colors$c4, pt.size = 0, log=F, ncol = 2) & RestoreLegend()

FetchData(subset(combined, idents = "Mg"), vars = c("Ifitm3", "group"), layer = "data") %>% 
  ggplot(aes(y = Ifitm3, x = group, col = group, lower = 0.5, upper = 0.5)) + geom_boxplot(notch =)


VlnPlot(combined, idents =, c("Ccr5", "Ccr2", "Cd8a", "Tbx21", "Il2rb","Ifng"), cols = custom_colors$c15, pt.size = 0, log=F, ncol = 3)


# Tan subclustering -------------------------------------------------

Idents(combined) <- "celltype"
combined$orig.ident
tan <- subset(combined, idents = c("Tan"))
tan <- FindNeighbors(tan, dims = 1:35, reduction = "integrated")
tan <- FindClusters(tan, resolution = 0.4, algorithm = 4, graph.name = "RNA_nn", random.seed = 42)
tan <- RunUMAP(tan, umap.method = "umap-learn", graph = "RNA_nn", n.epoch = 200, min.dist = , seed.use = 42)

tan.umap1 <- DimPlot_scCustom(tan, colors_use = custom_colors$c8, aspect_ratio = 1, label.box = F)
tan.umap2 <- DimPlot_scCustom(tan, colors_use = custom_colors$c8, aspect_ratio = 1, label = F) 
ggarrange(tan.umap1, tan.umap2)
ggsave("./images/main/tan_dimplot.pdf", width = 7, height = 3, device = cairo_pdf)

FeaturePlot(tan ,c("Scn7a", "Col25a1", "Ptn", "Slc17a8", "Mafb", "Pcp4", "nFeature_RNA"))

Proportion_Plot(tan, split.by = "group", plot_type = "pie", colors_use = custom_colors$c8)
DimPlot_scCustom(tan, group.by = "orig.ident", colors_use = custom_colors$c4, pt.size = 1, alpha = 1, aspect_ratio = 1, reduction = "umap") & NoAxes()
ggsave("./images/main/tan_dimplot_compare.pdf", width = 4.1, height = 3)


Cluster_Stats_All_Samples(seurat_object = tan, group_by_var = "group")
labels <- forcats::fct_collapse(tan$seurat_clusters,
                                "TDN" = 1,
                                "α2" = 2,
                                "α1" = 3,
                                "β1" = 4,
                                "β2" = 5)

tan <- AddMetaData(tan, labels, col.name = "subtype")
tan$subtype <- factor(tan$subtype, levels = c("β2","β1", "α2", "α1", "TDN"))
Idents(tan) <- "subtype"


all.markers.tan <- FindAllMarkers(tan, latent.vars = , only.pos = T)
View(all.markers.tan %>% group_by(cluster) %>% slice_head(n=50))


Proportion_Plot(seurat_object = tan, plot_type = "pie", split.by = "group", colors_use = custom_colors$c8)
ggsave("./images/main/tan_cellnumber.pdf", width = 5, height = 4, device = cairo_pdf)


#comparison with brain aging spatial seq
VlnPlot_scCustom(combined, idents = "Tan", c("Oasl2", "Ifit1", "H2-K1", "Ifi27", "Ccnd2", "Ctnna2"), colors_use = custom_colors$c4, 
                 split.by = "group", pt.size = 0)
ggsave("./images/main/tan_comparison_Jin.pdf", width = 5, height = 4, device = cairo_pdf)


#https://www.biorxiv.org/content/10.1101/2023.07.06.547914v1.full.pdf
#https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.1129414/full
#https://www.biorxiv.org/content/10.1101/2020.11.02.359992v3.full
#tanycyte classifications
#pan-tan = Rax, Col23a1; pan-Epen = Ccdc153
#b2-tan = Col25a1, Scn7a, Ptn; b1-tan = Frzb, Ptn; a2-tan = Vcan, Frzb, Mafb; a1-tan = Slc17a8
#allen brain atlas IHC: Col25a1, Ptn, Frzb, Vcan, Gpr50, Mafb?, Slc17a8

VlnPlot(tan, c("Rax", "Gpr50", "Scn7a", "Adm", "Col25a1", "Ptn", "Frzb", "Pdzph1", "Vcan", "Rspo3", "Mafb", "Slc17a8", "Pcp4", "Tenm4"), 
        ncol = 4, cols = custom_colors$c8, pt.size = 0)
ggsave("./images/supp/tan_vp.pdf", scale = 1, height = 8.2, width = 9)

Stacked_VlnPlot(tan, c("Rax", "Gpr50", "Scn7a","Adm", "Col25a1", "Ptn", "Frzb", "Pdzph1", "Vcan", "Rspo3", "Mafb", "Slc17a8", "Pcp4", "Tenm4"), 
                x_lab_rotate = TRUE,
                colors_use = custom_colors$c8)
ggsave("./images/main/tan_vp_stack.pdf", scale = 1, height = 5.8, width = 3.5, device = cairo_pdf)


#https://pubmed.ncbi.nlm.nih.gov/12766484/ apoptosis inducing interferon response genes
ifn.apo <- c("Fas", "Xaf1","Eif2ak2", "Dapk1", "Casp8", "Oasl2","Plscr1", "Dido1")
Stacked_VlnPlot(tan, ifn.apo, x_lab_rotate = TRUE, colors_use = custom_colors$c3, split.by = "group") & RestoreLegend()
ggsave("./images/main/tan_IEG_fp.pdf", scale = 1, height = 6, width = 6, device = cairo_pdf)

#https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/REACTOME_APOPTOSIS.html REACTOME APOPTOSIS gene set
line <- readLines("REACTOME_APOPTOSIS.v2025.1.Mm.gmt")
gene_line <- strsplit(line[1], "\t")[[1]]
genes_reactome_apoptosis <- gene_line[-c(1,2)]
genes_filtered <- genes_reactome_apoptosis[genes_reactome_apoptosis %in% rownames(tan)]

tan <- AddModuleScore(tan, features = list(genes_filtered), name = "Apoptosis_REACTOME")
VlnPlot_scCustom(
  tan, features = "Apoptosis_REACTOME1", plot_box = TRUE, group.by = "group"
) + stat_compare_means()
ggsave("./images/main/tan_IEG_fp.pdf", scale = 1, height = 6, width = 6, device = cairo_pdf)



# saveRDS(tan, "tan.rds")

# Tan DEG -----------------------------------------------------------------
Idents(tan) <- "subtype"
a <- tan@meta.data$group
b <- Idents(tan)
c <- paste(b,a,sep = "_")
tan@meta.data$deg.compare <- as.factor(c)
table <- data.frame("young" = paste(levels(tan), "young", sep = "_"), "old" = paste(levels(tan), "old", sep="_") )
Idents(tan) <- "deg.compare"

deg.per.cluster.tan <- list()
for(i in 1:nrow(table)) {
  deg.per.cluster.tan[[i]] <-  FindMarkers(tan, ident.2 = table[i,1], ident.1 = table[i,2], test.use = ,
                                           logfc.threshold = , min.cells.group = ,
                                           pseudocount.use = , only.pos = F, latent.vars = )}


deg.per.cluster.tan.sig.up <- lapply(deg.per.cluster.tan, function(x) {
  return(x %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.585))})
deg.per.cluster.tan.sig.down <- lapply(deg.per.cluster.tan, function(x) {
  return(x %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < -0.585))})

a <- unlist(lapply(deg.per.cluster.tan.sig.up, nrow))
b <- -1*unlist(lapply(deg.per.cluster.tan.sig.down, nrow))
c <- data.frame(count = c(a,b),
                cell = rep(as.character(levels(tan$subtype))[1:5],2),
                dir = rep(c('up','down'), each = 5))


ggplot(c, aes(x=factor(cell, levels = cell[1:5]), y=count, fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + 
  scale_x_discrete(name = "Cluster", guide = guide_axis(angle = 90)) + scale_fill_jama() + geom_hline(yintercept = 0)
ggsave("./images/main/tan_deg_count_mast.pdf", width = 3.3, height = 3, device = cairo_pdf)

ggplot(c, aes(x=factor(cell, levels = cell[1:5]), y=abs(count), fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + 
  scale_x_discrete(name = "Cluster") + scale_fill_jama() + geom_hline(yintercept = 0)
ggsave("./images/main/tan_deg_count_abs_mast.pdf", width = 5.3, device = cairo_pdf)

names(deg.per.cluster.tan.sig.up) <- levels(tan$subtype)
names(deg.per.cluster.tan.sig.down) <- levels(tan$subtype)
openxlsx::write.xlsx(deg.per.cluster.tan.sig.up, "./supp/tan_DEG_MAST_up.xlsx", rowNames = T)
openxlsx::write.xlsx(deg.per.cluster.tan.sig.down, "./supp/tan_DEG_MAST_down.xlsx", rowNames = T)


# Tan GO ------------------------------------------------------------------
Idents(tan) <- "subtype"
ego.tan.up.list <- future_map(1:5, ~ run_enrich_rrvgo(
  deg_obj = deg.per.cluster.tan.sig.up[[.x]],
  universe = universe_genes
), .progress = TRUE)

ego.tan.down.list <- future_map(1:5, ~ run_enrich_rrvgo(
  deg_obj = deg.per.cluster.tan.sig.down[[.x]],
  universe = universe_genes
), .progress = TRUE)

ego.up.reduced <- future_map(ego.tan.up.list, extract_reduced)
ego.down.reduced <- future_map(ego.tan.down.list, extract_reduced)

ego.up.reduced <- map(ego.up.reduced, ~ {
  if (nrow(.x) == 0) return(.x)
  
  .x %>%
    filter(termDispensability < 0.5) %>%
    distinct(parentTerm, .keep_all = TRUE)
})
ego.down.reduced <- map(ego.down.reduced, ~ {
  if (nrow(.x) == 0) return(.x)
  
  .x %>%
    filter(termDispensability < 0.5) %>%
    distinct(parentTerm, .keep_all = TRUE)
})

names(ego.up.reduced) <- levels(tan)[seq_along(1:length(ego.up.reduced))]
names(ego.down.reduced) <- levels(tan)[seq_along(1:length(ego.down.reduced))]

openxlsx::write.xlsx(ego.up.reduced, "./supp/tan_GO_up.xlsx")
openxlsx::write.xlsx(ego.down.reduced, "./supp/tan_GO_down.xlsx")

color_lookup <- data.frame(cluster = levels(tan$subtype), col = custom_colors$c8[1:5]) 
dataup <- imap_dfr(ego.up.reduced, ~ mutate(.x, subtype = .y)) %>%  mutate(term2 = paste(parentTerm, subtype, sep = "_")) 
datadown <- imap_dfr(ego.down.reduced, ~ mutate(.x, subtype = .y)) %>%  mutate(term2 = paste(parentTerm, subtype, sep = "_")) 

data1 <- dataup %>%
  mutate(
    subtype = factor(subtype, levels = levels(tan$subtype)),
    term2 = stringr::str_wrap(term2, width = 70),
    term2 = forcats::fct_inorder(term2)
  ) %>%
  group_by(subtype) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(term2 = forcats::fct_rev(term2))

data2 <- datadown %>%
  mutate(
    subtype = factor(subtype, levels = levels(tan$subtype)),
    term2 = stringr::str_wrap(term2, width = 70),
    term2 = forcats::fct_inorder(term2)
  ) %>%
  group_by(subtype) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(term2 = forcats::fct_rev(term2))


legend_col <- custom_colors$c8[1:5]
names(legend_col) <- levels(tan$subtype)

ggplot(data1, aes(x = fct_reorder(term2, -p.adjust), y = -log10(p.adjust))) +
  geom_segment(aes(xend=term2, y=0, yend=-log10(p.adjust), color = subtype), lwd = 1) +
  geom_point(aes(size= Count, color = subtype)) + scale_color_manual(values = legend_col) +
  scale_x_discrete(labels = setNames(sub("_(.*)", "", data1$term2), data1$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") + coord_flip() +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2)) +
  facet_grid(subtype~., scales = "free", space = "free") + theme_classic()
ggsave("./images/supp/tan_go_deg_up_style1.pdf", width = 10, height = 8, device = cairo_pdf)


ggplot(data2, aes(x = fct_reorder(term2, -p.adjust), y = -log10(p.adjust))) +
  geom_segment(aes(xend=term2, y=0, yend=-log10(p.adjust), color = subtype), lwd = 1) +
  geom_point(aes(size= Count, color = subtype)) + scale_color_manual(values = legend_col) +
  scale_x_discrete(labels = setNames(sub("_(.*)", "", data2$term2), data2$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") + coord_flip() +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2)) +
  facet_grid(subtype~., scales = "free", space = "free") + theme_classic()
ggsave("./images/supp/tan_go_deg_down_style1.pdf", width = 10, height = 7, device = cairo_pdf)


all_ego_up <- map_dfr(seq_along(ego.tan.up.list), function(i) {
  res <- ego.tan.up.list[[i]]
  if (is.null(res$result$ego)) return(tibble())
  df <- as.data.frame(res$result$ego)
  if (nrow(df) == 0) return(tibble())
  df$cluster_id <- i
  df$direction <- "UP"
  return(df)
})
all_ego_down <- map_dfr(seq_along(ego.tan.down.list), function(i) {
  res <- ego.tan.down.list[[i]]
  if (is.null(res$result$ego)) return(tibble())
  df <- as.data.frame(res$result$ego)
  if (nrow(df) == 0) return(tibble())
  df$cluster_id <- i
  df$direction <- "DOWN"
  return(df)
})
all_ego_df <- bind_rows(all_ego_up, all_ego_down)
all_go_terms <- unique(all_ego_df$ID)

simMatrix_global <- calculateSimMatrix(
  all_go_terms,
  orgdb = "org.Mm.eg.db",
  ont = "BP",
  method = "Wang"
)

global_scores <- all_ego_df %>%
  group_by(ID) %>%
  summarise(score = -log10(min(p.adjust))) %>%
  deframe()

reducedTerms_global <- reduceSimMatrix(
  simMatrix_global,
  scores = global_scores,
  threshold = 0.93,
  orgdb = "org.Mm.eg.db"
) %>%
  mutate(p.adjust = 10^(-score))



dist_mat <- as.dist(1 - simMatrix_global[reducedTerms_global$go,reducedTerms_global$go])
hc <- hclust(dist_mat, method = "complete")
go_order <- hc$labels[hc$order]

id_to_parent <- reducedTerms_global %>%
  dplyr::select(go, parentTerm) %>%
  distinct()

parent_order <- id_to_parent %>%
  filter(go %in% go_order) %>%
  arrange(match(go, go_order)) %>%
  pull(parentTerm) %>%
  make.unique(sep = " ")

plot_data <- reducedTerms_global %>%
  left_join(
    all_ego_df %>%
      dplyr::select(ID, cluster_id, direction, Count) %>%
      distinct(ID, cluster_id, direction, Count),
    by = c("go" = "ID")
  ) %>%
  mutate(
    parentTerm = factor(parentTerm, levels = parent_order),
    cluster = factor(cluster_id, 
                     levels = seq_along(levels(tan$subtype)), 
                     labels = levels(tan$subtype)),
    direction = factor(direction, levels = c("UP", "DOWN"))
  ) %>%
  group_by(parentTerm, cluster, direction) %>%
  summarise(
    p.adjust = min(p.adjust, na.rm = TRUE),
    Count    = max(Count, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(logFDR = pmin(-log10(p.adjust), 6))

ggplot(plot_data, aes(x = cluster, y = fct_inorder(stringr::str_wrap(parentTerm, width = 70)))) +
  geom_point(aes(size = Count^1.5, color = logFDR)) +
  scale_color_viridis_c(option = "viridis", direction = 1, name = "-log10(FDR)", limits = c(1, 4)) +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Cluster",
    y = "GO term (RRVGO parent, HC ordered)",
    title = "Global GO terms - Reduced by RRVGO - HC ordered Y"
  ) +
  facet_wrap(~ direction) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 12),
    panel.spacing = unit(1, "lines")
  )

ggsave("./images/main/tan_global_go.pdf", height = 6, width = 12, device = cairo_pdf)


# Mg subclustering ----------------------------------------------------

mg <- subset(combined, idents = "Mg")
mg@active.assay <- "RNA"
mg <- FindNeighbors(mg, reduction = "integrated" , dims = 1:50, k.param = 10)
mg <- FindClusters(mg, resolution = 0.2, algorithm = 4, graph.name = "RNA_nn", random.seed = 42)
mg <- RunUMAP(mg, umap.method = "umap-learn", graph = "RNA_nn", min.dist = , n.epochs = 200, seed.use = 42)


mg.umap1 <- DimPlot_scCustom(mg, colors_use =  custom_colors$c3, aspect_ratio = 1) & NoLegend() & NoAxes()
mg.umap2 <- DimPlot_scCustom(mg, colors_use =  custom_colors$c3, aspect_ratio = 1, label = F) & NoAxes()
ggarrange(mg.umap1, mg.umap2)
ggsave("./images/main/mg_dimplot.pdf", width = 8.5, height = 4)

DimPlot(mg, split.by = "group", cols = custom_colors$c3, pt.size = 1, alpha = 1) & NoLegend() & NoAxes()
ggsave("./images/main/mg_dimplot_compare.pdf", width = 5.7, height = 3)

labels <- forcats::fct_collapse(mg$seurat_clusters, "AMg" = 1, "HMg" = 2, "LDAM" = 3)
mg <- AddMetaData(mg, labels, col.name = "subtype")
mg$subtype <- factor(mg$subtype, levels = c("HMg", "AMg", "LDAM"))
Idents(mg) <- "subtype"

Cluster_Stats_All_Samples(seurat_object = mg, group_by_var = "group")

Proportion_Plot(seurat_object = mg, plot_type = "pie", split.by = "group", colors_use = custom_colors$c3, num_columns = 1)
ggsave("./images/main/mg_cellnumber1.pdf", width = 5, height = 8)

Proportion_Plot(seurat_object = mg, plot_type = "bar", split.by = "group", colors_use = custom_colors$c3, num_columns = 1)
ggsave("./images/main/mg_cellnumber2.pdf", width = 3, height = 3.5)

all.markers.mg <- FindAllMarkers(mg, only.pos = T, logfc.threshold =, min.diff.pct =)
View(print(all.markers.mg %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)))

mg.markers.top50 <- all.markers.mg %>% group_by(cluster) %>% slice_head(n = 100)
openxlsx::write.xlsx(mg.markers.top50, file = "./supp/mg_marker_top50.xlsx")

genes <- list(Homeostatic = c("Cx3cr1","P2ry12", "P2ry13", "Tgfb1", "Sirpa","Tmem119", "Sall1", 
                              "Siglech"),
              ISR = c("Jun", "Hspa1a", "Atf3",  "Atf4", "Ddit3", "Gadd45b", "Xbp1", "Ppp1r15a"),
              Lipid = c("Apoe", "Plin2", "Plin3", "Abca1", "Soat1", "Cd63", "Grn", "Cebpa", "Cebpb"),
              "DAM & Inflammation" = c("Tnf", "Tyrobp", "Axl", "Cst7", "C4b",
                                       "Il1b", "Ccl3", "Ccl4", "H2-D1", "B2m", "Stat1", "Ifitm3"),
              Senescence = c("Lgals3", "Cdkn1a", "Cdkn2a", "Lmna", "H2afx", "Hmgb1"))


DotPlot(mg, genes, dot.min = 0, dot.scale = 5, scale = T) + 
  scale_color_viridis(option = "plasma", direction = 1, limits = c(0,1.2), oob = scales::squish) + 
  scale_x_discrete(name = NULL, guide = guide_axis(angle = 90)) + scale_y_discrete(limits = rev(levels(mg$subtype)))
ggsave("./images/main/mg_ident_dot.pdf", width = 12, height = 3.6)


sca <- c("Apoe", "Scarb1", "Lrp1", "Lrp6", "Cd44")
Stacked_VlnPlot(mg, sca, log = T, colors_use = custom_colors$c3)
ggsave("./images/supp/mg_scavenger_receptor.pdf", width = 3.5, height = 7)

my_comparisons <- list( c("HMg", "AMg"), c("HMg", "LDAM"), c("AMg", "LDAM"))

plot_list <- VlnPlot(mg, unlist(genes), pt.size = 0, cols = custom_colors$c3, combine = F)

plot_vln <- lapply(plot_list, function(p) {
  y_max <- max(p$data[[1]], na.rm = TRUE)
  p + stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    step.increase = 0.2,
    bracket.size = 0.2,
    size = 3) + ylim(c(0, y_max * 1.3 + 1.3))
})

wrap_plots(plot_vln, ncol =10) & theme(axis.title.x = element_blank(), axis.title.y = element_blank()) &NoLegend()
ggsave("./images/main/mg_multiple_vp.pdf", width = 18, height = 15)


# saveRDS(mg, "mg.rds")

# Mg DEG ------------------------------------------------------------------
Idents(mg) <- "subtype"
a <- mg@meta.data$group
b <- Idents(mg)
c <- paste(b,a,sep = "_")
mg@meta.data$deg.compare <- as.factor(c)
table <- data.frame("young" = paste(levels(mg), "young", sep = "_"), "old" = paste(levels(mg), "old", sep="_") )
Idents(mg) <- "deg.compare"


deg.per.cluster.mg <- list()
for(i in 1:nrow(table)) {
  deg.per.cluster.mg[[i]] <-  FindMarkers(mg, ident.2 = table[i,1], ident.1 = table[i,2], test.use = ,
                                          logfc.threshold = , min.cells.group = ,
                                          pseudocount.use = , only.pos = F, latent.vars = )}


deg.per.cluster.mg.sig.up <- lapply(deg.per.cluster.mg, function(x) {
  return(x %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0.585))})
deg.per.cluster.mg.sig.down <- lapply(deg.per.cluster.mg, function(x) {
  return(x %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < -0.585))})

a <- unlist(lapply(deg.per.cluster.mg.sig.up, nrow))
b <- -1*unlist(lapply(deg.per.cluster.mg.sig.down, nrow))
c <- data.frame(count = c(a,b),
                cell = rep(as.character(levels(mg$subtype))[1:3],2),
                dir = rep(c('up','down'), each = 3))

ggplot(c, aes(x=factor(cell, levels = cell[1:3]), y=count, fill = dir)) + geom_col() + theme_classic() + theme(legend.position="none") + 
  scale_x_discrete(name = "Cluster", guide = guide_axis(angle = 90)) + scale_fill_jama() + geom_hline(yintercept = 0)
ggsave("./images/main/tan_deg_count_mast.pdf", width = 3.3, height = 3, device = cairo_pdf)


# Mg GO Identity --------------------------------------------------------------
Idents(mg) <- "subtype"


all.markers.mg.list <- all.markers.mg %>% dplyr::group_by(cluster, .add = T) %>% 
  top_n(100, wt =-p_val_adj) %>%   group_split()
all.markers.mg.list <- lapply(all.markers.mg.list, function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x$gene
  return(x)
})

all.go.mg.list <- future_map(1:3, ~ run_enrich_rrvgo(
  deg_obj = all.markers.mg.list[[.x]],
  universe = universe_genes
), .progress = TRUE)

all.go <- future_map(all.go.mg.list, possibly(function(x)
{x <- as.data.frame(x$result$ego)
return(x)}))
names(all.go) <- levels(mg)

all.go2 <- imap_dfr(all.go, ~ mutate(.x, subtype = .y)) %>% mutate(subtype = factor(subtype, levels = levels(mg))) %>%
  group_by(subtype) %>% slice_min(n = 4, p.adjust, with_ties = F)

data <- data.frame(
  term = all.go2$Description,
  significance = -log10(all.go2$p.adjust),
  count = all.go2$Count,
  col = factor(all.go2$subtype, levels = levels(mg$subtype)),
  term2 = paste(all.go2$Description, all.go2$subtype, sep = "_")
)
data <- data[nrow(data):1,]

legend_col <- custom_colors$c3
names(legend_col) <- levels(mg$subtype)

ggplot(data, aes(x= fct_inorder(term2), y=significance)) +
  geom_segment(aes(xend=term2, y=0, yend= significance), color = rep(custom_colors$c3[3:1], each = 4), lwd = 1) +
  geom_point(aes(size= count, color = col)) + scale_color_manual(values = legend_col) +
  coord_flip() + theme_bw() + scale_x_discrete(labels = data$term, name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2))
ggsave("./images/supp/go_all_mg_style1.pdf", width = 12, height = 8)

ggplot(data, aes(x=factor(term2, levels = term2), y=significance, color = term)) +
  geom_segment(aes(xend=term2, y=0, yend=significance), color = rep(custom_colors$c3, each = 4), lwd = 1.5) +
  geom_point(aes(size= count), color = rep(custom_colors$c3, each = 4)) +
  coord_flip() + theme_classic() +  facet_grid(col~., scales = "free", space = "free") +
  scale_x_discrete(labels=setNames(data$term, data$term2), name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") +
  labs(y="-log10(p-value)", x="GO Term") + 
  theme(legend.position = "bottom", strip.text = element_text(color = "black", face = "bold", size = 7), 
        strip.background = element_rect(fill = "gray90"))
ggsave("./images/supp/go_all_mg_style2.pdf", width = 10, height = 11)


# Lipid -------------------------------------------------------------------

lipidset.name <- list.files(path = "./genesets", pattern = "^GOBP", full.names = T)
lipidset <- lapply(lipidset.name, function(x) unlist(strsplit(read.delim(x)[17,2], ",")))
lipidset <- lapply(lipidset, function(x) intersect(x, rownames(combined)))
names(lipidset) <- list.files(path = "./genesets", pattern = "^GOBP", full.names = F)

combined2 <- AddModuleScore(combined, features = lipidset, name = "lipid_")
Idents(combined2) <- "deg.compare"
lvls <- levels(combined2)
comparisons <- lapply(seq(1, length(lvls) - 1, by = 2), function(i) c(lvls[i], lvls[i+1]))
VlnPlot(combined2, grep("^lipid", colnames(combined2@meta.data), value = T), split.by ="group", pt.size = 0) &
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif", label.y = 0,
    hide.ns = T
  )


Idents(combined2) <- "celltype"
VlnPlot(combined2, "lipid_3", split.by ="group", pt.size = 0, cols = custom_colors$c4) + ggtitle(stringr::str_sub(names(lipidset)[3], start = 1 , end = -16))
ggsave("./images/supp/mg_lipid_biosyntesis.pdf", width =8, height = 4)

p1 <- VlnPlot(combined2, idents = "Mg", "lipid_1", split.by = "group", pt.size = 0, cols = custom_colors$c4[3:4]) + 
  labs(title = stringr::str_sub(names(lipidset)[1], start = 1 , end = -16)) + theme(plot.title = element_text(size = 6)) +
  stat_compare_means(label = "p.signif")
p2 <- VlnPlot(combined2, idents = "Mg", "lipid_2", split.by = "group", pt.size = 0, cols = custom_colors$c4[3:4]) + 
  labs(title =stringr::str_sub(names(lipidset)[2], start = 1 , end = -16)) + theme(plot.title = element_text(size = 6)) +
  stat_compare_means(label = "p.signif")
ggarrange(plotlist = list(p1,p2))
ggsave("./images/supp/mg_response_ldl.pdf", width = 7, height = 4)


Stacked_VlnPlot(combined, c("Acsl1", "Srebf1", "Fasn", "Acaca", "Elovl6", "Gpam", "Dgat1", "Dgat2", "Mlxipl"), split.by = "group", idents = "Mg", pt.size = 0)
ggsave("./images/supp/mg_denovo_lipid.pdf", width = 3, height = 9)



# HYPOMAP comparison ------------------------------------------------------
hypomap <- readRDS("HYPOMAP_human.rds")

subset.mg <- subset(hypomap, subset = C2 == "C2-49")
subset.mg$age_years <- as.integer(subset.mg$age_years)
subset.mg$age65 <- if_else(as.numeric(as.character(subset.mg$age_years)) > 65,"yes", "no")
subset.mg <- subset(subset.mg, subset = Donor_ID != "3u5kk")
subset.mg <- FindNeighbors(subset.mg, reduction = "scvi", nn.method = "annoy", compute.SNN = T, dims = 1:80)
subset.mg <- FindClusters(subset.mg, resolution = 0.5, algorithm = 4, graph.name = "RNA_snn", random.seed = 42)
subset.mg <- RunUMAP(subset.mg, dims = 1:80, reduction = "scvi")
DimPlot(subset.mg, reduction = "umap", group.by = "Donor_ID", cols = rep(custom_colors$c15,6), shuffle = T)
FeaturePlot_scCustom(subset.mg, "age_years")
DimPlot(subset.mg, cols = custom_colors$c15)
   

FeatureScatter_scCustom(subset(subset.mg, PLIN2 > 0 | CCL3 > 0), feature1 =  "PLIN2", feature2 = "CCL3", group.by = "age_years", colors_use = custom_colors$c15)

FeatureScatter_scCustom(mg, feature1 =  "Apoe", feature2 = "Plin2", group.by = "ident", colors_use = custom_colors$c15)

Nebulosa::plot_density(subset.mg, c("PLIN2", "CCL3", "APOE"), joint = T)       
FeaturePlot(subset.mg, c("PLIN2", "CCL3", "APOE"))
FeaturePlot_scCustom(subset.mg, c("age_years"))

expr_matrix <- GetAssayData(subset.mg, slot = "data") # "data" slot = log-normalized

dec <- correlatePairs(expr_matrix, pairings = list("CCL3", rownames(expr_matrix)))
dec <- dec[order(-dec$rho), ] # sort by correlation
head(dec)                        
as.data.frame(dec) %>% filter(gene2 == "PLIN2")

subset(subset.mg, subset = age65 == "yes") %>% Nebulosa::plot_density(c("PLIN2", "CCL3", "APOE"), joint = T)
ggsave("./images/main/hMg_old.pdf", width = 7.5, height = 6)

subset(subset.mg, subset = age65 == "no") %>% Nebulosa::plot_density(c("PLIN2", "CCL3", "APOE"), joint = T)  
ggsave("./images/main/hMg_young.pdf", width = 7.5, height = 6)




# cd8 ---------------------------------------------------------------------
cd8 <-  subset(combined, idents = "CD8")
cd8 <- FindVariableFeatures(cd8, nfeatures = 10000) %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
colnames(cd8@reductions$umap@cell.embeddings) <-  c("UMAP_1", "UMAP_2")

library(ProjecTILs)
# list_of_ref.maps <- get.reference.maps(collection = "mouse")
list_of_ref.maps <- get.reference.maps(reference = "TILs")
reference_map1 <- list_of_ref.maps[["mouse"]][["TILs"]]


proj.tils <- Run.ProjecTILs(cd8, ref = reference_map1, k = 20, split.by = "group", filter.cells = F)
p1 <- plot.projection(reference_map1, subset(proj.tils, subset = group == "young"), linesize = 0.5, pointsize = 0.5)
p2 <- plot.projection(reference_map1, subset(proj.tils, subset = group == "old"), linesize = 0.5, pointsize = 0.5)
p1 + p2
ggsave("./images/main/cd8_project_tils.pdf", width = 10, height = 6)

p1 <- plot.statepred.composition(reference_map1, subset(proj.tils, subset = group == "young"), metric = "Percent")
p2 <- plot.statepred.composition(reference_map1, subset(proj.tils, subset = group == "old"), metric = "Percent")
p1 + p2
ggsave("./images/main/cd8_proportion_tils.pdf", width = 10, height = 3)


# https://www.tandfonline.com/doi/full/10.1080/2162402X.2020.1737369#d1e806
plot.states.radar(reference_map1, list(young = subset(proj.tils, subset = group == "young"), old = subset(proj.tils, subset = group == "old")), 
                  genes4radar = c("Cd3e", "Cd3d", "Cd8a", "Cd4", "Foxp3", "Lef1", "Tcf7", "S1pr1", "Il7r", "Id3", "Tbx21", "Gzmk", 
                                  "Cd69", "Havcr2", "Pdcd1", "Lag3", "Mki67"), cols = custom_colors$c4[c(1,4)],
                  min.cells = 1, labels.col = "functional.cluster")
ggsave("./images/main/cd8_radar.pdf", width = 15, height = 15)

VlnPlot(cd8, c("Tbx21", "Gzmk", "Fasl", "Prf1", "Ifng"), split.by = "group", pt.size = 0, cols = custom_colors$c4[c(1,4)])
ggsave("./images/main/cd8_activation.pdf", width = 8, height = 6)
# https://www.sciencedirect.com/science/article/pii/S2211124720313176


tcr_old <- read.table("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/18-month__TCR/airr_rearrangement.tsv", sep = "\t", header = T) 
tcr_young <- read.table("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/11-week__TCR/airr_rearrangement.tsv", sep = "\t", header = T) 
cdr_old <- read.table("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/18-month__TCR/filtered_contig_annotations.csv", sep = ",", header = T) 
cdr_young <- read.table("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/11-week__TCR/filtered_contig_annotations.csv", sep = ",", header = T) 
clone_old <- read.table("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/18-month__TCR/clonotypes.csv", sep = ",", header = T) 
clone_young <- read.table("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/11-week__TCR/clonotypes.csv", sep = ",", header = T) 


process_trb <- function(df, condition_label) {
  df %>%
    filter(chain == "TRB", high_confidence == "true", productive == "true") %>%
    transmute(
      CDR3b   = cdr3,
      TRBV    = v_gene,
      TRBJ    = j_gene,
      CDR3a   = "NA",
      clonotype = paste0(raw_clonotype_id, ":", condition_label),
      count   = 1
    )
}

out <- bind_rows(
  process_trb(cdr_old,  "old"),
  process_trb(cdr_young,"young")
)

write.table(out, file = "trb_only_old_young.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#clonotype young
classified <- clone_young %>%
  mutate(
    category = case_when(
      (!is.na(inkt_evidence) & inkt_evidence != "") &
        (!is.na(mait_evidence) & mait_evidence != "") ~ "Dual iNKT/MAIT",
      !is.na(inkt_evidence) & inkt_evidence != "" ~ "iNKT",
      !is.na(mait_evidence) & mait_evidence != "" ~ "MAIT",
      TRUE ~ "Conventional"
    )
  )
use_prop <- "proportion" %in% names(classified) && any(!is.na(classified$proportion))

summary_cat <- classified %>%
  group_by(category) %>%
  summarise(
    weight = if (use_prop) sum(proportion, na.rm = TRUE) else sum(frequency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    frac = weight / sum(weight),
    label = paste0(category, " (", percent(frac, accuracy = 0.1), ")")
  ) %>%
  arrange(desc(frac))

summary_cat$category <- fct_inorder(summary_cat$category)

p1 <- ggplot(summary_cat, aes(x = "", y = frac, fill = category)) +
  geom_col(width = 1, color = "grey20", alpha = 0.5) +
  coord_polar(theta = "y") +
  geom_text(aes(label = ifelse(frac >= 0.03, percent(frac, accuracy = 0.1), "")),
            position = position_stack(vjust = 0.5), size = 4) +
  labs(title = "Clonotype Categories", fill = "Category") +
  theme_void(base_size = 12) + scale_fill_cosmic()

#clonotype old
classified <- clone_old %>%
  mutate(
    category = case_when(
      (!is.na(inkt_evidence) & inkt_evidence != "") &
        (!is.na(mait_evidence) & mait_evidence != "") ~ "Dual iNKT/MAIT",
      !is.na(inkt_evidence) & inkt_evidence != "" ~ "iNKT",
      !is.na(mait_evidence) & mait_evidence != "" ~ "MAIT",
      TRUE ~ "Conventional"
    )
  )
use_prop <- "proportion" %in% names(classified) && any(!is.na(classified$proportion))

summary_cat <- classified %>%
  group_by(category) %>%
  summarise(
    weight = if (use_prop) sum(proportion, na.rm = TRUE) else sum(frequency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    frac = weight / sum(weight),
    label = paste0(category, " (", percent(frac, accuracy = 0.1), ")")
  ) %>%
  arrange(desc(frac))

summary_cat$category <- fct_inorder(summary_cat$category)

p2 <- ggplot(summary_cat, aes(x = "", y = frac, fill = category)) +
  geom_col(width = 1, color = "grey20", alpha = 0.5) +
  coord_polar(theta = "y") +
  geom_text(aes(label = ifelse(frac >= 0.05, percent(frac, accuracy = 0.1), "")),
            position = position_stack(vjust = 0.5), size = 4) +
  labs(title = "Clonotype Categories", fill = "Category") +
  theme_void(base_size = 12) + scale_fill_cosmic()

p1 + p2

ggsave("./images/main/cd8_mait_prop.pdf", width = 7, height = 4)


#cloniotype freq
summarize_clone_freq <- function(df, label) {
  df %>%
    mutate(
      freq_bin = case_when(
        frequency > 2 ~ ">2",
        frequency == 2 ~ "=2",
        frequency == 1 ~ "=1",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(freq_bin) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      dataset = label,
      prop = n / sum(n),
      group = paste0(label, "_", freq_bin)  # unique combo
    )
}

old_summary   <- summarize_clone_freq(clone_old, "Old")
young_summary <- summarize_clone_freq(clone_young, "Young")

summary_all <- bind_rows(young_summary, old_summary)

pal_old   <- colorRampPalette(c("#e3a1bb", "#d11778"))(3)
pal_young <- colorRampPalette(c("#a6dbe2", "#1289A7"))(3)

names(pal_old)   <- paste0("Old_",   c("=1", "=2", ">2"))
names(pal_young) <- paste0("Young_", c("=1", "=2", ">2"))

mycols <- c(pal_old, pal_young)

summary_all <- summary_all %>%
  mutate(pct_label = percent(prop, accuracy = 0.1))

summary_all$dataset <- fct_inorder(summary_all$dataset)

# plot
ggplot(summary_all, aes(x = "", y = prop, fill = group)) +
  geom_col(width = 1, color = "gray20") +
  coord_polar(theta = "y") +
  facet_wrap(~dataset) +
  scale_fill_manual(
    values = mycols,
    labels = c("=1", "=2", ">2"),
    breaks = c("Old_=1", "Old_=2", "Old_>2"),
    guide = guide_legend(title = "Frequency bin")
  ) +
  geom_text(
    aes(label = ifelse(prop > 0.05, pct_label, "")),
    position = position_stack(vjust = 0.5),
    color = "gray10", size = 3
  ) +
  labs(title = "Clonotype Frequency Distribution") +
  theme_void(base_size = 12)

ggsave("./images/main/cd8_clonotype_prop.pdf", width = 7, height = 4)


# scRepertoire ------------------------------------------------------------

contig_list <- list(read.csv("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/11-week__TCR/filtered_contig_annotations.csv"), 
                    read.csv("/disk1/aging_scRNAseq2/HN00190425/HN00190425_result_10X/18-month__TCR/filtered_contig_annotations.csv"))

combined.TCR <- combineTCR(contig_list, 
                           samples = c("young", "old"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

clonalScatter(combined.TCR, cloneCall ="aa", x.axis = "young", y.axis = "old",
              dot.size = "total", graph = "proportion", chain = "both")
ggsave("./images/supp/cd8_overlap.pdf", width = 5, height = 3.5)


percentAA(combined.TCR, chain = "TRA", aa.length = 20) + percentAA(combined.TCR, chain = "TRB", aa.length = 20)


clonalProportion(combined.TCR, 
                 cloneCall = "aa",
                 clonalSplit = c(5, 10, 20, 50, 100)) 
ggsave("./images/main/cd8_clonal_domination.pdf", width = 3, height = 3)


