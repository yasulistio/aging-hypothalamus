{library(DESeq2)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(ggsci)
  library(viridis)
  library(pheatmap)
  library(BayesPrism)
  library(tibble)
  library(cowplot)
  library(unikn)}


cts <- read.delim("Expression_Profile.mm10.gene.txt", row.names = 3)
cts <- cts[,c(11:13,8:10,14:16)] #selecting column with readcounts
colnames(cts) <- gsub("_Read_Count","" ,colnames(cts))
colnames(cts) <- gsub("M.H", "H", colnames(cts))

# QC ----------------------------------------------------------------------
cts.m <- reshape2::melt(cts)
cts.m$value <- log2(cts.m$value+1)
ggdensity(cts.m, x='value', fill='variable', alpha = 0.5) + ylim(c(0,0.07))
ggsave("histogram.pdf")


a <-  prcomp(cts)
b <-  summary(a)
coldata <- data.frame(sample = colnames(cts), age = as.factor(rep(c("2mo","12mo","26mo"), each =3)))
ggplot(data= as.data.frame(a$rotation), aes(x=PC1, y=PC2, col= coldata$age )) + geom_point(size = 2.5) + 
  geom_text_repel(aes(label= colnames(cts)), box.padding = 0.2) + theme_bw() + theme(legend.position = "none") + 
  scale_color_jama() + labs(x= paste0("PC1 (",b$importance[2]*100,"%)"), y= paste0("PC2 (",b$importance[5]*100,"%)"))

ggsave("./images/pca.pdf", scale = 1)


# DESEq2 ------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(cts, coldata, ~age)
dds <- DESeq(dds)
res <- results(dds, contrast = c('age', '26mo', '2mo'))
res.sig <- na.omit(res) %>% .[.$padj < 0.01,]
View(data.frame(res.sig) %>% arrange(padj))
deg <- rownames(res.sig)

as.data.frame(res.sig) %>% filter(log2FoldChange < -1) %>% nrow
as.data.frame(res.sig) %>% filter(log2FoldChange > 1) %>% nrow

vsd <- vst(dds)
cts.norm <- assay(vsd)

ann_data <- data.frame(Age= rep(c("2mo", "12mo", "26mo"), each =3))
rownames(ann_data) <- colnames(cts.norm)
ann_col <- list(Age = c("2mo" =  pal_karpfenblau[[1]], "12mo" = pal_karpfenblau[[3]], "26mo" = pal_karpfenblau[[5]]))
pheatmap(cts.norm[deg,], scale = "row", cluster_cols = F, show_rownames = F, color = colorRampPalette(viridis(10, option = "D", direction = 1))(10),
         annotation_col = ann_data, annotation_colors = ann_col, show_colnames = F, treeheight_row = 0, annotation_names_col = F, cellwidth = 14, cellheight = 0.4,
         breaks = seq(-2,2, length.out = 11))


res26v2 <- results(dds, contrast = c('age', '26mo', '2mo')) %>% na.omit() 
res12v2 <- results(dds, contrast = c('age', '12mo', '2mo')) %>% na.omit() 
res26v12 <- results(dds, contrast = c('age', '26mo', '12mo')) %>% na.omit() 
res.list <- list(res26v2, res12v2, res26v12)
names(res.list) <- c("2mo vs 26mo", "2mo vs 12mo", "12mo vs 26mo")


openxlsx::write.xlsx(lapply(res.list, function(x) as.data.frame(x) %>% filter(padj < 0.05) %>% arrange(padj)), "bulk_DEG.xlsx", rowNames = T)


plotlist <- lapply(res.list, function(x) {
  res <- as.data.frame(x)
  col <-ifelse(res$log2FoldChange > 1 & res$padj <0.05, "Up", 
               ifelse(res$log2FoldChange < -1 & res$padj <0.05, "Down", "Non"))
  stroke <-ifelse(res$log2FoldChange > 1 & res$padj <0.05, 0.1, 
                  ifelse(res$log2FoldChange < -1 & res$padj <0.05, 0.1, NA))
  z <- as.data.frame(res[c("C4b", "Ccl3", "Oasl2", "Lag3", "Clec7a"),])
  p <- ggplot(data= as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj), fill = col, stroke = stroke)) + geom_point(shape = 21, size =2) +
    geom_vline(xintercept = c(-1,1), color="grey10", linewidth=0.2, lty=2, alpha = 0.5)+ 
    geom_hline(yintercept = -log10(0.05), color="grey10", linewidth=0.2, lty=2, alpha = 0.5) +
    geom_label_repel(data=z, x= z$log2FoldChange, y=-log10(z$padj), aes(label=rownames(z)), inherit.aes = F, box.padding = unit(1, "lines"),
                     point.padding = unit(0.5, "lines"), max.overlaps = Inf, force = 5 ) + xlim(-15,15) + ylim(0,40) + 
    labs(title = stringr::str_sub(x@elementMetadata$description[4], 20,-1)) +
    theme_classic(base_size = 15) + theme(legend.position = "none", font.size = 5) + scale_fill_manual(values= c("#2d9eb3", "grey80", "#E64B36FF")) + 
    ylab("-log10(FDR)")
  return(p)})

ggarrange(plotlist = plotlist)
ggsave("./images/main/1A_volcano.pdf", width = 7, height = 7)  

# clusterprofiler ---------------------------------------------------------

res.sig <- res26v2 %>% as.data.frame() %>% filter(padj <0.05)

go.list <- lapply(list("BP", "MF", "CC"), function(cat){ 
  x <- enrichGO(gene          = rownames(res.sig %>% filter(log2FoldChange > 1)),
                universe      = rownames(res),
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = cat,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = F)
  go <- clusterProfiler::simplify(x, cutoff = 0.5, measure = "Wang")
  return(go)
})

go.result.list <- lapply(1:3, function(x) go.list[[x]]@result)
names(go.result.list) <- c("BP", "MF", "CC")
openxlsx::write.xlsx(go.result.list, "bulk_GO.xlsx")

all.go <- do.call(rbind,lapply(go.list, function(x) head(x,5)))
all.go <- all.go %>% mutate(type = factor(rep(c("BP", "MF", "CC"), each =5), levels = c("BP", "MF", "CC")))

data <- data.frame(
  term = stringr::str_trunc(rev(all.go$Description), width = 50),
  significance = rev(-log10(all.go$p.adjust)),
  count = rev(all.go$Count),
  term2 = rev(paste(all.go$Description, all.go$celltype, sep = "_")),
  type = rev(factor(all.go$type, levels = levels(all.go$type)))
)

ggplot(data, aes(x=factor(term2, levels = term2), y=significance)) +
  geom_segment(aes(xend=term2, y=0, yend=significance, color = type), lwd = 1) +
  geom_point(aes(size= count, color = type)) +
  coord_flip() + theme_bw() + scale_x_discrete(labels = data$term, name = "GO BP") + scale_y_continuous(name = "-log10(FDR)") +
  labs(y="-log10(p-value)", x="GO Term") + theme(legend.position = "right") + guides(color = guide_legend(ncol = 2)) + scale_color_npg()

ggsave("./images/main/bulk_DEG_up.pdf", width = 6.5, height = 5)



# graphs ------------------------------------------------------------------


make_cpm_graph <- function(data, genes) {
  df <- as.data.frame(data[genes,]) %>% rownames_to_column("gene") %>% pivot_longer(cols = -gene)
  df$name <- stringr::str_sub(df$name, 1, -3)
  df$name <- stringr::str_sub(df$name, 6, -1)
  df$name <- fct_inorder(df$name)
  
  if (is.null(genes)) {
    genes <- rownames(data)
  }
  ggplot(df, aes(x = name, y = value, col = name, fill = name)) +
    geom_point(size = 1)+
    stat_summary(fun = base::mean, 
                 geom = "bar", 
                 position = position_dodge(0.9),
                 alpha = 0.2) +
    stat_summary(fun.data = mean_sdl, 
                 fun.args = list(mult = 1), 
                 geom = "errorbar", 
                 position = position_dodge(0.9), 
                 width = 0.2) +
    labs(y = "CPM (mean ± SD)", x = "Timepoint", fill = "Gene") +
    theme_cowplot() + facet_wrap(~gene, scales = "free")}


cts.cpm <- edgeR::cpm(cts)
make_cpm_graph(cts.cpm, c("Ccl3","Ccl4"))
ggsave("./images/main/ccl34_cpm.pdf", width = 6, height = 4)  


df <- cts.cpm[as.data.frame(res.sig) %>% filter(log2FoldChange >1) %>% rownames(),]
df <- as.data.frame(df) %>% bind_rows(colSums(.)/nrow(.))
make_cpm_graph(df, rownames(df)[nrow(df)])
ggsave("./images/main/average_cpm.pdf", width = 3, height = 4) 


deg.names <- lapply(res.list, function(x){
  up <- as.data.frame(x) %>% filter(padj < 0.05, log2FoldChange > 1) %>% rownames()
  down <- as.data.frame(x) %>% filter(padj < 0.05, log2FoldChange < -1) %>% rownames()
  return(list(up = up, down = down))})



# Venn --------------------------------------------------------------------


library(ggVennDiagram)
gene_sets <- list(
  "2vs26" = deg.names$`2mo vs 26mo`$up,
  "2vs12" = deg.names$`2mo vs 12mo`$up,
  "12vs26" = deg.names$`12mo vs 26mo`$up
)

ggVennDiagram(gene_sets) + ggplot2::theme_void()
ggsave("./images/supp/venn_up.pdf", width = 4, height = 4)
lapply(res.list, function(x) as.data.frame(x) %>% filter(padj < 0.01, log2) %>% nrow())

gene_sets <- list(
  "2vs26" = deg.names$`2mo vs 26mo`$down,
  "2vs12" = deg.names$`2mo vs 12mo`$down,
  "12vs26" = deg.names$`12mo vs 26mo`$down
)
ggVennDiagram(gene_sets) + ggplot2::theme_void()
ggsave("./images/supp/venn_down.pdf", width = 4, height = 4)


p1 <- make_cpm_graph(cts.cpm, c("Serpina3n", "Serpina3c", "H2-K1", "H2-Q4", "H2-Q7", "H2-K2", "Clec7a", "Cst7", "Tnf")) + 
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p2 <- make_cpm_graph(cts.cpm, c("Ifi44", "Ifi207", "Isg15", "Oas1a", "Oas2", "Lag3", "Lilrb4a", "Mmp12", "Stat1" )) + 
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p3 <- make_cpm_graph(cts.cpm, c("C1qa", "C1qb", "C3", "C4a", "C4b", "Ly9", "Lyz2")) + 
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggarrange(plotlist = list(p1,p2,p3))
ggsave("./images/supp/bulkseq_gene_trends.pdf", width = 8, height = 11.5)
