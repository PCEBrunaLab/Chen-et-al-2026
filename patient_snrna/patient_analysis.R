## HB patients: snRNA-seq alone

## /////////////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

setwd('/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb_patients/snrna')

renv::load()

library(Seurat)
library(scran)
library(scater)
library(scuttle)
library(clusterProfiler)
library(msigdbr)
library(tidyverse)
library(patchwork)
library(reshape2)
library(pheatmap)

## [ Plot themes ] ----

## Custom ggplot theme
umap.theme <- function() {
  require(ggplot2)
  return(theme_bw() +
           theme(aspect.ratio = 1, 
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 axis.line = element_line(colour = "#16161D", linewidth = 0.8),
                 axis.ticks = element_line(colour = "#16161D", linewidth = 0.8),
                 legend.title = element_blank()))
}

condition.cols <- c("P1_F1" = "#F7C9E3",
                    "P1_F2" = "#E58ABD",
                    "P2_F1" = "#FFE0B8",
                    "P2_F2" = "#FFAA55")

## [ QC filtering ] ----

## Load in data
patient.sce <- readRDS("data/HB_patients_SCE-norm_relaxed.RDS")
length(unique(colnames(patient.sce))) ## 75897 cells

## Compute logcounts
patient.sce <- logNormCounts(patient.sce)

## ENSEMBL database
# library(EnsDb.Hsapiens.v86)
# library(biomaRt)
# ens.86.genes <- genes(EnsDb.Hsapiens.v86)
# linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype == "lincRNA"]
# ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
# ens.bm <- getBM(attributes = c("ensembl_gene_id",
#                                "chromosome_name",
#                                "start_position",
#                                "end_position",
#                                "external_gene_name",
#                                "hgnc_symbol",
#                                "description"),
#                 mart = ensembl,
#                 filters = "ensembl_gene_id",
#                 values = rownames(patient.sce))
# # Save ensembl.csv for reuse
# write.csv(ens.bm, "data/patient_ensembl_biomart.csv")
# gc()
ens.bm <- read.csv("data/patient_ensembl_biomart.csv", row.names = NULL)
ens.bm$X <- NULL

## Identify duplicated gene names
dup.ensg <- ens.bm$ensembl_gene_id[duplicated(ens.bm$ensembl_gene_id)] 
dup.ensg <- setNames(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% dup.ensg],
                     ens.bm$ensembl_gene_id[ens.bm$ensembl_gene_id %in% dup.ensg])
## Duplicated gene names: LINC00595

## Created filtered ens.bm and keep only useful information
ens.filt.bm <- ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)
## This filters out around 1000 genes
nrow(ens.filt.bm) ## 35456
nrow(patient.sce) ## 36601
ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]
## Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name == "", ]
## Around 10000 genes filtered out
nrow(ens.filt.bm) ## 25577
nrow(patient.sce) ## 36601
## Keep track of everything that didn't get annotated
ens.removed.bm <- ens.bm[!ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm) ## 9879
length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
## 9733 are novel transcripts or novel transcripts antisense to a gene
grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
## The rest are mostly either novel, mitochondrially encoded, or just empty description
## Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])
## "ELFN2"  "GOLGA8M" "LINC00595" "LINC01115" "LINC03021" "LINC03023" "LINC03025" "RAET1E-AS1"  "SPATA13"  

## Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm) ## 25568

## Identify mitochondrial, ribosomal and linc genes
mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]
linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]

## Add QC metrics to SCE object
patient.sce <- addPerCellQCMetrics(patient.sce, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))
patient.sce <- patient.sce[ens.filt.bm$ensembl_gene_id,]
rownames(patient.sce) <- ens.filt.bm$external_gene_name
rownames(ens.filt.bm) <- ens.filt.bm$external_gene_name

## Assign rowData to SCE object with filtered ens.bm
rowData(patient.sce) <- ens.filt.bm
saveRDS(patient.sce, "data/patient_sce.rds")

## Extract meta data from SCE object
patient.meta <- as.data.frame(colData(patient.sce))
patient.meta$batch <- NULL
## Generate Seurat object from SCE object
patient.seurat <- CreateSeuratObject(counts = assay(patient.sce, "counts"),
                                     assay = "RNA",
                                     meta.data = patient.meta)
saveRDS(patient.seurat, "data/patient_seurat.rds")

patient.sce$Condition <- paste(patient.sce$Patient_ID, patient.sce$Fraction, sep = "_")
patient.sce$Condition_Rep <- paste(patient.sce$Condition, patient.sce$Replicate, sep = "_")
patient.seurat$Condition <- paste(patient.seurat$Patient_ID, patient.seurat$Fraction, sep = "_")
patient.seurat$Condition_Rep <- paste(patient.seurat$Condition, patient.seurat$Replicate, sep = "_")
## Plot QC metrics
dir.create("plots/qc", recursive = TRUE)
det.sce.gg <- plotColData(patient.sce, y = "detected", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
umi.sce.gg <- plotColData(patient.sce, y = "total", x = "Condition_Rep", colour_by = "Condition_Rep") + 
  labs(x = element_blank(), y = "UMI / cell", title = "Detected UMI per cell") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(sides = "l", outside = TRUE) + coord_cartesian(clip = "off")
mt.sce.gg <- plotColData(patient.sce, y = "subsets_mt_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
rb.sce.gg <- plotColData(patient.sce, y = "subsets_ribo_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
det.sce.gg + umi.sce.gg + mt.sce.gg + rb.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.8,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/qc/patient_qc_plots.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/patient_qc_plots.png", width = 8.3, height = 5.8)

## Examine whether there are any genes that dominate expression in cells - large % of expression occupied by single gene
counts.assay <- counts(patient.sce)
counts.assay@x <- counts.assay@x/rep.int(colSums(counts.assay), diff(counts.assay@p))
top.expr <- order(Matrix::rowSums(counts.assay), decreasing = TRUE)[20:1]
top.expr <- as.matrix(t(counts.assay[top.expr, ]))
top.expr <- as.data.frame(top.expr)
ord.idx <- colnames(top.expr)
top.expr <- melt(top.expr, variable.name = "Gene", value.name = "Prop")
top.expr$Gene <- factor(top.expr$Gene, levels = ord.idx)

ggplot(top.expr, mapping = aes(x = Gene, y = Prop * 100)) +
  geom_boxplot() +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_line(linetype = "dotted")) +
  coord_flip() +
  labs(y = "% total expression", x = "Gene", title = "Percentage expression of a single gene per total cell expression")
ggsave("plots/qc/patient_qc_top_genes.pdf", width = 8.3, height = 5.8)
ggsave("plots/qc/patient_qc_top_genes.png", width = 8.3, height = 5.8)

rm(counts.assay)
gc()

## Visualise MALAT1 expression and lincRNA per sample
malat1.sce.gg <- plotExpression(patient.sce, 
                                features = "MALAT1", 
                                exprs_values = "logcounts",
                                x = "Condition_Rep", 
                                colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Log2 Expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
lincrna.sce.gg <- plotColData(patient.sce, y = "subsets_linc_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "lincRNA\nexpression % of total") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)
malat1.sce.gg + lincrna.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.8,
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/qc/patient_qc_malat1_lincrna.pdf", width = 5.8, height = 8.3)
ggsave("plots/qc/patient_qc_malat1_lincrna.png", width = 5.8, height = 8.3)

## Cell cycle scoring
patient.seurat <- NormalizeData(patient.seurat)
patient.seurat <- CellCycleScoring(patient.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
patient.sce$Seurat.Phase <- patient.seurat$Phase
patient.sce$Seurat.S <- patient.seurat$S.Score
patient.sce$Seurat.G2M <- patient.seurat$G2M.Score
saveRDS(patient.sce, "data/patient_sce_cell_cycle.rds")
saveRDS(patient.seurat, "data/patient_seurat_cell_cycle.rds")

# patient.sce <- readRDS("data/patient_sce_cell_cycle.rds")
# patient.seurat <- readRDS("data/patient_seurat_cell_cycle.rds")
## Check gene and cell meta data
colnames(colData(patient.sce))
colnames(rowData(patient.sce))
table(patient.sce$Condition)
table(patient.sce$Condition_Rep)

det.sce.gg <- plotColData(patient.sce, y = "detected", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.sce$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 1000, linetype = "dashed")
# geom_hline(yintercept = 200, linetype = "dotted")
mt.sce.gg <- plotColData(patient.sce, y = "subsets_mt_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.sce$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_hline(yintercept = 15, linetype = "dotted")
rb.sce.gg <- plotColData(patient.sce, y = "subsets_ribo_percent", x = "Condition_Rep", colour_by = "Condition_Rep") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.sce$Condition_Rep)))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed")
det.sce.gg + mt.sce.gg + rb.sce.gg & theme(aspect.ratio = 0.6,
                                           axis.text = element_text(size = 12),
                                           axis.text.x = element_text(size = 12))
ggsave("plots/qc/patient_qc_plots_filtering.pdf", width = 11.7, height = 5.8)
ggsave("plots/qc/patient_qc_plots_filtering.png", width = 11.7, height = 5.8)

## No need for expression filtering
detected.filt <- colnames(patient.sce)[patient.sce$detected > 1000] ## No cells removed
## Genes have to have at least 5 reads 
rowcounts.filt <- rownames(patient.sce)[Matrix::rowSums(counts(patient.sce)) > 5] ## 2375 genes removed
## Mitochondrial genes filtering
mt.genes.filt <- grep("^MT-", rownames(patient.sce), invert = TRUE, value = TRUE) ## Doesn't remove any genes
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)
patient.filt.sce <- patient.sce[genes.filt, detected.filt]
dim(patient.sce) - dim(patient.filt.sce)
## We lose 2375 genes and 0 cells

## Gene expression of each cell has to have less than 15% mitochondrial
mito.filt <- colnames(patient.filt.sce)[patient.filt.sce$subsets_mt_percent < 15]
## No ribosomal filtering - too harsh for snRNA-seq
#ribo.filt <- colnames(patient.filt.sce)[patient.filt.sce$subsets_ribo_percent > 5]
#qc.filt <- intersect(mito.filt, ribo.filt)

patient.filt.sce <- patient.filt.sce[, mito.filt]
dim(patient.sce) - dim(patient.filt.sce)
## We lose 195 cells in total

## Percentage of the library that we filter out
round(100 - (100 * (table(patient.filt.sce$Condition_Rep) / table(patient.sce$Condition_Rep))), 1)
table(patient.filt.sce$Condition_Rep)

## We can do the same filtering in Seurat
rowcounts.filt <- rownames(patient.seurat)[Matrix::rowSums(GetAssayData(patient.seurat, slot = "counts", assay = "RNA")) > 5]
mt.genes.filt <- grep("^MT-", rownames(patient.seurat), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

patient.filt.seurat <- subset(patient.seurat, 
                              cells = WhichCells(patient.seurat, 
                                                 expression = subsets_mt_percent < 15 &
                                                   detected > 1000),
                              features = genes.filt)

## The dimensions are equal
dim(patient.filt.seurat) == dim(patient.filt.sce)

saveRDS(patient.filt.sce, "data/patient_sce_filtered.rds")
saveRDS(patient.filt.seurat, "data/patient_seurat_filtered.rds")

## [ Filter HVGs ] ----

patient.filt.seurat <- readRDS("data/patient_seurat_filtered.rds")

## Yura's script to filter HVGs
patient.filt.seurat <- NormalizeData(patient.filt.seurat)
patient.filt.seurat <- FindVariableFeatures(patient.filt.seurat, selection.method = "vst", nfeatures = 1000)
LabelPoints(plot = VariableFeaturePlot(patient.filt.seurat, assay = "RNA"),
            points = head(VariableFeatures(patient.filt.seurat), 20), repel = TRUE)
ggsave("plots/filter_hvgs/patient_filter_hvgs.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/patient_filter_hvgs.png", width = 8.3, height = 5.8)

## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(patient.filt.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(patient.filt.seurat), value = TRUE),
              grep("[\\.]", rownames(patient.filt.seurat), value = TRUE),
              grep("^LINC", rownames(patient.filt.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
VariableFeatures(patient.filt.seurat) <- VariableFeatures(patient.filt.seurat)[which(!VariableFeatures(patient.filt.seurat) %in% rm.genes)]
## Scale data and score cell cycle
patient.filt.seurat <- ScaleData(patient.filt.seurat, features = VariableFeatures(patient.filt.seurat))

tmp.seurat <- CellCycleScoring(object = patient.filt.seurat, 
                               g2m.features = cc.genes$g2m.genes, 
                               s.features = cc.genes$s.genes)
tmp.seurat$Cycle.Score <- tmp.seurat$S.Score - tmp.seurat$G2M.Score

seurat.cycle.melt <- ggsankey::make_long(
  data.frame("Sample" = patient.filt.seurat$Sample, 
             "SCTransform" = tmp.seurat$Phase,
             "Normal" = patient.filt.seurat$Phase),
  SCTransform, Normal)

library(ggsankey)
ggplot(seurat.cycle.melt, 
       aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = factor(node),
           label = node)) +
  geom_sankey(flow.alpha = 0.7) +
  geom_sankey_label() +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none") +
  labs(x = element_blank(), title = "Seurat Cell Cycle Prediction") +
  scale_fill_manual(values = hues::iwanthue(3))
## Not a huge difference
rm(tmp.seurat)
gc()

DefaultAssay(patient.filt.seurat) <- "RNA"
patient.filt.seurat$Seurat.Cycle.Score <- patient.filt.seurat$S.Score - patient.filt.seurat$G2M.Score
saveRDS(patient.filt.seurat, "data/patient_seurat_hvgs.rds")

# patient.seurat <- readRDS("data/patient_seurat_hvgs.rds")
patient.seurat <- RunPCA(patient.seurat, verbose = FALSE, npcs = 100, 
                         features = VariableFeatures(patient.seurat))
ElbowPlot(object = patient.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate two separate UMAP dim reds +/-
harmony.seurat <- harmony::RunHarmony(patient.seurat, group.by.vars = "Condition_Rep",
                                      theta = 0.5, lambda = 1, sigma = 0.05,
                                      assay.use = "RNA", reduction = "pca",
                                      dims.use = 1:30, reduction.save = "Harmony",
                                      max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

umap.gg <- DimPlot(harmony.seurat, group.by = "Condition_Rep", order = TRUE) +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient.seurat$Condition_Rep)))) +
  umap.theme() + labs(title = "Conditions")
umap.gg
ggsave("plots/filter_hvgs/umap_harmony_t05_l1_s005.pdf", width = 8.3, height = 5.8)
ggsave("plots/filter_hvgs/umap_harmony_t05_l1_s005.png", width = 8.3, height = 5.8)

# patient.seurat <- RunUMAP(patient.seurat, reduction = "pca", dims = 1:30)
# patient.seurat <- FindNeighbors(patient.seurat, reduction = "pca", dims = 1:30)

cluster.search <- function(seurat, from = 0.2, to = 0.8, by = 0.2){
  tmp.res <- lapply(seq(from = from, to = to, by = by), function(x){
    tmp.clust <- FetchData(FindClusters(seurat, verbose = TRUE, resolution = x), 
                           c(paste0("RNA_snn_res.", x), "seurat_clusters"))
    colnames(tmp.clust) <- c(paste0("snn_res.", x), paste0("seurat_clusters.", x))
    return(tmp.clust)
  })
  tmp.res <- do.call("cbind", tmp.res)
  return(AddMetaData(seurat, metadata = tmp.res))
}
harmony.seurat <- cluster.search(harmony.seurat)
saveRDS(harmony.seurat, "data/patient_seurat_hvgs_harmony.rds")

patient.seurat <- readRDS("data/patient_seurat_hvgs_harmony.rds")

## [ Jaccard similarity ] ----

## Split data
patient.seurat <- readRDS("data/patient_seurat_hvgs.rds")
obj.list <- SplitObject(patient.seurat, split.by = "Patient_ID")

## Rerun variable feature selection
obj.list <- lapply(obj.list, FindVariableFeatures)
## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(patient.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(patient.seurat), value = TRUE),
              grep("[\\.]", rownames(patient.seurat), value = TRUE),
              grep("^LINC", rownames(patient.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
output <- list()
for (i in 1:length(obj.list)) {
  obj <- obj.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(obj.list)
rm(patient.seurat, obj.list, obj)

## Scale data
output <- lapply(output, function(obj) {
  ScaleData(obj, features = VariableFeatures(obj))
})
## Run PCA
output <- lapply(output, function(obj) {
  RunPCA(obj, verbose = FALSE, npcs = 100,
         features = VariableFeatures(obj))
})

## Run UMAP and compute nearest neighbours 
output <- lapply(output, RunUMAP, reduction = "pca", dims = 1:30)
output <- lapply(output, FindNeighbors, reduction = "pca", dims = 1:30)
## Clustering 
cluster.search <- function(seurat, from = 0.2, to = 0.8, by = 0.2){
  tmp.res <- lapply(seq(from = from, to = to, by = by), function(x){
    tmp.clust <- FetchData(FindClusters(seurat, verbose = TRUE, resolution = x), 
                           c(paste0("RNA_snn_res.", x), "seurat_clusters"))
    colnames(tmp.clust) <- c(paste0("snn_res.", x), paste0("seurat_clusters.", x))
    return(tmp.clust)
  })
  tmp.res <- do.call("cbind", tmp.res)
  return(AddMetaData(seurat, metadata = tmp.res))
}
output <- lapply(output, cluster.search)
saveRDS(output, "data/patient_seurat_split_id_list.rds")

patient.list <- readRDS("data/patient_seurat_split_id_list.rds")
patient.seurat <- readRDS("data/patient_seurat_hvgs_harmony.rds")
## Look through UMAPs and choose best cluster resolution
# DimPlot(patient.seurat, group.by = "seurat_clusters.0.4", order = TRUE) +
#   umap.theme()
Idents(patient.list$P1) <- patient.list$P1$seurat_clusters.0.2
Idents(patient.list$P2) <- patient.list$P2$seurat_clusters.0.4
Idents(patient.seurat) <- patient.seurat$seurat_clusters.0.4

## Get cluster markers - MAST takes longer to run
library(future)
parallel::detectCores() ## 10
plan(multisession, workers = 8) ## 6-8 recommended on M1 Pro
options(future.globals.maxSize = 8 * 1024^3) ## 16 GB RAM
p1.markers <- FindAllMarkers(patient.list$P1, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
p2.markers <- FindAllMarkers(patient.list$P2, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
all.markers <- FindAllMarkers(patient.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)

## Filter significant and higher expressed markers
p1.markers <- p1.markers[p1.markers$p_val_adj < 0.05, ]
p2.markers <- p2.markers[p2.markers$p_val_adj < 0.05, ]
p3.markers <- p3.markers[p3.markers$p_val_adj < 0.05, ]
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]

p1.markers <- p1.markers[p1.markers$avg_log2FC > 0.5, ]
p2.markers <- p2.markers[p2.markers$avg_log2FC > 0.5, ]
p3.markers <- p3.markers[p3.markers$avg_log2FC > 0.5, ]
all.markers <- all.markers[all.markers$avg_log2FC > 0.5, ]

## Split into lists of genes per cluster
p1.markers <- split(p1.markers, p1.markers$cluster)
p2.markers <- split(p2.markers, p2.markers$cluster)
p3.markers <- split(p3.markers, p3.markers$cluster)
all.markers <- split(all.markers, all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
## Yura uses 500 per cluster
lapply(p1.markers, nrow) ## 0 = 1009, 1 = 601, 2 = 1559, 3 = 443, 4 = 1117, 5 = 1158, 6 = 635, 7 = 1009, 8 = 1416
lapply(p2.markers, nrow) ## 0 = 214, 1 = 538, 2 = 631, 3 = 502, 4 = 304, 5 = 523, 6 = 292, 7 = 1336, 8 = 632, 9 = 1244
lapply(p3.markers, nrow) ## 0 = 214, 1 = 538, 2 = 631, 3 = 502, 4 = 304, 5 = 523, 6 = 292, 7 = 1336, 8 = 632, 9 = 1244
lapply(all.markers, nrow) ## 0 = 1127, 1 = 931, 2 = 1494, 3 = 1287, 4 = 1390, 5 = 1607, 6 = 1337, 7 = 803, 8 = 976, 9 = 1547, 10 = 1493,
## 11 = 897, 12 = 1275, 13 = 1228, 14 = 808, 15 = 1507, 16 = 1800, 17 = 760, 18 = 1201, 19 = 1415, 20 = 1166, 21 = 1104

p1.markers <- lapply(p1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 400)
})
p2.markers <- lapply(p2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
p3.markers <- lapply(p3.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 200)
})
all.markers <- lapply(all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 500)
})

## Give your clusters unique names
names(p1.markers) <- paste0("P1_", names(p1.markers))
names(p2.markers) <- paste0("P2_", names(p2.markers))
names(p3.markers) <- paste0("P3_", names(p3.markers))
names(all.markers) <- paste0("All_", names(all.markers))

## Make a big list of all the cluster markers
patient.markers.list <- c(p1.markers,
                          p2.markers,
                          p3.markers,
                          all.markers)

## Define your similarity function, you can use another but Jaccard works well for this
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}
## Set up the matrix you will populate with data
patient.jaccard.mat <- matrix(data = NA, nrow = length(patient.markers.list),
                              ncol = length(patient.markers.list),
                              dimnames = list(names(patient.markers.list),
                                              names(patient.markers.list)))
## Run your pairwise Jaccard similarity
for (i in rownames(patient.jaccard.mat)) {
  for (j in colnames(patient.jaccard.mat)) {
    patient.jaccard.mat[i,j] <- jaccard(patient.markers.list[[i]]$gene, patient.markers.list[[j]]$gene)
  }
}

anno.df <- data.frame("Sample" = gsub("_\\d+", "", colnames(patient.jaccard.mat)),
                      row.names = colnames(patient.jaccard.mat))
anno.col <- list("Sample" = c("All" = "black",
                              "P1" = as.character(condition.cols[1]),
                              "P2" = as.character(condition.cols[3]),
                              "P3" = as.character(condition.cols[5])))
dev.off()
gt <- pheatmap(patient.jaccard.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               annotation_row = anno.df,
               annotation_col = anno.df,
               annotation_colors = anno.col)$gtable
ggsave("plots/jaccard/jaccard_heatmap.png", plot = gt, height = 16.5, width = 16.5)
ggsave("plots/jaccard/jaccard_heatmap.pdf", plot = gt, height = 16.5, width = 16.5)

## Run the parameter search for PAM
patient.jaccard.dist <- as.dist(1 - patient.jaccard.mat)
silhouette.res <- numeric()
## Run PAM for 2-15 clusters and see what gives you the highest silhouette
for (k in 2:30) {  # Assuming you want to check from 2 to 30 clusters
  pam.fit <- pam(patient.jaccard.dist, k, diss = TRUE)
  silhouette.res[k] <- mean(silhouette(pam.fit)[,"sil_width"])
}
silhouette.res <- silhouette.res[-1]
names(silhouette.res) <- 2:30
names(which.max(silhouette.res)) ## gives you the k that is best for your data
## For this data it was 11
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:31))
## Draw the geom_vline at 10 for this data
ggplot(silhouette.df, mapping = aes(x = k, y = silhouette, group = 1)) +
  geom_line() +
  geom_vline(xintercept = 10, colour = "red", linetype = "solid", linewidth = 0.8) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor.y = element_line(linetype = "dotted"),
        panel.grid.major.y = element_line(linetype = "dotted")) +
  labs(x = "k", y = "Silhouette", title = "PAM clustering silhouette score")
ggsave("plots/jaccard/pam_clustering_silhouette_score.png", width = 8.3, height = 5.8)
ggsave("plots/jaccard/pam_clustering_silhouette_score.pdf", width = 8.3, height = 5.8)

## Run the PAM and get the cluster assignment
patient.k11.pam <- pam(patient.jaccard.dist, k = 11, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- patient.k11.pam
## Assign the clusters back to your Seurat objects
patient.clusters.df <- anno.df$PAM.Cluster
names(patient.clusters.df) <- rownames(anno.df)
p1.set <- paste0("P1_", as.character(patient.list$P1$seurat_clusters.0.2))
p2.set <- paste0("P2_", as.character(patient.list$P2$seurat_clusters.0.4))
all.set <- paste0("All_", as.character(patient.seurat$seurat_clusters.0.4))

patient.list$P1$PAM.Cluster <- factor(as.character(patient.clusters.df[p1.set]))
patient.list$P2$PAM.Cluster <- factor(as.character(patient.clusters.df[p2.set]))
patient.seurat$PAM.Cluster <- factor(as.character(patient.clusters.df[all.set]))
saveRDS(patient.list, "data/patient_seurat_split_list_pam.rds")
saveRDS(patient.seurat, "data/patient_seurat_hvgs_pam.rds")

## [ Cell annotation ] ----

patient.seurat <- readRDS("data/patient_seurat_hvgs_meta.rds")
patient.sce <- as.SingleCellExperiment(patient.seurat)
hpca.se <- celldex::HumanPrimaryCellAtlasData()

library(SingleR)
pred.patient <- SingleR(test = patient.sce, ref = hpca.se, assay.type.test = 1, labels = hpca.se$label.main)
saveRDS(pred.patient, "data/patient_singler_prediction.rds")

table(pred.patient$labels)
## Astrocyte               B_cell         Chondrocytes                   DC Embryonic_stem_cells    Endothelial_cells     Epithelial_cells 
##      6669                   26                 4891                  133                    8                 1781                   15 
## Erythroblast          Fibroblasts                  GMP          Hepatocytes           HSC_-G-CSF            HSC_CD34+            iPS_cells 
##            3                 1541                    3                44920                    1                    1                  294 
## Keratinocytes           Macrophage             Monocyte                  MSC Neuroepithelial_cell              Neurons              NK_cell 
##             1                 1295                  573                  400                  322                 3634                   65 
## Osteoblasts            Platelets  Smooth_muscle_cells              T_cells    Tissue_stem_cells 
##        2376                   41                 6012                  145                  552 

pdf("plots/anno/patient_singler_heatmap.pdf", width = 8.3, height = 8.3)
png("plots/anno/patient_singler_heatmap.png", units = "in", res = 200, width = 8.3, height = 8.3)
plotScoreHeatmap(pred.patient)
dev.off()

plotDeltaDistribution(pred.patient, ncol = 3)
ggsave("plots/anno/patient_singler_delta.pdf", width = 5.8, height = 11.7)
ggsave("plots/anno/patient_singler_delta.png", width = 5.8, height = 11.7)

colnames(pred.patient) <- paste("SingleR", colnames(pred.patient), sep = ".")
colData(patient.sce) <- merge(colData(patient.sce), pred.patient, by = "row.names", all.x = TRUE)
identical(colnames(patient.seurat), rownames(pred.patient)) ## TRUE
patient.seurat <- AddMetaData(object = patient.seurat,
                              metadata = as.data.frame(pred.patient))
saveRDS(patient.seurat, "data/patient_seurat_singler.rds")

anno.cols <- c("#D62728", "#FF9896", "#1F77B4", "#FFBB78",
               "#FF7F0E", "#17BECF", "#9467BD", "#2CA02C",
               "#E377C2", "#AEC7E8", "#843C39", "#BCBD22",
               "#9EDAE5", "#393B79", "#C5B0D5", "#8C564B",
               "#7F7F7F", "#F7B6D2", "#7B4173", "#C7C7C7",
               "#DBDB8D", "#C49C94", "#98DF8A", "#637939",
               "#8C6D31", "#525252")
names(anno.cols) <- unique(patient.seurat$SingleR.labels)
saveRDS(anno.cols, "data/singler_cols.rds")

DimPlot(patient.seurat, group.by = "SingleR.labels", order = TRUE) +
  scale_colour_manual(values = anno.cols) +
  umap.theme() + labs(title = "SingleR annotation") +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
ggsave("plots/umaps/patient_umap_singler.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/patient_umap_singler.png", width = 8.3, height = 5.8)

## [ InferCNV ] ----

## Gene reference
patient.seurat <- readRDS("data/patient_seurat_singler.rds")
ens.bm <- read.csv("data/patient_ensembl_biomart.csv", row.names = NULL)
ens.ref <- ens.bm[ens.bm$external_gene_name %in% rownames(patient.seurat), ]
ens.ref <- ens.ref[!duplicated(ens.ref$external_gene_name), ]
gene.order <- ens.ref[c("external_gene_name", "chromosome_name", "start_position", "end_position")]
saveRDS(gene.order, "data/gene_order_ref.rds")
write.table(gene.order, "data/gene_order_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

anno.df <- patient.seurat@meta.data[c("CellID", "SingleR.labels")]
immune.cells <- c("B_cell", "NK_cell", "T_cells", "Monocyte", "Macrophage", "DC")
anno.df$SingleR.labels[anno.df$SingleR.labels %in% immune.cells] <- "Immune_ref"
saveRDS(anno.df, "data/cell_anno_ref.rds")
write.table(anno.df, "data/cell_anno_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

# ## Create InferCNV object
# library(infercnv)
# infercnv.obj <- CreateInfercnvObject(raw_counts_matrix = patient.seurat@assays$RNA$counts,
#                                      annotations_file = "data/cell_anno_ref.txt",
#                                      delim = "\t",
#                                      gene_order_file = "data/gene_order_ref.txt",
#                                      ref_group_names = "Immune_ref")
# saveRDS(infercnv.obj, "data/infercnv_input.rds")
# options(scipen = 999)
# ## Run InferCNV
# infercnv.obj <- infercnv::run(infercnv.obj,
#                               cutoff = 0.1, ## Use 1 for smart-seq, 0.1 for 10x-genomics
#                               out_dir = "infercnv",  
#                               cluster_by_groups = F,
#                               denoise = T,
#                               HMM = F) ## Was failing hspike modelling
# saveRDS(infercnv.obj, "data/infercnv_output.rds")

infercnv.obj <- readRDS("data/infercnv/15_tumor_subclusters.leiden.infercnv_obj")

expr.mat <- infercnv.obj@expr.data
nrow(expr.mat) ## 8761
ncol(expr.mat) ## 75702

## Calculate CNV scores
cnv.scores <- colMeans(abs(expr.mat))
cnv.scores.df <- data.frame(cell = names(cnv.scores),
                            cnv_score = as.numeric(cnv.scores),
                            stringsAsFactors = FALSE)
write.table(cnv.scores.df,
            file = "data/infercnv/cnv_scores_per_cell.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

gene.anno <- readRDS("data/gene_order_ref.rds")
gene.anno$gene_idx <- match(gene.anno$external_gene_name, rownames(expr.mat))
gene.anno <- gene.anno[!is.na(gene.anno$gene_idx), ]
gene.anno$chr <- paste0("chr", gene.anno$chr)

## Calculate scores for specific chromosomes
calculate_score <- function(chr_name, expr_matrix, gene_info) {
  ## Handle whole chromosomes (e.g., "chr7")
  if (!grepl("[pq]$", chr_name)) {
    gene_indices <- gene_info$gene_idx[gene_info$chr == chr_name]
  } else {
    ## Handle chromosome arms (e.g., "chr17q")
    chr_base <- gsub("[pq]$", "", chr_name)
    arm <- substr(chr_name, nchar(chr_name), nchar(chr_name))
    ## Get centromere position (median of gene positions)
    chr_genes <- gene_info[gene_info$chr == chr_base, ]
    if (nrow(chr_genes) == 0) {
      return(NULL)
    }
    centromere_pos <- median(chr_genes$start)
    if (arm == "p") {
      gene_indices <- chr_genes$gene_idx[chr_genes$start < centromere_pos]
    } else {
      gene_indices <- chr_genes$gene_idx[chr_genes$start >= centromere_pos]
    }
  }
  if (length(gene_indices) < 10) {
    message("  WARNING: Skipping ", chr_name, " (only ", length(gene_indices), " genes)")
    return(NULL)
  }
  ## Extract expression for this chromosome/arm
  chr_expr <- expr_matrix[gene_indices, , drop = FALSE]
  ## Calculate score per cell: mean absolute deviation
  chr_scores <- colMeans(abs(chr_expr))
  message(
    "  ", chr_name, ": ", length(gene_indices), " genes, mean score = ",
    round(mean(chr_scores), 4)
  )
  return(chr_scores)
}
chr.arms <- c("chr1q", "chr2", "chr8", "chr20", ## Frequent gains
              "chr1p", "chr4", "chr11q", "chr18q", "chr5p", "chr16q") ## Frequent losses
## Tomlinson and Kappler, (2012); doi:10.1002/pbc.24213
## Barros et al., (2021); doi: 10.3389/fonc.2021.741526
## 1p also gain? Wu et al., (2013); doi: 10.1007/s12072-012-9350-y
chr.scores.list <- list()
for (arm in chr.arms) {
  scores <- calculate_score(arm, expr.mat, gene.anno)
  if (!is.null(scores)) {
    chr.scores.list[[arm]] <- scores
  }
}

## Add CNV scores
patient.seurat <- readRDS("data/patient_seurat_singler.rds")
patient.seurat$infercnv_score <- NA
common.cells <- intersect(colnames(patient.seurat), cnv.scores.df$cell) ## 75702
patient.seurat@meta.data[common.cells, "infercnv_score"] <- cnv.scores.df$cnv_score[match(common.cells, cnv.scores.df$cell)]

## Add chromosome-specific scores
for (arm in names(chr.scores.list)) {
  col_name <- paste0(arm, "_score")
  patient.seurat@meta.data[[col_name]] <- NA
  chr.common.cells <- intersect(colnames(patient.seurat), names(chr.scores.list[[arm]]))
  patient.seurat@meta.data[chr.common.cells, col_name] <- chr.scores.list[[arm]][chr.common.cells]
}

saveRDS(patient.seurat, "data/patient_seurat_infercnv.rds")

## CNV score histograms
ggplot(patient.seurat@meta.data, aes(x = infercnv_score)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.001,
                 fill = "#1F77B4", colour = "#1F77B4", alpha = 0.7) +
  geom_density(colour = "black", size = 1.2) +
  labs(title = "InferCNV Score Distribution", x = "InferCNV score") +
  theme_bw() +
  theme(text = element_text(size = 15))
ggsave("plots/infercnv/infercnv_score_histogram.pdf", width = 8.3, height = 5.8)
ggsave("plots/infercnv/infercnv_score_histogram.png", width = 8.3, height = 5.8)

chr.score <- colnames(patient.seurat@meta.data)
chr.score <- chr.score[grep("^chr.*_score$", chr.score)]
make_title <- function(x) {
  x |> 
    gsub("_", " ", x = _) |> 
    tools::toTitleCase()
}
for (score in chr.score) {
  name <- score
  p <- ggplot(patient.seurat@meta.data, aes(x = .data[[score]])) +
    geom_histogram(aes(y = ..density..), binwidth = 0.01,
                   fill = "#1F77B4", colour = "#1F77B4", alpha = 0.7) +
    geom_density(colour = "black", size = 1.2) +
    labs(title = paste(make_title(score), "Distribution"), x = "InferCNV score") +
    theme_bw() +
    theme(text = element_text(size = 15))
  ggsave(plot = p, paste0("plots/infercnv/", score, "_histogram.pdf"), width = 8.3, height = 5.8)
  ggsave(plot = p, paste0("plots/infercnv/", score, "_histogram.png"), width = 8.3, height = 5.8)
}

## CNV score UMAPs
FeaturePlot(patient.seurat, features = "infercnv_score") +
  scale_colour_gradientn(colours = pals::coolwarm(100)) +
  labs(title = "InferCNV Score", colour = "") +
  umap.theme() +
  theme(text = element_text(size = 15),
        legend.key.size = unit(1, "cm"))
ggsave("plots/infercnv/infercnv_score_umap.pdf", width = 8.3, height = 5.8)
ggsave("plots/infercnv/infercnv_score_umap.png", width = 8.3, height = 5.8)

for (score in chr.score) {
  FeaturePlot(patient.seurat, features = score) +
    scale_colour_gradientn(colours = pals::coolwarm(100)) +
    labs(title = paste(make_title(score)), colour = "") +
    umap.theme() +
    theme(text = element_text(size = 15),
          legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/infercnv/", score, "_umap.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/infercnv/", score, "_umap.png"), width = 8.3, height = 5.8)
}
