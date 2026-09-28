## Patient 1 snRNA-seq analysis

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

## Colours for PAM clusters
clust.cols <- c("1" = "#d67673ff",
                "2" = "#de8f8cff",
                "3" = "#e9aaa8ff",
                "4" = "#f3c7c5ff",
                "5" = "#9fcf9aff",
                "6" = "#c0e3baff",
                "7" = "#e8f3dcff",
                "8" = "#8fc9b6ff",
                "9" = "#d9d1a6ff")

## [ Patient 1 prep ] ----

patient1.seurat <- RunPCA(patient1.seurat, verbose = FALSE, npcs = 100, 
                          features = VariableFeatures(patient1.seurat))
ElbowPlot(object = patient1.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate two separate UMAP dim reds +/-
harmony.seurat <- harmony::RunHarmony(patient1.seurat, group.by.vars = "Condition_Rep",
                                      theta = 0.5, lambda = 1, sigma = 0.01,
                                      assay.use = "RNA", reduction = "pca",
                                      dims.use = 1:30, reduction.save = "Harmony",
                                      max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

umap.gg <- DimPlot(harmony.seurat, group.by = "Condition_Rep", order = TRUE) +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient1.seurat$Condition_Rep)))) +
  umap.theme() + labs(title = "Conditions")
umap.gg
ggsave("patient1/plots/filter_hvgs/umap_harmony_t05_l1_s001.pdf", width = 8.3, height = 5.8)
ggsave("patient1/plots/filter_hvgs/umap_harmony_t05_l1_s001.png", width = 8.3, height = 5.8)

# patient1.seurat <- RunUMAP(patient1.seurat, reduction = "pca", dims = 1:30)
# patient1.seurat <- FindNeighbors(patient1.seurat, reduction = "pca", dims = 1:30)
# DimPlot(patient1.seurat, group.by = "Condition_Rep", order = TRUE) +
#   scale_colour_manual(values = hues::iwanthue(length(unique(patient1.seurat$Condition_Rep)))) +
#   umap.theme() + labs(title = "Conditions")

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
saveRDS(harmony.seurat, "patient1/data/patient1_seurat_hvgs_harmony.rds")

## [ Jaccard similarity ] ----

## Split data
patient1.seurat <- readRDS("patient1/data/patient1_seurat_hvgs_harmony.rds")

## Look through UMAPs and choose best cluster resolution
DimPlot(patient1.seurat, group.by = "seurat_clusters.0.2", order = TRUE) +
  umap.theme()
Idents(patient1.seurat) <- patient1.seurat$seurat_clusters.0.2

## Get cluster markers - MAST takes longer to run
library(future)
parallel::detectCores() ## 10
plan(multisession, workers = 8) ## 6-8 recommended on M1 Pro
options(future.globals.maxSize = 8 * 1024^3) ## 16 GB RAM
all.markers <- FindAllMarkers(patient1.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
saveRDS(all.markers, "patient1/data/patient1_all_markers.rds")

p1f1.markers <- readRDS("data/jaccard/P1_F1_markers.rds")
p1f2.markers <- readRDS("data/jaccard/P1_F2_markers.rds")

## Filter significant and higher expressed markers
p1f1.markers <- p1f1.markers[p1f1.markers$p_val_adj < 0.05, ]
p1f2.markers <- p1f2.markers[p1f2.markers$p_val_adj < 0.05, ]
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]

p1f1.markers <- p1f1.markers[p1f1.markers$avg_log2FC > 0.5, ]
p1f2.markers <- p1f2.markers[p1f2.markers$avg_log2FC > 0.5, ]
all.markers <- all.markers[all.markers$avg_log2FC > 0.5, ]

## Split into lists of genes per cluster
p1f1.markers <- split(p1f1.markers, p1f1.markers$cluster)
p1f2.markers <- split(p1f2.markers, p1f2.markers$cluster)
all.markers <- split(all.markers, all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
## Yura uses 500 per cluster
lapply(p1f1.markers, nrow) ## 0 = 1009, 1 = 601, 2 = 1559, 3 = 443, 4 = 1117, 5 = 1158, 6 = 635, 7 = 1009, 8 = 1416
lapply(p1f2.markers, nrow) ## No markers
lapply(all.markers, nrow) ## 0 = 1007, 1 = 603, 2 = 1556, 3 = 324, 4 = 843, 5 = 1102, 6 = 1009, 7 = 336, 8 = 1391

p1f1.markers <- lapply(p1f1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 400)
})
# p1f2.markers <- lapply(p1f2.markers, \(x) {
#   x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
#   head(x, 100)
# })
all.markers <- lapply(all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 400)
})

## Give your clusters unique names
names(p1f1.markers) <- paste0("P1_F1_", names(p1f1.markers))
# names(p1f2.markers) <- paste0("P1_F2_", names(p1f2.markers))
names(all.markers) <- paste0("P1_all_", names(all.markers))

## Make a big list of all the cluster markers
patient.markers.list <- c(p1f1.markers, all.markers)

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
anno.col <- list("Sample" = c("P1_all" = "black",
                              "P1_F1" = as.character(condition.cols[1]),
                              "P1_F2" = as.character(condition.cols[2])))
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
ggsave("patient1/plots/jaccard/jaccard_heatmap.png", plot = gt, height = 16.5, width = 16.5)
ggsave("patient1/plots/jaccard/jaccard_heatmap.pdf", plot = gt, height = 16.5, width = 16.5)

## Run the parameter search for PAM
patient.jaccard.dist <- as.dist(1 - patient.jaccard.mat)
silhouette.res <- numeric()
## Run PAM for 2-15 clusters and see what gives you the highest silhouette
library(cluster)
for (k in 2:15) {
  pam.fit <- pam(patient.jaccard.dist, k, diss = TRUE)
  silhouette.res[k] <- mean(silhouette(pam.fit)[,"sil_width"])
}
silhouette.res <- silhouette.res[-1]
names(silhouette.res) <- 2:15
names(which.max(silhouette.res)) ## gives you the k that is best for your data
## For this data it was 9
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:16))
## Draw the geom_vline at 8 for this data
ggplot(silhouette.df, mapping = aes(x = k, y = silhouette, group = 1)) +
  geom_line() +
  geom_vline(xintercept = 8, colour = "red", linetype = "solid", linewidth = 0.8) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor.y = element_line(linetype = "dotted"),
        panel.grid.major.y = element_line(linetype = "dotted")) +
  labs(x = "k", y = "Silhouette", title = "PAM clustering silhouette score")
ggsave("patient1/plots/jaccard/pam_clustering_silhouette_score.png", width = 8.3, height = 5.8)
ggsave("patient1/plots/jaccard/pam_clustering_silhouette_score.pdf", width = 8.3, height = 5.8)

## Run the PAM and get the cluster assignment
patient.k9.pam <- pam(patient.jaccard.dist, k = 9, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- patient.k9.pam
## Assign the clusters back to your Seurat objects
patient.clusters.df <- anno.df$PAM.Cluster
names(patient.clusters.df) <- rownames(anno.df)
all.set <- paste0("P1_all_", as.character(patient1.seurat$seurat_clusters.0.2))
patient1.seurat$PAM.Cluster <- factor(as.character(patient.clusters.df[all.set]))
saveRDS(patient1.seurat, "patient1/data/patient1_seurat_hvgs_pam.rds")

DimPlot(patient1.seurat, group.by = "Condition", order = TRUE) +
  scale_colour_manual(values = condition.cols) +
  umap.theme() + labs(title = "Conditions")
ggsave("patient1/plots/umaps/umap_harmony_condition.pdf", width = 8.3, height = 5.8)
ggsave("patient1/plots/umaps/umap_harmony_condition.png", width = 8.3, height = 5.8)

patient1.seurat$PAM.Cluster <- factor(patient1.seurat$PAM.Cluster, levels = 1:9)
pam.cols <- hues::iwanthue(length(unique(patient1.seurat$PAM.Cluster)))
names(pam.cols) <- 1:9
saveRDS(pam.cols, "patient1/data/pam_cluster_cols.rds")
DimPlot(patient1.seurat, group.by = "PAM.Cluster", order = TRUE) +
  scale_colour_manual(values = pam.cols) +
  umap.theme() + labs(title = "PAM Cluster")
ggsave("patient1/plots/umaps/umap_harmony_pam.pdf", width = 8.3, height = 5.8)
ggsave("patient1/plots/umaps/umap_harmony_pam.png", width = 8.3, height = 5.8)

## [ Cluster annotation ] ----

patient1.seurat <- readRDS("patient1/data/patient1_seurat_hvgs_pam.rds")

DimPlot(patient1.seurat, group.by = "Condition", order = TRUE) +
  umap.theme() + labs(title = "Condition") +
  scale_colour_manual(values = condition.cols) +
  theme(text = element_text(size = 15))
ggsave("patient1/plots/figure3/umap_harmony_conditions.pdf", width = 8.3, height = 5.8)
ggsave("patient1/plots/figure3/umap_harmony_conditions.png", width = 8.3, height = 5.8)

## Annotate PAM clusters
patient1.seurat$Exprs_clust <- recode(as.character(patient1.seurat$PAM.Cluster),
                                      "1" = "2",
                                      "2" = "4",
                                      "3" = "9",
                                      "4" = "3",
                                      "5" = "5",
                                      "6" = "8",
                                      "7" = "6",
                                      "8" = "7",
                                      "9" = "1")

DimPlot(patient1.seurat, group.by = "Exprs_clust", order = TRUE) +
  umap.theme() + labs(title = "Expression clusters") +
  scale_colour_manual(values = clust.cols) +
  theme(text = element_text(size = 15))
ggsave("patient1/plots/figure3/umap_harmony_pam_nums.pdf", width = 8.3, height = 5.8)
ggsave("patient1/plots/figure3/umap_harmony_pam_nums.png", width = 8.3, height = 5.8)

saveRDS(patient1.seurat, "patient1/data/patient1_seurat_hvgs_clust.rds")

anno.data <- patient1.seurat@meta.data[c("Condition", "Exprs_clust")]
anno.data <- anno.data[!(anno.data$Exprs_clust %in% 5:9), ]
anno.table <- as.data.frame(table(anno.data$Condition, anno.data$Exprs_clust))
ggplot(anno.table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("Cluster proportion") +
  guides(fill = guide_legend(title = "Expression clusters")) +
  scale_fill_manual(values = clust.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1.2)
ggsave("patient1/plots/figure3/exprs_clust_proportion.pdf", width = 5.8, height = 5.8)
ggsave("patient1/plots/figure3/exprs_clust_proportion.png", width = 5.8, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## Annotation heatmaps /////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

filt.seurat <- subset(patient1.seurat, subset = Exprs_clust %in% 1:4)

pseudo.pam <- AggregateExpression(filt.seurat, assays = "RNA", return.seurat = T,
                                  group.by = c("Exprs_clust"))
saveRDS(pseudo.pam, "patient1/data/patient1_pseudobulk_clust.rds")

## HB sigs
library(readxl)
sig.df <- as.data.frame(read_xlsx(path = "/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/hb_sigs_filt.xlsx", col_names = FALSE))
colnames(sig.df) <- c("Gene", "Signature")
sigs <- lapply(unique(sig.df$Signature), function(x){
  sig.df[sig.df$Signature==x, "Gene"]
})
names(sigs) <- unique(sig.df$Signature)
sigs
sig.names <- paste(names(sigs), ".Sig", sep = "")
sig.names <- gsub("_", "-", sig.names)

pseudo.pam <- AddModuleScore(pseudo.pam, features = sigs, assay = "RNA", seed = 12345, 
                             name = sig.names)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["HB_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.pam, vars = sig.names)))
pam.hb.mat <- as.matrix(pseudo.pam@assays$HB_sigs@data)
anno.df <- data.frame(Exprs_clust = 1:4)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.hb.mat)
anno.cols <- clust.cols[1:4]
gt <- pheatmap(pam.hb.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:4,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient1/plots/figure3/pseudo_pam_hb_pheatmap.pdf", plot = gt)
ggsave("patient1/plots/figure3/pseudo_pam_hb_pheatmap.png", plot = gt)

## Development markers (Wesley et al., 2022)
wesley.df <- as.data.frame(read_xlsx(path = "/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/welsey2022_sig.xlsx", col_names = FALSE))
colnames(wesley.df) <- c("Gene", "Signature")
wesley.sig <- lapply(unique(wesley.df$Signature), function(x){
  wesley.df[wesley.df$Signature==x, "Gene"]
})
names(wesley.sig) <- unique(wesley.df$Signature)
wesley.sig
wesley.names <- paste(names(wesley.sig), ".Sig", sep = "")

pseudo.pam <- AddModuleScore(pseudo.pam, features = wesley.sig, assay = "RNA", seed = 12345, 
                             name = wesley.names)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["Wesley_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.pam, vars = wesley.names)))
pam.wes.mat <- as.matrix(pseudo.pam@assays$Wesley_sigs@data)
anno.df <- data.frame(Exprs_clust = 1:4)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.wes.mat)
# anno.cols <- clust.cols[1:4]
gt <- pheatmap(pam.wes.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:4,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient1/plots/figure3/pseudo_pam_wesley_pheatmap.pdf", plot = gt)
ggsave("patient1/plots/figure3/pseudo_pam_wesley_pheatmap.png", plot = gt)

custom <- c("PROM1", "NCAM1", "EPCAM", "CLDN3", ## stemmy
            "ICAM1", "KRT8", "KRT18", "KRT19", ## progenitor
            "ALB", "AFP", "CEBPA", "HNF4A", "MAT1A") ## hepatocyte
custom %in% rownames(pseudo.pam)

pam.custom.mat <- as.matrix(pseudo.pam@assays$RNA$scale.data[custom, ])
anno.df <- data.frame(Exprs_clust = 1:4)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.custom.mat)
# anno.cols <- clust.cols[1:4]
hc <- hclust(dist(pam.custom.mat), method = "ward.D2")
library(dendextend)
dend <- as.dendrogram(hc)
dend <- rotate(dend, order = c("AFP", "ALB", "CEBPA", "HNF4A", "MAT1A", "EPCAM", "PROM1", "CLDN3", "NCAM1", "ICAM1", "KRT19", "KRT8", "KRT18"))
## Find best order that corresponds to marker programs for ease of visualisation
## Does not force order; works within constraints of tree topology
gt <- pheatmap(pam.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               cluster_rows = as.hclust(dend),
               color = Seurat:::SpatialColors(100),
               labels_col = 1:4,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient1/plots/figure3/pseudo_pam_custom_pheatmap.pdf", plot = gt)
ggsave("patient1/plots/figure3/pseudo_pam_custom_pheatmap.png", plot = gt)

saveRDS(pseudo.pam, "patient1/data/patient1_pseudobulk_clust_sigs.rds")

## HuH6-derived marker signature
unique.lists <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/huh6_temp_sig.rds")
pseudo.pam <- AddModuleScore(pseudo.pam, features = unique.lists, name = paste0(names(unique.lists), ".Sig"),
                             assay = "RNA", seed = 12345)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["Markers_sigs"]] <- CreateAssayObject(
  data = t(FetchData(object = pseudo.pam,
                     vars = c("Hep.Sig", "Prog.Sig", "Stem.Sig"))))
pam.marker.mat <- as.matrix(pseudo.pam@assays$Markers_sigs@data)
anno.df <- data.frame(Exprs_clust = 1:4)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.marker.mat)
# anno.cols <- clust.cols[1:4]
gt <- pheatmap(pam.marker.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:4,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient1/plots/figure3/pseudo_pam_marker_pheatmap.pdf", plot = gt)
ggsave("patient1/plots/figure3/pseudo_pam_marker_pheatmap.png", plot = gt)

## Add SingleR and InferCNV meta data
patient1.seurat <- readRDS("patient1/data/patient1_seurat_hvgs_clust.rds")
ref.seurat <- readRDS("data/patient_seurat_infercnv.rds")
all(colnames(patient1.seurat) %in% colnames(ref.seurat)) ## TRUE

cnv.scores <- grep("^(infercnv|chr).*_score$", colnames(ref.seurat@meta.data), value = TRUE)
meta.add <- ref.seurat@meta.data[colnames(patient1.seurat), c("SingleR.labels", cnv.scores)]
patient1.seurat <- AddMetaData(patient1.seurat, metadata = meta.add)
saveRDS(patient1.seurat, "patient1/data/patient1_seurat_hvgs_clust_meta.rds")

## SingleR plots
clust.df <- as.data.frame.matrix(table(patient1.seurat$Exprs_clust, patient1.seurat$SingleR.labels))
singler.cols <- readRDS("data/singler_cols.rds")
clust.anno.cols <- list(Exprs_clust = clust.cols,
                        CellType = singler.cols)
clust.anno.row <- data.frame(Exprs_clust = factor(rownames(clust.df)))
rownames(clust.anno.row) <- rownames(clust.df)
clust.anno.col <- data.frame(CellType = factor(colnames(clust.df)))
rownames(clust.anno.col) <- colnames(clust.df)
clust.log.df <- log1p(clust.df)
gt <- pheatmap(clust.log.df,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               main = "log(Cell Number + 1)",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               annotation_row = clust.anno.row,
               annotation_col = clust.anno.col,
               annotation_colors = clust.anno.cols)
ggsave(plot = gt, "patient1/plots/anno/heatmap_clust_cell.pdf", width = 11.7, height = 16.5)
ggsave(plot = gt, "patient1/plots/anno/heatmap_clust_cell.png", width = 11.7, height = 16.5)

meta.df <- patient1.seurat@meta.data[, c("SingleR.labels", "Exprs_clust")]
cluster.df <- meta.df %>%
  group_by(Exprs_clust, SingleR.labels) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Exprs_clust) %>%
  mutate(prop = n / sum(n))
ggplot(cluster.df, aes(x = Exprs_clust, y = prop, fill = SingleR.labels)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = singler.cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = 0.8)
ggsave("patient1/plots/anno/singler_bar_clust.pdf", width = 8.3, height = 5.8)
ggsave("patient1/plots/anno/singler_bar_clust.png", width = 8.3, height = 5.8)

## InferCNV plots
infercnv.meta <- patient1.seurat@meta.data[, c("Exprs_clust", cnv.scores), drop = FALSE]
patient1.agg <- aggregate(infercnv.meta[, cnv.scores],
                          by = list(cluster = infercnv.meta[["Exprs_clust"]]),
                          FUN = mean,
                          na.rm = TRUE)
rownames(patient1.agg) <- patient1.agg$cluster
patient1.agg$cluster <- NULL
patient1.mat <- as.matrix(patient1.agg)
gt <- pheatmap(patient1.mat,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE)
ggsave(plot = gt, "patient1/plots/anno/heatmap_clust_infercnv.pdf", width = 5.8, height = 5.8)
ggsave(plot = gt, "patient1/plots/anno/heatmap_clust_infercnv.png", width = 5.8, height = 5.8)
