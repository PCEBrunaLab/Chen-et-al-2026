## HB patient 2 snRNA-seq analysis

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
clust.cols <- c("1" = "#c85f5cff",
                "2" = "#d67673ff",
                "3" = "#de8f8cff",
                "4" = "#e9aaa8ff",
                "5" = "#f1bcbcff",
                "6" = "#f6d0ceff",
                "7" = "#fbe3e2ff",
                "8" = "#bfe0ecff",
                "9" = "#c0e3baff")

## [ Patient 2 prep ] ----

patient2.seurat <- RunPCA(patient2.seurat, verbose = FALSE, npcs = 100, 
                          features = VariableFeatures(patient2.seurat))
ElbowPlot(object = patient2.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate two separate UMAP dim reds +/-
harmony.seurat <- harmony::RunHarmony(patient2.seurat, group.by.vars = "Condition_Rep",
                                      theta = 1, lambda = 1, sigma = 0.07,
                                      assay.use = "RNA", reduction = "pca",
                                      dims.use = 1:30, reduction.save = "Harmony",
                                      max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

umap.gg <- DimPlot(harmony.seurat, group.by = "Condition_Rep", order = TRUE) +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient2.seurat$Condition_Rep)))) +
  umap.theme() + labs(title = "Conditions")
umap.gg
ggsave("patient2/plots/filter_hvgs/umap_harmony_t1_l1_s007.pdf", width = 8.3, height = 5.8)
ggsave("patient2/plots/filter_hvgs/umap_harmony_t1_l1_s007.png", width = 8.3, height = 5.8)

patient2.seurat <- RunUMAP(patient2.seurat, reduction = "pca", dims = 1:30)
patient2.seurat <- FindNeighbors(patient2.seurat, reduction = "pca", dims = 1:30)
DimPlot(patient2.seurat, group.by = "Condition_Rep", order = TRUE) +
  scale_colour_manual(values = hues::iwanthue(length(unique(patient2.seurat$Condition_Rep)))) +
  umap.theme() + labs(title = "Conditions")

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
saveRDS(harmony.seurat, "patient2/data/patient2_seurat_hvgs_harmony.rds")

## [ Jaccard similarity ] ----

## Split data
patient2.seurat <- readRDS("patient2/data/patient2_seurat_hvgs_harmony.rds")

## Look through UMAPs and choose best cluster resolution
DimPlot(patient2.seurat, group.by = "seurat_clusters.0.2", order = TRUE) +
  umap.theme()
Idents(patient2.seurat) <- patient2.seurat$seurat_clusters.0.2

## Get cluster markers - MAST takes longer to run
library(future)
parallel::detectCores() ## 10
plan(multisession, workers = 8) ## 6-8 recommended on M1 Pro
options(future.globals.maxSize = 8 * 1024^3) ## 16 GB RAM
all.markers <- FindAllMarkers(patient2.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
saveRDS(all.markers, "patient2/data/patient2_all_markers.rds")

p2f1.markers <- readRDS("data/jaccard/P2_F1_markers.rds")
p2f2.markers <- readRDS("data/jaccard/P2_F2_markers.rds")

## Filter significant and higher expressed markers
p2f1.markers <- p2f1.markers[p2f1.markers$p_val_adj < 0.05, ]
p2f2.markers <- p2f2.markers[p2f2.markers$p_val_adj < 0.05, ]
all.markers <- all.markers[all.markers$p_val_adj < 0.05, ]

p2f1.markers <- p2f1.markers[p2f1.markers$avg_log2FC > 0.5, ]
p2f2.markers <- p2f2.markers[p2f2.markers$avg_log2FC > 0.5, ]
all.markers <- all.markers[all.markers$avg_log2FC > 0.5, ]

## Split into lists of genes per cluster
p2f1.markers <- split(p2f1.markers, p2f1.markers$cluster)
p2f2.markers <- split(p2f2.markers, p2f2.markers$cluster)
all.markers <- split(all.markers, all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
## Yura uses 500 per cluster
lapply(p2f1.markers, nrow) ## 0 = 214, 1 = 538, 2 = 631, 3 = 502, 4 = 304, 5 = 523, 6 = 292, 7 = 1336, 8 = 632, 9 = 1244
lapply(p2f2.markers, nrow) ## 0 = 226, 1 = 118, 2 = 339, 3 = 394, 4 = 598, 5 = 438, 6 = 643, 7 = 572, 8 = 690
lapply(all.markers, nrow) ## 0 = 588, 1 = 458, 2 = 510, 3 = 609, 4 = 625, 5 = 954, 6 = 565, 7 = 319, 8 = 669, 9 = 619

p2f1.markers <- lapply(p2f1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
p2f2.markers <- lapply(p2f2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
all.markers <- lapply(all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})

## Give your clusters unique names
names(p2f1.markers) <- paste0("P2_F1_", names(p2f1.markers))
names(p2f2.markers) <- paste0("P2_F2_", names(p2f2.markers))
names(all.markers) <- paste0("P2_all_", names(all.markers))

## Make a big list of all the cluster markers
patient.markers.list <- c(p2f1.markers, p2f2.markers, all.markers)

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
anno.col <- list("Sample" = c("P2_all" = "black",
                              "P2_F1" = as.character(condition.cols[1]),
                              "P2_F2" = as.character(condition.cols[2])))
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
ggsave("patient2/plots/jaccard/jaccard_heatmap.png", plot = gt, height = 16.5, width = 16.5)
ggsave("patient2/plots/jaccard/jaccard_heatmap.pdf", plot = gt, height = 16.5, width = 16.5)

## Run the parameter search for PAM
patient.jaccard.dist <- as.dist(1 - patient.jaccard.mat)
silhouette.res <- numeric()
## Ran PAM for 2-15 clusters
## 15 highest silhouette so run for 2-20 clusters
library(cluster)
for (k in 2:20) {
  pam.fit <- pam(patient.jaccard.dist, k, diss = TRUE)
  silhouette.res[k] <- mean(silhouette(pam.fit)[,"sil_width"])
}
silhouette.res <- silhouette.res[-1]
names(silhouette.res) <- 2:20
names(which.max(silhouette.res)) ## gives you the k that is best for your data
## For this data it was 15
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:21))
## Draw the geom_vline at 14 for this data
ggplot(silhouette.df, mapping = aes(x = k, y = silhouette, group = 1)) +
  geom_line() +
  geom_vline(xintercept = 14, colour = "red", linetype = "solid", linewidth = 0.8) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor.y = element_line(linetype = "dotted"),
        panel.grid.major.y = element_line(linetype = "dotted")) +
  labs(x = "k", y = "Silhouette", title = "PAM clustering silhouette score")
ggsave("patient2/plots/jaccard/pam_clustering_silhouette_score.png", width = 8.3, height = 5.8)
ggsave("patient2/plots/jaccard/pam_clustering_silhouette_score.pdf", width = 8.3, height = 5.8)

## Run the PAM and get the cluster assignment
patient.k15.pam <- pam(patient.jaccard.dist, k = 15, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- patient.k15.pam
## Assign the clusters back to your Seurat objects
patient.clusters.df <- anno.df$PAM.Cluster
names(patient.clusters.df) <- rownames(anno.df)
all.set <- paste0("P2_all_", as.character(patient2.seurat$seurat_clusters.0.2))
patient2.seurat$PAM.Cluster <- factor(as.character(patient.clusters.df[all.set]))
saveRDS(patient2.seurat, "patient2/data/patient2_seurat_hvgs_pam.rds")

DimPlot(patient2.seurat, group.by = "Condition", order = TRUE) +
  scale_colour_manual(values = condition.cols) +
  umap.theme() + labs(title = "Conditions")
ggsave("patient2/plots/umaps/umap_harmony_condition.pdf", width = 8.3, height = 5.8)
ggsave("patient2/plots/umaps/umap_harmony_condition.png", width = 8.3, height = 5.8)

patient2.seurat$PAM.Cluster <- factor(patient2.seurat$PAM.Cluster, levels = 1:15)
pam.cols <- hues::iwanthue(length(unique(patient2.seurat$PAM.Cluster)))
names(pam.cols) <- unique(patient2.seurat$PAM.Cluster)
saveRDS(pam.cols, "patient2/data/pam_cluster_cols.rds")
DimPlot(patient2.seurat, group.by = "PAM.Cluster", order = TRUE) +
  scale_colour_manual(values = pam.cols) +
  umap.theme() + labs(title = "PAM Cluster")
ggsave("patient2/plots/umaps/umap_harmony_pam.pdf", width = 8.3, height = 5.8)
ggsave("patient2/plots/umaps/umap_harmony_pam.png", width = 8.3, height = 5.8)

## [ PAM annotation ] ----

patient2.seurat <- readRDS("patient2/data/patient2_seurat_hvgs_pam.rds")

pseudo.pam <- AggregateExpression(patient2.seurat, assays = "RNA", return.seurat = T,
                                  group.by = c("PAM.Cluster"))
saveRDS(pseudo.pam, "patient2/data/patient2_pseudobulk_pam.rds")

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
gt <- pheatmap(pam.hb.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100))$gtable
ggsave("patient2/plots/class/pseudo_pam_hb_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/class/pseudo_pam_hb_pheatmap.png", plot = gt)

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
gt <- pheatmap(pam.wes.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100))$gtable
ggsave("patient2/plots/class/pseudo_pam_wesley_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/class/pseudo_pam_wesley_pheatmap.png", plot = gt)

custom <- c("PROM1", "NCAM1", "EPCAM", "CLDN3", ## stemmy
            "ICAM1", "KRT8", "KRT18", "KRT19", ## progenitor
            "ALB", "AFP", "CEBPA", "HNF4A", "MAT1A") ## hepatocyte
custom %in% rownames(pseudo.pam)

pam.custom.mat <- as.matrix(pseudo.pam@assays$RNA$scale.data[custom, ])
gt <- pheatmap(pam.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100))$gtable
ggsave("patient2/plots/class/pseudo_pam_custom_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/class/pseudo_pam_custom_pheatmap.png", plot = gt)

custom.df <- as.data.frame(read_xlsx(path = "/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/custom_sig.xlsx", col_names = FALSE))
colnames(custom.df) <- c("Gene", "Signature")
custom.sig <- lapply(unique(custom.df$Signature), function(x){
  custom.df[custom.df$Signature==x, "Gene"]
})
names(custom.sig) <- unique(custom.df$Signature)
custom.sig
custom.names <- paste(names(custom.sig), ".Sig", sep = "")

pseudo.pam <- AddModuleScore(pseudo.pam, features = custom.sig, assay = "RNA", seed = 12345, 
                             name = custom.names)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["Custom_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.pam, vars = custom.names)))
pam.sig.mat <- as.matrix(pseudo.pam@assays$Custom_sigs@data)
gt <- pheatmap(pam.sig.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100))$gtable
ggsave("patient2/plots/class/pseudo_pam_custom_sig_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/class/pseudo_pam_custom_sig_pheatmap.png", plot = gt)

saveRDS(pseudo.pam, "patient2/data/patient2_pseudobulk_pam_sigs.rds")

## HuH6-derived marker signature
pseudo.pam <- readRDS("patient2/data/patient_pseudobulk_pam_sigs.rds")

unique.lists <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/huh6_temp_sig.rds")
pseudo.pam <- AddModuleScore(pseudo.pam, features = unique.lists, name = paste0(names(unique.lists), ".Sig"),
                             assay = "RNA", seed = 12345)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["Markers_sigs"]] <- CreateAssayObject(
  data = t(FetchData(object = pseudo.pam,
                     vars = c("Hep.Sig", "Prog.Sig", "Stem.Sig"))))
pam.marker.mat <- as.matrix(pseudo.pam@assays$Markers_sigs@data)
gt <- pheatmap(pam.marker.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100))$gtable
ggsave("patient2/plots/class/pseudo_pam_marker_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/class/pseudo_pam_marker_pheatmap.png", plot = gt)

## Add SingleR and InferCNV meta data
ref.seurat <- readRDS("data/patient_seurat_infercnv.rds")
all(colnames(patient2.seurat) %in% colnames(ref.seurat)) ## TRUE

cnv.scores <- grep("^(infercnv|chr).*_score$", colnames(ref.seurat@meta.data), value = TRUE)
meta.add <- ref.seurat@meta.data[colnames(patient2.seurat), c("SingleR.labels", cnv.scores)]
patient2.seurat <- AddMetaData(patient2.seurat, metadata = meta.add)
saveRDS(patient2.seurat, "patient2/data/patient2_seurat_hvgs_meta.rds")

## SingleR plots
pam.df <- as.data.frame.matrix(table(patient2.seurat$PAM.Cluster, patient2.seurat$SingleR.labels))
pam.cols <- readRDS("patient2/data/pam_cluster_cols.rds")
singler.cols <- readRDS("data/singler_cols.rds")
pam.anno.cols <- list(PAM = pam.cols,
                      CellType = singler.cols)
pam.anno.row <- data.frame(PAM = factor(rownames(pam.df)))
rownames(pam.anno.row) <- rownames(pam.df)
pam.anno.col <- data.frame(CellType = factor(colnames(pam.df)))
rownames(pam.anno.col) <- colnames(pam.df)
pam.log.df <- log1p(pam.df)
gt <- pheatmap(pam.log.df,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               main = "log(Cell Number + 1)",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               annotation_row = pam.anno.row,
               annotation_col = pam.anno.col,
               annotation_colors = pam.anno.cols)
ggsave(plot = gt, "patient2/plots/anno/heatmap_pam_cell.pdf", width = 11.7, height = 16.5)
ggsave(plot = gt, "patient2/plots/anno/heatmap_pam_cell.png", width = 11.7, height = 16.5)

meta.df <- patient2.seurat@meta.data[, c("SingleR.labels", "PAM.Cluster")]
cluster.df <- meta.df %>%
  group_by(PAM.Cluster, SingleR.labels) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PAM.Cluster) %>%
  mutate(prop = n / sum(n))
ggplot(cluster.df, aes(x = PAM.Cluster, y = prop, fill = SingleR.labels)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = singler.cols) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Proportion") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        aspect.ratio = 0.8)
ggsave("patient2/plots/anno/singler_bar_pam.pdf", width = 8.3, height = 5.8)
ggsave("patient2/plots/anno/singler_bar_pam.png", width = 8.3, height = 5.8)

## InferCNV plots
make_title <- function(x) {
  x |> 
    gsub("_", " ", x = _) |> 
    tools::toTitleCase()
}
for (score in cnv.scores) {
  FeaturePlot(patient2.seurat, features = score) +
    scale_colour_gradientn(colours = pals::coolwarm(100),
                           values = scales::rescale(c(min(patient2.seurat[[score]]),
                                                      1,
                                                      max(patient2.seurat[[score]])))) +
    labs(title = paste(make_title(score)), colour = "") +
    umap.theme() +
    theme(text = element_text(size = 15),
          legend.key.size = unit(1, "cm"))
  ggsave(paste0("patient2/plots/anno/infercnv/", score, "_umap.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("patient2/plots/anno/infercnv/", score, "_umap.png"), width = 8.3, height = 5.8)
}

for (score in cnv.scores) {
  Idents(patient2.seurat) <- patient2.seurat$PAM.Cluster
  VlnPlot(patient2.seurat, features = score,
          pt.size = 0, cols = pam.cols) &
    geom_hline(yintercept = 1, linetype = "dashed") &
    xlab("") &
    labs(title = paste(make_title(score))) &
    theme(legend.position = "none",
          text = element_text(size = 12),
          aspect.ratio = 0.5)
  ggsave(paste0("patient2/plots/anno/infercnv/", score, "_pam_violin.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("patient2/plots/anno/infercnv/", score, "_pam_violin.png"), width = 8.3, height = 5.8)
}

infercnv.meta <- patient2.seurat@meta.data[, c("PAM.Cluster", cnv.scores), drop = FALSE]
patient2.agg <- aggregate(infercnv.meta[, cnv.scores],
                          by = list(cluster = infercnv.meta[["PAM.Cluster"]]),
                          FUN = mean,
                          na.rm = TRUE)
rownames(patient2.agg) <- patient2.agg$cluster
patient2.agg$cluster <- NULL
patient2.mat <- as.matrix(patient2.agg)
gt <- pheatmap(patient2.mat,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE)
ggsave(plot = gt, "patient2/plots/anno/heatmap_pam_infercnv.pdf", width = 5.8, height = 5.8)
ggsave(plot = gt, "patient2/plots/anno/heatmap_pam_infercnv.png", width = 5.8, height = 5.8)

## [ Cluster annotation ] ----

patient2.seurat <- readRDS("patient2/data/patient2_seurat_hvgs_pam.rds")

DimPlot(patient2.seurat, group.by = "Condition", order = TRUE) +
  umap.theme() + labs(title = "Condition") +
  scale_colour_manual(values = condition.cols) +
  theme(text = element_text(size = 15))
ggsave("patient2/plots/figure3/umap_harmony_conditions.pdf", width = 8.3, height = 5.8)
ggsave("patient2/plots/figure3/umap_harmony_conditions.png", width = 8.3, height = 5.8)

## Annotate PAM clusters
patient2.seurat$Exprs_clust <- recode(as.character(patient2.seurat$PAM.Cluster),
                                      "1" = "5",
                                      "2" = "3",
                                      "3" = "4",
                                      "4" = "8",
                                      "6" = "9",
                                      "7" = "7",
                                      "8" = "6",
                                      "14" = "2",
                                      "15" = "1")

DimPlot(patient2.seurat, group.by = "Exprs_clust", order = TRUE) +
  umap.theme() + labs(title = "Expression clusters") +
  scale_colour_manual(values = clust.cols) +
  theme(text = element_text(size = 15))
ggsave("patient2/plots/figure3/umap_harmony_pam_nums.pdf", width = 8.3, height = 5.8)
ggsave("patient2/plots/figure3/umap_harmony_pam_nums.png", width = 8.3, height = 5.8)

saveRDS(patient2.seurat, "patient2/data/patient2_seurat_hvgs_clust.rds")

anno.data <- patient2.seurat@meta.data[c("Condition", "Exprs_clust")]
anno.data <- anno.data[!(anno.data$Exprs_clust == 9), ]
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
ggsave("patient2/plots/figure3/exprs_clust_proportion.pdf", width = 5.8, height = 5.8)
ggsave("patient2/plots/figure3/exprs_clust_proportion.png", width = 5.8, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## Annotation heatmaps /////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

filt.seurat <- subset(patient2.seurat, subset = Exprs_clust %in% 1:8)

pseudo.pam <- AggregateExpression(filt.seurat, assays = "RNA", return.seurat = T,
                                  group.by = c("Exprs_clust"))
saveRDS(pseudo.pam, "patient2/data/patient2_pseudobulk_clust.rds")

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
anno.df <- data.frame(Exprs_clust = 1:8)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.hb.mat)
anno.cols <- clust.cols[1:8]
gt <- pheatmap(pam.hb.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:8,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient2/plots/figure3/pseudo_pam_hb_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/figure3/pseudo_pam_hb_pheatmap.png", plot = gt)

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
anno.df <- data.frame(Exprs_clust = 1:8)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.wes.mat)
# anno.cols <- clust.cols[1:8]
gt <- pheatmap(pam.wes.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:8,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient2/plots/figure3/pseudo_pam_wesley_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/figure3/pseudo_pam_wesley_pheatmap.png", plot = gt)

custom <- c("PROM1", "NCAM1", "EPCAM", "CLDN3", ## stemmy
            "ICAM1", "KRT8", "KRT18", "KRT19", ## progenitor
            "ALB", "AFP", "CEBPA", "HNF4A", "MAT1A") ## hepatocyte
custom %in% rownames(pseudo.pam)

pam.custom.mat <- as.matrix(pseudo.pam@assays$RNA$scale.data[custom, ])
anno.df <- data.frame(Exprs_clust = 1:8)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.custom.mat)
# anno.cols <- clust.cols[1:8]
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
               labels_col = 1:8,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient2/plots/figure3/pseudo_pam_custom_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/figure3/pseudo_pam_custom_pheatmap.png", plot = gt)

saveRDS(pseudo.pam, "patient2/data/patient2_pseudobulk_clust_sigs.rds")

## HuH6-derived marker signature
unique.lists <- readRDS("/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/hb sc analysis/huh6/data/huh6_temp_sig.rds")
pseudo.pam <- AddModuleScore(pseudo.pam, features = unique.lists, name = paste0(names(unique.lists), ".Sig"),
                             assay = "RNA", seed = 12345)
colnames(pseudo.pam@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.pam@meta.data))
pseudo.pam[["Markers_sigs"]] <- CreateAssayObject(
  data = t(FetchData(object = pseudo.pam,
                     vars = c("Hep.Sig", "Prog.Sig", "Stem.Sig"))))
pam.marker.mat <- as.matrix(pseudo.pam@assays$Markers_sigs@data)
anno.df <- data.frame(Exprs_clust = 1:8)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(pam.marker.mat)
# anno.cols <- clust.cols[1:8]
gt <- pheatmap(pam.marker.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:8,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggsave("patient2/plots/figure3/pseudo_pam_marker_pheatmap.pdf", plot = gt)
ggsave("patient2/plots/figure3/pseudo_pam_marker_pheatmap.png", plot = gt)

## Add SingleR and InferCNV meta data
patient2.seurat <- readRDS("patient2/data/patient2_seurat_hvgs_clust.rds")
ref.seurat <- readRDS("data/patient_seurat_infercnv.rds")
all(colnames(patient2.seurat) %in% colnames(ref.seurat)) ## TRUE

cnv.scores <- grep("^(infercnv|chr).*_score$", colnames(ref.seurat@meta.data), value = TRUE)
meta.add <- ref.seurat@meta.data[colnames(patient2.seurat), c("SingleR.labels", cnv.scores)]
patient2.seurat <- AddMetaData(patient2.seurat, metadata = meta.add)
saveRDS(patient2.seurat, "patient2/data/patient2_seurat_hvgs_clust_meta.rds")

## SingleR plots
clust.df <- as.data.frame.matrix(table(patient2.seurat$Exprs_clust, patient2.seurat$SingleR.labels))
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
ggsave(plot = gt, "patient2/plots/anno/heatmap_clust_cell.pdf", width = 11.7, height = 16.5)
ggsave(plot = gt, "patient2/plots/anno/heatmap_clust_cell.png", width = 11.7, height = 16.5)

meta.df <- patient2.seurat@meta.data[, c("SingleR.labels", "Exprs_clust")]
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
ggsave("patient2/plots/anno/singler_bar_clust.pdf", width = 8.3, height = 5.8)
ggsave("patient2/plots/anno/singler_bar_clust.png", width = 8.3, height = 5.8)

## InferCNV plots
infercnv.meta <- patient2.seurat@meta.data[, c("Exprs_clust", cnv.scores), drop = FALSE]
patient2.agg <- aggregate(infercnv.meta[, cnv.scores],
                          by = list(cluster = infercnv.meta[["Exprs_clust"]]),
                          FUN = mean,
                          na.rm = TRUE)
rownames(patient2.agg) <- patient2.agg$cluster
patient2.agg$cluster <- NULL
patient2.mat <- as.matrix(patient2.agg)
gt <- pheatmap(patient2.mat,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 12, fontsize_col = 12,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE)
ggsave(plot = gt, "patient2/plots/anno/heatmap_clust_infercnv.pdf", width = 5.8, height = 5.8)
ggsave(plot = gt, "patient2/plots/anno/heatmap_clust_infercnv.png", width = 5.8, height = 5.8)
