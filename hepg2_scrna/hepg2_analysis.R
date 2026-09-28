## HepG2 scRNA-seq analysis

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(scater)
library(scran)
library(Seurat)
library(TSCAN)
library(slingshot)
library(biomaRt)
library(harmony)
library(tricycle)
library(cluster)
library(clusterProfiler)
library(msigdbr)
library(RcisTarget)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(ggbeeswarm)
library(ggsankey)
library(pheatmap)
library(reshape2)
library(patchwork)
library(viridis)
library(hues)
library(writexl)
library(readxl)

## Labels for experimental conditions
conditions.list <- c("POT", "Untreated", "Cisplatin", bquote("Recovery t"[1]), bquote("Recovery t"[2]))

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

## Colours for experimental condition
group.cols <- c("POT" = "#7f7f7fff", 
                "Untreated" = "#d2d2d2ff", 
                "Cisplatin_off" = "#e5f5daff", 
                "Cisplatin_recovery1" = "#b6d3a5ff", 
                "Cisplatin_recovery2" = "#a7daedff")

## Colours for PAM clusters
clust.cols <- c("1" = "#d67673ff",
                "2" = "#de8f8cff",
                "3" = "#e9aaa8ff",
                "4" = "#f3c7c5ff",
                "5" = "#f4c596ff",
                "6" = "#f8ddbbff",
                "7" = "#bfe0ecff")

## Colours for expression states
state.cols <- c("Hep-like"  = "#cf6f6dff",
                "Pro-like"  = "#f2b785ff",  
                "Stem-like" = "#9ecfe0ff")

## [ Data preparation ] ----

hepg2.sce <- readRDS("data/HepG2_SCE-norm.RDS")

## Demux samples
demux.list <- list.files("data/data prep/", pattern = "demuxed", full.names = TRUE)
d.list <- list()
for(x in seq_along(demux.list)){
  x.d <- read.table(demux.list[x], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x.name <- gsub(demux.list[x], pattern = "(.*/)(Multiplex_[0-9]*)(_demuxed\\.tsv)", replacement = "\\2")
  x.d$Sample <- x.name
  d.list[[x.name]] <- x.d
}

## Make demux data frame
demux.df <- do.call(rbind.data.frame, d.list)
demux.df$CellID <- paste(demux.df$Sample, demux.df$ID, sep = "_")

## Remove path format from metadata values
col.df <- as.data.frame(colData(hepg2.sce))
col.df$Sample <- gsub(pattern = "(\\S+)(scRNA-seq/)", "", col.df$Sample)
col.df$Sample <- gsub(pattern = "(/)", "", col.df$Sample)
col.df$CellID <- paste(col.df$Sample, col.df$Barcode, sep="_")

unique(col.df$Sample)

## Demux sample conditions
config.df <- read.table("data/data prep/HepG2_config.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
## Recovery column
config.df$Rec <- FALSE
rec.patterns <- c("off", "recovery")
config.df$Rec[grepl(config.df$description, pattern = paste(rec.patterns, collapse = "|"))] <- TRUE
## Condition column
config.df$Condition <- gsub(config.df$description, pattern = "(\\S+)(_[A-C])", replacement="\\1")

## Get size factors for QC
sf.df <- data.frame("CellID" = colnames(hepg2.sce), "Sample" = colData(hepg2.sce)$Sample, "sizeFactor" = sizeFactors(hepg2.sce))
sf.df$Sample <- gsub(pattern = "(\\S+)(scRNA-seq/)", "", sf.df$Sample)
sf.df$Sample <- gsub(pattern = "(/)", "", sf.df$Sample)

meta.df <- Reduce(x = list(col.df, demux.df, sf.df), f = function(x, y) merge(x, y, by = c("CellID", "Sample")))
meta.df <- merge(meta.df, config.df, by.x = c("Sample", "CMO"), by.y = c("Sample", "cmo_ids"))

## CMO and cell counts metadata
cmo.sum.files <- list.files("data/data prep/", pattern = "CMOsummary", full.names = TRUE)
cell.sum.files <- list.files("data/data prep/", pattern = "Cellsummary", full.names = TRUE)
cmo.list <- list()
cell.list <- list()
for(x in seq_along(cmo.sum.files)){
  x.cmo <- read.table(cmo.sum.files[x], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x.cell <- read.table(cell.sum.files[x], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cmo.list[[x]] <- x.cmo
  cell.list[[x]] <- x.cell
}

cmo.sum.df <- do.call(rbind.data.frame, cmo.list)
cell.sum.df <- do.call(rbind.data.frame, cell.list)

cell.sum.df$CellID <- gsub(pattern = "(\\S+)(scRNA-seq/)", "", cell.sum.df$CellID)
cell.sum.df$CellID <- gsub(pattern = "(/)", "", cell.sum.df$CellID)

cell.cmo.df <- merge(demux.df, cell.sum.df, by = c("CellID", "Sample"))

intersect(colnames(cell.cmo.df), colnames(meta.df))
test.df <- merge(cell.cmo.df, meta.df, by = c("CellID", "Sample", "ID", "Class", "CMO"))
test.df <- test.df %>%
  dplyr::select(-c(sizeFactor.y)) %>%
  mutate(sizeFactor = sizeFactor.x) %>%
  dplyr::select(-c(sizeFactor.x))
sum(is.na(test.df$CellID)) ## 0

saveRDS(test.df, "data/HepG2_metadata.rds")

## Add metadata to SCE object
final.meta <- readRDS("data/HepG2_metadata.rds")
colData(hepg2.sce)$CellID <- rownames(colData(hepg2.sce))
hepg2.sce$Sample <- NULL
## Remove cells not in metadata
hepg2.sce$CellID[!(hepg2.sce$CellID %in% final.meta$CellID)]

cell.nms <- colnames(hepg2.sce)[colnames(hepg2.sce) %in% final.meta$CellID]
hepg2.sce <- hepg2.sce[ , final.meta$CellID]
dim(hepg2.sce) ## 36613 48134
colData(hepg2.sce) <- merge(colData(hepg2.sce), final.meta, by = c("CellID", "Barcode", "sizeFactor"))
colnames(hepg2.sce) <- cell.nms 

saveRDS(hepg2.sce, "data/HepG2_sce_meta.rds")

## [ Initial filtering ] ----

hepg2.sce <- readRDS("data/HepG2_sce_meta.rds")

ensembl <- useMart("ensembl", "hsapiens_gene_ensembl", host = "http://useast.ensembl.org")
#ens.bm <- getBM(attributes = c("ensembl_gene_id", 
#                               "chromosome_name",
#                               "start_position",
#                               "end_position",
#                               "external_gene_name", 
#                               "hgnc_symbol", 
#                               "description"), 
#                mart = ensembl, 
#                filters = "ensembl_gene_id", 
#                values = rownames(hepg2.sce))
#write.csv(ens.bm, "data/ensembl_biomart.csv")
ens.bm <- read.csv("data/ensembl_biomart.csv", row.names = 1)

## Keep track of any duplicates
dup.ensg <- ens.bm$ensembl_gene_id[duplicated(ens.bm$ensembl_gene_id)]
## Make a vector and give elements names at the same time
dup.ensg <- setNames(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% dup.ensg],
                     ens.bm$ensembl_gene_id[ens.bm$ensembl_gene_id %in% dup.ensg])

## Two duplicates
## ENSG00000230417 = LINC00595
## ENSG00000276085 = CCL3L3

## Remove unnecessary annotation from the description
ens.filt.bm <- ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)

## This filters out around 200 genes
nrow(ens.filt.bm) ## 36392
nrow(hepg2.sce) ## 36613

## Keep only the simple chromosome information
table(ens.filt.bm$chromosome_name)
ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]

## Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name == "", ]

## Now we've filtered out around 10,000 genes
nrow(ens.filt.bm) ## 25310
nrow(hepg2.sce) ## 36613

## Keep genes that didn't get annotated so we know what gets dropped
ens.removed.bm <- ens.bm[!ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm) ## 11082

length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
## 10924 / 11082 are novel transcripts or novel transcripts antisense to a gene

grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
## The rest is a mixed bag: ~160 genes which are either novel, mitochondrially encoded, non-coding, or just empty description
## It seems unlikely to be particularly informative in our case

## Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])
## ARMCX5-GPRASP2
## CCL3L1 - part of chemokine cluster on chr17
## DNAJC9-AS1 - heatshock protein homolog -AS1 is readthrough
## ELFN2 - postsynaptic cell adhesion molecule in group III mGluRs
## GOLGA8M - something to do with the Golgi membrane
## GPR84-AS1 - a g-protein coupled receptor complex member readthrough
## MATR3 - nuclear matrix protein          
## MKKS - centrosomal shuttling protein
## NPIPA9 - nuclear pore complex related
## PRICKLE2-AS1 - homolog of prickle readthrough
## RAET1E-AS1 - MHC class I related genes member readthrough
## SPATA13 - spermatogenesis associated
## TNFRSF10A-DT - TNF-receptor superfamily member; DT = divergent transcript (?)

## LINC00484
## LINC00595
## LINC01115
## LINC01238
## LINC01605
## LINC03021
## LINC03023
## LINC03025

## Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm) ## 25288

mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]

ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype=="lincRNA"]

linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]

## Basic QC
hepg2.sce <- addPerCellQCMetrics(hepg2.sce, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))

hepg2.sce <- hepg2.sce[ens.filt.bm$ensembl_gene_id, ]
rownames(hepg2.sce) <- ens.filt.bm$external_gene_name
rownames(ens.filt.bm) <- ens.filt.bm$external_gene_name

metadata.df <-  as.data.frame(colData(hepg2.sce))

hepg2.seurat <- CreateSeuratObject(counts = assay(hepg2.sce, "counts"),
                                   assay = "RNA",
                                   meta.data = metadata.df)

saveRDS(hepg2.sce, "data/hepg2_sce_start.rds")
saveRDS(hepg2.seurat, "data/hepg2_seurat_start.rds")

## [ QC ] ----

hepg2.sce <- readRDS("data/hepg2_sce_start.rds")

dir.create("plots/qc/", recursive = TRUE)

table(hepg2.sce$Condition)
## Cisplatin_off    Cisplatin_recovery1   Cisplatin_recovery2   POT     Untreated 
## 13644            8649                  5559                  8870    11412 
hepg2.sce$Condition <- factor(hepg2.sce$Condition,
                              levels = c("POT", "Untreated", "Cisplatin_off", "Cisplatin_recovery1", "Cisplatin_recovery2"))

det.sce.gg <- plotColData(hepg2.sce, y = "detected", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

umi.sce.gg <- plotColData(hepg2.sce, y = "total", x = "Condition", colour_by = "Condition") + 
  labs(x = element_blank(), y = "UMI / cell", title = "Detected UMI per cell") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(sides = "l", outside = TRUE) + coord_cartesian(clip = "off")

mt.sce.gg <- plotColData(hepg2.sce, y = "subsets_mt_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

rb.sce.gg <- plotColData(hepg2.sce, y = "subsets_ribo_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

det.sce.gg + umi.sce.gg + mt.sce.gg + rb.sce.gg + plot_layout(nrow = 2) & theme(aspect.ratio = 0.5)

ggsave("plots/qc/hepg2_QC_plots.pdf", width = 8.7, height = 5.8)
ggsave("plots/qc/hepg2_QC_plots.png", width = 8.7, height = 5.8)

## Examine whether there are any genes that dominate expression in cells - large % of expression occupied by single gene
## Use sparse matrix operations, as if your dataset is large, doing matrix operations the regular way will take a long time
counts.assay <- counts(hepg2.sce)
counts.assay@x <- counts.assay@x/rep.int(colSums(counts.assay), diff(counts.assay@p))
top.expr <- order(Matrix::rowSums(counts.assay), decreasing = TRUE)[20:1]
top.expr <- as.matrix(t(counts.assay[top.expr, ]))
top.expr <- as.data.frame(top.expr)
ord.idx <- colnames(top.expr)
top.expr <- melt(top.expr, variable.name = "Gene", value.name = "Prop")
top.expr$Gene <- factor(top.expr$Gene, levels = ord.idx)

ggplot(top.expr,
       mapping = aes(x = Gene, y = Prop * 100)) +
  geom_boxplot() +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_line(linetype = "dotted"),
                     text = element_text(size = 15)) +
  coord_flip() +
  labs(y = "% total expression", x = "Gene", title = "Percentage expression of a single gene per total cell expression")

ggsave("plots/qc/hepg2_QC_topgenes.pdf", width = 8.7, height = 5.8)
ggsave("plots/qc/hepg2_QC_topgenes.png", width = 8.7, height = 5.8)

## Remove the counts as it's not needed.
rm(counts.assay)
gc()

## MALAT1 is a long non-coding RNA and typically comes up enriched in skim-seq like 10x
## It is routinely removed as part of preprocessing
## However, since we're doing drug treatments, it's worth checking to see if it might indicate stress response
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7099402/

malat1.sce.gg <- plotExpression(hepg2.sce, 
                                features = "MALAT1", 
                                exprs_values = "logcounts",
                                x = "Condition", 
                                colour_by = "Condition") +
  labs(x = element_blank(), y = expression("Log"[2]*" Expression")) +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

## During preprocessing we added a QC annotation for lincRNA

lincrna.sce.gg <- plotColData(hepg2.sce, y = "subsets_linc_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "lincRNA\nexpression % of total") +
  scale_colour_manual(values = group.cols) +
  scale_x_discrete(labels = conditions.list) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

malat1.sce.gg + lincrna.sce.gg + plot_layout(nrow = 2) &
  theme(aspect.ratio = 0.5, text = element_text(size = 15))

## Cisplatin is having some kind of increase in lincRNA, and this is increased in recovery 1 and 2

ggsave("plots/qc/hepg2_QC_lincrna_malat1.pdf", width = 8.7, height = 8.7)
ggsave("plots/qc/hepg2_QC_lincrna_malat1.png", width = 8.7, height = 8.7)

lincrna.sce.gg <- lincrna.sce.gg +
  labs(title = "Proportion of lincRNA\nin total cell expression") +
  theme(axis.text.x = element_text(size = 6))

mt.sce.gg + rb.sce.gg + lincrna.sce.gg + theme(aspect.ratio = 0.5)

ggsave("plots/qc/hepg2_QC_mt_rb_lincrna.pdf", width = 11.6, height = 5.8)
ggsave("plots/qc/hepg2_QC_mt_rb_lincrna.png", width = 11.6, height = 5.8)

gc()

## Cell cycle prediction
## Tricycle assigns cell cycle along a radial trajectory between 0 and 2pi
tricycle.sce <- project_cycle_space(hepg2.sce, gname.type = "SYMBOL", species = "human")
tricycle.sce <- estimate_cycle_position(tricycle.sce)

## TOP2a is a sanity check as it should peak at a specific part of the cell cycle that we can confirm is consistent with estimation
top2a.idx <- which(rowData(tricycle.sce)$Symbol == "TOP2A")
fit.l <- fit_periodic_loess(tricycle.sce$tricyclePosition,
                            assay(tricycle.sce, 'logcounts')[top2a.idx, ],
                            plot = TRUE,
                            x_lab = "Cell cycle position \u03b8", 
                            y_lab = "log2(TOP2A)",
                            fig.title = paste0("Expression of TOP2A along \u03b8 (n=", 
                                               ncol(tricycle.sce), 
                                               ")"))
fit.l$fig + theme_bw(base_size = 14)
## TOP2A should peak at pi - it does?

## From the tricycle guide: users can approximately relate:
## 0.5pi to be the start of S stage, 
## pi to be the start of G2M stage, 
## 1.5pi to be the middle of M stage
## and 1.75pi-0.25pi to be G1/G0 stage

## Make that into some kind of maths
score.tricycle <- function(x){
  ifelse(x >= 0.5*pi & x < pi, "S", 
         ifelse(x >= pi & x < 1.75*pi, "G2M", 
                ifelse(x >= 1.75*pi | x < 0.25*pi, "G1", "NA")))
}

hepg2.seurat <- readRDS("data/hepg2_seurat_start.rds")

hepg2.cycle.seurat <- NormalizeData(hepg2.seurat)
hepg2.cycle.seurat <- CellCycleScoring(hepg2.cycle.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

cell.cycle.df <- data.frame("Cell.ID" = colnames(hepg2.sce),
                            "Seurat" = hepg2.cycle.seurat$Phase,
                            "Tricycle" = score.tricycle(tricycle.sce$tricyclePosition))

cell.cycle.melt <- make_long(cell.cycle.df, Seurat, Tricycle)

ggplot(cell.cycle.melt, 
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
  labs(x = element_blank(), title = "Cell Cycle Prediction") +
  scale_fill_manual(values = hues::iwanthue(4))

ggsave("plots/qc/cell_cycle_comparison.pdf", width = 5.8, height = 5.8)
ggsave("plots/qc/cell_cycle_comparison.png", width = 5.8, height = 5.8)

## Both results roughly agree
## Tricycle is newest so maybe we go with that?

## Add cell cycle results back into the object
hepg2.sce$Tricycle.Phase <- score.tricycle(tricycle.sce$tricyclePosition)
hepg2.sce$Tricycle.Position <- tricycle.sce$tricyclePosition

hepg2.sce$Seurat.Phase <- hepg2.cycle.seurat$Phase
hepg2.sce$Seurat.S <- hepg2.cycle.seurat$S.Score
hepg2.sce$Seurat.G2M <- hepg2.cycle.seurat$G2M.Score

hepg2.seurat$Seurat.Phase <- hepg2.cycle.seurat$Phase
hepg2.seurat$Seurat.S <- hepg2.cycle.seurat$S.Score
hepg2.seurat$Seurat.G2M <- hepg2.cycle.seurat$G2M.Score

hepg2.seurat$Tricycle.Phase <- score.tricycle(tricycle.sce$tricyclePosition)
hepg2.seurat$Tricycle.Position <- tricycle.sce$tricyclePosition

rm(tricycle.sce, hepg2.cycle.seurat)
gc()

saveRDS(hepg2.sce, "data/hepg2_sce_working.rds")
saveRDS(hepg2.seurat, "data/hepg2_seurat_working.rds")

## [ QC Filtering ] ----

hepg2.sce <- readRDS("data/hepg2_sce_working.rds")
hepg2.seurat <- readRDS("data/hepg2_seurat_working.rds")

colnames(colData(hepg2.sce))
colnames(rowData(hepg2.sce))

table(hepg2.sce$Class)
## All singlet

table(hepg2.sce$sample_id, hepg2.sce$Condition)

sample.cols <- iwanthue(length(names(table(hepg2.sce$sample_id))))
names(sample.cols) <- names(table(hepg2.sce$sample_id))

det.sce.gg <- plotColData(hepg2.sce, y = "detected", x = "sample_id", colour_by = "sample_id") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = sample.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 1000, linetype = "dotted")

mt.sce.gg <- plotColData(hepg2.sce, y = "subsets_mt_percent", x = "sample_id", colour_by = "sample_id") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = sample.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 10, linetype = "dotted")

rb.sce.gg <- plotColData(hepg2.sce, y = "subsets_ribo_percent", x = "sample_id", colour_by = "sample_id") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = sample.cols) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  geom_hline(yintercept = 5, linetype = "dotted")

det.sce.gg + mt.sce.gg + rb.sce.gg & theme(aspect.ratio = 0.5)

ggsave("plots/qc/hepg2_QC_plots_sample.pdf", width = 8.7, height = 5.8)
ggsave("plots/qc/hepg2_QC_plots_sample.png", width = 8.7, height = 5.8)

## Set the expression threshold over 1000
detected.filt <- colnames(hepg2.sce)[hepg2.sce$detected > 1000] ## Doesn't remove any cells
## Genes have to have at least 5 reads 
rowcounts.filt <- rownames(hepg2.sce)[Matrix::rowSums(counts(hepg2.sce)) > 5]
## Mitochondrial genes filtering
mt.genes.filt <- grep("^MT-", rownames(hepg2.sce), invert = TRUE, value = TRUE) ## Doesn't remove any genes
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

hepg2.filt.sce <- hepg2.sce[genes.filt, detected.filt]

dim(hepg2.sce) - dim(hepg2.filt.sce)
## We lose ~ 3900 genes

mito.filt <- colnames(hepg2.filt.sce)[hepg2.filt.sce$subsets_mt_percent < 15]
ribo.filt <- colnames(hepg2.filt.sce)[hepg2.filt.sce$subsets_ribo_percent > 5]
qc.filt <- intersect(mito.filt, ribo.filt)

hepg2.filt.sce <- hepg2.filt.sce[, qc.filt]

dim(hepg2.sce) - dim(hepg2.filt.sce)
## We lose ~ 5000 cells

## Percentage of the library that we filter out
round(100 - (100 * (table(hepg2.filt.sce$sample_id) / table(hepg2.sce$sample_id))), 1)
## HepG2_sample1 HepG2_sample10 HepG2_sample11 HepG2_sample12 HepG2_sample13 HepG2_sample14  HepG2_sample2  HepG2_sample3  HepG2_sample4  HepG2_sample5 
## 4.1           18.3           20.6            6.0            7.6           18.1            6.0           32.3           22.5            6.1 
## HepG2_sample6  HepG2_sample7  HepG2_sample8  HepG2_sample9 
## 8.7            6.8           21.4            7.0 

table(hepg2.filt.sce$sample_id)
## HepG2_sample1 HepG2_sample10 HepG2_sample11 HepG2_sample12 HepG2_sample13 HepG2_sample14  HepG2_sample2  HepG2_sample3  HepG2_sample4  HepG2_sample5 
## 6322           2018           1933           5267           1920           3067           4212           1094            922           2141 
## HepG2_sample6  HepG2_sample7  HepG2_sample8  HepG2_sample9 
## 2203           4208           2166           5544 

## The cisplatin recovery 1 and 2 samples got hit the hardest
## But relatively evenly hit looking at the percentages
round(100 - (100 * (table(hepg2.filt.sce$Condition) / table(hepg2.sce$Condition))), 1)
## POT           Untreated       Cisplatin_off Cisplatin_recovery1 Cisplatin_recovery2 
## 4.6                 6.9                 6.7                18.9                24.8 

table(hepg2.filt.sce$Condition)
## POT           Untreated       Cisplatin_off Cisplatin_recovery1 Cisplatin_recovery2 
## 8463               10623               12731                7018                4182 

## We can do the same filtering in Seurat
rowcounts.filt <- rownames(hepg2.seurat)[Matrix::rowSums(GetAssayData(hepg2.seurat, slot = "counts", assay = "RNA")) > 5]
mt.genes.filt <- grep("^MT-", rownames(hepg2.seurat), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

## Seurat and Scater have their own methods for calling # of detected genes so these will disagree slightly
## To be consistent we'll filter on $detected rather than $nFeature_RNA
hepg2.seurat$nFeature_RNA == hepg2.seurat$detected

hepg2.filt.seurat <- subset(hepg2.seurat, 
                            cells = WhichCells(hepg2.seurat, 
                                               expression = subsets_mt_percent < 15 &
                                                 subsets_ribo_percent > 5 &
                                                 detected > 1000),
                            features = genes.filt)

## The dimensions are equal
dim(hepg2.filt.seurat) == dim(hepg2.filt.sce)

saveRDS(hepg2.filt.sce, "data/hepg2_sce_filtered.rds")
saveRDS(hepg2.filt.seurat, "data/hepg2_seurat_filtered.rds")

## [ Filter HVGs ] ----

hepg2.filt.seurat <- readRDS("data/hepg2_seurat_filtered.rds")

## Filter HVGs
hepg2.filt.seurat <- NormalizeData(hepg2.filt.seurat)
hepg2.filt.seurat <- FindVariableFeatures(hepg2.filt.seurat, selection.method = "vst", nfeatures = 1000)
LabelPoints(plot = VariableFeaturePlot(hepg2.filt.seurat, assay = "RNA"),
            points = head(VariableFeatures(hepg2.filt.seurat), 20), repel = TRUE)
ggsave("plots/norm/filter_hvgs.pdf", width = 8.3, height = 5.8)
ggsave("plots/norm/filter_hvgs.png", width = 8.3, height = 5.8)

## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(hepg2.filt.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(hepg2.filt.seurat), value = TRUE),
              grep("[\\.]", rownames(hepg2.filt.seurat), value = TRUE),
              grep("^LINC", rownames(hepg2.filt.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
VariableFeatures(hepg2.filt.seurat) <- VariableFeatures(hepg2.filt.seurat)[which(!VariableFeatures(hepg2.filt.seurat) %in% rm.genes)]
## Scale data and score cell cycle
hepg2.filt.seurat <- ScaleData(hepg2.filt.seurat, features = VariableFeatures(hepg2.filt.seurat))

tmp.seurat <- CellCycleScoring(object = hepg2.filt.seurat, 
                               g2m.features = cc.genes$g2m.genes, 
                               s.features = cc.genes$s.genes)
tmp.seurat$Cycle.Score <- tmp.seurat$S.Score - tmp.seurat$G2M.Score

seurat.cycle.melt <- make_long(
  data.frame("Sample" = hepg2.filt.seurat$Sample, 
             "SCTransform" = tmp.seurat$Phase,
             "Normal" = hepg2.filt.seurat$Seurat.Phase),
  SCTransform, Normal)

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

DefaultAssay(hepg2.filt.seurat) <- "RNA"
hepg2.filt.seurat$Seurat.Cycle.Score <- hepg2.filt.seurat$Seurat.S - hepg2.filt.seurat$Seurat.G2M

saveRDS(hepg2.filt.seurat, "data/hepg2_seurat_hvgs.rds")

hepg2.seurat <- readRDS("data/hepg2_seurat_hvgs.rds")

hepg2.seurat <- RunPCA(hepg2.seurat, verbose = FALSE, npcs = 100, 
                       features = VariableFeatures(hepg2.seurat))
ElbowPlot(object = hepg2.seurat, ndims = 50, reduction = "pca")

## Integrate with Harmony and generate two separate UMAP dim reds +/-
harmony.seurat <- RunHarmony(hepg2.seurat, group.by.vars = "Sample",
                             theta = 1, lambda = 1, sigma = 0.07,
                             assay.use = "RNA", reduction = "pca",
                             dims.use = 1:30, reduction.save = "Harmony",
                             max_iter = 10, plot_convergence = FALSE)

harmony.seurat <- RunUMAP(harmony.seurat, reduction = "Harmony", dims = 1:30)
harmony.seurat <- FindNeighbors(harmony.seurat, reduction = "Harmony", dims = 1:30)

harmony.seurat$Condition <- factor(harmony.seurat$Condition,
                                   levels = c("POT",
                                              "Untreated",
                                              "Cisplatin_off",
                                              "Cisplatin_recovery1",
                                              "Cisplatin_recovery2"))

umap.gg <- DimPlot(harmony.seurat, group.by = "Condition", order = TRUE) +
  scale_colour_manual(values = group.cols) +
  umap.theme() + labs(title = "Conditions")
umap.gg
ggsave("plots/norm/umap_harmony_t1_l1_s007.png", width = 8.3, height = 5.8)
ggsave("plots/norm/umap_harmony_t1_l1_s007.pdf", width = 8.3, height = 5.8)

# hepg2.seurat <- RunUMAP(hepg2.seurat, reduction = "pca", dims = 1:30)
# hepg2.seurat <- FindNeighbors(hepg2.seurat, reduction = "pca", dims = 1:30)

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
saveRDS(harmony.seurat, "data/hepg2_seurat_hvgs_harmony.rds")

## [ Jaccard Similarity ] ----

## Split data
hepg2.seurat <- readRDS("data/hepg2_seurat_hvgs.rds")
obj.list <- SplitObject(hepg2.seurat, split.by = "Condition")

## Rerun variable feature selection
obj.list <- lapply(obj.list, FindVariableFeatures)
## Pick which genes to remove from HVG
rm.genes <- c(grep("^MT-", rownames(hepg2.seurat), value = TRUE),
              grep("^M?RP[SL]", rownames(hepg2.seurat), value = TRUE),
              grep("[\\.]", rownames(hepg2.seurat), value = TRUE),
              grep("^LINC", rownames(hepg2.seurat), value = TRUE),
              c("MALAT1"))
## Set back the genes you want to keep
output <- list()
for (i in 1:length(obj.list)) {
  obj <- obj.list[[i]]
  VariableFeatures(obj) <- VariableFeatures(obj)[which(!VariableFeatures(obj) %in% rm.genes)]
  output[[i]] <- obj
}
names(output) <- names(obj.list)
rm(hepg2.seurat, obj.list, obj)

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
saveRDS(output, "data/hepg2_seurat_split_list.rds")

hepg2.list <- readRDS("data/hepg2_seurat_split_list.rds")
hepg2.seurat <- readRDS("data/hepg2_seurat_hvgs_harmony.rds")
## Look through UMAPs and choose best cluster resolution
Idents(hepg2.list$POT) <- hepg2.list$POT$seurat_clusters.0.4
Idents(hepg2.list$Untreated) <- hepg2.list$Untreated$seurat_clusters.0.4
Idents(hepg2.list$Cisplatin_off) <- hepg2.list$Cisplatin_off$seurat_clusters.0.2
Idents(hepg2.list$Cisplatin_recovery1) <- hepg2.list$Cisplatin_recovery1$seurat_clusters.0.2
Idents(hepg2.list$Cisplatin_recovery2) <- hepg2.list$Cisplatin_recovery2$seurat_clusters.0.2
Idents(hepg2.seurat) <- hepg2.seurat$seurat_clusters.0.4

## Get cluster markers - MAST takes longer to run
pot.markers <- FindAllMarkers(hepg2.list$POT, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
ut.markers <- FindAllMarkers(hepg2.list$Untreated, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
cis.markers <- FindAllMarkers(hepg2.list$Cisplatin_off, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
rec1.markers <- FindAllMarkers(hepg2.list$Cisplatin_recovery1, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
rec2.markers <- FindAllMarkers(hepg2.list$Cisplatin_recovery2, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
hepg2.all.markers <- FindAllMarkers(hepg2.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)

dir.create("data/jaccard/", recursive = TRUE)
saveRDS(pot.markers, "data/jaccard/pot_markers.rds")
saveRDS(ut.markers, "data/jaccard/ut_markers.rds")
saveRDS(cis.markers, "data/jaccard/cis_markers.rds")
saveRDS(rec1.markers, "data/jaccard/rec1_markers.rds")
saveRDS(rec2.markers, "data/jaccard/rec2_markers.rds")
saveRDS(hepg2.all.markers, "data/jaccard/hepg2_all_markers.rds")

## Filter significant and higher expressed markers
pot.markers <- pot.markers[pot.markers$p_val_adj < 0.05, ]
ut.markers <- ut.markers[ut.markers$p_val_adj < 0.05, ]
cis.markers <- cis.markers[cis.markers$p_val_adj < 0.05, ]
rec1.markers <- rec1.markers[rec1.markers$p_val_adj < 0.05, ]
rec2.markers <- rec2.markers[rec2.markers$p_val_adj < 0.05, ]
hepg2.all.markers <- hepg2.all.markers[hepg2.all.markers$p_val_adj < 0.05, ]

pot.markers <- pot.markers[pot.markers$avg_log2FC > 0.5, ]
ut.markers <- ut.markers[ut.markers$avg_log2FC > 0.5, ]
cis.markers <- cis.markers[cis.markers$avg_log2FC > 0.5, ]
rec1.markers <- rec1.markers[rec1.markers$avg_log2FC > 0.5, ]
rec2.markers <- rec2.markers[rec2.markers$avg_log2FC > 0.5, ]
hepg2.all.markers <- hepg2.all.markers[hepg2.all.markers$avg_log2FC > 0.5, ]
## Split into lists of genes per cluster
pot.markers <- split(pot.markers, pot.markers$cluster)
ut.markers <- split(ut.markers, ut.markers$cluster)
cis.markers <- split(cis.markers, cis.markers$cluster)
rec1.markers <- split(rec1.markers, rec1.markers$cluster)
rec2.markers <- split(rec2.markers, rec2.markers$cluster)
hepg2.all.markers <- split(hepg2.all.markers, hepg2.all.markers$cluster)

## Check how many markers you get per cluster and change the number you input to the comparison
lapply(pot.markers, nrow) ## 0 = 176, 1 = 425, 2 = 994, 3 = 71, 4 = 232, 5 = 551
lapply(ut.markers, nrow) ## 0 = 28, 1 = 465, 2 = 238, 3 = 254, 4 = 486, 5 = 1012, 6 = 349
lapply(cis.markers, nrow) ## 0 = 523, 1 = 596, 2 = 734, 3 = 2013
lapply(rec1.markers, nrow) ## 0 = 313, 1 = 619, 2 = 437
lapply(rec2.markers, nrow) ## 0 = 591, 1 = 1204, 2 = 370, 3 = 791
lapply(hepg2.all.markers, nrow) ## 0 = 378, 1 = 305, 2 = 1353, 3 = 634, 4 = 588, 5 = 94, 6 = 1044, 7 = 804

pot.markers <- lapply(pot.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 100)
})
ut.markers <- lapply(ut.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 100)
})
cis.markers <- lapply(cis.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 500)
})
rec1.markers <- lapply(rec1.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
rec2.markers <- lapply(rec2.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 300)
})
hepg2.all.markers <- lapply(hepg2.all.markers, \(x) {
  x <- x[order(x$avg_log2FC, decreasing = TRUE), ]
  head(x, 100)
})

## Give your clusters unique names
names(pot.markers) <- paste0("POT_", names(pot.markers))
names(ut.markers) <- paste0("UT_", names(ut.markers))
names(cis.markers) <- paste0("Cis_", names(cis.markers))
names(rec1.markers) <- paste0("Rec1_", names(rec1.markers))
names(rec2.markers) <- paste0("Rec2_", names(rec2.markers))
names(hepg2.all.markers) <- paste0("hepg2_all_", names(hepg2.all.markers))
## Make a big list of all the cluster markers
hepg2.markers.list <- c(pot.markers,
                        ut.markers,
                        cis.markers,
                        rec1.markers,
                        rec2.markers,
                        hepg2.all.markers)

## Define your similarity function, you can use another but Jaccard works well for this
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection/union)
}
## Set up the matrix you will populate with data
hepg2.jaccard.mat <- matrix(data = NA, nrow = length(hepg2.markers.list),
                            ncol = length(hepg2.markers.list),
                            dimnames = list(names(hepg2.markers.list),
                                            names(hepg2.markers.list)))
## Run your pairwise Jaccard similarity
for (i in rownames(hepg2.jaccard.mat)) {
  for (j in colnames(hepg2.jaccard.mat)) {
    hepg2.jaccard.mat[i,j] <- jaccard(hepg2.markers.list[[i]]$gene, hepg2.markers.list[[j]]$gene)
  }
}

rownames(hepg2.jaccard.mat) <- gsub("hepg2_", "", rownames(hepg2.jaccard.mat))
colnames(hepg2.jaccard.mat) <- gsub("hepg2_", "", colnames(hepg2.jaccard.mat))
rownames(hepg2.jaccard.mat) <- gsub("all", "All", rownames(hepg2.jaccard.mat))
colnames(hepg2.jaccard.mat) <- gsub("all", "All", colnames(hepg2.jaccard.mat))
anno.df <- data.frame("Timepoint" = gsub("_.+", "", colnames(hepg2.jaccard.mat)),
                      row.names = colnames(hepg2.jaccard.mat))
anno.col <- list("Timepoint" = c("All" = "black",
                                 "POT" = as.character(group.cols[1]),
                                 "UT" = as.character(group.cols[2]),
                                 "Cis" = as.character(group.cols[3]),
                                 "Rec1" = as.character(group.cols[4]),
                                 "Rec2" = as.character(group.cols[5])))
dev.off()
gt <- pheatmap(hepg2.jaccard.mat,
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
ggsave("plots/jaccard/jaccard_heatmap.png", plot = gt, width = 8.3, height = 8.3)
ggsave("plots/jaccard/jaccard_heatmap.pdf", plot = gt, width = 8.3, height = 8.3)

## /////////////////////////////////////////////////////////////////////////////
## PAM clustering //////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

## Run the parameter search for PAM
hepg2.jaccard.dist <- as.dist(1 - hepg2.jaccard.mat)
silhouette.res <- numeric()
## Run PAM for 2-15 clusters and see what gives you the highest silhouette
for (k in 2:15) {  # Assuming you want to check from 2 to 15 clusters
  pam.fit <- pam(hepg2.jaccard.dist, k, diss = TRUE)
  silhouette.res[k] <- mean(silhouette(pam.fit)[,"sil_width"])
}
silhouette.res <- silhouette.res[-1]
names(silhouette.res) <- 2:15
names(which.max(silhouette.res)) ## gives you the k that is best for your data
## For this data it was 15
silhouette.df <- as.data.frame(silhouette.res)
silhouette.df$k <- rownames(silhouette.df)
colnames(silhouette.df)[1] <- "silhouette"
silhouette.df$k <- factor(silhouette.df$k, levels = c(2:16))
## Draw the geom_vline at 11 for this data
ggplot(silhouette.df, mapping = aes(x = k, y = silhouette, group = 1)) +
  geom_line() +
  geom_vline(xintercept = 14, colour = "red", linetype = "solid", linewidth = 0.8) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major.x = element_line(linetype = "dotted"),
        panel.grid.minor.y = element_line(linetype = "dotted"),
        panel.grid.major.y = element_line(linetype = "dotted")) +
  labs(x = "k", y = "Silhouette", title = "PAM clustering silhouette score")
ggsave("plots/jaccard/pam_clustering_silhouette_score.png", width = 8.3, height = 5.8)
ggsave("plots/jaccard/pam_clustering_silhouette_score.pdf", width = 8.3, height = 5.8)

## Run the PAM and get the cluster assignment
hepg2.k15.pam <- pam(hepg2.jaccard.dist, k = 15, diss = TRUE, cluster.only = TRUE)
anno.df$PAM.Cluster <- hepg2.k15.pam
## Assign the clusters back to your Seurat objects
hepg2.clusters.df <- anno.df$PAM.Cluster
names(hepg2.clusters.df) <- rownames(anno.df)
pot.set <- paste0("POT_", as.character(hepg2.list$POT$seurat_clusters.0.4))
ut.set <- paste0("UT_", as.character(hepg2.list$Untreated$seurat_clusters.0.4))
cis.set <- paste0("Cis_", as.character(hepg2.list$Cisplatin_off$seurat_clusters.0.2))
rec1.set <- paste0("Rec1_", as.character(hepg2.list$Cisplatin_recovery1$seurat_clusters.0.2))
rec2.set <- paste0("Rec2_", as.character(hepg2.list$Cisplatin_recovery2$seurat_clusters.0.2))
all.set <- paste0("All_", as.character(hepg2.seurat$seurat_clusters.0.4))

hepg2.list$POT$PAM.Cluster <- factor(as.character(hepg2.clusters.df[pot.set]))
hepg2.list$Untreated$PAM.Cluster <- factor(as.character(hepg2.clusters.df[ut.set]))
hepg2.list$Cisplatin_off$PAM.Cluster <- factor(as.character(hepg2.clusters.df[cis.set]))
hepg2.list$Cisplatin_recovery1$PAM.Cluster <- factor(as.character(hepg2.clusters.df[rec1.set]))
hepg2.list$Cisplatin_recovery2$PAM.Cluster <- factor(as.character(hepg2.clusters.df[rec2.set]))
hepg2.seurat$PAM.Cluster <- factor(as.character(hepg2.clusters.df[all.set]))
saveRDS(hepg2.list, "data/hepg2_seurat_split_list_pam.rds")
saveRDS(hepg2.seurat, "data/hepg2_seurat_hvgs_pam.rds")

## [ Markers ] ----

## PAM cluster markers
Idents(hepg2.seurat) <- hepg2.seurat$PAM.Cluster
hepg2.pam.markers <- FindAllMarkers(hepg2.seurat, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, densify = TRUE)
saveRDS(hepg2.pam.markers, "data/jaccard/hepg2_pam_markers.rds")

table(hepg2.pam.markers$cluster)
## 1   14   15    2    3    4    6 
## 4833 4854 3716 4685 5949 3892 4425 

hepg2.pam.list <- lapply(unique(hepg2.pam.markers$cluster), function(x){
  tmp.df <- hepg2.pam.markers[hepg2.pam.markers$cluster==x & hepg2.pam.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})
names(hepg2.pam.list) <- unique(hepg2.pam.markers$cluster)

saveRDS(hepg2.pam.list, "data/jaccard/hepg2_pam_markers_list.rds")

hepg2.pam.GO.BP.list <- lapply(hepg2.pam.list, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

hepg2.pam.C2.list <- lapply(hepg2.pam.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

hepg2.pam.H.list <- lapply(hepg2.pam.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

saveRDS(hepg2.pam.GO.BP.list, "data/jaccard/hepg2_pam_markers_BP_ontology.rds")
saveRDS(hepg2.pam.C2.list, "data/jaccard/hepg2_pam_markers_C2_geneset.rds")
saveRDS(hepg2.pam.H.list, "data/jaccard/hepg2_pam_markers_Hallmarks.rds")

hepg2.pam.lists <- list("hepg2.pam.GO.BP.list" = hepg2.pam.GO.BP.list,
                        "hepg2.pam.C2.list" = hepg2.pam.C2.list,
                        "hepg2.pam.H.list" = hepg2.pam.H.list,
                        "hepg2.pam.list" = hepg2.pam.list)

hepg2.pam.names <- sub("list", "sheets", names(hepg2.pam.lists))

for(i in 1:length(hepg2.pam.lists)){
  hepg2.pam.list <- hepg2.pam.lists[[i]]
  list.sheets <- lapply(hepg2.pam.list, as.data.frame)
  assign(hepg2.pam.names[[i]], list.sheets)
}

write_xlsx(hepg2.pam.GO.BP.sheets, "data/jaccard/hepg2_pam_markers_BP_ontology.xlsx")
write_xlsx(hepg2.pam.C2.sheets, "data/jaccard/hepg2_pam_markers_C2_geneset.xlsx")
write_xlsx(hepg2.pam.H.sheets, "data/jaccard/hepg2_pam_markers_Hallmarks.xlsx")
write_xlsx(hepg2.pam.sheets, "data/jaccard/hepg2_pam_markers.xlsx")

## [ Cluster annotation ] ----

hepg2.seurat <- readRDS("data/hepg2_seurat_hvgs_pam_anno.rds")

DimPlot(hepg2.seurat, group.by = "Condition", order = TRUE) +
  umap.theme() + labs(title = "Condition") +
  scale_colour_manual(values = group.cols, labels = conditions.list) +
  theme(text = element_text(size = 15))
ggsave("plots/umaps/umap_harmony_conditions.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/umap_harmony_conditions.png", width = 8.3, height = 5.8)

## Annotate PAM clusters
hepg2.seurat$Exprs_clust <- recode(as.character(hepg2.seurat$PAM.Cluster),
                                   "1" = "2",
                                   "2" = "3",
                                   "3" = "6",
                                   "4" = "7",
                                   "6" = "1",
                                   "14" = "5",
                                   "15" = "4")

DimPlot(hepg2.seurat, group.by = "Exprs_clust", order = TRUE) +
  umap.theme() + labs(title = "Expression clusters") +
  scale_colour_manual(values = clust.cols) +
  theme(text = element_text(size = 15))
ggsave("plots/umaps/umap_harmony_pam_nums.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/umap_harmony_pam_nums.png", width = 8.3, height = 5.8)

anno.data <- hepg2.seurat@meta.data[c("Condition", "Exprs_clust")]
anno.table <- as.data.frame(table(anno.data$Condition, anno.data$Exprs_clust))
ggplot(anno.table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("Cluster proportion") +
  guides(fill = guide_legend(title = "Expression clusters")) +
  scale_x_discrete(labels = conditions.list) +
  scale_fill_manual(values = clust.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/exprs_clust_proportion.pdf", width = 8.3, height = 5.8)
ggsave("plots/exprs_clust_proportion.png", width = 8.3, height = 5.8)

hepg2.seurat$Exprs_state <- recode(as.character(hepg2.seurat$Exprs_clust),
                                   "1" = "Hep-like",
                                   "2" = "Hep-like",
                                   "3" = "Hep-like",
                                   "4" = "Hep-like",
                                   "5" = "Pro-like",
                                   "6" = "Pro-like",
                                   "7" = "Stem-like")
saveRDS(hepg2.seurat, "data/hepg2_seurat_hvgs_pam_states.rds")

DimPlot(hepg2.seurat, group.by = "Exprs_state", order = TRUE) +
  umap.theme() + labs(title = "Expression states") +
  scale_colour_manual(values = state.cols) +
  theme(text = element_text(size = 15))
ggsave("plots/umaps/umap_harmony_states.pdf", width = 8.3, height = 5.8)
ggsave("plots/umaps/umap_harmony_states.png", width = 8.3, height = 5.8)

anno.data <- hepg2.seurat@meta.data[c("description", "Condition", "Exprs_state")]
anno.table <- as.data.frame(table(anno.data$Condition, anno.data$Exprs_state))
ggplot(anno.table, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("Cluster proportion") +
  guides(fill = guide_legend(title = "Expression state")) +
  scale_x_discrete(labels = conditions.list) +
  scale_fill_manual(values = state.cols) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/exprs_state_proportion.pdf", width = 8.3, height = 5.8)
ggsave("plots/exprs_state_proportion.png", width = 8.3, height = 5.8)

## /////////////////////////////////////////////////////////////////////////////
## Annotation heatmaps /////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

hepg2.seurat <- readRDS("data/hepg2_seurat_hvgs_pam_states.rds")
pseudo.hepg2 <- AggregateExpression(hepg2.seurat, assays = "RNA", return.seurat = T,
                                    group.by = c("Exprs_clust"))
Idents(pseudo.hepg2) <- "Exprs_clust"

## HB sigs
sig.df <- as.data.frame(read_xlsx(path = "../huh6/data/hb_sigs_filt.xlsx", col_names = FALSE))
colnames(sig.df) <- c("Gene", "Signature")
sigs <- lapply(unique(sig.df$Signature), function(x){
  sig.df[sig.df$Signature==x, "Gene"]
})
names(sigs) <- unique(sig.df$Signature)
sigs
sig.names <- paste(names(sigs), ".Sig", sep = "")
sig.names <- gsub("_", "-", sig.names)

pseudo.hepg2 <- AddModuleScore(pseudo.hepg2, features = sigs, assay = "RNA", seed = 12345, 
                               name = sig.names)
colnames(pseudo.hepg2@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.hepg2@meta.data))
pseudo.hepg2[["HB_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.hepg2, vars = sig.names)))

hepg2.mat <- as.matrix(pseudo.hepg2@assays$HB_sigs@data)
anno.df <- data.frame(Exprs_clust = 1:7)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(hepg2.mat)
anno.cols <- clust.cols
gt <- pheatmap(hepg2.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:7,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggplot2::ggsave("plots/heatmaps/pseudo_pam_hb_pheatmap.pdf", plot = gt)
ggplot2::ggsave("plots/heatmaps/pseudo_pam_hb_pheatmap.png", plot = gt)

## Development markers (Wesley et al., 2022)
wesley.df <- as.data.frame(read_xlsx(path = "../huh6/data/welsey2022_sig.xlsx", col_names = FALSE))
colnames(wesley.df) <- c("Gene", "Signature")
wesley.sig <- lapply(unique(wesley.df$Signature), function(x){
  wesley.df[wesley.df$Signature==x, "Gene"]
})
names(wesley.sig) <- unique(wesley.df$Signature)
wesley.sig
wesley.names <- paste(names(wesley.sig), ".Sig", sep = "")

pseudo.hepg2 <- AddModuleScore(pseudo.hepg2, features = wesley.sig, assay = "RNA", seed = 12345, 
                               name = wesley.names)
colnames(pseudo.hepg2@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.hepg2@meta.data))
pseudo.hepg2[["Wesley_sigs"]] <- CreateAssayObject(data = t(FetchData(object = pseudo.hepg2, vars = wesley.names)))

hepg2.wes.mat <- as.matrix(pseudo.hepg2@assays$Wesley_sigs@data)
anno.df <- data.frame(Exprs_clust = 1:7)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(hepg2.wes.mat)
# anno.cols <- clust.cols
gt <- pheatmap(hepg2.wes.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:7,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggplot2::ggsave("plots/heatmaps/pseudo_pam_wesley_pheatmap.pdf", plot = gt)
ggplot2::ggsave("plots/heatmaps/pseudo_pam_wesley_pheatmap.png", plot = gt)

## Liver development markers
custom <- c("PROM1", "NCAM1", "EPCAM", "CLDN3", ## stemmy
            "ICAM1", "KRT8", "KRT18", "KRT19", ## progenitor
            "ALB", "AFP", "CEBPA", "HNF4A", "MAT1A") ## hepatocyte
custom %in% rownames(pseudo.hepg2)

hepg2.custom.mat <- as.matrix(pseudo.hepg2@assays$RNA$scale.data[custom, ])
anno.df <- data.frame(Exprs_clust = 1:7)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(hepg2.custom.mat)
# anno.cols <- clust.cols
hc <- hclust(dist(hepg2.custom.mat), method = "ward.D2")
library(dendextend)
dend <- as.dendrogram(hc)
dend <- rotate(dend, order = c("AFP", "ALB", "CEBPA", "HNF4A", "MAT1A", "EPCAM", "PROM1", "CLDN3", "NCAM1", "ICAM1", "KRT19", "KRT8", "KRT18"))
## Find best order that corresponds to marker programs for ease of visualisation
## Does not force order; works within constraints of tree topology
gt <- pheatmap(hepg2.custom.mat,
               border_color = NA,
               cellwidth = 8, cellheight = 8,
               fontsize_row = 8, fontsize_col = 8,
               cluster_rows = as.hclust(dend),
               color = Seurat:::SpatialColors(100),
               labels_col = 1:7,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggplot2::ggsave("plots/heatmaps/pseudo_pam_custom_pheatmap.pdf", plot = gt)
ggplot2::ggsave("plots/heatmaps/pseudo_pam_custom_pheatmap.png", plot = gt)

saveRDS(pseudo.hepg2, "data/hepg2_seurat_pseudobulk_exprs_clust_sigs.rds")

## HuH6-derived signature
unique.lists <- readRDS("../huh6/data/huh6_temp_sig.rds")
pseudo.hepg2 <- AddModuleScore(pseudo.hepg2, features = unique.lists, name = paste0(names(unique.lists), ".Sig"),
                               assay = "RNA", seed = 12345)
colnames(pseudo.hepg2@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(pseudo.hepg2@meta.data))
pseudo.hepg2[["Markers_sigs"]] <- CreateAssayObject(
  data = t(FetchData(object = pseudo.hepg2,
                     vars = c("Hep.Sig", "Prog.Sig", "Stem.Sig"))))

hepg2.marker.mat <- as.matrix(pseudo.hepg2@assays$Markers_sigs@data)
anno.df <- data.frame(Exprs_clust = 1:7)
anno.df$Exprs_clust <- factor(anno.df$Exprs_clust)
rownames(anno.df) <- colnames(hepg2.marker.mat)
# anno.cols <- clust.cols
gt <- pheatmap(hepg2.marker.mat,
               border_color = NA,
               cellwidth = 12, cellheight = 12,
               fontsize_row = 10, fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2",
               color = Seurat:::SpatialColors(100),
               labels_col = 1:7,
               annotation_col = anno.df,
               annotation_colors = list(Exprs_clust = anno.cols))$gtable
ggplot2::ggsave("plots/heatmaps/pseudo_pam_marker_pheatmap.pdf", plot = gt)
ggplot2::ggsave("plots/heatmaps/pseudo_pam_marker_pheatmap.png", plot = gt)

## [ Entropy ] ----

hepg2.seurat <- readRDS("data/hepg2_seurat_hvgs_pam_anno.rds")
hepg2.sce <- as.SingleCellExperiment(hepg2.seurat)
saveRDS(hepg2.sce, "data/hepg2_sce_hvgs_pam_anno.rds")

dir.create("plots/entropy/", recursive = TRUE)

entropy <- perCellEntropy(hepg2.sce)
entropy.df <- data.frame(cluster = hepg2.sce$PAM_Name, entropy = entropy)

ggplot(entropy.df, aes(x = cluster, y = entropy, fill = cluster)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size = 12)) +
  scale_fill_manual(values = pam.cols) +
  labs(x = "Cluster", y = "Entropy")
ggsave("plots/entropy/entropy_cluster_violin.pdf", width = 8.3, height = 5.8)
ggsave("plots/entropy/entropy_cluster_violin.png", width = 8.3, height = 5.8)

## [ Trajectory inference ] ----

## Trajectory inference was performed using Slingshot
hepg2.slingshot.sce <- slingshot(hepg2.sce, reducedDim = "PCA",
                                 clusterLabels = "PAM_Name", start.clus = "Late progenitor",
                                 stretch = 0)
saveRDS(hepg2.slingshot.sce, "data/hepg2_sce_slingshot.rds")

## Extract matrix of pseudotime values along each lineage
slingshot.lineages <- slingPseudotime(hepg2.slingshot.sce)
head(slingshot.lineages) ## 3 lineages

## Single pseudotime for all cells
shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

umap.gg <- plotUMAP(hepg2.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(hepg2.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                 aes(x = umap_1, y=umap_2), linewidth = 1.2)
}
umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))
ggsave("plots/slingshot/umap_slingshot.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/umap_slingshot.png", width = 8.3, height = 5.8)

saveRDS(slingshot.embedded.all, "data/slingshot_embedded_all.rds")

## Save as SDS
slingshot.embedded.all.sds <- embedCurves(hepg2.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
saveRDS(slingshot.embedded.all.sds, "data/slingshot_sds_embedded_all.rds")

## Plot pseudotime values for cells in each lineage
no.cols <- 2
pseudotime <- slingPseudotime(hepg2.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
pal <- viridis(100)
pdf("plots/slingshot//slingshot_pseudotime_lineages.pdf", width = 8.3, height = 4.8)
png("plots/slingshot/slingshot_pseudotime_lineages.png", width = 8.3, height = 4.8,
    units = "in", res = 200)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(hepg2.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}
dev.off()

## Embed lineages one by one
slingshot.embedded <- embedCurves(hepg2.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(hepg2.slingshot.sce, colour_by = "PAM_Name") +
  umap.theme() + scale_colour_manual(values = pam.cols) +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

dir.create("plots/slingshot/by_lineage/", recursive = TRUE)

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("data/slingshot_embedded_", lineage.names[[i]], ".rds"))
  
  plotUMAP(hepg2.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), size = 1.2) +
    labs(title = lineage.names[[i]]) + umap.theme() +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], ".pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], ".png"), width = 8.7, height = 5.8)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), size = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_cluster.pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_cluster.png"), width = 8.7, height = 5.8)
}

## [ Gene dynamics ] ----

slingshot.sce <- readRDS("data/slingshot_late_progenitor/hepg2_sce_slingshot.rds")
slingshot.cols <- grep("^slingPseudotime_", colnames(colData(slingshot.sce)), value = TRUE)
cluster.cols <- pam.cols
names(cluster.cols) <- names(pam.cols)

for (pseudotime in slingshot.cols) {
  lineage <- TSCAN::testPseudotime(slingshot.sce, pseudotime = slingshot.sce[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/hepg2_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/hepg2_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(slingshot.sce, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "PAM_Name") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(slingshot.sce[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(slingshot.sce[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "PAM_Name",
                                  column_annotation_colours = list(PAM_Name = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/hepg2_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/hepg2_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(slingshot.sce, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "PAM_Name") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(slingshot.sce[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "PAM_Name",
                                column_annotation_colours = list(PAM_Name = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/slingshot/expression/hepg2_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## [ Gene annotation ] ----

file.paths <- list.files(path = here::here("data/slingshot_late_progenitor/trajectory_genes/"), pattern = "\\.rds", full.names = TRUE)
file.paths <- file.paths[grep("slingpseudotime1", file.paths, ignore.case = TRUE)]
file.names <-  gsub("_genes.rds$", "",
                    gsub("hepg2_sling", "", x = basename(file.paths)))
dynamic.genes <- lapply(file.paths, readRDS)
names(dynamic.genes) <- file.names

## GSEA of dynamic genes for each lineage
dynamic.genes.GO.BP.list <- lapply(dynamic.genes, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

dynamic.genes.C2.list <- lapply(dynamic.genes, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

dynamic.genes.H.list <- lapply(dynamic.genes, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})
saveRDS(dynamic.genes.GO.BP.list, "data/slingshot_late_progenitor/trajectory_genes/hepg2_dynamic_genes_BP_ontology.rds")
saveRDS(dynamic.genes.C2.list, "data/slingshot_late_progenitor/trajectory_genes/hepg2_dynamic_genes_C2_geneset.rds")
saveRDS(dynamic.genes.H.list, "data/slingshot_late_progenitor/trajectory_genes/hepg2_dynamic_genes_Hallmarks.rds")

dynamic.genes.lists <- list("dynamic.genes.GO.BP.list" = dynamic.genes.GO.BP.list,
                            "dynamic.genes.C2.list" = dynamic.genes.C2.list,
                            "dynamic.genes.H.list" = dynamic.genes.H.list)
dynamic.genes.names <- sub("list", "sheets", names(dynamic.genes.lists))

for(i in 1:length(dynamic.genes.lists)){
  gs.list <- dynamic.genes.lists[[i]]
  list.sheets <- lapply(gs.list, as.data.frame)
  assign(dynamic.genes.names[[i]], list.sheets)
}
write_xlsx(dynamic.genes.GO.BP.sheets, "data/slingshot_late_progenitor/trajectory_genes/hepg2_dynamic_genes_BP_ontology.xlsx")
write_xlsx(dynamic.genes.C2.sheets, "data/slingshot_late_progenitor/trajectory_genes/hepg2_dynamic_genes_C2_geneset.xlsx")
write_xlsx(dynamic.genes.H.sheets, "data/slingshot_late_progenitor/trajectory_genes/hepg2_dynamic_genes_Hallmarks.xlsx")

## Plot GSEA results for hallmark and gene ontology biological processes gene sets
for (i in 1:length(dynamic.genes.GO.BP.list)) {
  dotplot(dynamic.genes.GO.BP.list[[i]], showCategory = 20) +
    ggtitle(paste0("GSEA ", gsub("_", " ", names(dynamic.genes.GO.BP.list)[i])))
  ggsave(paste0("plots/slingshot_late_progenitor/expression/hepg2_gobp_", names(dynamic.genes.GO.BP.list)[i], ".pdf"), width = 8.3, height = 8.3)
  ggsave(paste0("plots/slingshot_late_progenitor/expression/hepg2_gobp_", names(dynamic.genes.GO.BP.list)[i], ".png"), width = 8.3, height = 8.3)
}

for (i in 1:length(dynamic.genes.H.list)) {
  dotplot(dynamic.genes.H.list[[i]], showCategory = 20) +
    ggtitle(paste0("GSEA ", gsub("_", " ", names(dynamic.genes.H.list)[i])))
  ggsave(paste0("plots/slingshot_late_progenitor/expression/hepg2_hallmarks_", names(dynamic.genes.H.list)[i], ".pdf"), width = 8.3, height = 8.3)
  ggsave(paste0("plots/slingshot_late_progenitor/expression/hepg2_hallmarks_", names(dynamic.genes.H.list)[i], ".png"), width = 8.3, height = 8.3)
}

## /////////////////////////////////////////////////////////////////////////////
## TF enrichment ///////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

file.paths <- list.files(path = here::here("data/slingshot_late_progenitor/trajectory_genes/"), pattern = "\\.rds", full.names = TRUE)
file.paths <- file.paths[grep("slingpseudotime1", file.paths, ignore.case = TRUE)]
file.names <-  gsub("_genes.rds$", "",
                    gsub("huh6_sling", "", x = basename(file.paths)))
dynamic.genes <- lapply(file.paths, readRDS)
names(dynamic.genes) <- file.names
gene.lists <- lapply(dynamic.genes, `[[`, "gene")

data(motifAnnotations_hgnc)
# motifRankings <- importRankings("~/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
motifRankings <- importRankings("~/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

motif.enrichment <- cisTarget(gene.lists, motifRankings,
                              motifAnnot = motifAnnotations)
saveRDS(motif.enrichment, "data/slingshot_late_progenitor/trajectory_genes/hepg2_motif_enrichment_500bp.rds")
split.list <- split(motif.enrichment, motif.enrichment$geneSet)
write_xlsx(split.list, "data/slingshot_late_progenitor/trajectory_genes/hepg2_motif_enrichment_500bp.xlsx")

# results.subset <- motif.enrichment[motif.enrichment$geneSet == "pseudotime2_down", ]
# showLogo(results.subset)

