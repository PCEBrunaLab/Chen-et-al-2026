# Integrated Lineage Tracing in Hepatoblastoma Finds a Plasticity-Proliferation Axis That Drives Post-Treatment Adaptation

## Abstract
Paediatric cancers, with their developmental context and low mutational burden, are powerful models for studying plasticity in cancer evolution. Hepatoblastoma (HB) has among the lowest mutational burdens of childhood cancers, yet ~20% relapse, implicating non-genetic mechanisms in treatment adaptation. We integrate expressed DNA barcodes with longitudinal single-cell multiomics to trace clonal and phenotypic dynamics following cisplatin treatment in HB cell lines. Cisplatin imposes a phenotypic bottleneck, selecting for progenitor-like cells across clones without a clonal sweep. Exit from persistence involves asynchronous, partly fitness-biased re-expansion and re-diversification towards hepatocytic-like states, also captured in two longitudinal patient cases. Pseudotime and dynamical-landscape analyses reveal a plasticity-proliferation axis underlying recovery, where persister cells occupy an arrested-plastic state, exited to restart proliferation. Chromatin and regulatory analyses identify the E2F target BIRC5 (survivin) as a key node in this axis, epigenetically primed post-recovery. BIRC5 maintains hepatocytic-like states, and its knockdown shifts cells toward progenitor-like phenotypes, while pharmacological inhibition selectively eliminates cisplatin-tolerant persisters. These findings position BIRC5 as a regulator of plasticity-driven adaptation and a candidate vulnerability in the minimal-residual disease window.

## Analysis scripts

**HuH6 scRNA-seq analysis** 
   - Standard `Seurat` pre-processing, QC filtering, normalisation, scaling and dimensionality reduction
   - `Seurat` and partitioning around medoids (PAM) clustering
   - Trajectory inference using `Slingshot`
   - Transcritpion factor motif analysis of pseudotime-derived genes with `RcisTarget`
   - Integration with barcode expression data
   - `MuTrans` landscape analysis in `Python`

**HuH6 DNA barcode analysis**
   - Filtering and plotting of barcode abundance by counts in the DNA sequencing data
   - Calculation of beta diversity using `vegan`  

**HuH6 single-cell multiomics analysis**
   - snRNA-seq analysis with `Seurat`, analogous to as outlined above for the analysis of HuH6 scRNA-seq data
   - snATAC-seq analysis using `Signac` with MACS3 peak calling

**HepG2 scRNA-seq analysis**
   - (As outlined above for HuH6 scRNA-seq analysis)

**Patient snRNA-seq analysis**
   - Standard `Seurat` pre-processing, QC filtering, normalisation, scaling and dimensionality reduction
   - `Seurat` clustering
   - Cell type inference with `SingleR` annotation
   - Analysis of malignant cells using `InferCNV`

**BIRC5 siRNA HB303, HepG2, HuH6 scRNA-seq analysis**
   - Standard `Seurat` pre-processing, QC filtering, normalisation, scaling and dimensionality reduction
   - `Seurat` clustering
