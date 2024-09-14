library(Seurat)
library(qs)
library(STdeconvolve)
library(scuttle)
library(dplyr)
library(Matrix)

source("./src/util.R")
source("./config/config.R")

set.seed(1)

# load data as spatialExperiment object
spe.path <- get.std.path()
spe <- load.data(spe.path)

cd <- assay(spe, "counts")
dim(cd)

# QC: Filter out poor pixels and genes
counts <- cleanCounts(cd, min.lib.size = 1, min.reads = 1, verbose = TRUE)
spe.clean <- spe[,colnames(counts)]

# compute colsums to check lib size
spe.clean$lib.size <- colSums(counts(spe.clean))

# sum(lib.size <= 100)
hist(spe.clean$lib.size, breaks = 200)

# compute umap and  PCA reductions
## seurat
seurat_obj <- CreateSeuratObject(counts = counts) %>% 
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(features = VariableFeatures(.), npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:50)

umap_seurat <- Embeddings(seurat_obj, reduction = "umap")
pca_seurat <- Embeddings(seurat_obj, reduction = "pca")

reducedDim(spe.clean, "umap_seurat") <- umap_seurat
reducedDim(spe.clean, "pca_seurat") <- pca_seurat

## Huber paper method
logcounts(spe.clean) <- transformGamPoi::shifted_log_transform(spe.clean)
umap_huber <- scater::calculateUMAP(spe.clean, pca = 10)

reducedDim(spe.clean, "umap_huber", withDimnames = FALSE) <- umap_huber
reducedDim(spe.clean, "umap_sel", withDimnames = FALSE) <- umap_huber


# var.genes from seurat
all_genes <- rownames(spe.clean) 
variable_genes <- VariableFeatures(seurat_obj)

metadata(spe.clean)$var.genes <- variable_genes

# non var.genes for testing
nonzero.gene <- rownames(spe.clean)[rowSums(counts(spe.clean)) > 0]
gene.var <- apply(counts(spe.clean)[nonzero.gene,], 1, var) %>% sort()
# gene.var[1]
# sum(counts(spe.clean)["LDLRAD2",])
metadata(spe.clean)$nonvar.genes <- names(gene.var)[1:3] 

# select spatial features STdeconvolve
corpus <- restrictCorpus(counts, 
                         removeAbove=0.95, 
                         removeBelow = 0.05, 
                         nTopOD = NA,
                         plot = F,
                         alpha = 0.1) # 1e-8

corpus <- as.matrix(corpus)
dim(corpus)

# remove mitochondria genes => no  mitochondria genes
mask.mt <- startsWith(rownames(corpus), "MT-") | startsWith(rownames(corpus), "mt-")
mt_genes <- corpus[mask.mt, ]
length(mt_genes)
corpus <- corpus[!mask.mt, ]
dim(corpus)

# load reference data set
dir.data <- get.data.dir()
sce.path <- file.path(dir.data, "scRNA-seq", "sce.qs")
# qsave(sce, sce.path)
sce <- qread(sce.path)

shared.genes <- intersect(rownames(sce), rownames(spe))
length(shared.genes)
sce.subset <- sce[shared.genes,]

# QC
seurat_obj <- as.Seurat(sce.subset, 
                        counts = "counts", data = "logcounts", 
                        min.cells = 3, min.features = 200)

names(sce.subset$cell.group.7) <- colnames(sce.subset)
seurat_obj@active.ident <- sce.subset$cell.group.7

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 & percent.mt < 5)

# marker genes
marker.path <- "./results/markers.qs"
if (file.exists(marker.path)) {
  markers <- qread(marker.path)
} else {
  markers <- FindAllMarkers(seurat_obj)
  qsave(markers, marker.path)
}

# dim(markers)
# head(markers)
# markers <- unique(markers$gene)
# length(markers)

# select top marker genes for each cell type
top.markers <- markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.25 & pct.2 >= 0) %>%
  top_n(wt = avg_log2FC, n = 30) %>%
  ungroup() %>% 
  select(gene) 

# top.markers <- markers %>% 
#   group_by(cluster) %>%
#   top_n(n = 50, wt = avg_log2FC) %>% 
#   filter(cluster == "Endothelial Cells" & pct.1 > 0.25 & pct.2 >= 0)

# top.markers
# head(top.markers)
# "VWF" %in%  top.markers$gene
# union of spatial and the single cell features

genes.sel <- union(rownames(corpus), intersect(rownames(counts), top.markers$gene))

length(genes.sel)
corpus <- counts[genes.sel,]
spe.subset <- spe[genes.sel, colnames(counts)]
sce.subset <- sce[genes.sel, ]
dim(spe.subset)

# stop()
# store processed Data
sp_data_qc <- file.path(dir.data, "SHK166_RA_Knee","sp_counts_SHK166_RA_Knee.mtx")
sp_data_qc_berglund <- file.path(dir.data, "SHK166_RA_Knee","sp_counts_SHK166_RA_Knee.tsv")
bc_qc <- file.path(dir.data, "SHK166_RA_Knee","bc_qc.tsv")

dir.data.processed <- file.path(dir.data, "processed")

cellcount <- file.path(dir.data.processed, "sc_counts_SHK166_RA_Knee.mtx")
sp_data <- file.path(dir.data.processed, "sp_counts_SHK166_RA_Knee.mtx")
sp_data_featsel <- file.path(dir.data.processed, "sp_counts_SHK166_RA_Knee.tsv")

celltype <- file.path(dir.data.processed, "ctypes_SHK166_RA_Knee.csv")
celltype.7 <- file.path(dir.data.processed, "ctypes7_SHK166_RA_Knee.csv")

counts.sc <- counts(sce)
ctypes <- sce$ctypes

# fix typo => "Myeloid"
ctypes.7 <- sce$cell.group.7 %>%
  as.character() %>%
  replace(., . == "Meyloid", "Myeloid") %>%
  factor()

unique(ctypes.7)
sce$cell.group.7 <- ctypes.7

qsave(spe.subset, file.path(dir.data.processed, "spe_clean_var_feats_and_markers.qs"))
qsave(spe.clean, file.path(dir.data.processed, "spe_clean.qs"))

qsave(sce.subset, file.path(dir.data.processed, "sce_subset.qs"))

writeMM(counts.sc, file = cellcount, row.names = TRUE)
writeMM(as(counts, "sparseMatrix"), file = sp_data_qc, row.names = TRUE)

write.table(counts, file = sp_data_qc_berglund, sep = "\t", row.names = TRUE)
write.table(as.matrix(corpus), file = sp_data_featsel, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(ctypes, file = celltype, sep = ",", row.names = TRUE)
write.table(ctypes.7, file = celltype.7, sep = ",", row.names = TRUE)


