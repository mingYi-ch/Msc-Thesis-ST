library(qs)
library(STdeconvolve)
library(scuttle)
library(Seurat)

source("./src/util.R")
source("./config/config.R")

# load data as spatialExperiment object
data_dir <- get.std.path()
spe <- load.data(data_dir)

cd <- assay(spe, "counts")

# QC: Filter out poor pixels and genes
counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10, verbose = TRUE)

# Select spatial features
corpus <- restrictCorpus(counts, 
                         removeAbove=0.95, 
                         removeBelow = 0.05, 
                         nTopOD = NA,
                         plot = FALSE,
                         alpha = 0.5) # 1e-8
# remove mitochondria genes
mask.mt <- startsWith(rownames(corpus), "MT-") | startsWith(rownames(corpus), "mt-")
mt_genes <- corpus[mask.mt]
corpus <- corpus[!mask.mt]


# load reference data set
dir.res <- get.res.dir()
sce <- qread(file.path(dir.res, "SCE/sce.qs"))

shared.genes <- intersect(rownames(sce), rownames(spe))
sce.subset <- sce[shared.genes,]

seurat_obj <- as.Seurat(sce.subset, counts = "counts", data = "logcounts")
names(sce.subset$cell.group.7) <- colnames(sce.subset)
seurat_obj@active.ident <- sce.subset$cell.group.7

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

# select top 200 genes for each cell type
top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 200, with_ties = FALSE) %>%
  ungroup() %>% 
  select(gene) %>% 
  distinct()

dim(top.markers)

# union of spatial and the single cell features
genes.sel <- union(rownames(corpus), intersect(rownames(corpus), top.markers$gene))

corpus <- counts[genes.sel,]
spe.subset <- spe[genes.sel, colnames(spe)]
sce.subset <- sce[genes.sel, ]


# stop()
# store processed Data
dir.data <- get.data.dir()
dir.data <- file.path(dir.data, "processed")

cellcount <- file.path(dir.data, "sc_counts_SHK166_RA_Knee.mtx")
sp_data <- file.path(dir.data, "sp_counts_SHK166_RA_Knee.mtx")
sp_data_berglund <- file.path(dir.data, "sp_counts_SHK166_RA_Knee.tsv")

celltype <- file.path(dir.data, "ctypes_SHK166_RA_Knee.csv")
celltype7 <- file.path(dir.data, "ctypes7_SHK166_RA_Knee.csv")


counts.sc <- counts(sce)
ctypes <- sce$ctypes
ctypes.7 <- sce$cell.group.7

qsave(spe.subset, file.path(dir.data, "spe_subset.qs"))
qsave(sce.subset, file.path(dir.data, "sce_subset.qs"))

writeMM(counts.sc, file = cellcount, row.names = TRUE)
writeMM(as(corpus, "sparseMatrix"), file = sp_data, row.names = TRUE)

write.table(as.matrix(corpus), file = sp_data_berglund, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(ctypes, file = celltype, sep = ",", row.names = TRUE)
write.table(ctypes.7, file = celltype.7, sep = ",", row.names = TRUE)


