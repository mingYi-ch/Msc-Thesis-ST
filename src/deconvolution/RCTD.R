library(qs)
library(spacexr)
library(SpatialExperiment)

source("./src/util.R")
source("./config/config.R")

# load data as spatialExperiment object
data.dir <- get.data.dir()
spe.path <- file.path(data.dir, "processed", "spe_subset.qs")
spe <- qread(spe.path)

# load reference data set
dir.res <- get.res.dir()
ref <- qread(file.path(dir.res, "SCE/sce.qs"))
counts.sc <- counts(ref)

ctypes <- ref$cell.group.7
nUMI <- colSums(counts.sc)
names(nUMI) <- colnames(counts.sc) # create nUMI named list
names(ctypes) <- colnames(counts.sc) # create nUMI named list
ref <- Reference(counts.sc, ctypes, nUMI)
counts <- assay(spe, "counts") # gene by barcode

# Count occurrences of each column name
colname_counts <- table(rownames(counts))

# Identify and count duplicates: TODO CHECK DUPLICATES
duplicate_counts <- colname_counts[colname_counts > 1]
counts <- counts[!duplicated(rownames(counts)), ] # duplicates "TBCE"   "HSPA14"

coords <- spatialCoords(spe)
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

coords_df <- as.data.frame(coords)
puck <- SpatialRNA(coords_df, counts, nUMI)

## Examine SpatialRNA object (optional)
# print(dim(puck@counts)) # observe Digital Gene Expression matrix
# hist(log(puck@nUMI,2)) # histogram of log_2 nUMI

myRCTD <- create.RCTD(puck, ref, max_cores = 16, CELL_MIN_INSTANCE = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# normalization per pixel
results <- myRCTD@results
norm_weights <- normalize_weights(results$weights)
myRCTD@results$norm_weights <- norm_weights

# save results, TODO: save dataset name in file names
res.dir <- get.res.dir()
res.dir <- file.path(res.dir, "deconvolution", "RCTD")
if (!dir.exists(res.dir)) {
  dir.create(res.dir)
}
qs::qsave(myRCTD, file.path(res.dir, "myRCTD.qs"))
