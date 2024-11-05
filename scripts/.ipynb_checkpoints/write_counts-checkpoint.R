library(qs)
library(spacexr)
library(SpatialExperiment)

source("./src/util.R")
source("./config/config.R")

# path
cellcount  <-  "~/data/scRNA-seq/sc_counts_SHK166_RA_Knee.mtx"
sp_data <- "~/data/spatial-trancriptomics/st_counts_SHK166_RA_Knee.mtx"
sp_data_berglund <- "~/data/spatial-trancriptomics/st_counts_SHK166_RA_Knee.tsv"

celltype <- "~/data/scRNA-seq/ctypes_counts_SHK166_RA_Knee.csv"
celltype.7 <- "~/data/scRNA-seq/ctypes7_counts_SHK166_RA_Knee.csv"

# load data as spatialExperiment object
data_dir <- get.std.path()
spe <- load.data(data_dir)

# load reference scRNA-seq
# ref_path <- get.refd.path()
# ref <- readRDS(ref_path)
res.dir <- get.res.dir()
sce <- qread(file.path(res.dir, "SCE/sce.qs"))

counts.sc <- counts(sce)
ctypes <- sce$ctypes
ctypes.7 <- sce$cell.group.7

writeMM(counts.sc, file = cellcount, row.names = TRUE)
writeMM(counts(spe), file = sp_data, row.names = TRUE)

write.table(counts(spe), file = sp_data_berglund, sep = "/t", row.names = TRUE, col.names = TRUE)

write.table(ctypes, file = celltype, sep = ",", row.names = TRUE)
write.table(ctypes.7, file = celltype.7, sep = ",", row.names = TRUE)

