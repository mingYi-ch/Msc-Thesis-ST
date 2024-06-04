library(qs)
library(spacexr)
library(SpatialExperiment)

source("./src/util.R")
source("./config/config.R")

# path

cellcount  <-  "~/data/scRNA-seq/sc_counts_SHK166_RA_Knee.mtx"
sp_data <- "~/data/spatial-trancriptomics/st_counts_SHK166_RA_Knee.mtx"
celltype <- "~/data/scRNA-seq/ctypes_counts_SHK166_RA_Knee.csv"

# load data as spatialExperiment object
data_dir <- get.std.path()
se <- load.data(data_dir)

# load reference scRNA-seq
ref_path <- get.refd.path()
ref <- readRDS(ref_path)
counts.sc <- ref@counts
ctypes <- ref@cell_types

writeMM(counts.sc, file = cellcount, row.names = TRUE)
writeMM(counts(se), file = sp_data, row.names = TRUE)
write.table(ctypes, file = celltype, sep = ",", row.names = TRUE)

