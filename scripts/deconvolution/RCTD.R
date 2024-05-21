library(qs)
source("./src/util.R")

# load data as spatialExperiment object
# TODO config
data_dir <- '/shares/rheumatologie.usz/caroline/spatial/SHK166_RA_Knee'
data_dir <- './data/SHK166_RA_Knee'
se <- load.data(data_dir)

cd <- assay(se, "counts")

# load reference data set
dir.cluster <- "~/data/"
dir.local <- "../../Spatial/"
ref.name <- "ref_from_merged_subpop_down_main_clusters_for_CCI.rds"

ref.dir <- dir.cluster
if (!dir.exists(ref.dir)) {
  ref.dir <- dir.local
}

ref <- readRDS(paste0(ref.dir, ref.name))
