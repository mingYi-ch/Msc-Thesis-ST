# set project dir for local and S3IT
host <- "local"
if (dir.exists("/home/myi")) {
  host <- "cluster"
}
proj.dir.local <- "/Users/Yi/Library/Mobile\ Documents/com~apple~CloudDocs/Master-thesis/Msc-Thesis-ST"
proj.dir.cluster <- "/data/myi/mt-spatial-transcriptomics"

# define spatial transcriptomic data paths
get.std.path <- function() {
  path <- '/shares/rheumatologie.usz/caroline/spatial/SHK166_RA_Knee'
  if (host == "local") {
    path <- paste(proj.dir.local, 'data/SHK166_RA_Knee', sep = '/')
  }
  return(path)
}


# reference scRNA-seq data path
get.refd.path <- function() {
  path <- "~/data/mt-spatial-transcriptomics/data/scRNA-seq/ref_from_merged_subpop_down_main_clusters_for_CCI.rds"
  if (host == "local") {
    path <- paste(proj.dir.local, 'data/reference/ref_from_merged_subpop_down_main_clusters_for_CCI.rds', sep = '/')
  }
  return(path)
}


# results dir
get.res.dir <- function() {
  path <- file.path(proj.dir.cluster, "results")
  if (host == "local") {
    path <- file.path(proj.dir.local, "results")
  }
  
  # Check if the path exists, and if not, create it
  path.plot <- file.path(path, "plots")
  if (!dir.exists(path.plot)) {
    dir.create(path.plot, recursive = TRUE)  # 'recursive = TRUE' creates any necessary parent directories
    message("Directory created at: ", path.plot)
  }
  return(path)
}

# data dir
get.data.dir <- function() {
  path <- file.path(proj.dir.cluster, "data")
  if (host == "local") {
    path <- file.path(proj.dir.local, "data")
  }
  return(path)
}