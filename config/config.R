# check host
host <- "local"
if (dir.exists("/home/myi")) {
  host <- "cluster"
}
proj.dir.local <- "/Users/Yi/Library/Mobile\ Documents/com~apple~CloudDocs/Master-thesis/Msc-Thesis-ST"
proj.dir.cluster <- "/home/myi/mt-spatial-transcriptomics/"

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
  path <- "~/data/ref_from_merged_subpop_down_main_clusters_for_CCI.rds"
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
  return(path)
}
