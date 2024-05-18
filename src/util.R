library(SpatialExperiment)
library(igraph)


load.data  <- function(dir) {
  # load counts and coords
  path.pos <- paste(dir, 'spatial/tissue_positions.csv', sep = '/')
  pos.info <- read.csv(path.pos)
  
  ## filt out out-tissue pixels
  pos.info <- subset(pos.info, in_tissue == 1)
  pos <- pos.info[, c('pxl_row_in_fullres', 'pxl_col_in_fullres')]
  rownames(pos) <- pos.info$barcode
  
  path.gexp <- paste(dir, 'filtered_feature_bc_matrix/matrix.mtx.gz', sep = '/')
  gexp <- Matrix::readMM(path.gexp)
  
  path.bc <- paste(dir, 'filtered_feature_bc_matrix/barcodes.tsv.gz', sep = '/')
  barcodes <- read.csv(path.bc, sep="\t", header=FALSE)
  colnames(gexp) <- barcodes$V1
  
  path.feature <- paste(dir, 'filtered_feature_bc_matrix/features.tsv.gz', sep = '/')
  features <- read.csv(path.feature, sep="\t", header=FALSE)
  dim(features)
  rownames(gexp) <- features$V2
  
  # build se object
  pos <- as.matrix(pos)
  ## "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
  colnames(pos) <- c("x", "y")
  
  se <- SpatialExperiment(assay = gexp, spatialCoords = pos)
  
  assayNames(se) <- 'counts'
  
  return(se)
}


get.dirs <- function(parent.dir) {
  
  fn <- list.files(data_dir)  # 5  biopsies
  biopsy_dirs <- character(0)
  biopsy_names <- character(0)
  
  for (name in fn) {
    path <- paste(data_dir, name, sep = '/')
    biopsy_dirs <- c(biopsy_dirs, path)
    
    file_path <-  paste(path, "filtered_feature_bc_matrix.h5", sep = '/')
    if (!file.exists(file_path)) {
      path <-  paste(path, "outs", sep = '/')
    }
    
    biopsy_names <- c(biopsy_names, name)
    
  }
  names(biopsy_dirs) <-  biopsy_names
  
  # print(object.size(biopsy_objs), units = 'Mb')
  return(biopsy_dirs)
  
}


# TODO: write test 
# each point has at most 8 neighbours, return row idxs of neighbours including itself
get.adjacency.matrix <- function(coords) {
  dim.bc <- dim(coords)[1]
  adjacency.matrix <-  matrix(0, nrow = dim.bc, ncol = dim.bc)
  
  for (idx.1 in 1:dim.bc) {
    for (idx.2 in idx.1 + seq_len(dim.bc - idx.1)) {
      
      x1 <- coords[idx.1, ][1]
      y1 <- coords[idx.1, ][2]
      
      x2 <- coords[idx.2, ][1]
      y2 <- coords[idx.2, ][2]
      
      # nei.diff <- c(-1, 0, 1)
      nei.diff <- seq(-3, 3)
      is.nei <- ((x1 - x2) %in% nei.diff) && ((y1 - y2) %in% nei.diff)
      if (is.nei) {
        adjacency.matrix[idx.1, idx.2] <- 1
        adjacency.matrix[idx.2, idx.1] <- 1
      }
    }
  }
  return(adjacency.matrix)
}


get.min.dists.mat <- function(g, from, to) {
  # get components
  membership <- components(g)$membership
  # find shortest dist, DP 
  min.dists.mat <- matrix(data = NA, nrow = length(from), ncol = length(to))
  
  for (i in from) {
    for (j in seq_len(length(to))) {
      target.idx <- to[j]
      target.idx
      if (!is.na(min.dists.mat[i, j])) {
        next
      }else if(membership[i] != membership[target.idx]){
        # no path between them
        cnt.v <- Inf
      }else if(i == target.idx) {
        cnt.v <- 0
      }else{
        shortest_path <- shortest_paths(g, from = i, to = target.idx, output = "vpath")
        vpath <- shortest_path$vpath 
        vpath <- vpath[[1]]
        len.path <- length(vpath)
        if (len.path == 0) {
          cnt.v <- Inf
        }else{
          # if a vertex is on the shortest path, we can calculate the shortest distance to the target vertex
          for (k in seq_len(len.path)) {
            idx <- vpath[k]
            min.dists.mat[idx, target.idx] <- len - k
          }
        }
      }
    }
    # break
  }
  return(min.dists.mat)
}


