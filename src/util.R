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
get.adjacency.matrix <- function(coords, nei.diff = seq(0, 3), dist.sqr.threshold = 9000 ) {
  dim.bc <- dim(coords)[1]
  adjacency.matrix <-  matrix(0, nrow = dim.bc, ncol = dim.bc)
  
  for (idx.1 in 1:dim.bc) {
    for (idx.2 in idx.1 + seq_len(dim.bc - idx.1)) {
      
      pt1 <- coords[idx.1, ]
      pt2 <- coords[idx.2, ]
      
      is.nei <- sum((pt1 - pt2)^2) < dist.sqr.threshold      
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
  len.from <- length(from)
  len.to <- length(to)
  
  from.map <- as.list(seq_len(len.from))
  names(from.map) <- from
  
  min.dists.mat <- matrix(data = NA, nrow = len.from, ncol = len.to)
  
  for (i in seq_len(len.from)) {
    for (j in seq_len(len.to)) {
      from.idx <- from[i]
      to.idx <- to[j]
      if (!is.na(min.dists.mat[i, j])) {
        next
      }else if(membership[from.idx] != membership[to.idx]){
        # no path between them
        cnt.v <- Inf
      }else if(from.idx == to.idx) {
        cnt.v <- 0
      }else{
        shortest_path <- shortest_paths(g, from = from.idx, to = to.idx, output = "vpath")
        vpath <- shortest_path$vpath 
        vpath <- vpath[[1]]
        len.path <- length(vpath)
        # if a vertex is on the shortest path, we can calculate the shortest distance to the target vertex
        for (k in seq_len(len.path)) {
          idx <- vpath[k]
          cat(idx, to.idx)
          min.dists.mat[from.map$idx, as.character(j)] <- len.path - k
      }
      }
    }
    # break
  }
  return(min.dists.mat)
}


