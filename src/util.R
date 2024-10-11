library(SpatialExperiment)
library(igraph)
library(ggplot2)
library(jsonlite)
library(scales)
# Function to calculate the average distance for each cluster
calculate.avg.dist <- function(x) {
  # Remove square brackets and parentheses
  cleaned_str <- str_replace_all(x, "[\\[\\]()]", "")
  
  # Convert the cleaned string to numeric values
  numeric_values <- strsplit(cleaned_str, ",")[[1]] %>% as.numeric()
  
  avg <- mean(numeric_values)
  if (!is.finite(avg)) {
    avg <- numeric_values[1]
  }
  
  return(avg)
}

stat_func_mean <- function(data, indices) {
    return(mean(data[indices], na.rm = TRUE))
}

stat_func_median <- function(data, indices) {
    return(median(data[indices], na.rm = TRUE))
}

# Function to perform bootstrapping and add orig and expr columns
# perform_bootstrapping_func <- function(origin, gene, sce.pseudo, stat_func) {
#   tibble(counts = as.data.frame(assay(sce.pseudo, origin)[1,])) %>% 
#     mutate(dist.cluster = sce.pseudo$dist.cluster) %>% 
#     # filter(complete.cases(counts)) %>% 
#     group_by(dist.cluster) %>% 
#     do({
#       data <- .$counts[!is.na(.$counts)]
#       
#       # in case no data to boot
#       if (length(data) == 0) {
#         expr <-  NA; lower <- NA; upper <- NA
#       }else{
#         boot.res <- boot(data, statistic = stat_func, R = 1000)
#         ci <- boot.ci(boot.res, conf = c(0.95), type = c('perc'))
#         expr <-  boot.res$t; lower <- ci$percent[4]; upper <- ci$percent[5]
#       }
#       
#       data.frame(expr = expr, n = nrow(.), lower = lower, upper = upper)
#     }) %>% 
#     ungroup() %>%
#     mutate(origin = origin, symbol = gene)
# }

perform_bootstrapping_func <- function(origin, genes_of_interest, sce.pseudo, stat_func) {
  # # mask for if a spot is in neighborhood
  # mask <- matrix(NA, nrow = length(genes_of_interest), ncol = ncol(fit_all), 
  #                dimnames = list(genes_of_interest, colnames(fit_all)))
  # 
  # mask2 <- matrix(1, nrow = length(genes_of_interest), ncol = ncol(fit_all), 
  #                 dimnames = list(genes_of_interest, colnames(fit_all)))
  # 
  # for(id in genes_of_interest){
  #   mask[id, filter(nei, name == id)$neighborhood[[1]]] <- 1
  #   mask2[id, filter(nei, name == id)$neighborhood[[1]]] <- NA
  # }
  # 
  # 
  # # construct two assays: inside and outside
  # sce.pseudo <- SingleCellExperiment(
  #   list(inside = as.matrix(logcounts(fit_all)[genes_of_interest,,drop=FALSE] * mask),
  #        outside = as.matrix(logcounts(fit_all)[genes_of_interest,,drop=FALSE] * mask2)
  #   ),
  #   colData = as.data.frame(colData(fit_all))
  # )
  
  res <- tibble(expr = as.data.frame(assay(sce.pseudo, origin)[genes_of_interest,])) %>% 
    mutate(dist.cluster = sce.pseudo$dist.cluster) %>% 
    group_by(dist.cluster) %>% 
    do({
      cluster <- unique(.$dist.cluster) %>% as.character()
      R <-  500
      # remove NA before sampling
      data <- .$expr[!is.na(.$expr)]
      
      # no data in a cluster
      if (length(data) == 0) {
        data.frame(
          expr = rep(0, R), 
          n = nrow(.), 
          lower = NA, 
          upper = NA, 
          SampleName = rep(glue("{origin}_{genes_of_interest}_{cluster}"), R),
          dist.cluster.1 = cluster,
          Replicate = seq(R)
        )
      }else{
        boot.res <- boot(data, statistic = mean_stat, R = R)
        ci <- boot.ci(boot.res, conf = c(0.95), type = c('perc'))
        lower <-  ci$percent[4]
        upper <-  ci$percent[5]
        
        if (length(lower) == 0) {
          lower <-  NA
          upper <-  NA
        }
        data.frame(
          expr = boot.res$t, 
          n = nrow(.), 
          lower = lower, 
          upper = upper, 
          SampleName = rep(glue("{origin}_{genes_of_interest}_{cluster}"), R),
          dist.cluster.1 = cluster,
          Replicate = seq(R)
        )
      }
    }) %>% 
    ungroup() %>%
    mutate(origin = origin, symbol = genes_of_interest)
  return(res)
}

lemur_nei_plot <- function(spe, nei, fit, gene_of_interest) {
  # refactor to a func
  is_inside <- tibble(symbol = gene_of_interest, cell_id = list(colnames(fit))) %>%
    left_join(dplyr::select(nei, name, neighborhood), by = c("symbol"= "name")) %>%
    mutate(inside = map2(cell_id, neighborhood, \(ref, nei_pixels) ref %in% nei_pixels)) %>%
    dplyr::select(-neighborhood) %>%
    unnest(c(inside, cell_id))
  
  # indicate if a spot is inside the neighborhood
  inside <- rep("excluded", ncol(spe))
  names(inside) <- colnames(spe)
  nei_id <- is_inside %>% 
    filter(inside == TRUE) %>% 
    pull(cell_id)
  
  nonnei_id <- is_inside %>% 
    filter(inside == FALSE) %>% 
    pull(cell_id)
  
  inside[nei_id] <- "nei"
  inside[nonnei_id] <- "non-nei"
  inside.factor <- factor(inside, levels = c("nei", "non-nei", "excluded"))
  specified_colors <- c("nei" = "blue", "non-nei" = "lightblue", "excluded" = "lightgrey")
  
  info_of_interest <- "neibourhood spots"
  cluster <- target.cell
  
  groups <- rep("others", ncol(spe))
  groups[unlist(is.ec.spot)] <- cluster
  
  strokes <- ifelse(groups == cluster, 0.5, 0)
  
  size.axis.text <- 8 * 2
  size.axis.title <- 8.5 * 2
  size.axis <- 7 * 2
  
  plots_ECs_nei <- as_tibble(spatialCoords(spe)) %>%
    mutate(groups = groups, strokes = strokes, inside.col = inside.factor) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(stroke = strokes, shape = groups, fill = inside.col), size = 2, alpha = 0.5, color = "black") + # here is the border color
    scale_shape_manual(values = c(24, 21)) +
    scale_fill_manual(values = specified_colors) +
    scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +  # Format y-axis labels
    scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +  # Format y-axis labels
    labs(subtitle = glue("{gene_of_interest} neighborhood cluster"),
         x = "row pixel",
         y = "col pixel",
         fill = "spots category by color: ",
         shape = "cell type group by shape:") + 
    
    theme(plot.subtitle = element_text(size = size.axis.title),
          legend.text = element_text(size = size.axis.text),
          legend.title = element_text(size = size.axis.text),
          axis.title.x = element_text(size = size.axis.text),
          axis.title.y = element_text(size = size.axis.text),
          axis.text.x = element_text(size = size.axis),  
          axis.text.y = element_text(size = size.axis),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2, override.aes = list(shape = 22, color = c("blue", "lightblue", "lightgrey"), size = 5)),
           shape = guide_legend(nrow = 2, override.aes = list(size = 5))  # Set legend as lines
    )
  return(plots_ECs_nei)
}

# retactor to multiple funcs
lemur_curve_plot <- function(nei, fit, gene_of_interest, spe, stat_func) {
  set.seed(1)
  mask <- matrix(NA, nrow = length(gene_of_interest), ncol = ncol(fit), 
                 dimnames = list(gene_of_interest, colnames(fit)))
  
  mask2 <- matrix(1, nrow = length(gene_of_interest), ncol = ncol(fit), 
                  dimnames = list(gene_of_interest, colnames(fit)))
  
  for(id in gene_of_interest){
    mask[id, filter(nei, name == id)$neighborhood[[1]]] <- 1
    mask2[id, filter(nei, name == id)$neighborhood[[1]]] <- NA
  }
  
  sce.pseudo <- SingleCellExperiment(
    list(inside = as.matrix(logcounts(fit)[gene_of_interest,,drop=FALSE] * mask),
         outside = as.matrix(logcounts(fit)[gene_of_interest,,drop=FALSE] * mask2)
    ),
    colData = as.data.frame(colData(fit))
  )
  
  # Perform bootstrapping for both "inside" and "outside"
  comparison_data <- bind_rows(
    perform_bootstrapping_func("inside", gene_of_interest, sce.pseudo, stat_func_mean),
    perform_bootstrapping_func("outside", gene_of_interest, sce.pseudo, stat_func_mean)
  )
  
  comparison_data <- tibble(comparison_data) %>%
    mutate(dist.cluster = ordered(dist.cluster, levels = levels(fit$colData$dist.cluster))) %>%
    arrange(dist.cluster)
  
  comparison_data$avg.dist <- lapply(as.character(comparison_data$dist.cluster), calculate.avg.dist) %>% unlist()
  
  # plot curves
  ## set y limits
  y_min <- min(comparison_data$expr)
  y_max <- max(comparison_data$expr)
  y_expansion <- 0.6
  
  size.axis.text <- 7.5 * 2
  size.axis.title <- 8 * 2
  size.axis <- 7 * 2
  
  ## degree of freedom of natural cubic spline 
  df <- 3
  
  expr_comparison_pl <- as_tibble(comparison_data) %>%
    mutate(avg.dist = as.factor(comparison_data$avg.dist), origin = as.factor(comparison_data$origin)) %>% 
    ggplot(aes(x = avg.dist, y = expr)) +
    geom_point(aes(color = origin), size = 1, alpha = 0) +
    geom_boxplot(aes(color = origin), width = 0.3, outlier.shape = NA) +  # Add boxplot
    geom_smooth(aes(color = origin, x = as.numeric(avg.dist), y = expr), method = 'lm', formula = y ~ ns(x, df = df), se = FALSE, linewidth = 0.5) +
    # geom_errorbar(aes(ymin = lower, ymax = upper, color = origin), width = 0.8, alpha = 0.5) +
    scale_color_manual(values = c("inside" = "black", "outside" = "lightgrey"), labels = c("inside" = "spots in neighborhood", "outside" = "All other spots")) +
    scale_x_discrete(breaks = unique(comparison_data$avg.dist)) +  # Show all x values
    scale_y_continuous(limits = c(y_min, y_max), expand = expansion(add = y_expansion)) +
    guides(x = guide_axis(angle = 45)) +
    
    labs(color = "",
         x = "avg.dist / spots",
         y = "Log Expr.",
         subtitle = paste0(gene_of_interest, " expr. vs. Dist to endothelial cell"))  +
    
    theme(plot.subtitle = element_text(size = size.axis.title),
          legend.text = element_text(size = size.axis.text),
          axis.title.x = element_text(size = size.axis.text),
          axis.title.y = element_text(size = size.axis.text),
          axis.text.x = element_text(size = size.axis),  
          axis.text.y = element_text(size = size.axis),
          legend.position = "bottom") +
    # annotate("text", x = Inf, y = Inf, label = glue("P.Val = {P.Val}"), hjust = 1.1, vjust = 1.5, size = 3) + 
    guides(color = guide_legend(override.aes = list(linetype = 1, size = 5)))  # Set legend as lines
  
  # Save the plot with rectangular dimensions
  return(expr_comparison_pl)
}

# add Inf, n + 1 groups
bin.dists <- function(x, n = 5, from = min(x), to = max(x[is.finite(x)]), add_label = TRUE, label_fmt = "%.2f"){
  groups <- seq(from, to, length.out = n+1)
  
  res <- rep(NA_integer_, length(x))
  labels <- rep(NA_character_, n + 1)
  
  groups <- c(groups, Inf)
  res[is.infinite(x)] <- n + 1
  labels[n + 1] <- "Inf"
  
  for(idx in seq_len(n)){
    if(idx == n){
      labels[idx] <- paste0("[", sprintf(label_fmt, groups[idx]), ",", sprintf(label_fmt, groups[idx+1]), "]")
      res[x >= groups[idx] & x <= groups[idx+1]] <- idx
    }else{
      labels[idx] <- paste0("[", sprintf(label_fmt, groups[idx]), ",", sprintf(label_fmt, groups[idx+1]), ")")
      res[x >= groups[idx] & x < groups[idx+1]] <- idx
    }
  }
  if(add_label){
    factor(res, levels = seq_len(n + 1), labels = labels)
  }else{
    res
  }
}

# add Inf, n + 1 groups
# !! Values in groups need to be changed due to different value range
bin.dists2 <- function(x, groups, add_label = TRUE, label_fmt = "%.3f"){
  # groups <- seq(from, to, length.out = n+1)
  n <- length(groups) -1
  
  res <- rep(NA_integer_, length(x))
  labels <- rep(NA_character_, n + 1)
  
  groups <- c(groups, Inf)
  res[is.infinite(x)] <- n + 1
  labels[n + 1] <- "Inf"
  
  for(idx in seq_len(n)){
    if(idx == n){
      labels[idx] <- paste0("[", sprintf(label_fmt, groups[idx]), ",", sprintf(label_fmt, groups[idx+1]), "]")
      res[x >= groups[idx] & x <= groups[idx+1]] <- idx
    }else{
      labels[idx] <- paste0("[", sprintf(label_fmt, groups[idx]), ",", sprintf(label_fmt, groups[idx+1]), ")")
      res[x >= groups[idx] & x < groups[idx+1]] <- idx
    }
  }
  if(add_label){
    factor(res, levels = seq_len(n + 1), labels = labels)
  }else{
    res
  }
}

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
  rownames(gexp) <- make.unique(features$V2) # TODO: check with Prof. if right
  
  # build spe object
  pos <- as.matrix(pos)
  ## "pxl_col_in_fullres" = "y" coordinates, and "pxl_row_in_fullres" = "x" coordinates
  colnames(pos) <- c("x", "y")
  
  # make bc the same order
  gexp <- gexp[,rownames(pos)]
  spe <- SpatialExperiment(assay = gexp, spatialCoords = pos)
  
  assayNames(spe) <- 'counts'
  
  # load image
  scalefactors <- fromJSON(file.path(dir, "spatial/scalefactors_json.json"))
  img.path <- file.path(dir, "spatial/tissue_hires_image.png")
  spe <- addImg(spe, 
                sample_id = "sample01", 
                image_id = "tissue_hires",
                imageSource = img.path, 
                scaleFactor = scalefactors$tissue_hires_scalef, 
                load = TRUE)
  return(spe)
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
# each point has at most 6 neighbours, return row idxs of neighbours including itself
get.adjacency.matrix <- function(coords, nei.diff = seq(0, 3), squared.radius = 9000 ) {
  dim.bc <- dim(coords)[1]
  adjacency.matrix <-  matrix(0, nrow = dim.bc, ncol = dim.bc)
  
  for (idx.1 in 1:dim.bc) {
    for (idx.2 in idx.1 + seq_len(dim.bc - idx.1)) {
      
      pt1 <- coords[idx.1, ]
      pt2 <- coords[idx.2, ]
      
      is.nei <- sum((pt1 - pt2)^2) < squared.radius      
      if (is.nei) {
        adjacency.matrix[idx.1, idx.2] <- 1
        adjacency.matrix[idx.2, idx.1] <- 1
      }
    }
  }
  return(adjacency.matrix)
}


get.all.dists <- function(g, from, to) {
  # get components
  membership <- igraph::components(g)$membership
  
  len.from <- length(from)
  len.to <- length(to)
  
  # idx in from may be not continuous number
  from.map <- as.list(seq_len(len.from))
  names(from.map) <- from
  
  dists.mat <- matrix(
    data = NA, 
    nrow = len.from, 
    ncol = len.to,
    dimnames = list(from, to)
    )
  
  for (i in seq_len(len.from)) {
    for (j in seq_len(len.to)) {
      from.idx <- from[i]
      to.idx <- to[j]
      
      if (!is.na(dists.mat[i, j])) {
        # print(to.idx)
        next
      }else if(membership[from.idx] != membership[to.idx]){
        # no path between them
        dists.mat[i, j] <- Inf
      }else if(from.idx == to.idx) {
        dists.mat[i, j] <- 0
      }else{
        shortest_path <- shortest_paths(g, from = from.idx, to = to.idx, output = "vpath")
        vpath <- shortest_path$vpath[[1]]
        
        # if a vertex is on the shortest path, we can calculate the shortest distance to the target vertex
        for (k in seq_along(vpath)) {
          idx <- vpath[k]
          idx <-  as.character(idx)
          dists.mat[from.map[[idx]], j] <- length(vpath) - k
        }
      }
    }
  }
  return(dists.mat)
}

# get.min.dists <- function(dists.mat, weights) {
#   # scale to the range of 0 to 1
#   dists.mat.finite <- dists.mat[is.finite(dists.mat)]
#   
#   spl.norm <- (dists.mat - min(dists.mat)) / (max(dists.mat.finite) - min(dists.mat))
#   spl.norm.prop <- sweep(spl.norm, MARGIN = 2, STATS = weights, FUN="/") # divide each row by target vector
#   
#   # get min dist to a certain target
#   min.dists <- rowMins(spl.norm.prop)
#   
#   return(min.dists)
# }

get.min.dists <- function(dists.mat, weights, normalize = FALSE) {
  # Scale by weights
  dists.mat <- sweep(dists.mat, MARGIN = 2, STATS = weights, FUN = "/")
  
  # Get minimum distances
  min.dists <- rowMins(dists.mat)
  
  # Ensure non-Inf values are considered
  min.mat.finite <- min.dists[is.finite(min.dists)]
  
  if (normalize) {
    # Normalize distances
    min.val <- min(min.mat.finite)
    max.val <- max(min.mat.finite)
    
    spl.norm <- (min.dists - min.val) / (max.val - min.val)
    # spl.norm[!is.finite(spl.norm)] <- 1  # Assign 1 to any Inf values after normalization
    
    # spl.norm[is.na(spl.norm)] <- Inf
    return(spl.norm)
  }
  return(min.dists)
}

######### Custom ggplot2 theme, code from https://github.com/const-ae/lemur-Paper/blob/master/notebooks/util.R #########

font_size <- 8
font_size_small <- 6
font_size_tiny <- 5
font_size_large <- 10
publ_theme <- cowplot::theme_cowplot(font_size = font_size, rel_small = font_size_small / font_size,
                                     rel_tiny = font_size_tiny / font_size, rel_large = font_size_large / font_size,
                                     line_size = 0.3) +
  theme(plot.title = element_text(size = font_size),
        axis.title = element_text(size = font_size_small),
        legend.title = element_text(size = font_size_small),
        strip.background = element_blank(),
        strip.text = element_text(size = font_size_small))
theme_set(publ_theme)

small_axis <- function(label = NULL, fontsize = 7, arrow_length = 10, label_offset = 1, fix_coord = TRUE, remove_axes = TRUE,
                       arrow_spec = grid::arrow(ends = "both", type = "closed", angle = 20, length = unit(arrow_length / 7, units)),
                       units = "mm", ...){
  coord <- if(fix_coord){
    coord_fixed(clip = "off", ...)
  }else{
    NULL
  }
  axis_theme <- if(remove_axes){
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  }else{
    NULL
  }
  lines <- annotation_custom(grid::polylineGrob(x = unit(c(0, 0, arrow_length), units), y = unit(c(arrow_length, 0, 0), units), 
                                                gp = grid::gpar(fill = "black"),
                                                arrow = arrow_spec))
  text <- if(! is.null(label)){
    annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                     x = unit(label_offset, units), y = unit(label_offset, units), hjust = 0, vjust = 0))
  }else{
    NULL
  }
  list(coord, axis_theme, lines, text)
}

small_arrow <- function(position = c(0.8, 0.95), offset = 0.01, label = NULL, direction = c("x", "y"), 
                        fontsize = 7, arrow_length = unit(10 / 7, "mm"), label_offset = 0, label_hjust = NULL, label_vjust = NULL,
                        arrow_spec = grid::arrow(ends = "last", type = "closed", angle = 20, length = arrow_length),
                        units = "npc"){
  direction <- match.arg(direction)
  if(!grid::is.unit(position)){
    position <- grid::unit(position, units = units)
  }
  if(!grid::is.unit(offset)){
    offset <- grid::unit(offset, units = units)
  }
  if(!grid::is.unit(label_offset)){
    label_offset <- grid::unit(label_offset, units = units)
  }
  if(direction == "x"){
    arrow <- annotation_custom(grid::polylineGrob(x = position, y = c(offset, offset), 
                                                  gp = grid::gpar(fill = "black"),
                                                  arrow = arrow_spec))
    text <- if(! is.null(label)){
      annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                       x = (position[1] + position[2]) / 2, y = offset + label_offset,
                                       hjust = label_hjust, vjust = label_vjust))
    }
  }else{
    arrow <- annotation_custom(grid::polylineGrob(y = position, x = c(offset, offset), 
                                                  gp = grid::gpar(fill = "black"),
                                                  arrow = arrow_spec))
    text <- if(! is.null(label)){
      annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                       y = (position[1] + position[2]) / 2, x = offset + label_offset,
                                       hjust = label_hjust, vjust = label_vjust, rot = 90))
    }
  }
  list(arrow, text)
}

# scale_colour_gradient2_rev <- function(..., low = scales::muted("blue"), mid = "white", high = scales::muted("red")){
#   scale_color_gradient2(..., low = low, mid = mid, high = high)
# }
signif_to_zero <- function(x, digits = 6){
  n_signif_digits <- digits - ceiling(log10(abs(x)))
  sign(x) * floor(abs(x) * 10^n_signif_digits) / 10^n_signif_digits
}
scale_color_de_gradient <- function(abs_max, mid_width = 0.1, ..., oob = scales::squish, limits = c(-1, 1) * abs_max, breaks = c(-1, 0, 1) * signif_to_zero(abs_max, 1)){
  colors <- c(scales::muted("blue"), "lightgrey", "lightgrey", scales::muted("red"))
  values <- c(0, 0.5 - mid_width/2, 0.5 + mid_width/2, 1)
  scale_color_gradientn(oob = oob, limits = limits, breaks = breaks, colors = colors, values = values, ...)
}

######### Custom plotting functions #########

convert_dims <- function( width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1){
  units <- match.arg(units)
  if(units == "inches"){
    units <- "in"
  }
  to_inches <- function(x) x/c(`in` = 1, cm = 2.54, mm = 2.54 * 
                                 10, px = dpi)[units]
  to_inches(c(width, height)) * scale
}

my_pdf <- function(filename, width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1, ...){
  dim <- convert_dims(width, height, units, dpi, scale)
  grDevices::pdf(filename, width = dim[1], height = dim[2], useDingbats = FALSE, ...)
}


my_tikz <- function(filename, width, height, units = c("inches", "in", "cm", "mm", "px"), dpi = 300, scale = 1, stand_alone = TRUE, ...){
  dim <- convert_dims(width, height, units, dpi, scale)
  tikzDevice::tikz(filename, width = dim[1], height = dim[2], standAlone = stand_alone, 
                   documentDeclaration = c(getOption("tikzDocumentDeclaration"), r"(\renewcommand{\familydefault}{\sfdefault})", r"(\usepackage{helvet})"),  # Use sans serif font Helvetica
                   packages = c(options("tikzLatexPackages")$tikzLatexPackages, "\\usepackage{amssymb}", "\\usepackage{amsmath}", "\\usepackage{bm}"),  ...)
}

save_plot <- function(filename, plot = ggplot2::last_plot(), width = 6.2328, height = 3.71, units = c("inches", "cm", "mm", "px"), dpi = 300, scale = 1, latex_support = FALSE, ...){
  
  old_dev <- grDevices::dev.cur()
  if(latex_support){
    filename <- if(stringr::str_ends(filename, "\\.pdf")){
      paste0(stringr::str_sub(filename, end  = -5L), ".tex")
    }
    my_tikz(filename, width = width, height = height, units = units, dpi = dpi, scale = scale, stand_alone = TRUE)
  }else{
    dim <- convert_dims(width, height, units, dpi, scale)
    dev <- ggplot2:::plot_dev(NULL, filename, dpi = dpi)
    dev(filename = filename, width = dim[1], height = dim[2], ...)
    on.exit(utils::capture.output({
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
    }))
  }
  
  grid::grid.draw(plot)
  
  if(latex_support){
    grDevices::dev.off()
    if (old_dev > 1) grDevices::dev.set(old_dev)
    
    withr::with_dir(dirname(filename), {
      tools::texi2pdf(basename(filename), clean = TRUE)
      # Remove .tex file
      file.remove(basename(filename)) 
      raw_file_name <- tools::file_path_sans_ext(basename(filename))
      # rastered images
      ras_files <- list.files(pattern = paste0("^", raw_file_name, "_ras\\d*.png"))
      if(length(ras_files) > 0){
        file.remove(ras_files)
      }
    })
  }
  
  invisible(filename)  
}

plot_assemble <- function(..., .plot_objs = NULL, width = 6.2328, height = 3.71, units = c("inches", "cm", "mm", "px"), 
                          latex_support = FALSE, show_grid_lines = TRUE, filename = NULL){
  units <- match.arg(units)
  
  plots <- if(is.null(.plot_objs)){
    list(...)
  }else{
    .plot_objs
  }
  
  if(show_grid_lines){
    x_breaks <- scales::breaks_pretty(n = 10)(seq(0, width, length.out = 100))
    y_breaks <- scales::breaks_pretty(n = 10)(seq(0, height, length.out = 100))
  }else{
    x_breaks <- c(0,Inf)
    y_breaks <- c(0,Inf)
  }
  
  if(! is.null(filename)){
    old_dev <- grDevices::dev.cur()
    if(latex_support){
      filename <- if(stringr::str_ends(filename, "\\.pdf")){
        paste0(stringr::str_sub(filename, end  = -5L), ".tex")
      }
      my_tikz(filename, width = width, height = height, units = units, stand_alone = TRUE)
    }else{
      my_pdf(filename, width = width, height = height, units = units)
      on.exit(utils::capture.output({
        grDevices::dev.off()
        if (old_dev > 1) grDevices::dev.set(old_dev)
      }))
    }
  }
  
  
  plotgardener::pageCreate(width = width, height = height, default.units = units, xgrid = diff(x_breaks)[1], ygrid = diff(y_breaks)[1], showGuides = show_grid_lines)
  for(obj in plots){
    if(is.ggplot(obj)){
      plotgardener::plotGG(obj, x = 0, y = 0, width = width, height = height, default.units = units)
    }else if(inherits(obj, "tikz_graphic_fun")){
      grid::grid.draw(obj$FUN(height))
    }else if(grid::is.grob(obj)){
      grid::grid.draw(obj)
    }else if(is.list(obj)){
      stopifnot(! is.null(names(obj)))
      stopifnot("plot" %in% names(obj))
      .x <- obj$x %||% 0
      .y <- obj$y %||% 0
      .width <- obj$width %||% width
      .height <- obj$height %||% height
      .units <- obj$units %||% units
      plotgardener::plotGG(obj$plot, x = .x, y = .y, width = .width, height = .height, default.units = .units)
    }else{
      warning("Cannot handle object of class: ", toString(class(obj)))
    }
  }
  
  if(! is.null(filename)){
    if(latex_support){
      grDevices::dev.off()
      if (old_dev > 1) grDevices::dev.set(old_dev)
      withr::with_dir(dirname(filename), {
        tools::texi2pdf(basename(filename), clean = TRUE)
        # Remove .tex file
        file.remove(basename(filename)) 
        raw_file_name <- tools::file_path_sans_ext(basename(filename))
        # rastered images
        ras_files <- list.files(pattern = paste0("^", raw_file_name, "_ras\\d*.png"))
        if(length(ras_files) > 0){
          file.remove(ras_files)
        }
      })
    }
  }
}

add_plot <- function(plot, x = 0, y = 0, width = NULL, height = NULL, units = NULL){
  list(plot = plot, x = x, y = y, width = width, height = height, units = units)
}

add_text <- function(label, x = 0, y = 0, fontsize = 12, hjust = 0, vjust = 0, ...){
  list(plot = cowplot::ggdraw() + cowplot::draw_label(label, size = fontsize, hjust = hjust, vjust = vjust, ...), x = x, y = y, width =0, height = 0)
}

# Note that x and y are from the lower left corner (instead of upper left :/)
add_graphic <- function(filename, x = 0, y = 0, width = NULL, height = NULL,
                        units = c("inches", "cm", "mm", "px", "user"),
                        anchor = c("north west", "south west", "base")){
  stopifnot(file.exists(filename))
  units <- match.arg(units)
  anchor <- anchor[1]
  abs_filepath <- tools::file_path_as_absolute(filename)
  size_spec <- if(!is.null(height) && !is.null(width)){
    paste0("[width=", width, units, ", height=", height, units, "]")
  }else if(!is.null(height)){
    paste0("[height=", height, units, "]")
  }else if(! is.null(width)){
    paste0("[width=", width, units, "]")
  }else{
    ""
  }
  content <- paste0(r"(\includegraphics)", size_spec, r"({")", abs_filepath, r"("})")
  res <- list(FUN = (\(figure_height){
    tikzDevice::grid.tikzNode(x = x, y = figure_height - y, units = units, 
                              opts = paste0("draw=none,fill=none,anchor=", anchor),  
                              content = content, draw = FALSE)
  }))
  class(res) <- "tikz_graphic_fun"
  res
}


correlationPlot2 <- function(mat, colTicks = NA, colLabs = NA, rowLabs = NA, title = NA, annotation = FALSE, annotation_size = 0.5 ){
  
  correlation_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209)
  correlation_breaks <- c(seq(-1,-0.01,length=100),
                          seq(-0.009,0.009,length=10),
                          seq(0.01,1,length=100))
  
  dat <- reshape2::melt(mat)
  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill=value)) +
    
    ggplot2::scale_y_discrete(breaks = as.character(dat$Var2), labels = as.character(dat$Var2))
  
  ## if correlation values are to be plotted in squares:
  if(annotation){
    plt <- plt + ggplot2::geom_text(ggplot2::aes(x = as.character(Var1), y = as.character(Var2), label = format(round(value, 2), nsmall = 2)), size = annotation_size)
  }
  
  plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, color = "black", hjust = 0, vjust = 0.5),
                              axis.text.y = ggplot2::element_text(size=12, color = "black"),
                              axis.title.y = ggplot2::element_text(size=13),
                              axis.title.x = ggplot2::element_text(size=13),
                              plot.title = ggplot2::element_text(size=15),
                              legend.text = ggplot2::element_text(size = 15, colour = "black"),
                              legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                              panel.background = ggplot2::element_blank(),
                              ## border around plot
                              panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                              plot.background = ggplot2::element_blank()
                              # legend.position="none"
  ) +
    ## fix up colorbar legend
    ggplot2::scale_fill_gradientn(limits = c(-1,1),
                                  breaks = c(-1,0,1),
                                  colors=(grDevices::colorRampPalette(c("blue","white","red")))(n = 209)
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Correlation",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    )) +
    
    scale_y_discrete(labels = colnames(mat)) +
    
    ggplot2::coord_fixed() 
  
  if (!is.na(colLabs)){
    plt <- plt + ggplot2::xlab(colLabs)
  }
  if (!is.na(rowLabs)){
    plt <- plt + ggplot2::ylab(rowLabs)
  }
  if (!is.na(title)){
    plt <- plt + ggplot2::ggtitle(title)
  }
  
  return(plt)
  
}

vizTopic2 <- function(theta, pos, topic,
                      groups = NA,
                      group_cols = NA,
                      size = 2,
                      stroke = 0.3,
                      alpha = 1,
                      low = "white",
                      high = "red",
                      plotTitle = NA,
                      showLegend = TRUE) {
  
  ## check that theta and pos are either data.frames or matrices
  if( !is.matrix(theta) & !is.data.frame(theta) == FALSE ){
    stop("`theta` must be a matrix or data.frame.")
  }
  if( !is.matrix(pos) == FALSE & !is.data.frame(pos) == FALSE ){
    stop("`pos` must be a matrix or data.frame with exactly 2 columns named `x` and `y`.")
  }
  
  if( (any(!colnames(pos) %in% c("x", "y")) == TRUE) | (dim(pos)[2] != 2) ){
    stop("`pos` must have exactly 2 columns named `x` and `y`.")
  }
  
  # ensure that `theta` and `pos` pixel rownames maintain same order
  # after the merge so as to not mess up the order of `groups`
  # if provided
  # make sure only using the shared pixels
  pixels <- intersect(rownames(theta), rownames(pos))
  pixels <- rownames(theta)[which(rownames(theta) %in% pixels)]
  
  proportion <- theta[,topic]
  dat <- merge(data.frame(proportion),
               data.frame(pos), by=0)
  
  rownames(dat) <- dat[,"Row.names"]
  ## make sure pixels in the original order before the merge
  dat <- dat[pixels,]
  
  # color spots by group:
  if (is.na(groups[1])) {
    Groups <- " "
  } else {
    Groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c(" " = "black")
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = dat, ggplot2::aes(x=x, y=y, fill=proportion, color = Groups),
                        shape = 21,
                        stroke = stroke, size = size, 
                        alpha = alpha) +
    ggplot2::scale_color_manual(values = group_cols)
  
  p <- p +
    ggplot2::theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8, colour = "black"),
      legend.title = ggplot2::element_text(size = 9, colour = "black")
    ) +
    
    ggplot2::scale_fill_gradientn(limits = c(0, 1.0),
                                  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                                  colors=(grDevices::colorRampPalette(c(low, high)))(n = 209)
    ) +
    
    ggplot2::guides(color = "none", fill = ggplot2::guide_colorbar(title = "Propotions",
                                                                   title.position = "left",
                                                                   title.hjust = 0.5,
                                                                   ticks.colour = "black",
                                                                   ticks.linewidth = 2,
                                                                   frame.colour= "black",
                                                                   frame.linewidth = 1.5,
                                                                   label.hjust = 0,
                                                                   title.theme = ggplot2::element_text(angle = 90)
    ))
  
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  p <- p + ggplot2::coord_equal()
  
  return(p)
  
}
