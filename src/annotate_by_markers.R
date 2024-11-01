library(Matrix)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

library(SpatialExperiment)
library(qs)
library(glue)

# source("./config/config.R")

# load data
dir <- get.data.dir()
res.dir <- get.res.dir()
# annotate Endothelial Cells
dir.data.processed <- file.path(dir, "processed")
spe <- qread(file.path(dir.data.processed, "spe_clean.qs"))

target.cell <- "Endothelial Cells"
markers <- c("VWF", "PECAM1")
thresholds <- c(2, 2)
n <- length(markers)
is.ec.spots <- list()

for (idx in seq(n)) {
  marker.sel <- markers[[idx]]
  threshold <- thresholds[[idx]]
  
  gexp.target <- counts(spe)[marker.sel, ]
  is.ec.spot <- (gexp.target > threshold)
  is.ec.spots <- c(is.ec.spots, list(is.ec.spot))
}
is.ec.spot <- do.call("&", is.ec.spots)

# save 
filename <- glue("ECs_manually_marker_final.qs")
fp <- file.path(res.dir, "deconvolution", filename)

# Create the directory if it doesn't exist
dir.create(dirname(fp), showWarnings = F, recursive = TRUE)
qsave(is.ec.spot, fp)


