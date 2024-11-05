library(Seurat)
library(patchwork)
library(tidyverse)
# library(dplyr)

data_dir <- '/shares/rheumatologie.usz/caroline/spatial/SHK166_RA_Knee'
# data_dir <- '/shares/rheumatologie.usz/caroline/spatial'
list.files(data_dir) 
obj <- Load10X_Spatial(
    data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = NULL
)

view(obj@meta.data)
view(obj)

plot1 <- VlnPlot(obj, features = "nCount_Spatial", pt.size = 0.01) + NoLegend()
plot2 <- SpatialFeaturePlot(obj, features = "nCount_Spatial") + theme(legend.position = "right")

wrap_plots(plot1, plot2)

{r}
# SDY data
path <- '/home/myi/data'
obj_test <- Load10X_Spatial(
  path,
  filename = "adata_spatial.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL
)
dim(obj_test)

