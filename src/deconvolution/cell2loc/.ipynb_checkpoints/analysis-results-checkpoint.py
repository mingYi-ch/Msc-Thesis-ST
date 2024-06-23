import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import cell2location
from cell2location.utils import select_slide


adata_file = "/home/myi/mt-spatial-transcriptomics/src/deconvolution/cell2loc/cell2loc/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
print(adata_vis.uns.keys())

slide = select_slide(adata_vis, 'SHK')
print(type(slide))