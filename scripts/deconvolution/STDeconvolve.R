library(STdeconvolve)
library(qs)
source("./src/util.R")

# load data as spatialExperiment object
data_dir <- '/shares/rheumatologie.usz/caroline/spatial'
se <- load.data(dir)

cd <- assay(se, "counts")

# QC: Filter out poor pixels and genes
counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10, verbose = TRUE)

# Select features
corpus <- restrictCorpus(counts, 
                         removeAbove=0.95, 
                         removeBelow = 0.05, 
                         nTopOD = NA,
                         plot = FALSE,
                         alpha = 0.05) # 1e-8

ldas <- fitLDA(t(as.matrix(corpus)), Ks = 8:20)

opt.K <- 10
# select the LDA model of interest and get the beta (cell-type transcriptional profiles) and theta (cell-type barcode proportions) matrices.
optLDA <- optimalModel(models = ldas, opt = opt.K)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)

qs::qsave(ldas, "./results/deconvolution/STDeconvolve/ldas.qs")
qs::qsave(results, "./results/deconvolution/STDeconvolve/deconRes.qs")

