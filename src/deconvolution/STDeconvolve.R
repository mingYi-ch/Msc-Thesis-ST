library(STdeconvolve)
library(qs)

source("./src/util.R")
source("./config/config.R")

# # load data as spatialExperiment object
# data_dir <- get.std.path()
# spe <- load.data(data_dir)
# 
# cd <- assay(spe, "counts")
# 
# # QC: Filter out poor pixels and genes
# counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10, verbose = TRUE)
# 
# # Select features
# corpus <- restrictCorpus(counts, 
#                          removeAbove=0.95, 
#                          removeBelow = 0.05, 
#                          nTopOD = NA,
#                          plot = FALSE,
#                          alpha = 0.1) # 1e-8

data_dir <- file.path(get.data.dir(), "processed", "spe_subset.qs")
spe <- qread(data_dir)
corpus <- counts(spe)
dim(corpus)
# ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(6, 10))
ldas <- fitLDA(t(as.matrix(corpus)), Ks = c(7))

# # opt.K <- ldas$kneedOptK
# 
# # select the LDA model of interest and get the beta (cell-type transcriptional profiles) and theta (cell-type barcode proportions) matrices.
# optLDA <- optimalModel(models = ldas, opt = opt.K)
# results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)

res.dir <- get.res.dir()
res.dir
qs::qsave(ldas, file.path(res.dir, "deconvolution/STDeconvolve/ldas.qs"))
# qs::qsave(results, paste0("./results/deconvolution/STDeconvolve/deconRes_K", opt.K, ".qs"))
