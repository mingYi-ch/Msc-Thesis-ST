# Counterfactual-prediction based analysis of spatial transcriptomics

This study explores the application of a novel method called LEMUR (*Ahlmann-Eltze and Huber (2023a)*) on spatial transcriptomics, emphasizing counterfactual prediction and differential gene expression analysis without relying on predefined discrete clusters. Utilizing this approach allows for the analysis of data-driven neighborhood clusters through counterfactual predictions of the top 200 marker genes of synovial fibroblast, revealing likely four synovial fibroblast subtypes that correspond well with known studies from scRNA-seq data (*Micheroli et al. (2022)*). Employing a data transformation technique validated by (*Ahlmann-Eltze and Huber (2023b)*) with excellent performance, on our reference scRNA-seq data set, we also demonstrate the dynamic changes in UMAP visualizations with different number of principal component inputs, elucidating the relationships among various cell and sub-cell types. The detailed report can be found on [overleaf](https://www.overleaf.com/read/yczfttmfrxfz#c1f5db).

## Environment setup

### Analysis

We use R (version 4.4) for the main analysis with the `renv` setup. Please refer to the [documentation](https://rstudio.github.io/renv/articles/renv.html) about how to restore the environment from `renv.lock` file.

### Deconvolution

With respect to the environments for the deconvolution methods, please refer to the respective git repository: [std-poisson,](https://github.com/SpatialTranscriptomicsResearch/std-poisson/issues) [RCTD](https://rdrr.io/github/dmcable/RCTD/man/spacexr.html), [STDeconvolve](https://www.bioconductor.org/packages/release/bioc/html/STdeconvolve.html), [cell2loc](https://cell2location.readthedocs.io/en/latest/) and [Celloscope](https://github.com/szczurek-lab/Celloscope)

## Run code

### Prepare data

The directory for loading data and outputing results can be set in `config/config.R`. Please [contact the author](yiming.civi@gmail.com) for the data used in this project.

Run the RMarkdown `analysis-notebook/scRNA-seq-analysis/analysis-scRNA-seq-2.Rmd` to process the reference data set and analyze it.

To pre-process spatial transcriptomics data, run

```         
src/preprocess.R
```

It will store pre-processed the data under the folder `data/preprocessed`.

### Run analysis

-   Firstly, run `analysis-notebook/prepare-LEMUR.Rmd` , this will annotate the spots containing endocidial cells using the marker genes by default, and perform counterfactual prediction as well as the neighborhood clustering.

-   To examine the neighborhood cluster of a specific gene, run `analysis-notebook/synovial-fibroblast-analysis/analysis-synFib-LEMUR.Rmd.`

-   Run `analysis-notebook/synovial-fibroblast-analysis/analysis-curve-diff.Rmd` to do the clustering on marker genes of synovial fibroblast.

Results are stored in the `results` folder.

## Reference

1.  [LEMUR Repository](https://github.com/const-ae/lemur)
>>>>>>> develop
