#! /bin/bash

# # build
# cd build
# ./gen_build.sh -DCMAKE_INSTALL_PREFIX=/data/myi/.local
# make release
# make install

# set env
source activate std
export PATH=/data/myi/.local/bin:$PATH
export LD_LIBRARY_PATH=/data/myi/.local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH


# run
rm -r ~/data/mt-spatial-transcriptomics/results/Berglund_Result/
cd ~/data/mt-spatial-transcriptomics

/data/myi/.local/bin/std --file ~/data/mt-spatial-transcriptomics/data/processed/sp_counts_SHK166_RA_Knee.tsv --iter 2000  --types 7 --top 2000  --output results/Berglund_Result/

# store results
zcat results/Berglund_Result/experiment0000-contributions_spot_type.tsv.gz > results/deconvolution/std-poisson/deconProp.tsv
zcat results/Berglund_Result/contributions_gene_type.tsv.gz > results/deconvolution/std-poisson/contributions_gene_type.tsv
zcat results/Berglund_Result/expected-features.tsv.gz > results/deconvolution/std-poisson/expected-features.tsv
