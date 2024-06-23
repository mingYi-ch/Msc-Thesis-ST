#! /bin/bash

cd build
./gen_build.sh -DCMAKE_INSTALL_PREFIX=/data/myi/.local
make release
make install

source activate std
export PATH=/data/myi/.local/bin:$PATH
export LD_LIBRARY_PATH=/data/myi/.local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH



/data/myi/.local/bin/std --file ~/data/spatial-trancriptomics/st_counts_SHK166_RA_Knee.tsv --iter 1000  --types 7 --top 3000 --output ./Berglund_Result/
zcat Berglund_Result/experiment0000-contributions_spot_type.tsv.gz > results/deconvolution/std-poisson/deconProp.tsv
zcat Berglund_Result/contributions_gene_type.tsv.gz > results/deconvolution/std-poisson/contributions_gene_type.tsv
zcat Berglund_Result/expected-features.tsv.gz > results/deconvolution/std-poisson/expected-features.tsv
