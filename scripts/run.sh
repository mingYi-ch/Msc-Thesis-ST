#!/bin/bash
#SBATCH --job-name=generate_matrix    # Job name
#SBATCH --output=generate_matrix.out  # Standard output log
#SBATCH --error=generate_matrix.err   # Standard error log
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --time=00:10:00               # Time limit hrs:min:sec
#SBATCH --mem=1G                      # Memory limit
#SBATCH --partition=short             # Partition name

# Load the R module
module load mamba
# mamba activate .renv-s3it
R

# Run the R script
Rscript /home/myi/mt-spatial-transcriptomics/scripts/explore_data.R
