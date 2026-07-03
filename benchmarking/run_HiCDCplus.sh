#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=./PBS/run_HiCDCPlus_%x.out
#SBATCH --error=./PBS/run_HiCDCPlus_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction_benchmarking
    echo "[INFO] Using conda environment diffraction_benchmarking. R 4.5.3"

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

HiCDCPlus_dir=$1
chrom=$2
resolution=$3
kspike=$4
p_val_threshold=$5


Rscript run_HiCDCPlus.R ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_A1.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_A2.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_B1.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_B2.table \
     ${HiCDCPlus_dir}/spikein_and_neighbors_coordinates_k${kspike}.txt ${p_val_threshold}