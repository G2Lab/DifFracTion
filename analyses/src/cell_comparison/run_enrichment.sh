#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=./PBS/statistical_significance_%x.out
#SBATCH --error=./PBS/statistical_significance_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction
    echo "[INFO] Using conda environment diffraction."

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

diffraction_type=$1
resolution=$2
log_fold=$3
cell_line=$4
marker=$5
bed_file=$6

cd ..

python -m analyses.src.cell_comparison.enrichment --diffraction_type ${diffraction_type} \
        --resolution ${resolution} \
        --log_fold ${log_fold} --tissue ${cell_line} \
        --marker ${marker} --bed_file ${bed_file}