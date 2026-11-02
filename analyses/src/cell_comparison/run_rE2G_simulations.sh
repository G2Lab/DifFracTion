#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=./PBS/rE2G_simulations_%x.out
#SBATCH --error=./PBS/rE2G_simulations_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction
    echo "[INFO] Using conda environment diffraction."

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

while getopts "d:r:l:c:" opt; do
  case $opt in
    d) diffraction_type="$OPTARG" ;;
    r) resolution="$OPTARG" ;;
    l) log_fold="$OPTARG" ;;
    c) cell_line="$OPTARG" ;;
    *) echo "Usage: $0 -d <diffraction_type> -r <resolution> -l <log_fold> -c <cell_line>"
       exit 1 ;;
  esac
done    


cd ..
python -m analyses.src.cell_comparison.rE2G_simulations --diffraction_type $diffraction_type --resolution $resolution --log_fold $log_fold --cell_line $cell_line