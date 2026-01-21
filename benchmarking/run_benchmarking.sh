#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --output=./benchmarking_%x.out
#SBATCH --error=./benchmarking_%x.err


ml miniforge3
ml R/4.5.2-mba

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate DifFracTion
    echo "[INFO] Using conda environment DifFracTion"
else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

usage() {
    echo ""
    echo "Usage: $0 --hic <HiC_file> --chrom <chromosome> --resolution <bp> --k <kspike> [options]"
    echo ""
    echo "Options:"
    echo "  -h,  --hic         Input Hi-C file (.hic format)"
    echo "  -c,  --chrom       Chromosome to process (e.g., 1)"
    echo "  -r,  --resolution  Resolution (bin size) in bp"
    echo "  -p,  --pval        P-value threshold for significance."
    echo "  -k,  --kspike      Spike-in value for benchmarking. Default: 2"
    echo ""
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--hic) hic="$2"; shift 2 ;;
        -c|--chrom) chrom="$2"; shift 2 ;;
        -r|--resolution) resolution="$2"; shift 2 ;;
        -k|--kspike) kspike="$2"; shift 2 ;;
        -p|--pval) p_val_threshold="$2"; shift 2 ;;
        *) echo "[DifFracTion] [ERROR] Unknown parameter: $1"; usage ;;
    esac
done
# Check required parameters
if [[ -z "$hic" || -z "$chrom" || -z "$resolution" || -z "$p_val_threshold" ]]; then
    echo "[DifFracTion] [ERROR] Missing required parameters."
    usage
fi  

# Set default value for kspike if not provided   
if [[ -z "$kspike" ]]; then
    kspike=2
fi
echo "[DifFracTion] [INFO] Starting benchmarking with the following parameters:"
echo "  Hi-C file: $hic"
echo "  Chromosome: $chrom"
echo "  Resolution: $resolution"
echo "  K-spike: $kspike"
echo "  P-value threshold: $p_val_threshold"
echo "  Input files generation for tools will be stored in:"
echo "  Directory: ./results/{tool_name}/${chrom}_${resolution}_${kspike}/input_files/"

echo "[DifFracTion] [INFO] Running benchmarking..."
echo "----------------------------------------"
echo "[DifFracTion] [INFO] Generating input files for benchmarking..."
echo "----------------------------------------"

python generate_inputs.py --hic $hic --chrom $chrom --resolution $resolution --k_spike $kspike


echo "[DifFracTion] [INFO] Input files generation completed."
echo "----------------------------------------"

run_all_tools=true

if $run_all_tools; then
     echo "[DifFracTion] [INFO] Evaluating diffHiC..."
     echo "[DifFracTion] [INFO] ./results/diffHiC/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"

     diffHiC_dir="./results/diffHic/${chrom}_${resolution}_${kspike}/input_files/"
     Rscript run_diffHiC.R ${diffHiC_dir}/diffHic_input_chr${chrom}_res${resolution}_k${kspike}.table \
     ${diffHiC_dir}/spikein_and_neighbors_coordinates_k${kspike}.txt ${p_val_threshold}

     echo "[DifFracTion] [INFO] diffHiC evaluation completed."
     echo "[DifFracTion] [INFO] diffHiC Output ./results/diffHic/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating HiCcompare."
     echo "[DifFracTion] [INFO] ./results/HiCcompare/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"
     HiCcompare_dir="./results/HiCcompare/${chrom}_${resolution}_${kspike}/input_files/"
     Rscript run_HiCcompare.R ${HiCcompare_dir}/HiCcompare_input_chr${chrom}_res${resolution}_k${kspike}.table \
     ${HiCcompare_dir}/spikein_and_neighbors_coordinates_k${kspike}.txt ${p_val_threshold}
     echo "[DifFracTion] [INFO] HiCcompare evaluation completed."
     echo "[DifFracTion] [INFO] HiCcompare Output ./results/HiCcompare/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating HiCDCPlus."
     echo "[DifFracTion] [INFO] ./results/HiCDCPlus/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"

     HiCDCPlus_dir="./results/HiCDCPlus/${chrom}_${resolution}_${kspike}/input_files/"
     Rscript run_HiCDCPlus.R ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_A1.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_A2.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_B1.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_k${kspike}_B2.table \
     ${HiCDCPlus_dir}/spikein_and_neighbors_coordinates_k${kspike}.txt ${p_val_threshold}
     echo "[DifFracTion] [INFO] HiCDCPlus evaluation completed."
     echo "[DifFracTion] [INFO] HiCDCPlus Output ./results/HiCDCPlus/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating multiHiCcompare."
     echo "[DifFracTion] [INFO] ./results/multiHiCcompare/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"

     multiHiCcompare_dir="./results/multiHiCcompare/${chrom}_${resolution}_${kspike}/input_files/"
     Rscript runmultiHiCcompare.R ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_k${kspike}_A1.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_k${kspike}_A2.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_k${kspike}_B1.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_k${kspike}_B2.table \
          ${multiHiCcompare_dir}/spikein_and_neighbors_coordinates_k${kspike}.txt ${p_val_threshold}
     echo "[DifFracTion] [INFO] multiHiCcompare evaluation completed."
     echo "[DifFracTion] [INFO] multiHiCcompare Output ./results/multiHiCcompare/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"
fi
