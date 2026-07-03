#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=./PBS/benchmarking_%x.out
#SBATCH --error=./PBS/benchmarking_%x.err

### This script runs the benchmarking of DifFracTion against other tools (diffHiC, HiCcompare, HiCDCPlus, multiHiCcompare) using the generated input files. It evaluates the performance of each tool based on the generated spike-ins and their neighbors.
### The main use is the test for spike in detection

#mamba env create -f ../environment_benchmarking.yml

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"  # Initialize Conda properly
    conda activate diffraction_benchmarking
    echo "[INFO] Using conda environment diffraction_benchmarking. R 4.5.3"

else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

if ! command -v Rscript &> /dev/null; then
    echo "[ERROR] Rscript is not available. Please ensure R is installed and Rscript is in your PATH."
    exit 1
fi


check_packs=false

if $check_packs; then
    missing_pkg=false
    set +H  # disable bash history expansion so ! in R expressions works
    # Check required R packages are installed in their lib dirs
    declare -A R_PACKAGES=(
        ["diffHic"]="diffHic_lib"
        ["HiCcompare"]="HiCcompare_lib"
        ["HiCDCPlus"]="HiCDCPlus_lib"
        ["multiHiCcompare"]="multiHiCcompare_lib"
    )
    for pkg in "${!R_PACKAGES[@]}"; do
        if ! Rscript -e "if (!requireNamespace('$pkg', quietly=TRUE)) quit(status=1)" &> /dev/null; then
            echo "[ERROR] R package '$pkg' not found. Make sure to activate the correct conda environment ."
            missing_pkg=true
        else
            echo "[INFO] R package '$pkg' is available."
        fi
    done
    if $missing_pkg; then exit 1; fi
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
hic=$(realpath "$hic")
echo "[DifFracTion] [INFO] Starting benchmarking (spike in detection) with the following parameters:"
echo "  Hi-C file: $(realpath "$hic")"
echo "  Chromosome: $chrom"
echo "  Resolution: $resolution"
echo "  K-spike: $kspike"
echo "  P-value threshold: $p_val_threshold"
echo "  Input files generation for tools will be stored in:"
echo "  Directory: ../results_spikeins/{tool_name}/${chrom}_${resolution}_${kspike}/input_files/"

echo "[DifFracTion] [INFO] Running benchmarking (spike in detection) ..."
echo "----------------------------------------"
echo "[DifFracTion] [INFO] Generating input files for benchmarking (spike in detection) ..."
echo "----------------------------------------"

generate_inputs=true
run_all_tools=false

if $generate_inputs; then
    echo "[DifFracTion] [INFO] Generating input files for DifFracTion..."
    echo "[DifFracTion] [INFO] ../results_spikeins/DifFracTion/${chrom}_${resolution}_${kspike}/input_files/ "
    echo "----------------------------------------"
    cd ..
    python -m benchmarking.generate_inputs --hic $hic --chrom $chrom --resolution $resolution --k_spike $kspike
    echo "[DifFracTion] [INFO] Input files for DifFracTion generated."
    cd benchmarking
    echo "----------------------------------------"
else
    echo "[DifFracTion] [INFO] Skipping input file generation for DifFracTion. Assuming they are already generated in ../results_spikeins/DifFracTion/${chrom}_${resolution}_${kspike}/input_files/ "
    echo "----------------------------------------"
fi


echo "[DifFracTion] [INFO] Input files generation completed."
echo "----------------------------------------"



if $run_all_tools; then
     echo "[DifFracTion] [INFO] Starting evaluation of all tools (DifFraction, diffHiC, HiCcompare, HiCDCPlus, multiHiCcompare) for spike in detection..."
     echo "----------------------------------------"
     echo "[DifFracTion] [INFO] Evaluating DifFracTion..."
     echo "[DifFracTion] [INFO] ../results_spikeins/DifFracTion/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"
     diffraction_dir=$(realpath "../results_spikeins/DifFracTion/${chrom}_${resolution}_${kspike}/input_files/")
     mkdir -p $(dirname ${diffraction_dir})/performance/
     matrix_A=$(find ${diffraction_dir} -name "*spike*.npz" | head -1)
     matrix_B=$(find ${diffraction_dir} -name "*kb.npz" | head -1)
     sbatch run_DifFracTion.sh --matrixA ${matrix_A}  \
                               --matrixB ${matrix_B} \
                               --chrom $chrom \
                               --resolution $resolution \
                               --pval $p_val_threshold \
                               --type_norm alpha \
                               --adjusted_pvalues_method distance \
                               --output_dir $(dirname ${diffraction_dir})/performance/
     sbatch run_DifFracTion.sh --matrixA ${matrix_A}  \
                               --matrixB ${matrix_B} \
                               --chrom $chrom \
                               --resolution $resolution \
                               --pval $p_val_threshold \
                               --type_norm iterative \
                               --adjusted_pvalues_method distance \
                               --output_dir $(dirname ${diffraction_dir})/performance/
     echo "[DifFracTion] [INFO] DifFracTion jobs submitted."
     echo "[DifFracTion] [INFO] DifFracTion Output $(dirname ${diffraction_dir})/performance/"
     echo "----------------------------------------"
     echo "[DifFracTion] [INFO] Evaluating diffHiC..."
     echo "[DifFracTion] [INFO] ../results_spikeins/diffHiC/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"

     diffHiC_dir=$(realpath "../results_spikeins/diffHic/${chrom}_${resolution}_${kspike}/input_files/")
     
     sbatch run_diffHiC.sh ${diffHiC_dir} ${chrom} ${resolution} ${kspike} ${p_val_threshold}

     echo "[DifFracTion] [INFO] diffHiC evaluation completed."
     echo "[DifFracTion] [INFO] diffHiC Output ${diffHiC_dir}/results_spikeins/diffHic/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating HiCcompare."
     echo "[DifFracTion] [INFO] ../results_spikeins/HiCcompare/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"

     HiCcompare_dir=$(realpath "../results_spikeins/HiCcompare/${chrom}_${resolution}_${kspike}/input_files/")
     #sbatch run_HiCcompare.sh ${HiCcompare_dir} ${chrom} ${resolution} ${kspike} ${p_val_threshold}

     echo "[DifFracTion] [INFO] HiCcompare Output ${HiCcompare_dir}/results_spikeins/HiCcompare/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating HiCDCPlus."
     echo "[DifFracTion] [INFO] ../results_spikeins/HiCDCPlus/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"

     HiCDCPlus_dir=$(realpath "../results_spikeins/HiCDCPlus/${chrom}_${resolution}_${kspike}/input_files/")
     sbatch run_HiCDCplus.sh ${HiCDCPlus_dir} ${chrom} ${resolution} ${kspike} ${p_val_threshold}

     echo "[DifFracTion] [INFO] HiCDCPlus evaluation completed."
     echo "[DifFracTion] [INFO] HiCDCPlus Output ${HiCDCPlus_dir}/results_spikeins/HiCDCPlus/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating multiHiCcompare."
     echo "[DifFracTion] [INFO] ../results_spikeins/multiHiCcompare/${chrom}_${resolution}_${kspike}/input_files/ "
     echo "----------------------------------------"

     multiHiCcompare_dir=$(realpath "../results_spikeins/multiHiCcompare/${chrom}_${resolution}_${kspike}/input_files/")
     sbatch run_multiHiCcompare.sh ${multiHiCcompare_dir} ${chrom} ${resolution} ${kspike} ${p_val_threshold}
     echo "[DifFracTion] [INFO] multiHiCcompare evaluation completed."
     echo "[DifFracTion] [INFO] multiHiCcompare Output ${multiHiCcompare_dir}/results_spikeins/multiHiCcompare/${chrom}_${resolution}_${kspike}/performance/ "
     echo "----------------------------------------"

fi

# Some tools filter out spikes in their filtering steps, which leads to have negative values on TN, because we are calculating it as TN=n_tests - TP -FP - FN. And TP and FN are calculated based on the known spike-ns.
