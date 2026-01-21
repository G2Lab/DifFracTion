#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --output=./ds_benchmarking_%x.out
#SBATCH --error=./ds_benchmarking_%x.err


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
    echo "Usage: $0 --hic <HiC_file> --chrom <chromosome> --resolution <bp> --downsample_factor <factor> [options]"
    echo ""
    echo "Options:"
    echo "  -h,  --hic         Input Hi-C file (.hic format)"
    echo "  -c,  --chrom       Chromosome to process (e.g., 1)"
    echo "  -r,  --resolution  Resolution (bin size) in bp"
    echo "  -d,  --downsample_factor  Downsample factor. Default: 1.0"
    echo "  -p,  --pval        P-value threshold for significance."
    echo ""
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--hic) hic="$2"; shift 2 ;;
        -c|--chrom) chrom="$2"; shift 2 ;;
        -r|--resolution) resolution="$2"; shift 2 ;;
        -d|--downsample_factor) downsample_factor="$2"; shift 2 ;;
        -p|--pval) p_val_threshold="$2"; shift 2 ;;
        *) echo "[DifFracTion] [ERROR] Unknown parameter: $1"; usage ;;
    esac
done
# Check required parameters
if [[ -z "$hic" || -z "$chrom" || -z "$resolution" || -z "$p_val_threshold" ]]; then
    echo "[DifFracTion] [ERROR] Missing required parameters."
    usage
fi  

# Set default value for downsample_factor if not provided   
if [[ -z "$downsample_factor" ]]; then
    downsample_factor=1.0
fi
echo "[DifFracTion] [INFO] Starting benchmarking with the following parameters:"
echo "  Hi-C file: $hic"
echo "  Chromosome: $chrom"
echo "  Resolution: $resolution"
echo "  Downsample factor: $downsample_factor"
echo "  P-value threshold: $p_val_threshold"
echo "  Input files generation for tools will be stored in:"
echo "  Directory: ./results_downsample/{tool_name}/${chrom}_${resolution}_${downsample_factor}/input_files/"

echo "[DifFracTion] [INFO] Running benchmarking..."
echo "----------------------------------------"
echo "[DifFracTion] [INFO] Generating input files for benchmarking..."
echo "----------------------------------------"

python generate_inputs_downsample.py --hic $hic --chrom $chrom --resolution $resolution --downsample_factor $downsample_factor


echo "[DifFracTion] [INFO] Input files generation completed."
echo "----------------------------------------"

run_all_tools=true

if $run_all_tools; then
     echo "[DifFracTion] [INFO] Evaluating diffHiC..."
     echo "[DifFracTion] [INFO] ./results_downsample/diffHiC/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"

     diffHiC_dir="./results_downsample/diffHic/${chrom}_${resolution}_${downsample_factor}/input_files/"
     Rscript run_diffHiC_ds.R ${diffHiC_dir}/diffHic_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table ${p_val_threshold}

     echo "[DifFracTion] [INFO] diffHiC evaluation completed."
     echo "[DifFracTion] [INFO] diffHiC Output ./results_downsample/diffHic/${chrom}_${resolution}_${downsample_factor}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating HiCcompare."
     echo "[DifFracTion] [INFO] ./results_downsample/HiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"

     HiCcompare_dir="./results_downsample/HiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/"
     Rscript run_HiCcompare_ds.R ${HiCcompare_dir}/HiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}.table ${p_val_threshold}

     echo "[DifFracTion] [INFO] HiCcompare evaluation completed."
     echo "[DifFracTion] [INFO] HiCcompare Output ./results_downsample/HiCcompare/${chrom}_${resolution}_${downsample_factor}/performance/ "
     echo "----------------------------------------"



     echo "[DifFracTion] [INFO] Evaluating HiCDCPlus."
     echo "[DifFracTion] [INFO] ./results_downsample/HiCDCPlus/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"

     HiCDCPlus_dir="./results_downsample/HiCDCPlus/${chrom}_${resolution}_${downsample_factor}/input_files/"
     Rscript run_HiCDCPlus_ds.R ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A1.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A2.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B1.table \
     ${HiCDCPlus_dir}/HiCDCPlus_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B2.table \
     ${p_val_threshold}

     echo "[DifFracTion] [INFO] HiCDCPlus evaluation completed."
     echo "[DifFracTion] [INFO] HiCDCPlus Output ./results_downsample/HiCDCPlus/${chrom}_${resolution}_${downsample_factor}/performance/ "
     echo "----------------------------------------"

     echo "[DifFracTion] [INFO] Evaluating multiHiCcompare."
     echo "[DifFracTion] [INFO] ./results_downsample/multiHiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/ "
     echo "----------------------------------------"

     multiHiCcompare_dir="./results_downsample/multiHiCcompare/${chrom}_${resolution}_${downsample_factor}/input_files/"
     Rscript runmultiHiCcompare_ds.R ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A1.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_A2.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B1.table \
          ${multiHiCcompare_dir}/multiHiCcompare_input_chr${chrom}_res${resolution}_ds${downsample_factor}_B2.table \
          ${p_val_threshold}
     echo "[DifFracTion] [INFO] multiHiCcompare evaluation completed."
     echo "[DifFracTion] [INFO] multiHiCcompare Output ./results_downsample/multiHiCcompare/${chrom}_${resolution}_${downsample_factor}/performance/ "
     echo "----------------------------------------"
fi
