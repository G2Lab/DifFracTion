#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --mem=80G
#SBATCH --output=./DifFracTion_validation_%x.out
#SBATCH --error=./DifFracTion_validation_%x.err

ml miniforge3
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
    echo "Usage: $0 --hic <HiC_file> --chrom <chromosome> --resolution <bp> [options]"
    echo ""
    echo "Options:"
    echo "  -h,  --hic         Input Hi-C file (.hic format)"
    echo "  -n , --normalization Type of normalization (alpha_based or cycle)"
    echo "  -c,  --chrom       Chromosome to process (e.g., 1)"
    echo "  -r,  --resolution  Resolution (bin size) in bp"
    echo "  -m,  --metric      Metric used to fit (mean or median). Default: mean"
    echo "  -u,  --upper       Upper distance limit to consider (in bp). Default: 7e6"
    echo "  -l,  --lower       Lower distance limit to consider (in bp). Default: 0"
    echo "  -p,  --pvalue      P-adjusted-value threshold for positive contacts. Default: 0.01"
    echo "  -o,  --output      Output directory for results. Default: ./synthetic_validation/"
    echo "  -d,  --downsample  Downsample factor (e.g., 0.1 for 10%). Default:0.8 "
    echo ""
    exit 1
}

# Default values
metric="median"
upper=7000000
lower=0
pvalue=0.01
output="./results/" # Default output directory after cd to synthetic_resolution
mkdir -p $output
downsample=0.8
normalization="alpha_based"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--hic) hic="$2"; shift 2 ;;
        -n|--normalization) normalization="$2"; shift 2 ;;
        -c|--chrom) chrom="$2"; shift 2 ;;
        -r|--resolution) resolution="$2"; shift 2 ;;
        -m|--metric) metric="$2"; shift 2 ;;
        -u|--upper) upper="$2"; shift 2 ;;
        -l|--lower) lower="$2"; shift 2 ;;
        -p|--pvalue) pvalue="$2"; shift 2 ;;
        -d|--downsample) downsample="$2"; shift 2 ;;
        -o|--output) output="$2"; shift 2 ;;
        *) echo "[DifFracTion] [ERROR] Unknown parameter: $1"; usage ;;
    esac
done

# Check required arguments
if [[ -z "$hic" || -z "$chrom" || -z "$resolution" ]]; then
    echo "[DifFracTion] [ERROR] Missing required arguments."
    usage
fi

output_file="${output}/${chrom}_${pvalue}_${resolution}_${downsample}_${normalization}/synthetic_validation.log"
mkdir -p $(dirname ${output_file})

echo "[DifFracTion] Parameters:" > $output_file
echo "  Hi-C file:       $hic" >> $output_file
echo "  Normalization:   $normalization" >> $output_file
echo "  Chromosome:      $chrom" >> $output_file
echo "  Resolution:      $resolution" >> $output_file
echo "  Metric:          $metric" >> $output_file
echo "  Upper limit:     $upper" >> $output_file
echo "  Lower limit:     $lower" >> $output_file
echo "  P-value:         $pvalue" >> $output_file
echo "  Downsample:      $downsample" >> $output_file
echo "  Output:          $output/${chrom}_${pvalue}_${resolution}_${downsample}_${normalization}"  >> $output_file
echo "" >> $output_file

python run_synthetic_validation_resolution.py \
    --hic $hic \
    --norm_option $normalization \
    --chrom $chrom \
    --resolution $resolution \
    --metric $metric \
    --upper $upper \
    --lower $lower \
    --pvalue $pvalue \
    --output $output \
    --downsample $downsample
