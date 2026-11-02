#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --output=./PBS/compartment_%x.out
#SBATCH --error=./PBS/compartment_%x.err

if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate diffraction
    echo "[INFO] Using conda environment diffraction."
else
    echo "[ERROR] Conda is not installed or not in PATH. Please install Conda."
    exit 1
fi

hic_file=$1
resolution=$2
type=$3
genome_version=${4:-"hg19"}  # default to hg19 if not provided
flip=${5:-"auto"}   # pass "flip" to force flip, "keep" to force keep, default "auto" uses GC

name=$(basename ${hic_file} .hic)
if [ $type == "tissue" ]; then
    results_path="./tissue_comparison/compartments"
    TADS_path="./tissue_comparison/TADs"

elif [ $type == "cell_line" ]; then
    results_path="./cell_comparison/compartments"
    TADS_path="./cell_comparison/TADs"
else
    echo "[ERROR] Unknown type: ${type}. Must be 'tissue' or 'cell_line'."
    exit 1
fi
output_path=${results_path}

mkdir -p ${results_path}



COOL_OUT="${output_path}/${name}_${resolution}.cool"
potential_cool="$(realpath ${TADS_path}/${name}_${resolution}.cool)"

if [ -f ${potential_cool} ]; then
    ln -s ${potential_cool} ${COOL_OUT}
    echo "[INFO] Cooler file ${COOL_OUT} already exists. Skipping conversion."
fi 

if [ ! -f "${COOL_OUT}" ]; then
    echo "[INFO] Converting ${hic_file} to ${COOL_OUT} at resolution ${resolution}..."
HIC_FILE="${hic_file}" COOL_OUT="${output_path}/${name}_${resolution}.cool" RESOLUTION="${resolution}" python3 - <<'PYEOF'
import hicstraw, os

hic_file   = os.environ['HIC_FILE']
cool_file  = os.environ['COOL_OUT']
resolution = int(os.environ['RESOLUTION'])

hic = hicstraw.HiCFile(hic_file)

assert resolution in hic.getResolutions(), \
    f"{resolution} not in available resolutions: {hic.getResolutions()}"

chromosomes = [c for c in hic.getChromosomes() if c.name != "All"]

size_file = cool_file.replace('.cool', '.chrom.sizes')
txt_file  = cool_file.replace('.cool', '.txt')

with open(size_file, 'w') as fsize:
    for chrom in chromosomes:
        fsize.write(f"{chrom.name}\t{chrom.length}\n")

with open(txt_file, 'w') as fo:
    for i in range(len(chromosomes)):
        for j in range(i, len(chromosomes)):
            chrom1 = chromosomes[i].name
            chrom2 = chromosomes[j].name
            result = hicstraw.straw('observed', 'NONE', hic_file, chrom1, chrom2, 'BP', resolution)
            for r in result:
                fo.write(f"{chrom1}\t{r.binX}\t{r.binX + resolution}\t"
                         f"{chrom2}\t{r.binY}\t{r.binY + resolution}\t{r.counts}\n")

os.system(f"cooler load -f bg2 {size_file}:{resolution} {txt_file} {cool_file}")
os.remove(txt_file)
os.remove(size_file)
PYEOF
else
    echo "[INFO] Cooler file ${COOL_OUT} already exists. Skipping conversion."
fi

cooler balance ${output_path}/${name}_${resolution}.cool --force
##Compute GC content track for phasing using bedtools nuc
#echo "[INFO] Extracting bins from cooler file..."
cooler dump -t bins ${output_path}/${name}_${resolution}.cool | tail -n +2 | cut -f1-3 | awk 'OFS="\t" {print ($1 ~ /^chr/) ? $1 : "chr"$1, $2, $3}' > ${results_path}/${name}_${resolution}_bins.bed


#
if [ ${genome_version} == "hg19" ]; then
    genome_fasta="${output_path}/hg19.fa"
    if [ ! -f "${genome_fasta}" ]; then
        echo "[INFO] Downloading hg19 genome for GC content calculation."
        wget -O ${output_path}/hg19.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
        gunzip -f ${output_path}/hg19.fa.gz
    fi
elif [ ${genome_version} == "hg38" ]; then
    genome_fasta="${output_path}/hg38.fa"
    if [ ! -f "${genome_fasta}" ]; then
        echo "[INFO] Downloading hg38 genome for GC content calculation."
        wget -O ${output_path}/hg38.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        gunzip -f ${output_path}/hg38.fa.gz
    fi
else
    echo "[ERROR] Unsupported genome version: ${genome_version}. Use 'hg19' or 'hg38'."
    exit 1
fi

##$5 is the GC content column in the bedtools nuc output
#echo "[INFO] Computing GC content with bedtools nuc..."
bedtools nuc -fi ${genome_fasta} -bed ${results_path}/${name}_${resolution}_bins.bed | \
    tail -n +2 | awk 'OFS="\t" {print $1, $2, $3, $5}' > ${results_path}/${name}_${resolution}_gc.bedgraph
#
echo "[INFO] Running eigs-cis..."
cooltools eigs-cis ${output_path}/${name}_${resolution}.cool \
    -o ${results_path}/${name}_${resolution}_compartments

# Flip E1 sign so that positive = A compartment (correlated with GC)
#echo "[INFO] Orienting E1 by GC content..."
echo "[INFO] Orienting E1 by GC content (top vs bottom decile)..."
python3 - <<EOF
import pandas as pd

eigs = pd.read_csv("${results_path}/${name}_${resolution}_compartments.cis.vecs.tsv", sep="\t").dropna()
if not eigs['chrom'].iloc[0].startswith('chr'):
    eigs['chrom'] = 'chr' + eigs['chrom']

gc = pd.read_csv("${results_path}/${name}_${resolution}_gc.bedgraph", sep="\t", header=None, names=["chrom","start","end","GC"])
merged = eigs.merge(gc, on=["chrom","start","end"], how="left").dropna()

flip_arg = "${flip}"
if flip_arg == "flip":
    eigs["E1"] = -eigs["E1"]
    print("[INFO] Flipped E1 globally: forced by user")
elif flip_arg == "keep":
    print("[INFO] E1 kept as-is: forced by user")
else:
    for chrom, grp in merged.groupby('chrom'):
        top_gc    = grp["GC"].quantile(0.90)
        bottom_gc = grp["GC"].quantile(0.10)
        median_e1_high_gc = grp.loc[grp["GC"] >= top_gc, "E1"].median()
        median_e1_low_gc  = grp.loc[grp["GC"] <= bottom_gc, "E1"].median()
        print(f"[INFO] {chrom} — Median E1 top-10% GC: {median_e1_high_gc:.3f}, bottom-10% GC: {median_e1_low_gc:.3f}")
        if pd.isna(median_e1_high_gc) or pd.isna(median_e1_low_gc):
            print(f"[WARN] {chrom} — not enough GC data to orient, skipping")
            continue
        if median_e1_high_gc < median_e1_low_gc:
            eigs.loc[eigs['chrom'] == chrom, "E1"] *= -1
            print(f"[INFO] {chrom} — E1 flipped")
        else:
            print(f"[INFO] {chrom} — E1 orientation OK")

eigs.to_csv("${results_path}/${name}_${resolution}_compartments.cis.vecs.tsv", sep="\t", index=False)
EOF

if [ $? -ne 0 ]; then
    echo "Some error ocurred during the compartment calling process. Please check the logs for details."
    exit 1
elif [ ! -f "${results_path}/${name}_${resolution}_compartments.cis.vecs.tsv" ]; then
    echo "Compartment calling failed: output file not found."
    exit 1
else
    echo "[INFO] Compartment calling completed successfully."
    rm -f ${results_path}/${name}_${resolution}_bins.bed
    rm -f ${results_path}/${name}_${resolution}_gc.bedgraph
    #rm -f ${genome_fasta}
fi