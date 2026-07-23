#!/bin/bash
## Submits only the missing DifFracTion performance files.
## Missing both alpha+iterative: all combinations below except the 1.0 iterative-only cases.
## Missing iterative only: chr1@50k@1.0, chr2@50k@1.0

hic_path=$(realpath "../test_data/GM12878-HRC.hic")
pval=0.01

submit() {
    local chrom=$1
    local resolution=$2
    local ds=$3
    local norm=$4
    echo "Submitting: chr${chrom} res${resolution} ds${ds} norm=${norm}"
    sbatch --job-name=ds_Benchmarking_DifFracTion_chr${chrom}_res${resolution}_ds${ds}_${norm} \
        run_benchmarking_downsample.sh \
        -f ${hic_path} -c ${chrom} -r ${resolution} -d ${ds} -p ${pval} -s diffraction -t ${norm}
}

# --- Missing iterative only (alpha already exists) ---
submit 1  50000  1.0  iterative
submit 2  50000  1.0  iterative

# --- Missing both alpha + iterative ---
submit 1  250000 0.5  both
submit 2  250000 0.25 both
submit 5  250000 0.25 both
submit 5  250000 0.5  both
submit 5  250000 0.75 both
submit 6  250000 0.1  both
submit 6  250000 0.25 both
submit 7  100000 0.25 both
submit 8  100000 0.1  both
submit 8  100000 0.25 both
submit 10 100000 1.0  both
submit 14 100000 0.5  both
submit 17 50000  0.5  both
submit 17 50000  0.75 both