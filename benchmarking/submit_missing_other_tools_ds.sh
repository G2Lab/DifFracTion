#!/bin/bash
## Submits missing combinations for HiCcompare, diffHiC, HiCDCPlus, multiHiCcompare.
## Each job runs only the specified tool via -s flag.

hic_path=$(realpath "../test_data/GM12878-HRC.hic")
pval=0.01

submit() {
    local tool=$1 chrom=$2 resolution=$3 ds=$4
    echo "Submitting: ${tool} chr${chrom} res${resolution} ds${ds}"
    sbatch --job-name=ds_Benchmarking_${tool}_chr${chrom}_res${resolution}_ds${ds} \
        run_benchmarking_downsample.sh \
        -f ${hic_path} -c ${chrom} -r ${resolution} -d ${ds} -p ${pval} -s ${tool}
}

#  HiCcompare 
submit HiCcompare 1  250000 0.5
submit HiCcompare 2  250000 0.25
submit HiCcompare 5  250000 0.25
submit HiCcompare 5  250000 0.5
submit HiCcompare 5  250000 0.75
submit HiCcompare 6  250000 0.1
submit HiCcompare 6  250000 0.25
submit HiCcompare 7  100000 0.25
submit HiCcompare 8  100000 0.1
submit HiCcompare 8  100000 0.25
submit HiCcompare 10 100000 1.0
submit HiCcompare 14 100000 0.5
submit HiCcompare 17 50000  0.5
submit HiCcompare 17 50000  0.75

#  diffHic 
submit diffHic 1  250000 0.5
submit diffHic 2  250000 0.25
submit diffHic 5  250000 0.25
submit diffHic 5  250000 0.5
submit diffHic 5  250000 0.75
submit diffHic 6  250000 0.1
submit diffHic 6  250000 0.25
submit diffHic 7  100000 0.25
submit diffHic 8  100000 0.1
submit diffHic 8  100000 0.25
submit diffHic 10 100000 1.0
submit diffHic 14 100000 0.5
submit diffHic 17 50000  0.5
submit diffHic 17 50000  0.75

#  multiHiCcompare 
submit multiHiCcompare 1  250000 0.5
submit multiHiCcompare 2  250000 0.25
submit multiHiCcompare 5  250000 0.25
submit multiHiCcompare 5  250000 0.5
submit multiHiCcompare 5  250000 0.75
submit multiHiCcompare 6  250000 0.1
submit multiHiCcompare 6  250000 0.25
submit multiHiCcompare 7  100000 0.25
submit multiHiCcompare 8  100000 0.1
submit multiHiCcompare 8  100000 0.25
submit multiHiCcompare 10 100000 1.0
submit multiHiCcompare 14 100000 0.5
submit multiHiCcompare 17 50000  0.5
submit multiHiCcompare 17 50000  0.75

#  HiCDCPlus 
submit HiCDCPlus 1  250000 0.5
submit HiCDCPlus 2  50000  0.75
submit HiCDCPlus 2  100000 0.1
submit HiCDCPlus 2  250000 0.25
submit HiCDCPlus 3  50000  0.75
submit HiCDCPlus 3  50000  1.0
submit HiCDCPlus 3  250000 0.1
submit HiCDCPlus 3  250000 0.25
submit HiCDCPlus 3  250000 0.5
submit HiCDCPlus 3  250000 0.75
submit HiCDCPlus 3  250000 1.0
submit HiCDCPlus 4  100000 0.1
submit HiCDCPlus 4  250000 0.1
submit HiCDCPlus 4  250000 0.25
submit HiCDCPlus 4  250000 0.5
submit HiCDCPlus 4  250000 0.75
submit HiCDCPlus 4  250000 1.0
submit HiCDCPlus 5  100000 0.25
submit HiCDCPlus 5  100000 0.75
submit HiCDCPlus 5  250000 0.25
submit HiCDCPlus 5  250000 0.5
submit HiCDCPlus 5  250000 0.75
submit HiCDCPlus 6  50000  0.25
submit HiCDCPlus 6  50000  0.5
submit HiCDCPlus 6  50000  0.75
submit HiCDCPlus 6  50000  1.0
submit HiCDCPlus 6  100000 0.75
submit HiCDCPlus 6  250000 0.1
submit HiCDCPlus 6  250000 0.25
submit HiCDCPlus 7  100000 0.25
submit HiCDCPlus 7  250000 0.1
submit HiCDCPlus 7  250000 0.25
submit HiCDCPlus 7  250000 0.5
submit HiCDCPlus 7  250000 0.75
submit HiCDCPlus 7  250000 1.0
submit HiCDCPlus 8  100000 0.1
submit HiCDCPlus 8  100000 0.25
submit HiCDCPlus 10 100000 0.1
submit HiCDCPlus 10 100000 0.75
submit HiCDCPlus 10 100000 1.0
submit HiCDCPlus 12 50000  0.5
submit HiCDCPlus 12 50000  0.75
submit HiCDCPlus 12 50000  1.0
submit HiCDCPlus 14 100000 0.1
submit HiCDCPlus 14 100000 0.5
submit HiCDCPlus 17 50000  0.5
submit HiCDCPlus 17 50000  0.75
submit HiCDCPlus 22 50000  0.1
submit HiCDCPlus 22 50000  0.25
submit HiCDCPlus 22 50000  0.5
submit HiCDCPlus 22 50000  0.75
submit HiCDCPlus 22 50000  1.0
submit HiCDCPlus 22 100000 0.1
submit HiCDCPlus 22 100000 0.25
submit HiCDCPlus 22 100000 0.5
submit HiCDCPlus 22 100000 0.75
submit HiCDCPlus 22 100000 1.0