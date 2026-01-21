# Test on different resolutions / Downsamples
cd synthetic_resolution

Folder names are chromosome_pvalue_resolution_downsample_normalization

     chromosomes=(20 21 22)
     resolutions=(100000 250000 500000)
     downsample=(1.0 0.75 0.5 0.25 0.1)
     pval=(0.01)
     method="cycle" #alpha_based

     for p in ${pval[@]}; do
          for c in ${chromosomes[@]}; do
               for r in ${resolutions[@]}; do
                    for d in ${downsample[@]}; do
                         sbatch --job-name=HiCDip_chr${c}_res${r}_ds${d}_p${p}_${method} \
                         run_synthetic_validation_resolution.sh \
                         -f "../../test_data/GM12878-HRC.hic" -c ${c} -r ${r} -d ${d} -p ${p} -n ${method}

                         echo "Chromosome ${c}, Resolution ${r}, Downsample ${d}, p-val ${p}"
                    done
               done
          done
     done


### Test trial 
c=21
r=500000
d=0.8
p=0.01
method='cycle'

sbatch --job-name="DiFracTion_chr${c}_res${r}_ds${d}_p${p}_${method}" run_synthetic_validation_resolution.sh \
          -f "../../../data/HiC-Dip/test_data/GM12878-HRC.hic" -c ${c} -r ${r} -d ${d} -p ${p} -n ${method}
