# Execution
```

chromosomes=(20 21 22)
resolutions=(100000)
downsample=(1.0 0.75 0.5 0.25 0.1)
pval=(0.01 0.05)
k_spikes=(2 3 4)
method="cycle" #alpha_based

for k in ${k_spikes[@]}; do
     for p in ${pval[@]}; do
          for c in ${chromosomes[@]}; do
               for r in ${resolutions[@]}; do
                    for d in ${downsample[@]}; do
                         sbatch --job-name=DifFracTion_chr${c}_res${r}_ds${d}_p${p}_k${k}_${method} \
                         run_spikein_validation.sh \
                         -f "../../test_data/GM12878-HRC.hic" -c ${c} -r ${r} -d ${d} -p ${p} -k ${k} -n ${method}

                         echo "Chromosome ${c}, Resolution ${r}, Downsample ${d}, p-val ${p}"
                    done
               done
          done
     done
done
```


# Rs for high resolution
```

chromosomes=(20 21 22)
resolutions=(50000)
downsample=(1.0)
pval=(0.01)
k_spikes=(2 3 4)
method="cycle"

for k in ${k_spikes[@]}; do
     for p in ${pval[@]}; do
          for c in ${chromosomes[@]}; do
               for r in ${resolutions[@]}; do
                    for d in ${downsample[@]}; do
                         sbatch --job-name=DifFracTion_chr${c}_res${r}_ds${d}_p${p}_k${k}_${method} \
                         run_spikein_validation_rs.sh \
                         -f "../../test_data/GM12878-HRC.hic" -c ${c} -r ${r} -d ${d} -p ${p} -k ${k} -n ${method}

                    done
               done
          done
     done
done
```
