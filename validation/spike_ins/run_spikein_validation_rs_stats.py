import hicstraw
import sys,os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from matplotlib.patches import Circle
sys.stdout.flush()

import diffraction
from diffraction import utils as DifFracTion_utils
from diffraction import spikein as DifFracTion_spikes
from diffraction import src as DifFracTion

def parse_args():
     
     parser = argparse.ArgumentParser(description="Plot the distance distribution of Hi-C contacts")
     
     parser.add_argument("-f", "--hic", required=True, help="Input Hi-C file (hic format)")
     parser.add_argument("-c", "--chrom", required=True, help="Chromosome to process (e.g., 1)")
     parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution (bin size) in bp")
     parser.add_argument("-m", "--metric", default="median", required=False,help="Metric used to fit (mean or median). Default: median")
     parser.add_argument("-u","--upper", type=int, default=7e6, required=False, help="Upper distance limit to consider (in bp). Default: 7e6")
     parser.add_argument("-l","--lower", type=int, default=0, required=False, help="Lower distance limit to consider (in bp). Default: 0")
     parser.add_argument("-p","--pvalue", type=float, default=0.01, required=False, help="P-adjusted-value threshold to consider a contact as positive. Default: 0.01")
     parser.add_argument("-d","--downsample", type=float, default=0.5, required=False, help="Downsampling factor for the second dataset. Default: 0.8")
     parser.add_argument("-o","--output", required=False, default='./spike_ins/',help="Output prefix for saving results. Default: .")
     parser.add_argument("-n","--option", type=str, choices=['alpha_based','cycle'], default='alpha_based', required=False, help="Normalization option to use: 'alpha_based' or 'cycle'. Default: alpha_based")
     parser.add_argument("-k","--k_spike", type=int, default=5, required=True, help="Spike-in fold change factor k. Default: 5")
     parser.add_argument("-e","--rust_output", required=False, help="Path to Rust output CSV file from permutation test. Default: None")
     return parser.parse_args()

args = parse_args()

##### Inputs 
##### This version removes 'option' as we don't do normalization, just differential analysis after spike-ins
##### User inputs
hic_path  = args.hic
chromosome = str(args.chrom)
resolution = args.resolution
metric = args.metric
p_value_threshold = args.pvalue
downsampling_factor = args.downsample
k = args.k_spike
option= args.option  # [alpha_based, cycle]

upper_limit = args.upper
lower_limit = args.lower
output_path = args.output 
neighbor_degrees = 2

### Create output directory and files
os.makedirs(output_path, exist_ok=True)
os.makedirs(f"{output_path}/{chromosome}_{p_value_threshold}_{resolution}_{downsampling_factor}_{k}_{option}", exist_ok=True)
     
output_path = output_path + f'/{chromosome}_{p_value_threshold}_{resolution}_{downsampling_factor}_{k}_{option}/'
coordinate_file =  f"{output_path}/spikein_and_neighbors_coordinates_k{k}_{p_value_threshold}_{resolution}.txt"
plot_spikes_file = f"{output_path}/spikein_and_neighbors_k{k}_{p_value_threshold}_{resolution}.png"

#######
# cool files
cooler_file_1 = output_path + f'chr{chromosome}_M1.cool'
cooler_file_2 = output_path + f'chr{chromosome}_M2.cool'
cooler_file_ref = output_path + f'chr{chromosome}_REF.cool'
#######

# Rust output
chromosome_length = DifFracTion_utils.get_chromosome_length(hic_path, chromosome)
rust_output = args.rust_output


print("Starting stats task...",flush=True) 
permutation_df = pd.read_csv(rust_output)
print(permutation_df.head(),flush=True)
permutation_df = DifFracTion_utils.reformat_permutation_df(permutation_df,resolution,upper_limit,lower_limit)  
permutation_df=DifFracTion_utils.adjust_p_values(permutation_df)

# Since we only care about recovering the spike-ins, we do not use Bayesian approach here
significant_contacts=DifFracTion_utils.significant_interactions_no_bayesian(permutation_df,p_value_threshold,upper_limit)
significant_contacts.to_csv(f"{output_path}/significant_contacts_{p_value_threshold}_{resolution}.tsv",sep='\t',index=False)

if len(significant_contacts) == 0:
      print("No significant contacts found with the given p-value threshold.")
      significant_contacts.to_csv(f"{output_path}/supported_significant_contacts_{p_value_threshold}_{resolution}.tsv",sep='\t',index=False)
else:
      significant_contacts=DifFracTion_utils.neighbor_support(significant_contacts) #matrix is just for coordinates, we don't modify values here
      significant_contacts_neighbors = significant_contacts[significant_contacts['neighbor_support'] == True]
      significant_contacts_neighbors.to_csv(f"{output_path}/supported_significant_contacts_{p_value_threshold}_{resolution}.tsv",sep='\t',index=False)
print("Stats task completed...")