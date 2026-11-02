import hicstraw
import sys,os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from matplotlib.patches import Circle

import diffraction
from diffraction import utils as DifFracTion_utils
from diffraction import spikein as DifFracTion_spikes
from diffraction import src as DifFracTion

# Main scipt to generate DifFracTion results on synthetic datasets with spike-ins

def parse_args():
     
     parser = argparse.ArgumentParser(description="Plot the distance distribution of Hi-C contacts")
     
     parser.add_argument("--task", required=True, choices=["spikes"])

     parser.add_argument("-f", "--hic", required=True, help="Input Hi-C file (hic format)")
     parser.add_argument("-c", "--chrom", required=True, help="Chromosome to process (e.g., 1)")
     parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution (bin size) in bp")
     parser.add_argument("-m", "--metric", default="median", required=False,help="Metric used to fit (mean or median). Default: median")
     parser.add_argument("-u","--upper", type=int, default=7e6, required=False, help="Upper distance limit to consider (in bp). Default: 7e6")
     parser.add_argument("-l","--lower", type=int, default=0, required=False, help="Lower distance limit to consider (in bp). Default: 0")
     parser.add_argument("-p","--pvalue", type=float, default=0.01, required=False, help="P-adjusted-value threshold to consider a contact as positive. Default: 0.01")
     parser.add_argument("-d","--downsample", type=float, default=0.5, required=False, help="Downsampling factor for the second dataset. Default: 0.8")
     parser.add_argument("-n","--option", type=str, choices=['alpha_based','cycle'], default='alpha_based', required=False, help="Normalization option to use: 'alpha_based' or 'cycle'. Default: alpha_based")
     parser.add_argument("-o","--output", required=False, default='./spike_ins/',help="Output prefix for saving results. Default: .")
     parser.add_argument("-k","--k_spike", type=int, default=5, required=True, help="Spike-in fold change factor k. Default: 5")
     return parser.parse_args()

args = parse_args()

##### Inputs 
##### User inputs
hic_path  = args.hic
chromosome = str(args.chrom)
resolution = args.resolution
metric = args.metric
p_value_threshold = args.pvalue
downsampling_factor = args.downsample
k = args.k_spike
option = args.option #[alpha_based, cycle]

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
desired_alpha = -1.08

def task_spikes():
    print("Starting spikes task...")
    matrixZoom = DifFracTion_utils.prepare_MatrixZoomData(hic_path, chromosome, resolution)
    raw_matrix = DifFracTion_utils.rawMatrix(matrixZoom, chromosome_length)
    raw_matrix_sub = raw_matrix.copy()

    # Apply spike-ins
    _, spikein_coordinates = DifFracTion_spikes.select_spikeins(
        raw_matrix_sub,
        resolution,
        n_spikes=50,
        fold_change=k
    )

    altered_matrix, spikein_all_lists, all_neighbors = DifFracTion_spikes.apply_spikeins(
        raw_matrix_sub,
        spikein_coordinates,
        neighbor_degrees,
        k
    )

    # Save plot + spike metadata
    DifFracTion_spikes.save_spike_and_neighbors_coordinates(
        spikein_all_lists,
        all_neighbors,
        coordinate_file
    )

    plot_mat = DifFracTion_spikes.plot_original_spiked(raw_matrix_sub, altered_matrix)
    plot_mat.savefig(plot_spikes_file, bbox_inches="tight", dpi=350)

    # Dataset A
    matrix_1 = altered_matrix
    # Dataset B (downsampled from raw)
    matrix_2 = DifFracTion_utils.synthetic_datasets(raw_matrix, downsampling_factor)

    #Normalization

    matrix_1_ds,matrix_2_ds=DifFracTion_utils.get_counts_by_distance_both_densematrices(matrix_1,matrix_2,
    chromosome_length,resolution)
    d_1, y_1, m_1, b_1, alpha_1 = DifFracTion_utils.calculate_distance_decay(matrix_1_ds, metric=metric,
                             lower_limit=lower_limit,upper_limit=upper_limit)
    d_2, y_2, m_2, b_2, alpha_2 = DifFracTion_utils.calculate_distance_decay(matrix_2_ds, metric=metric,
                             lower_limit=lower_limit,upper_limit=upper_limit) 
    
    if option == 'alpha_based':
        matrix_1_scaled,matrix_2_scaled,scaled_counts_1,scaled_counts_2=DifFracTion.alpha_normalization(matrix_1_ds,matrix_2_ds,
                                                                 chromosome_length,resolution,
                                                                 desired_alpha=desired_alpha,
                                                                 upper_limit=upper_limit,
                                                                 lower_limit=lower_limit,
                                                                 metric=metric)
          
        maplot_refit=DifFracTion_utils.plot_MA(matrix_1_scaled,matrix_2_scaled,resolution,log2_fc_cutoff=0,upper_limit=upper_limit)
        maplot_refit.savefig(f"{output_path}/MA_plot_refit_k{k}_{p_value_threshold}_{resolution}.png",dpi=300)
     
    if option == 'cycle':
        matrix_1_scaled,matrix_2_scaled,scaled_counts_1,scaled_counts_2=DifFracTion.iterative_normalization(matrix_1,matrix_2,
                                                                 chromosome_length,resolution)
          
        maplot_refit=DifFracTion_utils.plot_MA(matrix_1_scaled,matrix_2_scaled,resolution,log2_fc_cutoff=0,upper_limit=upper_limit)
        maplot_refit.savefig(f"{output_path}/MA_plot_refit_k{k}_{p_value_threshold}_{resolution}.png",dpi=300)
          
    # Reference
    matrix_ref = matrix_1_scaled + matrix_2_scaled

    # Generate cool/mcool
    DifFracTion_utils.generate_cool_matrix(matrix_1_scaled, chromosome, resolution, cooler_file_1)
    DifFracTion_utils.generate_cool_matrix(matrix_2_scaled, chromosome, resolution, cooler_file_2)
    DifFracTion_utils.generate_cool_matrix(matrix_ref, chromosome, resolution, cooler_file_ref)

    print("Spikes task completed...")



def main():
    args = parse_args()

    if args.task == "spikes":
        task_spikes()

if __name__ == "__main__":
    main()