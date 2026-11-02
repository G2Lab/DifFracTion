import hicstraw
import diffraction
import sys,os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import numpy as np
import pandas as pd
import argparse


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
     parser.add_argument("-n","--norm_option", required=False, default='alpha_based', help="Normalization option: 'alpha_based' or 'cycle'. Default: alpha_based")
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
upper_limit = args.upper
lower_limit = args.lower
output_path = args.output #cd to spike_ins before running, ./results/ is created on shell
downsampling_factor = args.downsample
option = args.norm_option  # 'alpha_based' or 'cycle'
k = args.k_spike

neighbor_degrees = 2

### Create output directory and files
os.makedirs(output_path, exist_ok=True)
os.makedirs(f"{output_path}/{chromosome}_{p_value_threshold}_{resolution}_{downsampling_factor}_{k}_{option}", exist_ok=True)
     
output_path = output_path + f'/{chromosome}_{p_value_threshold}_{resolution}_{downsampling_factor}_{k}_{option}/'
coordinate_file =  f"{output_path}/spikein_and_neighbors_coordinates_k{k}_{p_value_threshold}_{resolution}.txt"
plot_spikes_file = f"{output_path}/spikein_and_neighbors_k{k}_{p_value_threshold}_{resolution}.png"
#######

desired_alpha = -1.08
chromosome_length = DifFracTion_utils.get_chromosome_length(hic_path, chromosome)
likelihood=False


def main():
    matrixZoom=DifFracTion_utils.prepare_MatrixZoomData(hic_path, chromosome, resolution)
    raw_matrix = DifFracTion_utils.rawMatrix(matrixZoom,chromosome_length)
    raw_matrix_sub = raw_matrix.copy()

    #Apply spike-ins
    _, spikein_coordinates = DifFracTion_spikes.select_spikeins(raw_matrix_sub,resolution,n_spikes=50,fold_change=k)
    altered_matrix,spikein_all_lists,all_neighbors = DifFracTion_spikes.apply_spikeins(raw_matrix_sub, spikein_coordinates, neighbor_degrees, k)

    #Save plot and coordinates
    DifFracTion_spikes.save_spike_and_neighbors_coordinates(spikein_all_lists, all_neighbors, coordinate_file)    
    plot_mat=DifFracTion_spikes.plot_original_spiked(raw_matrix_sub, altered_matrix)
    plot_mat.savefig(plot_spikes_file,bbox_inches="tight",dpi=350)
    
    ###################################################################
    #Dataset A
    matrix_1 = altered_matrix
    #Dataset B with downsampling
    matrix_2 = DifFracTion_utils.synthetic_datasets(raw_matrix,downsampling_factor)# factor or 1

    #Save MA plot
    maplot=DifFracTion_utils.plot_MA(matrix_1,matrix_2,resolution,log2_fc_cutoff=0,upper_limit=upper_limit)
    maplot_file = f"{output_path}/MA_plot_initial_k{k}_{p_value_threshold}_{resolution}.png"
    maplot.savefig(maplot_file,dpi=300)
    ###################################################################
    # For initial visualization and distance decay calculation

    #Counts on synthetic datasets
    synthetic_counts1,synthetic_counts2=DifFracTion_utils.get_counts_by_distance_both_densematrices(matrix_1,matrix_2,
                                                                                         chromosome_length,resolution)
    
    #Initial distance decay
    d_1, y_1, m_1, b_1, alpha_1 = DifFracTion_utils.calculate_distance_decay(synthetic_counts1, metric=metric,
                             lower_limit=lower_limit,upper_limit=upper_limit)
    d_2, y_2, m_2, b_2, alpha_2 = DifFracTion_utils.calculate_distance_decay(synthetic_counts2, metric=metric,
                             lower_limit=lower_limit,upper_limit=upper_limit) 
    
    #Save initial distance decay plot
    nofit_plt=DifFracTion_utils.decay_boxPlot(d_1, synthetic_counts1, m_1, b_1,
		   d_2, synthetic_counts2, m_2, b_2, upper_limit=upper_limit)
    nofit_plt.savefig(f"{output_path}/decay_boxplot_k{k}_{p_value_threshold}_{resolution}.png",dpi=300)
     ###################################################################

    if option == 'alpha_based':
          adjusted_matrix_1_v2,adjusted_matrix_2_v2,scaled_counts_1,scaled_counts_2=DifFracTion.alpha_normalization(matrix_1,matrix_2,
                                                                 chromosome_length,resolution,
                                                                 desired_alpha=desired_alpha,
                                                                 upper_limit=upper_limit,
                                                                 lower_limit=lower_limit,
                                                                 metric=metric)
          
          maplot_refit=DifFracTion_utils.plot_MA(adjusted_matrix_1_v2,adjusted_matrix_2_v2,resolution,log2_fc_cutoff=0,upper_limit=upper_limit)
          maplot_refit.savefig(f"{output_path}/MA_plot_refit_k{k}_{p_value_threshold}_{resolution}.png",dpi=300)
     
    elif option == 'cycle':
          adjusted_matrix_1_v2,adjusted_matrix_2_v2,scaled_counts_1,scaled_counts_2=DifFracTion.iterative_normalization(matrix_1,matrix_2,
                                                                 chromosome_length,resolution)
          
          maplot_refit=DifFracTion_utils.plot_MA(adjusted_matrix_1_v2,adjusted_matrix_2_v2,resolution,log2_fc_cutoff=0,upper_limit=upper_limit)
          maplot_refit.savefig(f"{output_path}/MA_plot_refit_k{k}_{p_value_threshold}_{resolution}.png",dpi=300)
          
    
    ##### Recalculate distance decay 
    d_1, y_1, m_1, b_1, alpha_1 = DifFracTion_utils.calculate_distance_decay(scaled_counts_1, upper_limit=upper_limit,
                                                                 lower_limit=lower_limit,metric=metric)
    d_2, y_2, m_2, b_2, alpha_2 = DifFracTion_utils.calculate_distance_decay(scaled_counts_2,upper_limit=upper_limit,
                                                                 lower_limit=lower_limit, metric=metric)

    fit_plot=DifFracTion_utils.decay_boxPlot(d_1, scaled_counts_1, m_1, b_1,
                   d_2, scaled_counts_2, m_2, b_2)
    
    fit_plot.savefig(f"{output_path}/decay_boxplot_refit_k{k}_{p_value_threshold}_{resolution}.png",dpi=300)
 

    if likelihood:
          dcis_bayesian=DifFracTion.identify_differential_interactions(adjusted_matrix_1_v2,adjusted_matrix_2_v2,resolution,
                                                                      method='log2fc', subsample_n=1000,
                                                                      log2_fc_cutoff=0,p_value_threshold=p_value_threshold,bayesian=likelihood)
          dcis_bayesian.to_csv(f"{output_path}/significant_contacts_{p_value_threshold}_{resolution}.tsv",sep='\t',index=False)
          significant_contacts=dcis_bayesian
    else:
          dcis_no_bayesian=DifFracTion.identify_differential_interactions(adjusted_matrix_1_v2,adjusted_matrix_2_v2,resolution,
                                                                      method='log2fc', subsample_n=1000,
                                                                      log2_fc_cutoff=0,p_value_threshold=p_value_threshold,bayesian=likelihood)
          dcis_no_bayesian.to_csv(f"{output_path}/significant_contacts_{p_value_threshold}_{resolution}.tsv",sep='\t',index=False)
          significant_contacts=dcis_no_bayesian
     
    if len(significant_contacts) == 0:
          print("No significant contacts found with the given p-value threshold.")
          significant_contacts.to_csv(f"{output_path}/supported_significant_contacts_{p_value_threshold}_{resolution}.tsv",sep='\t',index=False)

    else:
          significant_contacts_neighbors = significant_contacts[significant_contacts['neighbor_support'] == True]
          significant_contacts_neighbors.to_csv(f"{output_path}/supported_significant_contacts_{p_value_threshold}_{resolution}.tsv",sep='\t',index=False)

if __name__ == "__main__":
    main()