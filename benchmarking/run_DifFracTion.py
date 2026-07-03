import hicstraw
import scipy
from scipy.sparse import coo_matrix
from scipy import ndimage
from scipy.stats import linregress
from statsmodels.stats.multitest import multipletests

import numpy as np
import pandas as pd
import os
import glob

import time
import tracemalloc
import argparse

import matplotlib.pyplot as plt
from matplotlib import colors

from diffraction_edits import utils as DifFracTion_utils
from diffraction_edits import src as DifFracTion

def parse_args():
    parser = argparse.ArgumentParser(description='Run DifFracTion on a given dataset')
    parser.add_argument('--matrixA', type=str, required=True, help='Path to the dataset (matrix or cooler file)')
    parser.add_argument('--matrixB', type=str, required=True, help='Path to the second dataset (matrix or cooler file)')
    parser.add_argument('--chrom', type=str, required=True, help='Chromosome to analyze') 
    parser.add_argument('--resolution', type=int, default=10000, help='Resolution for the analysis (default: 10000)')
    parser.add_argument('--pval', type=float, default=0.05, help='FDR threshold for significant interactions (default: 0.05)')
    parser.add_argument('--type_norm', type=str, default='alpha', help='Type of normalization to apply alpha or iterative (default: alpha)')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the results')
    parser.add_argument('--adjusted_pvalues_method', type=str, default='distance', help='Method for adjusting p-values (default: distance)')
    
    return parser.parse_args()


def get_confusion_matrix(spike_in_files, alpha_padjusted_interactions, n_bins):
    spike_in_coordinates = spike_in_files[0] if spike_in_files else None
    print(f"Spike-in coordinates file: {spike_in_coordinates}")
    spikes_df = pd.read_csv(spike_in_coordinates, sep='\t',header=None) if spike_in_coordinates else None
    if spikes_df is not None:
        merged = pd.merge(alpha_padjusted_interactions, spikes_df, left_on=['row', 'col'], right_on=[0, 1], how='outer', indicator=True)
        TP = len(merged[merged['_merge'] == 'both'])
        FP = len(merged[merged['_merge'] == 'left_only']) # if it was only identified in the interactions but not in the spike-in coordinates, it's a false positive
        FN = len(merged[merged['_merge'] == 'right_only']) # if it was only in the spike-in coordinates but not identified in the interactions, it's a false negative
        TN = n_bins - TP - FP - FN # the rest of the tested interactions that were not identified and are not in the spike-in coordinates are true negatives
        confusion_m = pd.DataFrame([{'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN}])
        return confusion_m
    else:
        print("No spike-in coordinates file found. Cannot compute confusion matrix.")
        return None
    
if __name__ == "__main__":
     args = parse_args()
     print(f"Running DifFracTion with the following parameters:\nMatrix A: {args.matrixA}\nMatrix B: {args.matrixB}\nChromosome: {args.chrom}\nResolution: {args.resolution}\nP-value threshold: {args.pval}\nNormalization type: {args.type_norm}\nAdjusted p-values method: {args.adjusted_pvalues_method}\nOutput directory: {args.output_dir}")
     matrix_A=DifFracTion_utils.load_dense_matrix(args.matrixA)
     matrix_B=DifFracTion_utils.load_dense_matrix(args.matrixB)

     

     if args.type_norm == 'alpha':
        tracemalloc.start()
        start_time = time.time()

        # Is always True on filter 
        norm_matrix_A, norm_matrix_B = DifFracTion.alpha_normalization(matrix_A,matrix_B,
                                                                         args.resolution,metric='median',
                                                                         desired_alpha=1.08,filter=True)
        normalization_time = time.time() - start_time
        _, norm_mem_peak = tracemalloc.get_traced_memory()
        mem_mb_norm = norm_mem_peak / 1024 / 1024

        alpha_padjusted_interactions,n_bins = DifFracTion.identify_differential_interactions(norm_matrix_A, norm_matrix_B,
                                                                      args.resolution, filter_by='padjusted', pvalue_cutoff=args.pval,
                                                                      log2fc_cutoff=0, neighbor_support=True, adjusted_pvalues_method=args.adjusted_pvalues_method)
        fully_elapsed_sec = time.time() - start_time
        _, mem_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        mem_mb = mem_peak / 1024 / 1024

        perf_log = pd.DataFrame([{'tool': 'DifFracTion_alpha', 'fully_elapsed_sec': fully_elapsed_sec, 'fully_mem_mb': mem_mb,'normalization_time': normalization_time, 'norm_mem_mb': mem_mb_norm}])
        os.makedirs(args.output_dir, exist_ok=True) #/performance
        perf_log.to_csv(f"{args.output_dir}/runtime_DifFracTion_alpha.txt", sep='\t', index=False)
        #alpha_padjusted_interactions.to_csv(f"{args.output_dir}/alpha_padjusted_interactions.txt", sep='\t', index=False)
        spike_in_files = glob.glob(os.path.dirname(args.matrixA) + "/spikein_and_neighbors_coordinates_k*.txt")
        confusion_m = get_confusion_matrix(spike_in_files, alpha_padjusted_interactions, n_bins)
        confusion_m.to_csv(f"{args.output_dir}/confusion_matrix_DifFracTion_alpha.txt", sep='\t', index=False)

        # Only True in neighbor support
        alpha_padjusted_interactions = alpha_padjusted_interactions[alpha_padjusted_interactions['neighbor_support'] == True]
        confusion_m_neighbors = get_confusion_matrix(spike_in_files, alpha_padjusted_interactions, n_bins)
        confusion_m_neighbors.to_csv(f"{args.output_dir}/confusion_matrix_DifFracTion_alpha_neighbors.txt", sep='\t', index=False)
     
     elif args.type_norm == 'iterative':
        tracemalloc.start()
        start_time = time.time()

        norm_matrix_A_iterative, norm_matrix_B_iterative = DifFracTion.iterative_normalization(matrix_A, matrix_B,
                                                                                            args.resolution, metric='median',
                                                                                            weight_factor=-1.08,filter=True)
        
        normalization_time = time.time() - start_time
        _, norm_mem_peak = tracemalloc.get_traced_memory()
        mem_mb_norm = norm_mem_peak / 1024 / 1024

        iterative_padjusted_interactions,n_bins= DifFracTion.identify_differential_interactions(norm_matrix_A_iterative, norm_matrix_B_iterative,
                                                                      args.resolution, filter_by='padjusted', pvalue_cutoff=args.pval,
                                                                      log2fc_cutoff=0, neighbor_support=True)
        fully_elapsed_sec = time.time() - start_time
        _, mem_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        mem_mb = mem_peak / 1024 / 1024

        perf_log = pd.DataFrame([{'tool': 'DifFracTion_iterative', 'fully_elapsed_sec': fully_elapsed_sec, 'fully_mem_mb': mem_mb, 'normalization_time': normalization_time, 'norm_mem_mb': mem_mb_norm}])
        os.makedirs(args.output_dir, exist_ok=True)
        perf_log.to_csv(f"{args.output_dir}/runtime_DifFracTion_iterative.txt", sep='\t', index=False)
        #iterative_padjusted_interactions.to_csv(f"{args.output_dir}/iterative_padjusted_interactions.txt", sep='\t', index=False)
        spike_in_files = glob.glob(os.path.dirname(args.matrixA) + "/spikein_and_neighbors_coordinates_k*.txt")
        confusion_m = get_confusion_matrix(spike_in_files, iterative_padjusted_interactions, n_bins)
        confusion_m.to_csv(f"{args.output_dir}/confusion_matrix_DifFracTion_iterative.txt", sep='\t', index=False)

        # Only True in neighbor support
        iterative_padjusted_interactions = iterative_padjusted_interactions[iterative_padjusted_interactions['neighbor_support'] == True]
        confusion_m_neighbors = get_confusion_matrix(spike_in_files, iterative_padjusted_interactions, n_bins)
        confusion_m_neighbors.to_csv(f"{args.output_dir}/confusion_matrix_DifFracTion_iterative_neighbors.txt", sep='\t', index=False)
     else:
        print(f"Normalization type {args.type_norm} not recognized. Please choose 'alpha' or 'iterative'.")
    #Differential analysis                                                                              