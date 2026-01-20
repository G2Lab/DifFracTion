## Main function to normalize
import os
import hicstraw 
import cooler
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from statistics import mean, median
from scipy.stats import linregress
from scipy.sparse import coo_matrix
import math
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from . import utils as DifFracTion_utils

#### DifFracTion functions ####

def alpha_normalization(matrixA,matrixB,chromosome_length,resolution,desired_alpha=-1.08,upper_limit=7e6,lower_limit=0,metric='median'):
	'''
	Docstring for alpha_normalization
	
	:param matrixA: nxn numpy array representing the first contact matrix
	:param matrixB: nxn numpy array representing the second contact matrix
	:param chromosome_length: Length of the chromosome in base pairs
	:param resolution: Bin resolution in base pairs
	:param desired_alpha: Desired alpha value for normalization
	:param upper_limit: Upper distance limit for normalization, default is 7e6
	:param lower_limit: Lower distance limit for normalization, default is 0
	:param metric: Metric to use for distance decay calculation (e.g., 'median' or 'mean')
	'''

	matrixA_counts,matrixB_counts=DifFracTion_utils.get_counts_by_distance_both_densematrices(matrixA,matrixB,chromosome_length,resolution)
	d_1, y_1, m_1, b_1, alpha_1 = DifFracTion_utils.calculate_distance_decay(matrixA_counts, metric=metric,
                              lower_limit=lower_limit,upper_limit=upper_limit)
	d_2, y_2, m_2, b_2, alpha_2 = DifFracTion_utils.calculate_distance_decay(matrixB_counts, metric=metric,
                              lower_limit=lower_limit,upper_limit=upper_limit)
     																

	scaled_counts_1,scaled_counts_2,full_factors_r,full_factors_s=DifFracTion_utils.pre_counts(matrixA_counts,matrixB_counts,d_1,d_2,y_1,y_2,desired_alpha,resolution)
     
	#### Apply factors to dense matrices
	matrixA_adj=DifFracTion_utils.apply_dict_factors_to_dense(matrixA,full_factors_r,resolution)
	matrixB_adj=DifFracTion_utils.apply_dict_factors_to_dense(matrixB,full_factors_s,resolution)
	return matrixA_adj,matrixB_adj,scaled_counts_1,scaled_counts_2

def iterative_normalization(matrixA,matrixB,chromosome_length,resolution,upper_limit=7e6,lower_limit=0,metric='median', verbose=True):
     '''
     Docstring for iterative_normalization

     :param matrixA: nxn numpy array representing the first contact matrix
     :param matrixB: nxn numpy array representing the second contact matrix
     :param chromosome_length: Length of the chromosome in base pairs
     :param resolution: Bin resolution in base pairs
     :param upper_limit: Upper distance limit for normalization, default is 7e6
     :param lower_limit: Lower distance limit for normalization, default is 0
     :param metric: Metric to use for distance decay calculation (e.g., 'median' or 'mean')
     :param verbose: -
     '''
     matrixA_counts,matrixB_counts=DifFracTion_utils.get_counts_by_distance_both_densematrices(matrixA,matrixB,chromosome_length,resolution)
     matrixB_counts_scaled, best_error, updated_factor = DifFracTion_utils.cycle_correction(matrixA_counts, matrixB_counts,metric=metric,
											error_type='mse',upper_limit=upper_limit,lower_limit=lower_limit)
     matrixA_counts_scaled=matrixA_counts
     #### Apply factors to dense matrices
     matrixA_adj=DifFracTion_utils.apply_factors_to_dense(matrixA,1,resolution)
     matrixB_adj=DifFracTion_utils.apply_factors_to_dense(matrixB,updated_factor,resolution)
     if verbose:
          print(f"Final scaling factor applied to matrix B: {updated_factor:.4f}")
          print(f"Final error after iterative normalization: {best_error:.4f}")
     return matrixA_adj,matrixB_adj,matrixA_counts_scaled,matrixB_counts_scaled

def identify_differential_interactions(matrixA_adj,matrixB_adj,resolution,method='log2fc',log2_fc_cutoff=0,p_value_threshold=0.01,subsample_n=1000,bayesian=False):
     '''
     Docstring for call_differential_interactions

     :param matrixA_adj: nxn numpy array representing the first normalized contact matrix
     :param matrixB_adj: nxn numpy array representing the second normalized contact matrix
     :param method: Method to use for calling differential interactions ('log2fc' or 'difference')
     :param p_value_threshold: P-value threshold for significance
     :param bayesian: Whether to perform Bayesian fold change estimation
     :param log2_fc_cutoff: Log2 fold change cutoff for significance
     :param permutations_n: Number of permutations for the permutation test
     '''
     
     diff_interactions=DifFracTion_utils.perform_permutation_test(matrixA_adj,matrixB_adj,resolution,method=method,permutations_n=subsample_n)
     diff_interactions=DifFracTion_utils.get_p_values(diff_interactions,method=method)
     diff_interactions=DifFracTion_utils.adjust_p_values(diff_interactions)
     if bayesian:
          diff_interactions=DifFracTion_utils.bayesian_fc_estimation(diff_interactions)
          significant_contacts=DifFracTion_utils.significant_interactions(diff_interactions,p_value_threshold,log2_fc_cutoff)
          significant_contacts=DifFracTion_utils.neighbor_support(significant_contacts)
     else:
          significant_contacts=DifFracTion_utils.significant_interactions_no_bayesian(diff_interactions,p_value_threshold,log2_fc_cutoff) 
          significant_contacts=DifFracTion_utils.neighbor_support(significant_contacts)
     return significant_contacts 