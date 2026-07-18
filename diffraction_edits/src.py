import hicstraw
import scipy
from scipy.sparse import coo_matrix
from scipy import ndimage
from scipy.stats import linregress
from statsmodels.stats.multitest import multipletests

import numpy as np
import pandas as pd
import os
from typing import Any,Literal

import matplotlib.pyplot as plt
from matplotlib import colors

from . import utils as DifFracTion_utils

#### DifFracTion functions ####

def alpha_normalization(
	matrix_A: np.ndarray, 
	matrix_B: np.ndarray, 
	resolution: int, 
	metric: Literal['median', 'mean'] = 'median',
	desired_alpha: float = 1.08, filter: bool = True
	) -> tuple[np.ndarray, np.ndarray]:
	"""
	Performs alpha normalization on two Hi-C contact matrices.
	
	Parameters:
	:param matrix_A (numpy.ndarray): First Hi-C contact matrix.
	:param matrix_B (numpy.ndarray): Second Hi-C contact matrix.
	:param resolution (int): Resolution of the Hi-C data in base pairs.
	:param metric (str): Metric to use for calculating the scaling factor ('median' or 'mean').
	:param desired_alpha (float): The desired alpha value for normalization.
	
	Returns:
	norm_matrix_A (numpy.ndarray): Alpha-normalized first matrix.
	norm_matrix_B (numpy.ndarray): Alpha-normalized second matrix.
	"""

	# These have all the bins within the distance decay limits, the three objects are the same size
	# Consider that the alpha values are contingent on the pseudocount added to the real counts, we added 0.01 to minimize the error in estimation of the decay
	# Adding 1 changed the alpha considerably
	#distances,counts_A,counts_B=get_distances_and_counts(matrix_A, matrix_B, resolution)


	# We remove overdispersed counts from our normalization as they can bias the estimation of the decay;
	# Overdispersion is defined as the mean counts of bins outside the distance decay limits. This comes from the fact that as genomic distance increases, counts become more unreliabre/hard to capture and therefore we cannot trust the values. 
	# All values on any distance that are smaller than the mean of the counts outside the distance decay limits are considered overdispersed and are filtered out.
	# Each matrix is filtered independently
	
	#  Because the filtering occurs inside this function, we just have the distances and counts of the bins that are not overdispersed, which is what we want for the normalization
	distances_A, counts_A, counts_B = DifFracTion_utils.get_distances_and_counts_paired(matrix_A, matrix_B, resolution, filter=filter)
	distances_B = distances_A  # Since we're dealing with paired matrices, the distances should be the same

	# These are the values of: 
	# - [0] unique distances that lay between the limits of the distance decay,
	# - [1] the observed metric values (e.g., median) for each distance [one per distance],
	# - [2] the expected values from the power-law fit for each distance,
	# - [3] the alpha (slope) of the decay,
	# - [4] the intercept of the decay,
	# - [5] the log of the distances (for log-log plots),
	# - [6] the log of the observed metric values (for log-log plots)
	#A_values,B_values=calculate_distance_decay(distances,counts_A,counts_B,metric='median')
	# Here it does not matter because all distances were already filtered
	# Here we get a print of the current alpha
	A_values=DifFracTion_utils.calculate_distance_decay_individual(distances_A, counts_A, metric=metric,verbose=False)
	B_values=DifFracTion_utils.calculate_distance_decay_individual(distances_B, counts_B, metric=metric,verbose=False)
	#Arrays here have elements per distance, not per bin
	
	new_intercept_A = DifFracTion_utils.refit_intercept(A_values[6], A_values[5], desired_alpha, form='log')
	new_intercept_B = DifFracTion_utils.refit_intercept(B_values[6], B_values[5], desired_alpha, form='log')

	b_shared = DifFracTion_utils.shared_intercept(A_values, B_values, desired_alpha)

	fit_A_new = DifFracTion_utils.refit_decay(b_shared,desired_alpha,A_values[5],form='log')
	fit_B_new = DifFracTion_utils.refit_decay(b_shared,desired_alpha,B_values[5],form='log')

	#Now we have the factors by which we need to multiply the original counts to get the fitted counts with the desired alpha
	#However these calculations are based on the range of distances of the distance decay, we need to apply a factor to the distances that do not lay on the range
	# A factor per distance inside the range
	factors_A = DifFracTion_utils.get_factors(A_values[1], fit_A_new)
	factors_B = DifFracTion_utils.get_factors(B_values[1], fit_B_new)

	norm_matrix_A = DifFracTion_utils.apply_decay_factors(matrix_A,factors_A,A_values[0],resolution,outside_range='interp')
	norm_matrix_B = DifFracTion_utils.apply_decay_factors(matrix_B,factors_B,B_values[0],resolution,outside_range='interp')

	return norm_matrix_A,norm_matrix_B

def iterative_normalization(
	matrix_A: np.ndarray, 
	matrix_B: np.ndarray, 
	resolution: int, 
	metric: Literal['median', 'mean'] = 'median', 
	weight_factor: float = -1.08, 
	eta: float = 0.15,filter: bool = True
	) -> tuple[np.ndarray, np.ndarray]:
	"""
	Performs iterative normalization on two Hi-C contact matrices.
	
	Parameters:
	matrix_A (numpy.ndarray): First Hi-C contact matrix.
	matrix_B (numpy.ndarray): Second Hi-C contact matrix.
	resolution (int): Resolution of the Hi-C data in base pairs.
	metric (str): Metric to use for calculating the scaling factor ('median' or 'mean').
	weight_factor (float): Weight factor for the iterative normalization. Default is -1.08.
	eta (float): Convergence threshold for the iterative normalization. Default is 0.15.
	filter (bool): Whether to filter out overdispersed counts based on CI for normalization. Default is True.
	Returns:
	norm_matrix_A (numpy.ndarray): Iteratively normalized first matrix.
	norm_matrix_B (numpy.ndarray): Iteratively normalized second matrix.
	"""
	norm_matrix_A_iterative, norm_matrix_B_iterative = DifFracTion_utils.iterative_correction(matrix_A, matrix_B, 
															resolution, metric=metric, 
															weight_factor=weight_factor, 
															eta=eta, tolerance=1e-6, 
															filter=True, plot=False)
	return norm_matrix_A_iterative, norm_matrix_B_iterative	

def identify_differential_interactions(
    norm_matrix_A: np.ndarray,
    norm_matrix_B: np.ndarray,
    resolution: int,
    filter_by: Literal['padjusted', 'pvalue'] = 'padjusted',
    pvalue_cutoff: float = 0.05,
    log2fc_cutoff: float = 1,
    neighbor_support: bool = True,
    bayesian: bool = True,
    permutations: int = 1000,
    upper_limit: float = 7e6,
    lower_limit: float = 5e4,
    adjusted_pvalues_method:Literal['global','distance'] = 'distance'
    ) -> pd.DataFrame | Any | None:
	"""Identifies differential interactions between two normalized Hi-C contact matrices.
	Arguments:
	norm_matrix_A (numpy.ndarray): First normalized Hi-C contact matrix.
	norm_matrix_B (numpy.ndarray): Second normalized Hi-C contact matrix.
	resolution (int): Resolution of the Hi-C data in base pairs.
	filter_by (str): Metric to filter interactions by ('pvalue' or 'padjusted').
	pvalue_cutoff (float): P-value cutoff for significance.
	log2fc_cutoff (float): Log2 fold change cutoff for significance.
	neighbor_support (bool): Whether to require significant interactions to have support from neighboring interactions.
	For neighbor support all interactions are reported, the flag is only used for visualization of significant interactions on the MA plot, where only those with neighbor support are highlighted.
	Returns:
     pandas.DataFrame: DataFrame containing significant interactions and their associated statistics.
	"""
	interactions,n_bins=DifFracTion_utils.differential_interactions(norm_matrix_A, norm_matrix_B, resolution, permutations_n=permutations,
                                                    filter_by=filter_by, pvalue_cutoff=pvalue_cutoff, log2fc_cutoff=log2fc_cutoff, 
										  neighbor_support=neighbor_support, bayesian=bayesian, upper_limit=upper_limit, lower_limit=lower_limit,
										  adjusted_pvalues_method=adjusted_pvalues_method,do_cutoff=True)
	
	return interactions,n_bins