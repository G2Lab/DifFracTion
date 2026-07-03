import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


# -------- Module: Spike-in generation and analysis --------
def select_spikeins(matrix,resolution,
					n_spikes=1000,fold_change=2,
					upper_limit=None, lower_limit=5e4):
	'''Args:
		matrix: dense matrix to add spike ins to
		resolution: resolution of the matrix in bp
		n_spikes: number of spike ins to add
		fold_change: fold change to add to the spike ins
		upper_limit: upper distance limit for spike ins in bp
		lower_limit: lower distance limit for spike ins in bp
	Returns:
		- matrix_A: matrix with spike ins added
		- spikes_coordinates: list of coordinates where spike ins were added'''
	
	matrix_A = matrix.copy()
	n_bins = matrix_A.shape[0]
 
	spikes_record= 0

	# Save the coordinates of the spike ins
	spikes_coordinates = []
	
	# While we haven't added the desired number of spike ins, keep trying random coordinates
	while spikes_record < n_spikes:
		# Randomly select two bins (x,y) in the matrix
		row = np.random.randint(0,n_bins) # row index
		col = np.random.randint(0,n_bins) # column index
		
		if row==col or row > col : # restrict to upper triangle
			continue
          # here
		distance = abs(row-col)*resolution

		if upper_limit is None:
			upper_limit = n_bins * resolution # maximum distance in the matrix
		
		if distance < lower_limit or distance > upper_limit:
			continue

		base_counts = matrix_A[row,col]

		if base_counts < 20:
			continue
		 
		# This is the form 
		counts_to_fold = (base_counts * (2 ** fold_change))
		matrix_A[row,col] = counts_to_fold
		matrix_A[col,row] = counts_to_fold
  
		spikes_record+= 1
		spikes_coordinates.append((row,col))
	print(f"Added {spikes_record} spikes to matrix")

	return matrix_A, spikes_coordinates

def k_neighbor_extraction(matrix,spikes_coordinates,k=2):
	#previously named analyze_neighbors and get_neighbors_coordinates 
	'''Args:
		matrix: Unmodified matrix 
		spikes_coordinates: (tuple of coordinates) where spike ins were added
		k: number of neighbors to extract
	Returns:A list per degree of neighbors that can be accessed by indexing the output with the degree of interest'''

	#List of lists (k)
	#One list inside neighbors_per_degree[degree] 
	neighbors_per_degree = [[] for _ in range(k)]
	seen = set() # to keep track of seen coordinates and avoid duplicates
	
	first_degree_set = []
	#Core of each spike
	for (i,j) in spikes_coordinates:
		seen.add((i,j))
		neighbors = extract_bin_neighbors(matrix,i,j)
		for i_prime,j_prime, value in neighbors:
			if (i_prime,j_prime) not in seen:
				seen.add((i_prime,j_prime))
				first_degree_set.append((i_prime,j_prime,value)) # Set with all first degree neighbors
				neighbors_per_degree[0].append((i_prime,j_prime,value)) #A list inside neighbors_per_degree[0] with only the values of the first degree neighbors, not the coordinates
	
	# now we have the 1st neighbors, repeat iteratively
	for degree in range(1,k):
		degree_frontier = []
		for (i_first,j_first,value_first) in first_degree_set:
			new_neighbors = extract_bin_neighbors(matrix,i_first,j_first)
			for i_first,j_first, value in new_neighbors:
				if (i_first,j_first) not in seen:
					seen.add((i_first,j_first)) # tuple
					degree_frontier.append((i_first,j_first,value)) # Set with all neighbors of the current degree
					neighbors_per_degree[degree].append((i_first,j_first,value))
		#A list inside neighbors_per_degree[degree] with only the values of the neighbors of the current degree, not the coordinates
		first_degree_set = degree_frontier # update the frontier for the next iteration

		if len(degree_frontier) == 0:
			break # if there are no more neighbors to explore, stop the loop
	
	#the * is to unpack the list of lists into separate lists for each degree, so that we can access them by indexing the output with the degree of interest
	return (*neighbors_per_degree,)

def extract_bin_neighbors(matrix,i,j):
	#Previously named get_neighbors 
	'''Arguments:
	matrix: original matrix to retrieve values from
	i,j: coordiantes of a bin'''

	neighbors = []
	rows, cols = matrix.shape

	if i + 1 < rows: #down
		neighbors.append((i+1,j,matrix[i+1,j]))
	if i - 1 >= 0: #up
		neighbors.append((i-1,j,matrix[i-1,j]))
	if j + 1 < cols: #right
		neighbors.append((i,j+1,matrix[i,j+1]))
	if j - 1 >= 0: #left
		neighbors.append((i,j-1,matrix[i,j-1]))

	return neighbors

def apply_kernel_based_on_degree(matrix,i,j,d0,d1,fold_change,sigma=1.0):
	'''Applies a kernel-based fold change to the neighbors of a spike in, with the fold change decaying based on the degree of the neighbor.'''
	'''Arguments:
	matrix: original matrix
	i,j: spike in coordinates
	d0: list of first degree neighbors (list of tuples with coordinates )
	d1: list of second degree neighbors (list of tuples with coordinates )
	FC: fold change to apply to the spike in'''

	matrix_copy = matrix.copy()
	rows, cols = matrix_copy.shape

	c_old = matrix_copy[i,j]
	c_new = c_old * (2 ** fold_change)
	D=c_new - c_old #delta

	matrix_copy[i,j] = c_new
	matrix_copy[j,i] = c_new # ensure symmetry

	#Gaussian kernel
	K = lambda d: np.exp(- (d**2) / (2 * sigma**2))
	w_d0 = K(1) # weight for first degree neighbors
	w_d1 = K(2) # weight for second degree neighbors

	# Iterate over the neighbors that are on the form [(i,j,value),...]
	for i_d0, j_d0, value_d0 in d0:
		#They have to be within the bounds of the matrix
		if 0 <= i_d0 < rows and 0 <= j_d0 < cols:
			matrix_copy[i_d0, j_d0] += D * w_d0
			matrix_copy[j_d0, i_d0] += D * w_d0 # ensure symmetry
	for i_d1, j_d1, value_d1 in d1:
		#They have to be within the bounds of the matrix
		if 0 <= i_d1 < rows and 0 <= j_d1 < cols:
			matrix_copy[i_d1, j_d1] += D * w_d1
			matrix_copy[j_d1, i_d1] += D * w_d1 # ensure symmetry
	
	return matrix_copy

# -------- Executors --------
def run_spike_ins(matrix_A,resolution,n_spikes=1000,fold_change=2,neighbor_degrees=2,upper_limit=7e6, lower_limit=5e4):
	# Previous name: def apply_spikeins(raw_matrix_sub, spikein_coordinates, neighbor_degrees, k):
	# spikes_coordinates is a list of tuples
	_,spikes_coordinates=select_spikeins(matrix_A,resolution,n_spikes=n_spikes,fold_change=fold_change,upper_limit=upper_limit, lower_limit=lower_limit)


	non_modified_raw_matrix = matrix_A.copy()
	altered_matrix = matrix_A.copy() # is not matrix_a_prime because on apply_kernel_based_on_degree we will apply the FC to the spike in and the neighbors

	
	all_neighbors = {deg: [] for deg in range(neighbor_degrees)} 
	all_spikes_list = [] # to keep track of the spike ins and their values for downstream analysis

	for spikein in spikes_coordinates:
		spike_i = spikein[0]
		spike_j = spikein[1] 

		# ----- This whole chunck is enough for the process
		#Gives a list of lists
		neighbors_per_degree = k_neighbor_extraction(non_modified_raw_matrix,[(spike_i,spike_j)],k=neighbor_degrees)
		
		#Updates the matrix per each spikein
		altered_matrix = apply_kernel_based_on_degree(altered_matrix,spike_i,spike_j,
								neighbors_per_degree[0],neighbors_per_degree[1],fold_change=fold_change,sigma=1.0)
		# ------------------------------------------

		individual_spike_values = []

		# we append a tuple with i,j, value to a list
		individual_spike_values.append((spike_i, spike_j, non_modified_raw_matrix[spike_i, spike_j]))
		#By extending it avoid dealing with nested lists, and all spikes will be in a tuple format in the all_spikes_list
		all_spikes_list.extend(individual_spike_values)
		
		
		# For each degree, take the neighbors and add it to the all neighbors list of lists
		for degree in range(neighbor_degrees):
			d_neighbors = neighbors_per_degree[degree]
			all_neighbors[degree].extend(d_neighbors) # we keep track of all the neighbors for all the spike ins, separated by degree

	altered_matrix = np.triu(altered_matrix)

	return altered_matrix, all_spikes_list, all_neighbors
	
def save_spike_and_neighbors_coordinates(
		spikes: list[tuple[int, int, float]], 
		all_neighbors: dict[int, list[tuple[int, int, float]]], 
		output_file: str):
	'''Save spike coordinates to a text file
	Args:
		spikes: list of tuples containing spike coordinates (i,j)
		all_neighbors: dictionary of neighbor lists per degree, with degree as key and list of (i,j,value) tuples as values
		output_file: path to the output text file'''
	
	with open(output_file, 'w') as f:
		for (i, j,_) in spikes:
			f.write(f"{i}\t{j}\tspike\n")
		for deg, neighbors in all_neighbors.items():
			for (i, j,_) in neighbors:
				f.write(f"{i}\t{j}\tneighbor_{deg+1}\n")

def plot_original_spiked(non_modified_raw_matrix, altered_matrix):	
    fig, axes = plt.subplots(
        1, 2,
        figsize=(9,5),
        constrained_layout=True,
        dpi=350
    )

    vmax = np.percentile(altered_matrix, 99)
    vmin = 0

    im0 = axes[0].imshow(non_modified_raw_matrix, vmin=vmin, vmax=vmax, cmap='RdBu_r')
    axes[0].set_title("Original Matrix")

    im1 = axes[1].imshow(altered_matrix, vmin=vmin, vmax=vmax, cmap='RdBu_r')
    axes[1].set_title("Spiked-In Matrix")

    fig.colorbar(im1, ax=axes, location='bottom', fraction=0.05, pad=0.08)
    return fig