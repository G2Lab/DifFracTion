import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from . import utils as DifFracTion_utils


def apply_spikeins(raw_matrix_sub, spikein_coordinates, neighbor_degrees, k):
	non_modified_raw_matrix=raw_matrix_sub.copy()
	altered_matrix = non_modified_raw_matrix.copy()

	all_neighbors = {deg: [] for deg in range(neighbor_degrees)} 
	spikein_all_lists = []

	for spikein in spikein_coordinates:
		#print(f"Processing spikein: {spikein}")	
		signal_neigh=analyze_neighbors(non_modified_raw_matrix,spikes=[spikein],K=neighbor_degrees)
		coords_neigh=get_neighbors_coordinates(non_modified_raw_matrix,spikes=[spikein],K=neighbor_degrees)

		altered_matrix=spike_with_degree_kernel(altered_matrix, spikein[0], spikein[1], 
												deg1_coords=coords_neigh[0], deg2_coords=coords_neigh[1], log2_fc=k)
		# it actually gets the  value
		spikes_list=get_spikes_coordinates(non_modified_raw_matrix,spikein, verbose=False)

		#We have to return the newest value, our current list contain the old value
		for deg in range(neighbor_degrees):
				all_neighbors[deg].extend(coords_neigh[deg])
		spikein_all_lists.extend(spikes_list)

	return altered_matrix, spikein_all_lists, all_neighbors

def select_spikeins(matrix1,resolution,n_spikes=1000,fold_change=2,upper_limit=7e6):
	'''Select spike ins in a Hi-C matrix'''
	lower_limit=resolution*5
	matrix_1 = matrix1.copy()
	n_bins = matrix_1.shape[0]
 
	spikes1 = 0

	# Save the coordinates of the spike ins
	spikes1_coords = []
  
	while spikes1 < n_spikes:
		x = np.random.randint(0,n_bins)
		y = np.random.randint(0,n_bins)
		
		if x==y or x > y : # restrict to upper triangle
			continue
            
		d = abs(x-y)*resolution

		if d < lower_limit or d > upper_limit:
			continue

		base_counts = matrix_1[x,y]
		if base_counts < 20:
			continue
		base_counts = matrix_1[x,y] 
		counts_fold_1 = (base_counts * (2 ** fold_change))
		matrix_1[x,y] = counts_fold_1
		matrix_1[y,x] = counts_fold_1
  
		spikes1 += 1
		spikes1_coords.append((x,y))
	print(f"Added {spikes1} spikes to matrix 1")

	return matrix_1, spikes1_coords

def spike_with_degree_kernel(M, i, j, deg1_coords, deg2_coords, log2_fc=4, sigma=1.0):
    '''Apply a spike-in at position (i,j) and propagate it to neighbors with a Gaussian kernel over degree'''
    M = M.copy()
    n = M.shape[0]

    c_old = M[i, j]
    c_new = c_old * (2 ** log2_fc)
    delta = c_new - c_old

    # set center
    M[i, j] = c_new
    M[j, i] = c_new  # keep symmetry

    # kernel weights by distance (degree)
    K = lambda d: np.exp(-(d**2) / (2 * sigma**2))  # Gaussian over degree
    w1 = K(1)
    w2 = K(2)

    for (u, v,_) in deg1_coords:
        if 0 <= u < n and 0 <= v < n:
            M[u, v] += delta * w1
            M[v, u] = M[u, v]

    for (u, v,_) in deg2_coords:
        if 0 <= u < n and 0 <= v < n:
            M[u, v] += delta * w2
            M[v, u] = M[u, v]

    return M

def get_neighbors(M, i, j):
    
    """Return direct (4-connected) neighbors of a matrix element.
	d is diagonal of the spike in
	row, column,value"""
    neighbors = []
    n_rows, n_cols = M.shape

    if i + 1 < n_rows:
        neighbors.append((i + 1, j, M[i + 1][j]))
    if i - 1 >= 0:
        neighbors.append((i - 1, j, M[i - 1][j]))
    if j + 1 < n_cols:
        neighbors.append((i, j + 1, M[i][j + 1]))
    if j - 1 >= 0:
        neighbors.append((i, j - 1, M[i][j - 1]))

    return neighbors

def analyze_neighbors(M, spikes, K=5, verbose=False):
    """
    General k-degree neighbor extraction.
    Returns: list of lists (values_per_degree)
        values_per_degree[d] = list of values for degree d+1
    """
    seen = set()
    values_per_degree = [[] for _ in range(K)]

    # Initialize frontier with all first-degree neighbors for each spike
    frontier = []  # list of (i,j,val)

    # Collect all first-degree neighbors from all spikes
    for (i0, j0) in spikes:
        seen.add((i0, j0))
        first_neighbors = get_neighbors(M, i0, j0)
        for i, j, val in first_neighbors:
            if (i, j) not in seen:
                frontier.append((i, j, val))
                seen.add((i, j))
        # Save first-degree values
        values_per_degree[0].extend([val for _, _, val in first_neighbors])

    if verbose:
        print(f"Degree 1 neighbors: {len(values_per_degree[0])}")

    # Iterate through higher degrees
    for degree in range(1, K):
        next_frontier = []
        for (i, j, val) in frontier:
            new_neighbors = get_neighbors(M, i, j)
            for ni, nj, nval in new_neighbors:
                if (ni, nj) not in seen:
                    next_frontier.append((ni, nj, nval))
                    seen.add((ni, nj))
                    values_per_degree[degree].append(nval)

        if verbose:
            print(f"Degree {degree+1} neighbors: {len(values_per_degree[degree])}")

        frontier = next_frontier  # Move frontier outward

        # If no new neighbors at this degree, stop early
        if len(frontier) == 0:
            break

    return (*values_per_degree,) # unpack the list of lists into separate lists

def get_neighbors_coordinates(M, spikes, K=5, verbose=False):
    """
    General K-degree neighbor expansion.
    Returns:
        A list of lists:
            neighbors_per_degree[d] = list of (i, j, val) for degree d+1
    """
    seen = set()
    neighbors_per_degree = [[] for _ in range(K)]

    # First degree neighbors 
    frontier = []  # list of (i, j, val)

    for (i0, j0) in spikes:
        seen.add((i0, j0))
        first_neighbors = get_neighbors(M, i0, j0)

        for i, j, val in first_neighbors:
            if (i, j) not in seen:
                frontier.append((i, j, val))
                neighbors_per_degree[0].append((i, j, val))
                seen.add((i, j))

    if verbose:
        print(f"Degree 1 neighbors: {len(neighbors_per_degree[0])}")

    # Higher degree neighbors
    for degree in range(1, K):
        next_frontier = []
        for (i, j, val) in frontier:
            new_neighbors = get_neighbors(M, i, j)
            for ni, nj, nval in new_neighbors:
                if (ni, nj) not in seen:
                    neighbors_per_degree[degree].append((ni, nj, nval))
                    next_frontier.append((ni, nj, nval))
                    seen.add((ni, nj))

        if verbose:
            print(f"Degree {degree+1} neighbors: {len(neighbors_per_degree[degree])}")

        frontier = next_frontier

        # No more expansion possible
        if len(frontier) == 0:
            break

    return (*neighbors_per_degree,)

def get_spikes_coordinates(M, spike,verbose=False):
	'''Return a list of lists containing the spike values from the matrix M at the given spike coordinates.'''
	
	
	individual_spike_values = [] #Add a list for each spike in so we can keep track of the respective values and neighbors
	individual_spike_values.append((spike[0], spike[1], M[spike[0], spike[1]]))
	
	if verbose:
		print("Spike values:", individual_spike_values)
	return individual_spike_values

def save_spike_and_neighbors_coordinates(spikes, all_neighbors, output_file):
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