from . import utils as DifFracTion_utils
from . import spikein as DifFracTion_spikes
import numpy as np
import pandas as pd
from functools import reduce
import os


# This is the main function to generate spike-ins and save coordinates for benchmarking
# [x] This function is used in benchmarking/generate_inputs.py to generate spike-ins and save coordinates for benchmarking
# [x] Done
main_path = DifFracTion_utils.main_path
##

def matrix2longdf(matrix,resolution=None):
        '''Args:
        matrix: 2D numpy array resultant from getRecordsAsMatrix 
        Returns:
        DataFrame with five columns: start1, end1, start2, end2, count [0,1,2,3,4]'''
        if resolution is None:
            resolution = 1

        matrix = np.triu(matrix,k=1)
        i,j = np.nonzero(matrix)# Two arrays of the same lenght
        return pd.DataFrame({
            'start1': i * resolution,
            'end1': (i + 1) * resolution,
            'start2': j * resolution,
            'end2': (j + 1) * resolution,
            'count': matrix[i, j]
        })


#---- For downsampling benchmarking

# One replicate per condition
def generate_HiCcompare_input_ds(
    raw_matrix,
    chromosome: str,
    resolution: int,
    downsample_factor: float,
    output_path: str
):
    def make_replicate(depth):
        return matrix2longdf(DifFracTion_utils.synthetic_datasets(raw_matrix, depth), resolution)

    coord_cols = ['start1', 'end1', 'start2', 'end2']

    # Condition A: full depth | Condition B: downsampled
    IF1 = make_replicate(0.9999999)
    IF2 = make_replicate(downsample_factor)

    merged = IF1.merge(IF2, on=coord_cols, how='inner').fillna(0)
    merged.columns = coord_cols + ['IF1', 'IF2']

    merged['chr1'] = chromosome
    merged['chr2'] = chromosome
    merged['D']    = abs(merged['start1'] - merged['start2']) // resolution
    merged['M']    = np.log2((merged['IF2'] + 0.000001) / (merged['IF1'] + 0.000001))

    merged = merged[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'IF1', 'IF2', 'D', 'M']]

    out_path = f"{output_path}/HiCcompare_input_chr{chromosome}_res{resolution}_ds{downsample_factor}.table"
    merged.to_csv(out_path, sep='\t', index=False, header=True)

    return merged, out_path

## ------- More than 1 replicate
def generate_HiCDCPlus_input_ds(
    raw_matrix,
    chromosome: str,
    resolution: int,
    downsample_factor: float,
    output_path: str
):
    def make_frame(matrix, resolution):
        matrix = np.triu(matrix, k=1)
        i, j = np.nonzero(matrix)
        return pd.DataFrame({
            'chr':    f'chr{chromosome}',
            'startI': i * resolution,
            'startJ': j * resolution,
            'counts': matrix[i, j]
        })

    def make_replicate(depth,resolution, condition):
        frame = make_frame(DifFracTion_utils.synthetic_datasets(raw_matrix, depth), resolution)
        out_path = f"{output_path}/HiCDCPlus_input_chr{chromosome}_res{resolution}_ds{downsample_factor}_{condition}.table"
        frame.to_csv(out_path, sep='\t', index=False, header=True)
        return out_path

    # Condition A: full depth | Condition B: downsampled
    save_path_A1 = make_replicate(0.9999999, resolution, "A1")
    save_path_A2 = make_replicate(0.9999999, resolution, "A2")
    save_path_B1 = make_replicate(downsample_factor, resolution, "B1")
    save_path_B2 = make_replicate(downsample_factor, resolution, "B2")

    return save_path_A1, save_path_A2, save_path_B1, save_path_B2

def generate_diffHic_input_ds(
	raw_matrix,
	chromosome: str,
	resolution: int,
	downsample_factor: float,
	save_path: str
):
	def make_replicate(matrix,resolution, col_name):
		return matrix2longdf(matrix,resolution).rename(columns={'count': col_name})

	# Condition A: two full-depth replicates
	A1 = make_replicate(DifFracTion_utils.synthetic_datasets(raw_matrix, 0.999999999999), resolution, 'sampleA_rep1')
	A2 = make_replicate(DifFracTion_utils.synthetic_datasets(raw_matrix, 0.999999999999), resolution, 'sampleA_rep2')

	# Condition B: two downsampled replicates
	B1 = make_replicate(DifFracTion_utils.synthetic_datasets(raw_matrix, downsample_factor), resolution, 'sampleB_rep1')
	B2 = make_replicate(DifFracTion_utils.synthetic_datasets(raw_matrix, downsample_factor), resolution, 'sampleB_rep2')
	# Merge all replicates on genomic coordinates
	coord_cols = ['start1', 'end1', 'start2', 'end2']
     
	def merge_replicates(left, right):
		return pd.merge(left, right, on=coord_cols, how='inner')

	merged = reduce(merge_replicates, [A1, A2, B1, B2]).fillna(0)

	merged['chr1'] = chromosome
	merged['chr2'] = chromosome
	merged = merged[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2',
	                  'sampleA_rep1', 'sampleA_rep2', 'sampleB_rep1', 'sampleB_rep2']]

	out_path = f'{save_path}/diffHic_input_chr{chromosome}_res{resolution}_ds{downsample_factor}.table'
	merged.to_csv(out_path, sep='\t', index=False, header=True)
	return merged, out_path

def generate_multiHiCcompare_input_ds(
    raw_matrix,
    chromosome: str,
    resolution: int,
    downsample_factor: float,
    output_path: str
):
    

    def generate_individual_frame(matrix, resolution, col_name):
        matrix_df = matrix2longdf(matrix, resolution)
        matrix_df.columns = ['start1', 'end1', 'start2', 'end2', col_name]
        matrix_df['chr1'] = chromosome
        matrix_df['chr2'] = chromosome
        return matrix_df[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', col_name]]

    
    def make_replicate(depth,condition):
        frame = generate_individual_frame(DifFracTion_utils.synthetic_datasets(raw_matrix, depth), resolution, condition)
        save_path = f"{output_path}/multiHiCcompare_input_chr{chromosome}_res{resolution}_ds{downsample_factor}_{condition}.table"
        frame.to_csv(save_path, sep='\t', index=False, header=True)
        return save_path
    
    # Condition A: Full Depth
    save_path_A1 = make_replicate(0.9999999, "IF_A1")
    save_path_A2 = make_replicate(0.9999999, "IF_A2")

    # Condition B: Downsampled
    save_path_B1 = make_replicate(downsample_factor, "IF_B1")
    save_path_B2 = make_replicate(downsample_factor, "IF_B2")


    return save_path_A1, save_path_A2, save_path_B1, save_path_B2


#---- For spike-in benchmarking
## --- Main 
def generate_and_save_spike_ins(
    hic_path: str,
    chromosome: str,
    resolution: int,
    fold_change: int,
    neighbor_degrees: int,
    matrix_save_path: str,
    coordinate_file_to_save: str
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int, float]]]:
    #matrix_save_path=f"{main_path}/generated_data/"
    # This function is only used in benchmarking against other tools
    # find its use on DifFracTion/benchmarking/generate_inputs.py
	n_spikes = 100
	resolution_kb = resolution // 1000
	lower_limit = 5e4
	#has to be on mainpath/generated_data
	cell_type=os.path.basename(hic_path).split("-")[0]
	name_matrix=f"{matrix_save_path}/{cell_type}_chr{chromosome}_{resolution_kb}kb.npz"

	#if not os.path.exists(f"{name_matrix}"):
	DifFracTion_utils.save_sparse_matrix(hic_path,chromosome,resolution,f"{name_matrix}")
	matrix_A = DifFracTion_utils.load_dense_matrix(f"{name_matrix}")


	altered_matrix_A,spike_coordinates,neighbors= DifFracTion_spikes.run_spike_ins(
		matrix_A,resolution,n_spikes=n_spikes,fold_change=fold_change,neighbor_degrees=neighbor_degrees,
		upper_limit=None, lower_limit=lower_limit)


	DifFracTion_spikes.save_spike_and_neighbors_coordinates(
			spike_coordinates,
			neighbors,
			coordinate_file_to_save
		)

	return altered_matrix_A, matrix_A, spike_coordinates,name_matrix

## ===== More than 1 replicate per condition
def generate_HiCDCPlus_input(
    altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str
):  
    #Made based on requirements from HiCDCPlus
    def make_frame(matrix):
        matrix = np.triu(matrix, k=1)
        i, j = np.nonzero(matrix)
        return pd.DataFrame({
            'chr':    f'chr{chromosome}',
            'startI': i * resolution,
            'startJ': j * resolution,
            'counts': matrix[i, j]
        })
    
    
    
    def make_replicate(m, resolution, spikes_k, condition):
        frame=make_frame(DifFracTion_utils.synthetic_datasets(m, 0.9999999))
        save_path = (
            f"{output_path}/HiCDCPlus_input_chr{chromosome}_res{resolution}_k{spikes_k}_{condition}.table"
        )
        frame.to_csv(save_path, sep='\t', index=False, header=True)
        return save_path
    

    # Condition A: full depth with spike-ins
    save_path_A1 = make_replicate(altered_matrix, resolution, spikes_k, "A1")
    save_path_A2 = make_replicate(altered_matrix, resolution, spikes_k, "A2")

    # Condition B: raw matrix with no spike-ins but same depth
    save_path_B1 = make_replicate(raw_matrix, resolution, spikes_k, "B1")
    save_path_B2 = make_replicate(raw_matrix, resolution, spikes_k, "B2")

    return save_path_A1, save_path_A2, save_path_B1, save_path_B2
    
def generate_diffHic_input(
    altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str):

    def make_replicate(matrix, col_name):
        return matrix2longdf(matrix, resolution).rename(columns={'count': col_name})
    # Condition A: full depth with spike-ins
    A1 = make_replicate(DifFracTion_utils.synthetic_datasets(altered_matrix, 0.999999999999), "sampleA_rep1")
    A2 = make_replicate(DifFracTion_utils.synthetic_datasets(altered_matrix, 0.999999999999), "sampleA_rep2")
    # Condition B: raw matrix with no spike-ins but same depth
    B1 = make_replicate(DifFracTion_utils.synthetic_datasets(raw_matrix, 0.999999999999), "sampleB_rep1")
    B2 = make_replicate(DifFracTion_utils.synthetic_datasets(raw_matrix, 0.999999999999), "sampleB_rep2")

    coord_cols = ['start1', 'end1', 'start2', 'end2']

    # Merge  all replicates

    def merge_replicates(left, right):
        return pd.merge(left, right, on=coord_cols, how='inner')

    merged = reduce(merge_replicates, [A1, A2, B1, B2]).fillna(0)
    merged['chr1'] = chromosome
    merged['chr2'] = chromosome
    merged = merged[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2',
                      'sampleA_rep1', 'sampleA_rep2', 'sampleB_rep1', 'sampleB_rep2']]
    

    save_path = (
        f"{output_path}/diffHic_input_chr{chromosome}_res{resolution}_k{spikes_k}.table"
    )

    merged.to_csv(save_path, sep='\t', index=False, header=True)
    return merged, save_path

def generate_multiHiCcompare_input(
    altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str
):
    

    def generate_individual_frame(matrix, resolution, col_name):
        matrix_df = matrix2longdf(matrix, resolution)  # upper triangle already
        
        matrix_df.columns = ['start1','end1','start2','end2',col_name]
        matrix_df['chr1'] = chromosome
        matrix_df['chr2'] = chromosome
        return  matrix_df[['chr1','start1','end1','chr2','start2','end2',col_name]]
    
    def make_replicate(m, resolution, spikes_k, condition_rep):
        matrix = DifFracTion_utils.synthetic_datasets(m, 0.9999999) # this is because we already have spike-ins
        frame = generate_individual_frame(matrix, resolution, f"{condition_rep}")
        save_path = (
            f"{output_path}/multiHiCcompare_input_chr{chromosome}_res{resolution}_k{spikes_k}_{condition_rep}.table"
        )
        frame.to_csv(save_path, sep='\t', index=False, header=True)
        return save_path
    
    # Condition A: full depth with spike-ins
    save_path_A1 = make_replicate(altered_matrix, resolution, spikes_k, "IF_A1")
    save_path_A2 = make_replicate(altered_matrix, resolution, spikes_k, "IF_A2")
    # Condition B: raw matrix with no spike-ins but same depth
    save_path_B1 = make_replicate(raw_matrix, resolution, spikes_k, "IF_B1")
    save_path_B2 = make_replicate(raw_matrix, resolution, spikes_k, "IF_B2")


    return save_path_A1, save_path_A2, save_path_B1, save_path_B2

# ===== One replicate per condition
def generate_HiCcompare_input(altered_matrix,
    raw_matrix,
    chromosome: str,
    resolution: int,
    spikes_k: float,
    output_path: str
):
    """
    altered_matrix : matrix with spike-ins applied
    raw_matrix     : original real matrix
    DifFracTion_utils  : module namespace (DifFracTion_utils) so function works in scripts too
    """
    def make_replicate(depth, resolution):
        return matrix2longdf(DifFracTion_utils.synthetic_datasets(raw_matrix, depth), resolution)

    
    coord_cols = ['start1', 'end1', 'start2', 'end2']

    #Condition A: full depth and altered with spike-ins | Condition B: raw matrix with no spike-ins but same depth
    IF1=matrix2longdf(altered_matrix, resolution)  # upper triangle already
    #previously m1
    IF2=make_replicate(0.999999, resolution)  # full depth 

    merged = IF1.merge(IF2, on=coord_cols, how='inner').fillna(0)
    merged.columns = coord_cols + ['IF1', 'IF2']

    merged['chr1'] = chromosome
    merged['chr2'] = chromosome

    merged['D'] = abs(merged['start1'] - merged['start2']) // resolution
    merged['M'] = np.log2((merged['IF2']) / (merged['IF1']))

    matrices_merged = merged[
        ['chr1','start1','end1','chr2','start2','end2','IF1','IF2','D','M']
    ]

    save_path = (
        f"{output_path}/HiCcompare_input_chr{chromosome}_res{resolution}_k{spikes_k}.table"
    )

    matrices_merged.to_csv(save_path, sep='\t', index=False, header=True)

    return matrices_merged, save_path

