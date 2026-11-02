import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
import pybedtools
import glob

from scipy.stats import fisher_exact


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from diffraction import utils as DifFracTion_utils


##### Core

def process_DifFracTion_results(results_path, chromosomes, resolution,log_fold,type="alpha"):
	main_frame = pd.DataFrame()
	for chro in chromosomes:
		file = glob.glob(results_path + f"chr{chro}_{resolution}_{type}_normalization_differential_interactions.csv")
		if file:
			df = pd.read_csv(file[0], sep='\t')
			df['chr'] = f'chr{chro}'
			main_frame = pd.concat([main_frame, df], ignore_index=True)
	main_frame = main_frame[main_frame['neighbor_support'] == True]
	main_frame = main_frame[main_frame['log2fc'].abs() >= log_fold]
	main_frame = main_frame[(main_frame['condition_A_counts'] >= 25) & (main_frame['condition_B_counts'] >= 25)]

	positive_interactions = main_frame[main_frame['log2fc'] > 0]
	negative_interactions = main_frame[main_frame['log2fc'] < 0]
	return  positive_interactions, negative_interactions


##### Data loading

def process_bed_pe(bedpe_file:str):
	bedpe = pd.read_csv(bedpe_file, sep='\t', header=None)
	bedpe = bedpe[[0, 1, 2, 3, 4, 5, 6]]  # Keep only the first seven columns (chrom1, start1, end1, chrom2, start2, end2, name)
	bedpe.columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name']  
	return bedpe

def process_bed(bed_file:str):
	bed = pd.read_csv(bed_file, sep='\t', header=None)
	bed = bed[[0, 1, 2]]  # Keep only the first three columns (chromosome, start, end)
	bed.columns = ['chr', 'start', 'end'] 
	return bed

##### Interaction processing

def process_interactions_both_anchors(interactions:pd.DataFrame,resolution:int):
	'''Process interactions to keep both anchors in the same row, with columns for chr1, start1, end1, chr2, start2, end2'''
	interactions = interactions[["row",'col','chr']]
	interactions = interactions.rename(columns={'row': 'start', 'col': 'end'})

	interactions = interactions[['chr', 'start', 'end']]
	interactions['start'] = interactions['start'] * resolution
	interactions['start_end'] = interactions['start'] + resolution
	interactions['end'] = interactions['end'] * resolution
	interactions['end_end'] = interactions['end'] + resolution
	interactions['chr2'] = interactions['chr']

	interactions = interactions[['chr', 'start', 'start_end', 'chr2', 'end', 'end_end']]  # Keep only the relevant columns
	interactions.columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
	interactions['ID'] = interactions['chr1'] + '_' + interactions['start1'].astype(str) + '_' + interactions['start2'].astype(str) #original
	return interactions


##### DEG analysis

def find_closest_gene(interactions, genes_coords):
    anchor1 = interactions[['chr1', 'start1', 'end1', 'ID']].copy()
    anchor2 = interactions[['chr2', 'start2', 'end2', 'ID']].copy()
    anchor1.columns = ['chr', 'start', 'end', 'ID']
    anchor2.columns = ['chr', 'start', 'end', 'ID']

    bed = genes_coords[['chr', 'start', 'end', 'GeneName']].copy()

    cl1 = pybedtools.BedTool.from_dataframe(anchor1).sort().closest(
              pybedtools.BedTool.from_dataframe(bed).sort(), d=True
          ).to_dataframe(names=['chr','start','end','ID','b_chr','b_start','b_end','gene_name_1','dist1'])

    cl2 = pybedtools.BedTool.from_dataframe(anchor2).sort().closest(
              pybedtools.BedTool.from_dataframe(bed).sort(), d=True
          ).to_dataframe(names=['chr','start','end','ID','b_chr','b_start','b_end','gene_name_2','dist2'])

    merged = pd.merge(
        cl1[['ID', 'chr', 'start', 'end', 'gene_name_1', 'dist1']],
        cl2[['ID', 'gene_name_2', 'dist2']],
        on='ID', how='outer'
    )
    merged['min_dist'] = merged[['dist1', 'dist2']].min(axis=1)
    merged = merged.loc[merged.groupby('ID')['min_dist'].idxmin()].reset_index(drop=True)
    merged = merged[['ID', 'gene_name_1', 'gene_name_2', 'dist1', 'dist2', 'min_dist']]
    return merged

##### Marks enrichment analysis

def get_overlap_both_anchors(interactions:pd.DataFrame, bed:pd.DataFrame):
    '''interactions: DataFrame with columns: chr1, start1, end1, chr2, start2, end2, ID
    bed: DataFrame with columns: chr, start, end
    Returns a DataFrame with columns: ID, b_chr, b_start, b_end for interactions that overlap the bed entries on both anchors
    '''

    anchor1 = interactions[['chr1', 'start1', 'end1', 'ID']].copy()
    anchor2 = interactions[['chr2', 'start2', 'end2', 'ID']].copy()
    anchor1.columns = ['chr', 'start', 'end', 'ID']
    anchor2.columns = ['chr', 'start', 'end', 'ID']

    bed = bed[['chr', 'start', 'end']].copy()
    
    #Overlap any of the anchors with the bedpe entries, keeping the ID of the bedpe entry that overlaps
    ov1 = pybedtools.BedTool.from_dataframe(anchor1).intersect(
              pybedtools.BedTool.from_dataframe(bed), wa=True
          ).to_dataframe(names=['chr','start','end','ID','b_chr','b_start','b_end'])

    ov2 = pybedtools.BedTool.from_dataframe(anchor2).intersect(
              pybedtools.BedTool.from_dataframe(bed), wa=True
          ).to_dataframe(names=['chr','start','end','ID','b_chr','b_start','b_end'])

    # Only keep hits where the same interaction
    merged = pd.merge(
        ov1[['ID', 'b_chr', 'b_start', 'b_end']],
        ov2[['ID', 'b_chr', 'b_start', 'b_end']],
        on=['ID']
    )
    #Return two separate objects, one for those interactions that overlap the bed entries on both anchors, and one for those that overlap on only one anchor
    both_interactions = merged[['ID']].drop_duplicates()
    single_interactions = pd.concat([ov1[['ID']], ov2[['ID']]]).drop_duplicates()
    single_interactions = single_interactions[~single_interactions['ID'].isin(both_interactions['ID'])]
    return both_interactions, single_interactions

def random_regions_both_anchors(interactions, resolution):
	'''
	Generate random interactions preserving chromosome and interaction distance.
	:param interactions: DataFrame from process_DifFracTion_results, with columns row, col, chr
	:param resolution: analysis resolution in base pairs
	:return: DataFrame with columns chr1, start1, end1, chr2, start2, end2, ID
	'''
	simulated_interactions = []
	chr_sizes = pd.read_csv(f"{DifFracTion_utils.main_path}/test_data/hg19_chrsizes.txt", sep='\t', header=None)

	for _, row in interactions.iterrows():
		distance_bin = row['col'] - row['row']
		chr_row = row['chr']
		chr_size = chr_sizes[chr_sizes[0] == chr_row][1].values[0]
		chr_size_bin = chr_size // resolution
		random_start = np.random.randint(1, chr_size_bin - distance_bin - 1)
		random_start_end = random_start + 1
		random_end = random_start + distance_bin
		random_end_end = random_end + 1
		simulated_interactions.append([chr_row, random_start * resolution, random_start_end * resolution,
		                               chr_row, random_end * resolution, random_end_end * resolution])

	simulated_df = pd.DataFrame(simulated_interactions, columns=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
	simulated_df['ID'] = (simulated_df['chr1'] + '_' + simulated_df['start1'].astype(str) + '_' +
	                      simulated_df['start2'].astype(str))
	return simulated_df

def fisher_enrichment(interactions, interactions_raw, bed, resolution, seed=42):
    """
    Fisher's exact test: observed DCI hits vs random background.
    Random regions preserve chromosome and interaction distance of each DCI.
    Rows: [observed DCIs, random DCIs]. Cols: [hits peaks, does not hit peaks].
    """
    random_interactions = random_regions_both_anchors(interactions_raw, resolution)

    obs_both, obs_single = get_overlap_both_anchors(interactions, bed)
    rnd_both, rnd_single = get_overlap_both_anchors(random_interactions, bed)

    x = len(obs_both) + len(obs_single)      # observed hits
    q = len(interactions) - x                # observed misses
    
    y = len(rnd_both) + len(rnd_single)      # random hits
    p = len(random_interactions) - y         # random misses

    odds_ratio, pvalue = fisher_exact([[x, q], [y, p]], alternative='greater')
    return odds_ratio, pvalue, x, q, y, p
##### rE2G analysis

def closest_bedpe_distance(interactions: pd.DataFrame, bedpe: pd.DataFrame):
    """
    For each interaction, find the closest paired bedpe entry and return
    the distance from each anchor to its matched bedpe anchor.
    Pairing is enforced: anchor1 and anchor2 must match the SAME bedpe entry.
    Returns a DataFrame with columns: ID, bedpe_ID, dist1, dist2, max_dist
    """
    bedpe = bedpe.copy()
    anchor1 = interactions[['chr1', 'start1', 'end1', 'ID']].copy()
    anchor2 = interactions[['chr2', 'start2', 'end2', 'ID']].copy()
    anchor1.columns = ['chr', 'start', 'end', 'ID']
    anchor2.columns = ['chr', 'start', 'end', 'ID']

    bed1 = bedpe[['chr1', 'start1', 'end1', 'name']].copy()
    bed2 = bedpe[['chr2', 'start2', 'end2', 'name']].copy()
    bed1.columns = ['chr', 'start', 'end', 'name']
    bed2.columns = ['chr', 'start', 'end', 'name']

    #Closest element to each anchor, with distance
    cl1 = pybedtools.BedTool.from_dataframe(anchor1).sort().closest(
              pybedtools.BedTool.from_dataframe(bed1).sort(), d=True
          ).to_dataframe(names=['chr','start','end','ID','b_chr','b_start','b_end','name','dist1'])

    cl2 = pybedtools.BedTool.from_dataframe(anchor2).sort().closest(
              pybedtools.BedTool.from_dataframe(bed2).sort(), d=True
          ).to_dataframe(names=['chr','start','end','ID','b_chr','b_start','b_end','name','dist2'])
    # Enforce same bedpe entry for both anchors, keep the closest pair per interaction
    merged = pd.merge(
        cl1[['ID', 'name', 'dist1']],
        cl2[['ID', 'name', 'dist2']],
        on=['ID', 'name'] # both have to match the same bedpe entry
    )
    # We might have multiple bedpe entries that match the same interaction, so from both anchors we take all those entries, and select the anchor that has the maximum distance (all cases)
    merged['max_dist'] = merged[['dist1', 'dist2']].max(axis=1)
    #for each interaction, keep the bedpe entry that minimizes max_dist
    merged = merged.loc[merged.groupby('ID')['max_dist'].idxmin()].reset_index(drop=True)
    return merged

def get_near_elements_count(df, resolution):
	threshold = resolution * 2
	near_anchor1 = set(df.loc[(df['dist1'] <= threshold), 'name'])
	near_anchor2 = set(df.loc[(df['dist2'] <= threshold), 'name'])

	return near_anchor1,near_anchor2


##### Compartment browser plot

def plot_compartments(compartment_dict, chrom, start, end, figsize=None, resolution=100000):
    """
    compartment_dict : dict  {label: DataFrame}
        Each DataFrame must have columns: chrom, start, end, E1.
        Example: {"GM12878": gm_df, "K562": k562_df}
    chrom : str   e.g. "chr1"
    start : int   region start (bp)
    end   : int   region end   (bp)
    figsize : tuple, optional
    resolution : int  bin size of compartment calls (default 50000)
    """
    n = len(compartment_dict)
    if figsize is None:
        figsize = (10, 1.8 * n)

    fig, axes = plt.subplots(n, 1, figsize=figsize, sharex=True)
    if n == 1:
        axes = [axes]

    colors_A = {"GM12878": "#D62728", "K562": "#D62728"}   # red  = A compartment
    colors_B = {"GM12878": "#1F77B4", "K562": "#1F77B4"}   # blue = B compartment

    for ax, (label, df) in zip(axes, compartment_dict.items()):
        region = df[(df['chrom'] == chrom) & (df['start'] >= start) & (df['end'] <= end)].copy()
        if region.empty:
            ax.text(0.5, 0.5, f"No data for {chrom}:{start}-{end}", transform=ax.transAxes, ha='center')
            ax.set_ylabel(label, fontsize=9, rotation=0, labelpad=40, va='center')
            continue

        x_mid = (region['start'].values + region['end'].values) / 2
        e1    = region['E1'].values

        # fill positive (A) and negative (B) separately
        col_A = colors_A.get(label, "#D62728")
        col_B = colors_B.get(label, "#1F77B4")
        ax.fill_between(x_mid, e1, 0, where=(e1 >= 0), interpolate=True, color=col_A, alpha=0.85, linewidth=0)
        ax.fill_between(x_mid, e1, 0, where=(e1 <  0), interpolate=True, color=col_B, alpha=0.85, linewidth=0)
        ax.plot(x_mid, e1, color='black', linewidth=0.4)
        ax.axhline(0, color='black', linewidth=0.5)

        ax.set_ylabel(label, fontsize=9, rotation=0, labelpad=40, va='center')
        ax.set_xlim(start, end)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(axis='y', labelsize=7)

        # A/B labels on right
        ymax = np.nanmax(np.abs(e1)) if len(e1) else 1
        ax.set_ylim(-ymax * 1.3, ymax * 1.3)
        ax.text(1.01, 0.75, 'A', transform=ax.transAxes, fontsize=7, color=col_A, va='center')
        ax.text(1.01, 0.25, 'B', transform=ax.transAxes, fontsize=7, color=col_B, va='center')

    # x-axis: Mb coordinates
    axes[-1].xaxis.set_major_formatter(
        plt.FuncFormatter(lambda val, _: f"{val/1e6:.1f} Mb")
    )
    axes[-1].tick_params(axis='x', labelsize=8)
    axes[-1].set_xlabel(chrom, fontsize=9)

    plt.suptitle(f"{chrom}:{start//1000}kb–{end//1000}kb  |  E1 compartment signal", fontsize=10, y=1.01)
    plt.tight_layout()
    return fig, axes


##### TAD analysis

def get_boundaries(tads):
	boundaries = tads[tads['is_boundary_250000'] == True].copy()
	print(f"Number of boundaries: {len(boundaries)}")
	tad_list = []
	for chrom, group in boundaries.groupby('chrom'):
		if chrom.startswith('chr'):
			chrom = chrom[3:]  # Remove 'chr' prefix
		group = group.sort_values('start').reset_index(drop=True)
		for i in range(len(group) - 1):
			tad_list.append({
				'chrom': f"chr{chrom}",
				'start': group.loc[i, 'end'],
				'end':   group.loc[i+1, 'start']
			})

	tad_df = pd.DataFrame(tad_list)
	return tad_df

def get_differential_tads(tads1, tads2):
    # Find TADs in tads1 with no exact match (chrom, start, end) in tads2
    merged = tads1.merge(tads2, on=['chrom', 'start', 'end'], how='left', indicator=True)
    differential_tads_df = merged[merged['_merge'] == 'left_only'][['chrom', 'start', 'end']].reset_index(drop=True)
    return differential_tads_df

def interactions_inside_tads(interactions, tad_df, resolution):
    interactions = interactions.reset_index(drop=True)

    a1 = pd.DataFrame({
        'chrom': interactions['chr'],
        'start': interactions['row'] * resolution,
        'end':   interactions['row'] * resolution + resolution,
        'idx':   interactions.index
    })
    a2 = pd.DataFrame({
        'chrom': interactions['chr'],
        'start': interactions['col'] * resolution,
        'end':   interactions['col'] * resolution + resolution,
        'idx':   interactions.index
    })

    tad_bt = pybedtools.BedTool.from_dataframe(tad_df[['chrom', 'start', 'end']])

    def get_tad_per_anchor(anchor_df):
        bt = pybedtools.BedTool.from_dataframe(anchor_df)
        result = bt.intersect(tad_bt, wa=True, wb=True, loj=True).to_dataframe(
            names=['chrom','start','end','idx','tad_chr','tad_start','tad_end'])
        result['tad_id'] = result['tad_chr'].astype(str) + '_' + result['tad_start'].astype(str) + '_' + result['tad_end'].astype(str)
        return result.drop_duplicates(subset='idx').set_index('idx')['tad_id']

    tad1 = get_tad_per_anchor(a1)
    tad2 = get_tad_per_anchor(a2)

    merged = pd.DataFrame({'tad_1': tad1, 'tad_2': tad2})
    no_match = '._-1_-1'  # bedtools loj no-match: tad_chr='.', tad_start=-1, tad_end=-1 → joined with '_'
    inside = ((merged['tad_1'] == merged['tad_2']) & (merged['tad_1'] != no_match)).mean()
    print(f"Proportion of interactions inside TADs: {inside:.4f}")
    return inside


#Comparments
def interactions_AB(interactions, AB_compartment, resolution):
    interactions = interactions.reset_index(drop=True)

    a1 = pd.DataFrame({
        'chrom': interactions['chr'],
        'start': interactions['row'] * resolution,
        'end':   interactions['row'] * resolution + resolution,
        'idx':   interactions.index
    })
    a2 = pd.DataFrame({
        'chrom': interactions['chr'],
        'start': interactions['col'] * resolution,
        'end':   interactions['col'] * resolution + resolution,
        'idx':   interactions.index
    })

    A_compartment = AB_compartment[AB_compartment['E1'] > 0]
    A_compartment = pybedtools.BedTool.from_dataframe(A_compartment[['chrom', 'start', 'end']])

    B_compartment = AB_compartment[AB_compartment['E1'] < 0]
    B_compartment = pybedtools.BedTool.from_dataframe(B_compartment[['chrom', 'start', 'end']])


    # Assign each anchor a compartment label (A or B) from a single intersect
    all_compartments = AB_compartment[['chrom', 'start', 'end', 'E1']].copy()
    all_compartments['label'] = all_compartments['E1'].apply(lambda x: 'A' if x > 0 else 'B')
    comp_bt = pybedtools.BedTool.from_dataframe(all_compartments[['chrom', 'start', 'end', 'label']])

    def get_label_per_anchor(anchor_df):
        bt = pybedtools.BedTool.from_dataframe(anchor_df)
        result = bt.intersect(comp_bt, wa=True, wb=True, loj=True).to_dataframe(
            names=['chrom','start','end','idx','comp_chr','comp_start','comp_end','label'])
        result = result.drop_duplicates(subset='idx').set_index('idx')['label']
        return result

    label1 = get_label_per_anchor(a1)
    label2 = get_label_per_anchor(a2)

    merged = pd.DataFrame({'label1': label1, 'label2': label2})
    insideA = ((merged['label1'] == 'A') & (merged['label2'] == 'A')).mean()
    insideB = ((merged['label1'] == 'B') & (merged['label2'] == 'B')).mean()
    insideAB = 1 - insideA - insideB

    print(f"Proportion of interactions inside A compartment: {insideA:.4f}")
    print(f"Proportion of interactions inside B compartment: {insideB:.4f}")
    print(f"Proportion of interactions A-B (mixed): {insideAB:.4f}")
    return insideA, insideB


##############################################################################
##### GRAVEYARD 
##############################################################################



def process_interactions(interactions:pd.DataFrame,resolution:int):
	interactions = interactions[["row",'col','chr']]
	interactions = interactions.rename(columns={'row': 'start', 'col': 'end'})

	interactions = interactions[['chr', 'start', 'end']]  # Keep only the relevant columns
	interactions['start'] = interactions['start'] * resolution
	interactions['start_end'] = interactions['start'] + resolution
	interactions['end'] = interactions['end'] * resolution
	interactions['end_end'] = interactions['end'] + resolution
	interactions = interactions[['chr', 'start', 'start_end', 'end', 'end_end']]  # Keep only the relevant columns
	interactions['ID'] = interactions['chr'] + '_' + interactions['start'].astype(str) + '_' + interactions['end'].astype(str)
	#print(interactions)
	interactions_start = interactions[['chr', 'start', 'start_end','ID']]
	interactions_end = interactions[['chr', 'end', 'end_end','ID']]
	interactions_end = interactions_end.rename(columns={'end': 'start', 'end_end': 'start_end'})

	final_interactions = pd.concat([interactions_start, interactions_end], ignore_index=True)
	final_interactions.columns = ['chr', 'start', 'end','ID']
	return final_interactions

def compute_overlap(interactions:pd.DataFrame, bed:pd.DataFrame):
	# Convert DataFrames to BedTool objects
	interactions_bed = pybedtools.BedTool.from_dataframe(interactions)
	bed_bed = pybedtools.BedTool.from_dataframe(bed)
	# Perform intersection
	overlap = interactions_bed.intersect(bed_bed, wa=True, wb=False)
	# Convert the result back to a DataFrame
	overlap_df = overlap.to_dataframe(names=interactions.columns)
	return overlap_df


def process_analysis(interactions:pd.DataFrame, files:list, resolution:int,mark:str):
	proportions = []
	for bed_file in files:
		if mark=="rE2G":
			bed = process_bed_pe(bed_file)
			interactions_processed = process_interactions_both_anchors(interactions,resolution)
		else:
			bed = process_bed(bed_file)
			interactions_processed = process_interactions(interactions,resolution)

		signal_overlap = compute_overlap(interactions_processed, bed)
		unique_overlap = signal_overlap[['chr', 'start', 'end', 'ID']].drop_duplicates() #This drops if we have multiple matches for an interaction
		unique_overlap=unique_overlap['ID'].value_counts()
		print(unique_overlap)
		if mark == "CTCF":
			unique_overlap = unique_overlap[unique_overlap > 0]  # Keep only interactions that have both anchors overlapping
		unique_overlap = unique_overlap[unique_overlap > 0]  # Keep only interactions that have both anchors overlapping
		print(f"Total Interactions: {len(interactions)}, Overlapping Interactions with {os.path.basename(bed_file)}: {len(unique_overlap)}")

		# Give proportions
		proportion = len(unique_overlap) / len(interactions) if len(interactions) > 0 else 0
		print(f"Proportion of Overlapping Interactions: {proportion:.2%}")
		proportions.append(proportion)
	return proportions,unique_overlap

def random_simulations(interactions, bed, resolution, mark,num_simulations=100):
	chr_sizes = {
		"chr1": 249250621,
		"chr2": 243199373,
		"chr3": 198022430,
		"chr4": 191154276,
		"chr5": 180915260,
		"chr6": 171115067,
        "chr7": 159138663,
        "chr8": 146364022,
        "chr9": 141213431,
        "chr10": 135534747,
        "chr11": 135006516,
        "chr12": 133851895,
        "chr13": 115169878,
        "chr14": 107349540,
        "chr15": 102531392,
        "chr16": 90354753,
        "chr17": 81195210,
        "chr18": 78077248,
        "chr19": 59128983,
        "chr20": 63025520,
        "chr21": 48129895,
        "chr22": 51304566
	}

	#select same number of interactions as in the original set
	np.random.seed(12)  # For reproducibility
	simulated_proportions = []
	for i in range(num_simulations):
		simulated_interactions = []
		for _, row in interactions.iterrows():
			distance_bin = row['col'] - row['row']
			distance = distance_bin * resolution
			chr_row = row['chr']
			chr_size = chr_sizes[chr_row]
			random_start = np.random.randint(1, chr_size - resolution)
			random_start_end = random_start + resolution

			random_end =random_start + distance
			random_end_end = random_end + resolution
			if random_start > random_end:
				random_start, random_end = random_end, random_start
				random_start_end, random_end_end = random_end_end, random_start_end

			simulated_interactions.append([chr_row, random_start, random_start_end, random_end, random_end_end])

		simulated_df = pd.DataFrame(simulated_interactions, columns=['chr', 'start', 'start_end', 'end', 'end_end'])
		simulated_df['ID'] = simulated_df['chr'] + '_' + simulated_df['start'].astype(str) + '_' + simulated_df['end'].astype(str)
		interactions_start = simulated_df[['chr', 'start', 'start_end','ID']]
		interactions_end = simulated_df[['chr', 'end', 'end_end','ID']]
		interactions_end = interactions_end.rename(columns={'end': 'start', 'end_end': 'start_end'})
		simulated_df = pd.concat([interactions_start, interactions_end], ignore_index=True)
		simulated_df.columns = ['chr', 'start', 'end','ID']
		simulated_overlap = compute_overlap(simulated_df, bed)
		unique_simulated_overlap = simulated_overlap[['chr', 'start', 'end','ID']].drop_duplicates() #This drops if we have multiple matches for an interaction
		unique_simulated_overlap=unique_simulated_overlap['ID'].value_counts()
		if mark == "CTCF":
			unique_simulated_overlap = unique_simulated_overlap[unique_simulated_overlap > 0]  # Keep only interactions that have both anchors overlapping

		unique_simulated_overlap = unique_simulated_overlap[unique_simulated_overlap > 0]  # Keep only interactions that have both anchors overlapping

		proportion = len(unique_simulated_overlap) / len(interactions) if len(interactions) > 0 else 0
		simulated_proportions.append(proportion)
		print(f"Simulation {i+1}/{num_simulations}: Proportion of Overlapping Interactions: {proportion:.2%}")
	return simulated_proportions

def random_regions(interactions,resolution):
	'''
	Docstring for random_regions

	:param interactions: data frae straight from process_DifFracTion_results, with columns 'row', 'col', and 'chr'
	:param resolution: analysis resolution, in base pairs
	:return: pd.DataFrame of random regions with the same number of interactions as in the input dataframe,
	'''
	simulated_interactions = []
	chr_sizes = {
		"chr1": 249250621,
		"chr2": 243199373,
		"chr3": 198022430,
		"chr4": 191154276,
		"chr5": 180915260,
		"chr6": 171115067,
        "chr7": 159138663,
        "chr8": 146364022,
        "chr9": 141213431,
        "chr10": 135534747,
        "chr11": 135006516,
        "chr12": 133851895,
        "chr13": 115169878,
        "chr14": 107349540,
        "chr15": 102531392,
        "chr16": 90354753,
        "chr17": 81195210,
        "chr18": 78077248,
        "chr19": 59128983,
        "chr20": 63025520,
        "chr21": 48129895,
        "chr22": 51304566
    }
	chromosomes = [f"chr{i}" for i in range(1, 23)]
	for _, row in interactions.iterrows():
			distance_bin = row['col'] - row['row']
			chr_row = row['chr']
			chr_size = chr_sizes[chr_row]
			chr_size_bin = chr_size // resolution
			random_start = np.random.randint(1, chr_size_bin)
			random_end = random_start + distance_bin
			simulated_interactions.append([chr_row, random_start, random_end])
	simulated_df = pd.DataFrame(simulated_interactions, columns=['chr', 'row', 'col'])

	return simulated_df

def random_regions_both_anchors(interactions,resolution):
	'''
	:param interactions: data frae straight from process_DifFracTion_results, with columns 'row', 'col', and 'chr'
	:param resolution: analysis resolution, in base pairs
	:return: pd.DataFrame of random regions with the same number of interactions as in the input dataframe,
	'''
	simulated_interactions = []
	chr_sizes = pd.read_csv(f"{DifFracTion_utils.main_path}/test_data/hg19_chrsizes.txt", sep='\t', header=None)
	chromosomes = [f"chr{i}" for i in range(1, 23)]

	for _, row in interactions.iterrows():
			distance_bin = row['col'] - row['row']
			chr_row = row['chr']
			chr_size = chr_sizes[chr_sizes[0] == chr_row][1].values[0]
			chr_size_bin = chr_size // resolution
			random_start = np.random.randint(1, chr_size_bin)
			random_start_end = random_start + 1
			random_end = random_start + distance_bin
			random_end_end = random_end + 1
			simulated_interactions.append([chr_row, random_start * resolution, random_start_end * resolution, chr_row,random_end * resolution, random_end_end * resolution])

	simulated_df = pd.DataFrame(simulated_interactions, columns=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
	simulated_df['ID'] = simulated_df['chr1'] + '_' + simulated_df['start1'].astype(str) + '_' + simulated_df['end1'].astype(str) + '_' + simulated_df['chr2'] + '_' + simulated_df['start2'].astype(str) + '_' + simulated_df['end2'].astype(str)

	return simulated_df