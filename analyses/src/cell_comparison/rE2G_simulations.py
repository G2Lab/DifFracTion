import argparse
import os,sys,glob
import pandas as pd
from pathlib import Path
from diffraction import utils as DifFracTion_utils
from diffraction import src as DifFracTion

sys.path.insert(0,str(Path(__file__).parent.parent))

from utils import utils


def parse_args():
     parser = argparse.ArgumentParser(description='TAD simulations enrichment.')
     parser.add_argument('--diffraction_type', type=str, required=True, help='Type of diffraction run to process. i.e. alpha or iterative.')
     parser.add_argument('--resolution', type=int, default=10000, help='Resolution of the Hi-C matrices (default: 10000).')
     parser.add_argument('--log_fold', type=float, default=2.0, help='Log fold change threshold for filtering interactions (default: 1.0).')
     parser.add_argument('--cell_line', type=str, required=True, help='Cell line name for the analysis. GM12878 or K562')
     return parser.parse_args()


def main():
     args = parse_args()

     chromosomes = [f"{i}" for i in range(1, 23)]
     main_path = DifFracTion_utils.main_path
     results_path=os.path.join(main_path, "analyses/cell_comparison/")


     positive_interactions,negative_interactions = utils.process_DifFracTion_results(
			results_path, chromosomes, args.resolution,log_fold=args.log_fold,type=args.diffraction_type)
	
     #### Load rE2G bedpe files and differential bedpe files
     bed_file_positive=f"{main_path}/test_data/rE2G/GM12878.hg19.bedpe"
     bed_file_negative=f"{main_path}/test_data/rE2G/K562.hg19.bedpe"
     bed_positive = utils.process_bed_pe(bed_file_positive)
     bed_negative = utils.process_bed_pe(bed_file_negative)
     ####### 



     if args.cell_line == 'GM12878':
          interactions = positive_interactions
          differential_bed = bed_positive[~bed_positive[['chr2', 'start2', 'end2']].isin(bed_negative[['chr2', 'start2', 'end2']]).any(axis=1)]
     elif args.cell_line == 'K562':
          interactions = negative_interactions
          differential_bed = bed_negative[~bed_negative[['chr2', 'start2', 'end2']].isin(bed_positive[['chr2', 'start2', 'end2']]).any(axis=1)]
     else:
          raise ValueError(f"Unknown cell line: {args.cell_line}. Must be GM12878 or K562.")

     random_props = []
     for i in range(5000):
          random_interactions = utils.random_regions_both_anchors(interactions, args.resolution)
          distances = utils.closest_bedpe_distance(random_interactions, differential_bed)
          near_anchor1, near_anchor2 = utils.get_near_elements_count(distances, args.resolution)
          random_near = len(near_anchor1 | near_anchor2)
          random_props.append((int(random_near) / len(differential_bed))*100)
     
     tad_path_sims = os.path.join(results_path, "rE2G_simulations")

     os.makedirs(tad_path_sims, exist_ok=True)
     output_path = os.path.join(tad_path_sims, f"rE2G_simulations_{args.diffraction_type}_{args.cell_line}.tsv")
     pd.DataFrame({'prop_w_EP': random_props}).to_csv(output_path, sep='\t', index=False)
     print(f"Saved results to {output_path}")

if __name__ == "__main__":
    main()