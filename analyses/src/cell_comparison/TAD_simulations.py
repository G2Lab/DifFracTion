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
     parser.add_argument('--type', type=str, required=True, help='Cell line or Tissue comparison. Options: cell_line or tissue')
     
     #Conditional argument based on the type of comparison
     known, _ = parser.parse_known_args()
     if known.type == 'cell_line':
          parser.add_argument('--cell_line', type=str, required=True, help='Cell line name for the analysis. GM12878 or K562')
     elif known.type == 'tissue':
          parser.add_argument('--tissue', type=str, required=True, help='Tissue name for the analysis. Options: liver, heart, brain, etc.')
     return parser.parse_args()


def main():
     args = parse_args()

     chromosomes = [f"{i}" for i in range(1, 23)]
     main_path = DifFracTion_utils.main_path
     if args.type == 'cell_line':
          results_path=os.path.join(main_path, "analyses/cell_comparison/")
          tad_path = os.path.join(results_path, "TADs")

          c1_tads = pd.read_csv(os.path.join(tad_path, "GM12878-HRC_TADs.tsv"), sep="\t", header=0)
          c2_tads = pd.read_csv(os.path.join(tad_path, "K562-HRC_TADs.tsv"), sep="\t", header=0)
     

     elif args.type == 'tissue':
          results_path=os.path.join(main_path, "analyses/tissue_comparison/")
          tad_path = os.path.join(results_path, "TADs")
          c1_tads = pd.read_csv(os.path.join(tad_path, "TransverseColon-ENCDO845WKR_TADs.tsv"), sep="\t", header=0)
          c2_tads = pd.read_csv(os.path.join(tad_path, "GastrocnemiusMedialis-ENCDO845WKR_TADs.tsv"), sep="\t", header=0)

     

     positive_interactions,negative_interactions = utils.process_DifFracTion_results(
			results_path, chromosomes, args.resolution,log_fold=args.log_fold,type=args.diffraction_type)
	
     c1_tads_boudaries = utils.get_boundaries(c1_tads)
     c2_tads_boudaries = utils.get_boundaries(c2_tads)

     diff_c1_tads = utils.get_differential_tads(c1_tads_boudaries, c2_tads_boudaries)
     diff_c2_tads = utils.get_differential_tads(c2_tads_boudaries, c1_tads_boudaries)

     if args.type == 'cell_line':
          if args.cell_line == 'GM12878':
               interactions = positive_interactions
               diff_tads = diff_c1_tads

          elif args.cell_line == 'K562':
               interactions = negative_interactions
               diff_tads = diff_c2_tads
          else:
               raise ValueError(f"Unknown cell line: {args.cell_line}. Must be GM12878 or K562.")
     
     elif args.type == 'tissue':
          if args.tissue == 'TransverseColon':
               interactions = positive_interactions
               diff_tads = diff_c1_tads

          elif args.tissue == 'GastrocnemiusMedialis':
               interactions = negative_interactions
               diff_tads = diff_c2_tads
          else:
               raise ValueError(f"Unknown tissue: {args.tissue}. Must be TransverseColon or GastrocnemiusMedialis.")

     random_props = []
     for i in range(5000):
          random_interactions = utils.random_regions(interactions, args.resolution)
          prop = utils.interactions_inside_tads(random_interactions, diff_tads, args.resolution)
          random_props.append(prop)
     
     tad_path_sims = os.path.join(results_path, "TAD_simulations")

     os.makedirs(tad_path_sims, exist_ok=True)
     if args.type == 'cell_line':
          output_path = os.path.join(tad_path_sims, f"TAD_simulations_{args.diffraction_type}_{args.cell_line}.tsv")
     elif args.type == 'tissue':
          output_path = os.path.join(tad_path_sims, f"TAD_simulations_{args.diffraction_type}_{args.tissue}.tsv")
     pd.DataFrame({'prop_inside_tad': random_props}).to_csv(output_path, sep='\t', index=False)
     print(f"Saved results to {output_path}")

if __name__ == "__main__":
    main()