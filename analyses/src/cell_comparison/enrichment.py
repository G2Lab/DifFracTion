import argparse
import os,sys,glob
import pandas as pd
from pathlib import Path
from diffraction import utils as DifFracTion_utils
from diffraction import src as DifFracTion

sys.path.insert(0,str(Path(__file__).parent.parent))

from utils import utils

def parse_args():
     parser = argparse.ArgumentParser(description='Process Hi-C data for cell line comparison.')
     parser.add_argument('--diffraction_type', type=str, required=True, help='Type of diffraction run to process. i.e. alpha or iterative.')
     parser.add_argument('--resolution', type=int, default=10000, help='Resolution of the Hi-C matrices (default: 50000).')
     parser.add_argument('--log_fold', type=float, default=1.0, help='Log fold change threshold for filtering interactions (default: 1.0).')
     parser.add_argument('--tissue', type=str, required=True, help='Tissue name for the analysis. TC or GM.')
     parser.add_argument('--bed_file', type=str, required=True, help='Path to the BED file for overlap analysis.')
     parser.add_argument('--marker', type=str, required=True, help='Marker name for the analysis.CTCF,H3K27ac,H3K4me3,H3K9me3,H3K36me3,H3K27me3 etc')
     return parser.parse_args()



def main():
    random_sim=True
    args = parse_args()
    main_path = DifFracTion_utils.main_path
    results_path=os.path.join(main_path, "analyses/tissue_comparison/")

    save_path=(os.path.join(results_path,"marks_enrichment") )
    os.makedirs(save_path, exist_ok=True)
    main_frame = pd.DataFrame()
    chromosomes = [f"{i}" for i in range(1, 23)]

    resolution = args.resolution

    types = ["alpha", "iterative"]
    

    if args.diffraction_type not in types:
       raise ValueError(f"Invalid diffraction_type. Must be one of {types}.")
       sys.exit(1)

    positive_interactions,negative_interactions = utils.process_DifFracTion_results(
          results_path, chromosomes, resolution,log_fold=args.log_fold,type=args.diffraction_type)
    
    if args.tissue == "TC":
        args.bed_file_tissue_up = args.bed_file

        proportions_positive,unique_overlap_pos = utils.process_analysis(positive_interactions,
                                                      [args.bed_file_tissue_up],
                                                      args.resolution,args.marker)

        save_file_unique=os.path.join(save_path,f"{args.tissue}_{args.marker}_{args.diffraction_type}_unique_counts.csv")
        unique_overlap_pos.to_csv(save_file_unique,sep="\t")

        if random_sim:
            random_simulations_positive = utils.random_simulations(positive_interactions,
                                                               utils.process_bed(args.bed_file_tissue_up),
                                                               args.resolution, args.marker,num_simulations=5000)
        
        save_file=os.path.join(save_path,f"{args.tissue}_{args.marker}_{args.diffraction_type}_statistical_significance.csv")

            #merge the two lists into a dataframe and save it
        if random_sim:
            final_list = proportions_positive + random_simulations_positive
        else:
            final_list = proportions_positive
            
        df_positive = pd.DataFrame({
            'proportions': final_list
        })
        df_positive.to_csv(save_file, index=False)

    elif args.tissue == "GM":
        args.bed_file_tissue_down = args.bed_file


        proportions_negative,unique_overlap_neg = utils.process_analysis(negative_interactions,
                                                      [args.bed_file_tissue_down],
                                                      args.resolution,args.marker)
       
        # Save the unique overlap counts      
        save_file_unique=os.path.join(save_path,f"{args.tissue}_{args.marker}_{args.diffraction_type}_unique_counts.csv")
        unique_overlap_neg.to_csv(save_file_unique,sep="\t")

        if random_sim:

            random_simulations_negative = utils.random_simulations(negative_interactions,
                                                               utils.process_bed(args.bed_file_tissue_down),
                                                               args.resolution, args.marker,num_simulations=5000)
        
        save_file=os.path.join(save_path,f"{args.tissue}_{args.marker}_{args.diffraction_type}_statistical_significance.csv")
        
        if random_sim:
            final_list = proportions_negative + random_simulations_negative
        else:
            final_list = proportions_negative   
            
        df_negative = pd.DataFrame({
            'proportions': final_list
        })
       
        df_negative.to_csv(save_file, index=False)

if __name__ == "__main__":
    main()