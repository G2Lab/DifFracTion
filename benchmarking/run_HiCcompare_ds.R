#ml R/4.5.2-mba
# HiC-Dip/benchmarking/run_HiCcompare.R

.libPaths("./HiCcompare_Rlibs")
library(HiCcompare)
library(data.table)
# Generate a table for comparisons
# We will generate the data from python using hic.straw so we dont have to load large matrices in R
# Previously we have saved the table with the spike ins and the neighbors to be identified as significant interactions
# We need to load this table too and then see if we can identify the spike ins


# Inputs 
sys.argv <- commandArgs(trailingOnly = TRUE)

table.file <- sys.argv[1] # file1.hic
pval_threshold <- as.numeric(sys.argv[2]) # 0.05

# Functions
return_to_bin_format <- function(significant_df) {
     # Get to bin format
     width <- significant_df$end1 - significant_df$start1 
     significant_df$start1 <- significant_df$start1 
     significant_df$start1 <- significant_df$start1 / width
     significant_df$end1 <- significant_df$end1 / width
     significant_df$start2 <- significant_df$start2 
     significant_df$start2 <- significant_df$start2 / width
     significant_df$end2 <- significant_df$end2 / width
     significant_df_bins <- significant_df[, c("chr1","start1","end1","chr2","start2","end2","M","p.value","p.adj")]
     return(significant_df_bins)
}

get_FPR <- function(significant_coords, results_df) {
     significant_coords <- significant_coords[, c("start1","start2")]
     significant_coords_dt <- as.data.table(significant_coords)
     colnames(significant_coords_dt) <- c("start1","end1")
     #n_test are the bins that were included in the test (after filtering)
     n_tests <- nrow(results_df)
     #What the tool found and is a spike in
     FP <- nrow(significant_coords_dt) # all significant coords are false positives since there are no true positives in this set up
     FPR <- if (n_tests > 0) FP / n_tests else 0
     tab_res <- list("FP"=FP, "FPR"=FPR, "n_tests"=n_tests)
     return(tab_res)
}

# chr1 start1 end1 chr2 start2 end2 IF1 IF2 D(distance in bins) M(log2 ratio IF2/IF1)

chr.table <- fread(table.file, header=TRUE)
loess.hic <- hic_loess(chr.table, Plot=TRUE)
results <- hic_compare(loess.hic) # Filters interactions based on a threshold
significant_results_pval <- results[results$p.value <= pval_threshold, ] # Filter significant interactions
significant_results_FDR <- results[results$p.adj <= pval_threshold, ] # Filter significant interactions

# Return to bin format
significant_results_pval_bins <- return_to_bin_format(significant_results_pval)
significant_results_FDR_bins <- return_to_bin_format(significant_results_FDR)


# Performance evaluation
FPR_pval <- get_FPR(significant_results_pval_bins, results)
FPR_FDR <- get_FPR(significant_results_FDR_bins, results)

# Save FPR
output.dir <- paste0(dirname(dirname(table.file)),"/performance/")
dir.create(output.dir,recursive=TRUE, showWarnings=FALSE)
data.table::fwrite(as.data.frame(FPR_pval),
            paste0(output.dir,'/FPR_pval_',pval_threshold,'_HiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)
data.table::fwrite(as.data.frame(FPR_FDR),
            paste0(output.dir,'/FPR_padj_',pval_threshold,'_HiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)