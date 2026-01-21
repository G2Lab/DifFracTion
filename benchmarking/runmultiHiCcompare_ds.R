.libPaths("./multiHiCcompare_Rlibs")
library(multiHiCcompare)
library(data.table)
library(dplyr)

sys.argv <- commandArgs(trailingOnly = TRUE)
sample_A1 <- sys.argv[1] # file1.hic 
sample_A2 <- sys.argv[2] # file2.hic
sample_B1 <- sys.argv[3] # file3.hic
sample_B2 <- sys.argv[4] # file4.hic
pval_threshold <- as.numeric(sys.argv[5]) # 0.05




# Functions
return_to_bin_format <- function(significant_df) {
     # Get to bin format
     width <- (significant_df$region2 - significant_df$region1) / significant_df$D

     significant_df$region1 <- significant_df$region1 / width
     significant_df$region2 <- significant_df$region2 / width

     significant_df_bins <- significant_df[, c("chr","region1","region2","logFC","p.value","p.adj")]
     return(significant_df_bins)
}

get_FPR <- function(significant_coords, results_df) {
     significant_coords <- significant_coords[, c("region1","region2")]
     significant_coords_dt <- as.data.table(significant_coords)
     colnames(significant_coords_dt) <- c("region1","region2")

     #n_tests are the bins that were included in the test (after filtering)
     n_tests <- nrow(results_df)
     FP <- nrow(significant_coords_dt)
     FPR <- FP / n_tests
     tab_res <- list("FP"=FP, "FPR"=FPR, "n_tests"=n_tests)
     return(tab_res)
}

A1 <- fread(sample_A1, header=TRUE)
A2 <- fread(sample_A2, header=TRUE)
B1 <- fread(sample_B1, header=TRUE)
B2 <- fread(sample_B2, header=TRUE)

hicexp <- make_hicexp(A1, A2, B1, B2, groups = c(1, 1, 2, 2))
hicexp <- cyclic_loess(hicexp)
hicexp <- hic_exactTest(hicexp)

results <- results(hicexp)
significant_results_pval <- results[results$p.value <= pval_threshold, ] # Filter significant interactions
significant_results_FDR <- results[results$p.adj <= pval_threshold, ] #

#Bin format
significant_results_pval_bins <- return_to_bin_format(significant_results_pval)
significant_results_FDR_bins <- return_to_bin_format(significant_results_FDR)

# Performance evaluation
FPR_pval <- get_FPR(significant_results_pval_bins, results)
FPR_FDR <- get_FPR(significant_results_FDR_bins, results)

# Save confusion matrices
output.dir <- paste0(dirname(dirname(sample_A1)),"/performance/")
dir.create(output.dir,recursive=TRUE, showWarnings=FALSE)
data.table::fwrite(as.data.frame(FPR_pval),
            paste0(output.dir,'/FPR_pval_',pval_threshold,'_multiHiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)
data.table::fwrite(as.data.frame(FPR_FDR),
            paste0(output.dir,'/FPR_padj_',pval_threshold,'_multiHiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)