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
pval_threshold <- as.numeric(sys.argv[3]) # 0.05
spike.in.file <- sys.argv[2] # spike_ins.bed

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

get_confusion_matrix <- function(significant_coords, spike_ins_coords, results_df) {
     significant_coords <- significant_coords[, c("start1","start2")]
     significant_coords_dt <- as.data.table(significant_coords)
     colnames(significant_coords_dt) <- c("start1","end1")
     #n_test are the bins that were included in the test (after filtering)
     n_tests <- nrow(results_df)
     #What the tool found and is a spike in
     TP <- merge(spike_ins_coords, significant_coords_dt, by=c("start1","end1"))
     #What the tool missed
     FN <- anti_join(spike_ins_coords, significant_coords_dt, by=c("start1","end1"))
     #What the tool found but is not a spike in
     FP <- anti_join(significant_coords_dt, spike_ins_coords, by=c("start1","end1"))
     TN = n_tests -  nrow(TP) - nrow(FP) - nrow(FN)
     confusion_matrix <- list("TP"=nrow(TP), "FP"=nrow(FP), "TN"=TN, "FN"=nrow(FN), "n_tests"=n_tests)
     return(confusion_matrix)
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


# Load spike ins
spike_ins <- fread(spike.in.file, header=FALSE)
colnames(spike_ins) <- c("start1","end1","type")
spike_ins_coords <- spike_ins[, c("start1","end1")]

# Performance evaluation
confusion_matrix_pval <- get_confusion_matrix(significant_results_pval_bins, spike_ins_coords, results)
confusion_matrix_FDR <- get_confusion_matrix(significant_results_FDR_bins, spike_ins_coords, results)

# Save confusion matrices
output.dir <- paste0(dirname(dirname(table.file)),"/performance/")
dir.create(output.dir,recursive=TRUE, showWarnings=FALSE)
data.table::fwrite(as.data.frame(confusion_matrix_pval),
            paste0(output.dir,'/confusion_matrix_pval_',pval_threshold,'_HiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)
data.table::fwrite(as.data.frame(confusion_matrix_FDR),
            paste0(output.dir,'/confusion_matrix_padj_',pval_threshold,'_HiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)