.libPaths("./multiHiCcompare_Rlibs")
library(multiHiCcompare)
library(data.table)
library(dplyr)

sys.argv <- commandArgs(trailingOnly = TRUE)
sample_A1 <- sys.argv[1] # file1.hic 
sample_A2 <- sys.argv[2] # file2.hic
sample_B1 <- sys.argv[3] # file3.hic
sample_B2 <- sys.argv[4] # file4.hic

pval_threshold <- as.numeric(sys.argv[6]) # 0.05
spike.in.file <- sys.argv[5] # spike_ins.bed



# Functions
return_to_bin_format <- function(significant_df) {
     # Get to bin format
     width <- (significant_df$region2 - significant_df$region1) / significant_df$D

     significant_df$region1 <- significant_df$region1 / width
     significant_df$region2 <- significant_df$region2 / width

     significant_df_bins <- significant_df[, c("chr","region1","region2","logFC","p.value","p.adj")]
     return(significant_df_bins)
}

get_confusion_matrix <- function(significant_coords, spike_ins_coords, results_df) {
     significant_coords <- significant_coords[, c("region1","region2")]
     significant_coords_dt <- as.data.table(significant_coords)
     colnames(significant_coords_dt) <- c("region1","region2")

     #n_test are the bins that were included in the test (after filtering)
     n_tests <- nrow(results_df)
     #What the tool found and is a spike in
     TP <- merge(spike_ins_coords, significant_coords_dt, by=c("region1","region2"))
     #What the tool missed
     FN <- anti_join(spike_ins_coords, significant_coords_dt, by=c("region1","region2"))
     #What the tool found but is not a spike in
     FP <- anti_join(significant_coords_dt, spike_ins_coords, by=c("region1","region2"))
     TN = n_tests -  nrow(TP) - nrow(FP) - nrow(FN)
     confusion_matrix <- list("TP"=nrow(TP), "FP"=nrow(FP), "TN"=TN, "FN"=nrow(FN), "n_tests"=n_tests)
     return(confusion_matrix)
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

#load spike in coordinates
spike_ins <- fread(spike.in.file, header=FALSE)
colnames(spike_ins) <- c("region1","region2","type")
spike_ins_coords <- spike_ins[, c("region1","region2")]


# Performance evaluation
confusion_matrix_pval <- get_confusion_matrix(significant_results_pval_bins, spike_ins_coords, results)
confusion_matrix_FDR <- get_confusion_matrix(significant_results_FDR_bins, spike_ins_coords, results)

# Save confusion matrices
output.dir <- paste0(dirname(dirname(sample_A1)),"/performance/")
dir.create(output.dir,recursive=TRUE, showWarnings=FALSE)
data.table::fwrite(as.data.frame(confusion_matrix_pval),
            paste0(output.dir,'/confusion_matrix_pval_',pval_threshold,'_multiHiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)
data.table::fwrite(as.data.frame(confusion_matrix_FDR),
            paste0(output.dir,'/confusion_matrix_padj_',pval_threshold,'_multiHiCcompare.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)