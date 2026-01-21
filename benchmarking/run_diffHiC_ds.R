#ml R/4.5.2-mba
# HiC-Dip/benchmarking/run_diffHiC_ds.R

.libPaths("./diffHiC_Rlibs")
library(diffHic)
library(edgeR)
library(InteractionSet)
library(GenomicRanges)
library(data.table)
library(dplyr)
# Inputs
sys.argv <- commandArgs(trailingOnly = TRUE)
table.file <- sys.argv[1] # file1.hic
pval_threshold <- as.numeric(sys.argv[2]) # 0.05



#Functions
generate_GInteractions <- function(chr.table) {
     anchor1 <- GRanges(chr.table$chr1,
                    IRanges(chr.table$start1+1, chr.table$end1))

     anchor2 <- GRanges(chr.table$chr2,
                    IRanges(chr.table$start2+1, chr.table$end2))
     gi <- GInteractions(anchor1, anchor2)

     # We assume that counts start from column 7 to the end
     count.mat <- as.matrix(chr.table[, 7:ncol(chr.table)]) # extract count columns and convert to matrix
     # it is actually not a matrix, is just a samplesxinteractions data.frame

     iset <- InteractionSet(assays = list(counts = count.mat),
                       interactions = gi)
     colData(iset)$totals <- colSums(count.mat)  # library size normalization

     #Regions are the unique bins
     #Rows are interactions
     return(iset)
}

return_to_bin_format <- function(significant_df) {
     # Get to bin format
     significant_df$start1 <- significant_df$start1 - 1
     significant_df$start1 <- significant_df$start1 / significant_df$width1
     significant_df$end1 <- significant_df$end1 / significant_df$width1
     significant_df$start2 <- significant_df$start2 - 1
     significant_df$start2 <- significant_df$start2 / significant_df$width2
     significant_df$end2 <- significant_df$end2 / significant_df$width2
     significant_df_bins <- significant_df[, c("seqnames1","start1","end1","seqnames2","start2","end2","logFC","PValue","FDR")]
     return(significant_df_bins)
}

get_significant_coordinates <- function(significant_results_bins) {
     coords <- significant_results_bins[, c("start1","start2")]
     return(coords)
}

get_FPR <- function(significant_coords, iset) {

     significant_coords_dt <- as.data.table(significant_coords)
     colnames(significant_coords_dt) <- c("start1","end1")

     #n_test are the bins that were included in the test (after filtering)
     n_tests <- dim(iset)[1]
     FP <- nrow(significant_coords_dt) # all significant coords are false positives since there are no true positives in this set up
     FPR <- if (n_tests > 0) FP / n_tests else 0
     tab_res <- list("FP"=FP, "FPR"=FPR, "n_tests"=n_tests)
     return(tab_res)
}

#Load table upper triangle
chr.table <- fread(table.file, header=TRUE)
iset <- generate_GInteractions(chr.table)

# This line literally swaps the interactions to ensure that the larger coordinate is first
interactions(iset) <- as(interactions(iset), "ReverseStrictGInteractions")
direct <- filterDirect(iset)
#abundance is log2 counts per million
#Retain interactions that are at least two fold above the log2 CPM threshold
# because the threshold is treated as background noise and adding 5 log2 CPM ensures
# we are above that noise (remember we are on log2 scale so adding log2(5) is multiplying by 5 in normal scale)
direct.keep <- direct$abundances > log2(2) + direct$threshold
iset <- iset[direct.keep, ]

y <- asDGEList(iset)
y <- csaw::normOffsets(y)
y <- estimateDisp(y)

# n_samples = 2 and design is rank 2, so theres no residual degrees of freedom
# We included replicates using the downsampling method with the respective metrics 
group <- factor(c("A","A","B","B"))
design <- model.matrix(~group)

fit <- glmQLFit(y, design)
res <- glmQLFTest(fit, coef=2)

res_table <- topTags(res, n= nrow(res$table), sort.by = "none")$table

# need to convert back to proper orientation
coords <- as.data.frame(as(interactions(iset), "StrictGInteractions"))
res_table <- cbind(coords, res_table)

significant_results <- res_table[res_table$PValue <= pval_threshold, ]
significant_results_FDR <- significant_results[significant_results$FDR <= pval_threshold, ]

significant_results_bins <- return_to_bin_format(significant_results)
significant_results_FDR_bins <- return_to_bin_format(significant_results_FDR)

significant_coords <- get_significant_coordinates(significant_results_bins)
significant_coords_FDR <- get_significant_coordinates(significant_results_FDR_bins)


#confusion_matrix_pval$TP,confusion_matrix_pval$TN,confusion_matrix_pval$FP,confusion_matrix_pval$FN,confusion_matrix_pval$n_tests
respval <- get_FPR(significant_coords, iset)
res_FDR <- get_FPR(significant_coords_FDR, iset)
# Save confusion matrices
output.dir <- paste0(dirname(dirname(table.file)),"/performance/")
dir.create(output.dir,recursive=TRUE, showWarnings=FALSE)
data.table::fwrite(as.data.frame(respval),
            paste0(output.dir,'/FPR_pval_',pval_threshold,'_diffHiC.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)
data.table::fwrite(as.data.frame(res_FDR),
            paste0(output.dir,'/FPR_padj_',pval_threshold,'_diffHiC.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)