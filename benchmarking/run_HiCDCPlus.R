library(HiCDCPlus)
library(data.table)
library(dplyr)
library(DESeq2)


#Only with p_adj because it requires multiple operations
start_time <- proc.time()

sys.argv <- commandArgs(trailingOnly = TRUE)
sample_A1 <- sys.argv[1]
sample_A2 <- sys.argv[2] # file2.hic  
sample_B1 <- sys.argv[3] # file2.hic  
sample_B2 <- sys.argv[4] # file2.hic 


arguments <-basename(dirname(dirname(sample_A1)))
arguments_split <- strsplit(arguments, "_")[[1]]
chr_number <- arguments_split[1]
resolution <- as.integer(arguments_split[2]) 
if (length(arguments_split) > 3) {
     spike_k <- arguments_split[4]
} else {
     spike_k <- arguments_split[3]
}

pval_threshold <- as.numeric(sys.argv[6]) #0.01
spike.in.file <- sys.argv[5] 

return_to_bin_format <- function(significant_df,resolution) {
     # Get to bin format
     significant_df$startI <- significant_df$startI / resolution
     significant_df$startJ <- significant_df$startJ / resolution
     significant_df_bins <- significant_df[, c("startI","startJ","log2FoldChange","pvalue","padj")]
     return(significant_df_bins)
}

get_confusion_matrix <- function(significant_coords, spike_ins_coords, n_tests) {
     significant_coords <- significant_coords[, c("startI","startJ")]
     significant_coords_dt <- as.data.table(significant_coords)
     colnames(significant_coords_dt) <- c("startI","startJ")


     #What the tool found and is a spike in
     TP <- merge(spike_ins_coords, significant_coords_dt, by=c("startI","startJ"))
     #What the tool missed
     FN <- anti_join(spike_ins_coords, significant_coords_dt, by=c("startI","startJ"))
     #What the tool found but is not a spike in
     FP <- anti_join(significant_coords_dt, spike_ins_coords, by=c("startI","startJ"))
     TN = n_tests -  nrow(TP) - nrow(FP) - nrow(FN)
     confusion_matrix <- list("TP"=nrow(TP), "FP"=nrow(FP), "TN"=TN, "FN"=nrow(FN), "n_tests"=n_tests)
     return(confusion_matrix)
}

A1 <- fread(sample_A1, header=TRUE)
A2 <- fread(sample_A2, header=TRUE)
B1 <- fread(sample_B1, header=TRUE)
B2 <- fread(sample_B2, header=TRUE)
chromosome <- unique(A1$chr)
#resolution <- A1$startJ[1] - A1$startI[1]

hicfile_paths <- c(sample_A1, sample_A2, sample_B1, sample_B2)

                    
directory_c <- dirname(hicfile_paths[1])
features <- construct_features(
     output_path = paste0(directory_c, "/", chromosome, "_", resolution, "_features"),
     chrs      = chromosome,
     binsize   = resolution,
     gen        = "Hsapiens",
     gen_ver    = "hg19")

indexfile<-data.frame()
counter = 1
for (df in list(A1, A2, B1, B2)) {

     current_hicfile <- hicfile_paths[counter]
     directory_c <- dirname(hicfile_paths[counter])
     output_path<-paste0(directory_c,"/",
                    gsub("^(.*[\\/])", "",gsub('.table','.txt.gz',current_hicfile)))

     df$counts <- round(df$counts)

     # Initialize GenomicInteractions object with bin-level 1D features (GC content, mappability) from reference genome
     # Is a list of all genomic intervals in the chromosome, with their features (gc, mappability) and will be updated with observed contact counts and fitted model results for each replicate
     gi_list<-generate_bintolen_gi_list(
     bintolen_path=features)

     # Add observed contact counts from this replicate's Hi-C matrix to the genomic intervals
     gi_list[[chromosome]] <- add_2D_features(gi_list[[chromosome]], as.data.frame(df))

     # Propagate 1D bin features (gc, len) to each interaction pair (bin_i x bin_j)
     gi_list<-expand_1D_features(gi_list)

     # Fit negative binomial model per interaction: expected counts ~ distance + GC + mappability
     # qvalue = whether observed counts significantly exceed expected background (interaction signal)
     set.seed(1010) #HiC-DC downsamples rows for modeling
     gi_list<-HiCDCPlus(gi_list,ssize=0.1)

     # Collect bin pairs with significant signal (qvalue <= threshold) across all replicates
     # This builds the whitelist of interactions to test for differential analysis
     for (i in seq(length(gi_list))) {
          indexfile<- unique(rbind(indexfile,
          as.data.frame(gi_list[[i]][gi_list[[i]]$qvalue <= pval_threshold] )[c('seqnames1','start1','start2')]))
     }
     # Save fitted gi_list for this replicate (used as input to hicdcdiff)
     gi_list_write(gi_list,
          fname = output_path)
     counter = counter + 1
}

# This is the filtered file list of interactions to test for differential analysis, which is the union of significant interactions across all replicates (signal in at least one replicate)
colnames(indexfile)<-c('chr','startI','startJ')

data.table::fwrite(indexfile,
            paste0(directory_c,'/analysis_indices.txt.gz'),
            sep='\t',row.names=FALSE,quote=FALSE)


current_hicfile <- hicfile_paths[1]
directory_c <- dirname(hicfile_paths[1])
output_path<-paste0(directory_c,"/",
                    gsub("^(.*[\\/])", "",gsub('.table','.txt.gz',current_hicfile)))


file_A1 <-  hicfile_paths[1]
gi_A_path <- paste0(dirname(file_A1),"/",
                    gsub("^(.*[\\/])", "",gsub('.table','.txt.gz',file_A1)))

file_A2 <-  hicfile_paths[2]
gi_A2_path <- paste0(dirname(file_A2),"/",
                    gsub("^(.*[\\/])", "",gsub('.table','.txt.gz',file_A2)))

file_B1 <-  hicfile_paths[3]
gi_B_path <- paste0(dirname(file_B1),"/",
                    gsub("^(.*[\\/])", "",gsub('.table','.txt.gz',file_B1)))

file_B2 <-  hicfile_paths[4]
gi_B2_path <- paste0(dirname(file_B2),"/",
                    gsub("^(.*[\\/])", "",gsub('.table','.txt.gz',file_B2)))

# Differential analysis
hicdcdiff(input_paths=list(A=c(gi_A_path, gi_A2_path),B=c(gi_B_path, gi_B2_path)),
          filter_file=paste0(directory_c,'/analysis_indices.txt.gz'),
          output_path=paste0(directory_c,'/HiCDCPlus_differential_results/'),
          fitType='mean',binsize=resolution)

# Load differential results
diff_results <- fread(paste0(directory_c,'/HiCDCPlus_differential_results/',"diff_resBoverA_",chromosome,".txt.gz"), header=TRUE)
diff_results <- diff_results[diff_results$padj <= pval_threshold, ]
diff_results <- return_to_bin_format(diff_results,resolution)
output.dir <- paste0(dirname(directory_c),"/performance/")
dir.create(output.dir,recursive=TRUE, showWarnings=FALSE)

# Read spike ins
spike_ins <- fread(spike.in.file, header=FALSE)
colnames(spike_ins) <- c("startI","startJ","type")
spike_ins_coords <- spike_ins[, c("startI","startJ")]

n_tests <- nrow(indexfile)

confusion_matrix <- get_confusion_matrix(diff_results,
                                        spike_ins_coords,
                                        n_tests)
confusion_matrix <- as.data.frame(confusion_matrix)
# Save confusion matrix
data.table::fwrite(confusion_matrix,
            paste0(output.dir,'/confusion_matrix_padj_',pval_threshold,'_HiCDCPlus.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)
elapsed  <- proc.time() - start_time
mem_mb   <- gc(reset = TRUE)[2, 6]  # max used heap in Mb
perf_log <- data.frame(  
    tool        = "HiCDCPlus",
    chr         = chr_number,
    resolution  = resolution,
    spike_k     = spike_k,
    elapsed_sec = elapsed["elapsed"],
    mem_mb      = mem_mb
)
data.table::fwrite(perf_log, paste0(output.dir, "/runtime_HiCDCPlus.txt"), sep = "\t", row.names = FALSE)