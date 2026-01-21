.libPaths("./HiCDPlus_Rlibs")
library(HiCDCPlus)
library(data.table)
library(dplyr)
library(DESeq2)


#Only with p_adj because it requires multiple operations

sys.argv <- commandArgs(trailingOnly = TRUE)
sample_A1 <- sys.argv[1] 
sample_A2 <- sys.argv[2] 
sample_B1 <- sys.argv[3] 
sample_B2 <- sys.argv[4] 

pval_threshold <- as.numeric(sys.argv[5]) #0.01


return_to_bin_format <- function(significant_df,resolution) {
     # Get to bin format
     significant_df$startI <- significant_df$startI / resolution
     significant_df$startJ <- significant_df$startJ / resolution
     significant_df_bins <- significant_df[, c("startI","startJ","log2FoldChange","pvalue","padj")]
     return(significant_df_bins)
}

get_FPR <- function(significant_coords, n_tests) {
     significant_coords <- significant_coords[, c("startI","startJ")]
     significant_coords_dt <- as.data.table(significant_coords)
     colnames(significant_coords_dt) <- c("startI","startJ")
     FP <- nrow(significant_coords_dt) # all significant coords are false positives since there are no true positives in this set up
     FPR <- if (n_tests > 0) FP / n_tests else 0
     tab_res <- list("FP"=FP, "FPR"=FPR, "n_tests"=n_tests)
     return(tab_res)
}

A1 <- fread(sample_A1, header=TRUE)
A2 <- fread(sample_A2, header=TRUE)
B1 <- fread(sample_B1, header=TRUE)
B2 <- fread(sample_B2, header=TRUE)
chromosome <- unique(A1$chr)
resolution <- A1$startJ[1] - A1$startI[1]

hicfile_paths <- c(sample_A1, sample_A2, sample_B1, sample_B2)

                    
indexfile<-data.frame()
counter = 1
for (df in list(A1, A2, B1, B2)) {
     
     current_hicfile <- hicfile_paths[counter]
     directory_c <- dirname(hicfile_paths[counter])
     output_path<-paste0(directory_c,"/",
                    gsub("^(.*[\\/])", "",gsub('.table','.txt.gz',current_hicfile)))

     features <- construct_features(
          output_path = paste0("./", chromosome, "_", resolution, "_features"),
          chrs      = chromosome,
          binsize   = resolution,
          gen        = "Hsapiens",
          gen_ver    = "hg19")

     df$counts <- round(df$counts)

     gi_list<-generate_bintolen_gi_list(
     bintolen_path=features)
     gi_list[[chromosome]] <- add_2D_features(gi_list[[chromosome]], as.data.frame(df))

     gi_list<-expand_1D_features(gi_list)
     gi_list <- lapply(gi_list, function(gi) {
     mcols(gi)$gc  <- mcols(gi)$gc  + rnorm(length(gi), sd = 1e-6)
     mcols(gi)$len <- mcols(gi)$len + rnorm(length(gi), sd = 1e-6)
     gi
     })
          
     set.seed(1010) #HiC-DC downsamples rows for modeling
     gi_list<-HiCDCPlus(gi_list, ssize = 1)

     for (i in seq(length(gi_list))) {
          indexfile<- unique(rbind(indexfile,
          as.data.frame(gi_list[[i]][gi_list[[i]]$qvalue <= pval_threshold] )[c('seqnames1','start1','start2')]))
     }
     gi_list_write(gi_list,
          fname = output_path)
     counter = counter + 1
}

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
diff_results <- return_to_bin_format(diff_results,resolution)
output.dir <- paste0(dirname(directory_c),"/performance/")
dir.create(output.dir,recursive=TRUE, showWarnings=FALSE)


n_tests <- length(gi_list[[chromosome]])
     
FPR <- get_FPR(diff_results,n_tests)
FPR <- as.data.frame(FPR)


# Save FPR
data.table::fwrite(FPR,
            paste0(output.dir,'/FPR_padj_',pval_threshold,'_HiCDCPlus.txt')
            ,sep='\t',row.names=FALSE,quote=FALSE)