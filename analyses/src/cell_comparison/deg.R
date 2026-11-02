library(DESeq2)
library(argparse)
library(org.Hs.eg.db)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicFeatures)

parser <- ArgumentParser()
parser$add_argument("--counts_a", required=TRUE, help="Raw counts file for condition A (e.g. GM12878.rawCounts)")
parser$add_argument("--counts_b", required=TRUE, help="Raw counts file for condition B (e.g. K562.rawCounts)")
parser$add_argument("--name_a",   required=TRUE, help="Name for condition A (e.g. GM12878)")
parser$add_argument("--name_b",   required=TRUE, help="Name for condition B (e.g. K562)")
parser$add_argument("--output",   required=TRUE, help="Output TSV path")
parser$add_argument("--lfc",      type="double", default=1.0, help="Log2 fold change cutoff (default: 1.0)")
parser$add_argument("--padj",     type="double", default=0.05, help="Adjusted p-value cutoff (default: 0.05)")
args <- parser$parse_args()

read_counts <- function(path) {
    df <- read.table(path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
    df
}

counts_a <- read_counts(args$counts_a)
counts_b <- read_counts(args$counts_b)

counts <- merge(counts_a, counts_b, by="row.names", all=FALSE)
rownames(counts) <- counts$Row.names
counts$Row.names <- NULL
counts <- round(counts)

n_a <- ncol(counts_a)
n_b <- ncol(counts_b)
condition <- factor(c(rep(args$name_a, n_a), rep(args$name_b, n_b)))
col_data <- data.frame(condition=condition, row.names=colnames(counts))

dds <- DESeqDataSetFromMatrix(countData=counts, colData=col_data, design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition", args$name_b, args$name_a),
               lfcThreshold=0, alpha=args$padj)
res <- lfcShrink(dds, contrast=c("condition", args$name_b, args$name_a),
                 res=res, type="ashr")

res_df <- as.data.frame(res)

# We need to actually identify what type of ID the rownames are. In this case, they are Entrez IDs. We will map them to gene symbols and coordinates.
tocheck <- rownames(res_df)
tocheck_clean <- sub("\\..*$", "", tocheck)  # strip version suffixes e.g. ENSG00000000003.14
head(tocheck_clean)
chek_names_fun <- function(x) {
    if (all(grepl("^[0-9]+$", x))) {
        return("ENTREZID")
    } else if (all(grepl("^ENSG[0-9]+$", x))) {
        return("ENSEMBL")
    } else {
        return("UNKNOWN")
    }
}
rownames(res_df) <- tocheck_clean  # use clean IDs from here on
res_df$GeneID <- rownames(res_df)
idtype <- chek_names_fun(tocheck_clean)
if (idtype == "UNKNOWN") {
    stop("Unknown gene ID type. Please provide Entrez IDs or Ensembl IDs.")
}
#Entrez ID 2 gene symbol
res_df$GeneName <- mapIds(org.Hs.eg.db,
                          keys=res_df$GeneID,
                          column="SYMBOL",
                          keytype=idtype,
                          multiVals="first")

# Map Entrez ID coordinates ()
# Use our GENEID (ENSEMBL in EnsDb.Hsapiens.v86) to get the coordinates from EnsDb.Hsapiens.v86.
# Keep the columns we need: GeneID, chr, start, end, strand. Filter out non-standard chromosomes and duplicates.

if (idtype == "ENSEMBL") {
    for_next <- "GENEID"
} else if (idtype == "ENTREZID") {
    for_next <- "ENTREZID"
}

coords <- select(EnsDb.Hsapiens.v86,
                 keys=res_df$GeneID,
                 columns=c(for_next, "SEQNAME", "GENESEQSTART", "GENESEQEND", "SEQSTRAND"),
                 keytype=for_next)

colnames(coords) <- c("GeneID", "chr", "start", "end", "strand")
coords <- coords[coords$chr %in% c(as.character(1:22), "X", "Y"), ]
coords <- coords[!duplicated(coords$GeneID), ]

# Merge coordinates with res_df (from DESeq2 results)
res_df <- merge(res_df, coords, by="GeneID", all.x=TRUE)


res_df <- res_df[, c("GeneID", "GeneName", "chr", "start", "end", "strand",
                     "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

# Save res_df here for null model

res_df <- res_df[order(res_df$padj, na.last=TRUE), ]
write.table(res_df, sub("\\.tsv$", "_all_results.tsv", args$output), sep="\t", quote=FALSE, row.names=FALSE)

sig <- subset(res_df, abs(log2FoldChange) >= args$lfc & padj < args$padj)
cat(sprintf("[INFO] Total genes tested: %d\n", nrow(res_df)))
cat(sprintf("[INFO] Significant DEGs (|LFC| >= %.1f, padj < %.2f): %d\n",
            args$lfc, args$padj, nrow(sig)))

write.table(sig, args$output, sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("[INFO] Results written to %s\n", args$output))

