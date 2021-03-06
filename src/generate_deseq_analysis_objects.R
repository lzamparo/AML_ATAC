require(data.table)
require(DESeq2)
require(plyr)
require(stringr)
require(Biobase)
require(GenomicRanges)
require(BiocParallel)
require(ggplot2)
require(readxl)
require(ggrepel)

register(MulticoreParam(4))
### compute median coverage for each peak in candidate set & experiment, to provide DESeq2 style raw-counts data.

# go to counts dir
setwd('/Users/zamparol/projects/AML_ATAC/data/counts')

# load all timepoints data into 
conditions <- list.files(path = ".", include.dirs = TRUE)
reps <- list.files(path = ".", recursive=TRUE, pattern = "*trimmed_001_read_counts.bed")
files <- ldply(strsplit(reps, '/', fixed = TRUE))
colnames(files) <- c("condition","replicate")
files$Filename = str_split(files$replicate, "_S[0-9]{1,2}", simplify=TRUE)[,1]

# Drop S3; it does not contain the SRSF2 mutation :/
# Also drop: A3 (ASXL1_34), SA1 (SfA_40), they show strange patterns from RNA-seq analyses
files = data.table(files)
files = files[!(Filename %in% c("SRSF2_c1_8", "ASXL1_34_1_repeat", "SfA_40")),]


# load data by condition & replicate into a count matrix
get_counts <- function(myrow){
  my_rep <- read.delim(file = paste(myrow[1], myrow[2], sep = "/"), header = FALSE)
  colnames(my_rep) <- c("chrom","start", "end", "count")
  my_rep <- data.table(my_rep)
  my_rep[chrom != "all", count]
}
count_matrix <- as.matrix(cbind(apply(files, 1, get_counts)))

# get the rownames for the count matrix
get_peak_names <- function(myrow){
  my_rep <- read.delim(file = paste(myrow$condition, myrow$replicate, sep = "/"), header = FALSE)
  colnames(my_rep) <- c("chrom","start", "end", "count")
  my_rep <- data.table(my_rep)
  my_rep[chrom != "all", paste(chrom,start,end, sep = "-")]
}

peaks <- get_peak_names(files[1,])
rownames(count_matrix) <- peaks

# prepare the col_data object 
translator <- read_excel("~/projects/AML_ATAC/data/filename_to_sample.xlsx")
col_data <- data.table(translator)

 
# ensure the batch is taken as a factor
col_data[, Batch := factor(Batch, levels = c(1,2))] 
col_data[, Condition := factor(Condition, levels = c("P", "A", "S", "SA", "SAR", "SARF", "SARN"))]
col_data[, Adapter := factor(c(rep("Good",10),"Bad", rep("Good",10)), levels=c("Good","Bad"))]  # Including adapter data from Tiansu

# drop S3 
col_data = col_data[!(Filename %in% c("SRSF2_c1_8", "ASXL1_34_1_repeat", "SfA_40")),]

# set colnames of the count_matrix to same
# *** ENSURE THAT THE COLNAMES OF THE COUNT MATRIX CORRESPOND WITH THE ORDER IN WHICH THE COLUMNS OF THE MATRIX WERE READ IN: i.e `files` ***
setkey(col_data, Filename)
reordered_col_data <- col_data[files$Filename, ]
rownames <- reordered_col_data[, paste(Condition, Replicate, sep="_")]
rownames(reordered_col_data) <- rownames
colnames(count_matrix) <- rownames

# Sanity check
if (!all(rownames(reordered_col_data) == colnames(count_matrix))){
  print("Mismatch in count matrix column names, rownames of reordered_col_data")
  stop()
}

# Set the DDS object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = reordered_col_data,
                              design = ~ Batch + Condition)

# load the flanks dds, transfer the estimates size factors
#flanks_dds <- readRDS("../../results/flank_dds.rds")
#sizeFactors(dds) <- sizeFactors(flanks_dds)
dds_peaks <- DESeq(dds, parallel=TRUE)
sig_alpha <- 0.05

### Get the results for each treatment condition versus P

# A vs P
res_A_P <- results(dds_peaks, contrast=c("Condition", "A", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_A_P)
idx <- res_A_P$padj < sig_alpha
A_P_diff_peaks <- rownames(res_A_P)[idx]

res_S_P <- results(dds_peaks, contrast=c("Condition", "S", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_S_P)
idx <- res_S_P$padj < sig_alpha
S_P_diff_peaks <- rownames(res_S_P)[idx]

# SA vs P
res_SA_P<- results(dds_peaks, contrast=c("Condition", "SA", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_P)
idx <- res_SA_P$padj < sig_alpha
SA_P_diff_peaks <- rownames(res_SA_P)[idx]

# SAR vs P
res_SAR_P <- results(dds_peaks, contrast=c("Condition", "SAR", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SAR_P)
idx <- res_SAR_P$padj < sig_alpha
SAR_P_diff_peaks <- rownames(res_SAR_P)[idx]

# SARF vs P
res_SARF_P <- results(dds_peaks, contrast=c("Condition", "SARF", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARF_P)
idx <- res_SARF_P$padj < sig_alpha
SARF_P_diff_peaks <- rownames(res_SARF_P)[idx]

# SARN vs P
res_SARN_P <- results(dds_peaks, contrast=c("Condition", "SARN", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARN_P)
idx <- res_SARN_P$padj < sig_alpha
SARN_P_diff_peaks <- rownames(res_SARN_P)[idx]

### Now look at SA vs S, SA vs A, SAR vs SA, SARN vs SAR , SARF vs SAR
res_SA_S <- results(dds_peaks, contrast=c("Condition", "SA", "S"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_S)
idx <- res_SA_S$padj < sig_alpha
SA_S_diff_peaks <- rownames(res_SA_S)[idx]

res_SA_A <- results(dds_peaks, contrast=c("Condition", "SA", "A"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_A)
idx <- res_SA_A$padj < sig_alpha
SA_A_diff_peaks <- rownames(res_SA_A)[idx]

res_SAR_SA <- results(dds_peaks, contrast=c("Condition", "SAR", "SA"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SAR_SA)
idx <- res_SAR_SA$padj < sig_alpha
SAR_SA_diff_peaks <- rownames(res_SAR_SA)[idx]

res_SARN_SAR <- results(dds_peaks, contrast=c("Condition", "SARN", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARN_SAR)
idx <- res_SARN_SAR$padj < sig_alpha
SARN_SAR_diff_peaks <- rownames(res_SARN_SAR)[idx]

res_SARF_SAR <- results(dds_peaks, contrast=c("Condition", "SARF", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARF_SAR)
idx <- res_SARF_SAR$padj < sig_alpha
SARF_SAR_diff_peaks <- rownames(res_SARF_SAR)[idx]

setwd('../../results/DESeq')

# write out the SAR vs P diff peaks
bed_matrix <- str_split_fixed(SAR_P_diff_peaks, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/SAR_vs_P_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those peaks that are differential among the A-P, SA-P, SAR-P, SARN-P, SARF-P pairwise comparisons
union_diff_peaks_conditions_vs_P <- unique(do.call(c, list(A_P_diff_peaks, SA_P_diff_peaks,SAR_P_diff_peaks, SARF_P_diff_peaks,SARN_P_diff_peaks)))
bed_matrix <- str_split_fixed(union_diff_peaks_conditions_vs_P, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/all_conditions_vs_P_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those peaks that are differential among the A-P, SA-A, SAR-SA, SARN-SAR, SARF-SAR pairwise comparisons
union_diff_peaks_A_first <- unique(do.call(c, list(A_P_diff_peaks, SA_A_diff_peaks, SAR_SA_diff_peaks, SARN_SAR_diff_peaks, SARF_SAR_diff_peaks)))
bed_matrix <- str_split_fixed(union_diff_peaks_A_first, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/P_A_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those peaks that are differential among the S-P, SA-S, SAR-SA, SARN-SAR, SARF-SAR pairwise comparisons
union_diff_peaks_S_first <- unique(do.call(c, list(S_P_diff_peaks, SA_S_diff_peaks, SAR_SA_diff_peaks, SARN_SAR_diff_peaks, SARF_SAR_diff_peaks)))
bed_matrix <- str_split_fixed(union_diff_peaks_S_first, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/P_S_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those peaks that are differential among any of the comparisons.  Rainbow product FTW.
# this is broken: rows [2,21] corresponding to A, S and SARF, SARN do not belong.  Also row 1 (A,P) should be inverted
# no longer includes: A vs S, SARF vs SARN
all_conditions <- t(combn(unique(reordered_col_data[,Condition]),2))
first_args = c(as.character(all_conditions[1,1]), as.character(all_conditions[3:20,2]))
second_args = c(as.character(all_conditions[1,2]), as.character(all_conditions[3:20,1]))
get_diff_peaks <- function(first, second){
  res_contrast <- results(dds_peaks, contrast=c("Condition", first, second), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
  idx <- res_contrast$padj < 0.05
  rownames(res_contrast)[idx]
}
diff_list = mapply(get_diff_peaks, first_args, second_args)
diff_list_set = unique(do.call(c, diff_list)) 
bed_matrix <- str_split_fixed(diff_list_set, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/all_pairwise_comparisons_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# save this set, plus the results objets for later retrieval.  
setwd('./objects/')
saveRDS(dds_peaks, file = "dds_object.rds")
saveRDS(res_A_P, file = "res_A_P.rds")
saveRDS(res_S_P, file = "res_S_P.rds")
saveRDS(res_SA_P, file = "res_SA_P.rds")
saveRDS(res_SAR_P, file = "res_SAR_P.rds")
saveRDS(res_SARF_P, file = "res_SARF_P.rds")
saveRDS(res_SARN_P, file = "res_SARN_P.rds")

saveRDS(res_SA_A, file = "res_SA_A.rds")
saveRDS(res_SA_S, file = "res_SA_S.rds")
saveRDS(res_SAR_SA, file = "res_SAR_SA.rds")
saveRDS(res_SARN_SAR, file = "res_SARN_SAR.rds")
saveRDS(res_SARF_SAR, file = "res_SARF_SAR.rds")
