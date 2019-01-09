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

rownames <- get_peak_names(files[1,])
rownames(count_matrix) <- rownames

# prepare the col_data object 
translator <- read_excel("~/projects/AML_ATAC/data/filename_to_sample.xlsx")
col_data <- data.table(translator)
 
# ensure the batch is taken as a factor
col_data[,Batch := factor(Batch, levels = c(1,2))] 
col_data[, Condition := factor(Condition, levels = c("P", "A", "S", "SA", "SAR", "SARF", "SARN"))]
col_data[, Adapter := factor(c(rep("Good",10),"Bad", rep("Good",10)), levels=c("Good","Bad"))]  # Including adapter data from Tiansu
rownames <- col_data[, paste(Condition, Replicate, sep="_")]
rownames(col_data) <- rownames

# set colnames of the count_matrix to same
colnames(count_matrix) <- rownames

# Set the DDS object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ Batch + Condition + Adapter)

# load the flanks dds, transfer the estimates size factors
#flanks_dds <- readRDS("../../results/flank_dds.rds")
#sizeFactors(dds) <- sizeFactors(flanks_dds)
dds_peaks <- DESeq(dds, parallel=TRUE)
sig_alpha <- 0.05


### Get the results for each treatment condition versus P

# A vs P
res_A_P <- results(dds_peaks, contrast=c("Condition", "A", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_A_P)
idx <- res_A_P$padj < 0.05
A_P_diff_peaks <- rownames(res_A_P)[idx]

# SA vs P
res_SA_P<- results(dds_peaks, contrast=c("Condition", "SA", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_P)
idx <- res_SA_P$padj < 0.05
SA_P_diff_peaks <- rownames(res_SA_P)[idx]

# SAR vs P
res_SAR_P <- results(dds_peaks, contrast=c("Condition", "SAR", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SAR_P)
idx <- res_SAR_P$padj < 0.05
SAR_P_diff_peaks <- rownames(res_SAR_P)[idx]

# SARF vs P
res_SARF_P <- results(dds_peaks, contrast=c("Condition", "SARF", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARF_P)
idx <- res_SARF_P$padj < 0.05
SARF_P_diff_peaks <- rownames(res_SARF_P)[idx]

# SARN vs P
res_SARN_P <- results(dds_peaks, contrast=c("Condition", "SARN", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARN_P)
idx <- res_SARN_P$padj < 0.05
SARN_P_diff_peaks <- rownames(res_SARN_P)[idx]

### Now look at S vs P, SA vs S, SA vs A, SAR vs SA
res_S_P <- results(dds_peaks, contrast=c("Condition", "S", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_S_P)
idx <- res_S_P$padj < 0.05
S_P_diff_peaks <- rownames(res_S_P)[idx]

res_SA_S <- results(dds_peaks, contrast=c("Condition", "SA", "S"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_S)
idx <- res_SA_S$padj < 0.05
SA_S_diff_peaks <- rownames(res_SA_S)[idx]

res_SA_A <- results(dds_peaks, contrast=c("Condition", "SA", "A"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_A)
idx <- res_SA_A$padj < 0.05
SA_A_diff_peaks <- rownames(res_SA_A)[idx]

res_SAR_SA <- results(dds_peaks, contrast=c("Condition", "SAR", "SA"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SAR_SA)
idx <- res_SAR_SA$padj < 0.05
SAR_SA_diff_peaks <- rownames(res_SAR_SA)[idx]

# Include SARN & SARF in the comparisons
# SARN vs SAR
res_SARN_SAR <- results(dds_peaks, contrast=c("Condition", "SARN", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARN_SAR)
idx <- res_SARN_SAR$padj < 0.05
SARN_SAR_diff_peaks <- rownames(res_SARN_SAR)[idx]

# SARF vs SAR
res_SARF_SAR <- results(dds_peaks, contrast=c("Condition", "SARF", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARF_SAR)
idx <- res_SARF_SAR$padj < 0.05
SARF_SAR_diff_peaks <- rownames(res_SARF_SAR)[idx]

setwd('../../results/DESeq')

# collect those peaks that are differential among the A-P, SA-P, SAR-P, SARN-P, SARF-P pairwise comparisons
union_diff_peaks_conditions_vs_P <- union(union(A_P_diff_peaks, SA_P_diff_peaks), union(union(SAR_P_diff_peaks, SARF_P_diff_peaks),SARN_P_diff_peaks))
bed_matrix <- str_split_fixed(union_diff_peaks_conditions_vs_P, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/all_conditions_vs_P_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those peaks that are differential among the A-P, SA-A, SAR-SA, SARN-SAR, SARF-SAR pairwise comparisons
union_diff_peaks_A_first <- union(union(A_P_diff_peaks, SA_A_diff_peaks), union(union(SAR_SA_diff_peaks, SARN_SAR_diff_peaks), union(SARF_SAR_diff_peaks, SAR_SA_diff_peaks)))
bed_matrix <- str_split_fixed(union_diff_peaks_A_first, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/P_A_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those peaks that are differential among the S-P, SA-S, SAR-SA, SARN-SAR, SARF-SAR pairwise comparisons
union_diff_peaks_S_first <- union(union(S_P_diff_peaks, SA_S_diff_peaks), union(union(SAR_SA_diff_peaks, SARN_SAR_diff_peaks), union(SARF_SAR_diff_peaks, SAR_SA_diff_peaks)))
bed_matrix <- str_split_fixed(union_diff_peaks_S_first, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../peaks/P_S_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those peaks that are differential among any of the comparisons.  Rainbow product FTW.
all_conditions <- t(combn(unique(col_data[,Condition]),2))
first_args = c(as.character(all_conditions[,2]))
second_args = c(as.character(all_conditions[,1]))
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
