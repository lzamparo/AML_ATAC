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

setwd("/Users/zamparol/projects/AML_ATAC/results/DESeq")
# load the experimental definitions
translator <- read_excel("../../data/filename_to_sample.xlsx")
col_data <- data.table(translator)
rm(translator)

col_data[,Batch := factor(Batch, levels = c(1,2))] 
col_data[, Condition := factor(Condition, levels = c("P", "A", "S", "SA", "SAR", "SARF", "SARN"))]
col_data[, Adapter := factor(c(rep("Good",10),"Bad", rep("Good",10)), levels=c("Good","Bad"))]  # Including adapter data from Tiansu
rownames <- col_data[, paste(Condition, Replicate, sep="_")]
rownames(col_data) <- rownames

# load the DESeq objects
dds_peaks = readRDS("objects/dds_object.rds")
res_A_P = readRDS("objects/res_A_P.rds")
res_SA_A = readRDS("objects/res_SA_A.rds")
res_SAR_SA = readRDS("objects/res_SAR_SA.rds")
res_S_P = readRDS("objects/res_S_P.rds")
res_SA_S = readRDS("objects/res_SA_S.rds")

res_A_P_dt = as.data.table(res_A_P)
res_S_P_dt = as.data.table(res_S_P)
res_SA_A_dt = as.data.table(res_SA_A)
res_SA_S_dt = as.data.table(res_SA_S)
res_SAR_SA_dt = as.data.table(res_SAR_SA)

# load the list of annotated peaks
atlas = read.csv("../peaks/annotated_atlas_reduced.csv")
atlas_dt = data.table(atlas)
rm(atlas)
setnames(atlas_dt, "seqnames", "chrom")
setnames(atlas_dt, "annot", "annotation")

sig_alpha <- 0.05

# Include SARN & SARF in the comparisons
# SARN vs SAR
res_SARN_SAR <- results(dds_peaks, contrast=c("Condition", "SARN", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARN_SAR)
plotMA(res_SARN_SAR, main="MA plot for gains, losses in SARN vs SAR", ylim=c(-2,2))
idx <- res_SARN_SAR$padj < 0.05
SARN_SAR_diff_peaks <- rownames(res_SARN_SAR)[idx]

# SARF vs SAR
res_SARF_SAR <- results(dds_peaks, contrast=c("Condition", "SARF", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SARF_SAR)
plotMA(res_SARF_SAR, main="MA plot for gains, losses in SARF vs SAR", ylim=c(-2,2))
idx <- res_SARF_SAR$padj < 0.05
SARF_SAR_diff_peaks <- rownames(res_SARF_SAR)[idx]

# Take the union of all peaks

# collect those list of peaks that are differential among the P-A, A-SA, SA-SAR, SAR-SARN, SAR-SARF pairwise comparisons
union_diff_peaks_A_first <- union(union(A_P_diff_peaks, SA_A_diff_peaks), union(union(SAR_SA_diff_peaks, SARN_SAR_diff_peaks), union(SARF_SAR_diff_peaks, SAR_SA_diff_peaks)))
bed_matrix <- str_split_fixed(union_diff_peaks_A_first, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../../results/peaks/P_A_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those list of peaks that are differential among the P-S, S-SA, SA-SAR, SAR-SARN, SAR-SARF pairwise comparisons
union_diff_peaks_S_first <- union(union(S_P_diff_peaks, SA_S_diff_peaks), union(union(SAR_SA_diff_peaks, SARN_SAR_diff_peaks), union(SARF_SAR_diff_peaks, SAR_SA_diff_peaks)))
bed_matrix <- str_split_fixed(union_diff_peaks_S_first, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "../../results/peaks/P_S_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
