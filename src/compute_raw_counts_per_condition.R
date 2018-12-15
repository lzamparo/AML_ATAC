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
translator <- read_excel("../filename_to_sample.xlsx")
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
flanks_dds <- readRDS("../../results/flank_dds.rds")
sizeFactors(dds) <- sizeFactors(flanks_dds)

dds_peaks <- DESeq(dds, parallel=TRUE)

# Check that we're getting replicates clustering together
rld <- rlog(dds_peaks, blind=FALSE)

# PCA of conditions where batch effects are removed
mat <- assay(rld)
mat_nobatch <- limma::removeBatchEffect(mat, design=model.matrix(~ col_data$Condition), batch=col_data$Batch, batch2=col_data$Adapter)
assay(rld) <- mat_nobatch
data <- plotPCA(rld, intgroup=c("Condition"), returnData=TRUE)  ### pick ntop based on number of up & down peaks.
percentVar <- round(100 * attr(data, "percentVar"))

# For ASH: need just P-A-SA-SAR conditions...
ggplot(subset(data, data$Condition %in% c("WT","ASXL1","SRSF2-ASXL1","SRSF2-ASXL1-NRAS")), aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  ggtitle("PCA of log transformed read counts") +
  geom_text_repel(
    aes(label = name),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) 


### Get the results for each treatment condition versus WT along the major 
### expected differentiation path

sig_alpha <- 0.05

# A vs P
res_A_P <- results(dds_peaks, contrast=c("Condition", "P", "A"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_A_P)
plotMA(res_A_P, main="MA plot for gains, losses in A vs P", ylim=c(-2,2))
idx <- res_A_P$padj < 0.05
A_P_diff_peaks <- rownames(res_A_P)[idx]

# SA vs P
res_SA_P<- results(dds_peaks, contrast=c("Condition", "P", "SA"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_P)
plotMA(res_SA_P, main="MA plot for gains, losses in SA vs P", ylim=c(-2,2))
idx <- res_SA_P$padj < 0.05
SA_P_diff_peaks <- rownames(res_SA_P)[idx]

# SAR vs P
res_SAR_P <- results(dds_peaks, contrast=c("Condition", "P", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SAR_P)
plotMA(res_SAR_P, main="MA plot for gains, losses in SAR vs P", ylim=c(-2,2))
idx <- res_SAR_P$padj < 0.05
SAR_P_diff_peaks <- rownames(res_SAR_P)[idx]

### Now look at P vs S, AS vs S, SA vs A, SAR vs SA
res_S_P <- results(dds_peaks, contrast=c("Condition", "P", "S"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_S_P)
plotMA(res_S_P, main="MA plot for gains, losses in S vs P", ylim=c(-2,2))
idx <- res_S_P$padj < 0.05
S_P_diff_peaks <- rownames(res_S_P)[idx]

res_SA_S <- results(dds_peaks, contrast=c("Condition", "S", "SA"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_S)
plotMA(res_SA_S, main="MA plot for gains, losses in SA vs S", ylim=c(-2,2))
idx <- res_SA_S$padj < 0.05
SA_S_diff_peaks <- rownames(res_SA_S)[idx]

res_SA_A <- results(dds_peaks, contrast=c("Condition", "A", "SA"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SA_A)
plotMA(res_SA_A, main="MA plot for gains, losses in SA vs A", ylim=c(-2,2))
idx <- res_SA_A$padj < 0.05
SA_A_diff_peaks <- rownames(res_SA_A)[idx]

res_SAR_SA <- results(dds_peaks, contrast=c("Condition", "SA", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
summary(res_SAR_SA)
plotMA(res_SAR_SA, main="MA plot for gains, losses in SAR vs SA", ylim=c(-2,2))
idx <- res_SAR_SA$padj < 0.05
SAR_SA_diff_peaks <- rownames(res_SAR_SA)[idx]

setwd('../results/DESeq')

# collect those list of peaks that are differential among the P-A, A-SA, SA-SAR pairwise comparisons
union_diff_peaks_main_trajectory <- union(union(A_WT_diff_peaks, SA_A_diff_peaks), SA_SAR_diff_peaks)
bed_matrix <- str_split_fixed(union_diff_peaks_main_trajectory, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "peak_lists/WT_A_SA_SAR_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# collect those list of peaks that are differential among the P-A, P-SA, P-SAR pairwise comparisons
union_diff_peaks_conditions_vs_WT <- union(union(A_WT_diff_peaks, SA_WT_diff_peaks), SAR_WT_diff_peaks)
bed_matrix <- str_split_fixed(union_diff_peaks_conditions_vs_WT, "-", 3)
bed_df <- as.data.frame(bed_matrix)
write.table(bed_df, "peak_lists/all_conditions_vs_WT_diff_peaks.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# save this set, plus the results objets for later retrieval.  
setwd('../results/DESeq')
saveRDS(dds_peaks, file = "dds_object.rds")
saveRDS(res_A_WT, file = "res_A_WT.rds")
saveRDS(res_S_WT, file = "res_S_WT.rds")
saveRDS(res_SA_WT, file = "res_SA_WT.rds")
saveRDS(res_SAR_WT, file = "res_SAR_WT.rds")

saveRDS(res_SA_A, file = "res_SA_A.rds")
saveRDS(res_SA_S, file = "res_SA_S.rds")
saveRDS(res_SA_SAR, file = "res_SA_SAR.rds")
