require(stats)
require(data.table)
require(plyr)
require(readxl)

### Read in a matrix of counts from the raw data.  Take the subsets we want to cluster (stage-wise differenes, and SAR vs P differential peaks)
### then get the clusters, cutree, send to computeMatrix / plotHeatmap for figures.

# load the data

# go to counts dir
setwd('/data/leslie/zamparol/AML_ATAC/data/counts')

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
translator <- read_excel("~/projects/AML_ATAC/data//filename_to_sample.xlsx")
col_data <- data.table(translator)

# ensure the batch is taken as a factor
col_data[,Batch := factor(Batch, levels = c(1,2))] 
col_data[, Condition := factor(Condition, levels = c("P", "A", "S", "SA", "SAR", "SARF", "SARN"))]
col_data[, Adapter := factor(c(rep("Good",10),"Bad", rep("Good",10)), levels=c("Good","Bad"))]  # Including adapter data from Tiansu
rownames <- col_data[, paste(Condition, Replicate, sep="_")]
rownames(col_data) <- rownames

# set colnames of the count_matrix to same
colnames(count_matrix) <- rownames

# load the peak subsets to cluster
all_pairwise_peaks_bed = read.table("/data/leslie/zamparol/AML_ATAC/data/results/differential_peak_lists/all_pairwise_comparisons_diff_peaks.bed", sep="\t")
colnames(all_pairwise_peaks_bed) = c("chrom", "start", "end")
all_pairwise_diff_peaks = paste(all_pairwise_peaks_bed$chrom, all_pairwise_peaks_bed$start, all_pairwise_peaks_bed$end, sep="-")
all_pairwise_diff_matrix = count_matrix[all_pairwise_diff_peaks,]

### calculate the distance objects, cluster, cut into distinct groups

spearman_all <- cor(t(all_pairwise_diff_matrix), method="spearman")
spearman_dist_all <- as.dist(1-spearman_all)
hc_all_spearman <- hclust(spearman_dist_all, method="ward.D")

# write out the hclust result to parse later
saveRDS(hc_all_spearman, "/data/leslie/zamparol/AML_ATAC/data/hc_all_spearman.rds")
saveRDS(spearman_dist_all, "/data/leslie/zamparol/AML_ATAC/data/spearman_dist_all.rds")


