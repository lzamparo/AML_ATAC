require(stats)
require(data.table)
require(tidyverse)
require(readxl)

### Read in a matrix of counts from the raw data.  Take the subsets we want to cluster (stage-wise differenes, and SAR vs P differential peaks)
### then get the clusters, cutree, send to computeMatrix / plotHeatmap for figures.

# load the data

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

# load the peak subsets to cluster
SAR_vs_P_peaks_bed = read.table("../../results/peaks/SAR_vs_P_diff_peaks.bed", sep="\t")
colnames(SAR_vs_P_peaks_bed) = c("chrom", "start", "end")
SAR_P_diff_peaks = paste(SAR_vs_P_peaks_bed$chrom,SAR_vs_P_peaks_bed$start,SAR_vs_P_peaks_bed$end,sep="-")
SAR_P_diff_matrix = count_matrix[SAR_P_diff_peaks,]

P_A_SA_SAR_SARN_SARF_diff_peaks_bed = read.table("../../results/peaks/P_A_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t")
colnames(P_A_SA_SAR_SARN_SARF_diff_peaks_bed) = c("chrom","start", "end")
P_A_SA_SAR_SARN_SARF_diff_peaks = paste(P_A_SA_SAR_SARN_SARF_diff_peaks_bed$chrom,P_A_SA_SAR_SARN_SARF_diff_peaks_bed$start,P_A_SA_SAR_SARN_SARF_diff_peaks_bed$end,sep="-")
P_A_SA_SAR_SARN_SARF_diff_matrix = count_matrix[P_A_SA_SAR_SARN_SARF_diff_peaks,]

P_S_SA_SAR_SARN_SARF_diff_peaks_bed = read.table("../../results/peaks/P_S_SA_SAR_SARN_SARF_diff_peaks.bed", sep="\t")
colnames(P_S_SA_SAR_SARN_SARF_diff_peaks_bed) = c("chrom","start", "end")
P_S_SA_SAR_SARN_SARF_diff_peaks = paste(P_S_SA_SAR_SARN_SARF_diff_peaks_bed$chrom,P_S_SA_SAR_SARN_SARF_diff_peaks_bed$start,P_S_SA_SAR_SARN_SARF_diff_peaks_bed$end,sep="-")
P_S_SA_SAR_SARN_SARF_diff_matrix = count_matrix[P_S_SA_SAR_SARN_SARF_diff_peaks,]

### calculate the distance objects, cluster, cut into distinct groups

# first, try 1 - spearman correlation as distance
# 1 - <.> because correlation is a similarity measure and we want distance
spearman_SAR_P <- cor(t(SAR_P_diff_matrix), method="spearman")
spearman_dist_SAR_P <- as.dist(1-spearman_SAR_P)  

spearman_A_first <- cor(t(P_A_SA_SAR_SARN_SARF_diff_matrix), method="spearman")
spearman_dist_A_first <- as.dist(1-spearman_A_first)

spearman_S_first <- cor(t(P_S_SA_SAR_SARN_SARF_diff_matrix), method="spearman")
spearman_dist_S_first <- as.dist(1-spearman_S_first)

hc_SAR_P_spearman <- hclust(spearman_dist_SAR_P, method="ward.D")
hc_A_first_spearman <- hclust(spearman_dist_A_first, method="ward.D")
hc_S_first_spearman <- hclust(spearman_dist_S_first, method="ward.D")

plot(hc_SAR_P_spearman, labels=FALSE)  # looks like k ~ 7 
plot(hc_A_first_spearman, labels=FALSE) # looks lke k ~ 6
plot(hc_S_first_spearman, labels=FALSE) # looks like k ~ 6

SAR_P_spearman_clusters <- cutree(hc_SAR_P_spearman, k = 7)
A_first_spearman_clusters <- cutree(hc_A_first_spearman, k = 7)
S_first_spearman_clusters <- cutree(hc_S_first_spearman, k = 7)

rm(list = c("spearman_SAR_P", "spearman_S_first", "spearman_A_first", "spearman_dist_SAR_P", "spearman_dist_S_first", "spearman_dist_A_first"))
gc()

# try canberra as distance
canberra_SAR_P <- dist(SAR_P_diff_matrix, method="canberra")
canberra_A_first <- dist(P_A_SA_SAR_SARN_SARF_diff_matrix, method="canberra")
canberra_S_first <- dist(P_S_SA_SAR_SARN_SARF_diff_matrix, method="canberra")

hc_SAR_P_canberra <- hclust(canberra_SAR_P, method="ward.D")
hc_A_first_canberra <- hclust(canberra_A_first, method="ward.D")
hc_S_first_canberra <- hclust(canberra_S_first, method="ward.D")

plot(as.dendrogram(hc_SAR_P_canberra), labels=FALSE, ylim=c(0,5000))  # looks like k ~ 10 but this is tough to gauge
plot(as.dendrogram(hc_A_first_canberra), labels=FALSE, ylim=c(0,5000)) # looks like k ~ 7
plot(as.dendrogram(hc_S_first_canberra), labels=FALSE, ylim=c(0,5000)) # looks like k ~ 7

SAR_P_canberra_clusters <- cutree(hc_SAR_P_canberra, k = 7)
A_first_canberra_clusters <- cutree(hc_A_first_canberra, k = 7)
S_first_canberra_clusters <- cutree(hc_S_first_canberra, k = 7)

# apply the cut to the bedfiles.  I don't have a great way to evaluate Spearman vs Canberra
# right now, but I'll go with spearman on the grounds that it should capture broad changes in ordering
# without being sensitive to absolute values of the counts
SAR_vs_P_peaks_bed$cluster_ID <- SAR_P_spearman_clusters
P_A_SA_SAR_SARN_SARF_diff_peaks_bed$cluster_ID <- A_first_spearman_clusters
P_S_SA_SAR_SARN_SARF_diff_peaks_bed$cluster_ID <- S_first_spearman_clusters

# write out the SAR_vs_P, A_first, S_first clusters to bedfiles
SAR_vs_P_prefix = "SAR_vs_P_"
A_prefix = "P_A_SA_SAR_SARN_SARF_"
S_prefix = "P_S_SA_SAR_SARN_SARF_"

write_section  = function(DF, prefix) {
  write.table(DF[, c("chrom","start", "end")], paste0(prefix,unique(DF$cluster_ID),".bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
  return(DF)
}

setwd("../../results/peaks/tornado_cluster_bedfiles/SAR_vs_P")
SAR_vs_P_peaks_bed %>% group_by(cluster_ID) %>% do(write_section(.,SAR_vs_P_prefix))
setwd("../A_first")
P_A_SA_SAR_SARN_SARF_diff_peaks_bed %>% group_by(cluster_ID) %>% do(write_section(.,A_prefix))
setwd("../S_first")
P_S_SA_SAR_SARN_SARF_diff_peaks_bed %>% group_by(cluster_ID) %>% do(write_section(.,S_prefix))

