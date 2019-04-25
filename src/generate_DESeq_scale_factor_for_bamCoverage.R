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
drop_list = c("SRSF2_c1_8", "ASXL1_34_1_repeat", "SfA_40")
files = data.table(files)
files = files[!(Filename %in% drop_list),]

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

# drop S3, A1, SA3
col_data = col_data[!(Filename %in% drop_list),]

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
dds_peaks <- DESeq(dds, parallel=TRUE)

### Save the sizeFactors for making .bw tracks, but save the invese of each, 
### since DESeq scales by *dividing* by the calculated factor, but deepTools
### bamCoverage scales by *multiplying* by the provided value (cf. https://www.biostars.org/p/285318/#362709)
size_factors <- sizeFactors(dds_peaks)
scale_factors <- size_factors ^ (-1)
scale_factors_df <- data.frame(Filename=reordered_col_data[,Filename], scale_factor=scale_factors)

setwd('/Users/zamparol/projects/AML_ATAC/results')
write.csv(scale_factors_df, file="scale_factors_deeptools.csv", quote=FALSE, row.names=FALSE)
