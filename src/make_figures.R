require(data.table)
require(ggplot2)
require(scales)
require(ggrepel)
require(cowplot)

# plot pie chart of atlas, differentiable peaks atlas

# load atlas file
setwd("/Users/zamparol/projects/AML_ATAC/peaks")
atlas = read.csv("annotated_atlas_reduced.csv")
atlas_dt = data.table(atlas)
setnames(atlas_dt, "seqnames", "chrom")
setnames(atlas_dt, "annot", "annotation")

total = atlas_dt[,.N]
grouped_peaks = atlas_dt[,.N, by=annotation]
grouped_peaks[, pos := cumsum(N) - N / 4 ]
grouped_peaks[, percentage := N / total]
#grouped_peaks[annotation %in% c("intron", "promoter", "intergenic"), pos := pos - 6000]

# pie chart of atlas
atlas_pie = ggplot(atlas_dt, aes(x=factor(1), fill=annotation)) + 
  geom_bar(width=1) + 
  coord_polar(theta = "y") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        axis.title=element_blank()) +
  ggtitle("All peaks in atlas (by annotation)") + 
  geom_text_repel(data=grouped_peaks, aes(x=factor(1), y=pos, label=percent(percentage)), vjust= 1.5, size=6)

# pie chart of diff'ble peaks
setwd("/Users/zamparol/projects/AML_ATAC/results/DESeq/peak_lists")
diff_peaks = read.delim(file="WT_A_SA_SAR_diff_peaks.bed", sep="\t")
colnames(diff_peaks) = c("chrom", "start", "end")
diff_peaks_dt = data.table(diff_peaks)
key_cols = c("chrom", "start", "end")
setkeyv(diff_peaks_dt, key_cols)
setkeyv(atlas_dt, key_cols)
annotated_diff_peaks = atlas_dt[diff_peaks_dt,]

# pie chart of differential peak atlas
total = annotated_diff_peaks[,.N]
grouped_diff_peaks = annotated_diff_peaks[,.N, by=annotation]
grouped_diff_peaks[, pos := cumsum(N) - N / 4 ]
# swap percentage labels for intron, intergenic.  Pie chart displays them incorrectly.
grouped_diff_peaks[, percentage := N / total]
temp_intron = grouped_diff_peaks[annotation == "intron", percentage]
temp_intergenic = grouped_diff_peaks[annotation == "intergenic", percentage]
grouped_diff_peaks[annotation == "intron", percentage := temp_intergenic]
grouped_diff_peaks[annotation == "intergenic", percentage := temp_intron]

diff_pie = ggplot(annotated_diff_peaks, aes(x=factor(1), fill=annotation)) + 
  geom_bar(width=1) + 
  coord_polar(theta = "y") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        axis.title=element_blank()) +
  ggtitle("Differential peaks along P-A-SA-SAR (by annotation)") + 
  geom_text_repel(data=grouped_diff_peaks, aes(x=factor(1), y=pos, label=percent(percentage)), size=6)


# Bar chart of opening / closing peaks in each transition (P -> A -> SA -> SAR)
setwd('/Users/zamparol/projects/AML_ATAC/results/DESeq/objects')
dds_peaks = readRDS("dds_object.rds")
res_WT_A = readRDS("res_A_WT.rds")
res_A_SA = readRDS("res_SA_A.rds")
res_SA_SAR = readRDS("res_SA_SAR.rds")
res_WT_S = readRDS("res_S_WT.rds")
res_S_SA = readRDS("res_SA_S.rds")

res_WT_A_dt = as.data.table(res_WT_A)
res_WT_S_dt = as.data.table(res_WT_S)
res_A_SA_dt = as.data.table(res_A_SA)
res_S_SA_dt = as.data.table(res_S_SA)
res_SA_SAR_dt = as.data.table(res_SA_SAR)

peaks_up = c(res_WT_A_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_WT_S_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_A_SA_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_S_SA_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SA_SAR_dt[padj < 0.05 & log2FoldChange > 0, .N ])
peaks_down = c(res_WT_A_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_WT_S_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_A_SA_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_S_SA_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SA_SAR_dt[padj < 0.05 & log2FoldChange < 0, .N ])
conditions = rep(c("A vs P", "S vs P", "SA vs A", "SA vs S", "SAR vs SA"),2)
amount = c(peaks_up, peaks_down)
change = c(rep("gained", 5),rep("lost",5))
df = data.frame(conditions, change, amount)
colnames(df) = c("condition", "change", "amount")
# correct the factor levels
df$condition = ordered(rep(c("A vs P", "S vs P", "SA vs A", "SA vs S", "SAR vs SA"),2), levels = c("A vs P", "S vs P", "SA vs A", "SA vs S", "SAR vs SA"))
rm(list =c("conditions", "amount", "change", "peaks_down", "peaks_up"))
df$path = c("ASXL truncation first", "SRSF2 P95L first", "ASXL truncation first", "SRSF2 P95L first", "common", "ASXL truncation first", "SRSF2 P95L first", "ASXL truncation first", "SRSF2 P95L first", "common")

peak_changes = ggplot(df, aes(x=condition, fill=change)) + 
  geom_bar(data = subset(df, change == "gained" & path %in% c("ASXL truncation first", "common")), aes(y = amount), position="stack", stat="identity") + 
  geom_bar(data = subset(df, change == "lost" & path %in% c("ASXL truncation first", "common")), aes(y = -amount), position="stack", stat="identity") +
  ggtitle("Peaks gained and lost in P-A-SA-SAR transitions") + ylab("Peaks") + xlab("Transition") + 
  scale_fill_manual(name = "Change", values=c("red", "blue"))

peak_changes_alt_path = ggplot(df, aes(x=condition, fill=change)) + 
  geom_bar(data = subset(df, change == "gained" & path %in% c("SRSF2 P95L first", "common")), aes(y = amount), position="stack", stat="identity") + 
  geom_bar(data = subset(df, change == "lost" & path %in% c("SRSF2 P95L first", "common")), aes(y = -amount), position="stack", stat="identity") +
  ggtitle("Peaks gained and lost in P-S-SA-SAR transitions") + ylab("Peaks") + xlab("Transition") + 
  scale_fill_manual(name = "Change", values=c("red", "blue"))

# What is up with the lack of changes for A to SA?
A_SA_hist = ggplot(as.data.frame(res_A_SA), aes(x=pvalue)) + geom_histogram(bins = 100) + xlab("P values") + ggtitle("P values for log2 fold change test for A -> SA transition")

# Load the clustered peaks from plotHeatmap, look at the annotation breakdown, and the genes represented
setwd("~/projects/AML_ATAC/peaks")
annotated_atlas <- read.csv(file="annotated_atlas_reduced.csv")
atlas_dt = data.table(annotated_atlas)
diff_peaks = read.delim(file="P-A-SA-SAR_clusters.bed", sep="\t", header=TRUE)
diff_peaks = diff_peaks[,c("X.chrom", "start","end","deepTools_group")]
diff_peaks_dt = data.table(diff_peaks)
setnames(diff_peaks_dt, "X.chrom", "chrom")
setnames(atlas_dt, "seqnames", "chrom")
setkeyv(atlas_dt, c("chrom", "start", "end"))
setkeyv(diff_peaks_dt, c("chrom", "start", "end"))
annotated_diff_peaks = atlas_dt[diff_peaks_dt,]

annotated_diff_peaks[, group_total := .N, by=deepTools_group]
annotated_diff_peaks[, group_percentage := .N / group_total, by=.(deepTools_group, annot)]
grouped_percentages = unique(annotated_diff_peaks[,.(annot, group_percentage, deepTools_group)])
ggplot(grouped_percentages, aes(x=annot, fill=annot)) + 
  geom_col(aes(y=percent(group_percentage))) + 
  facet_wrap(~ deepTools_group) + 
  xlab("Annotation") + 
  ylab("Percentage of peaks") + 
  ggtitle("Differential peak annotation breakdown by plotHeatmap cluster")
