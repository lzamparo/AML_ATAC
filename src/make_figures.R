require(data.table)
require(ggplot2)
require(scales)
require(ggrepel)
require(cowplot)

# plot pie chart of atlas, differentiable peaks atlas

# load atlas file
setwd("~/projects/AML_ATAC/peaks")
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
diff_peaks = read.delim(file="P_A_SA_SAR_diff_peaks.bed", sep="\t")
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


# Load the DESeq objects used to make differential plots
setwd('/Users/zamparol/projects/AML_ATAC/results/DESeq/objects')
dds_peaks = readRDS("dds_object.rds")
res_A_P = readRDS("res_A_P.rds")
res_SA_A = readRDS("res_SA_A.rds")
res_SAR_SA = readRDS("res_SAR_SA.rds")
res_S_P = readRDS("res_S_P.rds")
res_SA_S = readRDS("res_SA_S.rds")
res_SARF_SAR = readRDS("res_SARF_SAR.rds")
res_SARN_SAR = readRDS("res_SARN_SAR.rds")

res_A_P_dt = as.data.table(res_A_P)
res_S_P_dt = as.data.table(res_S_P)
res_SA_A_dt = as.data.table(res_SA_A)
res_SA_S_dt = as.data.table(res_SA_S)
res_SAR_SA_dt = as.data.table(res_SAR_SA)
res_SARN_SAR_dt = as.data.table(res_SARN_SAR)
res_SARF_SAR_dt = as.data.table(res_SARF_SAR)

### PCA plot
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


# MA plots versus P
plotMA(res_A_P, main="MA plot for gains, losses in A vs P", ylim=c(-2,2))
plotMA(res_S_P, main="MA plot for gains, losses in S vs P", ylim=c(-2,2))
plotMA(res_SA_P, main="MA plot for gains, losses in SA vs P", ylim=c(-2,2))
plotMA(res_SAR_P, main="MA plot for gains, losses in SAR vs P", ylim=c(-2,2))
plotMA(res_SARF_P, main="MA plot for gains, losses in SARF vs P", ylim=c(-2,2))
plotMA(res_SARN_P, main="MA plot for gains, losses in SARN vs P", ylim=c(-2,2))

# Remaining MA plots in (P -> A / S -> SA -> SAR -> SARN / SARF)
plotMA(res_SA_S, main="MA plot for gains, losses in SA vs S", ylim=c(-2,2))
plotMA(res_SA_A, main="MA plot for gains, losses in SA vs A", ylim=c(-2,2))
plotMA(res_SAR_SA, main="MA plot for gains, losses in SAR vs SA", ylim=c(-2,2))
plotMA(res_SARN_SAR, main="MA plot for gains, losses in SARN vs SAR", ylim=c(-2,2))
plotMA(res_SARF_SAR, main="MA plot for gains, losses in SARF vs SAR", ylim=c(-2,2))

# Bar chart of opening / closing peaks in each transition 
# (P -> A -> SA -> SAR-> SARN / SARF)
peaks_up = c(res_A_P_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SA_A_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SAR_SA_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SARN_SAR_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SARF_SAR_dt[padj < 0.05 & log2FoldChange > 0, .N ])
peaks_down = c(res_A_P_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SA_A_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SAR_SA_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SARN_SAR_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SARF_SAR_dt[padj < 0.05 & log2FoldChange < 0, .N ])
conditions = rep(c("A vs P", "SA vs A", "SAR vs SA", "SARN vs SAR", "SARF vs SAR"),2)
amount = c(peaks_up, peaks_down)
change = c(rep("gained", 5),rep("lost",5))
df_A = data.frame(conditions, change, amount)
df_A$path = "ASXL truncation first"
colnames(df_A) = c("condition", "change", "amount", "path")

# (P -> S -> SA -> SAR-> SARN / SARF)
peaks_up = c(res_S_P_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SA_S_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SAR_SA_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SARN_SAR_dt[padj < 0.05 & log2FoldChange > 0, .N ], res_SARF_SAR_dt[padj < 0.05 & log2FoldChange > 0, .N ])
peaks_down = c(res_S_P_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SA_S_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SAR_SA_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SARN_SAR_dt[padj < 0.05 & log2FoldChange < 0, .N ], res_SARF_SAR_dt[padj < 0.05 & log2FoldChange < 0, .N ])
conditions = rep(c("S vs P", "SA vs S", "SAR vs SA", "SARN vs SAR", "SARF vs SAR"),2)
amount = c(peaks_up, peaks_down)
change = c(rep("gained", 5),rep("lost",5))
df_S = data.frame(conditions, change, amount)
df_S$path = "SRSF2 P95L first"
colnames(df_S) = c("condition", "change", "amount", "path")

# correct the factor levels
rm(list =c("conditions", "amount", "change", "peaks_down", "peaks_up"))

peak_changes_A = ggplot(df_A, aes(x=condition, fill=change)) + 
  geom_bar(data = subset(df_A, change == "gained"), aes(y = amount), position="stack", stat="identity") + 
  geom_bar(data = subset(df_A, change == "lost"), aes(y = -amount), position="stack", stat="identity") +
  ggtitle("Peaks gained and lost in transitions", subtitle="ASXL truncation first") + ylab("Peaks") + xlab("Transition") + 
  scale_fill_manual(name = "Change", values=c("red", "blue"))

peak_changes_S = ggplot(df_S, aes(x=condition, fill=change)) + 
  geom_bar(data = subset(df_S, change == "gained"), aes(y = amount), position="stack", stat="identity") + 
  geom_bar(data = subset(df_S, change == "lost"), aes(y = -amount), position="stack", stat="identity") +
  ggtitle("Peaks gained and lost in transitions", subtitle="SRSF2 P95L first") + ylab("Peaks") + xlab("Transition") + 
  scale_fill_manual(name = "Change", values=c("red", "blue"))

# What is up with the lack of changes for A to SA?
A_SA_hist = ggplot(as.data.frame(res_A_SA), aes(x=pvalue)) + geom_histogram(bins = 100) + xlab("P values") + ggtitle("P values for log2 fold change test for A -> SA transition")

# Load the clustered peaks from plotHeatmap, look at the annotation breakdown, and the genes represented
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
