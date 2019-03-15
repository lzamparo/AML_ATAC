require(data.table)
require(ggplot2)
require(scales)
require(ggrepel)
require(cowplot)
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(readxl)

# load atlas file
setwd("~/projects/AML_ATAC/results/peaks")
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
atlas_anno <- annotatePeak("all_conditions_peak_atlas.bed", tssRegion=c(-2000,2000), TxDb=txdb, level = "gene", annoDb="org.Hs.eg.db")
atlas = data.table(as.data.frame(atlas_anno))
rm(atlas_anno)
setnames(atlas, "seqnames", "chrom")
atlas[grepl("Exon",annotation), simple_annotation := "Exon"]
atlas[grepl("Intron",annotation), simple_annotation := "Intron"]
atlas[grepl("Promoter",annotation), simple_annotation := "Promoter"]
atlas[grepl("Intergenic",annotation), simple_annotation := "Intergenic"]
atlas[grepl("Downstream",annotation), simple_annotation := "Intergenic"]
atlas[grepl("5' UTR",annotation), simple_annotation := "5' UTR"]
atlas[grepl("3' UTR",annotation), simple_annotation := "3' UTR"]
setnames(atlas, c("annotation", "simple_annotation"), c("complex_annotation", "annotation"))

# Four panel plot of atlas peaks

# (a) Peaks by annotation
total = atlas[,.N]
grouped_peaks = atlas[,.N, by=annotation]
grouped_peaks[, percentage := N / total]
pba <- ggplot(grouped_peaks, aes(x=annotation, fill=annotation)) + 
  geom_col(aes(y=N)) + 
  geom_text(aes(x=annotation, y=N, label=percent(percentage)), vjust = -.50) +
  scale_y_continuous(breaks=c(0,10000,25000,50000,75000,100000,125000,150000)) +
  ggtitle("Number of peaks by annotation") +
  ylab("Number of peaks") + 
  xlab("Peak annotation")


# (b) Peak length by chromosome
plc <- ggplot(atlas[width < 2500,], aes(x = width, y = chrom)) +
  geom_density_ridges(stat = "binline",bins=50) + 
  xlim(c(0,2500)) + 
  xlab("Peak length (bp)") + ylab("Chromosome") + 
  theme_ridges(grid = FALSE)


# (c) Fine-grained histogram of peak lengths
fgpl <- ggplot(atlas[width < 2500,], aes(x = width, y = annotation)) +
  geom_density_ridges(aes(fill=annotation), stat = "binline",bins=80) +
  labs(title="Peak lengths by annotation")

# (d) Gene complexity plot: number of peaks / gene

# Count peaks / gene
gene_count = atlas[!is.na(SYMBOL), .N, by = SYMBOL]
setnames(gene_count, "N", "count")
gene_lengths = unique(atlas[!is.na(SYMBOL), geneLength, by = SYMBOL])
setkey(gene_count, SYMBOL)
setkey(gene_lengths, SYMBOL)
gene_length_dt = gene_lengths[gene_count]

gcp = ggplot(gene_length_dt[count < 100,], aes(x=geneLength, y=count)) +
  geom_point(aes(alpha=1/40)) + guides(alpha=FALSE) +
  geom_density2d() + 
  scale_x_log10() + 
  geom_text_repel(
    data = gene_length_dt[count > 35 & count < 100,],
    aes(label = SYMBOL),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) + 
  ggtitle("Peaks per gene") + 
  xlab("Gene length (log10 bp)") + 
  ylab("Number of peaks")

setwd("../figures/")
pdf(file = "atlas_diagnostic_plots_refseq.pdf", width = 15, height = 13)

# compile plots into a list
pltList <- list()
pltList[[1]] <- pba
pltList[[2]] <- plc
pltList[[3]] <- fgpl
pltList[[4]] <- gcp

# display the plots in a grid
grid.arrange(grobs=pltList, ncol=2)
dev.off()


# plot pie chart of atlas, differentiable peaks atlas
total = atlas[,.N]
grouped_peaks = atlas[,.N, by=annotation]
grouped_peaks[, annotation := as.factor(annotation)]
grouped_peaks[, percent := N / total]
grouped_peaks = grouped_peaks[order(c(1,3,2,4,5,6))]
grouped_peaks[, pie_label := cumsum(percent) - percent / 2]

# pie chart of atlas
atlas_pie = ggplot(grouped_peaks, aes(x="", y=percent, fill=annotation)) + 
  geom_bar(stat = "identity") + 
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  geom_label_repel(aes(label = percent(percent), y = pie_label)) +
  theme(axis.text=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title=element_blank()) +
  ggtitle("All peaks in atlas (by annotation)") + 
  coord_polar(theta = "y")


# pie chart of diff'ble peaks: SAR vs P
diff_peaks_anno <- annotatePeak("../peaks/SAR_vs_P_diff_peaks.bed", tssRegion=c(-2000,2000), TxDb=txdb, level = "gene", annoDb="org.Hs.eg.db")
diff_peaks_dt = data.table(as.data.frame(diff_peaks_anno))
rm(diff_peaks_anno)
setnames(diff_peaks_dt, "seqnames", "chrom")
diff_peaks_dt[grepl("Exon",annotation), simple_annotation := "Exon"]
diff_peaks_dt[grepl("Intron",annotation), simple_annotation := "Intron"]
diff_peaks_dt[grepl("Promoter",annotation), simple_annotation := "Promoter"]
diff_peaks_dt[grepl("Intergenic",annotation), simple_annotation := "Intergenic"]
diff_peaks_dt[grepl("Downstream",annotation), simple_annotation := "Intergenic"]
diff_peaks_dt[grepl("5' UTR",annotation), simple_annotation := "5' UTR"]
diff_peaks_dt[grepl("3' UTR",annotation), simple_annotation := "3' UTR"]
setnames(diff_peaks_dt, c("annotation", "simple_annotation"), c("complex_annotation", "annotation"))

total = diff_peaks_dt[,.N]
grouped_peaks =  diff_peaks_dt[,.N, by=annotation]
grouped_peaks[, annotation := as.factor(annotation)]
grouped_peaks[, percent := N / total]
grouped_peaks = grouped_peaks[order(c(4,3,2,1,5,6))]
grouped_peaks[, pie_label := cumsum(percent) - percent / 2]

SAR_vs_P_atlas_pie = ggplot(grouped_peaks, aes(x="", y=percent, fill=annotation)) + 
  geom_bar(stat = "identity") + 
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  geom_label_repel(aes(label = percent(percent), y = pie_label)) +
  theme(axis.text=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title=element_blank()) +
  ggtitle("All peaks in atlas (by annotation)") + 
  coord_polar(theta = "y")

# pie cahrt of diff'ble peaks: P -> A -> SA -> ...
diff_peaks_anno <- annotatePeak("../peaks/P_A_SA_SAR_SARN_SARF_diff_peaks.bed", tssRegion=c(-2000,2000), TxDb=txdb, level = "gene", annoDb="org.Hs.eg.db")
diff_peaks_dt = data.table(as.data.frame(diff_peaks_anno))
rm(diff_peaks_anno)

setnames(diff_peaks_dt, "seqnames", "chrom")
diff_peaks_dt[grepl("Exon",annotation), simple_annotation := "Exon"]
diff_peaks_dt[grepl("Intron",annotation), simple_annotation := "Intron"]
diff_peaks_dt[grepl("Promoter",annotation), simple_annotation := "Promoter"]
diff_peaks_dt[grepl("Intergenic",annotation), simple_annotation := "Intergenic"]
diff_peaks_dt[grepl("Downstream",annotation), simple_annotation := "Intergenic"]
diff_peaks_dt[grepl("5' UTR",annotation), simple_annotation := "5' UTR"]
diff_peaks_dt[grepl("3' UTR",annotation), simple_annotation := "3' UTR"]
setnames(diff_peaks_dt, c("annotation", "simple_annotation"), c("complex_annotation", "annotation"))

total = diff_peaks_dt[,.N]
grouped_peaks =  diff_peaks_dt[,.N, by=annotation]
grouped_peaks[, annotation := as.factor(annotation)]
grouped_peaks[, percent := N / total]
grouped_peaks = grouped_peaks[c(2,3,1,4,5,6)]
grouped_peaks[, pie_label := cumsum(percent) - percent / 2]

P_A_diff_atlas_pie = ggplot(grouped_peaks, aes(x="", y=percent, fill=annotation)) + 
  geom_bar(stat = "identity") + 
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  geom_label_repel(aes(label = percent(percent), y = pie_label)) +
  theme(axis.text=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title=element_blank()) +
  ggtitle("All peaks in atlas (by annotation)") + 
  coord_polar(theta = "y")


# pie chart of diff'ble peaks: P -> S -> SA -> ...
diff_peaks_anno <- annotatePeak("../peaks/P_S_SA_SAR_SARN_SARF_diff_peaks.bed", tssRegion=c(-2000,2000), TxDb=txdb, level = "gene", annoDb="org.Hs.eg.db")
diff_peaks_dt = data.table(as.data.frame(diff_peaks_anno))
rm(diff_peaks_anno)
setnames(diff_peaks_dt, "seqnames", "chrom")
diff_peaks_dt[grepl("Exon",annotation), simple_annotation := "Exon"]
diff_peaks_dt[grepl("Intron",annotation), simple_annotation := "Intron"]
diff_peaks_dt[grepl("Promoter",annotation), simple_annotation := "Promoter"]
diff_peaks_dt[grepl("Intergenic",annotation), simple_annotation := "Intergenic"]
diff_peaks_dt[grepl("Downstream",annotation), simple_annotation := "Intergenic"]
diff_peaks_dt[grepl("5' UTR",annotation), simple_annotation := "5' UTR"]
diff_peaks_dt[grepl("3' UTR",annotation), simple_annotation := "3' UTR"]
setnames(diff_peaks_dt, c("annotation", "simple_annotation"), c("complex_annotation", "annotation"))

total = diff_peaks_dt[,.N]
grouped_peaks =  diff_peaks_dt[,.N, by=annotation]
grouped_peaks[, annotation := as.factor(annotation)]
grouped_peaks[, percent := N / total]
grouped_peaks = grouped_peaks[c(4,5,1,2,3,6)]
grouped_peaks[, pie_label := cumsum(percent) - percent / 2]

P_S_diff_atlas_pie = ggplot(grouped_peaks, aes(x="", y=percent, fill=annotation)) + 
  geom_bar(stat = "identity") + 
  guides(fill=guide_legend(override.aes=list(colour=NA))) +
  geom_label_repel(aes(label = percent(percent), y = pie_label)) +
  theme(axis.text=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank(),
        axis.title=element_blank()) +
  ggtitle("All peaks in atlas (by annotation)") + 
  coord_polar(theta = "y")



### Load the DESeq objects used to make differential plots
setwd('/Users/zamparol/projects/AML_ATAC/results/DESeq/objects')
dds_peaks = readRDS("dds_object.rds")
res_A_P = readRDS("res_A_P.rds")
res_SA_A = readRDS("res_SA_A.rds")
res_SAR_SA = readRDS("res_SAR_SA.rds")
res_S_P = readRDS("res_S_P.rds")
res_SA_S = readRDS("res_SA_S.rds")
res_SARF_SAR = readRDS("res_SARF_SAR.rds")
res_SARN_SAR = readRDS("res_SARN_SAR.rds")
res_SA_P = readRDS("res_SA_P.rds")
res_SAR_P = readRDS("res_SAR_P.rds")
res_SARF_P = readRDS("res_SARF_P.rds")
res_SARN_P = readRDS("res_SARN_P.rds")

res_A_P_dt = as.data.table(res_A_P)
res_S_P_dt = as.data.table(res_S_P)
res_SA_A_dt = as.data.table(res_SA_A)
res_SA_S_dt = as.data.table(res_SA_S)
res_SAR_SA_dt = as.data.table(res_SAR_SA)
res_SARN_SAR_dt = as.data.table(res_SARN_SAR)
res_SARF_SAR_dt = as.data.table(res_SARF_SAR)
res_SA_P_dt = as.data.table(res_SA_P)
res_SAR_P_dt = as.data.table(res_SAR_P)
res_SARF_P_dt = as.data.table(res_SARF_P)
res_SARN_P_dt = as.data.table(res_SARN_P)


### PCA plot
# Check that we're getting replicates clustering together
rld <- rlog(dds_peaks, blind=FALSE)
col_data = as.data.table(read_excel("../../../data/filename_to_sample.xlsx"))
col_data = col_data[Filename != "SRSF2_c1_8",]

# PCA of conditions where batch effects are removed
mat <- assay(rld)
#mat_nobatch <- limma::removeBatchEffect(mat, design=model.matrix(~ col_data$Condition), batch=col_data$Batch)
mat_nobatch <- limma::removeBatchEffect(mat, design=model.matrix(~ col_data$Condition), batch=col_data$Batch, batch2=col_data$Adapter)
assay(rld) <- mat_nobatch
data <- plotPCA(rld, intgroup=c("Condition"), returnData=TRUE)  ### pick ntop based on number of up & down peaks.
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  ggtitle("PCA of log transformed read counts: with batch correction") +
  geom_text_repel(
    aes(label = name),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) 


# MA plots versus P
DESeq2::plotMA(res_A_P, main="MA plot for gains, losses in A vs P", ylim=c(-5,5))
DESeq2::plotMA(res_S_P, main="MA plot for gains, losses in S vs P", ylim=c(-5,5))
DESeq2::plotMA(res_SA_P, main="MA plot for gains, losses in SA vs P", ylim=c(-5,5))
DESeq2::plotMA(res_SAR_P, main="MA plot for gains, losses in SAR vs P", ylim=c(-5,5))
DESeq2::plotMA(res_SARF_P, main="MA plot for gains, losses in SARF vs P", ylim=c(-5,5))
DESeq2::plotMA(res_SARN_P, main="MA plot for gains, losses in SARN vs P", ylim=c(-5,5))

# Remaining MA plots in (P -> A / S -> SA -> SAR -> SARN / SARF)
DESeq2::plotMA(res_SA_S, main="MA plot for gains, losses in SA vs S", ylim=c(-5,5))
DESeq2::plotMA(res_SA_A, main="MA plot for gains, losses in SA vs A", ylim=c(-5,5))
DESeq2::plotMA(res_SAR_SA, main="MA plot for gains, losses in SAR vs SA", ylim=c(-5,5))
DESeq2::plotMA(res_SARN_SAR, main="MA plot for gains, losses in SARN vs SAR", ylim=c(-5,5))
DESeq2::plotMA(res_SARF_SAR, main="MA plot for gains, losses in SARF vs SAR", ylim=c(-5,5))

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
atlas = data.table(annotated_atlas)
diff_peaks = read.delim(file="P-A-SA-SAR_clusters.bed", sep="\t", header=TRUE)
diff_peaks = diff_peaks[,c("X.chrom", "start","end","deepTools_group")]
diff_peaks_dt = data.table(diff_peaks)
setnames(diff_peaks_dt, "X.chrom", "chrom")
setnames(atlas, "seqnames", "chrom")
setkeyv(atlas, c("chrom", "start", "end"))
setkeyv(diff_peaks_dt, c("chrom", "start", "end"))
annotated_diff_peaks = atlas[diff_peaks_dt,]

annotated_diff_peaks[, group_total := .N, by=deepTools_group]
annotated_diff_peaks[, group_percentage := .N / group_total, by=.(deepTools_group, annot)]
grouped_percentages = unique(annotated_diff_peaks[,.(annot, group_percentage, deepTools_group)])
ggplot(grouped_percentages, aes(x=annot, fill=annot)) + 
  geom_col(aes(y=percent(group_percentage))) + 
  facet_wrap(~ deepTools_group) + 
  xlab("Annotation") + 
  ylab("Percentage of peaks") + 
  ggtitle("Differential peak annotation breakdown by plotHeatmap cluster")
