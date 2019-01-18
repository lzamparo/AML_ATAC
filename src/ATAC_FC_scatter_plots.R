require(data.table)
require(DESeq2)
require(ggplot2)
require(ashr)
require(cowplot)
require(gridExtra)

setwd("~/projects/AML_ATAC/results/DESeq/objects")
dds_peaks <- readRDS("dds_object.rds")
sig_alpha = 0.05

### track the DA peaks from A versus P across the SA versus A, SA versus P contrasts

#d get the DA peaks for A versus P
res_A_P <- results(dds_peaks, contrast = c("Condition", "A", "P"), parallel=TRUE)
res_A_P_shrunk <- lfcShrink(dds_peaks, coef="Condition_A_vs_P", type="ashr", parallel=TRUE)

idx <- res_A_P$padj < 0.05
A_P_diff_peaks <- rownames(res_A_P)[idx]
dt_A = data.table(as.data.frame(res_A_P))
dt_A[, peaks := rownames(res_A_P)]
setnames(dt_A, c("log2FoldChange"),c("A_vs_P_log2FoldChange"))

# get the DA peaks for SA versus A, SA versus P
res_SA_A <- results(dds_peaks, contrast=c("Condition", "SA", "A"), parallel=TRUE)
res_SA_A_shrunk <- lfcShrink(dds_peaks, contrast=c("Condition","SA","A"), type="ashr", parallel=TRUE)
dt_A[, SA_vs_A_log2FoldChange := res_SA_A$log2FoldChange]
dt_A[, SA_vs_A_padj := res_SA_A$padj]

res_SA_P <- results(dds_peaks, contrast=c("Condition", "SA", "P"), parallel=TRUE)
res_SA_P_shrunk <- lfcShrink(dds_peaks, coef="Condition_SA_vs_P", type="ashr", parallel=TRUE)
dt_A[, SA_vs_P_log2FoldChange := res_SA_P$log2FoldChange]
dt_A[, SA_vs_P_padj := res_SA_P$padj]

# mark the DA points in A vs P, SA vs A, SA vs P
dt_A[, A_vs_P_DA := "unchanged"]
dt_A[A_vs_P_log2FoldChange > 0 & padj < 0.05, A_vs_P_DA := "opening"]
dt_A[A_vs_P_log2FoldChange < 0 & padj < 0.05, A_vs_P_DA := "closing"]

dt_A[, SA_vs_A_DA := "unchanged"]
dt_A[SA_vs_A_log2FoldChange > 0 & SA_vs_A_padj < 0.05, SA_vs_A_DA := "opening"]
dt_A[SA_vs_A_log2FoldChange < 0 & SA_vs_A_padj < 0.05, SA_vs_A_DA := "closing"]

dt_A[, SA_vs_P_DA := "unchanged"]
dt_A[SA_vs_P_log2FoldChange > 0 & SA_vs_P_padj < 0.05, SA_vs_P_DA := "opening"]
dt_A[SA_vs_P_log2FoldChange < 0 & SA_vs_P_padj < 0.05, SA_vs_P_DA := "closing"]

# plot the sig peaks as an un-ordered scatter plot: A vs P -> SA vs P
scatter_AP_SAP = ggplot(dt_A, aes(x=A_vs_P_log2FoldChange, y=SA_vs_P_log2FoldChange, color=A_vs_P_DA, shape=SA_vs_P_DA)) + 
  geom_point(alpha = 2/10) + 
  geom_vline(xintercept=1, linetype=3, alpha=2/10) + 
  geom_vline(xintercept=-1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=-1, linetype=3, alpha=2/10) + 
  scale_color_manual(name = "Status in A versus P", values = c("opening" = "red", "unchanged" = "gray", "closing" = "blue")) + 
  scale_shape_manual(name = "Status in SA versus P", values = c("opening" = "triangle", "unchanged" = "square", "closing" = "circle")) + 
  ggtitle("Log2FC of (A vs P) against (SA vs P)") + 
  xlab("A vs P LFC") + 
  ylab("SA vs P LFC")

# plot the sig peaks as an un-ordered scatter plot: A vs P -> SA vs A
scatter_AP_SAA = ggplot(dt_A, aes(x=A_vs_P_log2FoldChange, y=SA_vs_A_log2FoldChange, color=A_vs_P_DA, shape=SA_vs_A_DA)) + 
  geom_point(alpha = 2/10) +
  geom_vline(xintercept=1, linetype=3, alpha=2/10) + 
  geom_vline(xintercept=-1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=-1, linetype=3, alpha=2/10) + 
  scale_color_manual(name = "Status in A versus P", values = c("opening" = "red", "unchanged" = "gray", "closing" = "blue")) + 
  scale_shape_manual(name = "Status in SA versus A", values = c("opening" = "triangle", "unchanged" = "square", "closing" = "circle")) + 
  ggtitle("Log2FC of (A vs P) against (SA vs A)") + 
  xlab("A vs P Log2FC") + 
  ylab("SA vs A Log2FC")



### track the DA peaks from S versus P across the SA versus S, SA versus P contrasts
res_S_P <- lfcShrink(dds_peaks, coef="Condition_S_vs_P", type="ashr", parallel=TRUE)
idx <- res_S_P$padj < 0.05
S_P_diff_peaks <- rownames(res_S_P)[idx]
dt_S = data.table(as.data.frame(res_S_P))
dt_S[, peaks := rownames(res_S_P)]
setnames(dt_S, c("log2FoldChange"),c("S_vs_P_log2FoldChange"))

# get the DA peaks for SA versus S, SA versus P
res_SA_S <- lfcShrink(dds_peaks, contrast=c("Condition","SA","S"), type="ashr", parallel=TRUE)
dt_S[, SA_vs_S_log2FoldChange := res_SA_S$log2FoldChange]
dt_S[, SA_vs_S_padj := res_SA_S$padj]

dt_S[, SA_vs_P_log2FoldChange := res_SA_P$log2FoldChange]
dt_S[, SA_vs_P_padj := res_SA_P$padj]

# mark the DA points in S vs P, SA vs S, SA vs P
dt_S[, S_vs_P_DA := "unchanged"]
dt_S[S_vs_P_log2FoldChange > 0 & padj < 0.05, S_vs_P_DA := "opening"]
dt_S[S_vs_P_log2FoldChange < 0 & padj < 0.05, S_vs_P_DA := "closing"]

dt_S[, SA_vs_A_DA := "unchanged"]
dt_S[SA_vs_A_log2FoldChange > 0 & SA_vs_A_padj < 0.05, SA_vs_A_DA := "opening"]
dt_S[SA_vs_A_log2FoldChange < 0 & SA_vs_A_padj < 0.05, SA_vs_A_DA := "closing"]

dt_S[, SA_vs_P_DA := "unchanged"]
dt_S[SA_vs_P_log2FoldChange > 0 & SA_vs_P_padj < 0.05, SA_vs_P_DA := "opening"]
dt_S[SA_vs_P_log2FoldChange < 0 & SA_vs_P_padj < 0.05, SA_vs_P_DA := "closing"]

# plot the sig peaks as an un-ordered scatter plot: A vs P -> SA vs P
scatter_AP_SAP = ggplot(dt_S, aes(x=A_vs_P_log2FoldChange, y=SA_vs_P_log2FoldChange, color=A_vs_P_DA, shape=SA_vs_P_DA)) + 
  geom_point(alpha = 2/10) + 
  geom_vline(xintercept=1, linetype=3, alpha=2/10) + 
  geom_vline(xintercept=-1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=-1, linetype=3, alpha=2/10) + 
  scale_color_manual(name = "Status in A versus P", values = c("opening" = "red", "unchanged" = "gray", "closing" = "blue")) + 
  scale_shape_manual(name = "Status in SA versus P", values = c("opening" = "triangle", "unchanged" = "square", "closing" = "circle")) + 
  ggtitle("Log2FC of (A vs P) against (SA vs P)") + 
  xlab("A vs P LFC") + 
  ylab("SA vs P LFC")

# plot the sig peaks as an un-ordered scatter plot: A vs P -> SA vs A
scatter_AP_SAA = ggplot(dt_S, aes(x=A_vs_P_log2FoldChange, y=SA_vs_A_log2FoldChange, color=A_vs_P_DA, shape=SA_vs_A_DA)) + 
  geom_point(alpha = 2/10) +
  geom_vline(xintercept=1, linetype=3, alpha=2/10) + 
  geom_vline(xintercept=-1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=1, linetype=3, alpha=2/10) + 
  geom_hline(yintercept=-1, linetype=3, alpha=2/10) + 
  scale_color_manual(name = "Status in A versus P", values = c("opening" = "red", "unchanged" = "gray", "closing" = "blue")) + 
  scale_shape_manual(name = "Status in SA versus A", values = c("opening" = "triangle", "unchanged" = "square", "closing" = "circle")) + 
  ggtitle("Log2FC of (A vs P) against (SA vs A)") + 
  xlab("A vs P Log2FC") + 
  ylab("SA vs A Log2FC")


### arrange, save the plots
setwd("~/projects/AML_ATAC/results/figures/")
pdf(file = "peak_tracking_plots.pdf", width = 15, height = 13)

# compile plots into a list
pltList <- list()
pltList[[1]] <- scatter_AP_SAP
pltList[[2]] <- scatter_AP_SAA
pltList[[3]] <- fgpl
pltList[[4]] <- gcp

# display the plots in a grid
grid.arrange(grobs=pltList, ncol=2)
dev.off()
