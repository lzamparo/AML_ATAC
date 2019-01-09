require(data.table)
require(DESeq2)
require(ggplot2)
require(gridExtra)

# compare the log2FC distributions of each DDS object

setwd("~/projects/AML_ATAC/results/DESeq/objects")
flank_normalized_dds <- readRDS("dds_object.rds")
peak_normalized_dds <- readRDS("dds_object_peak_est_sf.rds")
sig_alpha = 0.05

# A versus P
flank_A_P <- results(flank_normalized_dds, contrast=c("Condition", "A", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
peak_A_P <- results(peak_normalized_dds, contrast=c("Condition", "A", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
flank_A_P_df = as.data.frame(flank_A_P)
flank_A_P_df$norm = "flank"
peak_A_P_df = as.data.frame(peak_A_P)
peak_A_P_df$norm = "peak"
A_P_df = rbind(flank_A_P_df, peak_A_P_df)
A_vs_P = ggplot(A_P_df, aes(x=log2FoldChange, colour=norm, alpha=0.4)) + geom_density() + guides(alpha="none") + geom_vline(xintercept = 0, linetype="dotted") + ggtitle("Smoothed Log2FC distibution: A vs P")
rm(list = c("flank_A_P", "peak_A_P", "flank_A_P_df", "peak_A_P_df", "A_P_df"))

# S versus P
flank_S_P <- results(flank_normalized_dds, contrast=c("Condition", "S", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
peak_S_P <- results(peak_normalized_dds, contrast=c("Condition", "S", "P"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
flank_S_P_df = as.data.frame(flank_S_P)
flank_S_P_df$norm = "flank"
peak_S_P_df = as.data.frame(peak_S_P)
peak_S_P_df$norm = "peak"
S_P_df = rbind(flank_S_P_df, peak_S_P_df)
S_vs_P = ggplot(S_P_df, aes(x=log2FoldChange, colour=norm, alpha=0.4)) + geom_density() + guides(alpha="none") + geom_vline(xintercept = 0, linetype="dotted") + ggtitle("Smoothed Log2FC distibutions: S vs P")
rm(list = c("flank_S_P", "peak_S_P", "flank_S_P_df", "peak_S_P_df", "S_P_df"))


# SA versus A
flank_SA_A <- results(flank_normalized_dds, contrast=c("Condition", "SA", "A"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
peak_SA_A <- results(peak_normalized_dds, contrast=c("Condition", "SA", "A"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
flank_SA_A_df = as.data.frame(flank_SA_A)
flank_SA_A_df$norm = "flank"
peak_SA_A_df = as.data.frame(peak_SA_A)
peak_SA_A_df$norm = "peak"
SA_A_df = rbind(flank_SA_A_df, peak_SA_A_df)
SA_vs_A = ggplot(SA_A_df, aes(x=log2FoldChange, colour=norm, alpha=0.4)) + geom_density() + guides(alpha="none") + geom_vline(xintercept = 0, linetype="dotted") + ggtitle("Smoothed Log2FC distibutions: SA vs A")
rm(list = c("flank_SA_A", "peak_SA_A", "flank_SA_A_df", "peak_SA_A_df", "SA_A_df"))


# SA versus S
flank_SA_S <- results(flank_normalized_dds, contrast=c("Condition", "SA", "S"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
peak_SA_S <- results(peak_normalized_dds, contrast=c("Condition", "SA", "S"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
flank_SA_S_df = as.data.frame(flank_SA_S)
flank_SA_S_df$norm = "flank"
peak_SA_S_df = as.data.frame(peak_SA_S)
peak_SA_S_df$norm = "peak"
SA_S_df = rbind(flank_SA_S_df, peak_SA_S_df)
SA_vs_S = ggplot(SA_S_df, aes(x=log2FoldChange, colour=norm, alpha=0.4)) + geom_density() + guides(alpha="none") + geom_vline(xintercept = 0, linetype="dotted") + ggtitle("Smoothed Log2FC distibutions: SA vs S")
rm(list = c("flank_SA_S", "peak_SA_S", "flank_SA_S_df", "peak_SA_S_df", "SA_S_df"))

# SAR versus SA
flank_SAR_SA <- results(flank_normalized_dds, contrast=c("Condition", "SAR", "SA"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
peak_SAR_SA <- results(peak_normalized_dds, contrast=c("Condition", "SAR", "SA"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
flank_SAR_SA_df = as.data.frame(flank_SAR_SA)
flank_SAR_SA_df$norm = "flank"
peak_SAR_SA_df = as.data.frame(peak_SAR_SA)
peak_SAR_SA_df$norm = "peak"
SAR_SA_df = rbind(flank_SAR_SA_df, peak_SAR_SA_df)
SAR_vs_SA = ggplot(SAR_SA_df, aes(x=log2FoldChange, colour=norm, alpha=0.4)) + geom_density() + guides(alpha="none") + geom_vline(xintercept = 0, linetype="dotted") + ggtitle("Smoothed Log2FC distibutions: SAR vs SA")
rm(list = c("flank_SAR_SA", "peak_SAR_SA", "flank_SAR_SA_df", "peak_SAR_SA_df", "SAR_SA_df"))

# SARN versus SAR
flank_SARN_SAR <- results(flank_normalized_dds, contrast=c("Condition", "SARN", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
peak_SARN_SAR <- results(peak_normalized_dds, contrast=c("Condition", "SARN", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
flank_SARN_SAR_df = as.data.frame(flank_SARN_SAR)
flank_SARN_SAR_df$norm = "flank"
peak_SARN_SAR_df = as.data.frame(peak_SARN_SAR)
peak_SARN_SAR_df$norm = "peak"
SARN_SAR_df = rbind(flank_SARN_SAR_df, peak_SARN_SAR_df)
SARN_vs_SAR = ggplot(SARN_SAR_df, aes(x=log2FoldChange, colour=norm, alpha=0.4)) + geom_density() + guides(alpha="none") + geom_vline(xintercept = 0, linetype="dotted") + ggtitle("Smoothed Log2FC distibutions: SARN vs SAR")
rm(list = c("flank_SARN_SAR", "peak_SARN_SAR", "flank_SARN_SAR_df", "peak_SARN_SAR_df", "SARN_SAR_df"))

# SARF versus SAR
flank_SARF_SAR <- results(flank_normalized_dds, contrast=c("Condition", "SARF", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
peak_SARF_SAR <- results(peak_normalized_dds, contrast=c("Condition", "SARF", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
flank_SARF_SAR_df = as.data.frame(flank_SARF_SAR)
flank_SARF_SAR_df$norm = "flank"
peak_SARF_SAR_df = as.data.frame(peak_SARF_SAR)
peak_SARF_SAR_df$norm = "peak"
SARF_SAR_df = rbind(flank_SARF_SAR_df, peak_SARF_SAR_df)
SARF_vs_SAR = ggplot(SARF_SAR_df, aes(x=log2FoldChange, colour=norm, alpha=0.4)) + geom_density() + guides(alpha="none") + geom_vline(xintercept = 0, linetype="dotted") + ggtitle("Smoothed Log2FC distibutions: SARF vs SAR")
rm(list = c("flank_SARF_SAR", "peak_SARF_SAR", "flank_SARF_SAR_df", "peak_SARF_SAR_df", "SARF_SAR_df"))

# get the plots together, write out to pdf
setwd("../../figures/")
pdf(file = "normalization_diagnostic_plots.pdf", width = 12, height = 24)

# compile plots into a list
pltList <- list()
pltList[[1]] <- A_vs_P
pltList[[2]] <- S_vs_P
pltList[[3]] <- SA_vs_A
pltList[[4]] <- SA_vs_S
pltList[[5]] <- SAR_vs_SA
pltList[[6]] <- SARN_vs_SAR
pltList[[7]] <- SARF_vs_SAR

# display the plots in a grid
grid.arrange(grobs=pltList, ncol=2)
dev.off()
