require(DESeq2)
require(data.table)
require(ggplot2)
require(latex2exp)

# Just look at condition specific changes in expression vs peak log FC?

# Read the RNA-seq DESeq results into data.table
rna_dds_obj = readRDS("~/projects/AML_ATAC/results/DESeq/objects/RNA/all_res.rds")
gene_names = row.names(rna_dds_obj)
rna_dds_dt = as.data.table(rna_dds_obj)
rna_dds_dt$Gene = gene_names
rm(list = c("gene_names", "rna_dds_obj"))

# Load the annotated peak atlas
annotated_peak_atlas = data.table(read.csv("~/projects/AML_ATAC/results/peaks/all_conditions_peak_atlas_annotated.csv"))

# load the chromatin accesibility data for each condition change 
A_vs_P = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_A_P.rds")
A_vs_P_atac_dt = as.data.table(A_vs_P)
A_vs_P_atac_dt$Gene = annotated_peak_atlas$SYMBOL
setnames(A_vs_P_atac_dt, "log2FoldChange", "log2FoldChangeATAC") 

S_vs_P = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_S_P.rds")
S_vs_P_atac_dt = as.data.table(S_vs_P)
S_vs_P_atac_dt$Gene = annotated_peak_atlas$SYMBOL
setnames(S_vs_P_atac_dt, "log2FoldChange", "log2FoldChangeATAC")

SA_vs_A = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SA_A.rds")
SA_vs_A_atac_dt = as.data.table(SA_vs_A)
SA_vs_A_atac_dt$Gene = annotated_peak_atlas$SYMBOL
setnames(SA_vs_A_atac_dt, "log2FoldChange", "log2FoldChangeATAC")

SA_vs_S = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SA_S.rds")
SA_vs_S_atac_dt = as.data.table(SA_vs_S)
SA_vs_S_atac_dt$Gene = annotated_peak_atlas$SYMBOL
setnames(SA_vs_S_atac_dt, "log2FoldChange", "log2FoldChangeATAC")

SAR_vs_SA = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SAR_SA.rds")
SAR_vs_SA_atac_dt = as.data.table(SAR_vs_SA)
SAR_vs_SA_atac_dt$Gene = annotated_peak_atlas$SYMBOL
setnames(SAR_vs_SA_atac_dt, "log2FoldChange", "log2FoldChangeATAC")

SAR_vs_P = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SAR_P.rds")
SAR_vs_P_atac_dt = as.data.table(SAR_vs_P)
SAR_vs_P_atac_dt$Gene = annotated_peak_atlas$SYMBOL
setnames(SAR_vs_P_atac_dt, "log2FoldChange", "log2FoldChangeATAC")

### functions to do Wilcoxon RST in both greater, less directions
wrst_up <- function(df, all_vals){
  res <- wilcox.test(df$log2FoldChangeATAC, all_vals, alternative = "greater")
  res$p.value
}

wrst_down <- function(df, all_vals){
  res <- wilcox.test(df$log2FoldChangeATAC, all_vals, alternative = "less")
  res$p.value
}

# need to join expression data on ATAC data, group by Gene, plot expn FC on X, peak FC on Y
setkey(rna_dds_dt, Gene)
setkey(A_vs_P_atac_dt, Gene)
A_P_joint = A_vs_P_atac_dt[rna_dds_dt]
A_P_joint[, log2FC_med_ATAC := median(log2FoldChangeATAC), by = Gene]
A_P_joint = unique(A_P_joint[,.(Gene, AvP_log2FC, log2FC_med_ATAC)])

# some colouring of points based on the being over-expressed 
AP_all_peaks_fc = A_vs_P_atac_dt$log2FoldChangeATAC
ups = unique(A_P_joint[log2FC_med_ATAC > 0, Gene])
downs = unique(A_P_joint[log2FC_med_ATAC < 0, Gene])
A_P_joint[Gene %in% ups, wrst := wrst_up(.SD, AP_all_peaks_fc), by = Gene]
A_P_joint[Gene %in% downs, wrst := wrst_down(.SD, AP_all_peaks_fc), by = Gene]
A_P_joint = unique(A_P_joint[,.(Gene, AvP_log2FC, log2FC_med_ATAC, wrst)])
A_P_spearman = cor(A_P_joint$AvP_log2FC, A_P_joint$log2FC_med_ATAC, use="pairwise.complete.obs", method="spearman")
s_label = paste("Spearman: ", format(A_P_spearman, digits=3))

A_P_scatter = ggplot(A_P_joint[wrst > 0.05,], aes(x=AvP_log2FC, y=log2FC_med_ATAC)) +
  geom_point() + 
  geom_point(data=A_P_joint[log2FC_med_ATAC > 1 & wrst < 0.05,], aes(x=AvP_log2FC, y=log2FC_med_ATAC), colour="red") +
  geom_point(data=A_P_joint[log2FC_med_ATAC < -1 & wrst < 0.05,], aes(x=AvP_log2FC, y=log2FC_med_ATAC), colour="green") + 
  ggtitle("A versus P: gene expression vs median accessibility") + 
  xlab(TeX("RNA $\\log_{2}(FC)$")) + 
  ylab(TeX("median ATAC $\\log_{2}(FC)$")) + 
  geom_text(aes(-5, 3, label=s_label))

### TODO: replicate in other conditions, esp SAR vs P

setkey(S_vs_P_dt, Gene)
setkey(S_vs_P_atac_dt, Gene)
S_P_joint = S_vs_P_atac_dt[S_vs_P_dt]
S_P_joint[, log2FC_med_ATAC := median(log2FoldChangeATAC), by = Gene]
S_P_joint = unique(S_P_joint[,.(Gene, log2FoldChange, log2FC_med_ATAC)])

S_P_scatter = ggplot(S_P_joint, aes(x=log2FoldChange, y=log2FC_med_ATAC)) +
  geom_point() + 
  ggtitle("S versus P: gene expression vs median accessibility") + 
  xlab(TeX("RNA $\\log_{2}(FC)$")) + 
  ylab(TeX("median ATAC $\\log_{2}(FC)$"))

setkey(SA_vs_A_dt, Gene)
setkey(SA_vs_A_atac_dt, Gene)
SA_A_joint = SA_vs_A_atac_dt[SA_vs_A_dt]
SA_A_joint[, log2FC_med_ATAC := median(log2FoldChangeATAC), by = Gene]
SA_A_joint = unique(SA_A_joint[,.(Gene, log2FoldChange, log2FC_med_ATAC)])

SA_A_scatter = ggplot(SA_A_joint, aes(x=log2FoldChange, y=log2FC_med_ATAC)) +
  geom_point() + 
  ggtitle("SA versus A: gene expression vs median accessibility") + 
  xlab(TeX("RNA $\\log_{2}(FC)$")) + 
  ylab(TeX("median ATAC $\\log_{2}(FC)$"))

setkey(SA_vs_S_dt, Gene)
setkey(SA_vs_S_atac_dt, Gene)
SA_S_joint = SA_vs_S_atac_dt[SA_vs_S_dt]
SA_S_joint[, log2FC_med_ATAC := median(log2FoldChangeATAC), by = Gene]
SA_S_joint = unique(SA_S_joint[,.(Gene, log2FoldChange, log2FC_med_ATAC)])

SA_S_scatter = ggplot(SA_S_joint, aes(x=log2FoldChange, y=log2FC_med_ATAC)) +
  geom_point() + 
  ggtitle("SA versus S: gene expression vs median accessibility") + 
  xlab(TeX("RNA $\\log_{2}(FC)$")) + 
  ylab(TeX("median ATAC $\\log_{2}(FC)$"))

setkey(SAR_vs_SA_dt, Gene)
setkey(SAR_vs_SA_atac_dt, Gene)
SAR_SA_joint = SAR_vs_SA_atac_dt[SAR_vs_SA_dt]
SAR_SA_joint[, log2FC_med_ATAC := median(log2FoldChangeATAC), by = Gene]
SAR_SA_joint = unique(SAR_SA_joint[,.(Gene, log2FoldChange, log2FC_med_ATAC)])

SAR_SA_scatter = ggplot(SAR_SA_joint, aes(x=log2FoldChange, y=log2FC_med_ATAC)) +
  geom_point() + 
  ggtitle("SAR versus SA: gene expression vs accessibility") + 
  xlab(TeX("RNA $\\log_{2}(FC)$")) + 
  ylab(TeX("median ATAC $\\log_{2}(FC)$"))

#all_SAR_vs_P_lfcs = SAR_vs_P_atac_dt$log2FoldChangeATAC

# need to join expression data on ATAC data, group by Gene, plot expn FC on X, peak FC on Y
setkey(rna_dds_dt, Gene)
setkey(SAR_vs_P_atac_dt, Gene)
SAR_P_joint = SAR_vs_P_atac_dt[rna_dds_dt]
SAR_P_joint[, log2FC_med_ATAC := median(log2FoldChangeATAC), by = Gene]
SAR_P_joint = unique(SAR_P_joint[,.(Gene, log2FoldChangeATAC, SARvP_log2FC, log2FC_med_ATAC)])

# some colouring of points based on the being over-expressed 
SARvP_all_peaks_fc = SAR_vs_P_atac_dt$log2FoldChangeATAC
ups = unique(SAR_P_joint[log2FC_med_ATAC > 0, Gene])
downs = unique(SAR_P_joint[log2FC_med_ATAC < 0, Gene])
SAR_P_joint[Gene %in% ups, wrst := wrst_up(.SD, SARvP_all_peaks_fc), by = Gene]
SAR_P_joint[Gene %in% downs, wrst := wrst_down(.SD, SARvP_all_peaks_fc), by = Gene]
SAR_P_joint = unique(SAR_P_joint[,.(Gene, SARvP_log2FC, log2FC_med_ATAC, wrst)])
SAR_P_spearman = cor(SAR_P_joint$SARvP_log2FC, SAR_P_joint$log2FC_med_ATAC, use="pairwise.complete.obs", method="spearman")
s_label = paste("Spearman: ", format(SAR_P_spearman, digits=3)) 


SAR_P_scatter = ggplot(SAR_P_joint[wrst > 0.05,], aes(x=SARvP_log2FC, y=log2FC_med_ATAC)) +
  geom_point() + 
  geom_point(data=SAR_P_joint[log2FC_med_ATAC > 1 & wrst < 0.05,], aes(x=SARvP_log2FC, y=log2FC_med_ATAC), colour="red") +
  geom_point(data=SAR_P_joint[log2FC_med_ATAC < -1 & wrst < 0.05,], aes(x=SARvP_log2FC, y=log2FC_med_ATAC), colour="green") + 
  ggtitle("SAR versus P: gene expression vs median accessibility") + 
  xlab(TeX("RNA $\\log_{2}(FC)$")) + 
  ylab(TeX("median ATAC $\\log_{2}(FC)$")) + 
  geom_text(aes(-5, 3, label=s_label))

