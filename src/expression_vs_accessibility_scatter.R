require(DESeq2)
require(data.table)
require(ggplot2)
require(latex2exp)

# Just look at condition specific changes in expression vs peak log FC?

# Read the RNA-seq DESeq results into tables

A_vs_P = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/A_vs_P.txt", sep="")
S_vs_P = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/S_vs_P.txt", sep="")
SA_vs_A = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SA_vs_A.txt", sep="")
SA_vs_S = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SA_vs_S.txt", sep="")
SAR_vs_SA = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SAR_vs_SA.txt", sep="")

# Accumulate into one data table
A_vs_P_dt = as.data.table(A_vs_P)
A_vs_P_dt[, condition := "A_vs_P"]
S_vs_P_dt = as.data.table(S_vs_P)
S_vs_P_dt[, condition := "S_vs_P"]

SA_vs_A_dt = as.data.table(SA_vs_A)
SA_vs_A_dt[, condition := "SA_vs_A"]
SA_vs_S_dt = as.data.table(SA_vs_S)
SA_vs_S_dt[, condition := "SA_vs_S"]

SAR_vs_SA_dt = as.data.table(SAR_vs_SA)
SAR_vs_SA_dt[, condition := "SAR_vs_SA"]

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

# need to join expression data on ATAC data, group by Gene, plot expn FC on X, peak FC on Y
setkey(A_vs_P_dt, Gene)
setkey(A_vs_P_atac_dt, Gene)
A_P_joint = A_vs_P_atac_dt[A_vs_P_dt]

# sample 1000 random genes
AP_gene_sample = sample(unique(A_P_joint[,Gene], 1000))

A_P_scatter = ggplot(A_P_joint[Gene %in% AP_gene_sample], aes(x=log2FoldChange, y=log2FoldChangeATAC)) +
  geom_point() + 
  ggtitle("A versus P: gene expression vs accessibility") + 
  xlab(TeX("$\\log_{2}(FC)$")) + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

setkey(S_vs_P_dt, Gene)
setkey(S_vs_P_atac_dt, Gene)
S_P_joint = S_vs_P_atac_dt[S_vs_P_dt]

# sample 1000 random genes
SP_gene_sample = sample(unique(S_P_joint[,Gene], 1000))

S_P_scatter = ggplot(S_P_joint[Gene %in% SP_gene_sample], aes(x=log2FoldChange, y=log2FoldChangeATAC)) +
  geom_point() + 
  ggtitle("S versus P: gene expression vs accessibility") + 
  xlab(TeX("$\\log_{2}(FC)$")) + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))


setkey(SA_vs_A_dt, Gene)
setkey(SA_vs_A_atac_dt, Gene)
SA_A_joint = SA_vs_A_atac_dt[SA_vs_A_dt]

# sample 1000 random genes
SAA_gene_sample = sample(unique(SA_A_joint[,Gene], 1000))

SA_A_scatter = ggplot(SA_A_joint[Gene %in% SAA_gene_sample], aes(x=log2FoldChange, y=log2FoldChangeATAC)) +
  geom_point() + 
  ggtitle("SA versus A: gene expression vs accessibility") + 
  xlab(TeX("$\\log_{2}(FC)$")) + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))


setkey(SA_vs_S_dt, Gene)
setkey(SA_vs_S_atac_dt, Gene)
SA_S_joint = SA_vs_S_atac_dt[SA_vs_S_dt]

# sample 1000 random genes
SAS_gene_sample = sample(unique(SA_S_joint[,Gene], 1000))

SA_S_scatter = ggplot(SA_S_joint[Gene %in% SAS_gene_sample], aes(x=log2FoldChange, y=log2FoldChangeATAC)) +
  geom_point() + 
  ggtitle("SA versus S: gene expression vs accessibility") + 
  xlab(TeX("$\\log_{2}(FC)$")) + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))


setkey(SAR_vs_SA_dt, Gene)
setkey(SAR_vs_SA_atac_dt, Gene)
SAR_SA_joint = SAR_vs_SA_atac_dt[SAR_vs_SA_dt]

# sample 1000 random genes
SAR_SA_gene_sample = sample(unique(SAR_SA_joint[,Gene], 1000))

SAR_SA_scatter = ggplot(SAR_SA_joint[Gene %in% SAR_SA_gene_sample], aes(x=log2FoldChange, y=log2FoldChangeATAC)) +
  geom_point() + 
  ggtitle("A versus P: gene expression vs accessibility") + 
  xlab(TeX("$\\log_{2}(FC)$")) + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

