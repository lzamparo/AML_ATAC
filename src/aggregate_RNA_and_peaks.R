require(DESeq2)
require(data.table)
require(ggplot2)
require(latex2exp)

# Read the RNA-seq DESeq results into tables

A_vs_P = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/A_vs_P.txt", sep="")
SA_vs_A = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SA_vs_A.txt", sep="")
SAR_vs_SA = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SAR_vs_SA.txt", sep="")

# Accumulate into one data table
A_vs_P_dt = as.data.table(A_vs_P)
A_vs_P_dt[, condition := "A_vs_P"]
SA_vs_A_dt = as.data.table(SA_vs_A)
SA_vs_A_dt[, condition := "SA_vs_A"]
SAR_vs_SA_dt = as.data.table(SAR_vs_SA)
SAR_vs_SA_dt[, condition := "SAR_vs_SA"]

# Load the annotated peak atlas
annotated_peak_atlas = data.table(read.csv("~/projects/AML_ATAC/results/peaks/all_conditions_peak_atlas_annotated.csv"))

# select the top 50 & bottom 50 DE genes for each condition
A_vs_P_MD_expression_dt = data.table(rbind(A_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:50], A_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:50]))
SA_vs_A_MD_expression_dt = data.table(rbind(SA_vs_A_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:50], SA_vs_A_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:50]))
# The top 50, bottom 50 significant up regulated & down regulated genes overlap!  Just very little differential expression here.  
SAR_vs_SA_MD_expression_dt = data.table(rbind(SAR_vs_SA_dt[padj < 0.05 & log2FoldChange > 0, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:50][!is.na(Gene)], SAR_vs_SA_dt[padj < 0.05 & log2FoldChange < 0, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:50][!is.na(Gene)]))

# load the chromatin accesibility data for each condition change 
A_vs_P = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_A_P.rds")
A_vs_P_dt = as.data.table(A_vs_P)
A_vs_P_dt$Gene = annotated_peak_atlas$SYMBOL
A_vs_P_dt$annotation = annotated_peak_atlas$annotation
A_vs_P_MD_accessability = A_vs_P_dt[Gene %in% A_vs_P_MD_expression_dt$Gene,]

SA_vs_A = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SA_A.rds")
SA_vs_A_dt = as.data.table(SA_vs_A)
SA_vs_A_dt$Gene = annotated_peak_atlas$SYMBOL
SA_vs_A_dt$annotation = annotated_peak_atlas$annotation
SA_vs_A_MD_accessability = SA_vs_A_dt[Gene %in% SA_vs_A_MD_expression_dt$Gene,]

SAR_vs_SA = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SAR_SA.rds")
SAR_vs_SA_dt = as.data.table(SAR_vs_SA)
SAR_vs_SA_dt$Gene = annotated_peak_atlas$SYMBOL
SAR_vs_SA_dt$annotation = annotated_peak_atlas$annotation
SAR_vs_SA_MD_accessability = SAR_vs_SA_dt[Gene %in% SAR_vs_SA_MD_expression_dt$Gene,]

# annotate as 'opening', 'closing', 'unchanged' for fill purposes later
A_vs_P_MD_accessability[, change := "not significant"]
A_vs_P_MD_accessability[padj < 0.05, change := "significant"]
A_vs_P_MD_accessability[log2FoldChange > 0, direction := "opening"]
A_vs_P_MD_accessability[log2FoldChange < 0, direction := "closing"]

SA_vs_A_MD_accessability[, change := "not significant"]
SA_vs_A_MD_accessability[padj < 0.05, change := "significant"]
SA_vs_A_MD_accessability[log2FoldChange > 0, direction := "opening"]
SA_vs_A_MD_accessability[log2FoldChange < 0, direction := "closing"]

SAR_vs_SA_MD_accessability[, change := "not significant"]
SAR_vs_SA_MD_accessability[padj < 0.05,change := "not significant"]
SAR_vs_SA_MD_accessability[log2FoldChange > 0, direction := "opening"] 
SAR_vs_SA_MD_accessability[log2FoldChange < 0, direction := "closing"]

# decorate by peaks gained / lost as shape (annotation), fill (opening / closing, vs unchanged), color (opening / closing)
# join the expression and accessibility data
A_vs_P_plot_data = merge(A_vs_P_MD_expression_dt, A_vs_P_MD_accessability[,.(Gene, annotation, change, direction)])
SA_vs_A_plot_data = merge(SA_vs_A_MD_expression_dt, SA_vs_A_MD_accessability[,.(Gene, annotation, change, direction)])
SAR_vs_SA_plot_data = merge(SAR_vs_SA_MD_expression_dt, SAR_vs_SA_MD_accessability[,.(Gene, annotation, change, direction)])

## going to have to do my own: geom_point, calculate my own y-axis value for each point by arb. adding to start-point
A_vs_P_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
A_vs_P_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
A_vs_P_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

SA_vs_A_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
SA_vs_A_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
SA_vs_A_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

SAR_vs_SA_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
SAR_vs_SA_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
SAR_vs_SA_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

diamond_A_vs_P = ggplot(data = A_vs_P_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(A_vs_P_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size=0.5, linetype='dotted') + 
  geom_point(aes(y=yval,shape=annotation, fill=change, color=direction)) + 
  ggtitle("A versus P: top 20 up-regulated & 20 down-regulated genes") + 
  guides(direction="none", fill="none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-15,15))

# Apparently the 'diamond' plots only have peak symbols on top of gene info
diamond_A_vs_P_top = ggplot(data = A_vs_P_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(A_vs_P_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size=0.5, linetype='dotted') + 
  geom_point(aes(y=top_yval, shape=annotation, color=direction)) + 
  ggtitle("A versus P: top 50 up-regulated & 50 down-regulated genes") + 
  guides(direction="none", fill="none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

diamond_SA_vs_A = ggplot(data = SA_vs_A_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SA_vs_A_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') + 
  geom_point(aes(y=yval,shape=annotation, fill=change, color=direction)) +
  ggtitle("SA versus A: top 50 up-regulated & 50 down-regulated genes") + 
  guides(fill="none", direction="none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-15,15))

diamond_SA_vs_A_top =  ggplot(data = SA_vs_A_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SA_vs_A_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') + 
  geom_point(aes(y=top_yval,shape=annotation, fill=change, color=direction)) +
  ggtitle("SA versus A: top 50 up-regulated & 50 down-regulated genes") + 
  guides(fill="none", direction="none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

diamond_SAR_vs_SA = ggplot(data = SAR_vs_SA_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SAR_vs_SA_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene, log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') +
  geom_point(aes(y=yval,shape=annotation, fill=change, color=direction)) + 
  ggtitle("SAR versus SA: top 50 up-regulated & 50 down-regulated genes") + 
  guides(size="none", direction="none", fill="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-15,15))

diamond_SAR_vs_SA_top = ggplot(data = SAR_vs_SA_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SAR_vs_SA_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene, log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') +
  geom_point(aes(y=top_yval,shape=annotation, fill=change, color=direction)) + 
  ggtitle("SAR versus SA: top 50 up-regulated & 50 down-regulated genes") + 
  guides(size="none", direction="none", fill="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))