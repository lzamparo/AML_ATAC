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

# select the top 20 & bottom 20 DE genes for each condition
A_vs_P_MD_expression_dt = data.table(rbind(A_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:20], A_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:20]))
SA_vs_A_MD_expression_dt = data.table(rbind(SA_vs_A_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:20], SA_vs_A_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:20]))
SAR_vs_SA_MD_expression_dt = data.table(rbind(SAR_vs_SA_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:20], SAR_vs_SA_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:20]))

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

# plot the A vs P genes: line plot (smoothed?) (y) log2FC vs (x) gene names aes(x=reorder(class,-amount,sum),y=amount
# decorate by peaks gained / lost as shape (annotation), fill (opening / closing, vs unchanged), color (opening / closing)
# join the expression and accessibility data
A_vs_P_plot_data = merge(A_vs_P_MD_expression_dt, A_vs_P_MD_accessability[,.(Gene, annotation, change, direction)])


## going to have to do my own: geom_point, calculate my own y-axis value for each point by arb. adding to start-point



# geom_dotplot(stackgroups = TRUE, binwidth = 1, method = "histodot")
diamond_A_vs_P = ggplot(data = A_vs_P_plot_data, aes(x=reorder(Gene,log2FoldChange),y=peaks_y)) + 
  geom_point(data = A_vs_P_plot_data[direction ==  "opening",], aes(color="red", shape=annotation, fill=change)) + 
  geom_point(data = A_vs_P_plot_data[direction ==  "closing",], aes(color="blue", shape=annotation, fill=change)) +
  ggtitle("Top 20 Upregulated & 20 Downregulated genes in A versus P") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,10))



