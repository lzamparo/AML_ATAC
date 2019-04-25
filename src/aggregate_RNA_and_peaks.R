require(DESeq2)
require(data.table)
require(ggplot2)
require(latex2exp)
require(gtools)

# Read the RNA-seq DESeq results into tables

A_vs_P = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/A_vs_P.txt", sep="")
S_vs_P = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/S_vs_P.txt", sep="")
SA_vs_A = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SA_vs_A.txt", sep="")
SA_vs_S = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SA_vs_S.txt", sep="")
SAR_vs_SA = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SAR_vs_SA.txt", sep="")
SAR_vs_P = read.csv("~/projects/AML_ATAC/results/DESeq/gene_lists/SAR_vs_P.txt", sep="") 

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

SAR_vs_P_dt = as.data.table(SAR_vs_P)
SAR_vs_P_dt[, condition := "SAR_vs_P"]

# Load the annotated peak atlas
annotated_peak_atlas = data.table(read.csv("~/projects/AML_ATAC/results/peaks/all_conditions_peak_atlas_annotated.csv"))

# select the top 20 & bottom 20 DE genes for each condition
top_n = 20
A_vs_P_MD_expression_dt = data.table(rbind(A_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:top_n], A_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:top_n]))
S_vs_P_MD_expression_dt = data.table(rbind(S_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:top_n], S_vs_P_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:top_n]))

SA_vs_A_MD_expression_dt = data.table(rbind(SA_vs_A_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:top_n], SA_vs_A_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:top_n]))
SA_vs_S_MD_expression_dt = data.table(rbind(SA_vs_S_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:top_n], SA_vs_S_dt[padj < 0.05, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:top_n]))

# The top top_n, bottom top_n significant up regulated & down regulated genes overlap!  Just very little differential expression here.  
SAR_vs_SA_MD_expression_dt = data.table(rbind(SAR_vs_SA_dt[padj < 0.05 & log2FoldChange > 0, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:top_n][!is.na(Gene)], SAR_vs_SA_dt[padj < 0.05 & log2FoldChange < 0, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:top_n][!is.na(Gene)]))

SAR_vs_P_MD_expression_dt = data.table(rbind(SAR_vs_P_dt[padj < 0.05 & log2FoldChange > 0, .(Gene,log2FoldChange, padj)][order(-log2FoldChange)][1:top_n][!is.na(Gene)], SAR_vs_P_dt[padj < 0.05 & log2FoldChange < 0, .(Gene,log2FoldChange, padj)][order(log2FoldChange)][1:top_n][!is.na(Gene)]))


# load the chromatin accesibility data for each condition change 
A_vs_P = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_A_P.rds")
A_vs_P_dt = as.data.table(A_vs_P)
A_vs_P_dt$Gene = annotated_peak_atlas$SYMBOL
A_vs_P_dt$annotation = annotated_peak_atlas$annotation
A_vs_P_MD_accessability = A_vs_P_dt[Gene %in% A_vs_P_MD_expression_dt$Gene,]

S_vs_P = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_S_P.rds")
S_vs_P_dt = as.data.table(S_vs_P)
S_vs_P_dt$Gene = annotated_peak_atlas$SYMBOL
S_vs_P_dt$annotation = annotated_peak_atlas$annotation
S_vs_P_MD_accessability = S_vs_P_dt[Gene %in% S_vs_P_MD_expression_dt$Gene,]

SA_vs_A = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SA_A.rds")
SA_vs_A_dt = as.data.table(SA_vs_A)
SA_vs_A_dt$Gene = annotated_peak_atlas$SYMBOL
SA_vs_A_dt$annotation = annotated_peak_atlas$annotation
SA_vs_A_MD_accessability = SA_vs_A_dt[Gene %in% SA_vs_A_MD_expression_dt$Gene,]

SA_vs_S = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SA_S.rds")
SA_vs_S_dt = as.data.table(SA_vs_S)
SA_vs_S_dt$Gene = annotated_peak_atlas$SYMBOL
SA_vs_S_dt$annotation = annotated_peak_atlas$annotation
SA_vs_S_MD_accessability = SA_vs_S_dt[Gene %in% SA_vs_S_MD_expression_dt$Gene,]

SAR_vs_SA = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SAR_SA.rds")
SAR_vs_SA_dt = as.data.table(SAR_vs_SA)
SAR_vs_SA_dt$Gene = annotated_peak_atlas$SYMBOL
SAR_vs_SA_dt$annotation = annotated_peak_atlas$annotation
SAR_vs_SA_MD_accessability = SAR_vs_SA_dt[Gene %in% SAR_vs_SA_MD_expression_dt$Gene,]

SAR_vs_P = readRDS(file = "~/projects/AML_ATAC/results/DESeq/objects/res_SAR_P.rds")
SAR_vs_P_dt = as.data.table(SAR_vs_P)
SAR_vs_P_dt$Gene = annotated_peak_atlas$SYMBOL
SAR_vs_P_dt$annotation = annotated_peak_atlas$annotation
SAR_vs_P_MD_accessability = SAR_vs_P_dt[Gene %in% SAR_vs_P_MD_expression_dt$Gene,]

# annotate as 'opening', 'closing', 'unchanged' for fill purposes later
# calculate peak change magnitude as level (abs(log2FC) / max(abs(log2FC))
A_vs_P_MD_accessability[, change := "not significant"]
A_vs_P_MD_accessability[padj < 0.05, change := "significant"]
A_vs_P_MD_accessability[log2FoldChange > 0, direction := "opening"]
A_vs_P_MD_accessability[log2FoldChange < 0, direction := "closing"]
setnames(A_vs_P_MD_accessability, "log2FoldChange", "log2FoldChangeATAC")

S_vs_P_MD_accessability[, change := "not significant"]
S_vs_P_MD_accessability[padj < 0.05, change := "significant"]
S_vs_P_MD_accessability[log2FoldChange > 0, direction := "opening"]
S_vs_P_MD_accessability[log2FoldChange < 0, direction := "closing"]
setnames(S_vs_P_MD_accessability, "log2FoldChange", "log2FoldChangeATAC")

SA_vs_A_MD_accessability[, change := "not significant"]
SA_vs_A_MD_accessability[padj < 0.05, change := "significant"]
SA_vs_A_MD_accessability[log2FoldChange > 0, direction := "opening"]
SA_vs_A_MD_accessability[log2FoldChange < 0, direction := "closing"]
setnames(SA_vs_A_MD_accessability, "log2FoldChange", "log2FoldChangeATAC")

SA_vs_S_MD_accessability[, change := "not significant"]
SA_vs_S_MD_accessability[padj < 0.05, change := "significant"]
SA_vs_S_MD_accessability[log2FoldChange > 0, direction := "opening"]
SA_vs_S_MD_accessability[log2FoldChange < 0, direction := "closing"]
setnames(SA_vs_S_MD_accessability, "log2FoldChange", "log2FoldChangeATAC")

SAR_vs_SA_MD_accessability[, change := "not significant"]
SAR_vs_SA_MD_accessability[padj < 0.05,change := "not significant"]
SAR_vs_SA_MD_accessability[log2FoldChange > 0, direction := "opening"] 
SAR_vs_SA_MD_accessability[log2FoldChange < 0, direction := "closing"]
setnames(SAR_vs_SA_MD_accessability, "log2FoldChange", "log2FoldChangeATAC")

SAR_vs_P_MD_accessability[, change := "not significant"]
SAR_vs_P_MD_accessability[padj < 0.05,change := "not significant"]
SAR_vs_P_MD_accessability[log2FoldChange > 0, direction := "opening"] 
SAR_vs_P_MD_accessability[log2FoldChange < 0, direction := "closing"]
setnames(SAR_vs_P_MD_accessability, "log2FoldChange", "log2FoldChangeATAC")

# decorate by peaks gained / lost as fill (opening / closing), color (opening / closing), level (abs(log2FC) / max(abs(log2FC)))
# join the expression and accessibility data
A_vs_P_plot_data = merge(A_vs_P_MD_expression_dt, A_vs_P_MD_accessability[,.(Gene, annotation, change, direction, log2FoldChangeATAC)])
S_vs_P_plot_data = merge(S_vs_P_MD_expression_dt, S_vs_P_MD_accessability[,.(Gene, annotation, change, direction, log2FoldChangeATAC)])

SA_vs_A_plot_data = merge(SA_vs_A_MD_expression_dt, SA_vs_A_MD_accessability[,.(Gene, annotation, change, direction, log2FoldChangeATAC)])
SA_vs_S_plot_data = merge(SA_vs_S_MD_expression_dt, SA_vs_S_MD_accessability[,.(Gene, annotation, change, direction, log2FoldChangeATAC)])

SAR_vs_SA_plot_data = merge(SAR_vs_SA_MD_expression_dt, SAR_vs_SA_MD_accessability[,.(Gene, annotation, change, direction, log2FoldChangeATAC)])

SAR_vs_P_plot_data = merge(SAR_vs_P_MD_expression_dt, SAR_vs_P_MD_accessability[,.(Gene, annotation, change, direction, log2FoldChangeATAC)])

# for each *_data DT, need to quantize the ATAC log2FC values, create a custom palette, assign colours to the quantiles, plot:

# ggplot(mtcars) + 
#   geom_histogram(aes(factor(hp), fill=factor(hp))) + 
#   scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(colourCount)) +
#   theme(legend.position="bottom") +
#   guides(fill=guide_legend(nrow=2))

diamondPalette <- function(values, cut.low = 0.01, cut.high = 0.99,
                           type = "linear") {
  num.levels <- 20
  values.pos <- values[values >= 0]
  values.pos <- values.pos[values.pos > quantile(values.pos, cut.low) &
                             values.pos < quantile(values.pos, cut.high)]
  values.neg <- values[values < 0]
  values.neg <- values.neg[values.neg > quantile(values.neg, 1 - cut.high) &
                             values.neg < quantile(values.neg, 1 - cut.low)]
  
  breaks.pos <- seq(min(values.pos), max(values.pos),
                      length.out = num.levels)
  breaks.neg <- seq(min(values.neg), max(values.neg),
                      length.out = num.levels)
  red <- "#DC143C"
  green <- "#32CD32"
  light.grey <- "#F8F8F8"
  color.palette <-
    colorRampPalette(c(green, light.grey,red))(2 * num.levels + 1)
  names(color.palette) <- c(breaks.neg, 0, breaks.pos)
  color.palette
}


## going to have to do my own: geom_point, calculate my own y-axis value for each point by arb. adding to start-point
A_vs_P_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
A_vs_P_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
A_vs_P_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]
d_pal = diamondPalette(A_vs_P_plot_data$log2FoldChangeATAC)
A_vs_P_plot_data[, d_factor := cut(log2FoldChangeATAC, breaks = as.numeric(names(d_pal)))]

S_vs_P_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
S_vs_P_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
S_vs_P_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

SA_vs_A_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
SA_vs_A_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
SA_vs_A_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

SA_vs_S_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
SA_vs_S_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
SA_vs_S_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

SAR_vs_SA_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
SAR_vs_SA_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
SAR_vs_SA_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

SAR_vs_P_plot_data[direction == "opening", yval := log2FoldChange + 1:.N, by = Gene]
SAR_vs_P_plot_data[direction == "closing", yval := log2FoldChange - 1:.N, by = Gene]
SAR_vs_P_plot_data[, top_yval := log2FoldChange + 1:.N, by = Gene]

### Ok, finally the diamond plots

# Apparently the 'diamond' plots only have peak symbols on top of gene info
diamond_A_vs_P_top = ggplot(data = A_vs_P_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(A_vs_P_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size=0.5, linetype='dotted') + 
  geom_point(aes(y=top_yval, color=d_factor, fill=d_factor)) + 
  ggtitle("A versus P: top 25 up-regulated & 25 down-regulated genes") + 
  scale_color_manual(values = d_pal) +
  scale_fill_manual(values = d_pal) +
  guides(direction="none", fill="none", color = "none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$"))

diamond_S_vs_P_top = ggplot(data = S_vs_P_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(S_vs_P_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size=0.5, linetype='dotted') + 
  geom_point(aes(y=top_yval, shape=annotation, color=direction)) + 
  ggtitle("S versus P: top 25 up-regulated & 25 down-regulated genes") + 
  guides(direction="none", fill="none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

diamond_SA_vs_A_top =  ggplot(data = SA_vs_A_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SA_vs_A_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') + 
  geom_point(aes(y=top_yval,shape=annotation, fill=change, color=direction)) +
  ggtitle("SA versus A: top 25 up-regulated & 25 down-regulated genes") + 
  guides(fill="none", direction="none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

diamond_SA_vs_S_top =  ggplot(data = SA_vs_S_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SA_vs_S_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene,log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') + 
  geom_point(aes(y=top_yval,shape=annotation, fill=change, color=direction)) +
  ggtitle("SA versus S: top 25 up-regulated & 25 down-regulated genes") + 
  guides(fill="none", direction="none", size="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

diamond_SAR_vs_SA_top = ggplot(data = SAR_vs_SA_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SAR_vs_SA_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene, log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') +
  geom_point(aes(y=top_yval,shape=annotation, fill=change, color=direction)) + 
  ggtitle("SAR versus SA: top 25 up-regulated & 25 down-regulated genes") + 
  guides(size="none", direction="none", fill="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))

diamond_SAR_vs_P_top = ggplot(data = SAR_vs_P_plot_data, aes(x=reorder(Gene,log2FoldChange))) + 
  geom_line(data = unique(SAR_vs_P_plot_data[,.(Gene, log2FoldChange)]), aes(x=reorder(Gene, log2FoldChange), y=log2FoldChange, group=1), size = 0.5, linetype='dotted') +
  geom_point(aes(y=top_yval,shape=annotation, fill=change, color=direction)) + 
  ggtitle("SAR versus P: top 25 up-regulated & 25 down-regulated genes") + 
  guides(size="none", direction="none", fill="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 7)) +
  xlab("Gene (symbol)") + 
  ylab(TeX("$\\log_{2}(FC)$")) + 
  ylim(c(-10,20))
