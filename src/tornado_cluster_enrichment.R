require(data.table)
require(stringr)
require(clusterProfiler)
require(DESeq2)

# load the different clusters in each development path
base_dir = "/Users/zamparol/projects/AML_ATAC/results/peaks/tornado_cluster_bedfiles"
dev_path = "SAR_vs_P"
dds_peaks = readRDS("/Users/zamparol/projects/AML_ATAC/results/DESeq/objects/dds_object.rds")
sig_alpha = 0.05

setwd(paste(base_dir, dev_path, sep="/"))

cluster_dirs = dir(".", pattern="SAR_vs_P_", no..=TRUE)
get_cluster = function(filename){
  df = read.delim(file = paste0(filename), sep="\t", header=FALSE)
  colnames(df) = c("chrom", "start", "end")
 return(list(name=filename, dt=data.table(df)))
}

cluster_bed_list <- lapply(cluster_dirs, get_cluster)

# load the atlas
atlas = read.table("/Users/zamparol/projects/AML_ATAC/results/peaks/all_conditions_peak_atlas_annotated.csv", sep=",", header=TRUE)
atlas_dt = data.table(atlas)
setkeyv(atlas_dt, c("chrom", "start", "end"))

# load the results tables for each development path, to get the commensurate
# log2FC values for each peak

all_peaks <- results(dds_peaks, contrast=c("Condition", "P", "SAR"), alpha=sig_alpha, independentFiltering=TRUE, parallel=TRUE)
all_peaks_dt = as.data.table(all_peaks)
peak_names = str_split_fixed(rownames(all_peaks),"-", n=3)
all_peaks_dt$chrom = peak_names[,1]
all_peaks_dt$start = as.integer(peak_names[,2])
all_peaks_dt$end = as.integer(peak_names[,3])
setkeyv(all_peaks_dt, c("chrom", "start", "end"))
all_peaks_dt$ENSEMBL = atlas_dt$ENSEMBL

# drop peaks with no ENSEMBL ID.  Will not be useful for enrichment tasks.
most_peaks_dt = all_peaks_dt[!is.na(ENSEMBL),]
eIDs = bitr(most_peaks_dt$ENSEMBL, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eIDs_dt = data.table(eIDs)
setkey(eIDs_dt, ENSEMBL)
setkey(most_peaks_dt, ENSEMBL)
entrez_named_peaks = merge(most_peaks_dt, eIDs_dt, all.x=TRUE)
setkeyv(entrez_named_peaks, c("chrom", "start", "end"))

write_out_named_clusters = function(cluster){
  # get the peaks for this cluster
  dt = cluster$dt
  setkeyv(dt, c("chrom", "end"))
  cluster_named_peaks = atlas_dt[dt]
  
  named_cluster_filename = gsub(".bed",".named.cluster.csv", cluster$name)
  write.table(cluster_named_peaks$SYMBOL, named_cluster_filename, quote=FALSE, row.names = FALSE, col.names=FALSE)
}

### the bed files in the clusters is off by 1!!! why on earth off by one ?????
setkeyv(atlas_dt, c("chrom", "end"))
cluster_enrichment_results = lapply(cluster_bed_list, write_out_named_clusters)

analyse_cluster = function(cluster){
  # get the peaks for this cluster
  dt = cluster$dt
  setkeyv(dt, c("chrom", "start", "end"))
  cluster_named_peaks = entrez_named_peaks[dt]
  
  # get the rownames for these peaks, get the gene names associated,
  # translate to ENTREZ id, get the log2FC
  entrez_names = cluster_named_peaks$ENTREZID
  fc.values = cluster_named_peaks$log2FoldChange
  names(fc.values) = entrez_names
  # drop those loci without an ENSEMBL ID
  fc.values <- subset(fc.values, !is.na(names(fc.values)))
  fc.values = sort(fc.values, decreasing = TRUE)
  
  # KEGG enrichment
  cluster_kegg = enrichKEGG(gene = names(fc.values),
                  universe = entrez_named_peaks[, ENTREZID],
                  pvalueCutoff= 0.1,
                  pAdjustMethod = "BH")
  
  # GSEA
  cluster_gsea <- gseGO(geneList = fc.values,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                nPerm = 1000,
                keyType = "ENTREZID",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.2,
                verbose = FALSE)
  return(list(kegg=as.data.frame(cluster_kegg), gsea=as.data.frame(cluster_gsea)))
}

cluster_enrichment_results = lapply(cluster_bed_list, analyse_cluster)

# write out the KEGG, GSEA enrichment data.frames to their own files
write_out_enrichment_dfs <- function(beds, enrichments){
  kegg_filename = gsub(".bed",".kegg.enrichment.csv", beds$name)
  gsea_filename = gsub(".bed",".gsea.enrichment.csv", beds$name)
  
  write.csv(enrichments$kegg, kegg_filename, row.names = FALSE)
  write.csv(enrichments$gsea, gsea_filename, row.names = FALSE)
}

mapply(write_out_enrichment_dfs, cluster_bed_list, cluster_enrichment_results)
