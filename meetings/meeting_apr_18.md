Meeting 11 April
================

- Han presented GATA ChIP overlap with ATAC-seq peaks results
	- Gist: the effect of expression changes for genes proximal to peaks with GATA1/2 motif-bearing peaks disappears when we use only those 
	with overlap of GATA2 ChIP in k562

- Eirini thinks two things can be presented as findings:
	- A has all these transient changes, most do not persist into later states (SA, SAR, ...)
	- Is there a consequential / quantifiable difference to order of mutation (P -> A -> SA) versus (P -> S -> SA)
	- Examine differences in stable versus transient changes

- Lee finds two bugs:
	- the scaling code for .bws is using the *inverse* of the DESeq value; DESeq divides by the calculated SF, while 
	deepTools multiplies
	- the scaling factors used were calculated from *flanknig* regions, rather than peak regions.  

So I need to re-calculate the scaling factors, and re-generate the bigWig tracks with the proper normalization applied.

Done since last time:
=====================

1. Re-did the scatter plot for SAR vs P, with each point should be colored by one-sided Wilcoxon of FCs for peaks assoc w the gene versus FCs in whole atlas


Still TODO from last time:
==========================

1. Re-work diamond plot: 

	- have color of open / close depend on the LFC of the ATAC seq peak
	- show all peaks associated w the gene, do not threshold on FC
		- impose cutoff on log2FC for gene expression of 2 (~~also get DEseq.rds file from Han for updated results~~)??


2. Bump all diff heatmap to 2kb window, experiment with differences in dynamic range of colour bar, "meta-peak" colour coded info 

	- beside the tornado plot?
	- have to figure out how to present it
	- have to look at the distribution of dynamic range for changes: should I be 'capping' some values to reveal more subtle changes?


Ideas to try to make sense of the contradictions in the data:
=============================================================

Try clustering on tornado plots: 

	- First, super sized of plots: all conditions & peaks including both A & S
	- in tandem, look for the A-specific peaks from the differential accessibility data and do GSEA && pathway enrichment on these
    - Would like to figure out a way to do a track-style validation on the AvsP peaks that do not open further in SA (or diminish a bit)
	   - what's a good way to visualize the trajectory of these peaks?  Should I be looking at the library size corrected counts around the summits?

