Meeting 11 April
================


Done since last time:
=====================

1. Re-do the scatter plot, but each point should be colored by one-sided Wilcoxon of FCs for peaks assoc w the gene versus FCs in whole atlas

	- C thinks plots look weird.  Many changes in ATAC without changes in expression.
		- ~~Diagnostic: check genes that change in expression a lot, make IGV figures for them~~
	- ~~Report Spearman correlation~~
	- ~~**Add SAR vs P comparison**~~

2. Re-work diamond plot: 

	- have color of open / close depend on the LFC of the ATAC seq peak
	- show all peaks associated w the gene, do not threshold on FC
		- impose cutoff on log2FC for gene expression of 2 (~~also get DEseq.rds file from Han for updated results~~)??


3. Bump all diff heatmap to 2kb window, experiment with differences in dynamic range of colour bar, "meta-peak" colour coded info 

	- beside the tornado plot?
	- have to figure out how to present it
	- have to look at the distribution of dynamic range for changes: should I be 'capping' some values to reveal more subtle changes?

4. Try clustering on tornado plots: 

	- First, super sized of plots: all conditions & peaks including both A & S
	- in tandem, look for the A-specific peaks from the differential accessibility data and do GSEA && pathway enrichment on these
    - Would like to figure out a way to do a track-style validation on the AvsP peaks that do not open further in SA (or diminish a bit)
	   - what's a good way to visualize the trajectory of these peaks?  Should I be looking at the library size corrected counts around the summits?

5. Some clustering & visualizaton that captures the transient peaks along transition from P -> A -> SA -> SAR -> ... && P -> S -> SA -> ..

6. Conversely, want to identify regions that are differential in A vs P and are maintained in SAR vs P

Still TODO from last time:
==========================