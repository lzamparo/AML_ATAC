Meeting 21 March
==================

Some confusion about whether the 'batch' covariate represents the days on which the runs were sequenced, or whether the libraries were prepared

Done since last time:
=====================

- median ATAC per gene versus expression change for A/S vs P, SA vs A/S, SAR vs SA
- sent .bed files to Han for HOMER motif analysis 

Still TODO from last time:
==========================

1. Try clustering on tornado plots: 

	- First, super sized of plots: all conditions & peaks including both A & S
	- in tandem, look for the A-specific peaks from the differential accessibility data and do GSEA && pathway enrichment on these
    - Would like to figure out a way to do a track-style validation on the AvsP peaks that do not open further in SA (or diminish a bit)
	   - what's a good way to visualize the trajectory of these peaks?  Should I be looking at the library size corrected counts around the summits?

2. Some clustering & visualizaton that captures the transient peaks along transition from P -> A -> SA -> SAR -> ... && P -> S -> SA -> ..

3. Conversely, want to identify regions that are differential in A vs P and are maintained in SAR vs P

4. Bump all diff heatmap to 2kb window, experiment with differences in dynamic range of colour bar, "meta-peak" colour coded info 

	- beside the tornado plot?
	- have to figure out how to present it
	- have to look at the distribution of dynamic range for changes: should I be 'capping' some values to reveal more subtle changes?
	
5. Re-work diamond plot: have color of open / close depend on the LFC of the ATAC seq peak

	- show all peaks associated w the gene, do not threshold on FC
		- impose log2FC for gene expression of 2 (~~also get DEseq.rds file from Han for updated results~~)??
	
6. Re-do the scatter plot, but each point should be the median of peak FC
		(Clarification: color by one-sided Wilcoxon of FCs for peaks assoc w the gene versus FCs in whole atlas for same comparison)

	- C thinks plots look weird.  Many changes in ATAC without changes in expression.
		- Diagnostic: check genes that change in expression a lot, make IGV figures for them
	- Report Spearman correlation
	- Add SAR vs P comparison

New TODO:
=========

- all things from last week
- A vs P FC versus SAR vs P FC to ID transient versus persistent changes 