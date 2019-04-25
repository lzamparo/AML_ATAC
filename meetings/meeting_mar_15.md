Meeting 15 March
==================

Some confusion about whether the 'batch' covariate represents the days on which the runs were sequenced, or whether the libraries were prepared


Did since last time: 
====================
(done) Diamond plot in SAR vs P
	-  still looks as if ATAC & expression are not correlated
(done) Scatter plots contrasting ATAC changes versus expression changes 
	- A vs P / S vs P, SA vs A / SA vs S, SAR vs SA 
	- reinforces the diamond plot results; there is not a correlation between ATAC FC and expression FC


Still TODO from last time:
==========================

(1) Try clustering on tornado plots: 
	- First, super sized of plots: all conditions & peaks including both A & S
	- in tandem, look for the A-specific peaks from the differential accessibility data and do GSEA && pathway enrichment on these
    - Would like to figure out a way to do a track-style validation on the AvsP peaks that do not open further in SA (or diminish a bit)
	   - what's a good way to visualize the trajectory of these peaks?  Should I be looking at the library size corrected counts around the summits?

(2) Some clustering & visualizaton that captures the transient peaks along transition from P -> A -> SA -> SAR -> ... && P -> S -> SA -> ..


New TODO:
=========
(3) Conversely, want to identify regions that are differential in A vs P and are maintained in SAR vs P
(4) Bump all diff heatmap to 2kb window, experiment with differences in dynamic range of colour bar, "meta-peak" colour coded info 
	- beside the tornado plot?
	- have to figure out how to present it
	
(5) Re-work diamond plot: have color of open / close depend on the LFC of the ATAC seq peak
	- show all peaks associated w the gene, do not threshold on FC
	
(6) Re-do the scatter plot, but each point should be the median of peak FC.  Just do what Yuri does.
(7) Send .bed files from all-diff peaks to Han for HOMER motif analysis