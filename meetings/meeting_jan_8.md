Meeting w Eirini & Tiansu & Christina:


TODOS:

- calculate all pairwise differential peak list (rainbow product of all comparisons differential peaks)
	- re-do the tornado plot analysis with this set, re-do the enrichments for the clusters found in this set
- re-do PCA plot only with differential peaks 
- check signs on diamond plot colors
- re-do library size normalization, just use peaks as per usual
- motif analysis in the peaks that open at each step (HOMER?)
	- Non-redundant set of motifs, one per TF (Hyunwoo + Sagar + Han have the list)
		- enrichment of motifs found in differntial peaks versus whole atlas
	- restrict to TFs expressed in the celltype

Priorities:

 [done] have a look at norm factors computed by flanks versus peaks 
 	- difference is modest, but present.
 [done] library size normalization by peaks, re-calculate the peaks that are differentially expressed in each condition
 [done] Re-produce all figures in the make_figues script, upload to drive
 	[done] including PCA figure with SARN / SARF 
 4) Re-do clustering of count matrices for all sets of peaks, generate bedfiles for same.
 	[done] cluster each of the three conditions
 	[done] run compute_matrix && plotHeatmap script for each of contrast sets
 	[done] do previous two steps for rainbow product peak set
 		[done] number of clusters is ambiguous from dendrogram.  Trying k = 11, 10, 9
 	[] re-run enrichment scripts for each
 		[] especially for different clusterings of all differentially accessible genes
 		[] try GeneMania, CytoScape E-Maps
 [done] Re-do diamond plots, just have all symbols above gene DE value, connect lines of DE plot
 
 6) motif analysis in the peaks that open at each step (FIMO)
	- Non-redundant set of motifs, one per TF (got list from Han)
		- enrichment of motifs found in differntial peaks versus whole atlas
	- restrict to TFs expressed in the celltype