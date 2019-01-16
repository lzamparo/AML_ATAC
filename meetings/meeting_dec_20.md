Meeting w Eirini & Tiansu & Han:


TODOS:

- Repeat tornado plot but starting from those peaks differential between SAR vs P
	- four panel plot for S or A, SARN or SARF
	- number of changes for SARN, SARF (part of bar plot figure)

- revisit Tornado plots, more types of clustering that are more structured than k-means (import matrix into R, cluster by Ward)
	- excel sheet of nearest gene to each peak in the clusters

- remove non-coding RNA detritus, reannotate the atlas with just RefSeq genes

- Diamond plots for A vs P, SA vs A, SAR vs SA
	- top 20, bottom 20 significant diff'ble expressed genes 

- motif analysis in the peaks that open at each step (HOMER?)
	- Non-redundant set of motifs, one per TF (Hyunwoo + Sagar + Han have the list)
		- enrichment of motifs found in differntial peaks versus whole atlas
	- restrict to TFs expressed in the celltype
	


Priorities:
	- Tornado plots x 3
		- with GSEA on clusters
		- annotation breakdown for each cluster
	- Diamond plots
	- motif analysis of sequence for peak


My ordering:

[done] reannotate atlas 
[done] regenerate tornado plots (union of stage-wise differential peaks (S first & A first) [done], set of peaks differential between SAR vs P [done])
[done] characterize clusters: GSEA on gene associated with peaks, KEGG enrichment with peaks. [SAR vs P done]
[done] diamond plots for: A vs P, S vs P, SA vs A, SA vs P, SAR vs SA, SARN vs SAR, SARF vs SAR
5) motif analysis of sequence for each peak