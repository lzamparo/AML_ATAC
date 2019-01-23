Meeting Jan 15th
================


Tornado plot changes:
	[done] finer scale for tornados (10bp bins)
	[done] cap scale for tonado plots (re-scale so we see 0 - 100), and specify colour scheme
	- increase width of region surrounding peaks to 2kb from 1kb
	- look at 'meta-peak' signal for each cluster and each condition

Transient vs stable DA peaks:
	[done] median fold change of expression versus median fold changes versus scatter plot
	[done] C suggests looking at the set of differential peaks in A vs P (or S vs P) and track these peaks through each condition
	[done] C suggests a scatter plot of (S vs P logFC) versus (S vs A logFC) and (S vs P) versus (SA vs P) to see which sets of peaks are still open 
		when you get to SAR
	- Find individual examples of peaks, look at IGV track visualizations of each:
		[done] appearing in A vs P, but closing / unchanged in SA vs P (dt_A[A_vs_P_DA == "opening" & SA_vs_P_DA %in% c("unchanged", "closing")][order(-A_vs_P_log2FoldChange)])
		[done] closed in A vs P, but opening in SA vs P (dt_A[A_vs_P_DA %in% c("unchanged", "closing") & SA_vs_P_DA == "opening" ][order(A_vs_P_log2FoldChange)])
		- appearing in S vs P, but closing / unchanged in SA vs P (dt_S[S_vs_P_DA == "opening" & SA_vs_P_DA %in% c("unchanged", "closing")][order(-S_vs_P_log2FoldChange)])
		- closed in S vs P, but opening in SA vs P (dt_A[S_vs_P_DA %in% c("unchanged", "closing") & SA_vs_P_DA == "opening" ][order(-S_vs_P_log2FoldChange)])
		 

Diamond plot changes:
	[done] 50 up, 50 down
	- grade the color of the peaks by log2FC

Bar plots for up/down:
	- consolidate these onto one scale, plot both as faceted barplot


Overall need a clearer picture of what peak are contributing at which stage:

- look at sets of common open peaks through stages.  Start with changes from parental through to SA, SAR
