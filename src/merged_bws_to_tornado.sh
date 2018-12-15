# run from cpumulticore to complete tornado plots
#!/bin/bash

condition="P_S_SA_SAR_SARN_SARF"
bedfile=$(echo "/data/leslie/zamparol/AML_ATAC/data/results/differential_peak_lists/$condition"_diff_peaks.bed)

matrix_name=$(echo "$condition"_matrix.gz)
outfile=$(echo "../../results/$condition"_tornado_plot.png)
title=$(echo "Differentially Accessible Peaks Along $condition")


bash compute_matrix.sh -b /data/leslie/zamparol/AML_ATAC -t tracks/deepTools -o $condition -r $bedfile -m matrices/$condition -p 7
cd /data/leslie/zamparol/AML_ATAC/matrices/$condition
plotHeatmap -m $matrix_name -o $outfile --kmeans 4 --whatToShow "heatmap and colorbar" -T $title --verbose --refPointLabel "Peak center" --outFileSortedRegions SortedOut_test.bed