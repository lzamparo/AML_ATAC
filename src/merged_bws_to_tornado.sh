#!/bin/bash

# run from cpumulticore to complete tornado plots

resultsdir="/data/leslie/zamparol/AML_ATAC/data/results/"
condition="all_diff_peaks"
bedfilesdir="/data/leslie/zamparol/AML_ATAC/data/tornado_cluster_bedfiles/$condition"

matrix_name=$(echo "$condition"_matrix.gz)
outfile_png=$(echo "$resultsdir/figures/$condition"_tornado_plot.png)
title=$(echo "Differentially Accessible Peaks for $condition")

echo "Computing matrix for plotting..."
pushd ~/projects/tools/atac-atlas/tools
bash compute_matrix_multibed.sh -b /data/leslie/zamparol/AML_ATAC -t tracks/deepTools -o $condition -r $bedfilesdir -m matrices/$condition -p 7
popd
echo "Done."
echo "Plotting heatmap..."
cd /data/leslie/zamparol/AML_ATAC/data/matrices/$condition
plotHeatmap -m $matrix_name -o $outfile_png --whatToShow "heatmap and colorbar" --plotTitle "$title" --verbose --refPointLabel "Peak center" --colorMap "Reds" -min 0 -max 100
echo "Done."
