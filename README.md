# Mallotus-marker-analysis

R script for evaluating DNA barcode marker performance across _Mallotus_ species (Euphorbiaceae). Computes intra- and interspecific pairwise p-distances, tests for a barcode gap, measures overlap between distance distributions, and ranks single and combined markers by how reliably they separate species.

# Required R packages
- R (>= 4.0)
- Packages: `ape`, `tidyverse`, `writexl`
- Install with: 
```r
install.packages(c("ape", "tidyverse", "writexl"))
```

# Usage
1. Place your aligned FASTA files in the working directory.
2. Set the file paths in the seqs list at the top of the script.
3. Adjust species_lookup to match your own species labels.
4. Run the script end to end.

# What it does
1. Reads FASTA alignments for five single markers (rbcL, matK, ITS, phyC, ncpGS) and three combined marker sets. With the combined marker sets, the script expects them as pre-concatenated alignments.
2. Computes raw pairwise P-distances with pairwise deletion.
3. Splits distances into intraspecific and interspecific values per species.
4. Summarises each marker: mean, SD, min/max, and barcode gap (mean inter − mean intra).
5. Wilcoxon rank-sum tests (inter > intra) per marker.
6. Overlap analysis: for each species, the fraction of interspecific distances strictly smaller than that species’ maximum intraspecific distance, reported per species and pooled per marker.
7. Species absent from some markers (here Mnan, missing in phyC and ncpGS) are excluded from the pairwise cross-marker tests. Change exclude_sp and the filter in per_species_common for your own data. They are still included in the per-marker summary and overlap statistics.
8. Pairwise Wilcoxon tests between markers (intra, inter, and barcode gap), matched on shared sequence pairs/species. Pair distances sharing a sequence are not independent, so p-values from the pair-level Wilcoxon tests are approximate. Use the composite score as a relative ranking.
10. Composite ranking: each marker scores one win per significant pairwise comparison in the correct direction. A tie (p ≥ 0.05) gives no points to either marker. With 8 marker sets, the maximum score is 7 + 7 + 7 = 21. Ties in the ranking are broken by smaller pooled overlap %, then larger barcode gap.
11. Boxplots of distance distributions, barcode gap, and overlap.

# Output
An excel file (barcode_results.xlsx) with nine sheets:
1. Sheet	Contents
2. Summary Statistics	Mean/SD/min/max intra & inter, barcode gap
3. Rank-Sum Inter vs Intra	Wilcoxon test per marker
4. Overlap by Species	Per species × marker overlap counts
5. Overlap Summary	Pooled overlap % per marker
6. Pairwise Intra / Inter / Barcode Gap	Pairwise Wilcoxon results
7. Ranking	Composite score and final rank
8. All Distances	Long-format distance values for plotting

# Data
Sequence alignments are not included. The script is written so that FASTA inputs and species labels can be swapped for any dataset with the same structure.
