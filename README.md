# Mallotus-marker-analysis

R script for evaluating DNA barcode marker performance across *Mallotus* species (Euphorbiaceae). It computes intra- and interspecific pairwise p-distances, tests for a barcode gap, measures overlap between distance distributions, and ranks single and combined markers by how well they separate species.

## Requirements

- R (>= 4.0)
- Packages: `ape`, `tidyverse`, `writexl`

```r
install.packages(c("ape", "tidyverse", "writexl"))
```

## Usage

1. Put your aligned FASTA files in a folder and set it as the working directory (`setwd()` or an RStudio project).
2. Fill in the file names in the `seqs` list at the top of the script, e.g. `rbcL = read.FASTA("rbcL.fas")`.
3. Check the species labelling (see below) so it matches your sequence names.
4. Run the script end to end.

### Sample naming convention

- Query samples: the first 4 characters of the sequence name are the species code (e.g. `Mthi01` becomes `Mthi`).
- Reference sequences: names in the form `Genus_species...` (e.g. `Mthi_ref`) keep their full name as their own label.

## What it does

1. Reads FASTA alignments for five single markers (_rbcL, matK, ITS, phyC, ncpGS_) and three combined marker sets (_ITS_phyC, matK_phyC, rbcL_matK_ITS_). Combined sets must be pre-concatenated alignments.
2. Computes raw pairwise p-distances with pairwise deletion.
3. Splits distances into intraspecific and interspecific values per species.
4. Summarises each marker: mean, SD, min/max, and barcode gap (mean inter − mean intra).
5. Runs a Wilcoxon rank-sum test (inter > intra) per marker and reports the two Mann–Whitney U values (`W_plus` = U for inter vs intra, `W_minus` = n_inter × n_intra − `W_plus`).
6. Overlap analysis: for each species, counts the interspecific distances strictly smaller than that species' maximum intraspecific distance. Values equal to the maximum are not counted. Results are reported per species and pooled per marker (total overlapping values / total interspecific values).
7. Pairwise Wilcoxon signed-rank tests between markers for intra, inter, and barcode gap. Intra and inter are matched on shared sequence pairs; barcode gap is matched by species.
8. Composite score: each marker scores one win per significant pairwise comparison in the correct direction (smaller intra, larger inter, larger gap). A tie (p ≥ 0.05) gives no points to either marker. With 8 marker sets the maximum score is 7 + 7 + 7 = 21.
9. Plots: distance boxplots, barcode gap, and overlap.

## Ranking

Markers are ordered by:

1. Fewest species with overlap
2. Fewest overlapping interspecific values
3. Largest barcode gap
4. Smallest mean intraspecific distance
5. Largest mean interspecific distance

The composite score and the rank-sum statistics are reported next to the rank as supporting evidence. They are not used for sorting.

## Output

An excel file `barcode_results.xlsx` with nine sheets:

| Sheet | Contents |
|---|---|
| Summary Statistics | Mean/SD/min/max intra & inter, barcode gap |
| Rank-Sum Inter vs Intra | Wilcoxon rank-sum test per marker |
| Overlap by Species | Per species × marker overlap counts |
| Overlap Summary | Pooled overlap % per marker |
| Pairwise Intra | Pairwise Wilcoxon results (smaller intra wins) |
| Pairwise Inter | Pairwise Wilcoxon results (larger inter wins) |
| Pairwise Barcode Gap | Pairwise Wilcoxon results (larger gap wins) |
| Ranking | Final rank with overlap, barcode gap, distance and test statistics, composite score |
| All Distances | Long-format distance values for plotting |

Plots: `boxplot_distances.svg`, `barcode_gap.svg`, `overlap_pct.svg`.

## Notes

- **Species missing from some markers:** _Mnan_ is absent from phyC and ncpGS, so it is excluded from the pairwise cross-marker tests. It is still included in the per-marker summary and overlap statistics. For your own data, change `exclude_sp` and the filter in `per_species_common`.
- **Non-independence:** pair distances that share a sequence are not independent, so p-values from the pair-level Wilcoxon tests are approximate. Treat the composite score as supporting evidence, not strict proof.

## Data

Sequence alignments are not included. FASTA inputs and species labels can be swapped for any dataset with the same structure.
