# Main R Script: Metagenomics, Antimicrobial Resistance, and Ecological Network Analysis

## Overview

This repository contains the R script used to perform the analyses
documented in the associated research study. The script brings together
data preparation, ARG assessment, read-level taxonomic profiling, and
MAG-based co-occurrence network analyses.

The main analysis file is:

-   `Main_R_Research_Analysis.R`

The script expects the input files to be stored in a directory named
`data/` located in the project working directory. The two network
summary tables are written to the working directory when the script is
run.

## Analysis workflow

The script is organized into the following sections:

1.  **Required libraries and custom functions**\
    Loads the R packages used in the analyses and defines functions for
    top-feature plots, alpha-diversity visualization, and Bray--Curtis
    PCoA.

2.  **Data loading and metadata preparation**\
    Imports the metadata, Kraken taxonomy reference, Bracken profiles,
    read-count data, RGI/Bowtie results, CheckM2 quality report, MAG
    coverage data, and ARG count matrix.

3.  **ARG audit and correlation analysis**\
    Calculates per-sample ARG richness and the sum of ARG counts, then
    evaluates their association using Spearman correlation.

4.  **ARG rarefaction sensitivity analysis**\
    Retains samples with at least 8,000 ARG reads and estimates rarefied
    ARG richness using 1,000 repeated multinomial draws per eligible
    sample. A Kruskal--Wallis test compares rarefied richness among
    cities.

5.  **Read-level metagenomics and taxonomic profiling**\
    Filters genus-level Bracken assignments to bacterial taxa,
    normalizes estimated reads per million bacterial reads, and
    generates alpha-diversity and Bray--Curtis PCoA plots.

6.  **MAG co-occurrence network analysis**\
    Calculates mean MAG coverage, retains MAGs with mean coverage
    greater than 0.5, estimates pairwise Pearson correlations within
    each target city, applies Benjamini--Hochberg adjustment, and
    constructs networks using the thresholds `r >= 0.70` and
    `FDR < 0.05`. Community detection is performed using the Louvain
    algorithm. Modularity is compared with a rewired-network null model
    using 1,000 repetitions.

7.  **Correlation-threshold sensitivity analysis**\
    Repeats network construction and community detection at correlation
    thresholds of 0.70, 0.80, and 0.90.

## Input files

Place the following files in the `data/` directory. File names must
match the names used in the R script.

  ------------------------------------------------------------------------------------------
  Input file                                             Role in the script
  ------------------------------------------------------ -----------------------------------
  `metadata.tsv`                                         Sample metadata, including sample
                                                         code, city, year, month, and month
                                                         label

  `k2_nt_20230502.inspect.species.krona.phylotree.txt`   Kraken taxonomy reference used to
                                                         identify bacterial taxa

  `bracken_phylum.tsv`                                   Phylum-level Bracken estimates

  `bracken_genus.tsv`                                    Genus-level Bracken estimates

  `bracken_species.tsv`                                  Species-level Bracken estimates
                                                         (loaded by the script; not used in
                                                         the current genus-level profiling
                                                         workflow)

  `bacterial_read_df.tsv`                                Total and bacterial read counts per
                                                         sample

  `rgi_bowtie.tsv`                                       RGI/Bowtie antimicrobial-resistance
                                                         mapping results

  `checkm2_quality_report.csv`                           CheckM2 MAG quality report; the
                                                         script reads this file using
                                                         tab-delimited input

  `sample_to_dereplicated_bin_coverage.tsv`              Per-sample MAG/contig coverage
                                                         information

  `ARG_counts_427_aligned.tsv`                           Sample-by-ARG count matrix used in
                                                         the ARG audit and rarefaction
                                                         analysis
  ------------------------------------------------------------------------------------------

## Software requirements

-   R (a version compatible with the package versions installed on the
    analysis system)
-   The following R packages:

``` r
install.packages(c(
  "tidyverse",
  "dplyr",
  "tibble",
  "janitor",
  "reshape2",
  "RColorBrewer",
  "paletteer",
  "ggpubr",
  "ggalt",
  "vegan",
  "ecodist",
  "igraph",
  "ggraph",
  "ggforce"
))
```

Install any package that is unavailable from the default CRAN repository
using its official distribution source. Package installation is a
one-time setup step and is not performed automatically by the analysis
script.

## How to run

1.  Arrange the project directory as shown below.
2.  Open R or RStudio.
3.  Set the working directory to the project root (the directory
    containing `Main_R_Research_Analysis.R` and `data/`).
4.  Run the script.

Example:

``` r
setwd("path/to/project")
source("Main_R_Research_Analysis.R")
```

The script uses relative file paths. Running it from a different working
directory may cause input-file errors.

## Output files

The script writes the following CSV files to the current working
directory:

-   `Supplementary_Table_S9_Network_Topology.csv` --- network topology
    and modularity summary for each target city.
-   `Network_Correlation_Threshold_Sensitivity.csv` --- network metrics
    calculated at the three correlation thresholds.

The diversity and ordination plotting functions return ggplot objects
for display in the active R plotting device. The current script does not
explicitly save these plots to image files.

## Reproducibility settings and analysis parameters

The script specifies the following settings:

  -----------------------------------------------------------------------
  Analysis component                  Setting
  ----------------------------------- -----------------------------------
  ARG identity filtering              Not applied in the current ARG
                                      count-matrix audit; the matrix is
                                      used as supplied

  Rarefaction target depth            8,000 ARG counts

  Rarefaction repetitions             1,000

  Rarefaction random seed             `20261013`

  Rarefaction group comparison        Kruskal--Wallis test

  ARG richness/read-count association Spearman correlation

  MAG mean-coverage filter            `> 0.5`

  Network correlation method          Pearson

  Primary network threshold           `r >= 0.70`

  Network significance threshold      Benjamini--Hochberg adjusted
                                      `FDR < 0.05`

  Community detection                 Louvain algorithm

  Modularity null-model repetitions   1,000

  Network random seed                 `20261102`

  Sensitivity thresholds              `0.70`, `0.80`, and `0.90`
  -----------------------------------------------------------------------

## Project directory structure

``` text
project/
├── Main_R_Research_Analysis.R
├── README.md
└── data/
    ├── metadata.tsv
    ├── k2_nt_20230502.inspect.species.krona.phylotree.txt
    ├── bracken_phylum.tsv
    ├── bracken_genus.tsv
    ├── bracken_species.tsv
    ├── bacterial_read_df.tsv
    ├── rgi_bowtie.tsv
    ├── checkm2_quality_report.csv
    ├── sample_to_dereplicated_bin_coverage.tsv
    └── ARG_counts_427_aligned.tsv
```

## Important notes

-   Input data are not included in this repository and must be supplied
    separately.
-   Preserve the input file names and tabular formats expected by the
    script.
-   The file `checkm2_quality_report.csv` is imported with a tab
    delimiter in the current script; confirm that the file is
    tab-delimited despite its `.csv` extension.
-   The script includes `bracken_phylum.tsv` and `bracken_species.tsv`
    as input imports, although the current profiling section performs
    its displayed diversity and ordination analyses at the genus level.
-   The ARG rarefaction procedure uses repeated multinomial sampling
    (`rmultinom`) as implemented in the analysis script. It should not
    be described as hypergeometric rarefaction.
-   The script has not been independently executed in a clean R
    environment as part of preparing this README. Successful execution
    depends on the availability and compatibility of the required
    packages and input files.
-   No statistical method, threshold, or analysis result should be
    changed without checking the study's finalized methods and results.

## Citation

If you use or adapt this script, cite the associated research article or
repository record once it is available. Add the final article citation,
DOI, or repository DOI here.

------------------------------------------------------------------------
