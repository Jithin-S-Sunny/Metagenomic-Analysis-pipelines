# Metagenomics Analysis
Microbiome Data Analysis Pipeline
This R script is designed for the analysis of microbiome data, focusing on mealworm samples. It includes multiple steps such as quality checking, data transformation, visualization, and statistical analyses including alpha diversity, PERMANOVA, and differential abundance testing using ANCOM-BC2.

Prerequisites
Ensure you have the following R packages installed before running the script:
dplyr
ggplot2
reshape2
scales
vegan
ggrepel
zCompositions
compositions
ANCOMBC
phyloseq
TreeSummarizedExperiment
SummarizedExperiment
pairwiseAdonis
ape
stringr
You can install them using the following commands:
install.packages(c("dplyr", "ggplot2", "reshape2", "scales", "vegan", "ggrepel", "stringr"))
BiocManager::install(c("zCompositions", "compositions", "ANCOMBC", "phyloseq", "TreeSummarizedExperiment", "SummarizedExperiment", "pairwiseAdonis", "ape"))

Usage
Input Files
Abundance Data: The abundance data in CSV format, with samples as columns and genera as rows (default: mealworm_genus.csv).
Metadata File: Metadata for samples in CSV format (default: metadata_meal.csv).
User Inputs
Modify the following variables in the script based on your dataset:

input_file: Path to the abundance data file.
output_prefix: Prefix for output files (e.g., "mealworm").
metadata_file: Path to the metadata file.
columns_of_interest: Specify the columns (samples) you want to analyze. Set to NULL to analyze all samples.
Quality Checking
The script includes a function to check for unwanted characters and leading/trailing spaces in column names and values. This step ensures clean input data before proceeding.

Data Filtering
Genera with low abundance (less than 0.5% of the total abundance) are grouped as "Others". The filtered abundance table is saved as <output_prefix>_Filtered_Abundance_genus.csv.

Visualization
The script generates the following plots:

Stacked Bar Plot: Displays the relative abundance of genera across samples. Saved as <output_prefix>_stacked_bar_plot.png.
Bubble Plot: Visualizes the abundance of genera as bubbles of varying sizes. Saved as <output_prefix>_bubble_plot.png.
Alpha Diversity Metrics
Calculates Shannon, Simpson, inverse Simpson, richness, and evenness indices for each sample. Results are saved in <output_prefix>_alpha_diversity_metrics.csv and visualized in a boxplot (<output_prefix>_alpha_diversity_plot.png).

PCoA (Principal Coordinates Analysis)
The script performs PCoA using Bray-Curtis distances and displays both untransformed and CLR (Centered Log-Ratio) transformed PCoA plots:

<output_prefix>_pcoa_untransformed.png
<output_prefix>_pcoa_clr.png
Statistical Analysis
The script performs PERMANOVA for group comparisons using pairwise.adonis2 and visualizes R-squared values and Tukey HSD results. Outputs include:

PERMANOVA results: <output_prefix>_combined_results_genus.txt
R² plot: <output_prefix>_permanova_r2_plot.png
Tukey HSD plot: <output_prefix>_tukey_hsd_plot.png
Differential Abundance Testing (ANCOM-BC2)
The script runs ANCOM-BC2 to identify differentially abundant genera across groups. The results are saved as:

Pairwise results: <output_prefix>_res_pair.csv
Differentially abundant genera: <output_prefix>_differentially_abundant_genera.csv
Notes
Ensure that all input files have consistent formatting and matching sample IDs between the abundance and metadata files.
Adjust the parameters in the script based on your experimental design and requirements.
Troubleshooting
Data Import Issues: Verify that the file paths and names match your local files. Ensure there are no leading/trailing spaces in column names or unwanted characters.
Missing Packages: Install missing packages using install.packages() or BiocManager::install() for Bioconductor packages.
