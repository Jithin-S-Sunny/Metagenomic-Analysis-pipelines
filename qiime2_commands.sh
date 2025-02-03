#!/bin/bash

# Activate conda environment
conda activate qiime2-amplicon-2024.2

# Convert files made in Windows and check for whitespaces
tr ',' '\t' < metadata.csv > metadata.tsv
awk -F',' 'NF != N {print NR, $0}' metadata.tsv
sed 's/"//g' metadata.tsv > metadata_cleaned.tsv

# Import FASTQ files and create a QIIME2 artifact. The input.tsv file will be in the format as given in the file input.tsv in this directory 
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path input.tsv --output-path sequences.qza --input-format PairedEndFastqManifestPhred33V2

# Summarize sequences and visualise all the .qzv results at https://view.qiime2.org/
qiime demux summarize --i-data sequences.qza --o-visualization demux.qzv

# DADA2: Quality filtering, denoising, dereplication, merging, and chimera detection
qiime dada2 denoise-paired --i-demultiplexed-seqs sequences.qza --p-trunc-len-f 240 --p-trunc-len-r 240 --output-dir dada2 --verbose

# Visualize denoising stats
qiime metadata tabulate --m-input-file dada2/denoising_stats.qza --o-visualization dada2/denoising-stats.qzv

# To build a phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences dada2/representative_sequences.qza --output-dir tree

qiime empress tree-plot --i-tree tree/rooted_tree.qza --o-visualization tree/empress.qzv

# Summarize table for understading depth and read frequency
qiime feature-table summarize --i-table 16S_bacteria_dada2/table.qza --o-visualization table.qzv

# Perform diversity analysis (set sampling depth based on the quartile values observed in previous step)
qiime diversity core-metrics-phylogenetic --i-table dada2/table.qza --i-phylogeny tree/rooted_tree.qza --p-sampling-depth 3600 --m-metadata-file metadata.tsv --output-dir diversity

qiime diversity alpha-group-significance --i-alpha-diversity diversity/shannon_vector.qza --m-metadata-file metadata.tsv --o-visualization diversity/alpha_groups.qzv

# Export alpha diversity
qiime tools export --input-path diversity/shannon_vector.qza --output-path exported_alpha_diversity

# Run PERMANOVA test (Adonis)
qiime diversity adonis --i-distance-matrix diversity/weighted_unifrac_distance_matrix.qza --m-metadata-file metadata.tsv --p-formula "site" --p-n-jobs 2 --o-visualization diversity/permanova.qzv

# Taxonomic classification using the updated classifier. 
qiime feature-classifier classify-sklearn --i-reads dada2/representative_sequences.qza --i-classifier silva-138-99-nb-classifier.qza --p-n-jobs 2 --o-classification taxa.qza

# Taxa barplot
qiime taxa barplot --i-table 16S_dada2/table.qza --i-taxonomy 16S_taxa/classification.qza --m-metadata-file metadata.tsv --o-visualization 16S_taxa/taxa_barplot.qzv

#Extract ASV table
qiime tools export --input-path table.qza --output-path bacteria_exported_feature_table

biom convert --input-fp bacteria_exported_feature_table/feature-table.biom --output-fp bacteria_exported_feature_table/feature_table.tsv --to-tsv
