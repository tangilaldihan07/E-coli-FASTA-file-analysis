# FASTA-based Sequence Analysis of *E. coli*

This repository contains a beginner-level bioinformatics project focused on analyzing DNA coding sequences (CDS) of *Escherichia coli* K-12 MG1655 using FASTA files and Python.

## Overview
The project performs basic descriptive sequence analysis on a subset of *E. coli* genes to understand nucleotide composition and sequence patterns.

## Analyses Performed
- Gene length calculation
- GC content analysis
- 2-mer and 3-mer frequency profiling
- Codon usage analysis
- Comparison with genome-wide GC content
- Principal Component Analysis (PCA) based on 2-mer composition
- Visualization using bar plots, heatmaps, and PCA plots

## Data Source
- RefSeq CDS FASTA (`cds_from_genomic.fna`)
- Assembly: GCF_000005845.2  
- Organism: *Escherichia coli* K-12 MG1655

## Tools Used
- Python
- pandas
- matplotlib
- scikit-learn

## Notes
This project is exploratory and descriptive, intended for learning FASTA handling, feature extraction, and visualization in bioinformatics. No predictive modeling or functional inference is performed.

## Repository Contents
- Python scripts for analysis
- Generated plots and figures
- FASTA input file (or instructions to download it)

