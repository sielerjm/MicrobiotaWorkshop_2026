# Code Directory

Created by: Michael Sieler  
Last updated: 2026-03-06

This directory contains the tutorial and helper scripts used to download,
organize, and analyze the workshop microbiome dataset.

## Structure

- `Microbiome_Workshop_Tutorial.Rmd`
  - Main end-to-end tutorial (Part 1: DADA2 processing; Part 2: microbiome analysis)
  - Writes figures and tables to `../Results/`
  - Stores checkpoint `.rds` files in `../Data/DADA2` and `../Data/MicrobiomeAnalysis`
- `Scripts/`
  - `Microbiome_Workshop_Tutorial.R`
    - Script export of the R Markdown tutorial for line-by-line execution/debugging
  - `download_bioproject_fastq.sh`
    - Downloads FASTQ files for a BioProject
    - Supports test and full runs
    - Organizes output into dataset folders and logs
  - `map_srr_to_sampleid_and_build_samplesheet.R`
    - Maps SRR accessions to metadata sample IDs
    - Creates FASTQ symlinks and `samplesheet__TutorialSubset.csv`
    - Generates a human-readable mapping report

## Typical workflow

1. Run `Scripts/download_bioproject_fastq.sh` to get raw FASTQs.
2. Run `Scripts/map_srr_to_sampleid_and_build_samplesheet.R` to build symlinks and `samplesheet__TutorialSubset.csv`.
3. Knit `Microbiome_Workshop_Tutorial.Rmd` for processing and downstream analyses.
