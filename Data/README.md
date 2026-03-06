# Data Directory

Created by: Michael Sieler  
Last updated: 2026-03-06

This directory contains all input data used by the workshop, plus generated
dataset-specific files needed to connect sequencing runs to sample metadata.

## Structure

- `Metadata/`
  - `GutMicrobiota_BeeTracking_metadata.csv` - full metadata table for sequencing runs and sample context
  - `GutMicrobiota_BeeTracking_metadata__TutorialSubset.csv` - workshop subset used for this tutorial workflow
- `PRJNA792398__06.03.2026/`
  - `FASTQs/`
    - `fastq_raw/` - downloaded raw `.fastq.gz` files from the BioProject
    - `fastq_symlinks__TutorialSubset/` - metadata-matched symlinks used by the tutorial subset
  - `samplesheet__TutorialSubset.csv` - DADA2 input table with sample IDs and paired FASTQ paths
  - `srr_to_sampleid_mapping__TutorialSubset.txt` - plain-text SRR to sample ID mapping report
  - `metadata_sample_ids.normalized.txt` - normalized sample IDs used during mapping
  - `run_accessions.txt` - list of run accessions for the BioProject
  - `logs/` - download and processing logs
- `Databases/`
  - `silva_nr99_v138.2_toSpecies_trainset.fa.gz` - SILVA training set used for taxonomy assignment
  - (optional) `silva_species_assignment_v138.2.fa.gz` for species-level exact matching
- `DADA2/`
  - intermediate `.rds` checkpoint objects from Part 1 of the tutorial
- `MicrobiomeAnalysis/`
  - intermediate `.rds` checkpoint objects from Part 2 of the tutorial

## Notes

- `samplesheet__TutorialSubset.csv` and `fastq_symlinks__TutorialSubset/` define which samples are included in the tutorial.
- Raw FASTQ counts can differ from symlink counts because not all runs map to workshop metadata.
