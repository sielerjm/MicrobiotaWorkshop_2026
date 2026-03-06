# MicrobiotaWorkshop_2026

Created by: Joanito Liberti and Michael Sieler  
Last updated: 2026-03-06

## Project Overview

This repository contains the data, scripts, and analysis workflow for a
microbiome methods tutorial focused on honeybee gut 16S rRNA amplicon data.
The project is designed to support both:

- a reproducible, end-to-end workshop pipeline (raw FASTQ to ecological inference), and
- transparent documentation suitable for scientific review and reuse.

The tutorial implements a two-part workflow:

1. **DADA2 processing**: quality control, denoising, chimera removal, taxonomy assignment, and phyloseq construction.
2. **Microbiome analysis**: composition, alpha/beta diversity, and differential abundance using metadata-linked samples.

## Scientific Context

The analysis uses sequencing runs mapped to metadata-defined biological samples,
including explicit handling of controls (kit blanks, water blanks, mock
community, inoculum controls) for quality control and contamination-aware steps.
Taxonomy assignment is performed with SILVA reference databases.

## Data Source and Attribution

Parts of the data and analysis context in this workshop are derived from the
open repository:

- https://github.com/JoanitoLiberti/The-gut-microbiota-affects-the-social-network-of-honeybees/tree/master

If you reuse this material, please cite the original publication:

- Liberti J., Kay T., Quinn A., Kesner L., Frank E.T., Cabirol A., Richardson
  T.O., Engel P., Keller L. The gut microbiota affects the social network of
  honeybees. Nature Ecology & Evolution 6, 1471–1479 (2022).
  https://doi.org/10.1038/s41559-022-01840-w

## Repository Structure

- `Code/`
  - main tutorial notebook (`Microbiome_Workshop_Tutorial.Rmd`)
  - helper scripts for download and sample mapping in `Code/Scripts/`
  - see `Code/README.md` for code-level details
- `Data/`
  - raw and mapped sequencing inputs
  - metadata files
  - reference databases (for taxonomy assignment)
  - intermediate checkpoint `.rds` objects for both DADA2 and downstream analysis
  - see `Data/README.md` for data organization details
- `Results/`
  - generated figures and tables organized by workflow part (`DADA2`, `MicrobiomeAnalysis`)
  - see `Results/README.md` for output details

## Reproducibility Notes

- The tutorial writes intermediate checkpoint objects to `Data/` so long-running
  steps can be resumed.
- Figures and tables are written to `Results/` for clean separation of analysis
  artifacts from source data and scripts.
- The workflow is metadata-aware and sample ID consistent from FASTQ mapping
  through downstream statistical analyses.

## Intended Use

This repository is intended for:

- workshop instruction and training in microbiome analysis workflows,
- reproducible analysis handoff between collaborators,
- supplementary methodological transparency for manuscript preparation and review.

## License

This repository includes a general educational-use license in `LICENSE`. It is
maintained for teaching and training purposes in a University of Geneva course.
Please retain attribution when reusing or adapting course materials.
