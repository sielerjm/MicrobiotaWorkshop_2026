#!/usr/bin/env Rscript

# Script filename: map_srr_to_sampleid_and_build_samplesheet.R
# Description: Map SRR FASTQ files to metadata Sample_ID values using ENA run metadata,
#              create sample-ID-based FASTQ symlinks, and export a DADA2-ready samplesheet.csv.
# Created by: Michael Sieler
# Last updated: 06.03.2026
# Expected input:
#   - A metadata CSV/TXT containing a Sample_ID column
#   - A directory with downloaded SRR FASTQ files named SRR*_1.fastq.gz and SRR*_2.fastq.gz
#   - A BioProject ID for ENA run mapping lookup
# Expected output:
#   - Symlinked FASTQ files named SRR__SampleID_1.fastq.gz and SRR__SampleID_2.fastq.gz
#   - samplesheet.csv containing sample_id, forward, reverse, run_accession, and metadata columns
#   - a mapping text report summarizing matched/excluded runs and per-sample file links

args <- commandArgs(trailingOnly = TRUE)

resolve_latest_dataset_dir <- function(data_root = "Data", bioproject = "PRJNA792398") {
  candidates <- list.dirs(data_root, full.names = FALSE, recursive = FALSE)
  pattern <- paste0("^", bioproject, "__[0-9]{2}\\.[0-9]{2}\\.[0-9]{4}$")
  candidates <- candidates[grepl(pattern, candidates)]

  if (length(candidates) == 0) {
    return(paste0(bioproject, "__05.03.2026"))
  }

  candidate_paths <- file.path(data_root, candidates)
  candidates[which.max(file.info(candidate_paths)$mtime)]
}

parse_args <- function(args_vec) {
  latest_dataset_dir <- resolve_latest_dataset_dir("Data", "PRJNA792398")

  # Defaults target the current workshop repository layout.
  opts <- list(
    bioproject = "PRJNA792398",
    metadata = "Data/Metadata/GutMicrobiota_BeeTracking_metadata__TutorialSubset.csv",
    raw_dir = file.path("Data", latest_dataset_dir, "FASTQs", "fastq_raw"),
    out_dir = file.path("Data", latest_dataset_dir, "FASTQs", "fastq_symlinks__TutorialSubset"),
    samplesheet = file.path("Data", latest_dataset_dir, "samplesheet__TutorialSubset.csv"),
    mapping_txt = file.path("Data", latest_dataset_dir, "srr_to_sampleid_mapping__TutorialSubset.txt")
  )

  if (length(args_vec) == 0) {
    return(opts)
  }

  i <- 1
  while (i <= length(args_vec)) {
    key <- args_vec[[i]]
    if (key %in% c("-h", "--help")) {
      cat(
        "Usage:\n",
        "  Rscript Scripts/map_srr_to_sampleid_and_build_samplesheet.R [options]\n\n",
        "Options:\n",
        "  --bioproject ID     BioProject accession (default: PRJNA792398)\n",
        "  --metadata PATH     Metadata CSV with Sample_ID column\n",
        "  --raw-dir PATH      Folder containing SRR FASTQ files\n",
        "  --out-dir PATH      Folder for Sample_ID-based symlinks\n",
        "  --samplesheet PATH  Output samplesheet CSV path\n",
        "  --mapping-txt PATH  Output mapping report text file path\n",
        "  -h, --help          Show this help message\n\n",
        "Example:\n",
        "  Rscript Scripts/map_srr_to_sampleid_and_build_samplesheet.R \\\n",
        "    --bioproject PRJNA792398 \\\n",
        "    --metadata Data/Metadata/GutMicrobiota_BeeTracking_metadata__TutorialSubset.csv \\\n",
        "    --raw-dir Data/PRJNA792398__DD.MM.YYYY/FASTQs/fastq_raw \\\n",
        "    --out-dir Data/PRJNA792398__DD.MM.YYYY/FASTQs/fastq_symlinks__TutorialSubset \\\n",
        "    --samplesheet Data/PRJNA792398__DD.MM.YYYY/samplesheet__TutorialSubset.csv \\\n",
        "    --mapping-txt Data/PRJNA792398__DD.MM.YYYY/srr_to_sampleid_mapping__TutorialSubset.txt\n",
        sep = ""
      )
      quit(status = 0)
    }

    if (!startsWith(key, "--") || i == length(args_vec)) {
      stop("Invalid argument syntax near: ", key, call. = FALSE)
    }

    value <- args_vec[[i + 1]]
    if (key == "--bioproject") opts$bioproject <- value
    if (key == "--metadata") opts$metadata <- value
    if (key == "--raw-dir") opts$raw_dir <- value
    if (key == "--out-dir") opts$out_dir <- value
    if (key == "--samplesheet") opts$samplesheet <- value
    if (key == "--mapping-txt") opts$mapping_txt <- value
    i <- i + 2
  }

  opts
}

opts <- parse_args(args)

if (!file.exists(opts$metadata)) {
  stop("Metadata file not found: ", opts$metadata, call. = FALSE)
}
if (!dir.exists(opts$raw_dir)) {
  stop("Raw FASTQ directory not found: ", opts$raw_dir, call. = FALSE)
}

# Resolve key directories to absolute paths for robust symlink creation.
raw_dir_abs <- normalizePath(opts$raw_dir, mustWork = TRUE)
out_dir_abs <- normalizePath(opts$out_dir, mustWork = FALSE)

# Required package for tidy data wrangling and I/O.
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  stop("Package 'tidyverse' is required. Install with install.packages('tidyverse').", call. = FALSE)
}
if (!requireNamespace("magrittr", quietly = TRUE)) {
  stop("Package 'magrittr' is required for %>% pipes.", call. = FALSE)
}
`%>%` <- magrittr::`%>%`

sanitize_for_filename <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+", "-", x)
  x <- gsub("-{2,}", "-", x)
  x <- gsub("(^-+)|(-+$)", "", x)
  x[x == ""] <- "NA"
  x
}

normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  is_numeric_only <- grepl("^[0-9]+$", x)
  # Normalize purely numeric IDs so "01" and "1" are treated as the same sample.
  x[is_numeric_only] <- as.character(as.integer(x[is_numeric_only]))
  x
}

read_metadata_table <- function(path) {
  is_tsv_like <- grepl("\\.(txt|tsv)$", tolower(path))
  delim <- if (is_tsv_like) "\t" else ","

  tbl <- readr::read_delim(
    file = path,
    delim = delim,
    col_types = readr::cols(
      Sample_ID = readr::col_character(),
      .default = readr::col_guess()
    ),
    show_col_types = FALSE
  )

  if (!"Sample_ID" %in% colnames(tbl)) {
    stop("Metadata must contain a 'Sample_ID' column.", call. = FALSE)
  }

  tbl <- tbl %>%
    dplyr::mutate(
      Sample_ID = trimws(as.character(.data$Sample_ID)),
      Sample_ID_norm = normalize_sample_id(.data$Sample_ID)
    )

  if (is_tsv_like) {
    metadata_csv <- sub("\\.(txt|tsv)$", ".csv", path, ignore.case = TRUE)
    readr::write_csv(tbl %>% dplyr::select(-"Sample_ID_norm"), metadata_csv)
    message("Converted metadata TSV/TXT to CSV: ", metadata_csv)
  }

  tbl
}

metadata_tbl <- read_metadata_table(opts$metadata)

if (!"Sample_ID" %in% colnames(metadata_tbl)) {
  stop("Metadata must contain a 'Sample_ID' column.", call. = FALSE)
}

if (any(duplicated(metadata_tbl$Sample_ID_norm))) {
  dup_ids <- metadata_tbl$Sample_ID[duplicated(metadata_tbl$Sample_ID_norm)]
  stop(
    "Duplicate normalized Sample_ID values in metadata. Example duplicate(s): ",
    paste(utils::head(unique(dup_ids), 8), collapse = ", "),
    call. = FALSE
  )
}

# ENA run mapping provides run_accession <-> sample_alias linkage for this BioProject.
ena_url <- paste0(
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", opts$bioproject,
  "&result=read_run&fields=run_accession,sample_alias,sample_accession,library_layout,fastq_ftp&format=tsv"
)

run_map_tbl <- readr::read_tsv(ena_url, show_col_types = FALSE) %>%
  dplyr::transmute(
    run_accession = as.character(.data$run_accession),
    Sample_ID_run = trimws(as.character(.data$sample_alias)),
    Sample_ID_norm = normalize_sample_id(.data$sample_alias),
    sample_accession = as.character(.data$sample_accession),
    library_layout = as.character(.data$library_layout),
    fastq_ftp = as.character(.data$fastq_ftp)
  ) %>%
  dplyr::filter(!is.na(.data$run_accession), .data$run_accession != "")

total_bioproject_runs <- nrow(run_map_tbl)

# Keep only entries that appear in the workshop metadata Sample_ID list.
mapped_tbl <- run_map_tbl %>%
  dplyr::inner_join(
    metadata_tbl,
    by = "Sample_ID_norm"
  ) %>%
  dplyr::distinct(.data$run_accession, .data$Sample_ID, .keep_all = TRUE)

matched_metadata_runs <- nrow(mapped_tbl)
excluded_runs <- total_bioproject_runs - matched_metadata_runs

if (nrow(mapped_tbl) == 0) {
  stop("No ENA run mappings matched metadata Sample_ID values.", call. = FALSE)
}

if (any(duplicated(mapped_tbl$Sample_ID))) {
  dup_ids <- mapped_tbl$Sample_ID[duplicated(mapped_tbl$Sample_ID)]
  stop(
    "Duplicate Sample_ID mappings detected in ENA map. Example duplicate(s): ",
    paste(utils::head(unique(dup_ids), 5), collapse = ", "),
    call. = FALSE
  )
}

# Resolve forward/reverse SRR FASTQ paths and validate file existence.
mapped_tbl <- mapped_tbl %>%
  dplyr::mutate(
    forward_srr = file.path(raw_dir_abs, paste0(.data$run_accession, "_1.fastq.gz")),
    reverse_srr = file.path(raw_dir_abs, paste0(.data$run_accession, "_2.fastq.gz")),
    forward_exists = base::file.exists(.data$forward_srr),
    reverse_exists = base::file.exists(.data$reverse_srr)
  )

missing_tbl <- mapped_tbl %>%
  dplyr::filter(!.data$forward_exists | !.data$reverse_exists)

if (nrow(missing_tbl) > 0) {
  stop(
    "Missing FASTQ file(s) for ", nrow(missing_tbl), " sample(s). Example run(s): ",
    paste(utils::head(missing_tbl$run_accession, 8), collapse = ", "),
    call. = FALSE
  )
}

dir.create(out_dir_abs, recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing symlinks/files from prior naming schemes to avoid stale entries.
old_links <- list.files(out_dir_abs, pattern = "\\.fastq\\.gz$", full.names = TRUE)
if (length(old_links) > 0) {
  unlink(old_links, force = TRUE)
}

# Create symlinks using a compact naming scheme: SRR__SampleID_1/_2.fastq.gz
link_tbl <- mapped_tbl %>%
  dplyr::mutate(
    file_stub = paste(
      sanitize_for_filename(.data$run_accession),
      sanitize_for_filename(.data$Sample_ID),
      sep = "__"
    )
  ) %>%
  dplyr::mutate(
    # Keep samplesheet paths as user-provided paths (relative-friendly),
    # and store absolute paths only for symlink creation.
    forward = file.path(opts$out_dir, paste0(.data$file_stub, "_1.fastq.gz")),
    reverse = file.path(opts$out_dir, paste0(.data$file_stub, "_2.fastq.gz")),
    forward_abs = file.path(out_dir_abs, paste0(.data$file_stub, "_1.fastq.gz")),
    reverse_abs = file.path(out_dir_abs, paste0(.data$file_stub, "_2.fastq.gz"))
  )

for (i in seq_len(nrow(link_tbl))) {
  row <- link_tbl[i, ]

  # Remove pre-existing files/symlinks so reruns can safely recreate links.
  if (file.exists(row$forward_abs) || nzchar(Sys.readlink(row$forward_abs))) {
    unlink(row$forward_abs, force = TRUE)
  }
  if (file.exists(row$reverse_abs) || nzchar(Sys.readlink(row$reverse_abs))) {
    unlink(row$reverse_abs, force = TRUE)
  }

  ok_f <- file.symlink(from = row$forward_srr, to = row$forward_abs)
  ok_r <- file.symlink(from = row$reverse_srr, to = row$reverse_abs)
  if (!isTRUE(ok_f) || !isTRUE(ok_r)) {
    stop("Failed to create symlink(s) for Sample_ID: ", row$Sample_ID, call. = FALSE)
  }
}

# Build final samplesheet with tutorial-friendly columns + metadata.
samplesheet_tbl <- link_tbl %>%
  dplyr::select(
    "Sample_ID", "run_accession", "sample_accession", "forward", "reverse",
    dplyr::everything()
  ) %>%
  dplyr::select(
    -"Sample_ID_norm", -"Sample_ID_run", -"library_layout", -"fastq_ftp",
    -"forward_srr", -"reverse_srr", -"forward_exists", -"reverse_exists",
    -"file_stub", -"forward_abs", -"reverse_abs"
  ) %>%
  dplyr::rename(sample_id = "Sample_ID") %>%
  dplyr::arrange(.data$sample_id)

readr::write_csv(samplesheet_tbl, opts$samplesheet)

# Write a plain-text mapping report for quick inspection/documentation.
mapping_tbl <- link_tbl %>%
  dplyr::transmute(
    sample_id = .data$Sample_ID,
    run_accession = .data$run_accession,
    sample_accession = .data$sample_accession,
    raw_r1 = .data$forward_srr,
    raw_r2 = .data$reverse_srr,
    symlink_r1 = .data$forward,
    symlink_r2 = .data$reverse
  ) %>%
  dplyr::arrange(.data$sample_id)

mapping_lines <- c(
  paste0("BioProject: ", opts$bioproject),
  paste0("Metadata file: ", opts$metadata),
  paste0("Raw FASTQ directory: ", opts$raw_dir),
  paste0("Symlink directory: ", opts$out_dir),
  paste0("Samplesheet: ", opts$samplesheet),
  paste0("Generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("Total BioProject runs: ", total_bioproject_runs),
  paste0("Runs matched to metadata Sample_ID: ", matched_metadata_runs),
  paste0("Runs excluded (not in metadata): ", excluded_runs),
  "",
  "Columns:",
  "sample_id\trun_accession\tsample_accession\traw_r1\traw_r2\tsymlink_r1\tsymlink_r2",
  apply(mapping_tbl, 1, function(row) paste(row, collapse = "\t"))
)

writeLines(mapping_lines, con = opts$mapping_txt)

message("Mapping complete.")
message("BioProject: ", opts$bioproject)
message("Total BioProject runs: ", total_bioproject_runs)
message("Runs matched to metadata Sample_ID: ", matched_metadata_runs)
message("Runs excluded (not in metadata): ", excluded_runs)
message("Symlinked FASTQ directory: ", opts$out_dir)
message("Samplesheet: ", opts$samplesheet)
message("Mapping report: ", opts$mapping_txt)
