#!/usr/bin/env bash
set -euo pipefail

# Script: download_bioproject_fastq.sh
# Description: Download FASTQ files for SRA runs linked to an NCBI BioProject ID,
#              optionally restricted to Sample_ID values in a metadata table.
# Created by: Michael Sieler
# Last updated: 06.03.2026
# Expected input:
#   - BioProject ID (e.g., PRJNA792398)
#   - Optional flags for metadata-based sample subsetting, test mode,
#     output directory, and thread count
# Expected output:
#   - FASTQ(.gz) files in Data/<RUN_DIR>/FASTQs/fastq_raw
#   - run_accessions.txt, logs, and supplementary files in Data/<RUN_DIR>

usage() {
  cat <<'EOF'
Usage:
  download_bioproject_fastq.sh --bioproject PRJNA792398 [options]

Required:
  -p, --bioproject ID      NCBI BioProject ID (example: PRJNA792398)

Options:
      --test               Download only a small test subset (default: first 1 run)
  -n, --n-samples N        Number of runs in test mode (default: 1)
      --metadata PATH      Optional metadata CSV/TSV/TXT with Sample_ID column;
                           limits downloads to matching Sample_ID/sample_alias entries
  -o, --output-root DIR    Root output folder (default: ../../Data relative to this script)
  -t, --threads N          Threads for fasterq-dump and pigz (default: 4)
  -h, --help               Show this help message
Notes:
  - Test mode replaces same-day TEST directory for a clean rerun.
  - Full mode reuses same-day project directory and downloads only missing files.
  - This script stores all downloaded runs in FASTQs/fastq_raw.
  - Metadata-matched symlinks (FASTQs/fastq_symlinks) are created by the mapping R script and may contain fewer files.

Examples:
  # Test on 1 run (default test size)
  bash Scripts/download_bioproject_fastq.sh -p PRJNA792398 --test

  # Full download
  bash Scripts/download_bioproject_fastq.sh -p PRJNA792398

  # Test run with explicit output root
  bash Scripts/download_bioproject_fastq.sh \
    -p PRJNA792398 \
    --test \
    -o "../Projects/MicrobiotaWorkshop_2026/Data"
  # Output directories created:
  # ../Projects/MicrobiotaWorkshop_2026/Data/TEST__DD.MM.YYYY
  # ../Projects/MicrobiotaWorkshop_2026/Data/TEST__DD.MM.YYYY/FASTQs/fastq_raw

  # Full run with explicit output root and custom thread count
  bash Scripts/download_bioproject_fastq.sh \
    -p PRJNA792398 \
    -o "../Projects/MicrobiotaWorkshop_2026/Data" \
    -t 8
  # Output directories created:
  # ../Projects/MicrobiotaWorkshop_2026/Data/PRJNA792398__DD.MM.YYYY
  # ../Projects/MicrobiotaWorkshop_2026/Data/PRJNA792398__DD.MM.YYYY/FASTQs/fastq_raw
EOF
}

check_dependencies() {
  # Validate required and optional external tools before downloads begin.
  local required_cmds=(
    awk
    curl
    gzip
    tee
  )
  local missing=()
  local cmd=""

  for cmd in "${required_cmds[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      missing+=("$cmd")
    fi
  done

  if (( ${#missing[@]} > 0 )); then
    echo "ERROR: Missing required software dependency/dependencies:" >&2
    printf '  - %s\n' "${missing[@]}" >&2
    echo "Install the missing tools, then re-run this script." >&2
    exit 1
  fi

  # Optional toolsets used for fallback paths.
  HAVE_EDIRECT=0
  HAVE_SRA_TOOLS=0

  if command -v esearch >/dev/null 2>&1 && command -v efetch >/dev/null 2>&1; then
    HAVE_EDIRECT=1
  fi

  if command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1; then
    HAVE_SRA_TOOLS=1
  fi
}

fetch_ena_field() {
  local accession="$1"
  local fields="$2"
  # Return empty output on request errors so caller can decide fallback behavior.
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=${fields}&format=tsv" 2>/dev/null || true
}

# Default argument values (can be overridden by CLI flags).
BIOPROJECT=""
TEST_MODE=0
TEST_N=1
THREADS=4
METADATA_PATH=""
# Resolve script location so defaults work from any current directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_ROOT="${SCRIPT_DIR}/../../Data"

# Parse command-line arguments.
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--bioproject)
      BIOPROJECT="${2:-}"
      shift 2
      ;;
    --test)
      TEST_MODE=1
      shift
      ;;
    -n|--n-samples)
      TEST_N="${2:-}"
      shift 2
      ;;
    --metadata)
      METADATA_PATH="${2:-}"
      shift 2
      ;;
    -o|--output-root)
      OUTPUT_ROOT="${2:-}"
      shift 2
      ;;
    -t|--threads)
      THREADS="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# Validate required and numeric arguments.
if [[ -z "$BIOPROJECT" ]]; then
  echo "ERROR: --bioproject is required." >&2
  usage
  exit 1
fi

if [[ ! "$TEST_N" =~ ^[0-9]+$ ]] || (( TEST_N < 1 )); then
  echo "ERROR: --n-samples must be a positive integer." >&2
  exit 1
fi

if [[ ! "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS < 1 )); then
  echo "ERROR: --threads must be a positive integer." >&2
  exit 1
fi

if [[ -n "$METADATA_PATH" && ! -f "$METADATA_PATH" ]]; then
  echo "ERROR: --metadata file not found: $METADATA_PATH" >&2
  exit 1
fi

normalize_sample_id() {
  local id="$1"
  id="${id//$'\r'/}"
  id="${id#\"}"
  id="${id%\"}"
  id="$(printf '%s' "$id" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')"

  # Normalize common workshop naming variants before matching.
  id="${id#Bee_}"
  id="${id#Bee}"
  if [[ "$id" =~ ^X[0-9]+$ ]]; then
    id="${id#X}"
  fi
  if [[ "$id" =~ ^[0-9]+$ ]]; then
    id="$(printf '%s' "$id" | sed -E 's/^0+//')"
    [[ -n "$id" ]] || id="0"
  fi
  printf '%s' "$id"
}

load_metadata_ids() {
  local metadata_file="$1"
  local out_file="$2"
  local delim=","
  local sample_col=""
  local id_count=0

  if [[ "$metadata_file" =~ \.(txt|tsv)$ ]]; then
    delim=$'\t'
  fi

  sample_col="$(
    awk -F"$delim" '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          gsub(/\r/, "", $i)
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
          gsub(/^"|"$/, "", $i)
          if ($i == "Sample_ID") {
            print i
            exit
          }
        }
      }
    ' "$metadata_file"
  )"

  if [[ -z "$sample_col" ]]; then
    echo "ERROR: --metadata file must contain a 'Sample_ID' column: $metadata_file" >&2
    exit 1
  fi

  : > "$out_file"
  while IFS= read -r raw_id; do
    norm_id="$(normalize_sample_id "$raw_id")"
    [[ -n "$norm_id" ]] && printf '%s\n' "$norm_id" >> "$out_file"
  done < <(awk -F"$delim" -v c="$sample_col" 'NR > 1 {print $c}' "$metadata_file")

  # Keep unique IDs only.
  if [[ -s "$out_file" ]]; then
    sort -u "$out_file" -o "$out_file"
    id_count="$(wc -l < "$out_file" | tr -d '[:space:]')"
  fi

  if [[ "$id_count" == "0" ]]; then
    echo "ERROR: --metadata contained no non-empty Sample_ID values: $metadata_file" >&2
    exit 1
  fi

  echo "Loaded ${id_count} unique normalized Sample_ID values from metadata."
}

# Build destination folder name based on mode and current date.
DATE_TAG="$(date '+%d.%m.%Y')"
RUN_TS="$(date '+%Y%m%d__%H%M%S')"
if (( TEST_MODE == 1 )); then
  RUN_DIR_NAME="TEST__${DATE_TAG}"
else
  RUN_DIR_NAME="${BIOPROJECT}__${DATE_TAG}"
fi

# Create final project folder under Data and dedicated FASTQ subfolders.
PROJECT_DIR="${OUTPUT_ROOT%/}/${RUN_DIR_NAME}"
FASTQ_ROOT_DIR="${PROJECT_DIR}/FASTQs"
FASTQ_DIR="${FASTQ_ROOT_DIR}/fastq_raw"
SRA_CACHE_DIR="${PROJECT_DIR}/_sra_cache"
LOG_DIR="${PROJECT_DIR}/logs"
LOG_FILE="${LOG_DIR}/${RUN_DIR_NAME}__${RUN_TS}.log"
METADATA_ID_FILE="${PROJECT_DIR}/metadata_sample_ids.normalized.txt"

# For test runs, always allow a clean same-day rerun by replacing the folder.
# For full runs, reuse existing folder to support resuming only missing files.
if [[ -d "$PROJECT_DIR" ]]; then
  if (( TEST_MODE == 1 )); then
    rm -rf "$PROJECT_DIR"
  else
    echo "Detected existing full-run directory; resuming missing files only."
  fi
fi

# Ensure root output exists (create Data if missing), then create run folders.
mkdir -p "$OUTPUT_ROOT"
mkdir -p "$PROJECT_DIR" "$FASTQ_ROOT_DIR" "$FASTQ_DIR" "$SRA_CACHE_DIR" "$LOG_DIR"

# Capture all script messages/errors into a dated log while still printing to terminal.
exec > >(tee -a "$LOG_FILE") 2>&1

# Check software dependencies after logging is configured so failures are logged too.
check_dependencies

echo "BioProject: $BIOPROJECT"
echo "Mode: $([[ "$TEST_MODE" -eq 1 ]] && echo "TEST (${TEST_N} runs)" || echo "FULL")"
echo "Project directory: $PROJECT_DIR"
echo "FASTQ root directory: $FASTQ_ROOT_DIR"
echo "FASTQ raw directory: $FASTQ_DIR"
echo "Log file: $LOG_FILE"
echo "FASTQ output format: compressed (.fastq.gz)"
if [[ -n "$METADATA_PATH" ]]; then
  echo "Metadata subset file: $METADATA_PATH"
fi
if (( HAVE_EDIRECT == 0 )); then
  echo "WARNING: EDirect tools (esearch/efetch) not found; NCBI run-list fallback is unavailable."
fi
if (( HAVE_SRA_TOOLS == 0 )); then
  echo "WARNING: SRA tools (prefetch/fasterq-dump) not found; SRA fallback is unavailable."
fi
echo
echo "Collecting run accessions (ENA first, then NCBI fallback)..."

# First try ENA Portal API for stable run accession lookup by BioProject.
# If ENA does not return SRR accessions, fall back to NCBI EDirect.
RUNS=()
if [[ -n "$METADATA_PATH" ]]; then
  load_metadata_ids "$METADATA_PATH" "$METADATA_ID_FILE"
  ENA_RUN_ALIAS_TSV="$(fetch_ena_field "$BIOPROJECT" "run_accession,sample_alias")"
  if [[ -z "$ENA_RUN_ALIAS_TSV" ]]; then
    echo "ERROR: Could not retrieve ENA run/sample_alias mapping for metadata subset filtering." >&2
    exit 1
  fi

  while IFS=$'\t' read -r run sample_alias; do
    [[ -n "$run" ]] || continue
    sample_alias_norm="$(normalize_sample_id "$sample_alias")"
    if [[ -n "$sample_alias_norm" ]] && grep -Fxq "$sample_alias_norm" "$METADATA_ID_FILE"; then
      RUNS+=("$run")
    fi
  done < <(printf '%s\n' "$ENA_RUN_ALIAS_TSV" | awk 'NR > 1 && $1 ~ /^SRR/ {print $1 "\t" $2}')
else
  while IFS= read -r run; do
    [[ -n "$run" ]] && RUNS+=("$run")
  done < <(
    fetch_ena_field "$BIOPROJECT" "run_accession" \
      | awk 'NR > 1 && $1 ~ /^SRR/ {print $1}'
  )
fi

if (( ${#RUNS[@]} == 0 )); then
  if [[ -n "$METADATA_PATH" ]]; then
    echo "ERROR: No SRR runs matched metadata Sample_ID values for ${BIOPROJECT}." >&2
    exit 1
  fi
  if (( HAVE_EDIRECT == 0 )); then
    echo "ERROR: ENA run lookup returned no results, and EDirect fallback is unavailable." >&2
    echo "Install esearch/efetch (Entrez Direct) to enable NCBI fallback." >&2
    exit 1
  fi
  while IFS= read -r run; do
    [[ -n "$run" ]] && RUNS+=("$run")
  done < <(
    esearch -db sra -query "${BIOPROJECT}[BioProject]" \
      | efetch -format runinfo \
      | awk -F',' 'NR > 1 && $1 ~ /^SRR/ {print $1}'
  )
fi

if (( ${#RUNS[@]} == 0 )); then
  echo "ERROR: No SRR run accessions found for ${BIOPROJECT}." >&2
  exit 1
fi

# In test mode keep only the first N runs; otherwise use all runs.
if (( TEST_MODE == 1 )); then
  RUNS_TO_DOWNLOAD=()
  count=0
  for run in "${RUNS[@]}"; do
    RUNS_TO_DOWNLOAD+=("$run")
    ((count += 1))
    if (( count >= TEST_N )); then
      break
    fi
  done
else
  # In full mode, detect already completed runs and only queue missing runs/files.
  RUNS_TO_DOWNLOAD=()
  SKIPPED_COMPLETE=0
  for run in "${RUNS[@]}"; do
    run_complete=0

    ENA_FASTQ_FIELD="$(fetch_ena_field "$run" "fastq_ftp" | awk 'NR == 2 {print $2}')"

    if [[ -n "$ENA_FASTQ_FIELD" && "$ENA_FASTQ_FIELD" != "-" ]]; then
      run_complete=1
      for fastq_path in ${ENA_FASTQ_FIELD//;/ }; do
        [[ -n "$fastq_path" ]] || continue
        fastq_file="${FASTQ_DIR}/$(basename "$fastq_path")"
        if [[ ! -f "$fastq_file" ]]; then
          run_complete=0
          break
        fi
      done
    else
      # Heuristic for runs without ENA fastq_ftp metadata.
      # Consider complete if single-end OR paired-end compressed files are present.
      if [[ -f "${FASTQ_DIR}/${run}.fastq.gz" || ( -f "${FASTQ_DIR}/${run}_1.fastq.gz" && -f "${FASTQ_DIR}/${run}_2.fastq.gz" ) ]]; then
        run_complete=1
      fi
    fi

    if (( run_complete == 1 )); then
      ((SKIPPED_COMPLETE += 1))
    else
      RUNS_TO_DOWNLOAD+=("$run")
    fi
  done
  echo "Full mode resume check: ${SKIPPED_COMPLETE} run(s) already complete, ${#RUNS_TO_DOWNLOAD[@]} run(s) missing."
fi

if (( ${#RUNS_TO_DOWNLOAD[@]} == 0 )); then
  echo "All runs already present for this full-mode directory. Nothing to download."
  echo "Project directory: $PROJECT_DIR"
  echo "FASTQ root directory: $FASTQ_ROOT_DIR"
  echo "FASTQ raw directory: $FASTQ_DIR"
  echo "Log file: $LOG_FILE"
  exit 0
fi

# Save a manifest so users can track exactly which runs were requested.
printf '%s\n' "${RUNS_TO_DOWNLOAD[@]}" > "${PROJECT_DIR}/run_accessions.txt"
echo "Selected ${#RUNS_TO_DOWNLOAD[@]} run(s). Manifest: ${PROJECT_DIR}/run_accessions.txt"

# Keep a list of any runs that fail during download or conversion.
FAILED_RUNS=()

# Download and convert each run independently.
for run in "${RUNS_TO_DOWNLOAD[@]}"; do
  echo
  # Try direct FASTQ download from ENA first (faster and avoids SRA conversion).
  ENA_FASTQ_FIELD=""
  ENA_FASTQ_FIELD="$(fetch_ena_field "$run" "fastq_ftp" | awk 'NR == 2 {print $2}')"

  if [[ -n "$ENA_FASTQ_FIELD" && "$ENA_FASTQ_FIELD" != "-" ]]; then
    echo "[$run] Downloading FASTQ from ENA..."
    run_failed=0
    # ENA returns one or two FASTQ paths separated by semicolons.
    for fastq_path in ${ENA_FASTQ_FIELD//;/ }; do
      [[ -n "$fastq_path" ]] || continue
      fastq_url="https://${fastq_path#ftp://}"
      fastq_file="${FASTQ_DIR}/$(basename "$fastq_path")"
      if [[ -f "$fastq_file" ]]; then
        echo "[$run] Found existing file, skipping: $(basename "$fastq_file")"
        continue
      fi
      if ! curl -fL --retry 3 --retry-delay 2 -o "$fastq_file" "$fastq_url"; then
        echo "[$run] ERROR: ENA FASTQ download failed: $fastq_url" >&2
        run_failed=1
        break
      fi
    done
    if (( run_failed == 0 )); then
      continue
    fi
    echo "[$run] Falling back to prefetch + fasterq-dump..."
  fi

  if (( HAVE_SRA_TOOLS == 0 )); then
    echo "[$run] ERROR: ENA download was unavailable/failed and SRA fallback tools are missing (prefetch/fasterq-dump)." >&2
    FAILED_RUNS+=("$run")
    continue
  fi

  echo "[$run] Downloading .sra with prefetch..."
  if ! prefetch --output-directory "$SRA_CACHE_DIR" "$run"; then
    echo "[$run] ERROR: prefetch failed." >&2
    FAILED_RUNS+=("$run")
    continue
  fi

  # prefetch stores runs as <cache>/<SRR>/<SRR>.sra
  SRA_PATH="${SRA_CACHE_DIR}/${run}/${run}.sra"
  if [[ ! -f "$SRA_PATH" ]]; then
    echo "[$run] ERROR: Could not locate ${SRA_PATH}." >&2
    FAILED_RUNS+=("$run")
    continue
  fi

  echo "[$run] Converting to FASTQ..."
  # --split-files handles paired-end data into *_1 and *_2 files.
  if ! fasterq-dump --threads "$THREADS" --split-files --outdir "$FASTQ_DIR" "$SRA_PATH"; then
    echo "[$run] ERROR: fasterq-dump failed." >&2
    FAILED_RUNS+=("$run")
    continue
  fi
done

echo
echo "Ensuring FASTQ outputs remain compressed (.fastq.gz)..."
# Prefer pigz for faster parallel compression when available.
if command -v pigz >/dev/null 2>&1; then
  for fq in "${FASTQ_DIR}"/*.fastq; do
    # If no FASTQ files are present, break without error.
    [[ -e "$fq" ]] || break
    pigz -p "$THREADS" "$fq"
  done
else
  # Fallback to gzip when pigz is not installed.
  for fq in "${FASTQ_DIR}"/*.fastq; do
    [[ -e "$fq" ]] || break
    gzip "$fq"
  done
fi

echo "Cleaning temporary SRA cache..."
# Remove intermediate .sra cache to save disk space.
rm -rf "$SRA_CACHE_DIR"

echo
# Return non-zero if any run failed so workflows can detect partial failures.
if (( ${#FAILED_RUNS[@]} > 0 )); then
  echo "Completed with failures (${#FAILED_RUNS[@]} run(s)):" >&2
  printf '  - %s\n' "${FAILED_RUNS[@]}" >&2
  echo "Project directory: $PROJECT_DIR" >&2
  echo "FASTQ root directory: $FASTQ_ROOT_DIR" >&2
  echo "FASTQ raw directory: $FASTQ_DIR" >&2
  exit 1
fi

echo "Download complete."
echo "Project directory: $PROJECT_DIR"
echo "FASTQ root directory: $FASTQ_ROOT_DIR"
echo "FASTQ raw directory: $FASTQ_DIR"
echo "Log file: $LOG_FILE"
