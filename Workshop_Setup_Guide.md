# Workshop Environment Setup Guide

**Workshop:** 2-Hour Microbiome Analysis Tutorial (16S rRNA gene analysis)
**Date:** April 14, 2026
**Time:** 2-hour slot (exact time TBD with course coordinator)
**Location:** Room 119, Pavillon Ansermet, University of Geneva
**Course:** MSc course coordinated by Athanasia Tzika
**Organizer:** Joanito Liberti
**Tutorial lead:** Michael Sieler
**IT contact:** Stephan (via Joanito's email)

**Created by:** Michael Sieler
**Last updated:** 2026-02-20

---

## 1. Goal

Students log in to pre-configured RStudio sessions with all R packages
installed and all data files in place. Zero setup time during the session.

---

## 2. Recommended Approach: RStudio Server or Posit Cloud

### Option A: University-hosted RStudio Server (preferred)

If the university or department has an existing RStudio Server instance (or can
provision one), this is the cleanest option. Each student gets a browser-based
RStudio session backed by a shared server.

**What we need from IT (Stephan):**

- An RStudio Server instance accessible from the room 119 computers (or
  students' own laptops via university network/VPN)
- Student accounts (one per student, or a shared guest login)
- R >= 4.3.0 with Bioconductor 3.18+ installed
- All R packages listed in Section 4 pre-installed system-wide
- The workshop project directory (Section 5) copied into each user's home
  directory, **or** a shared read-only location that users can copy from
- Sufficient RAM: at least 4 GB per concurrent user (DADA2's `filterAndTrim()`
  and `removeBimeraDenovo()` can use ~2 GB on the tutorial dataset)

**Advantages:** No dependency on external services, university controls the
environment, fast local network access.

### Option B: Posit Cloud (cloud.posit.co)

Posit Cloud (formerly RStudio Cloud) provides browser-based RStudio workspaces
with no local installation required.

**Setup steps:**

1. Create a free Posit Cloud account at https://posit.cloud
2. Create a new **Space** for the workshop (e.g., "Microbiome Workshop 2026")
3. Create a **Project** within the Space:
   - Upload all data files (see Section 5)
   - Install all R packages (see Section 4)
   - Upload `Microbiome_Workshop_Tutorial.Rmd`
4. Set the Project as an **Assignment** — each student gets their own copy when
   they open it
5. Share the Space with students via a join link (no individual invitations
   needed)

**Free tier limitations:**

- 25 project hours/month per user
- 1 GB RAM, 1 CPU per project
- This is tight for DADA2 steps, but the tutorial uses pre-computed `.rds`
  files for the computationally heavy steps, so it should be workable

**Paid option:** Posit Cloud Premium ($99/month) gives 4 GB RAM per project and
higher CPU, which is more comfortable. A one-month subscription for the
workshop may be worthwhile.

**Advantages:** Zero IT setup, works on any computer with a browser, students
can continue working after the session.

### Option C: Local RStudio on room 119 computers

If the computers in room 119 already have R and RStudio installed, packages and
data can be pre-installed on each machine. This requires either:

- Admin access to install R packages into the system library, or
- A shared network drive with a pre-built R library path that we set via
  `.Rprofile` or `.libPaths()`

**What we would need from IT:**

- Confirm R and RStudio are installed on room 119 computers (version numbers)
- Admin access (even temporarily) to install Bioconductor packages, or a
  writable shared library path
- The workshop project folder placed on each machine's desktop or a mapped
  network drive

**Disadvantages:** Requires touching each machine; version inconsistencies
across computers can cause problems.



---

## 3. Questions for Stephan

When coordinating with Stephan, the key questions are:

1. **Does the university have an RStudio Server** 
   - That students can access from room 119?
2. **If not, what is installed on the room 119 computers?** Specifically:
   - Is R installed? What version?
   - Is RStudio installed?
   - Can we install R packages in advance (do we have admin access)?
   - Is there a shared network drive accessible from all machines?
3. **How many students** are expected? (affects RAM/CPU planning)
4. **Is there reliable internet access** in room 119? (needed for Posit Cloud;
   not needed for local/server options)
5. **Can we schedule a test session** in room 119 before April 14 to verify
   everything works?

---

## 4. Required R Packages

The tutorial itself uses a core set of packages. We also include extended
pipeline packages so that students can explore further or so that the
environment is ready if we extend examples during Q&A.

### Bioconductor packages

| Package | Purpose | Used in |
|---------|---------|---------|
| `BiocManager` | Bioconductor package installer | Setup |
| `dada2` | DADA2 amplicon processing pipeline | Tutorial Part 1 |
| `phyloseq` | Microbiome data structure and analysis | Tutorial Parts 1 & 2 |
| `decontam` | Contaminant identification | Tutorial Part 1 (example) |
| `microViz` | Extended phyloseq utilities (tax_fix, comp_barplot, validation) | Pipeline / extra |
| `DESeq2` | Differential abundance via negative binomial models | Pipeline / extra |
| `ShortRead` | FastQ file reading and counting | Pipeline / extra |
| `Biostrings` | DNA sequence manipulation (reverse complement, etc.) | Pipeline / extra |
| `DECIPHER` | Multiple sequence alignment for tree building | Pipeline / extra |
| `genefilter` | Filter taxa by prevalence/abundance thresholds | Pipeline / extra |

### CRAN packages

| Package | Purpose | Used in |
|---------|---------|---------|
| `tidyverse` | Data manipulation & visualization (includes ggplot2, dplyr, tidyr, stringr, purrr, readr, forcats, tibble) | Tutorial + Pipeline |
| `vegan` | Community ecology: Bray-Curtis, PERMANOVA, betadisper, Mantel | Tutorial Part 2 + Pipeline |
| `data.table` | Efficient data manipulation | Tutorial Part 2 |
| `knitr` | R Markdown rendering | Tutorial |
| `rmarkdown` | R Markdown rendering | Tutorial |
| `ade4` | Multivariate analysis, Procrustes rotation | Pipeline / extra |
| `phangorn` | Phylogenetic tree construction (NJ, GTR+G+I) | Pipeline / extra |
| `lmerTest` | Linear mixed-effects models with Satterthwaite p-values | Pipeline / extra |
| `car` | Type II/III ANOVA tables for mixed models | Pipeline / extra |
| `RColorBrewer` | Color palette sets for plots | Pipeline / extra |
| `cowplot` | Plot composition helpers (plot_grid) | Pipeline / extra |
| `ggpubr` | Publication-ready ggplot extras | Pipeline / extra |
| `ggsignif` | Significance brackets on ggplots | Pipeline / extra |
| `patchwork` | Combine multiple ggplots into composite figures | Pipeline / extra |
| `scales` | Axis/legend scale helpers for ggplot2 | Pipeline / extra |
| `gt` | Publication-quality tables | Pipeline / extra |
| `furrr` | Parallel map using future back-ends | Pipeline / extra |
| `future` | Parallelization engine | Pipeline / extra |
| `jsonlite` | JSON serialization | Pipeline / extra |

### Installation script

Save as `install_packages.R` and run on the server/each machine before the
workshop:

```r
# install_packages.R
# Pre-install all R packages for the Microbiome Workshop
# Includes both tutorial-essential packages and extended pipeline packages
# Run this script ONCE before the workshop session
#
# Estimated time: 20-40 minutes (Bioconductor packages compile from source)
# Requires: internet access, R >= 4.3.0

cat("=== Microbiome Workshop: Installing R packages ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# --- Bioconductor packages ---
cat("--- Installing Bioconductor packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_packages <- c(
  "dada2", "phyloseq", "decontam",    # Tutorial core
  "DESeq2", "ShortRead", "Biostrings", # Extended pipeline
  "DECIPHER", "genefilter"             # Extended pipeline
)
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

# microViz from r-universe (not on Bioconductor/CRAN)
cat("\n--- Installing microViz from r-universe ---\n")
install.packages(
  "microViz",
  repos = c(
    davidbarnett = "https://david-barnett.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)

# --- CRAN packages ---
cat("\n--- Installing CRAN packages ---\n")
cran_packages <- c(
  # Tutorial core
  "tidyverse", "vegan", "data.table", "knitr", "rmarkdown",
  # Extended pipeline: statistics & modeling
  "ade4", "phangorn", "lmerTest", "car",
  # Extended pipeline: visualization
  "RColorBrewer", "cowplot", "ggpubr", "ggsignif",
  "patchwork", "scales", "gt",
  # Extended pipeline: utilities
  "furrr", "future", "jsonlite"
)
install.packages(cran_packages, repos = "https://cloud.r-project.org")

# --- Verify all packages load ---
cat("\n=== Verification ===\n")
all_packages <- c(
  # Bioconductor
  "dada2", "phyloseq", "decontam", "microViz",
  "DESeq2", "ShortRead", "Biostrings", "DECIPHER", "genefilter",
  # CRAN - tidyverse components
  "ggplot2", "dplyr", "tidyr", "stringr", "purrr",
  "readr", "forcats", "tibble",
  # CRAN - ecology & stats
  "vegan", "ade4", "lmerTest", "car",
  # CRAN - visualization
  "RColorBrewer", "cowplot", "ggpubr", "ggsignif",
  "patchwork", "scales", "gt",
  # CRAN - utilities
  "data.table", "phangorn", "furrr", "future",
  "jsonlite", "knitr", "rmarkdown"
)

results <- sapply(all_packages, function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  status <- ifelse(ok, "OK", "FAILED")
  cat(sprintf("  %-15s %s\n", pkg, status))
  ok
})

n_ok     <- sum(results)
n_fail   <- sum(!results)
cat(sprintf("\n--- Results: %d / %d packages OK ---\n", n_ok, length(results)))

if (n_fail > 0) {
  cat("FAILED packages:\n")
  cat(paste("  ", names(results[!results])), sep = "\n")
  cat("\nPlease install these manually before the workshop.\n")
} else {
  cat("All packages installed successfully.\n")
}

cat("\nFinished:", format(Sys.time()), "\n")
```

---

## 5. Required Data Files

The complete file listing is in `Data_Inventory.md`. Summary:

| Directory | Contents | Size (approx.) |
|-----------|----------|----------------|
| `raw/` | 20 paired-end FastQ files (10 samples) | ~120 MB |
| `rds/` | Pre-computed DADA2 intermediate results (12 `.rds` files) | ~50 MB |
| `ref/` | SILVA v138.1 reference databases (optional) | ~120 MB |
| [TODO: Part 2 data dir] | Pre-built phyloseq object for Part 2 | ~5 MB |
| `filtered/` | Empty directory (created by DADA2 during the session) | 0 |
| Root | `[TODO: metadata_file]`, `Microbiome_Workshop_Tutorial.Rmd` | < 1 MB |

**Total without SILVA:** ~175 MB
**Total with SILVA:** ~295 MB

### Data setup script

[TODO: Update this script once the workshop dataset is finalized.]

```bash
#!/usr/bin/env bash
# setup_workshop_data.sh
# Assembles the workshop data directory

SOURCE="[TODO: path to source data]"
TUTORIAL="/path/to/workshop_directory"

mkdir -p "$TUTORIAL/raw" "$TUTORIAL/rds" "$TUTORIAL/ref" \
         "$TUTORIAL/filtered"

cp "$SOURCE/raw/"*.fastq.gz "$TUTORIAL/raw/"
cp "$SOURCE/[TODO: metadata_file]" "$TUTORIAL/"
cp "$SOURCE/rds/"*.rds "$TUTORIAL/rds/"
cp "$SOURCE/ref/"*.gz "$TUTORIAL/ref/"
# cp "$SOURCE/[TODO: phyloseq_object.rds]" "$TUTORIAL/"

echo "Workshop data assembled in: $TUTORIAL"
du -sh "$TUTORIAL"
```

### Per-user data placement

Depending on the platform, data should be placed:

| Platform | Location |
|----------|----------|
| RStudio Server | `/home/<user>/workshop/` (copy per user, or symlink to shared read-only) |
| Posit Cloud | Uploaded to the project (automatically copied per student via Assignment) |
| Local machines | `C:\Users\<user>\Desktop\workshop\` or similar |
| Docker | Baked into the image at `/home/rstudio/workshop/` |

---

## 6. Pre-Workshop Checklist

### 2+ weeks before (now)

- [ ] Contact Stephan to determine available infrastructure (RStudio Server,
      local installs, or need for Posit Cloud)
- [ ] Confirm number of students expected
- [ ] Decide on platform (Option A, B, C, or D above)

### 1 week before

- [ ] R packages installed and verified on the target platform
- [ ] Data files copied to the correct location(s)
- [ ] Test the full tutorial end-to-end on the target platform
  - Open `Microbiome_Workshop_Tutorial.Rmd` in RStudio
  - Run all chunks in Part 1 (should complete in ~5 min with pre-computed data)
  - Run all chunks in Part 2 (should complete in ~3 min)
  - Verify all plots render correctly
- [ ] Confirm login credentials / access URLs for students
- [ ] Prepare a short "getting started" slide or handout with:
  - URL to access RStudio (if server/cloud)
  - Login credentials
  - Instructions to open the `.Rmd` file

### 1-2 days before

- [ ] Do a final smoke test from a room 119 computer (or equivalent network)
- [ ] Ensure `filtered/` directory is empty (or recreated) so `filterAndTrim()`
      can write to it fresh during the session
- [ ] Clear any cached `.Rmd` outputs so students start clean
- [ ] Prepare backup: have a pre-rendered HTML version of the tutorial available
      in case something goes wrong with the live environment

### Day of (30 min before session)

- [ ] Log into the server/cloud and verify it is responding
- [ ] Open the Rmd in one session and run the first chunk to confirm packages
      load
- [ ] Have the rendered HTML tutorial open in a browser tab as a fallback
- [ ] Ensure room 119 projector/screen is working

---

## 7. Contingency Plans

| Scenario | Fallback |
|----------|----------|
| Server is down | Switch to pre-rendered HTML walkthrough; students follow along and discuss |
| A package fails to load | Have a pre-rendered HTML available; fix the package during a break |
| Internet is down (Posit Cloud) | If using cloud: switch to local fallback USB drive with portable RStudio + data |
| `filterAndTrim()` is slow | Already handled: tutorial uses pre-computed `.rds` files for all slow steps |
| Students can't log in | Have 2-3 shared "guest" accounts ready; pair students |

---

## 8. Draft Message to Stephan

Below is a draft follow-up message that Michael can send (or that Joanito can
forward) with specific technical requirements:

---

> Dear Stephan,
>
> Thank you for helping us set up the computing environment for the microbiome
> analysis practical on April 14 in room 119. Below are the technical details
> for what we need.
>
> **Summary:** Students will run an R Markdown tutorial in RStudio that covers
> 16S rRNA gene microbiome analysis. The tutorial requires several R packages
> (including Bioconductor packages) and ~175 MB of data files per student.
>
> **Our ideal setup:** A browser-based RStudio environment (RStudio Server or
> similar) where students log in and everything is ready to go. If that is not
> available, we can work with local R/RStudio installations on the room
> computers.
>
> **Specific questions:**
>
> 1. Does the university have an RStudio Server instance that MSc students can
>    access from room 119? If so, can we pre-install Bioconductor packages
>    (dada2, phyloseq, decontam) and CRAN packages (tidyverse, vegan,
>    data.table) on it?
>
> 2. If not, what software is currently installed on the room 119 computers?
>    (R version, RStudio, etc.)
>
> 3. Can we get temporary admin access to install packages, or would you be able
>    to run an installation script for us?
>
> 4. How many student workstations are available in room 119?
>
> 5. Is there reliable internet access in the room? (We have a cloud-based
>    backup plan if needed.)
>
> 6. Would it be possible to schedule a 30-minute test session in room 119
>    before April 14 to verify the setup?
>
> We have a complete installation script (`install_packages.R`) and data setup
> script ready to share whenever convenient.
>
> Thank you very much for your help!
>
> Best regards,
> Michael

---

## 9. What Students Need to Know Before the Session

Depending on the platform, send students a brief email ~3 days before:

**If using RStudio Server / Posit Cloud:**

> For the microbiome analysis practical on April 14, you will use a
> browser-based RStudio environment. No software installation is needed.
>
> - Open your web browser and go to: [URL]
> - Log in with: [credentials or "use your university login"]
> - Open the file `Microbiome_Workshop_Tutorial.Rmd`
> - We will work through the tutorial together during the session
>
> Please bring a laptop if you prefer to follow along on your own screen
> (optional -- room computers will also be available).

**If using local installations:**

> For the microbiome analysis practical on April 14, R and RStudio are
> pre-installed on the room 119 computers with all required packages.
>
> - Log in to the computer with your university credentials
> - Open RStudio from the desktop / applications menu
> - Navigate to the `workshop` folder and open `Microbiome_Workshop_Tutorial.Rmd`
>
> No prior R experience is required, though familiarity with R basics is
> helpful.
