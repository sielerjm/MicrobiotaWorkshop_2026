# 2-Hour Microbiome Analysis Workshop: From Raw Reads to Biological Insight

**Instructor:** Michael Sieler
**Last updated:** 2026-02-20

---

## Workshop Overview

This workshop provides a hands-on, end-to-end introduction to 16S rRNA amplicon
microbiome analysis. Participants will walk through the complete analytical
pipeline—from raw Illumina paired-end FastQ files through the DADA2 denoising
pipeline to a phyloseq object, then perform standard community-level analyses
including alpha diversity, beta diversity with PERMANOVA, and differential
abundance testing.

**Target audience:** Mixed—some participants comfortable with R, others newer to
both R and microbiome analysis.

**Format:** Instructor-led live-coding tutorial in R Markdown. Long-running
computation steps use pre-computed `.rds` objects so participants can follow
along in real time. Narrative explanations precede each code block, covering the
biological rationale, computational logic, and analytical interpretation.

**Dataset:** Honeybee 16S V4 amplicon data (515F–806R primers, Illumina MiSeq).
[TODO: Describe the dataset — experimental design, sample groups, replicates,
and biological context.]

---

## Timing Summary

| Section | Topic | Minutes | Mode |
|---------|-------|---------|------|
| 1.1 | Introduction & pipeline overview | 5 | Lecture/discussion |
| 1.2 | Inspect raw data & quality profiles | 10 | Live code |
| 1.3 | Primer detection & quality filtering | 10 | Live code |
| 1.4 | Error learning & denoising | 10 | Pre-computed + explain |
| 1.5 | ASV table, chimeras, taxonomy | 10 | Mixed live/pre-computed |
| 1.6 | Filtering & building phyloseq | 10 | Live code |
| — | **Break** | **5** | — |
| 2.1 | Data exploration & composition plots | 8 | Live code |
| 2.2 | Alpha diversity | 15 | Live code |
| 2.3 | Beta diversity & PERMANOVA | 20 | Live code |
| 2.4 | Differential abundance | 15 | Live code |
| 2.5 | Wrap-up | 2 | Discussion |
| | **Total** | **120** | |

---

## Part 1: From Raw FastQ to ASV Table (~55 min)

### 1.1 Introduction: What Are We Doing and Why? (5 min)

**Key teaching points:**

- **16S rRNA gene** is a universal phylogenetic marker in Bacteria and Archaea.
  Conserved regions enable universal primer binding; hypervariable regions (V1–V9)
  provide taxonomic discrimination.
- **Pipeline at a glance:**
  Raw paired-end reads → Quality filtering & primer trimming → Error-model-based
  denoising → Paired-read merging → Chimera removal → ASV table → Taxonomy
  assignment → Phyloseq object.
- **ASV vs. OTU:**
  - OTUs cluster sequences at ~97% identity—a pragmatic but arbitrary threshold
    that masks real biological variation.
  - ASVs (Amplicon Sequence Variants) use statistical error models to resolve
    exact biological sequences. They are reproducible, comparable across studies,
    and provide strain-level resolution.
  - DADA2 is the most widely used ASV inference tool (Callahan et al., 2016).

### 1.2 Inspect Raw Data and Quality (10 min) — *live code*

**Key teaching points:**

- **FastQ format:** Each read has a header, sequence, and per-base quality score
  encoded as ASCII characters.
- **Paired-end reads:** R1 (forward) and R2 (reverse) are sequenced from opposite
  ends of the same DNA fragment. They must be merged later.
- **Quality score interpretation:**
  - Q30 = 1 error per 1,000 bases (99.9% accuracy)
  - Q20 = 1 error per 100 bases (99% accuracy)
  - Quality degrades toward the 3' end because Illumina's sequencing-by-synthesis
    chemistry accumulates phasing errors over cycles.
  - R2 reads are generally lower quality than R1.
- **What participants do:** Run `plotQualityProfile()` on aggregated forward and
  reverse reads. Discuss where quality drops below acceptable thresholds and how
  that informs `truncLen` choices.

**Code functions:** `list.files()`, `plotQualityProfile()`

### 1.3 Primer Detection and Trimming/Filtering (10 min) — *live code*

**Key teaching points:**

- **Why check for primers:** Depending on the library prep protocol (Illumina
  two-step vs. EMP one-step), primer sequences may or may not be present at the
  5' end of reads. If present, they must be removed because:
  - Primers are synthetic, not biological signal.
  - Degenerate bases in primers create artificial sequence variants.
  - Primer-target mismatches inflate apparent diversity.
- **`filterAndTrim()` parameters explained:**
  - `trimLeft = c(19, 20)` — Remove forward (19 bp) and reverse (20 bp) primer
    sequences from the 5' end.
  - `truncLen` — Trim reads at the 3' end to remove low-quality tails. Set based
    on quality profile inspection. Must ensure sufficient overlap for merging
    (~20 bp minimum). For V4 with 515F–806R: ~252 bp insert, so F+R after
    trimming must sum to >252.
  - `maxN = 0` — Discard any read containing an ambiguous base (N). DADA2
    requires this.
  - `maxEE = c(2, 2)` — Maximum expected errors per read (sum of error
    probabilities across all bases). A read with maxEE=2 has, on average, no
    more than 2 incorrect bases. Stringent but retains most data.
  - `rm.phix = TRUE` — Remove PhiX spike-in sequences. PhiX is a short viral
    genome added to increase sequence complexity on Illumina runs (important for
    low-diversity amplicon libraries).
- **Verification:** After filtering, re-check for primer sequences to confirm
  removal. Inspect the `filterAndTrim` output table to ensure a reasonable
  fraction of reads pass (typically >80%).

**Code functions:** `stringr::str_locate()`, `filterAndTrim()`, `getSequences()`

### 1.4 Error Rate Learning and Denoising (10 min) — *pre-computed + explain*

**Key teaching points:**

- **Error model (the heart of DADA2):**
  - `learnErrors()` estimates the error rate for every possible nucleotide
    transition (e.g., A→C, A→G, ...) as a function of quality score.
  - The model is learned empirically from the data itself (~100M+ bases).
  - `plotErrors()` visualizes the fit: black line = estimated rates, red line =
    expected rates under nominal Q-scores, points = observed rates. The black line
    should approximately follow the points and decrease with increasing quality.
  - This is what distinguishes ASV methods from OTU clustering: instead of
    arbitrarily collapsing similar sequences, DADA2 uses a principled statistical
    model to determine which sequence differences are real biology vs. sequencing
    error.
- **Denoising and merging:**
  - `derepFastq()` collapses identical sequences for computational efficiency
    (keeping track of abundance).
  - `dada()` applies the error model to identify true biological sequences. Each
    unique input sequence is evaluated against the most abundant sequences to
    determine if it is likely an error variant.
  - `mergePairs()` joins denoised forward and reverse reads based on their
    overlap region. Default minimum overlap = 12 bp. Mismatches in the overlap
    region cause reads to be discarded.
  - Processing is done per-sample in a loop to manage memory.
- **Why pre-computed:** `learnErrors()` takes ~30 min and `dada()` takes ~10 min
  on this dataset. Participants load saved `.rds` files but see the actual code.

**Code functions:** `learnErrors()`, `plotErrors()`, `derepFastq()`, `dada()`,
`mergePairs()`

### 1.5 ASV Table, Chimera Removal, and Taxonomy (10 min) — *mixed*

**Key teaching points:**

- **Sequence table:** `makeSequenceTable()` creates the samples-by-ASV count
  matrix. Column names are the actual DNA sequences.
- **Length filtering:** Inspect the distribution of merged sequence lengths.
  For V4 (515F–806R), expect ~251–255 bp. Sequences far outside this range
  likely represent non-specific amplification or merging artifacts. Filter to
  the expected range.
- **Chimera removal:**
  - Chimeras are hybrid sequences formed during PCR when an incomplete extension
    product re-anneals to a different template in a subsequent cycle.
  - Can constitute ~30% of unique sequences in raw data.
  - `removeBimeraDenovo()` with `method="consensus"` identifies chimeras by
    checking if each ASV can be reconstructed as a combination of two more
    abundant parent sequences.
  - Chimeras typically represent a small fraction of total reads (because each
    chimeric sequence is low-abundance), so read loss is modest even if many
    unique sequences are flagged.
- **Tracking reads:** Build a table tracing read counts through every pipeline
  step (input → filtered → denoised → merged → length-filtered → non-chimeric).
  Expect >60% of input reads to survive. Large drops at specific steps indicate
  problems (e.g., massive loss at merging = insufficient overlap or over-trimming).
- **Taxonomy assignment:**
  - `assignTaxonomy()` uses a naive Bayesian classifier (Wang et al., 2007) to
    assign taxonomy against a reference database (SILVA v138).
  - `addSpecies()` attempts exact matching for species-level assignment.
  - Some ASVs will lack genus- or species-level classification—this is normal.
    The 16S V4 region cannot resolve all taxa to species, and some sequences
    have no close reference in the database.

**Code functions:** `makeSequenceTable()`, `removeBimeraDenovo()`,
`assignTaxonomy()`, `addSpecies()`

### 1.6 Filtering and Building Phyloseq (10 min) — *live code*

**Key teaching points:**

- **Remove non-target sequences:**
  - 16S primers can amplify host mitochondrial and chloroplast 16S rRNA.
  - Filter out: Kingdom != Bacteria/Archaea, Order == Chloroplast,
    Family == Mitochondria.
  - Use `%in%` operator instead of `!=` to correctly handle NA values in
    taxonomy assignments.
- **Prevalence filtering:**
  - Remove ASVs found in fewer than N samples (e.g., N = minimum group size).
  - Rationale: Very rare ASVs are often sequencing artifacts, contaminants, or
    transient environmental DNA. They add noise and sparsity without informative
    biological signal.
  - This also reduces the multiple-testing burden in downstream differential
    abundance analysis.
- **Build the phyloseq object:**
  - Combine ASV table + taxonomy table + sample metadata (+ optional
    phylogenetic tree) into a single phyloseq object.
  - This is the standard data structure for microbiome analysis in R—most
    downstream tools accept phyloseq objects directly.
  - Quick tour: `otu_table()`, `tax_table()`, `sample_data()`, `taxa_names()`,
    `sample_names()`, `ntaxa()`, `nsamples()`.

**Code functions:** `phyloseq()`, `subset_taxa()`, `prune_taxa()`,
`prune_samples()`, `sample_data()`

---

### Break (5 min)

---

## Part 2: From Phyloseq to Biological Insights (~60 min)

Participants load a pre-built, normalized phyloseq object as the starting point
for Part 2. This bridges the gap for anyone who encountered issues in Part 1 and
ensures everyone starts from the same object.

### 2.1 Data Exploration: Composition Overview (8 min) — *live code*

**Key teaching points:**

- **Sequencing depth check:** Use `sample_sums()` to examine read counts per
  sample. Highly uneven depth can bias diversity estimates.
- **Normalization (brief conceptual overview):**
  - *Rarefaction:* Subsample all samples to the same depth. Simple and widely
    used; discards data. Appropriate for diversity analyses.
  - *Relative abundance:* Convert counts to proportions. Compositional—changes
    in one taxon affect all others.
  - *CSS / DESeq2 variance stabilization:* Statistical normalization methods
    that account for library size differences without discarding data. Better for
    differential abundance testing.
  - For this tutorial: we use rarefied data for diversity and composition, and
    raw counts (with DESeq2 handling normalization internally) for differential
    abundance.
- **Composition barplots:**
  - Agglomerate ASVs to Family level with `tax_glom()`.
  - Convert to relative abundance with `transform_sample_counts()`.
  - Stacked barplots with `plot_bar()` faceted by treatment group.
  - Biological interpretation: Identify dominant taxa; note obvious shifts
    between groups; recognize these are descriptive (not statistical).

**Code functions:** `sample_sums()`, `tax_glom()`, `transform_sample_counts()`,
`plot_bar()`, `rarefy_even_depth()`

### 2.2 Alpha Diversity: Within-Sample Diversity (15 min) — *live code*

**Key teaching points:**

- **What is alpha diversity?** Measures the diversity *within* a single sample.
  Different metrics capture different facets of community structure.
- **Metrics explained biologically:**
  - **Observed richness** — "How many different taxa are present?"
    Simply counts distinct ASVs. Sensitive to sequencing depth (deeper sequencing
    detects more rare taxa) and to rare/singleton ASVs. Best used on rarefied
    data.
  - **Shannon index (H')** — "How even is the community?"
    Accounts for both richness and evenness. H' = -sum(p_i * ln(p_i)), where p_i
    is the proportion of species i. A community of 5 equally abundant species has
    higher Shannon than one dominated by a single species. Typical range for
    microbial communities: 0–5.
  - **Inverse Simpson (1/D)** — "How dominated is the community?"
    D = sum(p_i^2), giving probability that two randomly chosen individuals
    belong to the same species. Inverse Simpson = 1/D. Higher values = more
    diverse. Less sensitive to rare taxa than Shannon; emphasizes dominant
    species.
- **Visualization:** Boxplots of each metric by treatment group, with individual
  data points overlaid (jitter). Use ggplot2 with `facet_wrap()` for side-by-side
  comparison of metrics.
- **Statistical testing:**
  - For two groups: Wilcoxon rank-sum test (non-parametric; no assumption of
    normality—appropriate for diversity indices which are often non-normal).
  - For 3+ groups: Kruskal-Wallis test, followed by pairwise Wilcoxon with
    p-value adjustment if significant.
- **Interpretation guidance:**
  - Reduced alpha diversity is commonly observed in disease states ("dysbiosis")
    but is not universal.
  - Shannon and Simpson can give different signals: a community could have high
    richness but low evenness (high Observed, lower Shannon).
  - Always report which metric, the test used, and the effect size—not just
    the p-value.

**Code functions:** `estimate_richness()`, `ggplot()`, `wilcox.test()`,
`kruskal.test()`

### 2.3 Beta Diversity: Between-Sample Diversity (20 min) — *live code*

**2.3a Distance Metrics (5 min conceptual)**

- **Bray-Curtis dissimilarity** — Abundance-weighted. Ranges 0 (identical) to 1
  (completely different). The workhorse metric for most amplicon studies.
  Considers *how many* of each taxon are shared between samples.
- **Jaccard distance** — Presence/absence only. Ignores abundances entirely.
  Useful when you want to ask whether communities share the same *membership*
  regardless of relative proportions.
- **UniFrac (weighted and unweighted)** — Incorporates phylogenetic tree
  information (fraction of shared branch length between samples).
  - *Unweighted UniFrac:* Qualitative—considers only presence/absence on the
    tree. Sensitive to rare lineages.
  - *Weighted UniFrac:* Quantitative—weights branches by abundance. Emphasizes
    shifts in dominant taxa.
  - Requires a phylogenetic tree in the phyloseq object.
- **Choosing a metric:**
  - Bray-Curtis: General-purpose default for most studies.
  - Jaccard: When rare taxa should count equally (e.g., pathogen detection).
  - UniFrac: When evolutionary relationships matter (e.g., testing if communities
    shift toward phylogenetically distinct taxa under treatment).

**2.3b PCoA Ordination (5 min) — *live code***

- **What is PCoA?** Principal Coordinates Analysis projects a high-dimensional
  distance matrix into a low-dimensional space (typically 2D) that preserves
  pairwise distances as well as possible.
- Each axis captures a fraction of total variation (reported as %). Axis 1
  captures the most, Axis 2 the second most, etc.
- Points = samples. Points closer together = more similar communities.
- Color by treatment group to visually assess whether groups cluster separately.
- **Caution:** Ordination is visualization, not a statistical test. Apparent
  separation must be confirmed with PERMANOVA.

**Code functions:** `ordinate()`, `plot_ordination()`

**2.3c Statistical Testing: PERMANOVA + betadisper (10 min) — *live code***

- **PERMANOVA (`vegan::adonis2()`):**
  - Tests whether the centroids (mean positions) of groups differ significantly
    in multivariate space.
  - Permutation-based: shuffles group labels thousands of times to build a null
    distribution, then asks how extreme the observed between-group vs.
    within-group variation ratio (pseudo-F) is.
  - Reports: R-squared (proportion of variation explained by grouping), pseudo-F
    statistic, and p-value.
  - Formula: `distance_matrix ~ grouping_variable`, with `permutations = 9999`.
- **Betadisper (`vegan::betadisper()` + `permutest()`):**
  - PERMANOVA is sensitive to *both* location (centroid differences) and
    dispersion (spread differences). A significant PERMANOVA could arise because
    groups have different *spread* rather than different *centers*.
  - `betadisper()` tests whether group dispersions are homogeneous.
  - If betadisper is non-significant (p > 0.05): dispersions are similar,
    so the PERMANOVA result likely reflects true centroid differences.
  - If betadisper is significant: interpret PERMANOVA with caution—the groups
    may simply be more variable, not compositionally different.
- **Interpretation example:**
  "Bee type significantly explains X% of the variation in community composition
  (PERMANOVA: R^2 = X, pseudo-F = Y, p < 0.001). Group dispersions are
  homogeneous (betadisper p = 0.45), confirming that different bee castes harbor
  compositionally distinct gut microbiomes."

**Code functions:** `phyloseq::distance()`, `vegan::adonis2()`,
`vegan::betadisper()`, `permutest()`

### 2.4 Differential Abundance: Which Taxa Differ? (15 min) — *live code*

**Key teaching points:**

- **Why DA analysis?** Alpha and beta diversity tell us communities differ
  *overall*. DA analysis identifies the *specific taxa* driving those
  differences—the biologically actionable output.
- **Prevalence filtering before DA:** Remove ASVs present in <10% of samples.
  Reduces noise and the multiple-testing burden.
- **Approach: Wilcoxon rank-sum test per taxon:**
  - For each taxon, test: does abundance differ between groups?
  - Non-parametric—appropriate for the skewed, zero-inflated distributions
    typical of microbiome count data.
  - Use `lapply()` to loop across all taxa.
  - **Multiple testing correction:** When testing hundreds of taxa simultaneously,
    some will appear significant by chance. Apply FDR (Benjamini-Hochberg)
    correction. More liberal than Bonferroni but controls the expected proportion
    of false discoveries.
- **Visualization:**
  - *Beeswarm/strip plots:* Show abundance distributions for significant taxa,
    split by group. Include error bars (95% CI of the mean).
  - *Heatmap:* Log-transformed mean abundances for significant taxa across
    groups. Compact summary when many taxa are significant.
- **Brief mention of alternative methods:**
  - DESeq2: Negative binomial model. Handles normalization internally. Better
    for complex designs with multiple covariates.
  - ANCOM-BC: Addresses compositionality explicitly.
  - Different methods can yield different significant taxa (Nearing et al., 2022,
    Nature Communications). Recommend running 2+ methods and focusing on
    concordant results.
- **Critical interpretation:**
  - Effect size matters as much as p-value.
  - Consider biological plausibility—does this taxon make sense in the context
    of the experimental system?
  - DA results are hypothesis-generating, not definitive—validate with
    independent methods or datasets when possible.

**Code functions:** `wilcox.test()`, `p.adjust()`, `data.table::melt()`,
`ggplot()` + `ggbeeswarm::geom_beeswarm()`, `geom_tile()` for heatmap

### 2.5 Wrap-Up (2 min)

- **Recap the full workflow:**
  Raw reads → Quality filtering → DADA2 denoising → ASV table → Taxonomy →
  Phyloseq → Composition → Alpha diversity → Beta diversity + PERMANOVA →
  Differential abundance
- **Each analysis answers a different question:**

  | Analysis | Question |
  |----------|----------|
  | Composition plots | What taxa are present and at what proportions? |
  | Alpha diversity | How diverse is each sample, and do groups differ in diversity? |
  | Beta diversity + PERMANOVA | Do overall community structures differ between groups? |
  | Differential abundance | Which specific taxa drive the differences? |

- **Resources for further learning:**
  - DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
  - phyloseq documentation: https://joey711.github.io/phyloseq/
  - microViz package: https://david-barnett.github.io/microViz/
  - Nearing et al. (2022) DA method comparison: https://doi.org/10.1038/s41467-022-28034-z
  - Schloss (2024) on rarefaction: relevant for normalization decisions

---

## Data Requirements

> **[TODO: Update file names and sources once workshop dataset is finalized.]**

| File | Description | Source |
|------|-------------|--------|
| Raw FastQ files (R1/R2) | Paired-end 16S V4 reads | [TODO: source] |
| `[TODO: metadata_file]` | Sample metadata | [TODO: source] |
| `errF_sampled.rds`, `errR_sampled.rds` | Pre-computed error models | Generated from workshop data |
| `dadaFs_sampled.rds`, `dadaRs_sampled.rds` | Pre-computed denoised reads | Generated from workshop data |
| `mergers_sampled.rds` | Pre-computed merged pairs | Generated from workshop data |
| `seqtab_sampled.rds` | Pre-computed sequence table | Generated from workshop data |
| `seqtabnochim_sampled.rds` | Pre-computed chimera-free table | Generated from workshop data |
| `taxa_sampled.rds` | Pre-computed taxonomy | Generated from workshop data |
| `filterAndTrim_sampled.rds` | Pre-computed filter stats | Generated from workshop data |
| SILVA ref databases | `silva_nr99_v138.1_train_set.fa.gz`, species assignment | [SILVA/DADA2 Zenodo](https://zenodo.org/records/14169026) |
| `[TODO: phyloseq_object.rds]` | Pre-built phyloseq for Part 2 | Generated from workshop data |

---

## Deliverables

1. **This outline** (`Workshop_Outline.md`) — for planning and review
2. **Tutorial R Markdown** (`Microbiome_Workshop_Tutorial.Rmd`) — the complete
   runnable tutorial with narrative explanations
3. **Data inventory** (`Data_Inventory.md`) — documents all required files and
   their source paths
