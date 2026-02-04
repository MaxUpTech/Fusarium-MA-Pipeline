# Fusarium MA Pipeline - Complete Workflow Guide

## Overview

This pipeline detects de novo mutations in Mutation Accumulation (MA) lines by comparing
evolved lines to their ancestor. It integrates statistical models from the muver pipeline
to improve accuracy and reduce false positives.

---

## Pipeline Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           INPUT FILES                                            │
│  - FASTQ files (paired-end reads)                                               │
│  - Reference genome (FASTA)                                                      │
│  - Sample mapping file (which sample is ancestor, which are MA lines)           │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│  STEP 1: PREPROCESSING (fastp)                                                   │
│  - Quality trimming                                                              │
│  - Adapter removal                                                               │
│  - Quality reports                                                               │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│  STEP 2: ALIGNMENT (BWA-MEM2)                                                    │
│  - Align reads to reference                                                      │
│  - Sort BAM files                                                                │
│  - Mark duplicates (GATK)                                                        │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│  STEP 3: VARIANT CALLING (GATK)                                                  │
│  - HaplotypeCaller per sample (GVCF mode)                                        │
│  - Joint genotyping across all samples                                           │
│  - Hard filtering (QD, FS, MQ)                                                   │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│  STEP 4: MUVER MODEL PREPARATION                                                 │
│  Uses ANCESTOR BAM to characterize:                                              │
│  - Depth distribution (normal distribution fitting)                              │
│  - Strand bias distribution (log-normal fitting)                                 │
│  - Repeat regions (identify & fit INDEL rates)                                   │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│  STEP 5: ENHANCED MUTATION CALLING                                               │
│  For each variant in each MA line:                                               │
│  - Check if absent in ancestor (de novo)                                         │
│  - Apply muver depth filter                                                      │
│  - Apply muver strand bias filter                                                │
│  - Apply repeat INDEL correction                                                 │
│  - Calculate composite significance                                              │
│  - Detect subclonal variants                                                     │
│  - FWER correction                                                               │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│  STEP 6: STATISTICAL ANALYSIS                                                    │
│  - Mutation rates (per-sample, per-genotype)                                     │
│  - Confidence intervals (Bootstrap, Bayesian)                                    │
│  - Mutation spectrum analysis                                                    │
│  - Effect size calculations                                                      │
│  - Advanced statistical tests                                                    │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│  STEP 7: VISUALIZATION & REPORTING                                               │
│  - Individual sample plots                                                       │
│  - Group comparison plots                                                        │
│  - Combined analysis plots                                                       │
│  - Statistical diagnostic plots                                                  │
│  - Interactive HTML report                                                       │
└─────────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           OUTPUT FILES                                           │
│  - mutations.tsv / mutations.vcf (called mutations)                              │
│  - mutation_rates.tsv (rates per sample)                                         │
│  - Plots (PNG files in plots/ subdirectories)                                    │
│  - analysis_report.html (comprehensive report)                                   │
└─────────────────────────────────────────────────────────────────────────────────┘
```

---

## Detailed Step-by-Step Breakdown

---

### STEP 1: Preprocessing

**Module:** `core/preprocessing.py`

**Purpose:** Clean raw sequencing reads to improve alignment quality.

| Input | Output |
|-------|--------|
| Raw FASTQ files (R1, R2) | Trimmed FASTQ files |
| | QC reports (HTML, JSON) |

**What happens:**
1. fastp reads paired-end FASTQ files
2. Removes adapter sequences (auto-detected)
3. Trims low-quality bases (Q < 20)
4. Removes reads shorter than 50bp
5. Generates quality reports

**Example:**
```
Input:  sample1_R1.fastq.gz, sample1_R2.fastq.gz
Output: sample1_trimmed_R1.fastq.gz, sample1_trimmed_R2.fastq.gz
        sample1_fastp.html, sample1_fastp.json
```

---

### STEP 2: Alignment

**Module:** `core/alignment.py`

**Purpose:** Map reads to reference genome and remove PCR duplicates.

| Input | Output |
|-------|--------|
| Trimmed FASTQ files | Sorted, deduplicated BAM files |
| Reference genome (FASTA) | BAM index files (.bai) |
| | Duplicate metrics |

**What happens:**
1. BWA-MEM2 indexes reference (if needed)
2. Aligns reads to reference with read groups
3. samtools sorts alignments by coordinate
4. GATK MarkDuplicates flags PCR duplicates
5. Creates BAM index

**Example:**
```
Input:  sample1_trimmed_R1.fastq.gz, sample1_trimmed_R2.fastq.gz
        Fusarium_graminearum.fa
Output: sample1.markdup.bam, sample1.markdup.bai
        sample1.metrics.txt
```

---

### STEP 3: Variant Calling

**Module:** `core/variant_calling.py`

**Purpose:** Identify all genomic variants across samples.

| Input | Output |
|-------|--------|
| BAM files (all samples) | Joint-called VCF |
| Reference genome | Filtered VCF |

**What happens:**
1. GATK HaplotypeCaller runs on each sample → GVCF files
2. GenomicsDBImport combines all GVCFs
3. GenotypeGVCFs performs joint genotyping
4. VariantFiltration applies hard filters:
   - SNPs: QD < 2, FS > 60, MQ < 40 → FILTERED
   - INDELs: QD < 2, FS > 200 → FILTERED
5. Merges SNP and INDEL calls

**Example:**
```
Input:  ancestor.markdup.bam, ma_line1.markdup.bam, ma_line2.markdup.bam, ...
Output: filtered.vcf.gz (contains all samples, all variants)
```

**VCF format (simplified):**
```
#CHROM  POS     REF  ALT  QUAL   Ancestor  MA_Line1  MA_Line2
chr1    12345   A    G    1500   0/0       0/1       0/0
chr1    23456   C    T    2000   0/0       0/0       1/1
chr2    34567   AT   A    800    0/0       0/1       0/1
```

---

### STEP 4: Muver Model Preparation

**Module:** `muver/` (multiple files)

**Purpose:** Characterize the sequencing data to set appropriate filtering thresholds.

This step uses the **ANCESTOR BAM** to learn what "normal" looks like, then applies
these learned parameters when filtering mutations in MA lines.

#### 4A: Depth Distribution Analysis

**File:** `muver/depth_distribution.py`

| Input | Output |
|-------|--------|
| Ancestor BAM file | Normal distribution parameters (μ, σ) |
| Reference genome | Filtered regions BED file |

**What happens:**
1. Calculates read depth at every position using pileup
2. Fits a normal (Gaussian) distribution to depth values
3. Identifies regions with abnormal depth using CDF:
   - Too low depth: P(depth) < 0.0001
   - Too high depth: P(depth) > 0.9999
4. Uses sliding window to smooth and merge filtered regions

**Why this matters:**
- Regions with abnormal depth often indicate:
  - Repetitive sequences (high depth from mismapping)
  - Deletions or difficult regions (low depth)
  - Copy number variants
- Variants in these regions are unreliable

**Parameters learned:**
```
μ (mean depth) = 50x (example)
σ (std dev) = 15 (example)
Filtered regions: chr1:1000-2000, chr2:5000-5500, ...
```

#### 4B: Strand Bias Distribution Analysis

**File:** `muver/bias_distribution.py`

| Input | Output |
|-------|--------|
| Ancestor BAM file | Log-normal distribution parameters (μ, σ) |
| Reference genome | |

**What happens:**
1. For each position, counts reads on forward vs reverse strand
2. Calculates log-ratio: ln(forward_count / reverse_count)
3. Fits a Gaussian to the log-ratios (making it log-normal)
4. Real variants should have balanced strand representation
5. Artifacts often show extreme strand bias

**Why this matters:**
- True variants: reads from both strands support the variant
- Artifacts: often only one strand shows the "variant" (sequencing error, mapping error)

**Parameters learned:**
```
μ (mean log-ratio) ≈ 0 (balanced strands)
σ (std dev) ≈ 0.3 (typical variation)
```

#### 4C: Repeat Region Analysis

**File:** `muver/repeat_indels.py`

| Input | Output |
|-------|--------|
| Reference genome | Repeat regions BED file |
| Ancestor BAM | Logistic fit parameters per repeat unit |

**What happens:**
1. Scans reference for tandem repeats (e.g., ATATAT, GCGCGC)
2. Records: repeat unit (AT), tract length (6 = 3 repeats of AT)
3. Counts INDELs in repeat regions from BAM
4. Fits logistic function: INDEL_rate = f(tract_length)
5. Longer repeats → higher INDEL rate (polymerase slippage)

**Why this matters:**
- INDELs in repeats are often sequencing/alignment artifacts
- Longer repeat tracts have higher error rates
- This correction prevents false positive INDEL calls

**Parameters learned (per repeat unit):**
```
Repeat unit "AT":
  x0 = 8 (inflection point)
  L = 0.1 (maximum rate)
  k = 0.5 (steepness)
```

---

### STEP 5: Enhanced Mutation Calling

**Module:** `core/mutation_calling.py`

**Purpose:** Identify TRUE de novo mutations by comparing MA lines to ancestor,
using muver models to filter artifacts.

| Input | Output |
|-------|--------|
| Filtered VCF | mutations.tsv |
| Muver parameters | mutations.vcf |
| Sample information | Mutation objects |

**What happens for EACH variant in EACH MA line:**

```
┌─────────────────────────────────────────────────────────────────┐
│  For variant at chr1:12345 A→G in MA_Line_1:                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. CHECK: Is variant PRESENT in MA line?                       │
│     - Genotype = 0/1 or 1/1? → YES, continue                    │
│     - Genotype = 0/0? → SKIP (no variant)                       │
│                                                                 │
│  2. CHECK: Is variant ABSENT in ancestor?                       │
│     - Ancestor genotype = 0/0? → YES, de novo candidate         │
│     - Ancestor genotype ≠ 0/0? → SKIP (inherited, not de novo)  │
│                                                                 │
│  3. MUVER FILTER: Depth distribution                            │
│     - Is position in filtered region? → FAIL                    │
│     - Is depth within normal range? → PASS                      │
│                                                                 │
│  4. MUVER FILTER: Strand bias                                   │
│     - Calculate log(forward/reverse) for variant reads          │
│     - P-value from learned distribution                         │
│     - P < 0.01? → FAIL (significant bias)                       │
│                                                                 │
│  5. MUVER CORRECTION: Repeat INDELs (if INDEL)                  │
│     - Is variant in repeat region?                              │
│     - Calculate expected error rate from logistic fit           │
│     - Adjust significance threshold                             │
│                                                                 │
│  6. STATISTICAL TESTS:                                          │
│     - Binomial test (forward strand)                            │
│     - Binomial test (reverse strand)                            │
│     - Chi-square test (allele distribution)                     │
│                                                                 │
│  7. COMPOSITE SIGNIFICANCE:                                     │
│     - Combine p-values: sqrt(log(p1)² + log(p2)² + log(p3)²)    │
│     - Higher score = more confident                             │
│                                                                 │
│  8. SUBCLONAL DETECTION:                                        │
│     - Test frequencies: 0.5, 0.25, 0.125                        │
│     - Maximum likelihood genotype calling                       │
│                                                                 │
│  9. FWER CORRECTION:                                            │
│     - Adjust for multiple testing                               │
│     - Threshold = 0.01 / genome_size                            │
│                                                                 │
│  → RESULT: PASS/FAIL + confidence metrics                       │
└─────────────────────────────────────────────────────────────────┘
```

**Output fields in mutations.tsv:**
```
sample          chromosome  position  ref  alt  type  quality  depth  ...
MA_Line_1       chr1        12345     A    G    SNP   1500     45     ...
MA_Line_1       chr2        23456     AT   A    DEL   800      52     ...
MA_Line_2       chr1        34567     C    T    SNP   2000     60     ...
```

---

### STEP 6: Statistical Analysis

**Module:** `core/statistical_analysis.py`

**Purpose:** Calculate mutation rates and perform statistical comparisons.

| Input | Output |
|-------|--------|
| Called mutations | mutation_rates.tsv |
| Sample information | spectrum_analysis.tsv |
| | effect_sizes.tsv |
| | statistical_tests.json |

#### 6A: Mutation Rate Calculation

**For each sample:**
```
Rate = (Number of mutations) / (Callable sites × Generations)

Example:
  Mutations = 50
  Callable sites = 36,000,000 bp
  Generations = 25

  Rate = 50 / (36,000,000 × 25) = 5.6 × 10⁻⁸ per site per generation
```

**Confidence intervals:**
1. **Poisson CI:** Based on Poisson distribution (count data)
2. **Bootstrap CI:** Resample mutations 10,000 times
3. **Bayesian CI:** Gamma-Poisson conjugate prior

#### 6B: Mutation Spectrum Analysis

**Categories (6 types due to strand symmetry):**
```
C→A (= G→T)    Transversion
C→G (= G→C)    Transversion
C→T (= G→A)    Transition (often from deamination)
T→A (= A→T)    Transversion
T→C (= A→G)    Transition
T→G (= A→C)    Transversion
```

**Metrics:**
- Ts/Tv ratio (transitions / transversions)
- GC→AT bias
- Spectrum distribution

#### 6C: Genotype Comparisons

**Tests performed:**
| Test | Purpose |
|------|---------|
| Mann-Whitney U | Compare rates between genotypes |
| Welch's t-test | Parametric rate comparison |
| Kruskal-Wallis | Multi-group comparison |
| Kolmogorov-Smirnov | Compare to Poisson distribution |
| Shapiro-Wilk | Test normality |
| Permutation test | Non-parametric significance |

#### 6D: Effect Size Calculations

**Metrics:**
- **Cohen's d:** Standardized difference between means
- **Hedges' g:** Bias-corrected Cohen's d
- **Glass's delta:** Uses control group SD only

---

### STEP 7: Visualization & Reporting

**Module:** `visualization/` and `reporting/`

#### Individual Sample Plots (`visualization/sample_plots.py`)

| Plot | Shows |
|------|-------|
| Mutation spectrum | Bar chart of 6 mutation types |
| Chromosome distribution | Mutations per chromosome |
| Allele frequency | Histogram of variant allele frequencies |
| Depth distribution | Histogram of depths at mutation sites |
| Quality scatter | QUAL vs depth for each mutation |
| Quality dashboard | Combined multi-panel view |

#### Group Plots (`visualization/group_plots.py`)

| Plot | Shows |
|------|-------|
| Genotype comparison | Box plots of rates by genotype |
| Spectrum by genotype | Stacked bars comparing spectra |
| Accumulation curves | Mutations vs generation |
| Effect size forest | Forest plot of effect sizes |
| Replicate consistency | Variation among replicates |

#### Combined Plots (`visualization/combined_plots.py`)

| Plot | Shows |
|------|-------|
| Rate heatmap | All samples × chromosomes |
| Spectrum clustering | PCA of mutation spectra |
| Sample similarity | Hierarchical clustering dendrogram |
| Summary dashboard | Overview of all samples |

#### Statistical Plots (`visualization/statistical_plots.py`)

| Plot | Shows |
|------|-------|
| Depth fit | Histogram + fitted normal distribution |
| Strand bias fit | Histogram + fitted log-normal |
| Poisson Q-Q | Observed vs expected quantiles |
| P-value distribution | Histogram of p-values |
| Bootstrap distribution | Distribution of bootstrapped rates |
| Repeat INDEL fits | Logistic fits for each repeat unit |

---

## Muver Models Summary

| Model | Input | Output | Purpose |
|-------|-------|--------|---------|
| **Depth Distribution** | BAM pileup depths | Normal(μ, σ) parameters; Filtered regions | Filter variants in abnormal coverage |
| **Strand Bias** | Per-position strand counts | LogNormal(μ, σ) parameters | Filter strand-biased artifacts |
| **Depth Correction** | Depths near chromosome ends | Correction factors | Adjust for end-of-chromosome bias |
| **Repeat INDELs** | Repeat regions + INDEL counts | Logistic fit per repeat unit | Correct for polymerase slippage |
| **Composite Significance** | Multiple p-values | Single composite score | Combine evidence for confidence |
| **Subclonal Detection** | Allele frequencies | Genotype + subclonal calls | Detect low-frequency variants |

---

## Data Flow Summary

```
FASTQ → [Preprocess] → Clean FASTQ
                           ↓
Reference + Clean FASTQ → [Align] → BAM
                                      ↓
BAM (all samples) → [Call Variants] → VCF
                          ↓
             Ancestor BAM → [Muver Prep] → Distribution Parameters
                                                    ↓
VCF + Parameters → [Mutation Calling] → Filtered Mutations
                                               ↓
Mutations + Samples → [Statistics] → Rates, Tests, CIs
                                          ↓
All Results → [Reporting] → Plots + HTML Report
```

---

## Output Directory Structure

```
results/
├── qc/                      # Preprocessing outputs
│   ├── sample1_trimmed_R1.fastq.gz
│   ├── sample1_fastp.html
│   └── ...
├── bams/                    # Alignment outputs
│   ├── sample1.markdup.bam
│   ├── sample1.markdup.bai
│   └── ...
├── vcfs/                    # Variant calling outputs
│   ├── filtered.vcf.gz
│   ├── sample1.g.vcf.gz
│   └── ...
├── muver_data/              # Muver model outputs
│   ├── depth_distribution.txt
│   ├── strand_bias_distribution.txt
│   ├── repeats.bed
│   └── repeat_indel_fits.txt
├── mutations/               # Mutation calling outputs
│   ├── mutations.tsv
│   ├── mutations.vcf
│   └── mutations_summary.json
├── statistics/              # Statistical analysis outputs
│   ├── mutation_rates.tsv
│   ├── spectrum_analysis.tsv
│   ├── effect_sizes.tsv
│   └── statistical_tests.json
├── plots/                   # Visualization outputs
│   ├── individual/          # Per-sample plots
│   ├── groups/              # Group comparison plots
│   ├── combined/            # All-samples plots
│   └── statistics/          # Statistical diagnostic plots
└── reports/                 # Final reports
    └── analysis_report.html
```

---

## Running the Pipeline

```bash
# Full pipeline
python run_pipeline.py -c config.yaml

# Skip early steps (if BAM files exist)
python run_pipeline.py -c config.yaml --skip-preprocessing --skip-alignment

# Resume from mutation calling (if VCF exists)
python run_pipeline.py -c config.yaml --resume-from mutation

# Only regenerate statistics and reports
python run_pipeline.py -c config.yaml --resume-from stats

# Only regenerate reports and plots
python run_pipeline.py -c config.yaml --resume-from report
```
