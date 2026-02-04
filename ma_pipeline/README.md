# Fusarium graminearum MA Pipeline v2.0

**Enhanced Mutation Accumulation Analysis Pipeline**

Integrates muver statistical models for improved SNP/INDEL accuracy with comprehensive visualization.

---

## Directory Structure

```
pipeline/
├── run_pipeline.py          # Main entry point - run this to start
├── config.yaml              # Configuration - edit paths here
├── requirements.txt         # Python dependencies
├── sample_mapping.md        # Your sample mapping file
│
├── core/                    # Core pipeline modules
│   ├── preprocessing.py     # FASTQ QC and trimming
│   ├── alignment.py         # BWA-MEM2 alignment
│   ├── variant_calling.py   # GATK variant calling
│   ├── mutation_calling.py  # Enhanced mutation detection
│   ├── statistical_analysis.py  # Comprehensive statistics
│   └── sample_mapper.py     # Sample name mapping
│
├── muver/                   # Muver statistical models
│   ├── fitting.py           # Gaussian/logistic fitting
│   ├── bias_distribution.py # Strand bias (log-normal)
│   ├── depth_distribution.py # Depth filtering
│   ├── depth_correction.py  # Chromosome end correction
│   ├── repeat_indels.py     # Repeat INDEL correction
│   ├── composite_significance.py # Composite scoring
│   └── subclonal_detection.py    # Subclonal variants
│
├── visualization/           # Plotting modules
│   ├── sample_plots.py      # Individual sample plots
│   ├── group_plots.py       # Group comparisons
│   ├── combined_plots.py    # All samples together
│   └── statistical_plots.py # Statistical analysis plots
│
└── reporting/               # Report generation
    └── html_report.py       # Interactive HTML report
```

---

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Edit Configuration

Edit `config.yaml` with your paths:
- Reference genome path
- FASTQ directory
- Sample mapping file
- Output directory

### 3. Edit Sample Mapping

Edit `sample_mapping.md` with your sample information.

### 4. Run Pipeline

```bash
# Full pipeline
python run_pipeline.py -c config.yaml

# Skip early steps (if BAM/VCF already exist)
python run_pipeline.py -c config.yaml --resume-from mutation

# Only regenerate reports
python run_pipeline.py -c config.yaml --resume-from report
```

---

## Output Files

After running, results are in your output directory:

| Directory | Contents |
|-----------|----------|
| `mutations/` | Mutation calls (TSV, VCF) |
| `statistics/` | Rates, spectrum, effect sizes |
| `plots/individual/` | Per-sample visualizations |
| `plots/groups/` | Genotype comparisons |
| `plots/combined/` | All-sample analysis |
| `plots/statistics/` | Statistical test plots |
| `reports/` | Interactive HTML report |

---

## Key Features

### Muver Statistical Models
- **Strand bias filtering** - Log-normal distribution characterization
- **Depth filtering** - Normal distribution with CDF-based thresholds
- **Repeat INDEL correction** - Logistic fitting for repeat regions
- **Composite significance** - Combined statistical testing
- **Subclonal detection** - Maximum likelihood genotyping

### Enhanced Statistics
- Poisson uniformity test
- Kolmogorov-Smirnov test
- Shapiro-Wilk normality test
- Permutation tests
- Effect sizes (Cohen's d, Hedges' g)
- Bootstrap confidence intervals
- Bayesian credible intervals

### Comprehensive Visualization
- Individual sample dashboards
- Genotype comparison plots
- Mutation spectrum clustering
- Effect size forest plots
- Statistical diagnostic plots

---

## Citation

If you use this pipeline, please cite:
- The muver pipeline for statistical models
- GATK for variant calling
- BWA-MEM2 for alignment
