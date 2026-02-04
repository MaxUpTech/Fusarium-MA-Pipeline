# Fusarium graminearum MA Pipeline

**A comprehensive, end-to-end pipeline for mutation accumulation analysis, customized for your *Fusarium graminearum* experiment.**

This pipeline is designed to handle your specific Illumina naming convention, paired-end HiSeq data, and G25 experiment. It combines the best features of the JGI Resequencing Pipeline and the Muver tool to provide a complete and robust analysis of de novo mutations.

---

## Features

| Feature | Description |
|---|---|
| **Flexible Sample Mapping** | Use an Excel or Markdown file to map your complex Illumina names to meaningful sample names. |
| **Automated Workflow** | From raw FASTQ to final report in one command. |
| **Dependency Checker** | A script to check for all required tools and help you install them. |
| **Comprehensive QC** | Read trimming and quality control with `fastp`. |
| **Robust Alignment** | `bwa-mem2` alignment with duplicate marking using `GATK MarkDuplicates`. |
| **Accurate Variant Calling** | `GATK HaplotypeCaller` with joint genotyping for high sensitivity. |
| **MA-Specific Mutation Calling** | Identifies true de novo mutations using statistical tests (binomial, Fisher's exact) and filters out sequencing errors. |
| **Complete Statistical Analysis** | Calculates mutation rates, mutation spectrum (6-category), Ts/Tv ratio, INDEL size distribution, and more. |
| **Publication-Ready Outputs** | Generates an interactive HTML report, high-quality plots, and detailed TSV files. |

---

## How to Use

### Step 1: Check Dependencies

Before you start, run the dependency checker to make sure you have all the required tools. It can also help you install them.

```bash
cd /path/to/Fusarium_MA_Pipeline
./scripts/check_and_install_dependencies.sh
```

### Step 2: Prepare Your Mapping File

Copy the template `templates/sample_mapping_template.md` (or create an Excel file) and fill it in with your sample information.

**Example `sample_mapping.md`:**

| illumina_name | sample_name | sample_type |
|---|---|---|
| 20087FL-23-01-01_S142 | Mlh1-0_parent | ancestor |
| 20087FL-23-01-02_S143 | mlh1.1_G25-D1 | ma_line |
| ... | ... | ... |

### Step 3: Configure Your Pipeline

Copy the configuration template `config/config_template.yaml` to your working directory and edit the paths:

- `reference.fasta`: Path to your *Fusarium graminearum* reference genome.
- `samples.fastq_dir`: Path to the directory containing your FASTQ files.
- `samples.mapping_file`: Path to your sample mapping file from Step 2.

### Step 4: Run the Pipeline

Execute the main pipeline script with your configuration file:

```bash
python3 scripts/run_pipeline.py -c /path/to/your_config.yaml
```

That's it! The pipeline will run the complete analysis and generate all outputs in the directory specified in your config file.

---

## Understanding Your Output

Your results will be organized in the output directory:

- `qc/`: Quality control reports from `fastp`.
- `bams/`: Final, analysis-ready BAM files.
- `vcfs/`: VCF files from variant calling.
- `mutations/`: Your final list of de novo mutations.
  - `mutations.tsv`: The main result file with all significant mutations.
  - `all_candidates.tsv`: All potential mutations before final filtering.
  - `mutations.vcf`: Significant mutations in VCF format.
- `statistics/`: All statistical analysis results.
  - `mutation_rates.tsv`: Mutation rates per sample.
  - `mutation_spectrum.tsv`: Mutation spectrum counts.
  - `summary_statistics.json`: A summary of all key metrics.
- `plots/`: Publication-ready plots.
- `reports/`: The final interactive HTML report.

### Key Files to Check First

1.  **`reports/analysis_report.html`**: The best place to start. A full overview of your results with plots and tables.
2.  **`mutations/mutations.tsv`**: The definitive list of all de novo mutations found in your MA lines.
3.  **`statistics/mutation_rates.tsv`**: The calculated mutation rate for each of your MA lines.

---

## Pipeline Workflow

1.  **Preprocessing**: Trim adapters and low-quality bases from FASTQ files using `fastp`.
2.  **Alignment**: Align trimmed reads to the reference genome using `bwa-mem2`.
3.  **BAM Processing**: Sort and mark duplicate reads using `samtools` and `GATK MarkDuplicates`.
4.  **Variant Calling**: Call variants for each sample using `GATK HaplotypeCaller` in GVCF mode.
5.  **Joint Genotyping**: Combine all samples to create a joint-called VCF file using `GATK GenomicsDBImport` and `GenotypeGVCFs`.
6.  **Filtering**: Apply hard filters to remove low-quality variant calls.
7.  **Mutation Calling**: Compare each MA line to the ancestor to identify de novo mutations, applying statistical tests to ensure high confidence.
8.  **Statistical Analysis**: Calculate mutation rates, spectrum, and other key metrics.
9.  **Reporting**: Generate all final output files, plots, and the HTML report.
