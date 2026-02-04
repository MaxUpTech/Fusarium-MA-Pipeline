# Sample Mapping Template for Fusarium MA Experiment

## Instructions

Fill in this table with your sample information. The pipeline will use this to map
your Illumina FASTQ file names to meaningful sample names.

**Required columns:**
- `illumina_name`: Part of the FASTQ filename that uniquely identifies the sample (e.g., S142, or the full prefix)
- `sample_name`: Your meaningful name for the sample
- `sample_type`: Either `ancestor` (for parent) or `ma_line` (for G25 samples)

**Optional columns:**
- `genotype`: The genetic background (e.g., Mlh1, WT, Msh2)
- `generation`: Number of generations (0 for parent, 25 for G25)
- `replicate`: Replicate identifier (e.g., D1, D2)

## Sample Mapping Table

| illumina_name | sample_name | sample_type | genotype | generation | replicate |
|---------------|-------------|-------------|----------|------------|-----------|
| 20087FL-23-01-01_S142_L006 | Mlh1-0_parent | ancestor | Mlh1 | 0 | |
| 20087FL-23-01-02_S143_L006 | mlh1.1_G25-D1 | ma_line | Mlh1 | 25 | D1 |
| 20087FL-23-01-03_S144_L006 | mlh1.2_G25-D2 | ma_line | Mlh1 | 25 | D2 |
| 20087FL-23-01-04_S145_L006 | mlh1.3_G25-D3 | ma_line | Mlh1 | 25 | D3 |

## Tips

1. You don't need to include the full filename - just enough to uniquely identify each sample
2. The `S###` part (e.g., S142) is often sufficient for matching
3. Make sure exactly ONE sample has `sample_type` = `ancestor`
4. All other samples should have `sample_type` = `ma_line`

## Example with Multiple Genotypes

If you have multiple genetic backgrounds in your experiment:

| illumina_name | sample_name | sample_type | genotype | generation | replicate |
|---------------|-------------|-------------|----------|------------|-----------|
| S142 | WT_parent | ancestor | WT | 0 | |
| S143 | WT_G25-1 | ma_line | WT | 25 | 1 |
| S144 | WT_G25-2 | ma_line | WT | 25 | 2 |
| S145 | Mlh1_parent | ancestor | Mlh1 | 0 | |
| S146 | Mlh1_G25-1 | ma_line | Mlh1 | 25 | 1 |
| S147 | Mlh1_G25-2 | ma_line | Mlh1 | 25 | 2 |
