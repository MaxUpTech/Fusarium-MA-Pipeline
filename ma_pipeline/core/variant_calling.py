#!/usr/bin/env python3
"""
Variant Calling Module for Fusarium MA Pipeline
GATK4 HaplotypeCaller with joint genotyping.
"""

import logging
import subprocess
import shutil
from pathlib import Path
from typing import List


class VariantCallingModule:
    """Variant calling using GATK4."""
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.VariantCalling')
        self.output_dir = Path(config['output']['directory']) / 'vcfs'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = config['resources'].get('max_threads', 8)
        self.reference = config['reference']['fasta']
    
    def create_dict(self):
        """Create sequence dictionary if needed."""
        dict_file = Path(self.reference).with_suffix('.dict')
        if not dict_file.exists():
            subprocess.run(['gatk', 'CreateSequenceDictionary',
                          '-R', self.reference], check=True)
    
    def call_variants_per_sample(self, sample) -> str:
        """Run HaplotypeCaller on a sample."""
        self.logger.info(f"Calling variants for {sample.sample_name}...")
        
        gvcf = self.output_dir / f"{sample.sample_name}.g.vcf.gz"
        subprocess.run([
            'gatk', 'HaplotypeCaller',
            '-R', self.reference,
            '-I', sample.bam_file,
            '-O', str(gvcf),
            '-ERC', 'GVCF',
            '--native-pair-hmm-threads', str(self.threads)
        ], check=True)
        
        return str(gvcf)
    
    def joint_genotype(self, samples: List) -> str:
        """Perform joint genotyping."""
        self.logger.info("Running joint genotyping...")
        
        # Create sample map
        sample_map = self.output_dir / 'sample_map.txt'
        with open(sample_map, 'w') as f:
            for sample in samples:
                f.write(f"{sample.sample_name}\t{sample.gvcf_file}\n")
        
        # Get intervals
        intervals = self._get_intervals()
        
        # GenomicsDB
        db_path = self.output_dir / 'genomicsdb'
        if db_path.exists():
            shutil.rmtree(db_path)
        
        subprocess.run([
            'gatk', 'GenomicsDBImport',
            '--sample-name-map', str(sample_map),
            '--genomicsdb-workspace-path', str(db_path),
            '-L', intervals
        ], check=True)
        
        # GenotypeGVCFs
        joint_vcf = self.output_dir / 'joint.vcf.gz'
        subprocess.run([
            'gatk', 'GenotypeGVCFs',
            '-R', self.reference,
            '-V', f'gendb://{db_path}',
            '-O', str(joint_vcf)
        ], check=True)
        
        return str(joint_vcf)
    
    def filter_variants(self, vcf: str) -> str:
        """Apply hard filters."""
        # Filter SNPs
        snp_vcf = self.output_dir / 'snps.vcf.gz'
        subprocess.run([
            'gatk', 'SelectVariants', '-R', self.reference, '-V', vcf,
            '--select-type-to-include', 'SNP', '-O', str(snp_vcf)
        ], check=True)
        
        snp_filtered = self.output_dir / 'snps.filtered.vcf.gz'
        subprocess.run([
            'gatk', 'VariantFiltration', '-R', self.reference, '-V', str(snp_vcf),
            '--filter-expression', 'QD < 2.0 || FS > 60.0 || MQ < 40.0',
            '--filter-name', 'SNP_FILTER', '-O', str(snp_filtered)
        ], check=True)
        
        # Filter INDELs
        indel_vcf = self.output_dir / 'indels.vcf.gz'
        subprocess.run([
            'gatk', 'SelectVariants', '-R', self.reference, '-V', vcf,
            '--select-type-to-include', 'INDEL', '-O', str(indel_vcf)
        ], check=True)
        
        indel_filtered = self.output_dir / 'indels.filtered.vcf.gz'
        subprocess.run([
            'gatk', 'VariantFiltration', '-R', self.reference, '-V', str(indel_vcf),
            '--filter-expression', 'QD < 2.0 || FS > 200.0',
            '--filter-name', 'INDEL_FILTER', '-O', str(indel_filtered)
        ], check=True)
        
        # Merge using bcftools concat (handles unsorted inputs better)
        merged = self.output_dir / 'filtered.vcf.gz'
        subprocess.run([
            'bcftools', 'concat', '-a', '-O', 'z',
            '-o', str(merged),
            str(snp_filtered), str(indel_filtered)
        ], check=True)
        
        # Index the merged VCF
        subprocess.run(['bcftools', 'index', '-t', str(merged)], check=True)
        
        return str(merged)
    
    def _get_intervals(self) -> str:
        """Get chromosome intervals from reference."""
        fai = f"{self.reference}.fai"
        intervals = []
        with open(fai, 'r') as f:
            for line in f:
                intervals.append(line.split('\t')[0])
        
        interval_file = self.output_dir / 'intervals.list'
        with open(interval_file, 'w') as f:
            f.write('\n'.join(intervals))
        return str(interval_file)
    
    def run(self, samples: List) -> str:
        """Run variant calling pipeline."""
        self.create_dict()
        
        # Call variants per sample
        for sample in samples:
            sample.gvcf_file = self.call_variants_per_sample(sample)
        
        # Joint genotype
        joint_vcf = self.joint_genotype(samples)
        
        # Filter
        filtered_vcf = self.filter_variants(joint_vcf)
        
        self.logger.info("Variant calling complete")
        return filtered_vcf
