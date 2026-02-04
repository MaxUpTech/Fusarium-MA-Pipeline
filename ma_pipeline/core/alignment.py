#!/usr/bin/env python3
"""
Alignment Module for Fusarium MA Pipeline
BWA-MEM2 alignment with duplicate marking.
"""

import logging
import subprocess
import os
from pathlib import Path
from typing import List


class AlignmentModule:
    """Alignment of reads to reference genome."""
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.Alignment')
        self.output_dir = Path(config['output']['directory']) / 'bams'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = config['resources'].get('max_threads', 8)
        self.reference = config['reference']['fasta']
    
    def index_reference(self):
        """Index reference genome if needed."""
        if not Path(f"{self.reference}.bwt.2bit.64").exists():
            self.logger.info("Indexing reference genome...")
            subprocess.run(['bwa-mem2', 'index', self.reference], check=True)
        
        if not Path(f"{self.reference}.fai").exists():
            subprocess.run(['samtools', 'faidx', self.reference], check=True)
    
    def align_sample(self, sample) -> str:
        """Align a single sample."""
        self.logger.info(f"Aligning {sample.sample_name}...")
        
        # Align
        sam_file = self.output_dir / f"{sample.sample_name}.sam"
        rg = f"@RG\\tID:{sample.sample_name}\\tSM:{sample.sample_name}\\tPL:ILLUMINA"
        
        cmd = ['bwa-mem2', 'mem', '-t', str(self.threads), '-R', rg,
               self.reference, sample.fastq_r1]
        if sample.fastq_r2:
            cmd.append(sample.fastq_r2)
        
        with open(sam_file, 'w') as out:
            subprocess.run(cmd, stdout=out, check=True)
        
        # Sort
        sorted_bam = self.output_dir / f"{sample.sample_name}.sorted.bam"
        subprocess.run(['samtools', 'sort', '-@', str(self.threads),
                       '-o', str(sorted_bam), str(sam_file)], check=True)
        os.remove(sam_file)
        
        # Mark duplicates
        dedup_bam = self.output_dir / f"{sample.sample_name}.dedup.bam"
        metrics = self.output_dir / f"{sample.sample_name}.metrics.txt"
        subprocess.run([
            'gatk', 'MarkDuplicates',
            '-I', str(sorted_bam),
            '-O', str(dedup_bam),
            '-M', str(metrics),
            '--CREATE_INDEX', 'true'
        ], check=True)
        os.remove(sorted_bam)
        
        return str(dedup_bam)
    
    def run(self, samples: List) -> List:
        """Run alignment on all samples."""
        self.index_reference()
        
        for sample in samples:
            sample.bam_file = self.align_sample(sample)
        
        self.logger.info("Alignment complete")
        return samples
