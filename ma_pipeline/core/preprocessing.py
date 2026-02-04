#!/usr/bin/env python3
"""
Preprocessing Module for Fusarium MA Pipeline
Quality control and read trimming using fastp.
"""

import logging
import subprocess
from pathlib import Path
from typing import List


class PreprocessingModule:
    """Quality control and preprocessing of FASTQ files."""
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.Preprocessing')
        self.output_dir = Path(config['output']['directory']) / 'qc'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.threads = config['resources'].get('max_threads', 8)
    
    def run_fastp(self, sample) -> tuple:
        """Run fastp on a sample."""
        self.logger.info(f"Running fastp on {sample.sample_name}...")
        
        r1_out = self.output_dir / f"{sample.sample_name}_trimmed_R1.fastq.gz"
        r2_out = self.output_dir / f"{sample.sample_name}_trimmed_R2.fastq.gz"
        html_report = self.output_dir / f"{sample.sample_name}_fastp.html"
        json_report = self.output_dir / f"{sample.sample_name}_fastp.json"
        
        cmd = [
            'fastp',
            '-i', sample.fastq_r1,
            '-o', str(r1_out),
            '-h', str(html_report),
            '-j', str(json_report),
            '-w', str(self.threads),
            '-q', '20',  # Quality threshold
            '-l', '50',  # Min length
            '--detect_adapter_for_pe'
        ]
        
        if sample.fastq_r2:
            cmd.extend(['-I', sample.fastq_r2, '-O', str(r2_out)])
        
        subprocess.run(cmd, check=True)
        
        return (str(r1_out), str(r2_out) if sample.fastq_r2 else '')
    
    def run(self, samples: List) -> List:
        """Run preprocessing on all samples."""
        self.logger.info(f"Preprocessing {len(samples)} samples...")
        
        for sample in samples:
            r1, r2 = self.run_fastp(sample)
            sample.fastq_r1 = r1
            sample.fastq_r2 = r2
        
        self.logger.info("Preprocessing complete")
        return samples
