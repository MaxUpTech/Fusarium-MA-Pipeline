#!/usr/bin/env python3
"""
Sample Mapper Module for Fusarium MA Pipeline

Handles mapping between Illumina sample names and meaningful experiment names.
Supports Excel (.xlsx) and Markdown (.md) mapping files.

Example mapping file formats:

Excel (sample_mapping.xlsx):
| illumina_name                              | sample_name        | sample_type |
|--------------------------------------------|--------------------| ------------|
| 20087FL-23-01-01_S142_L006_R1_001.fastq.gz | Mlh1-0_parent      | ancestor    |
| 20087FL-23-01-02_S143_L006_R1_001.fastq.gz | mlh1.1_G25-D1      | ma_line     |

Markdown (sample_mapping.md):
| illumina_name | sample_name | sample_type |
|---------------|-------------|-------------|
| 20087FL-23-01-01_S142_L006 | Mlh1-0_parent | ancestor |
| 20087FL-23-01-02_S143_L006 | mlh1.1_G25-D1 | ma_line |
"""

import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class Sample:
    """Data class representing a sample in the experiment."""
    original_name: str  # Illumina name pattern
    sample_name: str    # Meaningful name (e.g., Mlh1-0_parent)
    sample_type: str    # 'ancestor' or 'ma_line'
    fastq_r1: str = ""
    fastq_r2: str = ""
    genotype: str = ""  # e.g., 'Mlh1', 'WT'
    generation: int = 0  # 0 for parent, 25 for G25
    replicate: str = ""  # e.g., 'D1', 'D2'
    
    def is_ancestor(self) -> bool:
        return self.sample_type == 'ancestor'
    
    def is_ma_line(self) -> bool:
        return self.sample_type == 'ma_line'


class SampleMapper:
    """
    Maps Illumina sample names to meaningful experiment names.
    
    Supports:
    - Excel mapping files (.xlsx)
    - Markdown mapping files (.md)
    - Auto-detection of FASTQ pairs
    - Parsing of sample metadata from names
    """
    
    def __init__(self, fastq_dir: str, mapping_file: str = None):
        self.fastq_dir = Path(fastq_dir)
        self.mapping_file = Path(mapping_file) if mapping_file else None
        self.logger = logging.getLogger('FusariumMA.SampleMapper')
        self.samples: Dict[str, Sample] = {}
        
    def load_mapping_from_excel(self, excel_file: str) -> Dict[str, dict]:
        """
        Load sample mapping from Excel file.
        
        Expected columns:
        - illumina_name: Pattern to match in FASTQ filename
        - sample_name: Meaningful name for the sample
        - sample_type: 'ancestor' or 'ma_line'
        
        Optional columns:
        - genotype, generation, replicate
        """
        try:
            import pandas as pd
        except ImportError:
            self.logger.error("pandas not installed. Run: pip install pandas openpyxl")
            raise
            
        df = pd.read_excel(excel_file)
        
        # Normalize column names
        df.columns = df.columns.str.lower().str.strip().str.replace(' ', '_')
        
        mapping = {}
        for _, row in df.iterrows():
            illumina_pattern = str(row.get('illumina_name', row.get('original_name', '')))
            sample_name = str(row.get('sample_name', row.get('name', '')))
            sample_type = str(row.get('sample_type', row.get('type', 'ma_line'))).lower()
            
            if illumina_pattern and sample_name:
                mapping[illumina_pattern] = {
                    'sample_name': sample_name,
                    'sample_type': sample_type,
                    'genotype': str(row.get('genotype', '')),
                    'generation': int(row.get('generation', 25 if sample_type == 'ma_line' else 0)),
                    'replicate': str(row.get('replicate', ''))
                }
                
        self.logger.info(f"Loaded {len(mapping)} sample mappings from Excel")
        return mapping
    
    def load_mapping_from_markdown(self, md_file: str) -> Dict[str, dict]:
        """
        Load sample mapping from Markdown table.
        
        Expected format:
        | illumina_name | sample_name | sample_type |
        |---------------|-------------|-------------|
        | pattern1      | name1       | ancestor    |
        | pattern2      | name2       | ma_line     |
        """
        mapping = {}
        
        with open(md_file, 'r') as f:
            lines = f.readlines()
        
        # Find table header
        header_idx = None
        headers = []
        for i, line in enumerate(lines):
            if '|' in line and 'illumina' in line.lower():
                headers = [h.strip().lower().replace(' ', '_') 
                          for h in line.split('|') if h.strip()]
                header_idx = i
                break
        
        if header_idx is None:
            self.logger.error("Could not find table header in Markdown file")
            return mapping
        
        # Skip separator line (|---|---|---|)
        data_start = header_idx + 2
        
        # Parse data rows
        for line in lines[data_start:]:
            if '|' not in line or line.strip().startswith('|--'):
                continue
            
            # Split and strip, but keep empty values
            parts = line.split('|')
            # Remove first and last empty parts from leading/trailing |
            if parts and not parts[0].strip():
                parts = parts[1:]
            if parts and not parts[-1].strip():
                parts = parts[:-1]
            values = [v.strip() for v in parts]
            
            if len(values) >= 3:  # At minimum need illumina_name, sample_name, sample_type
                # Pad values with empty strings if fewer than headers
                while len(values) < len(headers):
                    values.append('')
                row = dict(zip(headers, values))
                
                illumina_pattern = row.get('illumina_name', row.get('original_name', ''))
                sample_name = row.get('sample_name', row.get('name', ''))
                sample_type = row.get('sample_type', row.get('type', 'ma_line')).lower()
                
                if illumina_pattern and sample_name:
                    mapping[illumina_pattern] = {
                        'sample_name': sample_name,
                        'sample_type': sample_type,
                        'genotype': row.get('genotype', ''),
                        'generation': int(row.get('generation', 25 if sample_type == 'ma_line' else 0)),
                        'replicate': row.get('replicate', '')
                    }
        
        self.logger.info(f"Loaded {len(mapping)} sample mappings from Markdown")
        return mapping
    
    def find_fastq_files(self) -> List[Tuple[str, str, str]]:
        """
        Find all FASTQ files in the directory.
        
        Returns:
            List of tuples: (base_name, r1_path, r2_path)
        """
        fastq_files = {}
        
        # Find all FASTQ files
        patterns = ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']
        all_files = []
        for pattern in patterns:
            all_files.extend(self.fastq_dir.glob(pattern))
        
        for fq_file in all_files:
            name = fq_file.name
            
            # Determine if R1 or R2
            is_r1 = bool(re.search(r'_R1[_.]|_1\.f', name, re.IGNORECASE))
            is_r2 = bool(re.search(r'_R2[_.]|_2\.f', name, re.IGNORECASE))
            
            # Extract base name (remove R1/R2 and extension)
            base = re.sub(r'_R[12][_.].*$|_[12]\.f.*$', '', name, flags=re.IGNORECASE)
            
            if base not in fastq_files:
                fastq_files[base] = {'r1': '', 'r2': ''}
            
            if is_r1:
                fastq_files[base]['r1'] = str(fq_file.absolute())
            elif is_r2:
                fastq_files[base]['r2'] = str(fq_file.absolute())
        
        result = []
        for base, paths in fastq_files.items():
            if paths['r1']:  # At least R1 must exist
                result.append((base, paths['r1'], paths['r2']))
        
        self.logger.info(f"Found {len(result)} FASTQ pairs/files")
        return result
    
    def match_sample_to_mapping(self, fastq_base: str, 
                                 mapping: Dict[str, dict]) -> Optional[dict]:
        """
        Match a FASTQ file base name to a mapping entry.
        
        Uses flexible matching to handle variations in naming.
        """
        # Try exact match first
        if fastq_base in mapping:
            return mapping[fastq_base]
        
        # Try partial matching
        for pattern, info in mapping.items():
            # Check if pattern is contained in fastq_base
            if pattern in fastq_base:
                return info
            # Check if fastq_base is contained in pattern
            if fastq_base in pattern:
                return info
            # Try matching key parts (e.g., S142, L006)
            pattern_parts = re.findall(r'S\d+|L\d+|\d{5}FL-\d+-\d+-\d+', pattern)
            for part in pattern_parts:
                if part in fastq_base:
                    return info
        
        return None
    
    def parse_sample_name(self, sample_name: str) -> dict:
        """
        Parse meaningful sample name to extract metadata.
        
        Examples:
        - "Mlh1-0_parent" -> genotype=Mlh1, generation=0, type=ancestor
        - "mlh1.1_G25-D1" -> genotype=mlh1, generation=25, replicate=D1
        """
        info = {
            'genotype': '',
            'generation': 0,
            'replicate': '',
            'sample_type': 'ma_line'
        }
        
        # Check for parent/ancestor
        if 'parent' in sample_name.lower() or 'ancestor' in sample_name.lower():
            info['sample_type'] = 'ancestor'
            info['generation'] = 0
        
        # Extract generation (G25, G10, etc.)
        gen_match = re.search(r'G(\d+)', sample_name, re.IGNORECASE)
        if gen_match:
            info['generation'] = int(gen_match.group(1))
        
        # Extract replicate (D1, D2, R1, R2, etc.)
        rep_match = re.search(r'[_-]([DR]\d+)', sample_name, re.IGNORECASE)
        if rep_match:
            info['replicate'] = rep_match.group(1)
        
        # Extract genotype (first part before _ or -)
        geno_match = re.match(r'^([A-Za-z0-9]+)', sample_name)
        if geno_match:
            info['genotype'] = geno_match.group(1)
        
        return info
    
    def create_samples(self, mapping_file: str = None) -> List[Sample]:
        """
        Create Sample objects by matching FASTQ files to mapping.
        
        Args:
            mapping_file: Path to Excel or Markdown mapping file
            
        Returns:
            List of Sample objects
        """
        # Load mapping
        if mapping_file:
            self.mapping_file = Path(mapping_file)
        
        if self.mapping_file is None:
            raise ValueError("No mapping file provided")
        
        if self.mapping_file.suffix == '.xlsx':
            mapping = self.load_mapping_from_excel(str(self.mapping_file))
        elif self.mapping_file.suffix == '.md':
            mapping = self.load_mapping_from_markdown(str(self.mapping_file))
        else:
            raise ValueError(f"Unsupported mapping file format: {self.mapping_file.suffix}")
        
        # Find FASTQ files
        fastq_files = self.find_fastq_files()
        
        # Match and create samples
        samples = []
        unmatched = []
        
        for base, r1, r2 in fastq_files:
            match = self.match_sample_to_mapping(base, mapping)
            
            if match:
                # Parse additional info from sample name
                parsed = self.parse_sample_name(match['sample_name'])
                
                sample = Sample(
                    original_name=base,
                    sample_name=match['sample_name'],
                    sample_type=match.get('sample_type', parsed['sample_type']),
                    fastq_r1=r1,
                    fastq_r2=r2,
                    genotype=match.get('genotype', parsed['genotype']),
                    generation=match.get('generation', parsed['generation']),
                    replicate=match.get('replicate', parsed['replicate'])
                )
                samples.append(sample)
                self.samples[sample.sample_name] = sample
            else:
                unmatched.append(base)
        
        if unmatched:
            self.logger.warning(f"Unmatched FASTQ files: {unmatched}")
        
        # Validate: must have at least one ancestor
        ancestors = [s for s in samples if s.is_ancestor()]
        ma_lines = [s for s in samples if s.is_ma_line()]
        
        self.logger.info(f"Created {len(samples)} samples: {len(ancestors)} ancestor(s), {len(ma_lines)} MA lines")
        
        if not ancestors:
            self.logger.error("No ancestor/parent sample found!")
        
        return samples
    
    def get_ancestor(self) -> Optional[Sample]:
        """Get the ancestor/parent sample."""
        for sample in self.samples.values():
            if sample.is_ancestor():
                return sample
        return None
    
    def get_ma_lines(self) -> List[Sample]:
        """Get all MA line samples."""
        return [s for s in self.samples.values() if s.is_ma_line()]
    
    def export_sample_sheet(self, output_file: str):
        """Export samples to TSV format for pipeline input."""
        with open(output_file, 'w') as f:
            headers = ['sample_name', 'original_name', 'fastq_r1', 'fastq_r2', 
                      'sample_type', 'genotype', 'generation', 'replicate']
            f.write('\t'.join(headers) + '\n')
            
            for sample in self.samples.values():
                values = [
                    sample.sample_name,
                    sample.original_name,
                    sample.fastq_r1,
                    sample.fastq_r2,
                    sample.sample_type,
                    sample.genotype,
                    str(sample.generation),
                    sample.replicate
                ]
                f.write('\t'.join(values) + '\n')
        
        self.logger.info(f"Exported sample sheet to {output_file}")
