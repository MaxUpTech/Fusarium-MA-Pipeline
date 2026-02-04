#!/usr/bin/env python3
"""
Fusarium graminearum MA Pipeline - Main Orchestrator

Complete pipeline for mutation accumulation analysis comparing
Parent strain to G25 MA lines.

Enhanced with muver statistical models for improved accuracy:
- Depth distribution filtering
- Strand bias distribution (log-normal)
- Depth correction for chromosome ends
- Repeat INDEL correction
- Composite significance scoring
- Subclonal variant detection

Comprehensive visualization:
- Individual sample plots
- Group comparison plots
- Combined analysis plots
- Statistical analysis plots

Author: Customized for Muhanad's Fusarium MA Experiment
"""

import argparse
import logging
import os
import sys
import yaml
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from core import (
    SampleMapper, Sample,
    PreprocessingModule,
    AlignmentModule,
    VariantCallingModule,
    EnhancedMutationCallingModule,
    EnhancedStatistics
)
from reporting import EnhancedReporting


class EnhancedFusariumMAPipeline:
    """
    Enhanced pipeline for Fusarium graminearum mutation accumulation analysis.
    
    Integrates muver statistical models for improved accuracy and
    comprehensive visualization capabilities.
    """
    
    def __init__(self, config_path: str):
        """Initialize pipeline with configuration."""
        self.config = self._load_config(config_path)
        self.setup_logging()
        self.setup_directories()
        
        # Store muver analysis data for reporting
        self.muver_data = {}
        
        self.logger.info("=" * 70)
        self.logger.info("Fusarium graminearum MA Pipeline - Enhanced")
        self.logger.info("Parent vs G25 Mutation Accumulation Analysis")
        self.logger.info("With Muver Statistical Models Integration")
        self.logger.info("=" * 70)
        
    def _load_config(self, config_path: str) -> dict:
        """Load configuration from YAML file."""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Set defaults for new configuration options
        config.setdefault('muver_models', {})
        config['muver_models'].setdefault('depth_distribution', {'enabled': True})
        config['muver_models'].setdefault('strand_bias', {'enabled': True})
        config['muver_models'].setdefault('depth_correction', {'enabled': True})
        config['muver_models'].setdefault('repeat_indels', {'enabled': True})
        config['muver_models'].setdefault('subclonal_detection', {'enabled': True})
        config['muver_models'].setdefault('composite_significance', {'enabled': True})
        
        config.setdefault('visualization', {})
        config['visualization'].setdefault('individual_plots', {'enabled': True})
        config['visualization'].setdefault('group_plots', {'enabled': True})
        config['visualization'].setdefault('combined_plots', {'enabled': True})
        config['visualization'].setdefault('statistical_plots', {'enabled': True})
        
        config.setdefault('statistics', {})
        config['statistics'].setdefault('additional_tests', {
            'kolmogorov_smirnov': True,
            'shapiro_wilk': True,
            'permutation_tests': True,
            'effect_size': True
        })
        
        return config
    
    def setup_logging(self):
        """Configure logging."""
        log_config = self.config.get('logging', {})
        level = getattr(logging, log_config.get('level', 'INFO'))
        
        # Create output directory if needed
        output_dir = Path(self.config['output']['directory'])
        output_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = output_dir / log_config.get('file', 'pipeline.log')
        
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger('FusariumMA')
    
    def setup_directories(self):
        """Create output directory structure."""
        output_dir = Path(self.config['output']['directory'])
        subdirs = [
            'qc', 'bams', 'vcfs', 'mutations', 'statistics',
            'plots/individual', 'plots/groups', 'plots/combined', 'plots/statistics',
            'reports', 'muver_data'
        ]
        for subdir in subdirs:
            (output_dir / subdir).mkdir(parents=True, exist_ok=True)
    
    def validate_inputs(self):
        """Validate that all input files exist."""
        # Check reference
        ref = self.config['reference']['fasta']
        if not Path(ref).exists():
            raise FileNotFoundError(f"Reference genome not found: {ref}")
        
        # Check FASTQ directory
        fastq_dir = self.config['samples']['fastq_dir']
        if not Path(fastq_dir).exists():
            raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")
        
        # Check mapping file
        mapping_file = self.config['samples']['mapping_file']
        if not Path(mapping_file).exists():
            raise FileNotFoundError(f"Mapping file not found: {mapping_file}")
        
        self.logger.info("Input validation passed")
    
    def load_samples(self) -> List:
        """Load and map samples from FASTQ files."""
        mapper = SampleMapper(
            fastq_dir=self.config['samples']['fastq_dir'],
            mapping_file=self.config['samples']['mapping_file']
        )
        
        samples = mapper.create_samples()
        
        # Export sample sheet for reference
        output_dir = Path(self.config['output']['directory'])
        mapper.export_sample_sheet(str(output_dir / 'sample_sheet.tsv'))
        
        return samples, mapper.get_ancestor(), mapper.get_ma_lines()
    
    def get_bam_files(self, samples: List) -> Dict[str, str]:
        """Get dictionary mapping sample names to BAM file paths."""
        bam_dir = Path(self.config['output']['directory']) / 'bams'
        bam_files = {}
        
        for sample in samples:
            # Try different naming conventions
            for suffix in ['.dedup.bam', '.markdup.bam', '.sorted.bam', '.bam']:
                bam_path = bam_dir / f"{sample.sample_name}{suffix}"
                if bam_path.exists():
                    bam_files[sample.sample_name] = str(bam_path)
                    break
        
        return bam_files
    
    def run_muver_preparation(self, samples: List, ancestor) -> Dict:
        """
        Run muver model preparation steps before mutation calling.
        
        This characterizes genome-wide distributions for:
        - Depth distribution
        - Strand bias distribution
        - Repeat regions
        """
        self.logger.info("=" * 70)
        self.logger.info("Muver Model Preparation")
        self.logger.info("=" * 70)
        
        muver_data = {}
        muver_config = self.config.get('muver_models', {})
        
        # Get BAM files
        bam_files = self.get_bam_files(samples)
        reference = self.config['reference']['fasta']
        
        if not bam_files:
            self.logger.warning("No BAM files found, skipping muver preparation")
            return muver_data
        
        # Use ancestor BAM for characterization (or first available)
        if ancestor and ancestor.sample_name in bam_files:
            char_bam = bam_files[ancestor.sample_name]
        else:
            char_bam = list(bam_files.values())[0]
        
        output_dir = Path(self.config['output']['directory']) / 'muver_data'
        
        # Depth distribution
        if muver_config.get('depth_distribution', {}).get('enabled', True):
            self.logger.info("Calculating depth distribution...")
            try:
                from muver_models.depth_distribution import DepthDistributionAnalyzer
                depth_analyzer = DepthDistributionAnalyzer(muver_config)
                muver_data['depth_distribution'] = depth_analyzer.calculate_from_bam(
                    char_bam, reference,
                    str(output_dir / 'depth_distribution.txt')
                )
                self.logger.info(f"  Depth: μ={muver_data['depth_distribution'].mu:.1f}, "
                               f"σ={muver_data['depth_distribution'].sigma:.1f}")
            except Exception as e:
                self.logger.warning(f"Depth distribution calculation failed: {e}")
        
        # Strand bias distribution
        if muver_config.get('strand_bias', {}).get('enabled', True):
            self.logger.info("Calculating strand bias distribution...")
            try:
                from muver_models.bias_distribution import BiasDistributionAnalyzer
                bias_analyzer = BiasDistributionAnalyzer(muver_config)
                muver_data['strand_bias'] = bias_analyzer.calculate_from_bam(
                    char_bam, reference,
                    str(output_dir / 'strand_bias_distribution.txt')
                )
                self.logger.info(f"  Strand bias: μ={muver_data['strand_bias'].mu:.3f}, "
                               f"σ={muver_data['strand_bias'].sigma:.3f}")
            except Exception as e:
                self.logger.warning(f"Strand bias calculation failed: {e}")
        
        # Repeat regions
        if muver_config.get('repeat_indels', {}).get('enabled', True):
            self.logger.info("Analyzing repeat regions...")
            try:
                from muver_models.repeat_indels import RepeatIndelAnalyzer
                repeat_analyzer = RepeatIndelAnalyzer(muver_config)
                
                # Check for existing repeat file or generate
                repeat_file = muver_config.get('repeat_indels', {}).get('repeat_file')
                if not repeat_file or not Path(repeat_file).exists():
                    repeat_file = str(output_dir / 'repeats.bed')
                    repeat_analyzer.generate_repeat_file(reference, repeat_file)
                else:
                    repeat_analyzer.load_repeat_file(repeat_file)
                
                muver_data['repeat_indels'] = repeat_analyzer.fit_repeat_indel_rates(
                    char_bam,
                    str(output_dir / 'repeat_indel_fits.txt'),
                    str(output_dir / 'repeat_indel')
                )
                self.logger.info("  Repeat INDEL fits calculated")
            except Exception as e:
                self.logger.warning(f"Repeat analysis failed: {e}")
        
        self.muver_data = muver_data
        return muver_data
    
    def run(self, skip_preprocessing: bool = False, skip_alignment: bool = False,
            skip_variant_calling: bool = False, resume_from: str = None):
        """
        Run the complete pipeline.
        
        Parameters
        ----------
        skip_preprocessing : bool
            Skip preprocessing step (use existing QC'd reads)
        skip_alignment : bool
            Skip alignment step (use existing BAM files)
        skip_variant_calling : bool
            Skip variant calling step (use existing VCF)
        resume_from : str
            Resume from a specific step ('muver', 'mutation', 'stats', 'report')
        """
        try:
            # Validate inputs
            self.validate_inputs()
            
            # Load samples
            self.logger.info("Loading samples...")
            samples, ancestor, ma_lines = self.load_samples()
            
            if not ancestor:
                raise ValueError("No ancestor sample found! Check your mapping file.")
            
            self.logger.info(f"Loaded {len(samples)} samples:")
            self.logger.info(f"  Ancestor: {ancestor.sample_name}")
            self.logger.info(f"  MA lines: {len(ma_lines)}")
            
            # Step 1: Preprocessing
            if not skip_preprocessing and resume_from not in ['muver', 'mutation', 'stats', 'report']:
                self.logger.info("=" * 70)
                self.logger.info("Step 1: Preprocessing (QC and trimming)")
                self.logger.info("=" * 70)
                preprocessor = PreprocessingModule(self.config)
                samples = preprocessor.run(samples)
            else:
                self.logger.info("Skipping preprocessing (using existing reads)")
            
            # Step 2: Alignment
            if not skip_alignment and resume_from not in ['muver', 'mutation', 'stats', 'report']:
                self.logger.info("=" * 70)
                self.logger.info("Step 2: Alignment")
                self.logger.info("=" * 70)
                aligner = AlignmentModule(self.config)
                samples = aligner.run(samples)
            else:
                self.logger.info("Skipping alignment (using existing BAM files)")
            
            # Step 3: Variant Calling
            if not skip_variant_calling and resume_from not in ['muver', 'mutation', 'stats', 'report']:
                self.logger.info("=" * 70)
                self.logger.info("Step 3: Variant Calling")
                self.logger.info("=" * 70)
                variant_caller = VariantCallingModule(self.config)
                filtered_vcf = variant_caller.run(samples)
            else:
                self.logger.info("Skipping variant calling (using existing VCF)")
                # Find existing VCF
                vcf_dir = Path(self.config['output']['directory']) / 'vcfs'
                vcf_files = list(vcf_dir.glob('*.filtered.vcf*'))
                if vcf_files:
                    filtered_vcf = str(vcf_files[0])
                else:
                    raise FileNotFoundError("No filtered VCF file found")
            
            # Step 4: Muver Model Preparation
            if resume_from not in ['mutation', 'stats', 'report']:
                muver_data = self.run_muver_preparation(samples, ancestor)
            else:
                self.logger.info("Skipping muver preparation")
                muver_data = {}
            
            # Step 5: Enhanced Mutation Calling
            if resume_from not in ['stats', 'report']:
                self.logger.info("=" * 70)
                self.logger.info("Step 5: Enhanced Mutation Calling (with Muver Models)")
                self.logger.info("=" * 70)
                mutation_caller = EnhancedMutationCallingModule(self.config)
                
                # Pass BAM files for muver initialization if not already done
                bam_files = self.get_bam_files(samples) if not muver_data else None
                mutations = mutation_caller.run(samples, ancestor, filtered_vcf, bam_files)
            else:
                self.logger.info("Skipping mutation calling")
                # Load existing mutations
                mutations_file = Path(self.config['output']['directory']) / 'mutations' / 'mutations.tsv'
                if mutations_file.exists():
                    # Would need to parse TSV back to Mutation objects
                    mutations = []  # Placeholder
                else:
                    raise FileNotFoundError("No mutations file found")
            
            # Step 6: Enhanced Statistical Analysis
            if resume_from not in ['report']:
                self.logger.info("=" * 70)
                self.logger.info("Step 6: Enhanced Statistical Analysis")
                self.logger.info("=" * 70)
                stats_module = EnhancedStatistics(self.config)
                stats = stats_module.run(mutations, samples, ancestor)
            else:
                self.logger.info("Skipping statistical analysis")
                stats = {}
            
            # Step 7: Comprehensive Reporting and Visualization
            self.logger.info("=" * 70)
            self.logger.info("Step 7: Generating Reports and Visualizations")
            self.logger.info("=" * 70)
            reporter = EnhancedReporting(self.config)
            reporter.run(mutations, stats, samples, ancestor, self.muver_data)
            
            # Summary
            self.logger.info("=" * 70)
            self.logger.info("Pipeline Complete!")
            self.logger.info("=" * 70)
            self.logger.info(f"Total mutations found: {len(mutations)}")
            self.logger.info(f"Output directory: {self.config['output']['directory']}")
            self.logger.info("")
            self.logger.info("Key output files:")
            self.logger.info("  Mutations:")
            self.logger.info("    - mutations/mutations.tsv")
            self.logger.info("    - mutations/mutations.vcf")
            self.logger.info("    - mutations/all_candidates.tsv")
            self.logger.info("  Statistics:")
            self.logger.info("    - statistics/mutation_rates.tsv")
            self.logger.info("    - statistics/mutation_spectrum.tsv")
            self.logger.info("    - statistics/effect_sizes.tsv")
            self.logger.info("    - statistics/summary_statistics.json")
            self.logger.info("  Plots:")
            self.logger.info("    - plots/individual/{sample}/")
            self.logger.info("    - plots/groups/")
            self.logger.info("    - plots/combined/")
            self.logger.info("    - plots/statistics/")
            self.logger.info("  Reports:")
            self.logger.info("    - reports/analysis_report.html")
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            raise


# Backward compatibility
FusariumMAPipeline = EnhancedFusariumMAPipeline


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Fusarium graminearum MA Pipeline - Enhanced Mutation Accumulation Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline
  python run_pipeline.py -c config.yaml
  
  # Resume from mutation calling (skip alignment, variant calling)
  python run_pipeline.py -c config.yaml --resume-from mutation
  
  # Skip preprocessing and alignment
  python run_pipeline.py -c config.yaml --skip-preprocessing --skip-alignment
  
  # Run only reporting (with existing results)
  python run_pipeline.py -c config.yaml --resume-from report

For help creating a config file, see config_template.yaml
        """
    )
    
    parser.add_argument(
        '-c', '--config',
        required=True,
        help='Path to configuration YAML file'
    )
    
    parser.add_argument(
        '--skip-preprocessing',
        action='store_true',
        help='Skip preprocessing step (use existing QC\'d reads)'
    )
    
    parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip alignment step (use existing BAM files)'
    )
    
    parser.add_argument(
        '--skip-variant-calling',
        action='store_true',
        help='Skip variant calling step (use existing VCF)'
    )
    
    parser.add_argument(
        '--resume-from',
        choices=['muver', 'mutation', 'stats', 'report'],
        help='Resume pipeline from a specific step'
    )
    
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='Fusarium MA Pipeline v2.0.0 (Enhanced with Muver Models)'
    )
    
    args = parser.parse_args()
    
    # Run pipeline
    pipeline = EnhancedFusariumMAPipeline(args.config)
    pipeline.run(
        skip_preprocessing=args.skip_preprocessing,
        skip_alignment=args.skip_alignment,
        skip_variant_calling=args.skip_variant_calling,
        resume_from=args.resume_from
    )


if __name__ == '__main__':
    main()
