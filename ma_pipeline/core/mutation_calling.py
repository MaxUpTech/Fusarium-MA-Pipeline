#!/usr/bin/env python3
"""
Enhanced Mutation Calling Module for Fusarium MA Pipeline

MA-specific mutation calling with comprehensive statistical validation.
Integrates all statistical models from the muver pipeline:
- Depth distribution filtering
- Strand bias distribution (log-normal)
- Depth correction for chromosome ends
- Repeat INDEL correction
- Composite significance scoring
- Subclonal variant detection
"""

import logging
import math
import gzip
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from scipy.stats import binom, chi2, fisher_exact
import numpy as np

# Import muver models
from muver import (
    BiasDistributionAnalyzer,
    BiasDistributionResult,
    DepthDistributionAnalyzer,
    DepthDistributionResult,
    DepthCorrector,
    RepeatIndelAnalyzer,
    RepeatIndelResult,
    CompositeSignificanceCalculator,
    SubclonalDetector,
    get_repeat_adjustment_value,
)


@dataclass
class Mutation:
    """Data class for a called mutation with enhanced muver metrics."""
    chromosome: str
    position: int
    ref_allele: str
    alt_allele: str
    sample_name: str
    mutation_type: str  # 'SNP', 'INS', 'DEL'
    
    # Genotypes
    ancestor_genotype: str = ""
    sample_genotype: str = ""
    
    # Depths (raw)
    ancestor_ref_depth: int = 0
    ancestor_alt_depth: int = 0
    sample_ref_depth: int = 0
    sample_alt_depth: int = 0
    
    # Corrected depths (after chromosome end correction)
    corrected_sample_depth: float = 0.0
    corrected_ancestor_depth: float = 0.0
    
    # Strand-specific depths
    sample_ref_forward: int = 0
    sample_ref_reverse: int = 0
    sample_alt_forward: int = 0
    sample_alt_reverse: int = 0
    ancestor_ref_forward: int = 0
    ancestor_ref_reverse: int = 0
    ancestor_alt_forward: int = 0
    ancestor_alt_reverse: int = 0
    
    # Quality metrics
    quality: float = 0.0
    
    # Basic statistical tests
    p_value_binomial: float = 1.0
    p_value_chi_square: float = 1.0
    p_value_fisher: float = 1.0
    
    # Muver enhanced tests
    p_value_binomial_forward: float = 1.0
    p_value_binomial_reverse: float = 1.0
    composite_score: float = 0.0
    
    # Strand bias (enhanced)
    strand_bias: float = 0.0
    strand_bias_pvalue: float = 1.0  # From log-normal distribution
    variant_distance_bias: float = 0.0
    
    # Repeat context
    in_repeat: bool = False
    repeat_unit: str = ""
    repeat_tract_length: int = 0
    repeat_adjustment: float = 0.0
    
    # Subclonal detection
    is_subclonal: bool = False
    subclonal_frequency: float = 0.0
    subclonal_allele: str = ""
    
    # Flags
    is_significant: bool = False
    passes_depth_filter: bool = True
    passes_strand_filter: bool = True
    passes_depth_distribution_filter: bool = True
    passes_composite_test: bool = False
    
    # Context
    trinucleotide_context: str = ""
    
    # Annotation (if available)
    gene: str = ""
    effect: str = ""
    impact: str = ""
    
    def to_dict(self) -> dict:
        return asdict(self)


class EnhancedMutationCallingModule:
    """
    Enhanced MA-specific mutation calling module with muver integration.
    
    Identifies de novo mutations by comparing MA lines to ancestor.
    Uses multiple statistical models from muver for improved accuracy:
    
    1. Depth distribution filtering - removes regions with abnormal depth
    2. Strand bias (log-normal) - filters variants with strand imbalance
    3. Depth correction - corrects for chromosome end bias
    4. Repeat INDEL correction - adjusts for repeat-associated errors
    5. Composite significance - combines multiple tests
    6. Subclonal detection - identifies subclonal variants
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.MutationCalling')
        self.output_dir = Path(config['output']['directory']) / 'mutations'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create muver_data output directory
        self.muver_dir = Path(config['output']['directory']) / 'muver_data'
        self.muver_dir.mkdir(parents=True, exist_ok=True)
        
        # Basic parameters
        mc_config = config.get('mutation_calling', {})
        self.min_depth = config.get('filtering', {}).get('min_depth', 20)
        self.min_af = mc_config.get('min_allele_frequency', 0.8)
        self.p_threshold = mc_config.get('p_value_threshold', 0.01)
        self.fwer = mc_config.get('fwer', 0.01)
        self.strand_bias_threshold = mc_config.get('strand_bias_threshold', 2.0)
        
        # Muver model parameters
        muver_config = config.get('muver_models', {})
        
        # Initialize muver analyzers
        self.bias_analyzer: Optional[BiasDistributionAnalyzer] = None
        self.depth_analyzer: Optional[DepthDistributionAnalyzer] = None
        self.depth_corrector: Optional[DepthCorrector] = None
        self.repeat_analyzer: Optional[RepeatIndelAnalyzer] = None
        self.composite_calculator: Optional[CompositeSignificanceCalculator] = None
        self.subclonal_detector: Optional[SubclonalDetector] = None
        
        # Muver results
        self.bias_result: Optional[BiasDistributionResult] = None
        self.depth_result: Optional[DepthDistributionResult] = None
        self.repeat_result: Optional[RepeatIndelResult] = None
        
        # Filtered positions
        self.filtered_positions: Set[Tuple[str, int]] = set()
        
        # Enable/disable muver features
        self.use_depth_distribution = muver_config.get('depth_distribution', {}).get('enabled', True)
        self.use_strand_bias = muver_config.get('strand_bias', {}).get('enabled', True)
        self.use_depth_correction = muver_config.get('depth_correction', {}).get('enabled', True)
        self.use_repeat_indels = muver_config.get('repeat_indels', {}).get('enabled', True)
        self.use_composite_significance = muver_config.get('composite_significance', {}).get('enabled', True)
        self.use_subclonal_detection = muver_config.get('subclonal_detection', {}).get('enabled', True)
        
    def initialize_muver_models(self, bam_files: Dict[str, str], reference: str):
        """
        Initialize and run muver model preparation steps.
        
        Parameters
        ----------
        bam_files : dict
            Dictionary mapping sample names to BAM file paths
        reference : str
            Path to reference FASTA file
        """
        self.logger.info("Initializing muver statistical models...")
        
        # Get first BAM file for characterization (typically use control/ancestor)
        first_bam = list(bam_files.values())[0] if bam_files else None
        
        if not first_bam:
            self.logger.warning("No BAM files provided, skipping muver model initialization")
            return
        
        # Initialize depth distribution analyzer
        if self.use_depth_distribution:
            self.logger.info("Calculating depth distribution...")
            self.depth_analyzer = DepthDistributionAnalyzer(self.config.get('muver_models', {}))
            self.depth_result = self.depth_analyzer.calculate_from_bam(
                first_bam, reference,
                str(self.muver_dir / 'depth_distribution.txt')
            )
        
        # Initialize strand bias analyzer
        if self.use_strand_bias:
            self.logger.info("Calculating strand bias distribution...")
            self.bias_analyzer = BiasDistributionAnalyzer(self.config.get('muver_models', {}))
            self.bias_result = self.bias_analyzer.calculate_from_bam(
                first_bam, reference,
                str(self.muver_dir / 'strand_bias_distribution.txt')
            )
        
        # Initialize depth corrector
        if self.use_depth_correction:
            self.logger.info("Initializing depth correction...")
            self.depth_corrector = DepthCorrector(self.config.get('muver_models', {}))
            # Load chromosome sizes from reference
            fai_file = reference + '.fai'
            if Path(fai_file).exists():
                self.depth_corrector.load_chromosome_sizes(fai_file)
        
        # Initialize repeat indel analyzer
        if self.use_repeat_indels:
            self.logger.info("Initializing repeat INDEL analyzer...")
            self.repeat_analyzer = RepeatIndelAnalyzer(self.config.get('muver_models', {}))
            
            # Check for existing repeat file or generate
            repeat_file = self.config.get('muver_models', {}).get('repeat_file')
            if repeat_file and Path(repeat_file).exists():
                self.repeat_analyzer.load_repeat_file(repeat_file)
            else:
                # Generate repeat file from reference
                repeat_file = str(self.muver_dir / 'repeats.bed')
                self.repeat_analyzer.generate_repeat_file(reference, repeat_file)
            
            # Calculate repeat indel rates from first BAM
            self.repeat_result = self.repeat_analyzer.fit_repeat_indel_rates(
                first_bam,
                str(self.muver_dir / 'repeat_indel_fits.txt'),
                str(self.muver_dir / 'repeat_indel_plots')
            )
        
        # Initialize composite significance calculator
        if self.use_composite_significance:
            self.composite_calculator = CompositeSignificanceCalculator(
                self.config.get('muver_models', {}).get('composite_significance', {})
            )
        
        # Initialize subclonal detector
        if self.use_subclonal_detection:
            self.subclonal_detector = SubclonalDetector(
                self.config.get('muver_models', {}).get('subclonal_detection', {})
            )
            if self.bias_result:
                self.subclonal_detector.set_strand_bias_std(self.bias_result.sigma)
            if self.repeat_result:
                self.subclonal_detector.set_repeat_fits(self.repeat_result.fits)
        
        self.logger.info("Muver model initialization complete")
    
    def parse_vcf(self, vcf_file: str) -> List[dict]:
        """Parse VCF file and extract variant information with strand counts."""
        variants = []
        open_func = gzip.open if vcf_file.endswith('.gz') else open
        mode = 'rt' if vcf_file.endswith('.gz') else 'r'
        
        with open_func(vcf_file, mode) as f:
            sample_names = []
            
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    sample_names = fields[9:]
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue
                    
                chrom, pos, _, ref, alt, qual, filt, info, fmt = fields[:9]
                sample_data = fields[9:]
                
                # Skip filtered variants (unless PASS or .)
                if filt not in ['PASS', '.']:
                    continue
                
                variant = {
                    'chromosome': chrom,
                    'position': int(pos),
                    'ref_allele': ref,
                    'alt_alleles': alt.split(','),
                    'alleles': [ref] + alt.split(','),
                    'quality': float(qual) if qual != '.' else 0,
                    'filter': filt,
                    'info': self._parse_info(info),
                    'format': fmt.split(':'),
                    'samples': {},
                    'strand_allele_counts': {}
                }
                
                for sample_name, data in zip(sample_names, sample_data):
                    parsed = self._parse_sample_data(variant['format'], data)
                    variant['samples'][sample_name] = parsed
                    
                    # Extract strand allele counts for muver
                    variant['strand_allele_counts'][sample_name] = \
                        self._extract_strand_allele_counts(parsed, variant['alleles'])
                
                # Check if in repeat region
                if self.repeat_analyzer:
                    repeat = self.repeat_analyzer.get_repeat_at_position(chrom, int(pos))
                    if repeat:
                        variant['intersected_repeat'] = {
                            'unit': repeat.unit,
                            'tract_length': repeat.tract_length,
                            'sequence': repeat.sequence
                        }
                
                variants.append(variant)
        
        self.logger.info(f"Parsed {len(variants)} variants from VCF")
        return variants
    
    def _parse_info(self, info_str: str) -> dict:
        """Parse INFO field."""
        info = {}
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info[key] = value
            else:
                info[item] = True
        return info
    
    def _parse_sample_data(self, format_fields: List[str], data: str) -> dict:
        """Parse sample genotype data."""
        values = data.split(':')
        return dict(zip(format_fields, values))
    
    def _extract_strand_allele_counts(self, sample_data: dict, 
                                       alleles: List[str]) -> Dict[str, Dict[str, int]]:
        """Extract strand-specific allele counts for muver composite test."""
        counts = {}
        
        # Try SAC (Strand Allele Counts) first
        if 'SAC' in sample_data:
            sac = sample_data['SAC'].split(',')
            for i, allele in enumerate(alleles):
                if i * 2 + 1 < len(sac):
                    counts[allele] = {
                        'forward': int(sac[i * 2]) if sac[i * 2] != '.' else 0,
                        'reverse': int(sac[i * 2 + 1]) if sac[i * 2 + 1] != '.' else 0
                    }
                else:
                    counts[allele] = {'forward': 0, 'reverse': 0}
        
        # Try AD (Allele Depths) - split evenly between strands
        elif 'AD' in sample_data:
            ad = sample_data['AD'].split(',')
            for i, allele in enumerate(alleles):
                if i < len(ad) and ad[i] != '.':
                    depth = int(ad[i])
                    counts[allele] = {
                        'forward': depth // 2,
                        'reverse': depth - depth // 2
                    }
                else:
                    counts[allele] = {'forward': 0, 'reverse': 0}
        
        else:
            for allele in alleles:
                counts[allele] = {'forward': 0, 'reverse': 0}
        
        return counts
    
    def get_allele_depths(self, sample_data: dict) -> Tuple[int, int, int, int]:
        """
        Extract allele depths from sample data.
        
        Returns: (ref_forward, ref_reverse, alt_forward, alt_reverse)
        """
        # Try SAC (Strand Allele Counts) first
        if 'SAC' in sample_data:
            sac = sample_data['SAC'].split(',')
            if len(sac) >= 4:
                return (int(sac[0]) if sac[0] != '.' else 0,
                        int(sac[1]) if sac[1] != '.' else 0,
                        int(sac[2]) if sac[2] != '.' else 0,
                        int(sac[3]) if sac[3] != '.' else 0)
        
        # Try AD (Allele Depths)
        if 'AD' in sample_data:
            ad = sample_data['AD'].split(',')
            ref_depth = int(ad[0]) if ad[0] != '.' else 0
            alt_depth = int(ad[1]) if len(ad) > 1 and ad[1] != '.' else 0
            # No strand info, split evenly
            return (ref_depth // 2, ref_depth - ref_depth // 2,
                    alt_depth // 2, alt_depth - alt_depth // 2)
        
        # Fall back to DP
        if 'DP' in sample_data:
            dp = int(sample_data['DP']) if sample_data['DP'] != '.' else 0
            return (dp // 2, dp - dp // 2, 0, 0)
        
        return (0, 0, 0, 0)
    
    def binomial_test_per_strand(self, sample_counts: Dict[str, int],
                                  control_counts: Dict[str, int]) -> Tuple[float, float]:
        """
        Perform binomial tests for each strand (muver method).
        
        For de novo mutations, we test if sample has significantly more
        alt reads than expected by chance (background error rate).
        
        Returns: (p_forward, p_reverse)
        """
        sample_total = sample_counts.get('forward', 0) + sample_counts.get('reverse', 0)
        control_total = control_counts.get('forward', 0) + control_counts.get('reverse', 0)
        
        # Background error rate when no control alt reads
        background_error = 0.01
        
        p_values = []
        
        for strand in ['forward', 'reverse']:
            sample_value = sample_counts.get(strand, 0)
            control_value = control_counts.get(strand, 0)
            
            if sample_value > 0 and sample_total > 0:
                # If control has alt reads, use that frequency
                # Otherwise use background error rate
                if control_total > 0 and control_value > 0:
                    expected_freq = float(control_value) / control_total
                else:
                    expected_freq = background_error
                
                # Test: probability of seeing this many or more alt reads by chance
                p_val = 1 - binom.cdf(sample_value - 1, sample_total, expected_freq)
            else:
                p_val = 1.0
            
            p_values.append(p_val)
        
        return tuple(p_values)
    
    def calculate_composite_significance(self, p_binom_fwd: float, p_binom_rev: float,
                                         p_chi: float) -> Tuple[float, bool]:
        """
        Calculate composite significance score (muver method).
        
        Returns: (composite_score, is_significant)
        """
        min_p = sys.float_info.min
        p_binom_fwd = max(p_binom_fwd, min_p)
        p_binom_rev = max(p_binom_rev, min_p)
        p_chi = max(p_chi, min_p)
        
        # Euclidean distance in -log10(p) space
        composite = math.sqrt(
            (-math.log10(p_binom_fwd)) ** 2 +
            (-math.log10(p_binom_rev)) ** 2 +
            (-math.log10(p_chi)) ** 2
        )
        
        log_threshold = -math.log10(self.p_threshold)
        individual_threshold = 0.1
        
        is_significant = (
            composite >= log_threshold and
            p_binom_fwd < individual_threshold and
            p_binom_rev < individual_threshold and
            p_chi < individual_threshold
        )
        
        return composite, is_significant
    
    def binomial_test(self, sample_alt: int, sample_total: int,
                      ancestor_alt: int, ancestor_total: int) -> float:
        """Binomial test comparing sample to ancestor allele frequency."""
        if sample_total == 0 or ancestor_total == 0:
            return 1.0
        
        expected_freq = (ancestor_alt + 0.5) / (ancestor_total + 1)
        
        if ancestor_alt == 0:
            expected_freq = 0.01  # Background error rate
        
        p_value = 1 - binom.cdf(sample_alt - 1, sample_total, expected_freq)
        
        return p_value
    
    def chi_square_test(self, sample_ref: int, sample_alt: int,
                        ancestor_ref: int, ancestor_alt: int) -> float:
        """Chi-square test comparing allele distributions."""
        observed = np.array([[sample_ref, sample_alt],
                            [ancestor_ref, ancestor_alt]])
        
        if observed.sum() < 10:
            return 1.0
        
        row_sums = observed.sum(axis=1)
        col_sums = observed.sum(axis=0)
        total = observed.sum()
        
        if total == 0:
            return 1.0
        
        expected = np.outer(row_sums, col_sums) / total
        
        with np.errstate(divide='ignore', invalid='ignore'):
            chi2_stat = np.sum((observed - expected) ** 2 / expected)
        
        if np.isnan(chi2_stat) or np.isinf(chi2_stat):
            return 1.0
        
        p_value = 1 - chi2.cdf(chi2_stat, 1)
        
        return p_value
    
    def fisher_exact_test(self, sample_ref: int, sample_alt: int,
                          ancestor_ref: int, ancestor_alt: int) -> float:
        """Fisher's exact test for small sample sizes."""
        table = [[sample_ref, sample_alt],
                 [ancestor_ref, ancestor_alt]]
        
        try:
            _, p_value = fisher_exact(table, alternative='two-sided')
            return p_value
        except:
            return 1.0
    
    def calculate_strand_bias(self, ref_f: int, ref_r: int,
                              alt_f: int, alt_r: int) -> Tuple[float, float]:
        """
        Calculate strand bias using both Fisher's test and log-normal distribution.
        
        Returns: (fisher_score, log_normal_pvalue)
        """
        # Fisher's exact test (traditional)
        fisher_score = 0.0
        if ref_f + ref_r > 0 and alt_f + alt_r > 0:
            table = [[ref_f, ref_r], [alt_f, alt_r]]
            try:
                _, p_value = fisher_exact(table)
                if p_value > 0:
                    fisher_score = -10 * math.log10(p_value)
            except:
                pass
        
        # Log-normal p-value (muver method)
        log_normal_pvalue = 1.0
        if self.bias_analyzer and self.bias_result:
            log_normal_pvalue = self.bias_analyzer.calculate_p_value(alt_f, alt_r)
        
        return fisher_score, log_normal_pvalue
    
    def call_mutation(self, variant: dict, ancestor_name: str,
                      sample_name: str) -> Optional[Mutation]:
        """
        Determine if a variant represents a de novo mutation.
        
        Uses enhanced muver statistical models for improved accuracy.
        """
        ancestor_data = variant['samples'].get(ancestor_name, {})
        sample_data = variant['samples'].get(sample_name, {})
        
        if not ancestor_data or not sample_data:
            return None
        
        # Get genotypes
        ancestor_gt = ancestor_data.get('GT', './.')
        sample_gt = sample_data.get('GT', './.')
        
        # Skip if genotypes are identical
        if ancestor_gt == sample_gt:
            return None
        
        # Get strand-specific depths
        anc_depths = self.get_allele_depths(ancestor_data)
        smp_depths = self.get_allele_depths(sample_data)
        
        anc_ref = anc_depths[0] + anc_depths[1]
        anc_alt = anc_depths[2] + anc_depths[3]
        smp_ref = smp_depths[0] + smp_depths[1]
        smp_alt = smp_depths[2] + smp_depths[3]
        
        anc_total = anc_ref + anc_alt
        smp_total = smp_ref + smp_alt
        
        # Depth filter
        passes_depth = anc_total >= self.min_depth and smp_total >= self.min_depth
        if not passes_depth:
            return None
        
        # Depth distribution filter (muver)
        passes_depth_dist = True
        if self.depth_analyzer:
            chrom = variant['chromosome']
            pos = variant['position']
            passes_depth_dist = not self.depth_analyzer.is_position_filtered(chrom, pos)
        
        # Allele frequency filter
        smp_af = smp_alt / smp_total if smp_total > 0 else 0
        anc_af = anc_alt / anc_total if anc_total > 0 else 0
        
        # For haploid: sample should have high AF, ancestor should have low AF
        if smp_af < self.min_af or anc_af > 0.1:
            return None
        
        # Determine mutation type
        ref = variant['ref_allele']
        alt = variant['alt_alleles'][0]
        
        if len(ref) == len(alt) == 1:
            mut_type = 'SNP'
        elif len(ref) < len(alt):
            mut_type = 'INS'
        else:
            mut_type = 'DEL'
        
        # Basic statistical tests
        p_binom = self.binomial_test(smp_alt, smp_total, anc_alt, anc_total)
        p_chi = self.chi_square_test(smp_ref, smp_alt, anc_ref, anc_alt)
        p_fisher = self.fisher_exact_test(smp_ref, smp_alt, anc_ref, anc_alt)
        
        # Per-strand binomial tests (muver)
        smp_strand_counts = variant['strand_allele_counts'].get(sample_name, {}).get(alt, {})
        anc_strand_counts = variant['strand_allele_counts'].get(ancestor_name, {}).get(alt, {})
        p_binom_fwd, p_binom_rev = self.binomial_test_per_strand(smp_strand_counts, anc_strand_counts)
        
        # Composite significance (muver)
        composite_score, passes_composite = self.calculate_composite_significance(
            p_binom_fwd, p_binom_rev, p_chi
        )
        
        # Strand bias (enhanced)
        strand_bias_fisher, strand_bias_pvalue = self.calculate_strand_bias(
            smp_depths[0], smp_depths[1], smp_depths[2], smp_depths[3]
        )
        passes_strand = strand_bias_fisher < self.strand_bias_threshold * 10
        if self.bias_result:
            passes_strand = passes_strand and strand_bias_pvalue >= self.p_threshold
        
        # Repeat INDEL context
        in_repeat = False
        repeat_unit = ""
        repeat_tract_length = 0
        repeat_adjustment = 0.0
        
        if self.repeat_analyzer and 'intersected_repeat' in variant:
            repeat_info = variant['intersected_repeat']
            in_repeat = True
            repeat_unit = repeat_info['unit']
            repeat_tract_length = repeat_info['tract_length']
            
            if self.repeat_result and mut_type in ['INS', 'DEL']:
                event_type = 'insertion' if mut_type == 'INS' else 'deletion'
                repeat_adjustment = get_repeat_adjustment_value(
                    len(repeat_unit), repeat_tract_length,
                    event_type, self.repeat_result.fits
                )
        
        # Corrected depths
        corrected_smp_depth = float(smp_total)
        corrected_anc_depth = float(anc_total)
        if self.depth_corrector:
            chrom = variant['chromosome']
            pos = variant['position']
            corrected_smp_depth = self.depth_corrector.correct_depth(chrom, pos, smp_total)
            corrected_anc_depth = self.depth_corrector.correct_depth(chrom, pos, anc_total)
        
        # Create mutation object
        mutation = Mutation(
            chromosome=variant['chromosome'],
            position=variant['position'],
            ref_allele=ref,
            alt_allele=alt,
            sample_name=sample_name,
            mutation_type=mut_type,
            ancestor_genotype=ancestor_gt,
            sample_genotype=sample_gt,
            ancestor_ref_depth=anc_ref,
            ancestor_alt_depth=anc_alt,
            sample_ref_depth=smp_ref,
            sample_alt_depth=smp_alt,
            corrected_sample_depth=corrected_smp_depth,
            corrected_ancestor_depth=corrected_anc_depth,
            sample_ref_forward=smp_depths[0],
            sample_ref_reverse=smp_depths[1],
            sample_alt_forward=smp_depths[2],
            sample_alt_reverse=smp_depths[3],
            ancestor_ref_forward=anc_depths[0],
            ancestor_ref_reverse=anc_depths[1],
            ancestor_alt_forward=anc_depths[2],
            ancestor_alt_reverse=anc_depths[3],
            quality=variant['quality'],
            p_value_binomial=p_binom,
            p_value_chi_square=p_chi,
            p_value_fisher=p_fisher,
            p_value_binomial_forward=p_binom_fwd,
            p_value_binomial_reverse=p_binom_rev,
            composite_score=composite_score,
            strand_bias=strand_bias_fisher,
            strand_bias_pvalue=strand_bias_pvalue,
            in_repeat=in_repeat,
            repeat_unit=repeat_unit,
            repeat_tract_length=repeat_tract_length,
            repeat_adjustment=repeat_adjustment,
            passes_depth_filter=passes_depth,
            passes_strand_filter=passes_strand,
            passes_depth_distribution_filter=passes_depth_dist,
            passes_composite_test=passes_composite,
            is_significant=(p_binom < self.p_threshold and p_fisher < self.p_threshold)
        )
        
        return mutation
    
    def apply_fwer_correction(self, mutations: List[Mutation],
                              num_tests: int, genome_size: int = None) -> List[Mutation]:
        """
        Apply FWER correction using Bonferroni method.
        
        For MA analysis, we use the number of candidate mutations as num_tests,
        not genome_size (which would be too conservative).
        """
        # Use number of candidate mutations for Bonferroni
        # This is appropriate for MA experiments
        adjusted_threshold = self.fwer / num_tests if num_tests > 0 else self.fwer
        
        self.logger.info(f"FWER-adjusted threshold: {adjusted_threshold:.2e} (based on {num_tests} tests)")
        
        significant = []
        for mut in mutations:
            # Use composite test if enabled, otherwise use traditional tests
            if self.use_composite_significance:
                passes = (
                    mut.passes_composite_test and
                    mut.passes_strand_filter and
                    mut.passes_depth_distribution_filter
                )
            else:
                min_p = min(mut.p_value_binomial, mut.p_value_fisher)
                passes = (
                    min_p < adjusted_threshold and
                    mut.passes_strand_filter and
                    mut.passes_depth_distribution_filter
                )
            
            if passes:
                mut.is_significant = True
                significant.append(mut)
        
        return significant
    
    def detect_subclonal_mutations(self, variants: List[dict],
                                    ancestor_name: str,
                                    sample_name: str) -> List[Mutation]:
        """
        Detect subclonal mutations using muver method.
        """
        if not self.subclonal_detector:
            return []
        
        subclonal_mutations = []
        
        for variant in variants:
            # Call genotypes for control and sample
            control_result = self.subclonal_detector.call_genotype(variant, ancestor_name)
            sample_result = self.subclonal_detector.call_genotype(
                variant, sample_name, control_genotype=control_result.called_genotype
            )
            
            if sample_result.subclonal and sample_result.subclonal.is_subclonal:
                # Create mutation object for subclonal
                ref = variant['ref_allele']
                alt = sample_result.subclonal.subclonal_allele
                
                if len(ref) == 1 and len(alt) == 1:
                    mut_type = 'SNP'
                elif len(ref) < len(alt):
                    mut_type = 'INS'
                else:
                    mut_type = 'DEL'
                
                mutation = Mutation(
                    chromosome=variant['chromosome'],
                    position=variant['position'],
                    ref_allele=ref,
                    alt_allele=alt,
                    sample_name=sample_name,
                    mutation_type=mut_type,
                    is_subclonal=True,
                    subclonal_frequency=sample_result.subclonal.subclonal_frequency or 0.0,
                    subclonal_allele=sample_result.subclonal.subclonal_allele or "",
                    quality=variant['quality']
                )
                
                subclonal_mutations.append(mutation)
        
        return subclonal_mutations
    
    def write_mutations(self, mutations: List[Mutation], output_file: str,
                        write_all: bool = False):
        """Write mutations to TSV file with enhanced columns."""
        with open(output_file, 'w') as f:
            # Enhanced header with muver columns
            headers = [
                'chromosome', 'position', 'ref_allele', 'alt_allele',
                'sample_name', 'mutation_type',
                'ancestor_genotype', 'sample_genotype',
                'ancestor_ref_depth', 'ancestor_alt_depth',
                'sample_ref_depth', 'sample_alt_depth',
                'sample_alt_forward', 'sample_alt_reverse',
                'corrected_sample_depth', 'corrected_ancestor_depth',
                'quality',
                'p_value_binomial', 'p_value_chi_square', 'p_value_fisher',
                'p_value_binomial_forward', 'p_value_binomial_reverse',
                'composite_score',
                'strand_bias', 'strand_bias_pvalue',
                'in_repeat', 'repeat_unit', 'repeat_tract_length', 'repeat_adjustment',
                'is_subclonal', 'subclonal_frequency',
                'passes_composite_test', 'passes_depth_distribution_filter',
                'is_significant',
                'trinucleotide_context', 'gene', 'effect', 'impact'
            ]
            f.write('\t'.join(headers) + '\n')
            
            for mut in mutations:
                if not write_all and not mut.is_significant:
                    continue
                    
                values = [
                    mut.chromosome, str(mut.position), mut.ref_allele, mut.alt_allele,
                    mut.sample_name, mut.mutation_type,
                    mut.ancestor_genotype, mut.sample_genotype,
                    str(mut.ancestor_ref_depth), str(mut.ancestor_alt_depth),
                    str(mut.sample_ref_depth), str(mut.sample_alt_depth),
                    str(mut.sample_alt_forward), str(mut.sample_alt_reverse),
                    f'{mut.corrected_sample_depth:.2f}', f'{mut.corrected_ancestor_depth:.2f}',
                    f'{mut.quality:.2f}',
                    f'{mut.p_value_binomial:.2e}', f'{mut.p_value_chi_square:.2e}',
                    f'{mut.p_value_fisher:.2e}',
                    f'{mut.p_value_binomial_forward:.2e}', f'{mut.p_value_binomial_reverse:.2e}',
                    f'{mut.composite_score:.4f}',
                    f'{mut.strand_bias:.2f}', f'{mut.strand_bias_pvalue:.4f}',
                    str(mut.in_repeat), mut.repeat_unit, str(mut.repeat_tract_length),
                    f'{mut.repeat_adjustment:.4f}',
                    str(mut.is_subclonal), f'{mut.subclonal_frequency:.4f}',
                    str(mut.passes_composite_test), str(mut.passes_depth_distribution_filter),
                    str(mut.is_significant),
                    mut.trinucleotide_context, mut.gene, mut.effect, mut.impact
                ]
                f.write('\t'.join(values) + '\n')
    
    def write_vcf(self, mutations: List[Mutation], output_file: str,
                  reference: str):
        """Write mutations to VCF format with enhanced INFO fields."""
        with open(output_file, 'w') as f:
            # Header
            f.write('##fileformat=VCFv4.2\n')
            f.write(f'##reference={reference}\n')
            f.write('##INFO=<ID=MT,Number=1,Type=String,Description="Mutation type">\n')
            f.write('##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name">\n')
            f.write('##INFO=<ID=PVAL,Number=1,Type=Float,Description="Fisher p-value">\n')
            f.write('##INFO=<ID=COMP,Number=1,Type=Float,Description="Composite significance score">\n')
            f.write('##INFO=<ID=SB,Number=1,Type=Float,Description="Strand bias score">\n')
            f.write('##INFO=<ID=REP,Number=0,Type=Flag,Description="In repeat region">\n')
            f.write('##INFO=<ID=SUB,Number=0,Type=Flag,Description="Subclonal variant">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            
            for mut in mutations:
                if not mut.is_significant:
                    continue
                
                info_parts = [
                    f'MT={mut.mutation_type}',
                    f'SAMPLE={mut.sample_name}',
                    f'PVAL={mut.p_value_fisher:.2e}',
                    f'COMP={mut.composite_score:.4f}',
                    f'SB={mut.strand_bias:.2f}'
                ]
                
                if mut.in_repeat:
                    info_parts.append('REP')
                if mut.is_subclonal:
                    info_parts.append('SUB')
                
                info = ';'.join(info_parts)
                line = f'{mut.chromosome}\t{mut.position}\t.\t{mut.ref_allele}\t{mut.alt_allele}\t{mut.quality:.2f}\tPASS\t{info}\n'
                f.write(line)
    
    def run(self, samples: List, ancestor, vcf_file: str,
            bam_files: Dict[str, str] = None) -> List[Mutation]:
        """
        Run enhanced mutation calling on all samples.
        
        Parameters
        ----------
        samples : list
            List of Sample objects
        ancestor : Sample
            Ancestor sample
        vcf_file : str
            Path to joint-called VCF file
        bam_files : dict, optional
            Dictionary mapping sample names to BAM files for muver initialization
        """
        self.logger.info("Starting enhanced mutation calling with muver models")
        
        # Initialize muver models if BAM files provided
        if bam_files:
            reference = self.config['reference']['fasta']
            self.initialize_muver_models(bam_files, reference)
        
        ancestor_name = ancestor.sample_name
        ma_lines = [s for s in samples if s.sample_name != ancestor_name]
        
        # Parse VCF
        variants = self.parse_vcf(vcf_file)
        
        # Estimate genome size for FWER correction
        genome_size = self.config.get('statistics', {}).get('callable_sites')
        if not genome_size:
            # Try to get from reference
            reference = self.config['reference']['fasta']
            try:
                genome_size = sum(1 for _ in open(reference) if not _.startswith('>'))
            except:
                genome_size = 36_000_000  # Default for F. graminearum
        
        # Call mutations for each MA line
        all_mutations = []
        subclonal_mutations = []
        num_tests = len(variants) * len(ma_lines)
        
        for sample in ma_lines:
            self.logger.info(f"Calling mutations for {sample.sample_name}...")
            
            sample_mutations = []
            for variant in variants:
                mutation = self.call_mutation(variant, ancestor_name, sample.sample_name)
                if mutation:
                    sample_mutations.append(mutation)
            
            all_mutations.extend(sample_mutations)
            self.logger.info(f"  Found {len(sample_mutations)} candidate mutations")
            
            # Detect subclonal mutations
            if self.use_subclonal_detection:
                subcl = self.detect_subclonal_mutations(variants, ancestor_name, sample.sample_name)
                subclonal_mutations.extend(subcl)
                if subcl:
                    self.logger.info(f"  Found {len(subcl)} subclonal variants")
        
        # Apply FWER correction
        significant = self.apply_fwer_correction(all_mutations, num_tests, genome_size)
        
        self.logger.info(f"Total significant mutations: {len(significant)}")
        
        # Write outputs
        self.write_mutations(significant, str(self.output_dir / 'mutations.tsv'))
        self.write_mutations(all_mutations, str(self.output_dir / 'all_candidates.tsv'), write_all=True)
        self.write_vcf(significant, str(self.output_dir / 'mutations.vcf'),
                      self.config['reference']['fasta'])
        
        # Write subclonal mutations
        if subclonal_mutations:
            self.write_mutations(subclonal_mutations, 
                               str(self.output_dir / 'subclonal_mutations.tsv'),
                               write_all=True)
        
        return significant


# Keep backward compatibility
MutationCallingModule = EnhancedMutationCallingModule
