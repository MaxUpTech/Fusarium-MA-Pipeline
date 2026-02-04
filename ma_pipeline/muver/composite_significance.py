#!/usr/bin/env python3
"""
Composite Significance Scoring for Fusarium MA Pipeline

Implements composite significance scoring from the muver pipeline.
Combines multiple statistical tests (binomial, chi-square) into a
single significance score using Euclidean distance in log-p-value space.
"""

import logging
import math
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import binom, chi2


@dataclass
class CompositeSignificanceResult:
    """Results from composite significance calculation."""
    composite_score: float
    p_binomial_forward: float
    p_binomial_reverse: float
    p_chi_square: float
    is_significant: bool
    passes_individual_threshold: bool


class CompositeSignificanceCalculator:
    """
    Calculates composite significance scores for variants.
    
    The composite significance score combines binomial tests (per strand)
    and chi-square tests into a single metric. This provides more robust
    detection of true mutations while reducing false positives from
    sequencing artifacts.
    """
    
    def __init__(self, config: dict = None):
        self.config = config or {}
        self.logger = logging.getLogger('FusariumMA.CompositeSignificance')
        
        # Thresholds
        self.p_threshold = self.config.get('p_threshold', 0.01)
        self.individual_threshold = self.config.get('individual_threshold', 0.1)
        self.fwer = self.config.get('fwer', 0.01)
        
    def set_genome_adjusted_threshold(self, genome_size: int):
        """
        Calculate FWER-adjusted p-value threshold.
        
        Parameters
        ----------
        genome_size : int
            Size of genome in bp
        """
        self.p_threshold = 1 - ((1 - self.fwer) ** (1 / genome_size))
        self.logger.info(f"FWER-adjusted threshold: {self.p_threshold:.2e}")
    
    def binomial_test_per_strand(self, sample_counts: Dict[str, int],
                                  control_counts: Dict[str, int],
                                  allele: str) -> Tuple[float, float]:
        """
        Perform binomial tests for each strand.
        
        Parameters
        ----------
        sample_counts : dict
            Sample allele counts with keys 'forward' and 'reverse'
        control_counts : dict
            Control allele counts with keys 'forward' and 'reverse'
        allele : str
            Allele to test
            
        Returns
        -------
        tuple
            (p_forward, p_reverse) p-values for each strand
        """
        sample_total = sum(sample_counts.values())
        control_total = sum(control_counts.values())
        
        p_values = {}
        
        for strand in ['forward', 'reverse']:
            sample_value = sample_counts.get(strand, 0)
            control_value = max(control_counts.get(strand, 0), 1)
            
            if sample_value > 0 and sample_total > 0 and control_total > 0:
                # Expected frequency based on control
                expected_freq = float(control_value) / control_total
                
                # P(observing <= (sample_total - sample_value) | expected_freq)
                # This tests if sample has MORE of this allele than expected
                p_values[strand] = binom.cdf(
                    sample_total - int(sample_value),
                    sample_total,
                    1.0 - expected_freq
                )
            else:
                p_values[strand] = 1.0
        
        return p_values.get('forward', 1.0), p_values.get('reverse', 1.0)
    
    def chi_square_test(self, sample_counts: Dict[str, Dict[str, int]],
                        control_counts: Dict[str, Dict[str, int]]) -> float:
        """
        Perform chi-square test comparing allele distributions.
        
        Parameters
        ----------
        sample_counts : dict
            Sample counts: allele -> strand -> count
        control_counts : dict
            Control counts: allele -> strand -> count
            
        Returns
        -------
        float
            Chi-square p-value
        """
        sample_total = sum(
            c for allele_counts in sample_counts.values()
            for c in allele_counts.values()
        )
        control_total = sum(
            c for allele_counts in control_counts.values()
            for c in allele_counts.values()
        )
        
        if sample_total == 0 or control_total == 0:
            return 1.0
        
        chi_sum = 0.0
        df = 0
        
        for allele in set(sample_counts.keys()) | set(control_counts.keys()):
            for strand in ['forward', 'reverse']:
                s = float(sample_counts.get(allele, {}).get(strand, 0))
                c = float(control_counts.get(allele, {}).get(strand, 0))
                
                if s > 0 or c > 0:
                    s_sum = float(sample_total)
                    c_sum = float(control_total)
                    
                    chi_sum += ((s * math.sqrt(c_sum / s_sum) - 
                                 c * math.sqrt(s_sum / c_sum)) ** 2) / (s + c)
                    df += 1
        
        if df == 0:
            return 1.0
        
        return 1.0 - chi2.cdf(chi_sum, df)
    
    def calculate_composite_score(self, p_binom_forward: float,
                                   p_binom_reverse: float,
                                   p_chi: float) -> float:
        """
        Calculate composite significance score.
        
        The composite score is the Euclidean distance in -log10(p) space.
        
        Parameters
        ----------
        p_binom_forward : float
            P-value from forward strand binomial test
        p_binom_reverse : float
            P-value from reverse strand binomial test
        p_chi : float
            P-value from chi-square test
            
        Returns
        -------
        float
            Composite significance score
        """
        # Ensure minimum p-values to avoid log(0)
        min_p = sys.float_info.min
        p_binom_forward = max(p_binom_forward, min_p)
        p_binom_reverse = max(p_binom_reverse, min_p)
        p_chi = max(p_chi, min_p)
        
        # Calculate composite score as Euclidean distance in -log10 space
        composite = math.sqrt(
            (-math.log10(p_binom_forward)) ** 2 +
            (-math.log10(p_binom_reverse)) ** 2 +
            (-math.log10(p_chi)) ** 2
        )
        
        return composite
    
    def check_significance(self, p_binom_forward: float,
                            p_binom_reverse: float,
                            p_chi: float) -> CompositeSignificanceResult:
        """
        Check if a variant is significant based on composite score.
        
        For a variant to be called significant:
        1. Composite score must exceed -log10(p_threshold)
        2. All individual p-values must be < individual_threshold (default 0.1)
        
        Parameters
        ----------
        p_binom_forward : float
            P-value from forward strand binomial test
        p_binom_reverse : float
            P-value from reverse strand binomial test
        p_chi : float
            P-value from chi-square test
            
        Returns
        -------
        CompositeSignificanceResult
            Significance results
        """
        composite = self.calculate_composite_score(
            p_binom_forward, p_binom_reverse, p_chi
        )
        
        log_threshold = -math.log10(self.p_threshold)
        
        # Check individual thresholds
        passes_individual = (
            p_binom_forward < self.individual_threshold and
            p_binom_reverse < self.individual_threshold and
            p_chi < self.individual_threshold
        )
        
        # Overall significance
        is_significant = (composite >= log_threshold) and passes_individual
        
        return CompositeSignificanceResult(
            composite_score=composite,
            p_binomial_forward=p_binom_forward,
            p_binomial_reverse=p_binom_reverse,
            p_chi_square=p_chi,
            is_significant=is_significant,
            passes_individual_threshold=passes_individual
        )
    
    def calculate_for_variant(self, variant: dict,
                              sample_name: str,
                              control_name: str) -> CompositeSignificanceResult:
        """
        Calculate composite significance for a variant.
        
        Parameters
        ----------
        variant : dict
            Variant data with strand allele counts
        sample_name : str
            Name of sample
        control_name : str
            Name of control
            
        Returns
        -------
        CompositeSignificanceResult
            Significance results
        """
        # Get strand allele counts
        sample_counts = variant.get('strand_allele_counts', {}).get(sample_name, {})
        control_counts = variant.get('strand_allele_counts', {}).get(control_name, {})
        
        # Get alleles
        alleles = variant.get('alleles', [])
        if not alleles:
            return CompositeSignificanceResult(
                composite_score=0.0,
                p_binomial_forward=1.0,
                p_binomial_reverse=1.0,
                p_chi_square=1.0,
                is_significant=False,
                passes_individual_threshold=False
            )
        
        # Calculate best composite score across alleles
        best_result = None
        
        for allele in alleles:
            allele_sample_counts = sample_counts.get(allele, {'forward': 0, 'reverse': 0})
            allele_control_counts = control_counts.get(allele, {'forward': 0, 'reverse': 0})
            
            p_forward, p_reverse = self.binomial_test_per_strand(
                allele_sample_counts, allele_control_counts, allele
            )
            
            p_chi = self.chi_square_test(sample_counts, control_counts)
            
            result = self.check_significance(p_forward, p_reverse, p_chi)
            
            if best_result is None or result.composite_score > best_result.composite_score:
                best_result = result
        
        return best_result


def calculate_composite_significance(variant: dict,
                                      sample_name: str,
                                      control_name: str,
                                      config: dict = None) -> CompositeSignificanceResult:
    """
    Convenience function to calculate composite significance.
    
    Parameters
    ----------
    variant : dict
        Variant data
    sample_name : str
        Sample name
    control_name : str
        Control name
    config : dict, optional
        Configuration
        
    Returns
    -------
    CompositeSignificanceResult
        Significance results
    """
    calculator = CompositeSignificanceCalculator(config)
    return calculator.calculate_for_variant(variant, sample_name, control_name)


def apply_composite_filter(variants: List[dict],
                           sample_name: str,
                           control_name: str,
                           p_threshold: float = 0.01,
                           individual_threshold: float = 0.1) -> List[dict]:
    """
    Apply composite significance filter to variants.
    
    Parameters
    ----------
    variants : list
        List of variant dictionaries
    sample_name : str
        Sample name
    control_name : str
        Control name
    p_threshold : float
        P-value threshold for composite score
    individual_threshold : float
        Threshold for individual tests
        
    Returns
    -------
    list
        Filtered variants that pass significance test
    """
    config = {
        'p_threshold': p_threshold,
        'individual_threshold': individual_threshold
    }
    
    calculator = CompositeSignificanceCalculator(config)
    
    filtered = []
    for var in variants:
        result = calculator.calculate_for_variant(var, sample_name, control_name)
        
        if result.is_significant:
            var['composite_score'] = result.composite_score
            var['composite_significant'] = True
            filtered.append(var)
    
    return filtered
