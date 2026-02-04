#!/usr/bin/env python3
"""
Subclonal Variant Detection for Fusarium MA Pipeline

Implements subclonal variant detection from the muver pipeline.
Detects variants present at subclonal frequencies (e.g., 0.5, 0.25, 0.125)
using maximum likelihood genotype calling.
"""

import logging
import math
import copy
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
from scipy.stats import binom

from .repeat_indels import RepeatIndelAnalyzer, get_repeat_adjustment_value


@dataclass
class SubclonalCall:
    """Result of subclonal variant detection."""
    genotype: Tuple[str, ...]  # Called genotype
    subclonal_genotype: Optional[Tuple[str, ...]] = None  # Subclonal genotype if present
    subclonal_allele: Optional[str] = None  # Subclonal allele
    subclonal_frequency: Optional[float] = None  # Estimated frequency
    log_ratio_sum: float = float('inf')  # Quality metric
    is_subclonal: bool = False
    strand_bias_pvalue: Optional[float] = None
    binomial_pvalue: Optional[float] = None


@dataclass
class GenotypeCallResult:
    """Result of genotype calling for a variant."""
    sample_name: str
    called_genotype: Tuple[str, ...]
    subclonal: Optional[SubclonalCall] = None
    expected_frequencies: Dict[str, float] = field(default_factory=dict)
    observed_frequencies: Dict[str, float] = field(default_factory=dict)


class SubclonalDetector:
    """
    Detects subclonal variants using maximum likelihood genotype calling.
    
    Subclonal variants are mutations present in only a fraction of cells,
    which can occur due to mutations arising during culture or sequencing
    artifacts. This class identifies such variants and provides quality
    metrics to distinguish true subclonals from artifacts.
    """
    
    def __init__(self, config: dict = None):
        self.config = config or {}
        self.logger = logging.getLogger('FusariumMA.SubclonalDetection')
        
        # Parameters
        self.ploidy = self.config.get('ploidy', 1)  # Haploid for Fusarium
        self.subclonal_frequencies = self.config.get(
            'subclonal_frequencies', [0.5, 0.25, 0.125]
        )
        self.p_threshold = self.config.get('p_threshold', 0.01)
        
        # For strand bias testing
        self.strand_bias_std = self.config.get('strand_bias_std', 1.0)
        
        # Repeat indel analyzer (optional)
        self.repeat_analyzer: Optional[RepeatIndelAnalyzer] = None
        self.repeat_fits: Optional[Dict] = None
    
    def set_strand_bias_std(self, std: float):
        """Set strand bias standard deviation from distribution fitting."""
        self.strand_bias_std = std
    
    def set_repeat_fits(self, fits: Dict):
        """Set repeat indel fits for frequency adjustment."""
        self.repeat_fits = fits
    
    def get_possible_genotypes(self, alleles: List[str]) -> List[Tuple[str, ...]]:
        """
        Generate all possible genotypes for given alleles and ploidy.
        
        Parameters
        ----------
        alleles : list
            List of allele strings
            
        Returns
        -------
        list
            List of possible genotype tuples
        """
        possible = []
        genotypes = [[]]
        
        for _ in range(self.ploidy):
            temp = []
            for allele in alleles:
                for geno in genotypes:
                    temp.append(geno + [allele])
            genotypes = temp
        
        for geno in genotypes:
            sorted_geno = tuple(sorted(geno, key=lambda x: alleles.index(x)))
            if sorted_geno not in possible:
                possible.append(sorted_geno)
        
        return possible
    
    def calculate_expected_frequencies(self, genotype: Tuple[str, ...],
                                        alleles: List[str],
                                        subclonal: Optional[Tuple] = None) -> Dict[str, float]:
        """
        Calculate expected allele frequencies for a genotype.
        
        Parameters
        ----------
        genotype : tuple
            Genotype tuple
        alleles : list
            List of all alleles
        subclonal : tuple, optional
            (subclonal_genotype, subclonal_allele, subclonal_frequency)
            
        Returns
        -------
        dict
            Allele -> expected frequency
        """
        frequencies = {}
        
        for allele in alleles:
            freq = float(genotype.count(allele)) / len(genotype)
            
            if subclonal is not None:
                sub_geno, sub_allele, sub_freq = subclonal
                sub_geno_freq = float(sub_geno.count(allele)) / len(sub_geno)
                freq = (1 - sub_freq) * freq + sub_freq * sub_geno_freq
            
            frequencies[allele] = freq
        
        return frequencies
    
    def calculate_observed_frequencies(self, strand_counts: Dict[str, Dict[str, int]],
                                        alleles: List[str]) -> Dict[str, float]:
        """
        Calculate observed allele frequencies from strand counts.
        
        Parameters
        ----------
        strand_counts : dict
            Allele -> strand -> count
        alleles : list
            List of alleles
            
        Returns
        -------
        dict
            Allele -> observed frequency
        """
        total = 0
        allele_totals = {}
        
        for allele in alleles:
            counts = strand_counts.get(allele, {'forward': 0, 'reverse': 0})
            allele_total = max(counts.get('forward', 0), 1) + max(counts.get('reverse', 0), 1)
            allele_totals[allele] = allele_total
            total += allele_total - 2  # Subtract the pseudocounts
        
        # Recalculate total including pseudocounts properly
        total = sum(allele_totals.values())
        
        frequencies = {}
        for allele in alleles:
            frequencies[allele] = float(allele_totals[allele]) / total
        
        return frequencies
    
    def call_genotype(self, variant: dict, sample_name: str,
                       control_genotype: Optional[Tuple[str, ...]] = None) -> GenotypeCallResult:
        """
        Call genotype for a variant using maximum likelihood.
        
        Parameters
        ----------
        variant : dict
            Variant data with strand allele counts
        sample_name : str
            Sample to genotype
        control_genotype : tuple, optional
            Control genotype for comparison
            
        Returns
        -------
        GenotypeCallResult
            Genotype call result
        """
        alleles = variant.get('alleles', [])
        strand_counts = variant.get('strand_allele_counts', {}).get(sample_name, {})
        
        if not alleles or not strand_counts:
            return GenotypeCallResult(
                sample_name=sample_name,
                called_genotype=tuple()
            )
        
        # Calculate observed frequencies
        observed_freq = self.calculate_observed_frequencies(strand_counts, alleles)
        
        # Generate possible genotypes
        possible_genotypes = self.get_possible_genotypes(alleles)
        
        # Build expected frequency dict for all genotypes and subclonals
        expected_frequencies = {}
        
        # Without subclonal
        for genotype in possible_genotypes:
            key = (genotype, None)
            expected_frequencies[key] = self.calculate_expected_frequencies(
                genotype, alleles, subclonal=None
            )
        
        # With subclonals
        for subclonal_allele in alleles:
            for genotype in possible_genotypes:
                for i, allele in enumerate(genotype):
                    if allele != subclonal_allele:
                        sub_geno = list(genotype)
                        sub_geno[i] = subclonal_allele
                        sub_geno = tuple(sorted(sub_geno, key=lambda x: alleles.index(x)))
                        
                        for sub_freq in self.subclonal_frequencies:
                            subclonal = (sub_geno, subclonal_allele, sub_freq)
                            key = (genotype, subclonal)
                            expected_frequencies[key] = self.calculate_expected_frequencies(
                                genotype, alleles, subclonal=subclonal
                            )
        
        # Apply repeat indel correction if available
        if self.repeat_fits and variant.get('intersected_repeat'):
            self._apply_repeat_correction(expected_frequencies, variant, alleles)
        
        # Find best genotype by minimizing log-ratio sum
        min_log_ratio = float('inf')
        best_key = None
        
        ref_allele = variant.get('ref_allele', alleles[0] if alleles else '')
        
        for key, exp_freq in expected_frequencies.items():
            genotype, subclonal = key
            log_ratio_sum = 0
            
            for allele in alleles:
                obs = observed_freq.get(allele, 1e-10)
                exp = max(exp_freq.get(allele, 1e-10), 1e-10)
                log_ratio_sum += abs(math.log(obs / exp))
            
            # Prefer genotypes sharing alleles with control
            shared_count = 0
            if control_genotype:
                for allele in genotype:
                    if allele in control_genotype:
                        shared_count += 1
            
            # Slight preference for simpler models (no subclonal)
            if subclonal is None:
                log_ratio_sum -= 0.001 * shared_count
            
            if log_ratio_sum < min_log_ratio - 1e-6:
                min_log_ratio = log_ratio_sum
                best_key = key
        
        if best_key is None:
            return GenotypeCallResult(
                sample_name=sample_name,
                called_genotype=tuple()
            )
        
        best_genotype, best_subclonal = best_key
        
        # Create subclonal call if applicable
        subclonal_call = None
        if best_subclonal is not None:
            sub_geno, sub_allele, sub_freq = best_subclonal
            
            # Calculate strand bias p-value for subclonal
            sub_counts = strand_counts.get(sub_allele, {'forward': 0, 'reverse': 0})
            strand_bias_pval = self._calculate_strand_bias_pvalue(
                sub_counts.get('forward', 0),
                sub_counts.get('reverse', 0)
            )
            
            # Calculate binomial p-value
            binom_pval = self._calculate_subclonal_binomial_pvalue(
                strand_counts, alleles, best_genotype, sub_allele, expected_frequencies[best_key]
            )
            
            subclonal_call = SubclonalCall(
                genotype=best_genotype,
                subclonal_genotype=sub_geno,
                subclonal_allele=sub_allele,
                subclonal_frequency=sub_freq,
                log_ratio_sum=min_log_ratio,
                is_subclonal=True,
                strand_bias_pvalue=strand_bias_pval,
                binomial_pvalue=binom_pval
            )
        
        return GenotypeCallResult(
            sample_name=sample_name,
            called_genotype=best_genotype,
            subclonal=subclonal_call,
            expected_frequencies=expected_frequencies.get(best_key, {}),
            observed_frequencies=observed_freq
        )
    
    def _apply_repeat_correction(self, expected_frequencies: Dict,
                                  variant: dict, alleles: List[str]):
        """Apply repeat indel correction to expected frequencies."""
        repeat = variant.get('intersected_repeat')
        if not repeat or not self.repeat_fits:
            return
        
        unit_length = len(repeat.get('unit', ''))
        tract_length = repeat.get('tract_length', 0)
        
        if unit_length < 1 or unit_length > 4:
            return
        
        for key, frequencies in expected_frequencies.items():
            adjustments = {allele: 0.0 for allele in alleles}
            
            for allele, freq in frequencies.items():
                # Calculate adjustment for insertion
                ins_adj = get_repeat_adjustment_value(
                    unit_length, tract_length, 'insertion', self.repeat_fits
                ) * freq
                
                # Calculate adjustment for deletion
                del_adj = get_repeat_adjustment_value(
                    unit_length, tract_length, 'deletion', self.repeat_fits
                ) * freq
                
                # Apply adjustments (simplified - full implementation would
                # check for neighboring alleles)
                adjustments[allele] -= (ins_adj + del_adj)
            
            for allele, adj in adjustments.items():
                frequencies[allele] = max(frequencies[allele] + adj, 0.01)
    
    def _calculate_strand_bias_pvalue(self, forward: int, reverse: int) -> float:
        """Calculate strand bias p-value using log-normal distribution."""
        if forward == 0 or reverse == 0:
            return 0.0
        
        log_ratio = -abs(math.log(float(forward) / reverse))
        pvalue = 2 * 0.5 * (1 + math.erf(log_ratio / (math.sqrt(2) * self.strand_bias_std)))
        
        return pvalue
    
    def _calculate_subclonal_binomial_pvalue(self, strand_counts: Dict,
                                              alleles: List[str],
                                              genotype: Tuple[str, ...],
                                              subclonal_allele: str,
                                              expected_frequencies: Dict[str, float]) -> float:
        """Calculate binomial p-value for subclonal allele."""
        sub_counts = strand_counts.get(subclonal_allele, {'forward': 0, 'reverse': 0})
        sub_total = sub_counts.get('forward', 0) + sub_counts.get('reverse', 0)
        
        total = sum(
            strand_counts.get(a, {'forward': 0, 'reverse': 0}).get('forward', 0) +
            strand_counts.get(a, {'forward': 0, 'reverse': 0}).get('reverse', 0)
            for a in alleles
        )
        
        if total == 0:
            return 1.0
        
        # Expected frequency without subclonal
        exp_freq = float(genotype.count(subclonal_allele)) / len(genotype)
        
        # Binomial test
        pvalue = binom.cdf(total - sub_total, total, 1.0 - exp_freq)
        
        return pvalue
    
    def is_valid_subclonal(self, subclonal_call: SubclonalCall,
                           is_filtered_site: bool = False,
                           is_excluded: bool = False) -> bool:
        """
        Check if a subclonal call is valid.
        
        Parameters
        ----------
        subclonal_call : SubclonalCall
            Subclonal call to validate
        is_filtered_site : bool
            Whether position is in a filtered region
        is_excluded : bool
            Whether position is excluded
            
        Returns
        -------
        bool
            True if subclonal is valid
        """
        if not subclonal_call.is_subclonal:
            return False
        
        if is_filtered_site or is_excluded:
            return False
        
        # Check strand bias
        if (subclonal_call.strand_bias_pvalue is not None and
            subclonal_call.strand_bias_pvalue < self.p_threshold):
            return False
        
        # Check binomial test
        if (subclonal_call.binomial_pvalue is not None and
            subclonal_call.binomial_pvalue >= self.p_threshold):
            return False
        
        return True


def detect_subclonal_variants(variants: List[dict],
                              sample_name: str,
                              control_name: str,
                              config: dict = None) -> List[GenotypeCallResult]:
    """
    Detect subclonal variants in a list of variants.
    
    Parameters
    ----------
    variants : list
        List of variant dictionaries
    sample_name : str
        Sample name
    control_name : str
        Control name
    config : dict, optional
        Configuration
        
    Returns
    -------
    list
        List of GenotypeCallResult for variants with subclonal calls
    """
    detector = SubclonalDetector(config)
    results = []
    
    # First call control genotypes
    control_genotypes = {}
    for var in variants:
        var_id = (var.get('chromosome'), var.get('position'))
        control_result = detector.call_genotype(var, control_name)
        control_genotypes[var_id] = control_result.called_genotype
    
    # Then call sample genotypes with control comparison
    for var in variants:
        var_id = (var.get('chromosome'), var.get('position'))
        control_geno = control_genotypes.get(var_id)
        
        result = detector.call_genotype(var, sample_name, control_genotype=control_geno)
        
        if result.subclonal is not None and result.subclonal.is_subclonal:
            results.append(result)
    
    return results
