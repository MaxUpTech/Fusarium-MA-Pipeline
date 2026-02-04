#!/usr/bin/env python3
"""
Strand Bias Distribution Model for Fusarium MA Pipeline

Implements log-normal strand bias characterization from the muver pipeline.
Calculates the distribution of strand bias across the genome and uses it
to filter variants with excessive strand imbalance.
"""

import logging
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit

from .fitting import gaussian


@dataclass
class BiasDistributionResult:
    """Results from strand bias distribution analysis."""
    mu: float  # Mean of log-ratios
    sigma: float  # Standard deviation of log-ratios
    log_ratios: List[float] = field(default_factory=list)
    histogram: Optional[np.ndarray] = None
    bin_edges: Optional[np.ndarray] = None


class BiasDistributionAnalyzer:
    """
    Analyzes strand bias distribution across the genome.
    
    The strand bias is characterized by fitting a log-normal distribution
    to the log-ratios of forward/reverse strand counts for each allele.
    This distribution is then used to identify variants with abnormal
    strand bias that may represent sequencing artifacts.
    """
    
    def __init__(self, config: dict = None):
        self.config = config or {}
        self.logger = logging.getLogger('FusariumMA.BiasDistribution')
        self.result: Optional[BiasDistributionResult] = None
        
    def calculate_from_bam(self, bam_file: str, reference: str,
                           output_file: Optional[str] = None) -> BiasDistributionResult:
        """
        Calculate strand bias distribution from a BAM file.
        
        Parameters
        ----------
        bam_file : str
            Path to BAM file
        reference : str
            Path to reference FASTA
        output_file : str, optional
            Path to output file for distribution data
            
        Returns
        -------
        BiasDistributionResult
            Fitted distribution parameters and data
        """
        import pysam
        
        self.logger.info(f"Calculating strand bias distribution from {bam_file}")
        
        log_ratios = []
        
        # Open BAM file and iterate through pileup
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for pileup_column in bam.pileup(min_base_quality=20, min_mapping_quality=30):
                plus_tally = defaultdict(int)
                minus_tally = defaultdict(int)
                
                for pileup_read in pileup_column.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    
                    base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                    
                    if pileup_read.alignment.is_reverse:
                        minus_tally[base.upper()] += 1
                    else:
                        plus_tally[base.upper()] += 1
                
                # Calculate log ratios for alleles present on both strands
                all_alleles = set(plus_tally.keys()) | set(minus_tally.keys())
                for allele in all_alleles:
                    if allele in plus_tally and allele in minus_tally:
                        if plus_tally[allele] > 0 and minus_tally[allele] > 0:
                            ratio = float(plus_tally[allele]) / minus_tally[allele]
                            log_ratios.append(math.log(ratio))
        
        return self._fit_distribution(log_ratios, output_file)
    
    def calculate_from_mpileup(self, mpileup_file: str, 
                                output_file: Optional[str] = None) -> BiasDistributionResult:
        """
        Calculate strand bias distribution from an mpileup file.
        
        Parameters
        ----------
        mpileup_file : str
            Path to mpileup file
        output_file : str, optional
            Path to output file
            
        Returns
        -------
        BiasDistributionResult
            Fitted distribution parameters
        """
        self.logger.info(f"Calculating strand bias distribution from {mpileup_file}")
        
        log_ratios = []
        
        with open(mpileup_file, 'r') as f:
            for line in f:
                line_split = line.strip().split()
                if len(line_split) < 5:
                    continue
                
                reference, coverage = line_split[2:4]
                bases = line_split[4] if int(coverage) > 0 else ''
                
                plus_tally = defaultdict(int)
                minus_tally = defaultdict(int)
                present_alleles = set()
                
                i = 0
                while i < len(bases):
                    if bases[i] == '.':
                        present_alleles.add(reference)
                        plus_tally[reference] += 1
                    elif bases[i] == ',':
                        present_alleles.add(reference)
                        minus_tally[reference] += 1
                    elif re.match('[ACGT]', bases[i]):
                        present_alleles.add(bases[i])
                        plus_tally[bases[i]] += 1
                    elif re.match('[acgt]', bases[i]):
                        present_alleles.add(bases[i].upper())
                        minus_tally[bases[i].upper()] += 1
                    elif re.match('[+-]', bases[i]):
                        # Handle indels
                        indel_type = bases[i]
                        i += 1
                        if i < len(bases) and bases[i].isdigit():
                            indel_length = int(bases[i])
                            i += 1
                            indel_seq = bases[i:i+indel_length].upper()
                            indel = indel_type + indel_seq
                            present_alleles.add(indel)
                            if indel_seq.isupper():
                                plus_tally[indel] += 1
                            else:
                                minus_tally[indel] += 1
                            i += indel_length - 1
                    elif bases[i] == '^':
                        i += 1  # Skip quality character
                    i += 1
                
                # Calculate log ratios
                for allele in present_alleles:
                    if allele in plus_tally and allele in minus_tally:
                        if plus_tally[allele] > 0 and minus_tally[allele] > 0:
                            ratio = float(plus_tally[allele]) / minus_tally[allele]
                            log_ratios.append(math.log(ratio))
        
        return self._fit_distribution(log_ratios, output_file)
    
    def _fit_distribution(self, log_ratios: List[float],
                          output_file: Optional[str] = None) -> BiasDistributionResult:
        """
        Fit a Gaussian distribution to log-ratios.
        
        Parameters
        ----------
        log_ratios : list
            List of log-ratios
        output_file : str, optional
            Path to output file
            
        Returns
        -------
        BiasDistributionResult
            Fitted distribution parameters
        """
        if not log_ratios:
            self.logger.warning("No log ratios calculated, using default parameters")
            return BiasDistributionResult(mu=0.0, sigma=1.0)
        
        # Initial fit using scipy.stats
        p0_mu, p0_sigma = norm.fit(log_ratios)
        if p0_sigma == 0:
            p0_sigma = 0.01
        
        # Create histogram
        bins = [float(x)/10 for x in range(-50, 51)]
        hist, bin_edges = np.histogram(log_ratios, bins=bins, density=True)
        
        # Fit Gaussian using curve_fit
        try:
            popt, pcov = curve_fit(gaussian, bin_edges[:-1], hist,
                                   p0=[p0_mu, p0_sigma], maxfev=100000)
            mu, sigma = popt
            sigma = abs(sigma)
        except:
            mu, sigma = p0_mu, abs(p0_sigma)
        
        result = BiasDistributionResult(
            mu=mu,
            sigma=sigma,
            log_ratios=log_ratios,
            histogram=hist,
            bin_edges=bin_edges
        )
        
        self.result = result
        
        # Write output file if requested
        if output_file:
            self._write_output(result, output_file)
        
        self.logger.info(f"Strand bias distribution: mu={mu:.4f}, sigma={sigma:.4f}")
        
        return result
    
    def _write_output(self, result: BiasDistributionResult, output_file: str):
        """Write distribution results to file."""
        with open(output_file, 'w') as f:
            f.write(f'Average log ratio: {result.mu}\n')
            f.write(f'Standard deviation of log ratios: {result.sigma}\n\n')
            f.write('Bias distribution:\n\n')
            f.write('\t'.join(['Strand log ratio', 'Frequency', 'Fit value']) + '\n')
            
            if result.histogram is not None and result.bin_edges is not None:
                for hist_val, bin_val in zip(result.histogram, result.bin_edges[:-1]):
                    fit_val = norm.pdf(bin_val, result.mu, result.sigma)
                    f.write(f'{round(bin_val, 1)}\t{hist_val}\t{fit_val}\n')
    
    def calculate_p_value(self, forward_count: int, reverse_count: int) -> float:
        """
        Calculate p-value for strand bias given the fitted distribution.
        
        Parameters
        ----------
        forward_count : int
            Count on forward strand
        reverse_count : int
            Count on reverse strand
            
        Returns
        -------
        float
            P-value for the observed strand bias
        """
        if self.result is None:
            raise ValueError("Must fit distribution before calculating p-values")
        
        if forward_count == 0 or reverse_count == 0:
            return 0.0
        
        log_ratio = math.log(float(forward_count) / reverse_count)
        # Two-tailed p-value using log-normal distribution
        abs_log_ratio = abs(log_ratio - self.result.mu)
        p_value = 2 * 0.5 * (1 + math.erf(-abs_log_ratio / (math.sqrt(2) * self.result.sigma)))
        
        return p_value
    
    def passes_filter(self, forward_count: int, reverse_count: int,
                      p_threshold: float = 0.01) -> bool:
        """
        Check if a variant passes the strand bias filter.
        
        Parameters
        ----------
        forward_count : int
            Count on forward strand
        reverse_count : int
            Count on reverse strand
        p_threshold : float
            P-value threshold
            
        Returns
        -------
        bool
            True if variant passes filter
        """
        p_value = self.calculate_p_value(forward_count, reverse_count)
        return p_value >= p_threshold


def calculate_strand_bias_distribution(bam_file: str, reference: str,
                                       output_file: Optional[str] = None,
                                       config: dict = None) -> BiasDistributionResult:
    """
    Convenience function to calculate strand bias distribution.
    
    Parameters
    ----------
    bam_file : str
        Path to BAM file
    reference : str
        Path to reference FASTA
    output_file : str, optional
        Path to output file
    config : dict, optional
        Configuration dictionary
        
    Returns
    -------
    BiasDistributionResult
        Fitted distribution parameters
    """
    analyzer = BiasDistributionAnalyzer(config)
    return analyzer.calculate_from_bam(bam_file, reference, output_file)


def apply_strand_bias_filter(variants: List[dict], 
                             bias_result: BiasDistributionResult,
                             p_threshold: float = 0.01) -> List[dict]:
    """
    Apply strand bias filter to a list of variants.
    
    Parameters
    ----------
    variants : list
        List of variant dictionaries with 'forward_count' and 'reverse_count'
    bias_result : BiasDistributionResult
        Fitted distribution parameters
    p_threshold : float
        P-value threshold for filtering
        
    Returns
    -------
    list
        Filtered variants that pass strand bias filter
    """
    analyzer = BiasDistributionAnalyzer()
    analyzer.result = bias_result
    
    filtered = []
    for var in variants:
        fwd = var.get('forward_count', 0)
        rev = var.get('reverse_count', 0)
        if analyzer.passes_filter(fwd, rev, p_threshold):
            var['strand_bias_pvalue'] = analyzer.calculate_p_value(fwd, rev)
            filtered.append(var)
    
    return filtered
