#!/usr/bin/env python3
"""
Depth Distribution Model for Fusarium MA Pipeline

Implements normal distribution-based depth filtering from the muver pipeline.
Fits depth values to a normal distribution and filters regions with abnormal
read depths that may indicate copy number variations or sequencing artifacts.
"""

import logging
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit

from .fitting import gaussian


@dataclass
class DepthDistributionResult:
    """Results from depth distribution analysis."""
    mu: float  # Mean depth per copy
    sigma: float  # Standard deviation
    depths: List[float] = field(default_factory=list)
    histogram: Optional[np.ndarray] = None
    bin_edges: Optional[np.ndarray] = None
    filtered_regions: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)


class DepthDistributionAnalyzer:
    """
    Analyzes read depth distribution across the genome.
    
    The depth distribution is characterized by fitting a normal distribution
    to per-position depths (normalized by ploidy). Regions with depths outside
    the expected distribution are flagged for filtering.
    """
    
    def __init__(self, config: dict = None):
        self.config = config or {}
        self.logger = logging.getLogger('FusariumMA.DepthDistribution')
        self.result: Optional[DepthDistributionResult] = None
        
        # Parameters
        self.ploidy = self.config.get('ploidy', 1)  # Fusarium is haploid
        self.p_threshold = self.config.get('p_threshold', 0.0001)
        self.window_size = self.config.get('window_size', 51)
        self.merge_window = self.config.get('merge_window', 1000)
    
    def calculate_from_bedgraph(self, bedgraph_file: str,
                                 output_file: Optional[str] = None,
                                 cnv_bedgraph: Optional[str] = None) -> DepthDistributionResult:
        """
        Calculate depth distribution from a bedGraph file.
        
        Parameters
        ----------
        bedgraph_file : str
            Path to bedGraph file with coverage data
        output_file : str, optional
            Path to output file for distribution data
        cnv_bedgraph : str, optional
            Path to CNV regions bedGraph (copy number per position)
            
        Returns
        -------
        DepthDistributionResult
            Fitted distribution parameters
        """
        self.logger.info(f"Calculating depth distribution from {bedgraph_file}")
        
        # Read CNV regions if provided
        cnv_regions = {}
        if cnv_bedgraph:
            cnv_regions = self._read_cnv_bedgraph(cnv_bedgraph)
        
        depths = []
        
        with open(bedgraph_file, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) < 4:
                    continue
                
                chromosome, start, end, coverage = fields[:4]
                start, end = int(start), int(end)
                coverage = float(coverage)
                
                for pos in range(start, end):
                    if (chromosome, pos) in cnv_regions:
                        copy_number = cnv_regions[(chromosome, pos)]
                    else:
                        copy_number = self.ploidy
                    
                    if copy_number > 0:
                        depths.append(int(coverage / copy_number))
        
        return self._fit_distribution(depths, output_file)
    
    def calculate_from_bam(self, bam_file: str, reference: str,
                           output_file: Optional[str] = None) -> DepthDistributionResult:
        """
        Calculate depth distribution from a BAM file.
        
        Parameters
        ----------
        bam_file : str
            Path to BAM file
        reference : str
            Path to reference FASTA
        output_file : str, optional
            Path to output file
            
        Returns
        -------
        DepthDistributionResult
            Fitted distribution parameters
        """
        import pysam
        
        self.logger.info(f"Calculating depth distribution from {bam_file}")
        
        depths = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for pileup_column in bam.pileup(min_base_quality=20, min_mapping_quality=30):
                depth = pileup_column.nsegments / self.ploidy
                depths.append(int(depth))
        
        return self._fit_distribution(depths, output_file)
    
    def _read_cnv_bedgraph(self, cnv_file: str) -> Dict[Tuple[str, int], int]:
        """Read CNV regions from bedGraph file."""
        cnv_regions = {}
        with open(cnv_file, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 4:
                    chrom, start, end, cn = fields[:4]
                    for pos in range(int(start), int(end)):
                        cnv_regions[(chrom, pos)] = int(float(cn))
        return cnv_regions
    
    def _fit_distribution(self, depths: List[float],
                          output_file: Optional[str] = None) -> DepthDistributionResult:
        """
        Fit a normal distribution to depth values.
        
        Parameters
        ----------
        depths : list
            List of depth values
        output_file : str, optional
            Path to output file
            
        Returns
        -------
        DepthDistributionResult
            Fitted distribution parameters
        """
        if not depths:
            self.logger.warning("No depth values, using default parameters")
            return DepthDistributionResult(mu=30.0, sigma=10.0)
        
        # Initial fit
        p0_mu, p0_sigma = norm.fit(depths)
        
        # Create histogram
        depth_max = max(depths)
        hist, bin_edges = np.histogram(depths, bins=range(1, depth_max + 2), density=True)
        
        # Fit Gaussian
        try:
            popt, pcov = curve_fit(gaussian, bin_edges[:-1], hist,
                                   p0=[p0_mu, p0_sigma])
            mu, sigma = popt
            sigma = abs(sigma)
        except:
            mu, sigma = p0_mu, abs(p0_sigma)
        
        result = DepthDistributionResult(
            mu=mu,
            sigma=sigma,
            depths=depths,
            histogram=hist,
            bin_edges=bin_edges
        )
        
        self.result = result
        
        # Write output if requested
        if output_file:
            self._write_output(result, output_file)
        
        self.logger.info(f"Depth distribution: mu={mu:.2f}, sigma={sigma:.2f}")
        
        return result
    
    def _write_output(self, result: DepthDistributionResult, output_file: str):
        """Write distribution results to file."""
        with open(output_file, 'w') as f:
            f.write(f'Average depth per copy: {result.mu}\n')
            f.write(f'Standard deviation of depths per copy: {result.sigma}\n\n')
            f.write('Depth distribution:\n\n')
            f.write('\t'.join(['Depth', 'Frequency', 'Fit value']) + '\n')
            
            if result.histogram is not None and result.bin_edges is not None:
                for hist_val, bin_val in zip(result.histogram, result.bin_edges[:-1]):
                    fit_val = norm.pdf(bin_val, result.mu, result.sigma)
                    f.write(f'{int(bin_val)}\t{hist_val}\t{fit_val}\n')
    
    def filter_regions_by_depth(self, depths_by_chrom: Dict[str, np.ndarray],
                                 chrom_sizes: Dict[str, int],
                                 output_file: Optional[str] = None) -> Dict[str, List[Tuple[int, int]]]:
        """
        Filter genomic regions with abnormal depth.
        
        Uses a sliding window approach to identify regions where the depth
        deviates significantly from the expected distribution.
        
        Parameters
        ----------
        depths_by_chrom : dict
            Dictionary mapping chromosome names to arrays of depth values
        chrom_sizes : dict
            Dictionary mapping chromosome names to sizes
        output_file : str, optional
            Path to BED file for filtered regions
            
        Returns
        -------
        dict
            Dictionary mapping chromosome names to lists of (start, end) tuples
        """
        if self.result is None:
            raise ValueError("Must fit distribution before filtering regions")
        
        filtered_regions = {}
        
        for chromosome in sorted(depths_by_chrom.keys()):
            chromosome_depths = depths_by_chrom[chromosome]
            regions = self._process_chromosome(
                chromosome, chromosome_depths,
                self.result.mu, self.result.sigma
            )
            if regions:
                filtered_regions[chromosome] = regions
        
        self.result.filtered_regions = filtered_regions
        
        # Write BED file if requested
        if output_file:
            with open(output_file, 'w') as f:
                for chrom, regions in sorted(filtered_regions.items()):
                    for start, end in regions:
                        f.write(f'{chrom}\t{start}\t{end}\n')
        
        total_filtered = sum(len(r) for r in filtered_regions.values())
        self.logger.info(f"Identified {total_filtered} filtered regions")
        
        return filtered_regions
    
    def _process_chromosome(self, chromosome: str, depths: np.ndarray,
                            mu: float, sigma: float) -> List[Tuple[int, int]]:
        """
        Process a single chromosome to identify filtered regions.
        
        Parameters
        ----------
        chromosome : str
            Chromosome name
        depths : np.ndarray
            Array of depth values
        mu : float
            Mean depth
        sigma : float
            Standard deviation of depth
            
        Returns
        -------
        list
            List of (start, end) tuples for filtered regions
        """
        d = int((self.window_size - 1) / 2)
        norm_dist = norm(mu, sigma)
        
        # Thresholds for keeping and filtering
        keep_threshold = [mu, mu]
        filter_threshold = [float('-inf'), float('inf')]
        
        intervals = []
        first = float('inf')
        last = float('-inf')
        side = 0
        last_side = 0
        
        max_pos = len(depths)
        
        for i in range(max_pos):
            # Calculate window bounds
            window_start = max(0, i - d)
            window_end = min(max_pos, i + d + 1)
            
            window_depth = np.mean(depths[window_start:window_end])
            
            # Check if outside acceptable range
            if not (keep_threshold[0] <= window_depth <= keep_threshold[1]):
                if window_depth <= filter_threshold[0] or window_depth >= filter_threshold[1]:
                    side = -1 if window_depth < mu else 1
                    
                    if i - last > self.merge_window or last_side * side == -1:
                        if last - first > 0:
                            intervals.append((first, last + 1))
                        first = i
                    
                    last = i
                    last_side = side
                else:
                    if window_depth < mu:
                        side = -1
                        p = norm_dist.cdf(window_depth)
                        
                        if p >= self.p_threshold:
                            keep_threshold[0] = window_depth
                        else:
                            filter_threshold[0] = window_depth
                            if i - last > self.merge_window or last_side * side == -1:
                                if last - first > 0:
                                    intervals.append((first, last + 1))
                                first = i
                            last = i
                            last_side = side
                    
                    elif window_depth > mu:
                        side = 1
                        p = 1.0 - norm_dist.cdf(window_depth)
                        
                        if p >= self.p_threshold:
                            keep_threshold[1] = window_depth
                        else:
                            filter_threshold[1] = window_depth
                            if i - last > self.merge_window or last_side * side == -1:
                                if last - first > 0:
                                    intervals.append((first, last + 1))
                                first = i
                            last = i
                            last_side = side
        
        # Final interval
        if last - first > 0:
            intervals.append((first, last + 1))
        
        return intervals
    
    def get_filtered_positions(self) -> Set[Tuple[str, int]]:
        """
        Get a set of all filtered positions.
        
        Returns
        -------
        set
            Set of (chromosome, position) tuples
        """
        if self.result is None or not self.result.filtered_regions:
            return set()
        
        positions = set()
        for chrom, regions in self.result.filtered_regions.items():
            for start, end in regions:
                for pos in range(start, end):
                    positions.add((chrom, pos))
        
        return positions
    
    def is_position_filtered(self, chromosome: str, position: int) -> bool:
        """
        Check if a position is in a filtered region.
        
        Parameters
        ----------
        chromosome : str
            Chromosome name
        position : int
            Position to check
            
        Returns
        -------
        bool
            True if position is in a filtered region
        """
        if self.result is None or not self.result.filtered_regions:
            return False
        
        regions = self.result.filtered_regions.get(chromosome, [])
        for start, end in regions:
            if start <= position < end:
                return True
        return False


def calculate_depth_distribution(bam_file: str, reference: str,
                                  output_file: Optional[str] = None,
                                  config: dict = None) -> DepthDistributionResult:
    """
    Convenience function to calculate depth distribution.
    
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
    DepthDistributionResult
        Fitted distribution parameters
    """
    analyzer = DepthDistributionAnalyzer(config)
    return analyzer.calculate_from_bam(bam_file, reference, output_file)


def filter_regions_by_depth(bedgraph_file: str, chrom_sizes: Dict[str, int],
                            mu: float, sigma: float,
                            output_file: str,
                            ploidy: int = 1,
                            p_threshold: float = 0.0001,
                            merge_window: int = 1000) -> Dict[str, List[Tuple[int, int]]]:
    """
    Filter regions by depth from a bedGraph file.
    
    Parameters
    ----------
    bedgraph_file : str
        Path to bedGraph file
    chrom_sizes : dict
        Chromosome sizes
    mu, sigma : float
        Distribution parameters
    output_file : str
        Path to output BED file
    ploidy : int
        Ploidy level
    p_threshold : float
        P-value threshold
    merge_window : int
        Window for merging adjacent regions
        
    Returns
    -------
    dict
        Filtered regions by chromosome
    """
    config = {
        'ploidy': ploidy,
        'p_threshold': p_threshold,
        'merge_window': merge_window
    }
    
    analyzer = DepthDistributionAnalyzer(config)
    analyzer.result = DepthDistributionResult(mu=mu, sigma=sigma)
    
    # Read depths from bedGraph
    depths_by_chrom = {}
    for chrom, size in chrom_sizes.items():
        depths_by_chrom[chrom] = np.zeros(size, dtype=np.int32)
    
    with open(bedgraph_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 4:
                chrom, start, end, coverage = fields[:4]
                start, end = int(start), int(end)
                coverage = float(coverage)
                
                if chrom in depths_by_chrom:
                    for pos in range(start, min(end, len(depths_by_chrom[chrom]))):
                        depths_by_chrom[chrom][pos] = int(coverage / ploidy)
    
    return analyzer.filter_regions_by_depth(depths_by_chrom, chrom_sizes, output_file)
