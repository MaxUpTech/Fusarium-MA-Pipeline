#!/usr/bin/env python3
"""
Depth Correction Model for Fusarium MA Pipeline

Implements chromosome end depth correction from the muver pipeline.
Corrects for systematic depth biases near chromosome ends using a
log-normal CDF combined with a linear function.
"""

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass
class DepthCorrectionParameters:
    """Parameters for depth correction model."""
    y_int: float  # Y-intercept of linear component
    scalar: float  # Scalar for log-normal component
    mean_log: float  # Mean of log-normal distribution
    sd_log: float  # Standard deviation of log-normal distribution
    slope: float  # Slope of linear component


class DepthCorrector:
    """
    Corrects read depth for systematic biases near chromosome ends.
    
    Uses a model combining a log-normal CDF and linear function to
    correct for position-dependent depth biases commonly observed
    near chromosome termini.
    """
    
    def __init__(self, config: dict = None):
        self.config = config or {}
        self.logger = logging.getLogger('FusariumMA.DepthCorrection')
        self.params: Optional[DepthCorrectionParameters] = None
        self.chrom_sizes: Dict[str, int] = {}
    
    def set_parameters(self, y_int: float, scalar: float, mean_log: float,
                       sd_log: float, slope: float):
        """
        Set depth correction parameters.
        
        Parameters
        ----------
        y_int : float
            Y-intercept of linear component
        scalar : float
            Scalar for log-normal CDF component
        mean_log : float
            Mean of log-normal distribution
        sd_log : float
            Standard deviation of log-normal distribution
        slope : float
            Slope of linear component
        """
        self.params = DepthCorrectionParameters(
            y_int=y_int,
            scalar=scalar,
            mean_log=mean_log,
            sd_log=sd_log,
            slope=slope
        )
        self.logger.info(f"Depth correction parameters set: {self.params}")
    
    def set_chromosome_sizes(self, chrom_sizes: Dict[str, int]):
        """
        Set chromosome sizes for position calculations.
        
        Parameters
        ----------
        chrom_sizes : dict
            Dictionary mapping chromosome names to sizes
        """
        self.chrom_sizes = chrom_sizes
    
    def load_chromosome_sizes(self, fai_file: str):
        """
        Load chromosome sizes from a FASTA index file.
        
        Parameters
        ----------
        fai_file : str
            Path to .fai index file
        """
        self.chrom_sizes = {}
        with open(fai_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    self.chrom_sizes[fields[0]] = int(fields[1])
        
        self.logger.info(f"Loaded {len(self.chrom_sizes)} chromosome sizes")
    
    def calculate_correction_factor(self, chromosome: str, position: int) -> float:
        """
        Calculate the depth correction factor for a genomic position.
        
        The correction factor is based on the distance to the nearest
        chromosome end, modeled using a log-normal CDF combined with
        a linear function.
        
        Parameters
        ----------
        chromosome : str
            Chromosome name
        position : int
            Genomic position (1-based)
            
        Returns
        -------
        float
            Correction factor (multiply depth by this to correct)
        """
        if self.params is None:
            return 1.0
        
        if chromosome not in self.chrom_sizes:
            return 1.0
        
        chrom_size = self.chrom_sizes[chromosome]
        
        # Calculate distance to nearest chromosome end
        relative_pos = min(position - 1, chrom_size - position)
        if relative_pos <= 0:
            relative_pos = 1
        
        # Calculate correction factor using log-normal CDF + linear
        try:
            log_term = (self.params.mean_log - math.log(relative_pos)) / \
                       (math.sqrt(2) * self.params.sd_log)
            lognorm_cdf = 0.5 + 0.5 * math.erf(log_term)
            
            correction = self.params.scalar * lognorm_cdf + \
                        self.params.y_int + \
                        self.params.slope * relative_pos
            
            return correction
        except (ValueError, ZeroDivisionError):
            return 1.0
    
    def correct_depth(self, chromosome: str, position: int, raw_depth: float) -> float:
        """
        Apply depth correction to a raw depth value.
        
        Parameters
        ----------
        chromosome : str
            Chromosome name
        position : int
            Genomic position
        raw_depth : float
            Raw read depth
            
        Returns
        -------
        float
            Corrected depth
        """
        correction = self.calculate_correction_factor(chromosome, position)
        if correction <= 0:
            return raw_depth
        return raw_depth / correction
    
    def correct_bedgraph(self, input_bedgraph: str, output_bedgraph: str):
        """
        Apply depth correction to an entire bedGraph file.
        
        Parameters
        ----------
        input_bedgraph : str
            Path to input bedGraph file
        output_bedgraph : str
            Path to output corrected bedGraph file
        """
        if self.params is None:
            raise ValueError("Must set parameters before correcting")
        
        if not self.chrom_sizes:
            raise ValueError("Must set chromosome sizes before correcting")
        
        self.logger.info(f"Correcting depths in {input_bedgraph}")
        
        last_chr = None
        last_val = None
        last_pos = None
        last_start = None
        
        with open(input_bedgraph, 'r') as fin, open(output_bedgraph, 'w') as fout:
            
            def write_interval(chrom, start, end, value):
                if value > 0:
                    fout.write(f'{chrom}\t{start - 1}\t{end}\t{value}\n')
            
            for line in fin:
                fields = line.strip().split()
                if len(fields) < 4:
                    continue
                
                chromosome, start, end, coverage = fields[:4]
                start = int(start) + 1  # Convert to 1-based
                end = int(end)
                coverage = float(coverage)
                
                for position in range(start, end + 1):
                    corrected = self.correct_depth(chromosome, position, coverage)
                    value = int(math.floor(corrected))
                    
                    # Initialize or continue interval
                    if last_chr is None:
                        last_chr = chromosome
                        last_val = value
                        last_pos = position
                        last_start = position
                    elif (chromosome != last_chr or value != last_val or 
                          position != last_pos + 1):
                        write_interval(last_chr, last_start, last_pos, last_val)
                        last_start = position
                    
                    last_chr = chromosome
                    last_val = value
                    last_pos = position
            
            # Write final interval
            if last_chr is not None:
                write_interval(last_chr, last_start, last_pos, last_val)
        
        self.logger.info(f"Corrected bedGraph written to {output_bedgraph}")


def calculate_depth_correction_parameters(depth_ratios: Dict[int, List[float]],
                                          max_distance: int = 50000) -> DepthCorrectionParameters:
    """
    Calculate depth correction parameters from observed depth ratios.
    
    This function fits the depth correction model to observed depth
    ratios as a function of distance from chromosome ends.
    
    Parameters
    ----------
    depth_ratios : dict
        Dictionary mapping distances to lists of depth ratios
    max_distance : int
        Maximum distance to consider for fitting
        
    Returns
    -------
    DepthCorrectionParameters
        Fitted model parameters
    """
    from scipy.optimize import curve_fit
    
    # Prepare data
    distances = []
    ratios = []
    
    for dist, ratio_list in depth_ratios.items():
        if dist <= max_distance and ratio_list:
            median_ratio = np.median(ratio_list)
            distances.append(dist)
            ratios.append(median_ratio)
    
    distances = np.array(distances)
    ratios = np.array(ratios)
    
    # Define model function
    def model(x, y_int, scalar, mean_log, sd_log, slope):
        result = np.zeros_like(x, dtype=float)
        for i, xi in enumerate(x):
            if xi > 0:
                log_term = (mean_log - np.log(xi)) / (np.sqrt(2) * sd_log)
                result[i] = scalar * (0.5 + 0.5 * math.erf(log_term)) + \
                           y_int + slope * xi
            else:
                result[i] = 1.0
        return result
    
    # Initial parameter estimates
    p0 = [1.0, 0.5, 8.0, 2.0, 0.0]
    
    try:
        popt, pcov = curve_fit(model, distances, ratios, p0=p0, maxfev=10000)
        y_int, scalar, mean_log, sd_log, slope = popt
        
        return DepthCorrectionParameters(
            y_int=y_int,
            scalar=scalar,
            mean_log=mean_log,
            sd_log=abs(sd_log),
            slope=slope
        )
    except:
        # Return default parameters if fitting fails
        return DepthCorrectionParameters(
            y_int=1.0,
            scalar=0.0,
            mean_log=8.0,
            sd_log=2.0,
            slope=0.0
        )


def apply_depth_correction(variants: List[dict], corrector: DepthCorrector) -> List[dict]:
    """
    Apply depth correction to variant calls.
    
    Parameters
    ----------
    variants : list
        List of variant dictionaries
    corrector : DepthCorrector
        Configured depth corrector
        
    Returns
    -------
    list
        Variants with corrected depth values
    """
    for var in variants:
        chrom = var.get('chromosome', var.get('chrom'))
        pos = var.get('position', var.get('pos'))
        
        if chrom and pos:
            # Correct depths
            if 'depth' in var:
                var['corrected_depth'] = corrector.correct_depth(chrom, pos, var['depth'])
            if 'sample_depth' in var:
                var['corrected_sample_depth'] = corrector.correct_depth(
                    chrom, pos, var['sample_depth']
                )
            if 'ancestor_depth' in var:
                var['corrected_ancestor_depth'] = corrector.correct_depth(
                    chrom, pos, var['ancestor_depth']
                )
    
    return variants
