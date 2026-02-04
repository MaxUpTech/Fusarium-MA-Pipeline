#!/usr/bin/env python3
"""
Repeat INDEL Correction Model for Fusarium MA Pipeline

Implements repeat-aware INDEL error rate correction from the muver pipeline.
Calculates expected indel rates within repeat regions and adjusts variant
calling thresholds accordingly.
"""

import csv
import logging
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import curve_fit

from .fitting import logistic, fit_logistic


@dataclass
class RepeatRegion:
    """Represents a repeat region in the genome."""
    chromosome: str
    start: int
    end: int
    unit: str  # Repeat unit sequence
    sequence: str  # Full repeat sequence
    
    @property
    def unit_length(self) -> int:
        return len(self.unit)
    
    @property
    def tract_length(self) -> int:
        return self.unit_length * self.sequence.count(self.unit)


@dataclass
class RepeatIndelFit:
    """Fitted parameters for repeat indel rates."""
    event_type: str  # 'insertion' or 'deletion'
    repeat_length: int  # Repeat unit length (1-4 bp)
    x0: float  # Midpoint
    L: float  # Amplitude
    M: float  # Baseline
    k: float  # Steepness


@dataclass
class RepeatIndelResult:
    """Results from repeat indel analysis."""
    fits: Dict[str, Dict[int, RepeatIndelFit]] = field(default_factory=dict)
    rates: Dict[str, Dict[int, Dict[int, float]]] = field(default_factory=dict)
    counts: Dict[str, Dict[int, Dict[int, int]]] = field(default_factory=dict)


class RepeatIndelAnalyzer:
    """
    Analyzes and corrects for INDEL errors in repeat regions.
    
    Repeat regions are prone to sequencing errors, particularly insertions
    and deletions of repeat units. This class models the error rates as a
    function of repeat unit length and tract length using logistic functions.
    """
    
    def __init__(self, config: dict = None):
        self.config = config or {}
        self.logger = logging.getLogger('FusariumMA.RepeatIndels')
        self.repeats: Dict[str, Dict[int, List[RepeatRegion]]] = {}
        self.result: Optional[RepeatIndelResult] = None
        
        # Minimum occurrences to report a rate
        self.occurrence_filter = self.config.get('occurrence_filter', 10)
    
    def load_repeat_file(self, repeat_file: str):
        """
        Load repeat regions from a file.
        
        Parameters
        ----------
        repeat_file : str
            Path to repeat file (BED-like format)
        """
        self.logger.info(f"Loading repeats from {repeat_file}")
        
        self.repeats = defaultdict(lambda: defaultdict(list))
        
        with open(repeat_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                
                chrom, start, end, unit, sequence = fields[:5]
                start, end = int(start), int(end)
                
                repeat = RepeatRegion(
                    chromosome=chrom,
                    start=start,
                    end=end,
                    unit=unit,
                    sequence=sequence
                )
                
                # Index by position for fast lookup
                for pos in range(start, end + 1):
                    self.repeats[chrom][pos].append(repeat)
        
        total = sum(len(pos_dict) for pos_dict in self.repeats.values())
        self.logger.info(f"Loaded {total} repeat positions")
    
    def generate_repeat_file(self, reference: str, output_file: str,
                             min_unit_length: int = 1,
                             max_unit_length: int = 4,
                             min_tract_length: int = 4):
        """
        Generate a repeat file from reference sequence.
        
        Parameters
        ----------
        reference : str
            Path to reference FASTA
        output_file : str
            Path to output repeat file
        min_unit_length : int
            Minimum repeat unit length
        max_unit_length : int
            Maximum repeat unit length
        min_tract_length : int
            Minimum tract length (in bp)
        """
        self.logger.info(f"Generating repeat file from {reference}")
        
        from Bio import SeqIO
        
        with open(output_file, 'w') as out:
            for record in SeqIO.parse(reference, 'fasta'):
                chrom = record.id
                seq = str(record.seq).upper()
                
                for unit_len in range(min_unit_length, max_unit_length + 1):
                    i = 0
                    while i < len(seq) - unit_len:
                        unit = seq[i:i + unit_len]
                        
                        # Skip units with N
                        if 'N' in unit:
                            i += 1
                            continue
                        
                        # Find repeat tract
                        j = i + unit_len
                        while j < len(seq) and seq[j:j + unit_len] == unit:
                            j += unit_len
                        
                        tract_len = j - i
                        if tract_len >= min_tract_length:
                            out.write(f'{chrom}\t{i}\t{j}\t{unit}\t{seq[i:j]}\n')
                            i = j
                        else:
                            i += 1
        
        self.logger.info(f"Repeat file written to {output_file}")
        self.load_repeat_file(output_file)
    
    def calculate_repeat_occurrences(self) -> Dict[Tuple[int, int], int]:
        """
        Count occurrences of repeats by unit length and tract length.
        
        Returns
        -------
        dict
            Dictionary mapping (unit_length, tract_length) to count
        """
        occurrences = defaultdict(int)
        
        for chrom_repeats in self.repeats.values():
            seen_repeats = set()
            for pos_repeats in chrom_repeats.values():
                for repeat in pos_repeats:
                    # Only count each repeat once (at its start)
                    repeat_id = (repeat.chromosome, repeat.start, repeat.unit)
                    if repeat_id not in seen_repeats:
                        seen_repeats.add(repeat_id)
                        occurrences[(repeat.unit_length, repeat.tract_length)] += 1
        
        return dict(occurrences)
    
    def calculate_indel_counts_from_bam(self, bam_file: str) -> Dict[str, Dict[int, Dict[int, int]]]:
        """
        Count indels in repeat regions from a BAM file.
        
        Parameters
        ----------
        bam_file : str
            Path to BAM file
            
        Returns
        -------
        dict
            Nested dictionary: event_type -> unit_length -> tract_length -> count
        """
        import pysam
        
        self.logger.info(f"Counting indels in repeats from {bam_file}")
        
        counts = {
            'depth': {i: defaultdict(int) for i in range(1, 5)},
            'insertion': {i: defaultdict(int) for i in range(1, 5)},
            'deletion': {i: defaultdict(int) for i in range(1, 5)}
        }
        
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            for read in bam.fetch():
                if read.is_unmapped or read.is_secondary:
                    continue
                
                chrom = read.reference_name
                ref_pos = read.reference_start
                
                # Parse CIGAR
                for op, length in read.cigartuples:
                    if op == 0:  # Match
                        for i in range(length):
                            pos = ref_pos + i
                            if chrom in self.repeats and pos in self.repeats[chrom]:
                                for repeat in self.repeats[chrom][pos]:
                                    if repeat.start == pos and 1 <= repeat.unit_length <= 4:
                                        counts['depth'][repeat.unit_length][repeat.tract_length] += 1
                        ref_pos += length
                    
                    elif op == 1:  # Insertion
                        if chrom in self.repeats and ref_pos in self.repeats[chrom]:
                            for repeat in self.repeats[chrom][ref_pos]:
                                if repeat.start == ref_pos and 1 <= repeat.unit_length <= 4:
                                    counts['insertion'][repeat.unit_length][repeat.tract_length] += 1
                    
                    elif op == 2:  # Deletion
                        if chrom in self.repeats and ref_pos in self.repeats[chrom]:
                            for repeat in self.repeats[chrom][ref_pos]:
                                if repeat.start == ref_pos and 1 <= repeat.unit_length <= 4:
                                    counts['deletion'][repeat.unit_length][repeat.tract_length] += 1
                        ref_pos += length
                    
                    elif op in [3, 7, 8]:  # Skip, match, mismatch
                        ref_pos += length
        
        return counts
    
    def calculate_indel_rates(self, indel_counts: Dict,
                              repeat_occurrences: Dict[Tuple[int, int], int]) -> Dict[str, Dict[int, Dict[int, float]]]:
        """
        Calculate indel rates from counts.
        
        Parameters
        ----------
        indel_counts : dict
            Counts from calculate_indel_counts_from_bam
        repeat_occurrences : dict
            Repeat occurrences from calculate_repeat_occurrences
            
        Returns
        -------
        dict
            Nested dictionary: event_type -> unit_length -> tract_length -> rate
        """
        rates = {
            'insertion': {i: {} for i in range(1, 5)},
            'deletion': {i: {} for i in range(1, 5)}
        }
        
        for event in ['insertion', 'deletion']:
            for unit_len in range(1, 5):
                for tract_len, count in indel_counts[event][unit_len].items():
                    if repeat_occurrences.get((unit_len, tract_len), 0) >= self.occurrence_filter:
                        depth = indel_counts['depth'][unit_len].get(tract_len, 0)
                        if depth > 0:
                            rates[event][unit_len][tract_len] = float(count) / depth
        
        return rates
    
    def fit_indel_rates(self, indel_rates: Dict) -> Dict[str, Dict[int, RepeatIndelFit]]:
        """
        Fit logistic functions to indel rates.
        
        Parameters
        ----------
        indel_rates : dict
            Rates from calculate_indel_rates
            
        Returns
        -------
        dict
            Nested dictionary: event_type -> unit_length -> RepeatIndelFit
        """
        fits = {'insertion': {}, 'deletion': {}}
        
        for event in ['insertion', 'deletion']:
            for unit_len, rates in indel_rates[event].items():
                if not rates:
                    continue
                
                tract_lengths = list(rates.keys())
                log_rates = [math.log10(r) for r in rates.values()]
                
                if len(tract_lengths) < 4:
                    continue
                
                fit_params = fit_logistic(tract_lengths, log_rates)
                
                if fit_params:
                    fits[event][unit_len] = RepeatIndelFit(
                        event_type=event,
                        repeat_length=unit_len,
                        **fit_params
                    )
        
        return fits
    
    def fit_repeat_indel_rates(self, bam_file: str,
                                output_file: Optional[str] = None,
                                plot_header: Optional[str] = None) -> RepeatIndelResult:
        """
        Complete analysis: calculate counts, rates, and fits.
        
        Parameters
        ----------
        bam_file : str
            Path to BAM file
        output_file : str, optional
            Path to output fits file
        plot_header : str, optional
            Header for output plots
            
        Returns
        -------
        RepeatIndelResult
            Complete analysis results
        """
        # Calculate occurrences
        occurrences = self.calculate_repeat_occurrences()
        
        # Calculate counts from BAM
        counts = self.calculate_indel_counts_from_bam(bam_file)
        
        # Calculate rates
        rates = self.calculate_indel_rates(counts, occurrences)
        
        # Fit rates
        fits = self.fit_indel_rates(rates)
        
        result = RepeatIndelResult(
            fits=fits,
            rates=rates,
            counts=counts
        )
        
        self.result = result
        
        # Write output file
        if output_file:
            self.write_fits(fits, output_file)
        
        # Generate plots
        if plot_header:
            self.plot_fits(rates, fits, plot_header)
        
        return result
    
    def write_fits(self, fits: Dict, output_file: str):
        """Write fit parameters to a TSV file."""
        with open(output_file, 'w') as f:
            writer = csv.DictWriter(f, 
                                    fieldnames=['Event', 'Repeat length', 'x0', 'L', 'M', 'k'],
                                    delimiter='\t')
            writer.writeheader()
            
            for event in ['insertion', 'deletion']:
                for unit_len in range(1, 5):
                    if unit_len in fits.get(event, {}):
                        fit = fits[event][unit_len]
                        writer.writerow({
                            'Event': event,
                            'Repeat length': unit_len,
                            'x0': fit.x0,
                            'L': fit.L,
                            'M': fit.M,
                            'k': fit.k
                        })
    
    def read_fits(self, fits_file: str) -> Dict[str, Dict[int, RepeatIndelFit]]:
        """Read fit parameters from a TSV file."""
        fits = {'insertion': {}, 'deletion': {}}
        
        with open(fits_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                event = row['Event']
                unit_len = int(row['Repeat length'])
                
                fits[event][unit_len] = RepeatIndelFit(
                    event_type=event,
                    repeat_length=unit_len,
                    x0=float(row['x0']),
                    L=float(row['L']),
                    M=float(row['M']),
                    k=float(row['k'])
                )
        
        return fits
    
    def plot_fits(self, rates: Dict, fits: Dict, output_header: str):
        """Generate plots comparing fits to observed rates."""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            self.logger.warning("matplotlib not available, skipping plots")
            return
        
        for event, event_fits in fits.items():
            for unit_len, fit in event_fits.items():
                event_rates = rates[event][unit_len]
                
                tract_lengths = sorted(event_rates.keys())
                observed = [math.log10(event_rates[t]) for t in tract_lengths]
                fitted = [logistic(t, fit.x0, fit.L, fit.M, fit.k) for t in tract_lengths]
                
                plt.figure(figsize=(8, 6))
                plt.plot(tract_lengths, observed, 'ko', label='Observed')
                plt.plot(tract_lengths, fitted, 'b-', label='Fitted')
                plt.xlabel('Repeat tract length (bp)')
                plt.ylabel('log10(indel rate)')
                plt.title(f'{event.capitalize()} rate vs tract length\n(repeat unit: {unit_len} bp)')
                plt.legend()
                plt.tight_layout()
                plt.savefig(f'{output_header}_{event}_{unit_len}.png', dpi=300)
                plt.close()
    
    def get_repeat_at_position(self, chromosome: str, position: int) -> Optional[RepeatRegion]:
        """
        Get the repeat region at a given position.
        
        Parameters
        ----------
        chromosome : str
            Chromosome name
        position : int
            Position to query
            
        Returns
        -------
        RepeatRegion or None
            Repeat region if found, else None
        """
        if chromosome not in self.repeats:
            return None
        
        if position not in self.repeats[chromosome]:
            return None
        
        repeats = self.repeats[chromosome][position]
        if not repeats:
            return None
        
        # Return the repeat with smallest unit length
        return min(repeats, key=lambda r: len(r.unit))


def get_repeat_adjustment_value(unit_length: int, tract_length: int,
                                 event_type: str,
                                 fits: Dict[str, Dict[int, RepeatIndelFit]]) -> float:
    """
    Get expected indel frequency for a repeat region.
    
    Parameters
    ----------
    unit_length : int
        Repeat unit length
    tract_length : int
        Repeat tract length
    event_type : str
        'insertion' or 'deletion'
    fits : dict
        Fitted parameters
        
    Returns
    -------
    float
        Expected frequency (0-1)
    """
    if event_type not in fits or unit_length not in fits[event_type]:
        return 0.0
    
    fit = fits[event_type][unit_length]
    
    log_rate = logistic(tract_length, fit.x0, fit.L, fit.M, fit.k)
    return 10 ** min(0, log_rate)


def calculate_repeat_indel_rates(bam_file: str, repeat_file: str,
                                  output_file: str,
                                  config: dict = None) -> RepeatIndelResult:
    """
    Convenience function to calculate and fit repeat indel rates.
    
    Parameters
    ----------
    bam_file : str
        Path to BAM file
    repeat_file : str
        Path to repeat file
    output_file : str
        Path to output fits file
    config : dict, optional
        Configuration dictionary
        
    Returns
    -------
    RepeatIndelResult
        Analysis results
    """
    analyzer = RepeatIndelAnalyzer(config)
    analyzer.load_repeat_file(repeat_file)
    return analyzer.fit_repeat_indel_rates(bam_file, output_file)
