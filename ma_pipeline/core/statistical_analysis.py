#!/usr/bin/env python3
"""
Enhanced Statistical Analysis Module for Fusarium MA Pipeline

Comprehensive statistical analysis combining methods from JGI pipeline and Muver
with additional tests for increased statistical power:
- Kolmogorov-Smirnov tests
- Shapiro-Wilk normality tests
- Permutation tests
- Effect size calculations (Cohen's d, Glass's delta)
- Bootstrap analysis for spectrum comparison
- Bayesian credible intervals
"""

import logging
import json
import math
import warnings
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict, field
import numpy as np
from scipy import stats


@dataclass
class MutationRateResult:
    """Mutation rate calculation result for a sample."""
    sample_name: str
    genotype: str
    generation: int
    total_mutations: int
    snp_count: int
    indel_count: int
    insertion_count: int
    deletion_count: int
    callable_sites: int
    rate_per_site_per_gen: float
    rate_per_kb_per_gen: float
    ci_lower: float  # 95% CI lower bound (per kb)
    ci_upper: float  # 95% CI upper bound (per kb)
    bootstrap_ci_lower: float = 0.0  # Bootstrap CI lower
    bootstrap_ci_upper: float = 0.0  # Bootstrap CI upper
    bayesian_ci_lower: float = 0.0  # Bayesian credible interval lower
    bayesian_ci_upper: float = 0.0  # Bayesian credible interval upper
    

@dataclass
class SpectrumResult:
    """Mutation spectrum analysis result."""
    sample_name: str
    # 6-category spectrum (collapsed from 12)
    C_to_A: int = 0
    C_to_G: int = 0
    C_to_T: int = 0
    T_to_A: int = 0
    T_to_C: int = 0
    T_to_G: int = 0
    # Derived metrics
    transitions: int = 0
    transversions: int = 0
    ts_tv_ratio: float = 0.0
    gc_to_at_bias: float = 0.0


@dataclass
class EffectSizeResult:
    """Effect size calculation result."""
    genotype1: str
    genotype2: str
    cohens_d: float
    cohens_d_ci_lower: float
    cohens_d_ci_upper: float
    glass_delta: float
    hedges_g: float
    interpretation: str  # small, medium, large


@dataclass
class PermutationTestResult:
    """Permutation test result."""
    test_name: str
    observed_statistic: float
    p_value: float
    n_permutations: int
    null_distribution_mean: float
    null_distribution_std: float


@dataclass
class NormalityTestResult:
    """Normality test results."""
    test_name: str
    statistic: float
    p_value: float
    is_normal: bool  # at alpha=0.05


@dataclass 
class DistributionTestResult:
    """Distribution comparison test result."""
    test_name: str
    statistic: float
    p_value: float
    conclusion: str


class EnhancedStatistics:
    """
    Enhanced statistical analysis for MA experiments.
    
    Outputs:
    - Mutation rates per sample and overall with multiple CI methods
    - Mutation spectrum (6-category and 96-category)
    - Ts/Tv ratio
    - GC→AT bias
    - INDEL size distribution
    - Per-chromosome distribution
    - Bootstrap confidence intervals
    - Poisson test for mutation rate variation
    - Kolmogorov-Smirnov tests
    - Shapiro-Wilk normality tests
    - Permutation tests
    - Effect size calculations
    - Bayesian analysis
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.Statistics')
        self.output_dir = Path(config['output']['directory']) / 'statistics'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.generations = config.get('statistics', {}).get('generations', 25)
        self.bootstrap_n = config.get('statistics', {}).get('bootstrap_iterations', 10000)
        self.permutation_n = config.get('statistics', {}).get('permutation_iterations', 10000)
        
        # Store results for plotting
        self.results_cache = {}
        
    def calculate_callable_sites(self, bam_file: str = None, 
                                  reference: str = None) -> int:
        """Calculate number of callable sites."""
        if 'callable_sites' in self.config.get('statistics', {}):
            return self.config['statistics']['callable_sites']
        
        reference = self.config['reference']['fasta']
        try:
            total = 0
            with open(reference, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        total += len(line.strip())
            return total
        except:
            return 36_000_000
    
    def calculate_mutation_rates(self, mutations: List, samples: List,
                                  ancestor) -> List[MutationRateResult]:
        """Calculate mutation rates with multiple CI methods."""
        self.logger.info("Calculating mutation rates...")
        
        callable_sites = self.calculate_callable_sites()
        self.logger.info(f"Callable sites: {callable_sites:,}")
        
        mutations_by_sample = defaultdict(list)
        for mut in mutations:
            mutations_by_sample[mut.sample_name].append(mut)
        
        results = []
        for sample in samples:
            if sample.sample_name == ancestor.sample_name:
                continue
            
            sample_muts = mutations_by_sample[sample.sample_name]
            
            snps = [m for m in sample_muts if m.mutation_type == 'SNP']
            indels = [m for m in sample_muts if m.mutation_type in ['INS', 'DEL']]
            insertions = [m for m in sample_muts if m.mutation_type == 'INS']
            deletions = [m for m in sample_muts if m.mutation_type == 'DEL']
            
            total = len(sample_muts)
            
            rate_per_site = total / (callable_sites * self.generations) if callable_sites > 0 else 0
            rate_per_kb = rate_per_site * 1000
            
            # Poisson-based CI
            ci_lower, ci_upper = self._poisson_ci(total, callable_sites, self.generations)
            
            # Bootstrap CI
            boot_lower, boot_upper = self._bootstrap_rate_ci(total, callable_sites, self.generations)
            
            # Bayesian credible interval (Gamma-Poisson conjugate)
            bayes_lower, bayes_upper = self._bayesian_credible_interval(total, callable_sites, self.generations)
            
            result = MutationRateResult(
                sample_name=sample.sample_name,
                genotype=sample.genotype,
                generation=sample.generation,
                total_mutations=total,
                snp_count=len(snps),
                indel_count=len(indels),
                insertion_count=len(insertions),
                deletion_count=len(deletions),
                callable_sites=callable_sites,
                rate_per_site_per_gen=rate_per_site,
                rate_per_kb_per_gen=rate_per_kb,
                ci_lower=ci_lower,
                ci_upper=ci_upper,
                bootstrap_ci_lower=boot_lower,
                bootstrap_ci_upper=boot_upper,
                bayesian_ci_lower=bayes_lower,
                bayesian_ci_upper=bayes_upper
            )
            results.append(result)
        
        return results
    
    def _poisson_ci(self, n_mutations: int, callable_sites: int,
                    generations: int, alpha: float = 0.05) -> Tuple[float, float]:
        """Calculate Poisson-based confidence interval."""
        if n_mutations == 0:
            return (0.0, 0.0)
        
        lower = stats.chi2.ppf(alpha/2, 2*n_mutations) / (2 * callable_sites * generations)
        upper = stats.chi2.ppf(1-alpha/2, 2*(n_mutations+1)) / (2 * callable_sites * generations)
        
        return (lower * 1000, upper * 1000)
    
    def _bootstrap_rate_ci(self, n_mutations: int, callable_sites: int,
                           generations: int, alpha: float = 0.05) -> Tuple[float, float]:
        """Calculate bootstrap confidence interval for mutation rate."""
        if n_mutations == 0:
            return (0.0, 0.0)
        
        # Bootstrap by resampling Poisson counts
        bootstrap_rates = []
        for _ in range(self.bootstrap_n):
            boot_count = np.random.poisson(n_mutations)
            boot_rate = boot_count / (callable_sites * generations) * 1000
            bootstrap_rates.append(boot_rate)
        
        lower = np.percentile(bootstrap_rates, 100 * alpha/2)
        upper = np.percentile(bootstrap_rates, 100 * (1 - alpha/2))
        
        return (lower, upper)
    
    def _bayesian_credible_interval(self, n_mutations: int, callable_sites: int,
                                     generations: int, alpha: float = 0.05) -> Tuple[float, float]:
        """
        Calculate Bayesian credible interval using Gamma-Poisson conjugate.
        
        Prior: Gamma(a=0.5, b=0.0001) - weakly informative
        Posterior: Gamma(a + n, b + exposure)
        """
        if n_mutations == 0:
            return (0.0, 0.0)
        
        # Prior parameters (weakly informative)
        prior_a = 0.5
        prior_b = 0.0001
        
        # Exposure (genome-generations)
        exposure = callable_sites * generations
        
        # Posterior parameters
        post_a = prior_a + n_mutations
        post_b = prior_b + exposure
        
        # Credible interval from Gamma posterior
        lower = stats.gamma.ppf(alpha/2, post_a, scale=1/post_b) * 1000
        upper = stats.gamma.ppf(1-alpha/2, post_a, scale=1/post_b) * 1000
        
        return (lower, upper)
    
    def bootstrap_ci(self, n_mutations: int, callable_sites: int,
                     generations: int, alpha: float = 0.05) -> Tuple[float, float]:
        """Public method for backward compatibility."""
        return self._poisson_ci(n_mutations, callable_sites, generations, alpha)
    
    def analyze_spectrum(self, mutations: List, by_sample: bool = True) -> Dict:
        """Analyze mutation spectrum."""
        self.logger.info("Analyzing mutation spectrum...")
        
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        def normalize_mutation(ref: str, alt: str) -> Tuple[str, str]:
            if ref in ['C', 'T']:
                return (ref, alt)
            else:
                return (complement.get(ref, ref), complement.get(alt, alt))
        
        if by_sample:
            results = {}
            mutations_by_sample = defaultdict(list)
            for mut in mutations:
                if mut.mutation_type == 'SNP':
                    mutations_by_sample[mut.sample_name].append(mut)
            
            for sample_name, sample_muts in mutations_by_sample.items():
                spectrum = Counter()
                for mut in sample_muts:
                    if len(mut.ref_allele) == 1 and len(mut.alt_allele) == 1:
                        ref, alt = normalize_mutation(mut.ref_allele, mut.alt_allele)
                        spectrum[f'{ref}>{alt}'] += 1
                
                result = SpectrumResult(
                    sample_name=sample_name,
                    C_to_A=spectrum.get('C>A', 0),
                    C_to_G=spectrum.get('C>G', 0),
                    C_to_T=spectrum.get('C>T', 0),
                    T_to_A=spectrum.get('T>A', 0),
                    T_to_C=spectrum.get('T>C', 0),
                    T_to_G=spectrum.get('T>G', 0)
                )
                
                result.transitions = result.C_to_T + result.T_to_C
                result.transversions = result.C_to_A + result.C_to_G + result.T_to_A + result.T_to_G
                result.ts_tv_ratio = (result.transitions / result.transversions 
                                      if result.transversions > 0 else 0)
                
                gc_to_at = result.C_to_T
                at_to_gc = result.T_to_C
                result.gc_to_at_bias = gc_to_at / at_to_gc if at_to_gc > 0 else float('inf')
                
                results[sample_name] = result
        else:
            spectrum = Counter()
            for mut in mutations:
                if mut.mutation_type == 'SNP' and len(mut.ref_allele) == 1 and len(mut.alt_allele) == 1:
                    ref, alt = normalize_mutation(mut.ref_allele, mut.alt_allele)
                    spectrum[f'{ref}>{alt}'] += 1
            
            results = {
                'C>A': spectrum.get('C>A', 0),
                'C>G': spectrum.get('C>G', 0),
                'C>T': spectrum.get('C>T', 0),
                'T>A': spectrum.get('T>A', 0),
                'T>C': spectrum.get('T>C', 0),
                'T>G': spectrum.get('T>G', 0)
            }
            
            ts = results['C>T'] + results['T>C']
            tv = results['C>A'] + results['C>G'] + results['T>A'] + results['T>G']
            results['transitions'] = ts
            results['transversions'] = tv
            results['ts_tv_ratio'] = ts / tv if tv > 0 else 0
        
        return results
    
    def analyze_indels(self, mutations: List) -> Dict:
        """Analyze INDEL characteristics."""
        self.logger.info("Analyzing INDELs...")
        
        insertions = [m for m in mutations if m.mutation_type == 'INS']
        deletions = [m for m in mutations if m.mutation_type == 'DEL']
        
        ins_sizes = [len(m.alt_allele) - len(m.ref_allele) for m in insertions]
        del_sizes = [len(m.ref_allele) - len(m.alt_allele) for m in deletions]
        
        # Repeat context analysis
        in_repeat_ins = sum(1 for m in insertions if hasattr(m, 'in_repeat') and m.in_repeat)
        in_repeat_del = sum(1 for m in deletions if hasattr(m, 'in_repeat') and m.in_repeat)
        
        results = {
            'total_indels': len(insertions) + len(deletions),
            'insertions': len(insertions),
            'deletions': len(deletions),
            'ins_del_ratio': len(insertions) / len(deletions) if len(deletions) > 0 else float('inf'),
            'insertion_sizes': ins_sizes,
            'deletion_sizes': del_sizes,
            'mean_insertion_size': np.mean(ins_sizes) if ins_sizes else 0,
            'mean_deletion_size': np.mean(del_sizes) if del_sizes else 0,
            'median_insertion_size': np.median(ins_sizes) if ins_sizes else 0,
            'median_deletion_size': np.median(del_sizes) if del_sizes else 0,
            'insertions_in_repeats': in_repeat_ins,
            'deletions_in_repeats': in_repeat_del,
            'repeat_fraction_ins': in_repeat_ins / len(insertions) if insertions else 0,
            'repeat_fraction_del': in_repeat_del / len(deletions) if deletions else 0,
            'size_distribution': {
                '1bp': sum(1 for s in ins_sizes + del_sizes if s == 1),
                '2-5bp': sum(1 for s in ins_sizes + del_sizes if 2 <= s <= 5),
                '6-10bp': sum(1 for s in ins_sizes + del_sizes if 6 <= s <= 10),
                '>10bp': sum(1 for s in ins_sizes + del_sizes if s > 10)
            }
        }
        
        return results
    
    def analyze_chromosome_distribution(self, mutations: List) -> Dict:
        """Analyze distribution of mutations across chromosomes."""
        self.logger.info("Analyzing chromosome distribution...")
        
        by_chrom = defaultdict(lambda: {'snps': 0, 'indels': 0, 'total': 0})
        
        for mut in mutations:
            by_chrom[mut.chromosome]['total'] += 1
            if mut.mutation_type == 'SNP':
                by_chrom[mut.chromosome]['snps'] += 1
            else:
                by_chrom[mut.chromosome]['indels'] += 1
        
        return dict(by_chrom)
    
    def poisson_test_uniformity(self, mutations: List, samples: List) -> Dict:
        """Test if mutation counts follow Poisson distribution."""
        self.logger.info("Testing Poisson uniformity...")
        
        counts = []
        for sample in samples:
            n = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            counts.append(n)
        
        if not counts or len(counts) < 2:
            return {'test': 'poisson', 'p_value': 1.0, 'conclusion': 'insufficient data'}
        
        mean_count = np.mean(counts)
        variance = np.var(counts, ddof=1)
        
        dispersion = variance / mean_count if mean_count > 0 else 0
        
        n = len(counts)
        chi2_stat = (n - 1) * variance / mean_count if mean_count > 0 else 0
        p_value = 1 - stats.chi2.cdf(chi2_stat, n - 1)
        
        conclusion = 'uniform' if p_value > 0.05 else 'non-uniform (possible selection or hotspots)'
        
        return {
            'test': 'poisson_dispersion',
            'mean_mutations': float(mean_count),
            'variance': float(variance),
            'dispersion_index': float(dispersion),
            'chi2_statistic': float(chi2_stat),
            'degrees_of_freedom': n - 1,
            'p_value': float(p_value),
            'conclusion': conclusion
        }
    
    # ============ NEW STATISTICAL TESTS ============
    
    def kolmogorov_smirnov_test(self, mutations: List, samples: List) -> DistributionTestResult:
        """
        Kolmogorov-Smirnov test comparing mutation counts to Poisson distribution.
        """
        self.logger.info("Performing Kolmogorov-Smirnov test...")
        
        counts = [sum(1 for m in mutations if m.sample_name == s.sample_name) for s in samples]
        
        if not counts:
            return DistributionTestResult(
                test_name='kolmogorov_smirnov',
                statistic=0.0,
                p_value=1.0,
                conclusion='insufficient data'
            )
        
        mean_count = np.mean(counts)
        
        # Compare to Poisson distribution with same mean
        statistic, p_value = stats.kstest(counts, 'poisson', args=(mean_count,))
        
        conclusion = 'consistent with Poisson' if p_value > 0.05 else 'deviates from Poisson'
        
        return DistributionTestResult(
            test_name='kolmogorov_smirnov_vs_poisson',
            statistic=float(statistic),
            p_value=float(p_value),
            conclusion=conclusion
        )
    
    def shapiro_wilk_test(self, mutations: List, samples: List) -> NormalityTestResult:
        """
        Shapiro-Wilk test for normality of mutation counts.
        """
        self.logger.info("Performing Shapiro-Wilk normality test...")
        
        counts = [sum(1 for m in mutations if m.sample_name == s.sample_name) for s in samples]
        
        if len(counts) < 3:
            return NormalityTestResult(
                test_name='shapiro_wilk',
                statistic=0.0,
                p_value=1.0,
                is_normal=True
            )
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            statistic, p_value = stats.shapiro(counts)
        
        return NormalityTestResult(
            test_name='shapiro_wilk',
            statistic=float(statistic),
            p_value=float(p_value),
            is_normal=p_value > 0.05
        )
    
    def permutation_test_genotypes(self, mutations: List, samples: List,
                                    n_permutations: int = None) -> PermutationTestResult:
        """
        Permutation test for difference between genotypes.
        """
        self.logger.info("Performing permutation test for genotype comparison...")
        
        if n_permutations is None:
            n_permutations = self.permutation_n
        
        # Group by genotype
        by_genotype = defaultdict(list)
        for sample in samples:
            n_muts = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            by_genotype[sample.genotype].append(n_muts)
        
        genotypes = list(by_genotype.keys())
        
        if len(genotypes) < 2:
            return PermutationTestResult(
                test_name='permutation_genotype',
                observed_statistic=0.0,
                p_value=1.0,
                n_permutations=0,
                null_distribution_mean=0.0,
                null_distribution_std=0.0
            )
        
        # Use first two genotypes
        g1, g2 = genotypes[0], genotypes[1]
        group1 = by_genotype[g1]
        group2 = by_genotype[g2]
        
        # Observed difference in means
        observed_diff = np.mean(group1) - np.mean(group2)
        
        # Pool all values
        pooled = group1 + group2
        n1 = len(group1)
        
        # Permutation test
        null_diffs = []
        for _ in range(n_permutations):
            np.random.shuffle(pooled)
            perm_diff = np.mean(pooled[:n1]) - np.mean(pooled[n1:])
            null_diffs.append(perm_diff)
        
        # Two-tailed p-value
        p_value = np.mean(np.abs(null_diffs) >= np.abs(observed_diff))
        
        return PermutationTestResult(
            test_name=f'permutation_{g1}_vs_{g2}',
            observed_statistic=float(observed_diff),
            p_value=float(p_value),
            n_permutations=n_permutations,
            null_distribution_mean=float(np.mean(null_diffs)),
            null_distribution_std=float(np.std(null_diffs))
        )
    
    def calculate_effect_sizes(self, mutations: List, samples: List) -> List[EffectSizeResult]:
        """
        Calculate effect sizes (Cohen's d, Hedges' g, Glass's delta) for genotype comparisons.
        """
        self.logger.info("Calculating effect sizes...")
        
        by_genotype = defaultdict(list)
        for sample in samples:
            n_muts = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            by_genotype[sample.genotype].append(n_muts)
        
        genotypes = list(by_genotype.keys())
        results = []
        
        for i, g1 in enumerate(genotypes):
            for g2 in genotypes[i+1:]:
                group1 = np.array(by_genotype[g1])
                group2 = np.array(by_genotype[g2])
                
                if len(group1) < 2 or len(group2) < 2:
                    continue
                
                # Cohen's d
                pooled_std = np.sqrt(((len(group1)-1)*np.var(group1, ddof=1) + 
                                      (len(group2)-1)*np.var(group2, ddof=1)) / 
                                     (len(group1) + len(group2) - 2))
                
                cohens_d = (np.mean(group1) - np.mean(group2)) / pooled_std if pooled_std > 0 else 0
                
                # Hedges' g (corrected for small sample sizes)
                correction = 1 - (3 / (4 * (len(group1) + len(group2)) - 9))
                hedges_g = cohens_d * correction
                
                # Glass's delta (using control group std)
                glass_delta = (np.mean(group1) - np.mean(group2)) / np.std(group2, ddof=1) if np.std(group2, ddof=1) > 0 else 0
                
                # Bootstrap CI for Cohen's d
                boot_ds = []
                for _ in range(1000):
                    boot1 = np.random.choice(group1, size=len(group1), replace=True)
                    boot2 = np.random.choice(group2, size=len(group2), replace=True)
                    boot_pooled_std = np.sqrt(((len(boot1)-1)*np.var(boot1, ddof=1) + 
                                               (len(boot2)-1)*np.var(boot2, ddof=1)) / 
                                              (len(boot1) + len(boot2) - 2))
                    if boot_pooled_std > 0:
                        boot_ds.append((np.mean(boot1) - np.mean(boot2)) / boot_pooled_std)
                
                ci_lower = np.percentile(boot_ds, 2.5) if boot_ds else 0
                ci_upper = np.percentile(boot_ds, 97.5) if boot_ds else 0
                
                # Interpretation
                abs_d = abs(cohens_d)
                if abs_d < 0.2:
                    interpretation = 'negligible'
                elif abs_d < 0.5:
                    interpretation = 'small'
                elif abs_d < 0.8:
                    interpretation = 'medium'
                else:
                    interpretation = 'large'
                
                results.append(EffectSizeResult(
                    genotype1=g1,
                    genotype2=g2,
                    cohens_d=float(cohens_d),
                    cohens_d_ci_lower=float(ci_lower),
                    cohens_d_ci_upper=float(ci_upper),
                    glass_delta=float(glass_delta),
                    hedges_g=float(hedges_g),
                    interpretation=interpretation
                ))
        
        return results
    
    def bootstrap_spectrum_comparison(self, mutations: List, samples: List,
                                       n_bootstrap: int = None) -> Dict:
        """
        Bootstrap analysis for comparing mutation spectra between genotypes.
        """
        self.logger.info("Bootstrap spectrum comparison...")
        
        if n_bootstrap is None:
            n_bootstrap = self.bootstrap_n
        
        # Get spectrum by genotype
        by_genotype = defaultdict(list)
        for sample in samples:
            sample_muts = [m for m in mutations if m.sample_name == sample.sample_name and m.mutation_type == 'SNP']
            by_genotype[sample.genotype].extend(sample_muts)
        
        genotypes = list(by_genotype.keys())
        
        if len(genotypes) < 2:
            return {'comparison': 'single_genotype'}
        
        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        def get_spectrum(muts):
            spectrum = Counter()
            for m in muts:
                if len(m.ref_allele) == 1 and len(m.alt_allele) == 1:
                    ref, alt = m.ref_allele, m.alt_allele
                    if ref in ['G', 'A']:
                        ref = complement.get(ref, ref)
                        alt = complement.get(alt, alt)
                    spectrum[f'{ref}>{alt}'] += 1
            return [spectrum.get(c, 0) for c in categories]
        
        results = {'comparisons': []}
        
        for i, g1 in enumerate(genotypes):
            for g2 in genotypes[i+1:]:
                muts1 = by_genotype[g1]
                muts2 = by_genotype[g2]
                
                spec1 = get_spectrum(muts1)
                spec2 = get_spectrum(muts2)
                
                total1 = sum(spec1)
                total2 = sum(spec2)
                
                if total1 == 0 or total2 == 0:
                    continue
                
                # Observed chi-square
                freq1 = np.array(spec1) / total1
                freq2 = np.array(spec2) / total2
                
                # Bootstrap difference
                boot_diffs = []
                all_muts = muts1 + muts2
                n1 = len(muts1)
                
                for _ in range(n_bootstrap):
                    boot_muts = np.random.choice(all_muts, size=len(all_muts), replace=True).tolist()
                    boot_spec1 = get_spectrum(boot_muts[:n1])
                    boot_spec2 = get_spectrum(boot_muts[n1:])
                    
                    boot_total1 = sum(boot_spec1)
                    boot_total2 = sum(boot_spec2)
                    
                    if boot_total1 > 0 and boot_total2 > 0:
                        boot_freq1 = np.array(boot_spec1) / boot_total1
                        boot_freq2 = np.array(boot_spec2) / boot_total2
                        boot_diffs.append(np.sum(np.abs(boot_freq1 - boot_freq2)))
                
                observed_diff = np.sum(np.abs(freq1 - freq2))
                
                # P-value from bootstrap distribution
                p_value = np.mean(np.array(boot_diffs) >= observed_diff) if boot_diffs else 1.0
                
                results['comparisons'].append({
                    'genotype1': g1,
                    'genotype2': g2,
                    'observed_difference': float(observed_diff),
                    'p_value': float(p_value),
                    'spectrum1': dict(zip(categories, spec1)),
                    'spectrum2': dict(zip(categories, spec2)),
                    'significant': p_value < 0.05
                })
        
        return results
    
    def kruskal_wallis_test(self, mutations: List, samples: List) -> Dict:
        """
        Kruskal-Wallis H-test for comparing more than 2 genotypes.
        """
        self.logger.info("Performing Kruskal-Wallis test...")
        
        by_genotype = defaultdict(list)
        for sample in samples:
            n_muts = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            by_genotype[sample.genotype].append(n_muts)
        
        groups = [v for v in by_genotype.values() if len(v) >= 2]
        
        if len(groups) < 2:
            return {
                'test': 'kruskal_wallis',
                'statistic': 0.0,
                'p_value': 1.0,
                'conclusion': 'insufficient groups'
            }
        
        statistic, p_value = stats.kruskal(*groups)
        
        return {
            'test': 'kruskal_wallis',
            'statistic': float(statistic),
            'p_value': float(p_value),
            'n_groups': len(groups),
            'conclusion': 'significant difference' if p_value < 0.05 else 'no significant difference'
        }
    
    def compare_genotypes(self, mutations: List, samples: List) -> Dict:
        """Enhanced genotype comparison with multiple tests."""
        self.logger.info("Comparing genotypes...")
        
        by_genotype = defaultdict(list)
        for sample in samples:
            n_muts = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            by_genotype[sample.genotype].append(n_muts)
        
        genotypes = list(by_genotype.keys())
        
        if len(genotypes) < 2:
            return {'comparison': 'single_genotype', 'genotypes': genotypes}
        
        results = {
            'genotypes': {},
            'pairwise_comparisons': [],
            'effect_sizes': [],
            'overall_test': None
        }
        
        for geno, counts in by_genotype.items():
            results['genotypes'][geno] = {
                'n_samples': len(counts),
                'mean_mutations': float(np.mean(counts)),
                'std_mutations': float(np.std(counts)),
                'median_mutations': float(np.median(counts)),
                'total_mutations': int(sum(counts))
            }
        
        # Pairwise comparisons
        for i, g1 in enumerate(genotypes):
            for g2 in genotypes[i+1:]:
                counts1 = by_genotype[g1]
                counts2 = by_genotype[g2]
                
                if len(counts1) >= 2 and len(counts2) >= 2:
                    # Mann-Whitney U test
                    u_stat, mw_p = stats.mannwhitneyu(counts1, counts2, alternative='two-sided')
                    
                    # Welch's t-test
                    t_stat, t_p = stats.ttest_ind(counts1, counts2, equal_var=False)
                    
                    results['pairwise_comparisons'].append({
                        'genotype1': g1,
                        'genotype2': g2,
                        'mann_whitney_u': float(u_stat),
                        'mann_whitney_p': float(mw_p),
                        'welch_t': float(t_stat),
                        'welch_p': float(t_p),
                        'significant': mw_p < 0.05
                    })
        
        # Overall test (Kruskal-Wallis)
        results['overall_test'] = self.kruskal_wallis_test(mutations, samples)
        
        # Effect sizes
        effect_sizes = self.calculate_effect_sizes(mutations, samples)
        results['effect_sizes'] = [asdict(es) for es in effect_sizes]
        
        return results
    
    def write_results(self, rates: List[MutationRateResult], spectrum: Dict,
                      indels: Dict, chrom_dist: Dict, poisson: Dict,
                      genotype_comparison: Dict, summary: Dict,
                      additional_tests: Dict = None):
        """Write all statistical results to files."""
        
        # Mutation rates TSV (enhanced)
        with open(self.output_dir / 'mutation_rates.tsv', 'w') as f:
            headers = ['sample_name', 'genotype', 'generation', 'total_mutations',
                      'snp_count', 'indel_count', 'insertion_count', 'deletion_count',
                      'callable_sites', 'rate_per_site_per_gen', 'rate_per_kb_per_gen',
                      'ci_lower', 'ci_upper', 'bootstrap_ci_lower', 'bootstrap_ci_upper',
                      'bayesian_ci_lower', 'bayesian_ci_upper']
            f.write('\t'.join(headers) + '\n')
            
            for r in rates:
                values = [r.sample_name, r.genotype, str(r.generation), str(r.total_mutations),
                         str(r.snp_count), str(r.indel_count), str(r.insertion_count),
                         str(r.deletion_count), str(r.callable_sites),
                         f'{r.rate_per_site_per_gen:.2e}', f'{r.rate_per_kb_per_gen:.2e}',
                         f'{r.ci_lower:.2e}', f'{r.ci_upper:.2e}',
                         f'{r.bootstrap_ci_lower:.2e}', f'{r.bootstrap_ci_upper:.2e}',
                         f'{r.bayesian_ci_lower:.2e}', f'{r.bayesian_ci_upper:.2e}']
                f.write('\t'.join(values) + '\n')
        
        # Spectrum TSV
        with open(self.output_dir / 'mutation_spectrum.tsv', 'w') as f:
            if spectrum and isinstance(list(spectrum.values())[0], SpectrumResult):
                headers = ['sample_name', 'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G',
                          'transitions', 'transversions', 'ts_tv_ratio', 'gc_to_at_bias']
                f.write('\t'.join(headers) + '\n')
                
                for name, s in spectrum.items():
                    values = [name, str(s.C_to_A), str(s.C_to_G), str(s.C_to_T),
                             str(s.T_to_A), str(s.T_to_C), str(s.T_to_G),
                             str(s.transitions), str(s.transversions),
                             f'{s.ts_tv_ratio:.3f}', 
                             f'{s.gc_to_at_bias:.3f}' if s.gc_to_at_bias != float('inf') else 'inf']
                    f.write('\t'.join(values) + '\n')
        
        # INDEL analysis TSV
        with open(self.output_dir / 'indel_analysis.tsv', 'w') as f:
            f.write('metric\tvalue\n')
            for key, value in indels.items():
                if not isinstance(value, (list, dict)):
                    f.write(f'{key}\t{value}\n')
        
        # Chromosome distribution TSV
        with open(self.output_dir / 'chromosome_distribution.tsv', 'w') as f:
            f.write('chromosome\tsnps\tindels\ttotal\n')
            for chrom, counts in sorted(chrom_dist.items()):
                f.write(f"{chrom}\t{counts['snps']}\t{counts['indels']}\t{counts['total']}\n")
        
        # Effect sizes TSV
        if genotype_comparison.get('effect_sizes'):
            with open(self.output_dir / 'effect_sizes.tsv', 'w') as f:
                headers = ['genotype1', 'genotype2', 'cohens_d', 'cohens_d_ci_lower',
                          'cohens_d_ci_upper', 'hedges_g', 'glass_delta', 'interpretation']
                f.write('\t'.join(headers) + '\n')
                for es in genotype_comparison['effect_sizes']:
                    values = [es['genotype1'], es['genotype2'],
                             f"{es['cohens_d']:.4f}", f"{es['cohens_d_ci_lower']:.4f}",
                             f"{es['cohens_d_ci_upper']:.4f}", f"{es['hedges_g']:.4f}",
                             f"{es['glass_delta']:.4f}", es['interpretation']]
                    f.write('\t'.join(values) + '\n')
        
        # Summary JSON (enhanced)
        summary_data = {
            'summary': summary,
            'poisson_test': poisson,
            'genotype_comparison': {k: v for k, v in genotype_comparison.items() 
                                    if k != 'effect_sizes'},
            'indel_analysis': {k: v for k, v in indels.items() if not isinstance(v, list)},
        }
        
        if additional_tests:
            summary_data['additional_tests'] = {
                'kolmogorov_smirnov': asdict(additional_tests.get('ks_test')) if additional_tests.get('ks_test') else None,
                'shapiro_wilk': asdict(additional_tests.get('shapiro_test')) if additional_tests.get('shapiro_test') else None,
                'permutation_test': asdict(additional_tests.get('permutation_test')) if additional_tests.get('permutation_test') else None,
                'kruskal_wallis': additional_tests.get('kruskal_wallis'),
                'bootstrap_spectrum': additional_tests.get('bootstrap_spectrum')
            }
        
        with open(self.output_dir / 'summary_statistics.json', 'w') as f:
            json.dump(summary_data, f, indent=2, default=str)
        
        # Statistical tests TSV
        with open(self.output_dir / 'statistical_tests.tsv', 'w') as f:
            f.write('test_name\tstatistic\tp_value\tconclusion\n')
            
            f.write(f"poisson_dispersion\t{poisson.get('chi2_statistic', 'NA')}\t{poisson.get('p_value', 'NA')}\t{poisson.get('conclusion', 'NA')}\n")
            
            if additional_tests:
                if additional_tests.get('ks_test'):
                    ks = additional_tests['ks_test']
                    f.write(f"{ks.test_name}\t{ks.statistic:.4f}\t{ks.p_value:.4f}\t{ks.conclusion}\n")
                
                if additional_tests.get('shapiro_test'):
                    sw = additional_tests['shapiro_test']
                    f.write(f"{sw.test_name}\t{sw.statistic:.4f}\t{sw.p_value:.4f}\t{'normal' if sw.is_normal else 'non-normal'}\n")
                
                if additional_tests.get('permutation_test'):
                    pt = additional_tests['permutation_test']
                    f.write(f"{pt.test_name}\t{pt.observed_statistic:.4f}\t{pt.p_value:.4f}\tpermutation_n={pt.n_permutations}\n")
        
        self.logger.info(f"Statistical results written to {self.output_dir}")
    
    def run(self, mutations: List, samples: List, ancestor) -> Dict:
        """Run all statistical analyses."""
        self.logger.info("Running comprehensive statistical analysis...")
        
        ma_lines = [s for s in samples if s.sample_name != ancestor.sample_name]
        
        # Basic analyses
        rates = self.calculate_mutation_rates(mutations, samples, ancestor)
        spectrum = self.analyze_spectrum(mutations, by_sample=True)
        overall_spectrum = self.analyze_spectrum(mutations, by_sample=False)
        indels = self.analyze_indels(mutations)
        chrom_dist = self.analyze_chromosome_distribution(mutations)
        poisson = self.poisson_test_uniformity(mutations, ma_lines)
        genotype_comp = self.compare_genotypes(mutations, ma_lines)
        
        # Additional tests
        additional_tests = {}
        
        if self.config.get('statistics', {}).get('additional_tests', {}).get('kolmogorov_smirnov', True):
            additional_tests['ks_test'] = self.kolmogorov_smirnov_test(mutations, ma_lines)
        
        if self.config.get('statistics', {}).get('additional_tests', {}).get('shapiro_wilk', True):
            additional_tests['shapiro_test'] = self.shapiro_wilk_test(mutations, ma_lines)
        
        if self.config.get('statistics', {}).get('additional_tests', {}).get('permutation_tests', True):
            additional_tests['permutation_test'] = self.permutation_test_genotypes(mutations, ma_lines)
        
        additional_tests['kruskal_wallis'] = self.kruskal_wallis_test(mutations, ma_lines)
        additional_tests['bootstrap_spectrum'] = self.bootstrap_spectrum_comparison(mutations, ma_lines)
        
        # Summary
        total_muts = len(mutations)
        total_snps = sum(1 for m in mutations if m.mutation_type == 'SNP')
        total_indels = total_muts - total_snps
        
        mean_rate = np.mean([r.rate_per_site_per_gen for r in rates]) if rates else 0
        
        summary = {
            'total_samples': len(ma_lines),
            'total_mutations': total_muts,
            'total_snps': total_snps,
            'total_indels': total_indels,
            'snp_fraction': total_snps / total_muts if total_muts > 0 else 0,
            'indel_fraction': total_indels / total_muts if total_muts > 0 else 0,
            'mean_mutations_per_sample': total_muts / len(ma_lines) if ma_lines else 0,
            'mean_mutation_rate_per_site_per_gen': float(mean_rate),
            'overall_ts_tv_ratio': overall_spectrum.get('ts_tv_ratio', 0),
            'generations': self.generations
        }
        
        # Store results for plotting
        self.results_cache = {
            'rates': rates,
            'spectrum': spectrum,
            'overall_spectrum': overall_spectrum,
            'indels': indels,
            'chromosome_distribution': chrom_dist,
            'poisson_test': poisson,
            'genotype_comparison': genotype_comp,
            'additional_tests': additional_tests,
            'summary': summary
        }
        
        # Write results
        self.write_results(rates, spectrum, indels, chrom_dist, poisson, 
                          genotype_comp, summary, additional_tests)
        
        return self.results_cache


# Backward compatibility
ComprehensiveStatistics = EnhancedStatistics
