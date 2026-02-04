#!/usr/bin/env python3
"""
Group Plots for Fusarium MA Pipeline

Generates visualizations for ancestor and generation groups:
- Mutation rate comparison by genotype
- Spectrum comparison by genotype
- Mutation accumulation curves
- Ancestor vs MA lines comparison
- Genotype effect size forest plots
- Replicate consistency plots
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional
from collections import defaultdict, Counter

HAS_PLOTTING = False
HAS_SEABORN = False

try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_PLOTTING = True
    try:
        import seaborn as sns
        HAS_SEABORN = True
    except ImportError:
        pass
except ImportError as e:
    logging.warning(f"Matplotlib/numpy not available: {e}")


class GroupPlotGenerator:
    """
    Generates group comparison plots for MA analysis.
    
    Creates visualizations comparing ancestor vs MA lines,
    different genotypes, and generation groups.
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.GroupPlots')
        self.output_dir = Path(config['output']['directory']) / 'plots' / 'groups'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        viz_config = config.get('visualization', {})
        self.dpi = viz_config.get('dpi', 300)
        
        if HAS_PLOTTING:
            try:
                plt.style.use('seaborn-v0_8-whitegrid')
            except:
                plt.style.use('seaborn-whitegrid')
            plt.rcParams['figure.figsize'] = (10, 6)
            plt.rcParams['font.size'] = 12
        
        self.spectrum_colors = {
            'C>A': '#1E88E5', 'C>G': '#000000', 'C>T': '#D32F2F',
            'T>A': '#757575', 'T>C': '#43A047', 'T>G': '#F9A825'
        }
    
    def generate_all_plots(self, mutations: List, samples: List,
                           rates: List, stats_results: Dict = None) -> List[str]:
        """Generate all group plots."""
        if not HAS_PLOTTING:
            self.logger.warning("Matplotlib not available, skipping plots")
            return []
        
        self.logger.info("Generating group plots...")
        
        plot_paths = []
        
        # Genotype comparison
        path = self.plot_genotype_comparison(mutations, samples, rates)
        if path:
            plot_paths.append(path)
        
        # Spectrum by genotype
        path = self.plot_spectrum_by_genotype(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # Mutation accumulation curve (if multiple generations)
        path = self.plot_accumulation_curve(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # Ancestor vs MA lines
        path = self.plot_ancestor_vs_ma(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # Effect size forest plot
        if stats_results and 'genotype_comparison' in stats_results:
            path = self.plot_effect_size_forest(stats_results['genotype_comparison'])
            if path:
                plot_paths.append(path)
        
        # Replicate consistency
        path = self.plot_replicate_consistency(mutations, samples, rates)
        if path:
            plot_paths.append(path)
        
        self.logger.info(f"Generated {len(plot_paths)} group plots")
        return plot_paths
    
    def plot_genotype_comparison(self, mutations: List, samples: List,
                                  rates: List) -> Optional[str]:
        """Create grouped bar chart comparing genotypes."""
        # Group rates by genotype
        by_genotype = defaultdict(list)
        for rate in rates:
            by_genotype[rate.genotype].append(rate.rate_per_kb_per_gen)
        
        if len(by_genotype) < 1:
            return None
        
        genotypes = list(by_genotype.keys())
        means = [np.mean(by_genotype[g]) for g in genotypes]
        stds = [np.std(by_genotype[g]) for g in genotypes]
        
        fig, ax = plt.subplots(figsize=(max(8, len(genotypes) * 2), 6))
        
        x = np.arange(len(genotypes))
        colors = plt.cm.Set2(np.linspace(0, 1, len(genotypes)))
        
        bars = ax.bar(x, means, yerr=stds, capsize=5, color=colors,
                     edgecolor='black', linewidth=1, error_kw={'elinewidth': 2})
        
        # Add sample counts
        for i, g in enumerate(genotypes):
            n = len(by_genotype[g])
            ax.text(i, means[i] + stds[i] + 0.001, f'n={n}', 
                   ha='center', va='bottom', fontsize=10)
        
        ax.set_xlabel('Genotype', fontsize=12)
        ax.set_ylabel('Mutation Rate (per kb per generation)', fontsize=12)
        ax.set_title('Mutation Rate Comparison by Genotype', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(genotypes, rotation=45, ha='right')
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'genotype_comparison.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_spectrum_by_genotype(self, mutations: List, samples: List) -> Optional[str]:
        """Create stacked/grouped spectrum comparison by genotype."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        
        # Get sample genotypes
        sample_genotype = {s.sample_name: s.genotype for s in samples}
        
        # Group mutations by genotype
        by_genotype = defaultdict(list)
        for mut in mutations:
            if mut.mutation_type == 'SNP' and len(mut.ref_allele) == 1 and len(mut.alt_allele) == 1:
                genotype = sample_genotype.get(mut.sample_name, 'Unknown')
                by_genotype[genotype].append(mut)
        
        if len(by_genotype) < 1:
            return None
        
        # Calculate spectra
        genotypes = list(by_genotype.keys())
        spectra = {}
        
        for geno, muts in by_genotype.items():
            spectrum = Counter()
            for mut in muts:
                ref, alt = mut.ref_allele, mut.alt_allele
                if ref in ['G', 'A']:
                    ref = complement.get(ref, ref)
                    alt = complement.get(alt, alt)
                spectrum[f'{ref}>{alt}'] += 1
            
            total = sum(spectrum.values())
            spectra[geno] = [spectrum.get(c, 0) / total if total > 0 else 0 for c in categories]
        
        # Create grouped bar plot
        fig, ax = plt.subplots(figsize=(12, 6))
        
        x = np.arange(len(categories))
        width = 0.8 / len(genotypes)
        colors = plt.cm.Set2(np.linspace(0, 1, len(genotypes)))
        
        for i, geno in enumerate(genotypes):
            offset = (i - len(genotypes)/2 + 0.5) * width
            ax.bar(x + offset, spectra[geno], width, label=geno, color=colors[i],
                  edgecolor='black', linewidth=0.5)
        
        ax.set_xlabel('Mutation Type', fontsize=12)
        ax.set_ylabel('Fraction', fontsize=12)
        ax.set_title('Mutation Spectrum Comparison by Genotype', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend(title='Genotype')
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'spectrum_by_genotype.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_accumulation_curve(self, mutations: List, samples: List) -> Optional[str]:
        """Create mutation accumulation curve by generation."""
        # Group by generation
        by_generation = defaultdict(list)
        sample_gen = {s.sample_name: s.generation for s in samples}
        
        for mut in mutations:
            gen = sample_gen.get(mut.sample_name, 0)
            by_generation[gen].append(mut)
        
        if len(by_generation) < 2:
            return None
        
        generations = sorted(by_generation.keys())
        counts = [len(by_generation[g]) for g in generations]
        cumulative = np.cumsum(counts)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Per-generation
        ax1.bar(generations, counts, color='steelblue', edgecolor='black')
        ax1.set_xlabel('Generation', fontsize=12)
        ax1.set_ylabel('Mutation Count', fontsize=12)
        ax1.set_title('Mutations per Generation', fontsize=12, fontweight='bold')
        
        # Cumulative
        ax2.plot(generations, cumulative, 'o-', color='steelblue', 
                linewidth=2, markersize=8)
        ax2.fill_between(generations, cumulative, alpha=0.3)
        ax2.set_xlabel('Generation', fontsize=12)
        ax2.set_ylabel('Cumulative Mutations', fontsize=12)
        ax2.set_title('Cumulative Mutation Accumulation', fontsize=12, fontweight='bold')
        
        plt.suptitle('Mutation Accumulation Over Generations', 
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'accumulation_curve.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_ancestor_vs_ma(self, mutations: List, samples: List) -> Optional[str]:
        """Create side-by-side comparison of ancestor vs MA lines."""
        # Identify ancestor
        ancestors = [s for s in samples if s.sample_type == 'ancestor']
        ma_lines = [s for s in samples if s.sample_type == 'ma_line']
        
        if not ancestors or not ma_lines:
            return None
        
        # Count mutations per sample type
        ancestor_muts = [m for m in mutations if any(m.sample_name == a.sample_name for a in ancestors)]
        ma_muts = [m for m in mutations if any(m.sample_name == ml.sample_name for ml in ma_lines)]
        
        # Create comparison
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # Panel 1: Total counts
        ax1 = axes[0]
        counts = [len(ancestor_muts), len(ma_muts)]
        labels = ['Ancestor', f'MA Lines (n={len(ma_lines)})']
        ax1.bar(labels, counts, color=['gray', 'steelblue'], edgecolor='black')
        ax1.set_ylabel('Mutation Count')
        ax1.set_title('Total Mutations', fontweight='bold')
        
        # Panel 2: Mutation types
        ax2 = axes[1]
        types = ['SNP', 'INS', 'DEL']
        anc_types = [sum(1 for m in ancestor_muts if m.mutation_type == t) for t in types]
        ma_types = [sum(1 for m in ma_muts if m.mutation_type == t) / len(ma_lines) if ma_lines else 0 for t in types]
        
        x = np.arange(len(types))
        width = 0.35
        ax2.bar(x - width/2, anc_types, width, label='Ancestor', color='gray')
        ax2.bar(x + width/2, ma_types, width, label='MA Lines (avg)', color='steelblue')
        ax2.set_xticks(x)
        ax2.set_xticklabels(types)
        ax2.set_ylabel('Count')
        ax2.set_title('Mutation Types', fontweight='bold')
        ax2.legend()
        
        # Panel 3: Quality comparison
        ax3 = axes[2]
        anc_qual = [m.quality for m in ancestor_muts] or [0]
        ma_qual = [m.quality for m in ma_muts] or [0]
        
        bp = ax3.boxplot([anc_qual, ma_qual], labels=['Ancestor', 'MA Lines'],
                        patch_artist=True)
        bp['boxes'][0].set_facecolor('gray')
        bp['boxes'][1].set_facecolor('steelblue')
        ax3.set_ylabel('Variant Quality')
        ax3.set_title('Quality Distribution', fontweight='bold')
        
        plt.suptitle('Ancestor vs MA Lines Comparison', fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'ancestor_vs_ma.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_effect_size_forest(self, genotype_comparison: Dict) -> Optional[str]:
        """Create forest plot of effect sizes."""
        effect_sizes = genotype_comparison.get('effect_sizes', [])
        
        if not effect_sizes:
            return None
        
        fig, ax = plt.subplots(figsize=(10, max(4, len(effect_sizes) * 0.5)))
        
        y_pos = np.arange(len(effect_sizes))
        
        for i, es in enumerate(effect_sizes):
            # Point estimate
            ax.scatter(es['cohens_d'], i, s=100, color='steelblue', zorder=3)
            
            # Confidence interval
            ax.hlines(i, es['cohens_d_ci_lower'], es['cohens_d_ci_upper'], 
                     color='steelblue', linewidth=2, zorder=2)
        
        # Reference line at 0
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, zorder=1)
        
        # Effect size guidelines
        for x, label in [(-0.8, 'Large'), (-0.5, 'Medium'), (-0.2, 'Small'),
                         (0.2, 'Small'), (0.5, 'Medium'), (0.8, 'Large')]:
            ax.axvline(x=x, color='lightgray', linestyle=':', linewidth=0.5, zorder=0)
        
        # Labels
        labels = [f"{es['genotype1']} vs {es['genotype2']}" for es in effect_sizes]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        
        ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12)
        ax.set_title('Genotype Effect Sizes\n(with 95% CI)', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'effect_size_forest.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_replicate_consistency(self, mutations: List, samples: List,
                                    rates: List) -> Optional[str]:
        """Create replicate consistency scatter plot."""
        # Group by genotype and replicate
        by_genotype = defaultdict(list)
        for rate in rates:
            by_genotype[rate.genotype].append((rate.sample_name, rate.rate_per_kb_per_gen))
        
        if len(by_genotype) < 1 or all(len(v) < 2 for v in by_genotype.values()):
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = plt.cm.Set2(np.linspace(0, 1, len(by_genotype)))
        
        for i, (geno, sample_rates) in enumerate(by_genotype.items()):
            if len(sample_rates) < 2:
                continue
            
            rates_only = [r[1] for r in sample_rates]
            mean_rate = np.mean(rates_only)
            
            # Plot individual points
            x_jitter = np.random.normal(i, 0.1, len(rates_only))
            ax.scatter(x_jitter, rates_only, c=[colors[i]], s=50, alpha=0.7)
            
            # Plot mean
            ax.scatter(i, mean_rate, c='black', s=100, marker='_', linewidths=3, zorder=5)
        
        ax.set_xlabel('Genotype', fontsize=12)
        ax.set_ylabel('Mutation Rate (per kb per generation)', fontsize=12)
        ax.set_title('Replicate Consistency by Genotype', fontsize=14, fontweight='bold')
        ax.set_xticks(range(len(by_genotype)))
        ax.set_xticklabels(list(by_genotype.keys()), rotation=45, ha='right')
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'replicate_consistency.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
