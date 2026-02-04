#!/usr/bin/env python3
"""
Individual Sample Plots for Fusarium MA Pipeline

Generates comprehensive visualizations for each sample:
- Mutation spectrum (6-category)
- Depth distribution histogram
- Strand bias distribution
- Allele frequency distribution
- Chromosome distribution
- Quality metrics dashboard
- Trinucleotide context (96-category)
- Variant quality scatter
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional
from collections import Counter, defaultdict

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


class SamplePlotGenerator:
    """
    Generates individual sample plots for mutation analysis.
    
    Creates a comprehensive set of visualizations for each sample,
    including mutation profiles, quality metrics, and genomic context.
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.SamplePlots')
        self.output_dir = Path(config['output']['directory']) / 'plots' / 'individual'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Visualization settings
        viz_config = config.get('visualization', {})
        self.dpi = viz_config.get('dpi', 300)
        self.formats = viz_config.get('individual_plots', {}).get('formats', ['png'])
        
        # Set style
        if HAS_PLOTTING:
            try:
                plt.style.use('seaborn-v0_8-whitegrid')
            except:
                plt.style.use('seaborn-whitegrid')
            plt.rcParams['figure.figsize'] = (10, 6)
            plt.rcParams['font.size'] = 12
        
        # Standard mutation spectrum colors
        self.spectrum_colors = {
            'C>A': '#1E88E5',
            'C>G': '#000000',
            'C>T': '#D32F2F',
            'T>A': '#757575',
            'T>C': '#43A047',
            'T>G': '#F9A825'
        }
    
    def generate_all_plots(self, mutations: List, samples: List, 
                           muver_data: Dict = None) -> Dict[str, List[str]]:
        """
        Generate all individual sample plots.
        
        Parameters
        ----------
        mutations : list
            List of Mutation objects
        samples : list
            List of Sample objects
        muver_data : dict, optional
            Muver analysis data (depth distribution, strand bias, etc.)
            
        Returns
        -------
        dict
            Dictionary mapping sample names to list of plot paths
        """
        if not HAS_PLOTTING:
            self.logger.warning("Matplotlib not available, skipping plots")
            return {}
        
        self.logger.info("Generating individual sample plots...")
        
        # Group mutations by sample
        mutations_by_sample = defaultdict(list)
        for mut in mutations:
            mutations_by_sample[mut.sample_name].append(mut)
        
        plot_paths = {}
        
        for sample in samples:
            sample_name = sample.sample_name
            sample_muts = mutations_by_sample[sample_name]
            
            # Create sample output directory
            sample_dir = self.output_dir / sample_name.replace(' ', '_').replace('/', '_')
            sample_dir.mkdir(parents=True, exist_ok=True)
            
            paths = []
            
            # Generate plots
            if sample_muts:
                paths.append(self.plot_mutation_spectrum(sample_muts, sample_name, sample_dir))
                paths.append(self.plot_chromosome_distribution(sample_muts, sample_name, sample_dir))
                paths.append(self.plot_allele_frequency_distribution(sample_muts, sample_name, sample_dir))
                paths.append(self.plot_variant_quality_scatter(sample_muts, sample_name, sample_dir))
                paths.append(self.plot_quality_dashboard(sample_muts, sample_name, sample_dir))
            
            # Muver-specific plots
            if muver_data:
                if 'depth_distribution' in muver_data:
                    paths.append(self.plot_depth_distribution(
                        muver_data['depth_distribution'], sample_name, sample_dir
                    ))
                if 'strand_bias' in muver_data:
                    paths.append(self.plot_strand_bias_distribution(
                        muver_data['strand_bias'], sample_name, sample_dir
                    ))
            
            plot_paths[sample_name] = [p for p in paths if p]
            self.logger.info(f"  Generated {len(plot_paths[sample_name])} plots for {sample_name}")
        
        return plot_paths
    
    def plot_mutation_spectrum(self, mutations: List, sample_name: str,
                                output_dir: Path) -> Optional[str]:
        """Create 6-category mutation spectrum plot for a sample."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        # Count mutations
        spectrum = Counter()
        for mut in mutations:
            if mut.mutation_type == 'SNP' and len(mut.ref_allele) == 1 and len(mut.alt_allele) == 1:
                ref, alt = mut.ref_allele, mut.alt_allele
                if ref in ['G', 'A']:
                    ref = complement.get(ref, ref)
                    alt = complement.get(alt, alt)
                spectrum[f'{ref}>{alt}'] += 1
        
        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        counts = [spectrum.get(c, 0) for c in categories]
        total = sum(counts)
        fractions = [c/total if total > 0 else 0 for c in counts]
        colors = [self.spectrum_colors[c] for c in categories]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Counts
        bars1 = ax1.bar(categories, counts, color=colors, edgecolor='black', linewidth=1)
        ax1.set_xlabel('Mutation Type', fontsize=12)
        ax1.set_ylabel('Count', fontsize=12)
        ax1.set_title(f'{sample_name}\nMutation Spectrum (Counts)', fontsize=14, fontweight='bold')
        
        for bar, count in zip(bars1, counts):
            if count > 0:
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                        str(count), ha='center', va='bottom', fontsize=10)
        
        # Fractions
        bars2 = ax2.bar(categories, fractions, color=colors, edgecolor='black', linewidth=1)
        ax2.set_xlabel('Mutation Type', fontsize=12)
        ax2.set_ylabel('Fraction', fontsize=12)
        ax2.set_title(f'{sample_name}\nMutation Spectrum (Fractions)', fontsize=14, fontweight='bold')
        ax2.set_ylim(0, max(fractions) * 1.2 if fractions and max(fractions) > 0 else 1)
        
        for bar, frac in zip(bars2, fractions):
            if frac > 0:
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{frac:.1%}', ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        output_file = str(output_dir / 'spectrum.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_chromosome_distribution(self, mutations: List, sample_name: str,
                                      output_dir: Path) -> Optional[str]:
        """Create chromosome distribution plot for a sample."""
        by_chrom = defaultdict(lambda: {'snps': 0, 'indels': 0})
        
        for mut in mutations:
            if mut.mutation_type == 'SNP':
                by_chrom[mut.chromosome]['snps'] += 1
            else:
                by_chrom[mut.chromosome]['indels'] += 1
        
        if not by_chrom:
            return None
        
        chroms = sorted(by_chrom.keys())
        snps = [by_chrom[c]['snps'] for c in chroms]
        indels = [by_chrom[c]['indels'] for c in chroms]
        
        fig, ax = plt.subplots(figsize=(max(10, len(chroms) * 0.8), 6))
        
        x = np.arange(len(chroms))
        width = 0.35
        
        ax.bar(x - width/2, snps, width, label='SNPs', color='steelblue', edgecolor='black')
        ax.bar(x + width/2, indels, width, label='INDELs', color='coral', edgecolor='black')
        
        ax.set_xlabel('Chromosome', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'{sample_name}\nMutations by Chromosome', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(chroms, rotation=45, ha='right')
        ax.legend()
        
        plt.tight_layout()
        output_file = str(output_dir / 'chromosome_distribution.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_allele_frequency_distribution(self, mutations: List, sample_name: str,
                                            output_dir: Path) -> Optional[str]:
        """Create allele frequency distribution plot."""
        # Calculate allele frequencies
        afs = []
        for mut in mutations:
            total = mut.sample_ref_depth + mut.sample_alt_depth
            if total > 0:
                af = mut.sample_alt_depth / total
                afs.append(af)
        
        if not afs:
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.hist(afs, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.mean(afs), color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {np.mean(afs):.2f}')
        ax.axvline(x=np.median(afs), color='green', linestyle='--', linewidth=2,
                   label=f'Median: {np.median(afs):.2f}')
        
        ax.set_xlabel('Allele Frequency', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'{sample_name}\nAlternate Allele Frequency Distribution', 
                     fontsize=14, fontweight='bold')
        ax.legend()
        ax.set_xlim(0, 1)
        
        plt.tight_layout()
        output_file = str(output_dir / 'allele_frequency.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_variant_quality_scatter(self, mutations: List, sample_name: str,
                                      output_dir: Path) -> Optional[str]:
        """Create variant quality vs depth scatter plot."""
        qualities = [mut.quality for mut in mutations]
        depths = [mut.sample_ref_depth + mut.sample_alt_depth for mut in mutations]
        types = [mut.mutation_type for mut in mutations]
        
        if not qualities:
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Color by mutation type
        colors = {'SNP': 'steelblue', 'INS': 'green', 'DEL': 'red'}
        
        for mut_type in ['SNP', 'INS', 'DEL']:
            mask = [t == mut_type for t in types]
            q = [qualities[i] for i in range(len(qualities)) if mask[i]]
            d = [depths[i] for i in range(len(depths)) if mask[i]]
            if q:
                ax.scatter(d, q, c=colors.get(mut_type, 'gray'), 
                          label=f'{mut_type} (n={len(q)})', alpha=0.6, s=30)
        
        ax.set_xlabel('Read Depth', fontsize=12)
        ax.set_ylabel('Variant Quality', fontsize=12)
        ax.set_title(f'{sample_name}\nVariant Quality vs Read Depth', 
                     fontsize=14, fontweight='bold')
        ax.legend()
        
        plt.tight_layout()
        output_file = str(output_dir / 'quality_scatter.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_depth_distribution(self, depth_result, sample_name: str,
                                 output_dir: Path) -> Optional[str]:
        """Plot depth distribution with fitted normal curve."""
        if depth_result is None or not hasattr(depth_result, 'depths'):
            return None
        
        depths = depth_result.depths
        if not depths:
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Histogram
        n, bins, patches = ax.hist(depths, bins=50, density=True, 
                                    color='steelblue', edgecolor='black', alpha=0.7)
        
        # Fitted curve
        if hasattr(depth_result, 'mu') and hasattr(depth_result, 'sigma'):
            from scipy.stats import norm
            x = np.linspace(min(depths), max(depths), 100)
            ax.plot(x, norm.pdf(x, depth_result.mu, depth_result.sigma), 
                   'r-', linewidth=2, label=f'Fitted (μ={depth_result.mu:.1f}, σ={depth_result.sigma:.1f})')
        
        ax.set_xlabel('Depth per Copy', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title(f'{sample_name}\nDepth Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        plt.tight_layout()
        output_file = str(output_dir / 'depth_distribution.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_strand_bias_distribution(self, bias_result, sample_name: str,
                                       output_dir: Path) -> Optional[str]:
        """Plot strand bias distribution with fitted log-normal curve."""
        if bias_result is None or not hasattr(bias_result, 'log_ratios'):
            return None
        
        log_ratios = bias_result.log_ratios
        if not log_ratios:
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Histogram
        ax.hist(log_ratios, bins=50, density=True, 
               color='steelblue', edgecolor='black', alpha=0.7)
        
        # Fitted curve
        if hasattr(bias_result, 'mu') and hasattr(bias_result, 'sigma'):
            from scipy.stats import norm
            x = np.linspace(min(log_ratios), max(log_ratios), 100)
            ax.plot(x, norm.pdf(x, bias_result.mu, bias_result.sigma), 
                   'r-', linewidth=2, label=f'Fitted (μ={bias_result.mu:.3f}, σ={bias_result.sigma:.3f})')
        
        ax.set_xlabel('Log(Forward/Reverse) Ratio', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title(f'{sample_name}\nStrand Bias Distribution', fontsize=14, fontweight='bold')
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax.legend()
        
        plt.tight_layout()
        output_file = str(output_dir / 'strand_bias.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_quality_dashboard(self, mutations: List, sample_name: str,
                                output_dir: Path) -> Optional[str]:
        """Create multi-panel quality metrics dashboard."""
        if not mutations:
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # Panel 1: Mutation type pie chart
        ax1 = axes[0, 0]
        type_counts = Counter(mut.mutation_type for mut in mutations)
        labels = list(type_counts.keys())
        sizes = list(type_counts.values())
        colors_pie = ['steelblue', 'green', 'red'][:len(labels)]
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors_pie, startangle=90)
        ax1.set_title('Mutation Types', fontsize=12, fontweight='bold')
        
        # Panel 2: Quality distribution
        ax2 = axes[0, 1]
        qualities = [mut.quality for mut in mutations]
        ax2.hist(qualities, bins=30, color='steelblue', edgecolor='black', alpha=0.7)
        ax2.set_xlabel('Quality Score')
        ax2.set_ylabel('Count')
        ax2.set_title('Quality Score Distribution', fontsize=12, fontweight='bold')
        
        # Panel 3: Depth distribution
        ax3 = axes[1, 0]
        depths = [mut.sample_ref_depth + mut.sample_alt_depth for mut in mutations]
        ax3.hist(depths, bins=30, color='coral', edgecolor='black', alpha=0.7)
        ax3.set_xlabel('Read Depth')
        ax3.set_ylabel('Count')
        ax3.set_title('Read Depth Distribution', fontsize=12, fontweight='bold')
        
        # Panel 4: Strand bias distribution
        ax4 = axes[1, 1]
        strand_biases = [mut.strand_bias for mut in mutations if hasattr(mut, 'strand_bias')]
        if strand_biases:
            ax4.hist(strand_biases, bins=30, color='green', edgecolor='black', alpha=0.7)
            ax4.set_xlabel('Strand Bias Score')
            ax4.set_ylabel('Count')
        ax4.set_title('Strand Bias Distribution', fontsize=12, fontweight='bold')
        
        plt.suptitle(f'{sample_name} - Quality Metrics Dashboard', 
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(output_dir / 'quality_dashboard.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
