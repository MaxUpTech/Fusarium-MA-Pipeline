#!/usr/bin/env python3
"""
Combined Plots for Fusarium MA Pipeline

Generates visualizations with all samples together:
- Mutation rate heatmap with clustering
- Spectrum clustering (PCA/hierarchical)
- Sample similarity matrix
- Combined distribution plots
- Summary dashboard
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


class CombinedPlotGenerator:
    """
    Generates combined visualizations for all samples.
    
    Creates heatmaps, clustering plots, and summary dashboards
    that show patterns across the entire dataset.
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.CombinedPlots')
        self.output_dir = Path(config['output']['directory']) / 'plots' / 'combined'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        viz_config = config.get('visualization', {})
        self.dpi = viz_config.get('dpi', 300)
        
        if HAS_PLOTTING:
            try:
                plt.style.use('seaborn-v0_8-whitegrid')
            except:
                plt.style.use('seaborn-whitegrid')
    
    def generate_all_plots(self, mutations: List, samples: List,
                           rates: List, stats_results: Dict = None) -> List[str]:
        """Generate all combined plots."""
        if not HAS_PLOTTING:
            self.logger.warning("Matplotlib not available, skipping plots")
            return []
        
        self.logger.info("Generating combined plots...")
        
        plot_paths = []
        
        # Rate heatmap
        path = self.plot_rate_heatmap(rates, samples)
        if path:
            plot_paths.append(path)
        
        # Spectrum clustering
        path = self.plot_spectrum_clustering(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # Sample similarity matrix
        path = self.plot_sample_similarity(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # Combined distributions
        path = self.plot_combined_distributions(mutations, samples, rates)
        if path:
            plot_paths.append(path)
        
        # Summary dashboard
        if stats_results:
            path = self.plot_summary_dashboard(stats_results, samples)
            if path:
                plot_paths.append(path)
        
        self.logger.info(f"Generated {len(plot_paths)} combined plots")
        return plot_paths
    
    def plot_rate_heatmap(self, rates: List, samples: List) -> Optional[str]:
        """Create heatmap of mutation rates with hierarchical clustering."""
        if not rates:
            return None
        
        # Prepare data matrix
        sample_names = [r.sample_name for r in rates]
        genotypes = [r.genotype for r in rates]
        rate_values = [r.rate_per_kb_per_gen for r in rates]
        snp_rates = [r.snp_count / (r.callable_sites * 25) * 1000 if r.callable_sites > 0 else 0 for r in rates]
        indel_rates = [r.indel_count / (r.callable_sites * 25) * 1000 if r.callable_sites > 0 else 0 for r in rates]
        
        # Create data matrix (samples x metrics)
        data = np.array([rate_values, snp_rates, indel_rates]).T
        
        if len(data) < 2:
            return None
        
        fig, ax = plt.subplots(figsize=(10, max(6, len(rates) * 0.3)))
        
        # Create heatmap with seaborn
        try:
            g = sns.clustermap(
                data,
                row_cluster=True if len(data) > 2 else False,
                col_cluster=False,
                yticklabels=sample_names,
                xticklabels=['Total Rate', 'SNP Rate', 'INDEL Rate'],
                cmap='YlOrRd',
                figsize=(10, max(6, len(rates) * 0.3)),
                dendrogram_ratio=(0.1, 0.1)
            )
            
            g.ax_heatmap.set_xlabel('Metric (per kb per generation)', fontsize=12)
            g.ax_heatmap.set_ylabel('Sample', fontsize=12)
            g.fig.suptitle('Mutation Rate Heatmap with Clustering', 
                          fontsize=14, fontweight='bold', y=1.02)
            
            output_file = str(self.output_dir / 'rate_heatmap.png')
            g.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            plt.close()
            
            return output_file
        except Exception as e:
            self.logger.warning(f"Could not create clustered heatmap: {e}")
            
            # Fallback to simple heatmap
            im = ax.imshow(data, aspect='auto', cmap='YlOrRd')
            ax.set_yticks(range(len(sample_names)))
            ax.set_yticklabels(sample_names)
            ax.set_xticks(range(3))
            ax.set_xticklabels(['Total Rate', 'SNP Rate', 'INDEL Rate'])
            plt.colorbar(im, ax=ax, label='Rate (per kb per gen)')
            ax.set_title('Mutation Rate Heatmap', fontsize=14, fontweight='bold')
            
            plt.tight_layout()
            output_file = str(self.output_dir / 'rate_heatmap.png')
            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            plt.close()
            
            return output_file
    
    def plot_spectrum_clustering(self, mutations: List, samples: List) -> Optional[str]:
        """Create PCA or hierarchical clustering of mutation spectra."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        
        # Calculate spectrum per sample
        spectra = {}
        for sample in samples:
            sample_muts = [m for m in mutations if m.sample_name == sample.sample_name 
                          and m.mutation_type == 'SNP']
            
            if not sample_muts:
                continue
            
            spectrum = Counter()
            for mut in sample_muts:
                if len(mut.ref_allele) == 1 and len(mut.alt_allele) == 1:
                    ref, alt = mut.ref_allele, mut.alt_allele
                    if ref in ['G', 'A']:
                        ref = complement.get(ref, ref)
                        alt = complement.get(alt, alt)
                    spectrum[f'{ref}>{alt}'] += 1
            
            total = sum(spectrum.values())
            spectra[sample.sample_name] = [spectrum.get(c, 0) / total if total > 0 else 0 
                                           for c in categories]
        
        if len(spectra) < 2:
            return None
        
        sample_names = list(spectra.keys())
        data = np.array([spectra[s] for s in sample_names])
        
        # Try PCA
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            
            # Standardize
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)
            
            # PCA
            pca = PCA(n_components=min(2, len(sample_names)-1, 6))
            pca_result = pca.fit_transform(data_scaled)
            
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Get genotypes for coloring
            sample_genotype = {s.sample_name: s.genotype for s in samples}
            genotypes = [sample_genotype.get(s, 'Unknown') for s in sample_names]
            unique_genotypes = list(set(genotypes))
            colors = plt.cm.Set2(np.linspace(0, 1, len(unique_genotypes)))
            color_map = {g: colors[i] for i, g in enumerate(unique_genotypes)}
            
            for i, (name, geno) in enumerate(zip(sample_names, genotypes)):
                ax.scatter(pca_result[i, 0], pca_result[i, 1], 
                          c=[color_map[geno]], s=100, label=geno if geno not in [genotypes[j] for j in range(i)] else '')
                ax.annotate(name, (pca_result[i, 0], pca_result[i, 1]), 
                           fontsize=8, alpha=0.7)
            
            ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)', fontsize=12)
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)' if pca.n_components_ > 1 else 'PC2', fontsize=12)
            ax.set_title('Mutation Spectrum PCA', fontsize=14, fontweight='bold')
            
            # Legend without duplicates
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), title='Genotype')
            
            plt.tight_layout()
            output_file = str(self.output_dir / 'spectrum_clustering.png')
            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            plt.close()
            
            return output_file
            
        except ImportError:
            # Fallback to simple bar comparison
            fig, ax = plt.subplots(figsize=(12, 6))
            
            x = np.arange(len(categories))
            width = 0.8 / len(sample_names)
            
            for i, name in enumerate(sample_names[:10]):  # Limit to 10 samples
                offset = (i - len(sample_names[:10])/2 + 0.5) * width
                ax.bar(x + offset, spectra[name], width, label=name)
            
            ax.set_xlabel('Mutation Type')
            ax.set_ylabel('Fraction')
            ax.set_title('Mutation Spectra Comparison', fontsize=14, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(categories)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            plt.tight_layout()
            output_file = str(self.output_dir / 'spectrum_clustering.png')
            plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
            plt.close()
            
            return output_file
    
    def plot_sample_similarity(self, mutations: List, samples: List) -> Optional[str]:
        """Create correlation heatmap of sample similarity."""
        # Count mutations per sample per chromosome
        sample_names = [s.sample_name for s in samples]
        
        # Build feature vectors (mutations per chromosome)
        chromosomes = sorted(set(m.chromosome for m in mutations))
        
        if not chromosomes:
            return None
        
        # Create count matrix
        counts = defaultdict(lambda: defaultdict(int))
        for mut in mutations:
            counts[mut.sample_name][mut.chromosome] += 1
        
        # Convert to matrix
        data = []
        valid_samples = []
        for sample in sample_names:
            if any(counts[sample][c] > 0 for c in chromosomes):
                data.append([counts[sample][c] for c in chromosomes])
                valid_samples.append(sample)
        
        if len(valid_samples) < 2:
            return None
        
        data = np.array(data)
        
        # Calculate correlation matrix
        try:
            corr_matrix = np.corrcoef(data)
        except:
            return None
        
        fig, ax = plt.subplots(figsize=(max(8, len(valid_samples) * 0.5), 
                                        max(6, len(valid_samples) * 0.4)))
        
        im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        
        ax.set_xticks(range(len(valid_samples)))
        ax.set_yticks(range(len(valid_samples)))
        ax.set_xticklabels(valid_samples, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(valid_samples, fontsize=8)
        
        plt.colorbar(im, ax=ax, label='Correlation')
        ax.set_title('Sample Similarity Matrix\n(Mutation Distribution Correlation)', 
                     fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'sample_similarity.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_combined_distributions(self, mutations: List, samples: List,
                                     rates: List) -> Optional[str]:
        """Create overlaid density plots for all samples."""
        if not mutations or not rates:
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Get genotypes for coloring
        sample_genotype = {s.sample_name: s.genotype for s in samples}
        unique_genotypes = list(set(sample_genotype.values()))
        colors = plt.cm.Set2(np.linspace(0, 1, len(unique_genotypes)))
        color_map = {g: colors[i] for i, g in enumerate(unique_genotypes)}
        
        # Panel 1: Mutation rate distribution
        ax1 = axes[0, 0]
        for geno in unique_genotypes:
            geno_rates = [r.rate_per_kb_per_gen for r in rates if r.genotype == geno]
            if geno_rates:
                ax1.hist(geno_rates, bins=10, alpha=0.5, label=geno, color=color_map[geno])
        ax1.set_xlabel('Mutation Rate (per kb per gen)')
        ax1.set_ylabel('Count')
        ax1.set_title('Mutation Rate Distribution', fontweight='bold')
        ax1.legend()
        
        # Panel 2: Quality score distribution
        ax2 = axes[0, 1]
        for geno in unique_genotypes:
            geno_quals = [m.quality for m in mutations if sample_genotype.get(m.sample_name) == geno]
            if geno_quals:
                ax2.hist(geno_quals, bins=30, alpha=0.5, label=geno, color=color_map[geno])
        ax2.set_xlabel('Variant Quality')
        ax2.set_ylabel('Count')
        ax2.set_title('Quality Score Distribution', fontweight='bold')
        ax2.legend()
        
        # Panel 3: Read depth distribution
        ax3 = axes[1, 0]
        for geno in unique_genotypes:
            geno_depths = [m.sample_ref_depth + m.sample_alt_depth for m in mutations 
                          if sample_genotype.get(m.sample_name) == geno]
            if geno_depths:
                ax3.hist(geno_depths, bins=30, alpha=0.5, label=geno, color=color_map[geno])
        ax3.set_xlabel('Read Depth')
        ax3.set_ylabel('Count')
        ax3.set_title('Read Depth Distribution', fontweight='bold')
        ax3.legend()
        
        # Panel 4: Allele frequency distribution
        ax4 = axes[1, 1]
        for geno in unique_genotypes:
            geno_afs = []
            for m in mutations:
                if sample_genotype.get(m.sample_name) == geno:
                    total = m.sample_ref_depth + m.sample_alt_depth
                    if total > 0:
                        geno_afs.append(m.sample_alt_depth / total)
            if geno_afs:
                ax4.hist(geno_afs, bins=20, alpha=0.5, label=geno, color=color_map[geno])
        ax4.set_xlabel('Allele Frequency')
        ax4.set_ylabel('Count')
        ax4.set_title('Allele Frequency Distribution', fontweight='bold')
        ax4.legend()
        
        plt.suptitle('Combined Distribution Plots by Genotype', 
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'combined_distributions.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_summary_dashboard(self, stats_results: Dict, samples: List) -> Optional[str]:
        """Create multi-panel summary dashboard."""
        summary = stats_results.get('summary', {})
        overall_spectrum = stats_results.get('overall_spectrum', {})
        indels = stats_results.get('indels', {})
        
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        
        # Panel 1: Summary stats
        ax1 = axes[0, 0]
        ax1.axis('off')
        stats_text = f"""
Total Samples: {summary.get('total_samples', 0)}
Total Mutations: {summary.get('total_mutations', 0)}
  - SNPs: {summary.get('total_snps', 0)} ({summary.get('snp_fraction', 0):.1%})
  - INDELs: {summary.get('total_indels', 0)} ({summary.get('indel_fraction', 0):.1%})
Mean per Sample: {summary.get('mean_mutations_per_sample', 0):.1f}
Mean Rate: {summary.get('mean_mutation_rate_per_site_per_gen', 0):.2e}
Ts/Tv Ratio: {summary.get('overall_ts_tv_ratio', 0):.2f}
Generations: {summary.get('generations', 0)}
"""
        ax1.text(0.1, 0.9, stats_text, transform=ax1.transAxes, fontsize=11,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
        ax1.set_title('Summary Statistics', fontsize=12, fontweight='bold')
        
        # Panel 2: Mutation spectrum
        ax2 = axes[0, 1]
        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        counts = [overall_spectrum.get(c, 0) for c in categories]
        colors = ['#1E88E5', '#000000', '#D32F2F', '#757575', '#43A047', '#F9A825']
        ax2.bar(categories, counts, color=colors, edgecolor='black')
        ax2.set_ylabel('Count')
        ax2.set_title('Overall Mutation Spectrum', fontsize=12, fontweight='bold')
        
        # Panel 3: Ts/Tv
        ax3 = axes[0, 2]
        ts = overall_spectrum.get('transitions', 0)
        tv = overall_spectrum.get('transversions', 0)
        ax3.pie([ts, tv], labels=['Transitions', 'Transversions'], autopct='%1.1f%%',
               colors=['steelblue', 'coral'], startangle=90)
        ax3.set_title('Transitions vs Transversions', fontsize=12, fontweight='bold')
        
        # Panel 4: INDEL sizes
        ax4 = axes[1, 0]
        size_dist = indels.get('size_distribution', {})
        if size_dist:
            sizes = list(size_dist.keys())
            counts = list(size_dist.values())
            ax4.bar(sizes, counts, color='steelblue', edgecolor='black')
            ax4.set_ylabel('Count')
            ax4.set_title('INDEL Size Distribution', fontsize=12, fontweight='bold')
        else:
            ax4.text(0.5, 0.5, 'No INDEL data', ha='center', va='center', fontsize=12)
            ax4.axis('off')
        
        # Panel 5: Ins/Del ratio
        ax5 = axes[1, 1]
        ins = indels.get('insertions', 0)
        dels = indels.get('deletions', 0)
        if ins + dels > 0:
            ax5.pie([ins, dels], labels=['Insertions', 'Deletions'], autopct='%1.1f%%',
                   colors=['green', 'red'], startangle=90)
            ax5.set_title('Insertion/Deletion Ratio', fontsize=12, fontweight='bold')
        else:
            ax5.text(0.5, 0.5, 'No INDEL data', ha='center', va='center', fontsize=12)
            ax5.axis('off')
        
        # Panel 6: Poisson test result
        ax6 = axes[1, 2]
        ax6.axis('off')
        poisson = stats_results.get('poisson_test', {})
        poisson_text = f"""
Poisson Uniformity Test
-----------------------
Mean mutations: {poisson.get('mean_mutations', 0):.2f}
Variance: {poisson.get('variance', 0):.2f}
Dispersion index: {poisson.get('dispersion_index', 0):.2f}
Chi2 statistic: {poisson.get('chi2_statistic', 0):.2f}
P-value: {poisson.get('p_value', 1.0):.4f}
Conclusion: {poisson.get('conclusion', 'N/A')}
"""
        ax6.text(0.1, 0.9, poisson_text, transform=ax6.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
        ax6.set_title('Statistical Tests', fontsize=12, fontweight='bold')
        
        plt.suptitle('Analysis Summary Dashboard', fontsize=16, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'summary_dashboard.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
