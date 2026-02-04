#!/usr/bin/env python3
"""
Statistical Analysis Plots for Fusarium MA Pipeline

Generates visualizations for each statistical analysis:
- Depth distribution fit
- Strand bias fit
- Poisson Q-Q plot
- Genotype comparison boxplot
- P-value distribution
- Effect size forest plot
- Dispersion analysis
- Bootstrap distribution
"""

import logging
from pathlib import Path
from typing import List, Dict, Optional
from collections import defaultdict

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


class StatisticalPlotGenerator:
    """
    Generates statistical analysis visualizations.
    
    Creates plots for each statistical test showing distributions,
    fits, and diagnostic information.
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.StatisticalPlots')
        self.output_dir = Path(config['output']['directory']) / 'plots' / 'statistics'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        viz_config = config.get('visualization', {})
        self.dpi = viz_config.get('dpi', 300)
        
        if HAS_PLOTTING:
            try:
                plt.style.use('seaborn-v0_8-whitegrid')
            except:
                plt.style.use('seaborn-whitegrid')
    
    def generate_all_plots(self, mutations: List, samples: List,
                           stats_results: Dict, muver_data: Dict = None) -> List[str]:
        """Generate all statistical analysis plots."""
        if not HAS_PLOTTING:
            self.logger.warning("Matplotlib not available, skipping plots")
            return []
        
        self.logger.info("Generating statistical analysis plots...")
        
        plot_paths = []
        
        # Poisson Q-Q plot
        path = self.plot_poisson_qq(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # Genotype comparison boxplot
        path = self.plot_genotype_boxplot(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # P-value distribution
        path = self.plot_pvalue_distribution(mutations)
        if path:
            plot_paths.append(path)
        
        # Dispersion analysis
        path = self.plot_dispersion_analysis(mutations, samples)
        if path:
            plot_paths.append(path)
        
        # Effect size forest plot (from stats results)
        if stats_results and 'genotype_comparison' in stats_results:
            path = self.plot_effect_size_forest(stats_results['genotype_comparison'])
            if path:
                plot_paths.append(path)
        
        # Bootstrap distribution
        if stats_results and 'rates' in stats_results:
            path = self.plot_bootstrap_distribution(stats_results['rates'])
            if path:
                plot_paths.append(path)
        
        # Muver-specific plots
        if muver_data:
            if 'depth_distribution' in muver_data:
                path = self.plot_depth_fit(muver_data['depth_distribution'])
                if path:
                    plot_paths.append(path)
            
            if 'strand_bias' in muver_data:
                path = self.plot_strand_bias_fit(muver_data['strand_bias'])
                if path:
                    plot_paths.append(path)
            
            if 'repeat_indels' in muver_data:
                path = self.plot_repeat_indel_fits(muver_data['repeat_indels'])
                if path:
                    plot_paths.append(path)
        
        self.logger.info(f"Generated {len(plot_paths)} statistical plots")
        return plot_paths
    
    def plot_depth_fit(self, depth_result) -> Optional[str]:
        """Plot depth distribution with fitted normal curve and CI bands."""
        if depth_result is None:
            return None
        
        depths = getattr(depth_result, 'depths', [])
        if not depths:
            return None
        
        from scipy.stats import norm
        
        mu = getattr(depth_result, 'mu', np.mean(depths))
        sigma = getattr(depth_result, 'sigma', np.std(depths))
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Histogram
        n, bins, patches = ax.hist(depths, bins=50, density=True, 
                                    color='steelblue', edgecolor='black', 
                                    alpha=0.7, label='Observed')
        
        # Fitted curve
        x = np.linspace(min(depths), max(depths), 100)
        y = norm.pdf(x, mu, sigma)
        ax.plot(x, y, 'r-', linewidth=2, label=f'Fitted Normal (μ={mu:.1f}, σ={sigma:.1f})')
        
        # CI bands (95%)
        ax.fill_between(x, norm.pdf(x, mu, sigma) * 0.95, 
                        norm.pdf(x, mu, sigma) * 1.05, 
                        color='red', alpha=0.2, label='95% CI')
        
        # Thresholds
        lower_thresh = norm.ppf(0.0001, mu, sigma)
        upper_thresh = norm.ppf(0.9999, mu, sigma)
        ax.axvline(x=lower_thresh, color='orange', linestyle='--', 
                   label=f'Filter thresholds (p=0.0001)')
        ax.axvline(x=upper_thresh, color='orange', linestyle='--')
        
        ax.set_xlabel('Depth per Copy', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title('Depth Distribution with Fitted Normal\nand Filter Thresholds', 
                     fontsize=14, fontweight='bold')
        ax.legend()
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'depth_fit.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_strand_bias_fit(self, bias_result) -> Optional[str]:
        """Plot strand bias distribution with log-normal fit."""
        if bias_result is None:
            return None
        
        log_ratios = getattr(bias_result, 'log_ratios', [])
        if not log_ratios:
            return None
        
        from scipy.stats import norm
        
        mu = getattr(bias_result, 'mu', np.mean(log_ratios))
        sigma = getattr(bias_result, 'sigma', np.std(log_ratios))
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Histogram
        ax.hist(log_ratios, bins=50, density=True, 
               color='steelblue', edgecolor='black', alpha=0.7, label='Observed')
        
        # Fitted curve
        x = np.linspace(min(log_ratios), max(log_ratios), 100)
        y = norm.pdf(x, mu, sigma)
        ax.plot(x, y, 'r-', linewidth=2, label=f'Fitted Gaussian (μ={mu:.3f}, σ={sigma:.3f})')
        
        # Reference line at 0 (balanced)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, label='Balanced (ratio=1)')
        
        # Significance thresholds
        z_thresh = norm.ppf(0.005)  # Two-tailed 1%
        ax.axvline(x=mu + z_thresh * sigma, color='orange', linestyle=':', 
                   label='1% significance threshold')
        ax.axvline(x=mu - z_thresh * sigma, color='orange', linestyle=':')
        
        ax.set_xlabel('Log(Forward/Reverse) Ratio', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title('Strand Bias Distribution with Gaussian Fit', 
                     fontsize=14, fontweight='bold')
        ax.legend()
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'strand_bias_fit.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_poisson_qq(self, mutations: List, samples: List) -> Optional[str]:
        """Create Q-Q plot comparing observed mutation counts to Poisson distribution."""
        from scipy import stats
        
        # Get mutation counts per sample
        counts = []
        for sample in samples:
            n = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            counts.append(n)
        
        if len(counts) < 3:
            return None
        
        counts = np.array(counts)
        mean_count = np.mean(counts)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Q-Q plot
        sorted_counts = np.sort(counts)
        n = len(sorted_counts)
        theoretical = stats.poisson.ppf(np.arange(1, n+1) / (n+1), mean_count)
        
        ax1.scatter(theoretical, sorted_counts, c='steelblue', s=50, alpha=0.7)
        
        # Reference line
        min_val = min(min(theoretical), min(sorted_counts))
        max_val = max(max(theoretical), max(sorted_counts))
        ax1.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label='y=x')
        
        ax1.set_xlabel('Theoretical Quantiles (Poisson)', fontsize=12)
        ax1.set_ylabel('Observed Quantiles', fontsize=12)
        ax1.set_title('Poisson Q-Q Plot', fontsize=12, fontweight='bold')
        ax1.legend()
        
        # Observed vs Expected histogram
        ax2.hist(counts, bins=max(10, len(counts)//3), density=True, 
                color='steelblue', edgecolor='black', alpha=0.7, label='Observed')
        
        x = np.arange(0, max(counts) + 5)
        ax2.plot(x, stats.poisson.pmf(x, mean_count), 'ro-', 
                linewidth=2, markersize=4, label=f'Poisson(λ={mean_count:.1f})')
        
        ax2.set_xlabel('Mutation Count', fontsize=12)
        ax2.set_ylabel('Density', fontsize=12)
        ax2.set_title('Observed vs Expected Distribution', fontsize=12, fontweight='bold')
        ax2.legend()
        
        plt.suptitle('Poisson Uniformity Test', fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'poisson_test.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_genotype_boxplot(self, mutations: List, samples: List) -> Optional[str]:
        """Create boxplot comparing mutation counts by genotype."""
        # Group by genotype
        by_genotype = defaultdict(list)
        for sample in samples:
            n_muts = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            by_genotype[sample.genotype].append(n_muts)
        
        if len(by_genotype) < 1:
            return None
        
        genotypes = list(by_genotype.keys())
        data = [by_genotype[g] for g in genotypes]
        
        fig, ax = plt.subplots(figsize=(max(8, len(genotypes) * 1.5), 6))
        
        # Box plot
        bp = ax.boxplot(data, labels=genotypes, patch_artist=True)
        
        # Color boxes
        colors = plt.cm.Set2(np.linspace(0, 1, len(genotypes)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Add individual points
        for i, (geno, counts) in enumerate(by_genotype.items()):
            x = np.random.normal(i + 1, 0.04, len(counts))
            ax.scatter(x, counts, c='black', s=30, alpha=0.5, zorder=3)
        
        # Add mean markers
        means = [np.mean(by_genotype[g]) for g in genotypes]
        ax.scatter(range(1, len(genotypes) + 1), means, c='red', s=100, 
                  marker='D', zorder=4, label='Mean')
        
        ax.set_xlabel('Genotype', fontsize=12)
        ax.set_ylabel('Mutation Count', fontsize=12)
        ax.set_title('Mutation Count Distribution by Genotype', 
                     fontsize=14, fontweight='bold')
        ax.legend()
        
        # Add Mann-Whitney p-value if 2 genotypes
        if len(genotypes) == 2:
            from scipy.stats import mannwhitneyu
            stat, pval = mannwhitneyu(data[0], data[1], alternative='two-sided')
            ax.text(0.5, 0.95, f'Mann-Whitney U p-value: {pval:.4f}',
                   transform=ax.transAxes, ha='center', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'genotype_boxplot.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_pvalue_distribution(self, mutations: List) -> Optional[str]:
        """Plot distribution of p-values from mutation calling."""
        p_values = []
        
        for mut in mutations:
            if hasattr(mut, 'p_value_fisher') and mut.p_value_fisher < 1.0:
                p_values.append(mut.p_value_fisher)
            elif hasattr(mut, 'p_value_binomial') and mut.p_value_binomial < 1.0:
                p_values.append(mut.p_value_binomial)
        
        if not p_values:
            return None
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Histogram of p-values
        ax1.hist(p_values, bins=50, color='steelblue', edgecolor='black', 
                alpha=0.7, density=True)
        ax1.axhline(y=1, color='red', linestyle='--', linewidth=2, 
                   label='Uniform expectation')
        ax1.set_xlabel('P-value', fontsize=12)
        ax1.set_ylabel('Density', fontsize=12)
        ax1.set_title('P-value Distribution', fontsize=12, fontweight='bold')
        ax1.legend()
        
        # -log10(p) histogram
        log_pvals = [-np.log10(p) for p in p_values if p > 0]
        ax2.hist(log_pvals, bins=50, color='coral', edgecolor='black', alpha=0.7)
        ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', 
                   linewidth=2, label='p=0.05')
        ax2.axvline(x=-np.log10(0.01), color='orange', linestyle='--', 
                   linewidth=2, label='p=0.01')
        ax2.set_xlabel('-log10(P-value)', fontsize=12)
        ax2.set_ylabel('Count', fontsize=12)
        ax2.set_title('-log10(P-value) Distribution', fontsize=12, fontweight='bold')
        ax2.legend()
        
        plt.suptitle('Statistical Significance Distribution', 
                     fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'pvalue_distribution.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_dispersion_analysis(self, mutations: List, samples: List) -> Optional[str]:
        """Plot variance/mean ratio analysis for overdispersion."""
        # Calculate counts per sample
        counts = []
        for sample in samples:
            n = sum(1 for m in mutations if m.sample_name == sample.sample_name)
            counts.append(n)
        
        if len(counts) < 3:
            return None
        
        counts = np.array(counts)
        mean_count = np.mean(counts)
        var_count = np.var(counts, ddof=1)
        dispersion = var_count / mean_count if mean_count > 0 else 0
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Scatter plot of mean vs variance
        ax1.scatter([mean_count], [var_count], s=200, c='steelblue', zorder=3)
        
        # Reference line (Poisson: var = mean)
        max_val = max(mean_count, var_count) * 1.2
        ax1.plot([0, max_val], [0, max_val], 'r--', linewidth=2, 
                label='Poisson (variance = mean)')
        
        ax1.set_xlabel('Mean Mutation Count', fontsize=12)
        ax1.set_ylabel('Variance', fontsize=12)
        ax1.set_title(f'Mean-Variance Relationship\nDispersion Index: {dispersion:.2f}', 
                     fontsize=12, fontweight='bold')
        ax1.legend()
        ax1.set_xlim(0, max_val)
        ax1.set_ylim(0, max_val)
        
        # Dispersion interpretation
        ax2.axis('off')
        interpretation = "Underdispersed" if dispersion < 1 else "Overdispersed" if dispersion > 1 else "Poisson"
        
        text = f"""
Dispersion Analysis Results
---------------------------
Number of samples: {len(counts)}
Mean mutation count: {mean_count:.2f}
Variance: {var_count:.2f}
Dispersion index (var/mean): {dispersion:.2f}

Interpretation: {interpretation}

If dispersion > 1: Overdispersion suggests
  - Mutation hotspots
  - Sample-specific rate variation
  - Possible selection effects

If dispersion < 1: Underdispersion suggests
  - More uniform than expected
  - Possible experimental artifacts
  
If dispersion ≈ 1: Consistent with Poisson
  - Mutations occur randomly
  - No evidence of clustering
"""
        ax2.text(0.1, 0.9, text, transform=ax2.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
        
        plt.suptitle('Dispersion Analysis', fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'dispersion_analysis.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_effect_size_forest(self, genotype_comparison: Dict) -> Optional[str]:
        """Create forest plot of effect sizes with CIs."""
        effect_sizes = genotype_comparison.get('effect_sizes', [])
        
        if not effect_sizes:
            return None
        
        fig, ax = plt.subplots(figsize=(10, max(4, len(effect_sizes) * 0.6)))
        
        y_pos = np.arange(len(effect_sizes))
        
        for i, es in enumerate(effect_sizes):
            d = es['cohens_d']
            ci_lower = es['cohens_d_ci_lower']
            ci_upper = es['cohens_d_ci_upper']
            
            # Point estimate
            color = 'green' if ci_lower > 0 or ci_upper < 0 else 'steelblue'
            ax.scatter(d, i, s=100, color=color, zorder=3)
            
            # Confidence interval
            ax.hlines(i, ci_lower, ci_upper, color=color, linewidth=2, zorder=2)
            
            # Add caps
            ax.scatter([ci_lower, ci_upper], [i, i], marker='|', s=50, color=color, zorder=3)
        
        # Reference lines
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1, zorder=1)
        
        # Effect size guidelines
        for x in [-0.8, -0.5, -0.2, 0.2, 0.5, 0.8]:
            ax.axvline(x=x, color='lightgray', linestyle=':', linewidth=0.5, zorder=0)
        
        # Labels
        labels = [f"{es['genotype1']} vs {es['genotype2']}\n({es['interpretation']})" 
                 for es in effect_sizes]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        
        ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12)
        ax.set_title("Genotype Effect Sizes with 95% CI\nForest Plot", 
                     fontsize=14, fontweight='bold')
        
        # Add effect size scale
        ax.text(0.02, 0.02, 'Small: 0.2 | Medium: 0.5 | Large: 0.8',
               transform=ax.transAxes, fontsize=9, color='gray')
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'effect_size_forest.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_bootstrap_distribution(self, rates: List) -> Optional[str]:
        """Plot bootstrap distribution of mutation rates."""
        if not rates:
            return None
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Get rate values and CIs
        sample_names = [r.sample_name for r in rates]
        rate_values = [r.rate_per_kb_per_gen for r in rates]
        ci_lower = [r.ci_lower for r in rates]
        ci_upper = [r.ci_upper for r in rates]
        boot_lower = [r.bootstrap_ci_lower for r in rates]
        boot_upper = [r.bootstrap_ci_upper for r in rates]
        
        x = np.arange(len(rates))
        
        # Plot bars with error bars
        ax.bar(x, rate_values, color='steelblue', edgecolor='black', alpha=0.7)
        
        # Poisson CI
        ax.errorbar(x - 0.1, rate_values, 
                   yerr=[np.array(rate_values) - np.array(ci_lower),
                         np.array(ci_upper) - np.array(rate_values)],
                   fmt='none', color='red', capsize=3, label='Poisson CI')
        
        # Bootstrap CI
        ax.errorbar(x + 0.1, rate_values,
                   yerr=[np.array(rate_values) - np.array(boot_lower),
                         np.array(boot_upper) - np.array(rate_values)],
                   fmt='none', color='green', capsize=3, label='Bootstrap CI')
        
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Mutation Rate (per kb per generation)', fontsize=12)
        ax.set_title('Mutation Rates with Confidence Intervals', 
                     fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([s[:15] + '...' if len(s) > 15 else s for s in sample_names], 
                          rotation=45, ha='right')
        ax.legend()
        
        plt.tight_layout()
        output_file = str(self.output_dir / 'bootstrap_distribution.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
    
    def plot_repeat_indel_fits(self, repeat_result) -> Optional[str]:
        """Plot repeat indel rate fits."""
        if repeat_result is None:
            return None
        
        rates = getattr(repeat_result, 'rates', {})
        fits = getattr(repeat_result, 'fits', {})
        
        if not rates or not fits:
            return None
        
        from .muver_models.fitting import logistic
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
        
        plot_idx = 0
        for event in ['insertion', 'deletion']:
            for repeat_len in [1, 2]:
                if plot_idx >= 4:
                    break
                
                ax = axes[plot_idx]
                event_rates = rates.get(event, {}).get(repeat_len, {})
                event_fits = fits.get(event, {}).get(repeat_len)
                
                if event_rates and event_fits:
                    tract_lengths = sorted(event_rates.keys())
                    observed = [np.log10(event_rates[t]) for t in tract_lengths]
                    
                    fit = event_fits
                    fitted = [logistic(t, fit.x0, fit.L, fit.M, fit.k) for t in tract_lengths]
                    
                    ax.scatter(tract_lengths, observed, c='black', s=30, label='Observed')
                    ax.plot(tract_lengths, fitted, 'b-', linewidth=2, label='Fitted')
                    
                    ax.set_xlabel('Repeat Tract Length (bp)')
                    ax.set_ylabel('log10(Indel Rate)')
                    ax.set_title(f'{event.capitalize()} - {repeat_len}bp repeat unit')
                    ax.legend()
                else:
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center')
                    ax.set_title(f'{event.capitalize()} - {repeat_len}bp repeat unit')
                
                plot_idx += 1
        
        plt.suptitle('Repeat INDEL Rate Fits', fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_file = str(self.output_dir / 'repeat_indel_fits.png')
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        return output_file
