#!/usr/bin/env python3
"""
Enhanced Reporting Module for Fusarium MA Pipeline

Generates comprehensive reports integrating all visualization modules:
- Individual sample plots
- Group comparison plots
- Combined analysis plots
- Statistical analysis plots
- Comprehensive HTML report with embedded visualizations
"""

import logging
import json
import base64
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False

from visualization import (
    SamplePlotGenerator,
    GroupPlotGenerator,
    CombinedPlotGenerator,
    StatisticalPlotGenerator
)


class EnhancedReporting:
    """
    Generate comprehensive reports for MA analysis.
    
    Integrates all visualization modules and generates:
    1. Individual sample plots
    2. Group comparison plots  
    3. Combined analysis plots
    4. Statistical analysis plots
    5. Comprehensive HTML report
    """
    
    def __init__(self, config: dict):
        self.config = config
        self.logger = logging.getLogger('FusariumMA.Reporting')
        self.output_dir = Path(config['output']['directory'])
        self.plots_dir = self.output_dir / 'plots'
        self.reports_dir = self.output_dir / 'reports'
        
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize plot generators
        self.sample_plotter = SamplePlotGenerator(config)
        self.group_plotter = GroupPlotGenerator(config)
        self.combined_plotter = CombinedPlotGenerator(config)
        self.stats_plotter = StatisticalPlotGenerator(config)
        
        # Store generated plot paths
        self.plot_paths = {
            'individual': {},
            'groups': [],
            'combined': [],
            'statistics': []
        }
    
    def generate_all_plots(self, mutations: List, samples: List, ancestor,
                           rates: List, stats_results: Dict,
                           muver_data: Dict = None) -> Dict[str, List[str]]:
        """
        Generate all plots from all visualization modules.
        
        Parameters
        ----------
        mutations : list
            List of Mutation objects
        samples : list
            List of Sample objects
        ancestor : Sample
            Ancestor sample
        rates : list
            List of MutationRateResult objects
        stats_results : dict
            Results from statistical analysis
        muver_data : dict, optional
            Muver analysis data (depth, strand bias, repeats)
            
        Returns
        -------
        dict
            Dictionary of all generated plot paths organized by category
        """
        self.logger.info("Generating all visualization plots...")
        
        # Individual sample plots
        self.logger.info("Generating individual sample plots...")
        self.plot_paths['individual'] = self.sample_plotter.generate_all_plots(
            mutations, samples, muver_data
        )
        
        # Group plots
        self.logger.info("Generating group comparison plots...")
        self.plot_paths['groups'] = self.group_plotter.generate_all_plots(
            mutations, samples, rates, stats_results
        )
        
        # Combined plots
        self.logger.info("Generating combined analysis plots...")
        self.plot_paths['combined'] = self.combined_plotter.generate_all_plots(
            mutations, samples, rates, stats_results
        )
        
        # Statistical plots
        self.logger.info("Generating statistical analysis plots...")
        self.plot_paths['statistics'] = self.stats_plotter.generate_all_plots(
            mutations, samples, stats_results, muver_data
        )
        
        total_plots = (
            sum(len(v) for v in self.plot_paths['individual'].values()) +
            len(self.plot_paths['groups']) +
            len(self.plot_paths['combined']) +
            len(self.plot_paths['statistics'])
        )
        
        self.logger.info(f"Generated {total_plots} total plots")
        
        return self.plot_paths
    
    def generate_html_report(self, mutations: List, stats: Dict,
                             samples: List, ancestor,
                             muver_data: Dict = None) -> str:
        """
        Generate comprehensive HTML report with embedded visualizations.
        
        Parameters
        ----------
        mutations : list
            List of Mutation objects
        stats : dict
            Statistical analysis results
        samples : list
            List of Sample objects
        ancestor : Sample
            Ancestor sample
        muver_data : dict, optional
            Muver analysis data
            
        Returns
        -------
        str
            Path to generated HTML report
        """
        self.logger.info("Generating comprehensive HTML report...")
        
        summary = stats.get('summary', {})
        rates = stats.get('rates', [])
        spectrum = stats.get('overall_spectrum', {})
        indels = stats.get('indels', {})
        poisson = stats.get('poisson_test', {})
        genotype_comp = stats.get('genotype_comparison', {})
        additional_tests = stats.get('additional_tests', {})
        
        # Generate HTML
        html = self._generate_html_header()
        html += self._generate_summary_section(summary, samples, ancestor)
        html += self._generate_mutation_spectrum_section(spectrum)
        html += self._generate_mutation_rates_section(rates)
        html += self._generate_indel_section(indels)
        html += self._generate_statistical_tests_section(poisson, genotype_comp, additional_tests)
        html += self._generate_plots_section()
        html += self._generate_sample_details_section(mutations, samples, ancestor)
        html += self._generate_html_footer()
        
        output_file = str(self.reports_dir / 'analysis_report.html')
        with open(output_file, 'w') as f:
            f.write(html)
        
        self.logger.info(f"HTML report saved to {output_file}")
        return output_file
    
    def _generate_html_header(self) -> str:
        """Generate HTML header with styles."""
        return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Fusarium MA Analysis Report</title>
    <style>
        * {{
            box-sizing: border-box;
        }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
            color: #333;
        }}
        .header {{
            background: linear-gradient(135deg, #1a5276, #2980b9);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
        }}
        .header p {{
            margin: 10px 0 0 0;
            opacity: 0.9;
        }}
        .section {{
            background: white;
            padding: 25px;
            margin-bottom: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #1a5276;
            border-bottom: 3px solid #2980b9;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        .section h3 {{
            color: #2980b9;
            margin-top: 20px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #2980b9;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .stat-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #f8f9fa, #e9ecef);
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            border-left: 4px solid #2980b9;
        }}
        .stat-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #1a5276;
        }}
        .stat-card .label {{
            color: #666;
            margin-top: 5px;
        }}
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .highlight {{
            background-color: #fff3cd;
            padding: 15px;
            border-radius: 5px;
            border-left: 4px solid #ffc107;
            margin: 15px 0;
        }}
        .success {{
            background-color: #d4edda;
            border-left-color: #28a745;
        }}
        .info {{
            background-color: #d1ecf1;
            border-left-color: #17a2b8;
        }}
        .warning {{
            background-color: #fff3cd;
            border-left-color: #ffc107;
        }}
        .collapsible {{
            background-color: #2980b9;
            color: white;
            cursor: pointer;
            padding: 15px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 16px;
            border-radius: 5px;
            margin-top: 10px;
        }}
        .collapsible:after {{
            content: '\\002B';
            color: white;
            font-weight: bold;
            float: right;
        }}
        .active:after {{
            content: "\\2212";
        }}
        .content {{
            padding: 0 18px;
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.2s ease-out;
            background-color: #f1f1f1;
            border-radius: 0 0 5px 5px;
        }}
        .footer {{
            text-align: center;
            color: #666;
            padding: 20px;
            font-size: 0.9em;
        }}
        .tab {{
            overflow: hidden;
            border: 1px solid #ccc;
            background-color: #f1f1f1;
            border-radius: 5px 5px 0 0;
        }}
        .tab button {{
            background-color: inherit;
            float: left;
            border: none;
            outline: none;
            cursor: pointer;
            padding: 14px 16px;
            transition: 0.3s;
        }}
        .tab button:hover {{
            background-color: #ddd;
        }}
        .tab button.active {{
            background-color: #2980b9;
            color: white;
        }}
        .tabcontent {{
            display: none;
            padding: 20px;
            border: 1px solid #ccc;
            border-top: none;
            border-radius: 0 0 5px 5px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Fusarium graminearum MA Analysis</h1>
        <p>Mutation Accumulation Analysis with Enhanced Statistical Models</p>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
"""
    
    def _generate_summary_section(self, summary: Dict, samples: List, ancestor) -> str:
        """Generate summary statistics section."""
        ma_lines = [s for s in samples if s.sample_name != ancestor.sample_name]
        
        return f"""
    <div class="section">
        <h2>Summary Statistics</h2>
        <div class="stat-grid">
            <div class="stat-card">
                <div class="value">{summary.get('total_samples', len(ma_lines))}</div>
                <div class="label">MA Lines Analyzed</div>
            </div>
            <div class="stat-card">
                <div class="value">{summary.get('total_mutations', 0)}</div>
                <div class="label">Total Mutations</div>
            </div>
            <div class="stat-card">
                <div class="value">{summary.get('total_snps', 0)}</div>
                <div class="label">SNPs</div>
            </div>
            <div class="stat-card">
                <div class="value">{summary.get('total_indels', 0)}</div>
                <div class="label">INDELs</div>
            </div>
            <div class="stat-card">
                <div class="value">{summary.get('mean_mutations_per_sample', 0):.1f}</div>
                <div class="label">Mean per Sample</div>
            </div>
            <div class="stat-card">
                <div class="value">{summary.get('mean_mutation_rate_per_site_per_gen', 0):.2e}</div>
                <div class="label">Mean Rate (per site/gen)</div>
            </div>
            <div class="stat-card">
                <div class="value">{summary.get('overall_ts_tv_ratio', 0):.2f}</div>
                <div class="label">Ts/Tv Ratio</div>
            </div>
            <div class="stat-card">
                <div class="value">{summary.get('generations', 25)}</div>
                <div class="label">Generations</div>
            </div>
        </div>
    </div>
"""
    
    def _generate_mutation_spectrum_section(self, spectrum: Dict) -> str:
        """Generate mutation spectrum section."""
        categories = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        
        rows = ""
        total = sum(spectrum.get(c, 0) for c in categories)
        for cat in categories:
            count = spectrum.get(cat, 0)
            frac = count / total * 100 if total > 0 else 0
            rows += f"<tr><td>{cat}</td><td>{count}</td><td>{frac:.1f}%</td></tr>"
        
        return f"""
    <div class="section">
        <h2>Mutation Spectrum</h2>
        <table>
            <tr><th>Mutation Type</th><th>Count</th><th>Fraction</th></tr>
            {rows}
            <tr style="font-weight: bold; background-color: #e9ecef;">
                <td>Total SNPs</td><td>{total}</td><td>100%</td>
            </tr>
        </table>
        <div class="highlight info">
            <strong>Transitions:</strong> {spectrum.get('transitions', 0)} | 
            <strong>Transversions:</strong> {spectrum.get('transversions', 0)} | 
            <strong>Ts/Tv Ratio:</strong> {spectrum.get('ts_tv_ratio', 0):.2f}
        </div>
    </div>
"""
    
    def _generate_mutation_rates_section(self, rates: List) -> str:
        """Generate mutation rates section."""
        if not rates:
            return ""
        
        rows = ""
        for r in rates:
            rows += f"""
            <tr>
                <td>{r.sample_name}</td>
                <td>{r.genotype}</td>
                <td>{r.total_mutations}</td>
                <td>{r.snp_count}</td>
                <td>{r.indel_count}</td>
                <td>{r.rate_per_kb_per_gen:.2e}</td>
                <td>{r.ci_lower:.2e} - {r.ci_upper:.2e}</td>
            </tr>
"""
        
        return f"""
    <div class="section">
        <h2>Mutation Rates by Sample</h2>
        <table>
            <tr>
                <th>Sample</th>
                <th>Genotype</th>
                <th>Total</th>
                <th>SNPs</th>
                <th>INDELs</th>
                <th>Rate (per kb/gen)</th>
                <th>95% CI</th>
            </tr>
            {rows}
        </table>
    </div>
"""
    
    def _generate_indel_section(self, indels: Dict) -> str:
        """Generate INDEL analysis section."""
        size_dist = indels.get('size_distribution', {})
        
        return f"""
    <div class="section">
        <h2>INDEL Analysis</h2>
        <div class="stat-grid">
            <div class="stat-card">
                <div class="value">{indels.get('insertions', 0)}</div>
                <div class="label">Insertions</div>
            </div>
            <div class="stat-card">
                <div class="value">{indels.get('deletions', 0)}</div>
                <div class="label">Deletions</div>
            </div>
            <div class="stat-card">
                <div class="value">{indels.get('ins_del_ratio', 0):.2f}</div>
                <div class="label">Ins/Del Ratio</div>
            </div>
            <div class="stat-card">
                <div class="value">{indels.get('mean_insertion_size', 0):.1f} bp</div>
                <div class="label">Mean Insertion Size</div>
            </div>
            <div class="stat-card">
                <div class="value">{indels.get('mean_deletion_size', 0):.1f} bp</div>
                <div class="label">Mean Deletion Size</div>
            </div>
        </div>
        <h3>Size Distribution</h3>
        <table>
            <tr><th>Size Category</th><th>Count</th></tr>
            <tr><td>1 bp</td><td>{size_dist.get('1bp', 0)}</td></tr>
            <tr><td>2-5 bp</td><td>{size_dist.get('2-5bp', 0)}</td></tr>
            <tr><td>6-10 bp</td><td>{size_dist.get('6-10bp', 0)}</td></tr>
            <tr><td>>10 bp</td><td>{size_dist.get('>10bp', 0)}</td></tr>
        </table>
    </div>
"""
    
    def _generate_statistical_tests_section(self, poisson: Dict, genotype_comp: Dict,
                                             additional_tests: Dict) -> str:
        """Generate statistical tests section."""
        # Poisson test
        poisson_conclusion = poisson.get('conclusion', 'N/A')
        poisson_class = 'success' if 'uniform' in poisson_conclusion.lower() else 'warning'
        
        # Additional tests
        additional_html = ""
        if additional_tests:
            ks_test = additional_tests.get('ks_test')
            shapiro_test = additional_tests.get('shapiro_test')
            perm_test = additional_tests.get('permutation_test')
            
            if ks_test:
                additional_html += f"""
                <tr>
                    <td>Kolmogorov-Smirnov</td>
                    <td>{getattr(ks_test, 'statistic', 'N/A'):.4f}</td>
                    <td>{getattr(ks_test, 'p_value', 'N/A'):.4f}</td>
                    <td>{getattr(ks_test, 'conclusion', 'N/A')}</td>
                </tr>
"""
            if shapiro_test:
                additional_html += f"""
                <tr>
                    <td>Shapiro-Wilk</td>
                    <td>{getattr(shapiro_test, 'statistic', 'N/A'):.4f}</td>
                    <td>{getattr(shapiro_test, 'p_value', 'N/A'):.4f}</td>
                    <td>{'Normal' if getattr(shapiro_test, 'is_normal', False) else 'Non-normal'}</td>
                </tr>
"""
            if perm_test:
                additional_html += f"""
                <tr>
                    <td>Permutation Test</td>
                    <td>{getattr(perm_test, 'observed_statistic', 'N/A'):.4f}</td>
                    <td>{getattr(perm_test, 'p_value', 'N/A'):.4f}</td>
                    <td>n={getattr(perm_test, 'n_permutations', 'N/A')}</td>
                </tr>
"""
        
        # Effect sizes
        effect_html = ""
        effect_sizes = genotype_comp.get('effect_sizes', [])
        for es in effect_sizes:
            effect_html += f"""
            <tr>
                <td>{es['genotype1']} vs {es['genotype2']}</td>
                <td>{es['cohens_d']:.3f}</td>
                <td>[{es['cohens_d_ci_lower']:.3f}, {es['cohens_d_ci_upper']:.3f}]</td>
                <td>{es['interpretation']}</td>
            </tr>
"""
        
        return f"""
    <div class="section">
        <h2>Statistical Tests</h2>
        
        <h3>Poisson Uniformity Test</h3>
        <div class="highlight {poisson_class}">
            <strong>Conclusion:</strong> {poisson_conclusion}<br>
            <strong>Dispersion Index:</strong> {poisson.get('dispersion_index', 0):.2f} | 
            <strong>Chi2 Statistic:</strong> {poisson.get('chi2_statistic', 0):.2f} | 
            <strong>P-value:</strong> {poisson.get('p_value', 1.0):.4f}
        </div>
        
        <h3>Additional Statistical Tests</h3>
        <table>
            <tr><th>Test</th><th>Statistic</th><th>P-value</th><th>Result</th></tr>
            <tr>
                <td>Poisson Dispersion</td>
                <td>{poisson.get('chi2_statistic', 0):.4f}</td>
                <td>{poisson.get('p_value', 1.0):.4f}</td>
                <td>{poisson_conclusion}</td>
            </tr>
            {additional_html}
        </table>
        
        <h3>Effect Sizes (Genotype Comparisons)</h3>
        <table>
            <tr><th>Comparison</th><th>Cohen's d</th><th>95% CI</th><th>Interpretation</th></tr>
            {effect_html if effect_html else '<tr><td colspan="4">No pairwise comparisons available</td></tr>'}
        </table>
    </div>
"""
    
    def _generate_plots_section(self) -> str:
        """Generate plots gallery section."""
        html = """
    <div class="section">
        <h2>Visualizations</h2>
        <div class="tab">
            <button class="tablinks active" onclick="openTab(event, 'GroupPlots')">Group Plots</button>
            <button class="tablinks" onclick="openTab(event, 'CombinedPlots')">Combined Plots</button>
            <button class="tablinks" onclick="openTab(event, 'StatisticsPlots')">Statistical Plots</button>
            <button class="tablinks" onclick="openTab(event, 'SamplePlots')">Individual Samples</button>
        </div>
"""
        
        # Group plots tab
        html += '<div id="GroupPlots" class="tabcontent" style="display: block;">'
        html += '<div class="plot-grid">'
        for path in self.plot_paths.get('groups', []):
            if path and Path(path).exists():
                rel_path = Path(path).relative_to(self.output_dir.parent) if self.output_dir.parent in Path(path).parents else path
                html += f'<div class="plot-container"><img src="../{rel_path}" alt="Group Plot"></div>'
        html += '</div></div>'
        
        # Combined plots tab
        html += '<div id="CombinedPlots" class="tabcontent">'
        html += '<div class="plot-grid">'
        for path in self.plot_paths.get('combined', []):
            if path and Path(path).exists():
                rel_path = Path(path).relative_to(self.output_dir.parent) if self.output_dir.parent in Path(path).parents else path
                html += f'<div class="plot-container"><img src="../{rel_path}" alt="Combined Plot"></div>'
        html += '</div></div>'
        
        # Statistics plots tab
        html += '<div id="StatisticsPlots" class="tabcontent">'
        html += '<div class="plot-grid">'
        for path in self.plot_paths.get('statistics', []):
            if path and Path(path).exists():
                rel_path = Path(path).relative_to(self.output_dir.parent) if self.output_dir.parent in Path(path).parents else path
                html += f'<div class="plot-container"><img src="../{rel_path}" alt="Statistics Plot"></div>'
        html += '</div></div>'
        
        # Individual sample plots tab
        html += '<div id="SamplePlots" class="tabcontent">'
        for sample_name, paths in self.plot_paths.get('individual', {}).items():
            html += f'<button class="collapsible">{sample_name}</button>'
            html += '<div class="content"><div class="plot-grid">'
            for path in paths:
                if path and Path(path).exists():
                    rel_path = Path(path).relative_to(self.output_dir.parent) if self.output_dir.parent in Path(path).parents else path
                    html += f'<div class="plot-container"><img src="../{rel_path}" alt="{sample_name}"></div>'
            html += '</div></div>'
        html += '</div>'
        
        html += '</div>'
        return html
    
    def _generate_sample_details_section(self, mutations: List, samples: List, ancestor) -> str:
        """Generate sample details section."""
        from collections import defaultdict
        
        mutations_by_sample = defaultdict(list)
        for mut in mutations:
            mutations_by_sample[mut.sample_name].append(mut)
        
        sample_rows = ""
        for sample in samples:
            if sample.sample_name == ancestor.sample_name:
                continue
            
            muts = mutations_by_sample[sample.sample_name]
            snps = sum(1 for m in muts if m.mutation_type == 'SNP')
            indels = sum(1 for m in muts if m.mutation_type != 'SNP')
            
            sample_rows += f"""
            <tr>
                <td>{sample.sample_name}</td>
                <td>{sample.genotype}</td>
                <td>{sample.generation}</td>
                <td>{getattr(sample, 'replicate', 'N/A')}</td>
                <td>{len(muts)}</td>
                <td>{snps}</td>
                <td>{indels}</td>
            </tr>
"""
        
        return f"""
    <div class="section">
        <h2>Sample Details</h2>
        <h3>Ancestor</h3>
        <table>
            <tr><th>Sample Name</th><th>Genotype</th><th>Generation</th></tr>
            <tr>
                <td>{ancestor.sample_name}</td>
                <td>{ancestor.genotype}</td>
                <td>{ancestor.generation}</td>
            </tr>
        </table>
        
        <h3>MA Lines</h3>
        <table>
            <tr>
                <th>Sample</th>
                <th>Genotype</th>
                <th>Generation</th>
                <th>Replicate</th>
                <th>Total Mutations</th>
                <th>SNPs</th>
                <th>INDELs</th>
            </tr>
            {sample_rows}
        </table>
    </div>
"""
    
    def _generate_html_footer(self) -> str:
        """Generate HTML footer with scripts."""
        return """
    <div class="footer">
        <p>Fusarium MA Pipeline - Enhanced Analysis Report</p>
        <p>Integrating muver statistical models for improved accuracy</p>
    </div>
    
    <script>
    function openTab(evt, tabName) {
        var i, tabcontent, tablinks;
        tabcontent = document.getElementsByClassName("tabcontent");
        for (i = 0; i < tabcontent.length; i++) {
            tabcontent[i].style.display = "none";
        }
        tablinks = document.getElementsByClassName("tablinks");
        for (i = 0; i < tablinks.length; i++) {
            tablinks[i].className = tablinks[i].className.replace(" active", "");
        }
        document.getElementById(tabName).style.display = "block";
        evt.currentTarget.className += " active";
    }
    
    var coll = document.getElementsByClassName("collapsible");
    for (var i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.maxHeight) {
                content.style.maxHeight = null;
            } else {
                content.style.maxHeight = content.scrollHeight + "px";
            }
        });
    }
    </script>
</body>
</html>
"""
    
    def run(self, mutations: List, stats: Dict, samples: List, ancestor,
            muver_data: Dict = None):
        """
        Generate all reports and plots.
        
        Parameters
        ----------
        mutations : list
            List of Mutation objects
        stats : dict
            Statistical analysis results
        samples : list
            List of Sample objects
        ancestor : Sample
            Ancestor sample
        muver_data : dict, optional
            Muver analysis data
        """
        self.logger.info("Generating comprehensive reports...")
        
        rates = stats.get('rates', [])
        
        # Generate all plots
        self.generate_all_plots(mutations, samples, ancestor, rates, stats, muver_data)
        
        # Generate HTML report
        self.generate_html_report(mutations, stats, samples, ancestor, muver_data)
        
        self.logger.info("Report generation complete")


# Backward compatibility
ComprehensiveReporting = EnhancedReporting
