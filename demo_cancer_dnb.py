"""
FRANCIS Demo: Cancer DNB Visualization

Demonstrates the Dynamical Network Biomarker analysis and survival stratification.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json

# Try to import analysis modules
try:
    from models.dnb_cancer import DNBAnalyzer, CANCER_PANELS
    from models.francis_survival import FRANCISSurvival, VALIDATED_TRIPLETS
except ImportError:
    print("Note: Run from repository root directory")
    VALIDATED_TRIPLETS = {
        'GBM': ['CXCL13', 'SOX2', 'LGALS1'],
        'SKCM': ['FABP4', 'CXCL13', 'CD8A'],
        'LUAD': ['CXCL13', 'GAPDH', 'KRAS'],
        'COADREAD': ['CXCL13', 'TOP2A', 'KRAS'],
        'BRCA': ['CXCL13', 'ARG1', 'CD8A'],
        'SARC': ['CXCL13', 'CD8A', 'PDCD1'],
        'OV': ['LDHA', 'MKI67', 'CXCL13'],
    }


def plot_validation_summary():
    """
    Plot summary of TCGA validation results.
    """
    # Validation data
    results = {
        'GBM': {'n': 134, 'deaths': 102, 'p': 1.6e-5},
        'SKCM': {'n': 276, 'deaths': 157, 'p': 0.0007},
        'LUAD': {'n': 392, 'deaths': 135, 'p': 0.0019},
        'COADREAD': {'n': 402, 'deaths': 84, 'p': 0.0046},
        'BRCA': {'n': 790, 'deaths': 114, 'p': 0.0053},
        'SARC': {'n': 169, 'deaths': 71, 'p': 0.0117},
        'OV': {'n': 252, 'deaths': 149, 'p': 0.0152},
    }
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    cancers = list(results.keys())
    
    # Panel 1: Sample sizes
    ax1 = axes[0]
    n_values = [results[c]['n'] for c in cancers]
    bars1 = ax1.barh(cancers, n_values, color='steelblue', edgecolor='black')
    ax1.set_xlabel('Number of Patients')
    ax1.set_title('Sample Size by Cancer Type')
    ax1.set_xlim(0, max(n_values) * 1.1)
    for bar, n in zip(bars1, n_values):
        ax1.text(bar.get_width() + 10, bar.get_y() + bar.get_height()/2, 
                 str(n), va='center', fontsize=9)
    
    # Panel 2: P-values (log scale)
    ax2 = axes[1]
    p_values = [results[c]['p'] for c in cancers]
    colors = ['darkgreen' if p < 0.01 else 'green' for p in p_values]
    bars2 = ax2.barh(cancers, [-np.log10(p) for p in p_values], 
                     color=colors, edgecolor='black')
    ax2.axvline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    ax2.axvline(-np.log10(0.01), color='orange', linestyle='--', label='p=0.01')
    ax2.set_xlabel('-log₁₀(p-value)')
    ax2.set_title('Statistical Significance')
    ax2.legend(loc='lower right')
    
    # Panel 3: CXCL13 universality
    ax3 = axes[2]
    cxcl13_presence = [1 if 'CXCL13' in VALIDATED_TRIPLETS[c] else 0 for c in cancers]
    bars3 = ax3.barh(cancers, cxcl13_presence, color='darkred', edgecolor='black')
    ax3.set_xlabel('CXCL13 in Triplet')
    ax3.set_title('CXCL13 Universality: 7/7 (100%)')
    ax3.set_xlim(0, 1.5)
    ax3.set_xticks([0, 1])
    ax3.set_xticklabels(['No', 'Yes'])
    
    plt.tight_layout()
    return fig


def plot_intermediate_instability_concept():
    """
    Conceptual illustration of the Intermediate Instability Trap.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Generate conceptual curve
    x = np.linspace(0, 35, 100)
    
    # Inverted U shape for intermediate regime
    y = np.zeros_like(x)
    y[x < 5] = 0.05 + 0.01 * x[x < 5]
    y[(x >= 5) & (x <= 25)] = 0.08 + 0.12 * np.sin((x[(x >= 5) & (x <= 25)] - 5) * np.pi / 20)
    y[x > 25] = 0.08 - 0.002 * (x[x > 25] - 25)
    
    # Add noise
    np.random.seed(42)
    y += np.random.normal(0, 0.01, len(y))
    y = np.clip(y, 0, None)
    
    # Plot
    ax.fill_between(x[(x >= 5) & (x <= 25)], 0, y[(x >= 5) & (x <= 25)], 
                    alpha=0.3, color='red', label='Intermediate Regime (5-25)')
    ax.plot(x, y, 'b-', linewidth=2)
    
    # Annotations
    ax.axvline(5, color='red', linestyle='--', alpha=0.5)
    ax.axvline(25, color='red', linestyle='--', alpha=0.5)
    
    ax.annotate('Low Instability\n(Stable)', xy=(2.5, 0.05), ha='center', fontsize=10)
    ax.annotate('INTERMEDIATE\nINSTABILITY TRAP', xy=(15, 0.18), ha='center', 
                fontsize=12, fontweight='bold', color='red')
    ax.annotate('High Instability\n(Catastrophic)', xy=(30, 0.05), ha='center', fontsize=10)
    
    ax.set_xlabel('Aneuploidy Score (Genomic Instability)', fontsize=12)
    ax.set_ylabel('Mortality Rate', fontsize=12)
    ax.set_title('The Intermediate Instability "Trap"\n'
                 'Paradoxically worse survival in middle genomic instability range', 
                 fontsize=14)
    ax.set_xlim(0, 35)
    ax.set_ylim(0, 0.25)
    ax.legend(loc='upper right')
    
    return fig


def plot_cxcl13_mechanism():
    """
    Diagram of CXCL13 mechanism in solid tumors vs B-cell malignancies.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel 1: Solid Tumors (FRANCIS scope)
    ax1 = axes[0]
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title('Solid Tumors (FRANCIS Scope)\nCXCL13 = Protective', fontsize=14, color='green')
    
    # Tumor
    tumor = plt.Circle((5, 5), 2, color='gray', alpha=0.5)
    ax1.add_patch(tumor)
    ax1.text(5, 5, 'TUMOR', ha='center', va='center', fontsize=10, fontweight='bold')
    
    # CXCL13 secretion
    ax1.annotate('CXCL13', xy=(5, 7.2), ha='center', fontsize=10, color='blue')
    ax1.annotate('', xy=(5, 8), xytext=(5, 7.2),
                 arrowprops=dict(arrowstyle='->', color='blue', lw=2))
    
    # B-cells
    for i, (bx, by) in enumerate([(2, 8), (5, 9), (8, 8)]):
        b_cell = plt.Circle((bx, by), 0.5, color='blue', alpha=0.7)
        ax1.add_patch(b_cell)
        ax1.text(bx, by, 'B', ha='center', va='center', color='white', fontweight='bold')
        ax1.annotate('', xy=(5, 7), xytext=(bx, by-0.5),
                     arrowprops=dict(arrowstyle='->', color='blue', lw=1.5))
    
    # TLS
    tls = plt.Circle((5, 2), 1.5, color='lightblue', alpha=0.5)
    ax1.add_patch(tls)
    ax1.text(5, 2, 'TLS\n(Anti-tumor)', ha='center', va='center', fontsize=9)
    ax1.annotate('', xy=(5, 3.5), xytext=(5, 3),
                 arrowprops=dict(arrowstyle='->', color='green', lw=2))
    ax1.text(5, 0.5, '↑ CXCL13 → ↑ B-cells → ↑ TLS → Better Survival', 
             ha='center', fontsize=10, color='green')
    
    # Panel 2: B-cell Malignancies (Excluded)
    ax2 = axes[1]
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.set_aspect('equal')
    ax2.axis('off')
    ax2.set_title('B-cell Malignancies (EXCLUDED)\nCXCL13 = Harmful', fontsize=14, color='red')
    
    # Malignant B-cell
    mal_b = plt.Circle((5, 5), 2, color='darkblue', alpha=0.7)
    ax2.add_patch(mal_b)
    ax2.text(5, 5, 'Malignant\nB-cell', ha='center', va='center', 
             color='white', fontsize=10, fontweight='bold')
    
    # Autocrine loop
    ax2.annotate('CXCL13', xy=(7.5, 6), fontsize=10, color='red')
    ax2.annotate('', xy=(7, 5), xytext=(7.5, 5.8),
                 arrowprops=dict(arrowstyle='->', color='red', lw=2,
                                connectionstyle='arc3,rad=-0.3'))
    ax2.annotate('', xy=(7, 5), xytext=(3, 5),
                 arrowprops=dict(arrowstyle='->', color='red', lw=2,
                                connectionstyle='arc3,rad=0.5'))
    
    ax2.text(5, 8, 'CXCR5', ha='center', fontsize=10)
    ax2.text(5, 1.5, 'Autocrine Loop\n↑ CXCL13 → ↑ Tumor Survival', 
             ha='center', fontsize=10, color='red')
    
    plt.tight_layout()
    return fig


def print_validation_table():
    """Print formatted validation results table."""
    print("\n" + "=" * 80)
    print("FRANCIS VALIDATION RESULTS - TCGA PanCancer Atlas")
    print("=" * 80)
    print(f"{'Cancer':<12} {'Triplet':<32} {'N':>6} {'Deaths':>8} {'P-Value':>12}")
    print("-" * 80)
    
    results = [
        ('GBM', 'CXCL13 + SOX2 + LGALS1', 134, 102, '1.6×10⁻⁵'),
        ('SKCM', 'FABP4 + CXCL13 + CD8A', 276, 157, '0.0007'),
        ('LUAD', 'CXCL13 + GAPDH + KRAS', 392, 135, '0.0019'),
        ('COADREAD', 'CXCL13 + TOP2A + KRAS', 402, 84, '0.0046'),
        ('BRCA', 'CXCL13 + ARG1 + CD8A', 790, 114, '0.0053'),
        ('SARC', 'CXCL13 + CD8A + PDCD1', 169, 71, '0.0117'),
        ('OV', 'LDHA + MKI67 + CXCL13', 252, 149, '0.0152'),
    ]
    
    total_n = 0
    total_deaths = 0
    
    for cancer, triplet, n, deaths, p in results:
        print(f"{cancer:<12} {triplet:<32} {n:>6} {deaths:>8} {p:>12}")
        total_n += n
        total_deaths += deaths
    
    print("-" * 80)
    print(f"{'TOTAL':<12} {'CXCL13 in 7/7 (100%)':<32} {total_n:>6} {total_deaths:>8} {'ALL <0.05':>12}")
    print("=" * 80)


def main():
    """Run demo visualizations."""
    print("\n" + "=" * 60)
    print("FRANCIS - Cancer Dynamics Model Demo")
    print("=" * 60)
    
    print_validation_table()
    
    print("\nGenerating visualizations...")
    
    # Generate plots
    fig1 = plot_validation_summary()
    fig1.savefig('demo_validation_summary.png', dpi=150, bbox_inches='tight')
    print("  Saved: demo_validation_summary.png")
    
    fig2 = plot_intermediate_instability_concept()
    fig2.savefig('demo_instability_trap.png', dpi=150, bbox_inches='tight')
    print("  Saved: demo_instability_trap.png")
    
    fig3 = plot_cxcl13_mechanism()
    fig3.savefig('demo_cxcl13_mechanism.png', dpi=150, bbox_inches='tight')
    print("  Saved: demo_cxcl13_mechanism.png")
    
    print("\nDemo complete!")
    print("\nKey findings:")
    print("  • CXCL13 present in 100% of optimal triplets (7/7 cancers)")
    print("  • All p-values < 0.05 by log-rank test")
    print("  • 2,415 patients, 812 deaths across 3 embryonic lineages")
    print("  • Intermediate instability regime (aneuploidy 5-25) is key")
    
    plt.show()


if __name__ == '__main__':
    main()
