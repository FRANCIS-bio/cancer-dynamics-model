"""
FRANCIS Demo - Cancer Dynamics Visualization

Demonstrates critical transition detection with plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from dnb_cancer import (
    detect_critical_transition, 
    compute_dii, 
    compute_dnb_scores,
    analyze_immunometabolic_panel
)


def simulate_cancer_progression(n_genes=500, n_timepoints=100, transition_point=70):
    """
    Simulate gene expression data with critical transition.
    
    Before transition_point: stable dynamics
    After transition_point: approaching bifurcation
    """
    np.random.seed(42)
    
    # Base expression with mild fluctuations
    data = np.random.randn(n_genes, n_timepoints) * 0.5
    
    # Select DNB candidate genes (will show critical slowing down)
    n_dnb = 30
    dnb_genes = np.random.choice(n_genes, n_dnb, replace=False)
    
    for t in range(n_timepoints):
        if t > transition_point:
            # After transition: increased variance in DNB genes
            progress = (t - transition_point) / (n_timepoints - transition_point)
            for g in dnb_genes:
                data[g, t] += np.random.randn() * (1 + 3 * progress)
                data[g, t] += progress * 2  # Trend
    
    # Cumulative to simulate expression levels
    data = np.cumsum(data * 0.1, axis=1) + 5
    
    return data, dnb_genes, transition_point


def plot_dnb_analysis(data, true_dnb_genes, transition_point):
    """
    Create visualization of DNB analysis results.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel A: Expression trajectories
    ax = axes[0, 0]
    # Plot non-DNB genes (gray)
    non_dnb = [i for i in range(min(100, data.shape[0])) if i not in true_dnb_genes[:10]]
    for g in non_dnb[:20]:
        ax.plot(data[g, :], color='gray', alpha=0.3, linewidth=0.5)
    # Plot DNB genes (red)
    for g in true_dnb_genes[:10]:
        ax.plot(data[g, :], color='red', alpha=0.7, linewidth=1.5)
    ax.axvline(x=transition_point, color='black', linestyle='--', label='Transition point')
    ax.set_xlabel('Time')
    ax.set_ylabel('Expression')
    ax.set_title('A. Gene Expression Trajectories\n(Red = DNB genes, Gray = Others)')
    ax.legend()
    
    # Panel B: Rolling DII
    ax = axes[0, 1]
    window_size = 15
    dii_values = []
    time_points = []
    
    for t in range(window_size, data.shape[1]):
        window_data = data[:100, t-window_size:t]  # Use subset for speed
        dii, _ = compute_dii(window_data)
        dii_values.append(dii)
        time_points.append(t)
    
    ax.plot(time_points, dii_values, 'b-', linewidth=2)
    ax.axvline(x=transition_point, color='black', linestyle='--', label='True transition')
    ax.axhline(y=0.5, color='red', linestyle=':', label='Critical threshold')
    ax.fill_between(time_points, dii_values, alpha=0.3)
    ax.set_xlabel('Time')
    ax.set_ylabel('DII')
    ax.set_title('B. Dynamical Instability Index Over Time')
    ax.legend()
    
    # Panel C: DNB scores distribution
    ax = axes[1, 0]
    scores = compute_dnb_scores(data)
    
    # Separate DNB and non-DNB
    dnb_scores = scores[true_dnb_genes]
    non_dnb_indices = [i for i in range(len(scores)) if i not in true_dnb_genes]
    non_dnb_scores = scores[non_dnb_indices]
    
    ax.hist(non_dnb_scores, bins=30, alpha=0.5, label='Non-DNB genes', color='gray')
    ax.hist(dnb_scores, bins=15, alpha=0.7, label='True DNB genes', color='red')
    ax.set_xlabel('DNB Score')
    ax.set_ylabel('Count')
    ax.set_title('C. DNB Score Distribution')
    ax.legend()
    
    # Panel D: Variance over time for DNB vs non-DNB
    ax = axes[1, 1]
    window = 10
    
    dnb_var = []
    non_dnb_var = []
    times = []
    
    for t in range(window, data.shape[1]):
        dnb_var.append(np.mean(np.var(data[true_dnb_genes, t-window:t], axis=1)))
        non_dnb_var.append(np.mean(np.var(data[non_dnb_indices[:100], t-window:t], axis=1)))
        times.append(t)
    
    ax.plot(times, dnb_var, 'r-', linewidth=2, label='DNB genes')
    ax.plot(times, non_dnb_var, 'gray', linewidth=2, label='Non-DNB genes')
    ax.axvline(x=transition_point, color='black', linestyle='--', label='Transition point')
    ax.set_xlabel('Time')
    ax.set_ylabel('Mean Variance')
    ax.set_title('D. Variance Amplification in DNB Genes')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('dnb_analysis_results.png', dpi=150, bbox_inches='tight')
    print("Saved: dnb_analysis_results.png")
    plt.show()


def main():
    print("=" * 60)
    print("FRANCIS - Cancer Dynamics Demo")
    print("=" * 60)
    
    # Simulate data
    print("\n[1] Simulating cancer progression data...")
    data, true_dnb_genes, transition_point = simulate_cancer_progression()
    print(f"    Data shape: {data.shape}")
    print(f"    True DNB genes: {len(true_dnb_genes)}")
    print(f"    Transition point: t={transition_point}")
    
    # Run detection
    print("\n[2] Running critical transition detection...")
    detected_genes, score, is_critical = detect_critical_transition(data)
    print(f"    Tipping score: {score:.3f}")
    print(f"    Critical: {is_critical}")
    
    # Check detection accuracy
    overlap = len(set(detected_genes) & set(true_dnb_genes))
    print(f"    Detection overlap: {overlap}/{len(true_dnb_genes)} true DNB genes found")
    
    # Compute final DII
    print("\n[3] Computing DII components...")
    dii, components = compute_dii(data[:100, :])
    print(f"    DII: {dii:.4f}")
    print(f"    - Entropy: {components['entropy']:.4f}")
    print(f"    - Divergence: {components['divergence']:.4f}")
    print(f"    - Trajectory: {components['trajectory_instability']:.4f}")
    
    # Generate plots
    print("\n[4] Generating visualization...")
    plot_dnb_analysis(data, true_dnb_genes, transition_point)
    
    print("\n" + "=" * 60)
    print("Demo complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
