"""
FRANCIS - Framework for Real-time Analysis of Nonlinear Critical Instability Signals

Dynamical Network Biomarkers (DNB) for Cancer Immunometabolic Collapse
Detects critical transition pre-cancer tipping points using nonlinear dynamics

Patent: U.S. Provisional Application No. 63/968,064
"""

import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import jensenshannon
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')


def compute_sample_entropy(signal, m=2, r=0.2):
    """
    Compute sample entropy of a time series.
    
    Parameters:
    -----------
    signal : array-like
        Input time series
    m : int
        Embedding dimension
    r : float
        Tolerance (fraction of std)
    
    Returns:
    --------
    float : Sample entropy value
    """
    N = len(signal)
    r_threshold = r * np.std(signal)
    
    def count_matches(template_length):
        count = 0
        for i in range(N - template_length):
            for j in range(i + 1, N - template_length):
                if np.max(np.abs(signal[i:i+template_length] - signal[j:j+template_length])) < r_threshold:
                    count += 1
        return count
    
    A = count_matches(m + 1)
    B = count_matches(m)
    
    if B == 0 or A == 0:
        return 0.0
    
    return -np.log(A / B)


def compute_multiscale_entropy(signal, scales=range(1, 6), m=2, r=0.2):
    """
    Compute multiscale entropy across multiple time scales.
    
    Parameters:
    -----------
    signal : array-like
        Input time series
    scales : range
        Time scales to analyze
    m : int
        Embedding dimension
    r : float
        Tolerance
    
    Returns:
    --------
    array : MSE values at each scale
    """
    mse = []
    for scale in scales:
        # Coarse-grain the signal
        n_segments = len(signal) // scale
        coarse = np.array([np.mean(signal[i*scale:(i+1)*scale]) for i in range(n_segments)])
        if len(coarse) > 10:
            mse.append(compute_sample_entropy(coarse, m, r))
        else:
            mse.append(0.0)
    return np.array(mse)


def compute_jensen_shannon_divergence(window1, window2, bins=20):
    """
    Compute Jensen-Shannon divergence between two time windows.
    
    Parameters:
    -----------
    window1 : array-like
        First time window
    window2 : array-like
        Second time window
    bins : int
        Number of histogram bins
    
    Returns:
    --------
    float : JSD value (0 = identical, 1 = maximally different)
    """
    # Create histograms
    all_data = np.concatenate([window1, window2])
    bin_edges = np.linspace(np.min(all_data), np.max(all_data), bins + 1)
    
    hist1, _ = np.histogram(window1, bins=bin_edges, density=True)
    hist2, _ = np.histogram(window2, bins=bin_edges, density=True)
    
    # Add small epsilon to avoid division by zero
    hist1 = hist1 + 1e-10
    hist2 = hist2 + 1e-10
    
    # Normalize
    hist1 = hist1 / np.sum(hist1)
    hist2 = hist2 / np.sum(hist2)
    
    return jensenshannon(hist1, hist2)


def estimate_lyapunov_exponent(signal, dt=1.0, min_sep=5):
    """
    Estimate largest Lyapunov exponent from time series.
    
    Parameters:
    -----------
    signal : array-like
        Input time series
    dt : float
        Time step
    min_sep : int
        Minimum temporal separation for neighbors
    
    Returns:
    --------
    float : Estimated Lyapunov exponent
    """
    N = len(signal)
    if N < 20:
        return 0.0
    
    divergences = []
    
    for i in range(N - min_sep - 1):
        # Find nearest neighbor (excluding temporal neighbors)
        distances = np.abs(signal - signal[i])
        distances[max(0, i-min_sep):min(N, i+min_sep+1)] = np.inf
        
        j = np.argmin(distances)
        
        if j < N - 1 and i < N - 1:
            initial_dist = np.abs(signal[i] - signal[j])
            final_dist = np.abs(signal[i+1] - signal[j+1])
            
            if initial_dist > 1e-10 and final_dist > 1e-10:
                divergences.append(np.log(final_dist / initial_dist))
    
    if len(divergences) > 0:
        return np.mean(divergences) / dt
    return 0.0


def compute_dnb_scores(timeseries, window=10):
    """
    Compute Dynamical Network Biomarker (DNB) scores.
    
    The DNB theory identifies a dominant group of molecules that:
    1. Show dramatically increased variance (SD)
    2. Show increased internal correlation (PCC_in)
    3. Show decreased external correlation (PCC_out)
    
    DNB score = SD * PCC_in / PCC_out
    
    Parameters:
    -----------
    timeseries : ndarray
        Gene expression matrix (n_genes x n_timepoints)
    window : int
        Analysis window size
    
    Returns:
    --------
    ndarray : DNB scores for each gene
    """
    n_genes, n_time = timeseries.shape
    
    if n_time < window:
        window = n_time
    
    # Standard deviation in recent window
    sd = np.std(timeseries[:, -window:], axis=1)
    
    # Pearson correlations
    pcc_matrix = np.zeros((n_genes, n_genes))
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            corr, _ = pearsonr(timeseries[i, -window:], timeseries[j, -window:])
            pcc_matrix[i, j] = corr
            pcc_matrix[j, i] = corr
    
    # Mean absolute correlation for each gene
    mean_corr = np.mean(np.abs(pcc_matrix), axis=1)
    
    # DNB score: high variance * high correlation
    # (in full implementation, would separate internal vs external correlation)
    dnb_score = sd * mean_corr
    
    return dnb_score


def compute_dii(timeseries, reference_window=None, current_window=10):
    """
    Compute Dynamical Instability Index (DII).
    
    DII(t) = w1 * Entropy(t) + w2 * Divergence(t) + w3 * TrajectoryInstability(t)
    
    Parameters:
    -----------
    timeseries : ndarray
        Multivariate time series (n_variables x n_timepoints)
    reference_window : tuple
        (start, end) indices for healthy reference
    current_window : int
        Size of current analysis window
    
    Returns:
    --------
    float : DII value
    dict : Component scores
    """
    n_vars, n_time = timeseries.shape
    
    if reference_window is None:
        reference_window = (0, min(20, n_time // 3))
    
    # Weights (can be tuned)
    w1, w2, w3 = 0.33, 0.33, 0.34
    
    # Component 1: Entropy change
    entropy_scores = []
    for i in range(n_vars):
        ref_entropy = compute_sample_entropy(timeseries[i, reference_window[0]:reference_window[1]])
        curr_entropy = compute_sample_entropy(timeseries[i, -current_window:])
        entropy_scores.append(np.abs(curr_entropy - ref_entropy))
    entropy_component = np.mean(entropy_scores)
    
    # Component 2: Jensen-Shannon divergence
    divergence_scores = []
    for i in range(n_vars):
        jsd = compute_jensen_shannon_divergence(
            timeseries[i, reference_window[0]:reference_window[1]],
            timeseries[i, -current_window:]
        )
        divergence_scores.append(jsd)
    divergence_component = np.mean(divergence_scores)
    
    # Component 3: Trajectory instability (Lyapunov-like)
    lyapunov_scores = []
    for i in range(n_vars):
        le = estimate_lyapunov_exponent(timeseries[i, -current_window:])
        lyapunov_scores.append(max(0, le))  # Positive = unstable
    trajectory_component = np.mean(lyapunov_scores)
    
    # Compute DII
    dii = w1 * entropy_component + w2 * divergence_component + w3 * trajectory_component
    
    components = {
        'entropy': entropy_component,
        'divergence': divergence_component,
        'trajectory_instability': trajectory_component,
        'weights': (w1, w2, w3)
    }
    
    return dii, components


def detect_critical_transition(gene_expr_data, threshold=3.0):
    """
    Main function: Detect critical transition and return DNB genes.
    
    Parameters:
    -----------
    gene_expr_data : ndarray
        Gene expression matrix (n_genes x n_timepoints)
    threshold : float
        DNB score threshold for tipping point detection
    
    Returns:
    --------
    dnb_genes : ndarray
        Indices of top DNB genes
    tipping_score : float
        Overall tipping point score
    is_critical : bool
        Whether system is approaching critical transition
    """
    scores = compute_dnb_scores(gene_expr_data)
    
    # Top 50 DNB genes (or 5% of genes, whichever is smaller)
    n_dnb = min(50, len(scores) // 20)
    dnb_genes = np.argsort(scores)[-n_dnb:]
    
    # Tipping score is mean of top DNB scores
    tipping_score = np.mean(scores[dnb_genes])
    
    # Critical if exceeds threshold
    is_critical = tipping_score > threshold
    
    return dnb_genes, tipping_score, is_critical


def analyze_immunometabolic_panel(expression_data, gene_names, marker_genes):
    """
    Analyze immunometabolic markers for critical transition.
    
    Parameters:
    -----------
    expression_data : ndarray
        Expression matrix (n_genes x n_timepoints)
    gene_names : list
        List of gene names corresponding to rows
    marker_genes : list
        List of immunometabolic marker genes to analyze
    
    Returns:
    --------
    dict : Analysis results including DII and DNB scores
    """
    # Find indices of marker genes
    gene_name_to_idx = {name: i for i, name in enumerate(gene_names)}
    marker_indices = [gene_name_to_idx[g] for g in marker_genes if g in gene_name_to_idx]
    
    if len(marker_indices) < 3:
        return {'error': 'Insufficient marker genes found'}
    
    # Extract marker expression
    marker_expr = expression_data[marker_indices, :]
    
    # Compute DII
    dii, dii_components = compute_dii(marker_expr)
    
    # Compute DNB scores
    dnb_scores = compute_dnb_scores(marker_expr)
    
    # Map back to gene names
    found_markers = [marker_genes[i] for i, g in enumerate(marker_genes) if g in gene_name_to_idx]
    
    results = {
        'dii': dii,
        'dii_components': dii_components,
        'dnb_scores': dict(zip(found_markers, dnb_scores)),
        'top_dnb_markers': sorted(zip(found_markers, dnb_scores), key=lambda x: x[1], reverse=True)[:10],
        'is_critical': dii > 0.5  # Threshold for marker panel
    }
    
    return results


# Example usage and demonstration
if __name__ == "__main__":
    print("=" * 60)
    print("FRANCIS - Cancer Dynamics Model")
    print("Dynamical Network Biomarker Analysis")
    print("=" * 60)
    
    # Simulated cancer progression data
    # (In real use, this would be TCGA or clinical data)
    np.random.seed(42)
    n_genes = 1000
    n_timepoints = 50
    
    # Simulate progression toward critical transition
    # Early timepoints: stable; Late timepoints: approaching bifurcation
    base_expr = np.random.randn(n_genes, n_timepoints)
    
    # Add trend toward instability in subset of genes (DNB candidates)
    dnb_candidate_genes = np.random.choice(n_genes, 50, replace=False)
    for g in dnb_candidate_genes:
        # Increasing variance over time
        base_expr[g, :] += np.linspace(0, 3, n_timepoints) * np.random.randn(n_timepoints)
        # Trend
        base_expr[g, :] += np.linspace(0, 2, n_timepoints)
    
    data = np.cumsum(base_expr * 0.1, axis=1)
    
    # Run DNB analysis
    print("\n[1] Running DNB Analysis...")
    genes, score, is_critical = detect_critical_transition(data)
    print(f"    Tipping score: {score:.3f}")
    print(f"    Critical transition detected: {is_critical}")
    print(f"    Number of DNB genes identified: {len(genes)}")
    
    # Compute DII
    print("\n[2] Computing Dynamical Instability Index (DII)...")
    dii, components = compute_dii(data[:100, :])  # Use subset for speed
    print(f"    DII: {dii:.4f}")
    print(f"    - Entropy component: {components['entropy']:.4f}")
    print(f"    - Divergence component: {components['divergence']:.4f}")
    print(f"    - Trajectory instability: {components['trajectory_instability']:.4f}")
    
    # Simulate immunometabolic markers
    print("\n[3] Immunometabolic Panel Analysis...")
    marker_genes = ['LDHA', 'HK2', 'GLUT1', 'PDK1', 'IL6', 'STAT3', 'CD274', 'IDO1', 'HIF1A', 'ARG1']
    fake_gene_names = marker_genes + [f'GENE_{i}' for i in range(n_genes - len(marker_genes))]
    
    results = analyze_immunometabolic_panel(data, fake_gene_names, marker_genes)
    print(f"    Immunometabolic DII: {results['dii']:.4f}")
    print(f"    Critical state: {results['is_critical']}")
    print(f"    Top DNB markers:")
    for gene, score in results['top_dnb_markers'][:5]:
        print(f"      - {gene}: {score:.3f}")
    
    print("\n" + "=" * 60)
    print("Analysis complete.")
    print("=" * 60)
