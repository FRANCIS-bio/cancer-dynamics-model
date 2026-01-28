"""
Dynamical Network Biomarker (DNB) Module

Implements DNB theory for detecting critical transitions in cancer.

Based on: Chen L, et al. (2012). Detecting early-warning signals for sudden 
deterioration of complex diseases by dynamical network biomarkers. Scientific Reports.

Patent: U.S. Provisional Application No. 63/968,064
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from typing import Dict, List, Tuple, Optional
import warnings


class DNBAnalyzer:
    """
    Dynamical Network Biomarker analyzer for detecting critical transitions.
    
    DNB theory predicts that before a system tips, a group of variables will show:
    1. Increased variance (critical slowing down)
    2. Increased correlation with each other (dominant group forms)
    3. Decreased correlation with the rest of the system
    
    The DNB score combines these three signals.
    """
    
    def __init__(self, marker_names: List[str] = None):
        """
        Initialize DNB analyzer.
        
        Parameters
        ----------
        marker_names : list of str, optional
            Names of markers to analyze. If None, uses all columns.
        """
        self.marker_names = marker_names
        self.results = None
        
    def compute_variance_score(self, data: pd.DataFrame) -> pd.Series:
        """Compute coefficient of variation for each marker."""
        cv = data.std() / data.mean().abs()
        cv = cv.replace([np.inf, -np.inf], np.nan).fillna(0)
        return cv
    
    def compute_intra_correlation(self, data: pd.DataFrame, 
                                   group: List[str]) -> float:
        """Compute mean absolute correlation within a group of markers."""
        if len(group) < 2:
            return 0.0
        
        group_data = data[group]
        corr_matrix = group_data.corr().abs()
        
        # Get upper triangle (excluding diagonal)
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
        upper_vals = corr_matrix.values[mask]
        
        return np.nanmean(upper_vals) if len(upper_vals) > 0 else 0.0
    
    def compute_inter_correlation(self, data: pd.DataFrame,
                                   group: List[str]) -> float:
        """Compute mean absolute correlation between group and rest."""
        rest = [c for c in data.columns if c not in group]
        
        if len(rest) == 0 or len(group) == 0:
            return 0.0
        
        correlations = []
        for g in group:
            for r in rest:
                corr = data[g].corr(data[r])
                if not np.isnan(corr):
                    correlations.append(abs(corr))
        
        return np.mean(correlations) if correlations else 0.0
    
    def compute_dnb_score(self, data: pd.DataFrame, 
                          candidate_group: List[str]) -> float:
        """
        Compute DNB score for a candidate marker group.
        
        DNB = (variance * intra_correlation) / inter_correlation
        
        High DNB indicates the group is behaving like a dynamical network
        biomarker, signaling approach to a critical transition.
        """
        # Get markers present in data
        group = [m for m in candidate_group if m in data.columns]
        
        if len(group) < 2:
            return 0.0
        
        # Component scores
        variance = self.compute_variance_score(data[group]).mean()
        intra_corr = self.compute_intra_correlation(data, group)
        inter_corr = self.compute_inter_correlation(data, group)
        
        # Avoid division by zero
        if inter_corr < 0.01:
            inter_corr = 0.01
        
        dnb = (variance * intra_corr) / inter_corr
        
        return dnb
    
    def compute_dii(self, data: pd.DataFrame, 
                    markers: List[str] = None) -> Dict:
        """
        Compute Dynamical Instability Index (DII).
        
        DII combines multiple indicators of critical transition:
        - Entropy (disorder in the system)
        - Divergence (spread of trajectories)
        - Instability (variance accumulation)
        
        Parameters
        ----------
        data : pd.DataFrame
            Time series or sample data (samples x markers)
        markers : list, optional
            Specific markers to use
            
        Returns
        -------
        dict
            DII components and total score
        """
        if markers is None:
            markers = list(data.columns)
        
        marker_data = data[[m for m in markers if m in data.columns]]
        
        if marker_data.empty:
            return {'dii': 0, 'entropy': 0, 'divergence': 0, 'instability': 0}
        
        # Entropy component (normalized)
        variance = marker_data.var()
        entropy = np.log(variance + 1).mean()
        entropy_norm = entropy / (np.log(marker_data.max().max() + 1) + 1e-10)
        
        # Divergence component (spread)
        if len(marker_data) > 1:
            distances = pdist(marker_data.values, metric='euclidean')
            divergence = np.mean(distances) / (np.std(distances) + 1e-10)
            divergence_norm = min(divergence / 10, 1.0)
        else:
            divergence_norm = 0
        
        # Instability component (coefficient of variation)
        cv = (marker_data.std() / marker_data.mean().abs()).replace([np.inf, -np.inf], 0)
        instability = cv.mean()
        instability_norm = min(instability, 1.0)
        
        # Combined DII
        dii = (entropy_norm + divergence_norm + instability_norm) / 3
        
        return {
            'dii': dii,
            'entropy': entropy_norm,
            'divergence': divergence_norm,
            'instability': instability_norm
        }
    
    def find_dnb_group(self, data: pd.DataFrame, 
                       min_size: int = 3,
                       max_size: int = 10) -> Tuple[List[str], float]:
        """
        Find the optimal DNB group from data.
        
        Uses a greedy approach to identify markers that collectively
        show DNB behavior (high intra-correlation, low inter-correlation,
        high variance).
        
        Parameters
        ----------
        data : pd.DataFrame
            Expression or marker data
        min_size : int
            Minimum group size
        max_size : int
            Maximum group size
            
        Returns
        -------
        best_group : list
            Optimal marker group
        best_score : float
            DNB score for the group
        """
        markers = list(data.columns)
        
        # Start with highest variance markers
        variances = data.var().sort_values(ascending=False)
        seed_markers = variances.head(max_size).index.tolist()
        
        best_group = seed_markers[:min_size]
        best_score = self.compute_dnb_score(data, best_group)
        
        # Greedy expansion
        for m in seed_markers[min_size:]:
            candidate = best_group + [m]
            score = self.compute_dnb_score(data, candidate)
            
            if score > best_score:
                best_group = candidate
                best_score = score
        
        return best_group, best_score
    
    def analyze(self, data: pd.DataFrame, 
                candidate_groups: List[List[str]] = None) -> Dict:
        """
        Run full DNB analysis.
        
        Parameters
        ----------
        data : pd.DataFrame
            Expression or marker data
        candidate_groups : list of lists, optional
            Pre-defined marker groups to test
            
        Returns
        -------
        dict
            Analysis results
        """
        if candidate_groups is None:
            # Find optimal group automatically
            best_group, dnb_score = self.find_dnb_group(data)
            candidate_groups = [best_group]
        
        results = []
        for group in candidate_groups:
            group_in_data = [m for m in group if m in data.columns]
            
            if len(group_in_data) < 2:
                continue
            
            dnb = self.compute_dnb_score(data, group_in_data)
            dii = self.compute_dii(data, group_in_data)
            
            results.append({
                'group': group_in_data,
                'dnb_score': dnb,
                'dii': dii['dii'],
                'entropy': dii['entropy'],
                'divergence': dii['divergence'],
                'instability': dii['instability']
            })
        
        # Sort by DNB score
        results = sorted(results, key=lambda x: x['dnb_score'], reverse=True)
        
        self.results = {
            'groups': results,
            'best_group': results[0] if results else None,
            'n_samples': len(data),
            'n_markers': len(data.columns)
        }
        
        return self.results
    
    def summary(self) -> str:
        """Return formatted summary of DNB analysis."""
        if self.results is None:
            return "No analysis run. Call analyze() first."
        
        r = self.results
        best = r['best_group']
        
        if best is None:
            return "No valid marker groups found."
        
        return f"""
DNB Analysis Results
====================
Samples: {r['n_samples']}
Markers: {r['n_markers']}

Best DNB Group: {', '.join(best['group'])}
DNB Score: {best['dnb_score']:.4f}

DII Components:
  Entropy:     {best['entropy']:.4f}
  Divergence:  {best['divergence']:.4f}
  Instability: {best['instability']:.4f}
  
Total DII: {best['dii']:.4f}

Interpretation:
  DNB > 1.0: Strong signal of approaching transition
  DNB 0.5-1.0: Moderate instability
  DNB < 0.5: System relatively stable
"""


# Pre-defined marker panels for cancer types
CANCER_PANELS = {
    'immunotherapy': ['CD274', 'PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT'],
    'metabolic': ['LDHA', 'HIF1A', 'GAPDH', 'PKM', 'SLC2A1', 'IDO1'],
    'immune': ['CD8A', 'CD4', 'FOXP3', 'CXCL13', 'CD19', 'MS4A1'],
    'proliferation': ['MKI67', 'TOP2A', 'PCNA', 'MCM2', 'CDK1', 'CCNB1'],
    'stemness': ['SOX2', 'NANOG', 'POU5F1', 'PROM1', 'ALDH1A1', 'CD44'],
}


if __name__ == '__main__':
    print("DNB Cancer Module")
    print("=" * 40)
    print("Implements Dynamical Network Biomarker theory")
    print("for detecting critical transitions in cancer.")
    print()
    print("Available marker panels:")
    for name, markers in CANCER_PANELS.items():
        print(f"  {name}: {', '.join(markers)}")
