"""
FRANCIS Survival Analysis Module

Implements Dynamical Immunometabolic Index (DII) computation and survival
stratification for cancer patients within intermediate genomic instability regimes.

Patent: U.S. Provisional Application No. 63/970,059
"""

import numpy as np
import pandas as pd
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from typing import Dict, List, Tuple, Optional
import json
import os

# Validated tissue-specific biomarker triplets
VALIDATED_TRIPLETS = {
    'GBM': ['CXCL13', 'SOX2', 'LGALS1'],
    'SKCM': ['FABP4', 'CXCL13', 'CD8A'],
    'LUAD': ['CXCL13', 'GAPDH', 'KRAS'],
    'COADREAD': ['CXCL13', 'TOP2A', 'KRAS'],
    'BRCA': ['CXCL13', 'ARG1', 'CD8A'],
    'SARC': ['CXCL13', 'CD8A', 'PDCD1'],
    'OV': ['LDHA', 'MKI67', 'CXCL13'],
}

# Intermediate instability regime bounds
ANEUPLOIDY_LOW = 5
ANEUPLOIDY_HIGH = 25


class FRANCISSurvival:
    """
    FRANCIS Survival Stratification System
    
    Identifies high-risk cancer patients within the intermediate genomic
    instability regime using CXCL13-anchored biomarker triplets.
    """
    
    def __init__(self, cancer_type: str):
        """
        Initialize FRANCIS model for a specific cancer type.
        
        Parameters
        ----------
        cancer_type : str
            TCGA cancer code (GBM, SKCM, LUAD, COADREAD, BRCA, SARC, OV)
        """
        if cancer_type not in VALIDATED_TRIPLETS:
            raise ValueError(
                f"Cancer type '{cancer_type}' not validated. "
                f"Supported types: {list(VALIDATED_TRIPLETS.keys())}"
            )
        
        self.cancer_type = cancer_type
        self.triplet = VALIDATED_TRIPLETS[cancer_type]
        self.results = None
        
    def compute_dii(self, expression_data: pd.DataFrame) -> pd.Series:
        """
        Compute Dynamical Immunometabolic Index (DII) from expression data.
        
        DII = mean(z1, z2, z3)
        
        Where z1, z2, z3 are z-score normalized expression values of the
        tissue-specific biomarker triplet.
        
        Parameters
        ----------
        expression_data : pd.DataFrame
            Gene expression matrix (samples x genes)
            
        Returns
        -------
        pd.Series
            DII values for each sample
        """
        # Check that all triplet genes are present
        missing = [g for g in self.triplet if g not in expression_data.columns]
        if missing:
            raise ValueError(f"Missing genes in expression data: {missing}")
        
        # Extract triplet expression
        triplet_expr = expression_data[self.triplet]
        
        # Z-score normalize each gene
        z_scores = (triplet_expr - triplet_expr.mean()) / triplet_expr.std()
        
        # DII = mean of z-scores
        dii = z_scores.mean(axis=1)
        dii.name = 'DII'
        
        return dii
    
    def filter_intermediate_instability(
        self,
        aneuploidy_scores: pd.Series,
        low: float = ANEUPLOIDY_LOW,
        high: float = ANEUPLOIDY_HIGH
    ) -> pd.Index:
        """
        Filter samples to intermediate genomic instability regime.
        
        Parameters
        ----------
        aneuploidy_scores : pd.Series
            Aneuploidy scores indexed by sample ID
        low : float
            Lower bound of intermediate regime (default: 5)
        high : float
            Upper bound of intermediate regime (default: 25)
            
        Returns
        -------
        pd.Index
            Sample IDs within intermediate regime
        """
        mask = (aneuploidy_scores >= low) & (aneuploidy_scores <= high)
        return aneuploidy_scores[mask].index
    
    def stratify(
        self,
        dii: pd.Series,
        method: str = 'median'
    ) -> pd.Series:
        """
        Stratify samples into risk groups based on DII.
        
        Parameters
        ----------
        dii : pd.Series
            DII values indexed by sample ID
        method : str
            Stratification method ('median', 'tertile', 'quartile')
            
        Returns
        -------
        pd.Series
            Risk group labels ('High Risk' or 'Low Risk')
        """
        if method == 'median':
            threshold = dii.median()
            groups = pd.Series(
                np.where(dii < threshold, 'High Risk', 'Low Risk'),
                index=dii.index
            )
        elif method == 'tertile':
            groups = pd.qcut(dii, 3, labels=['High Risk', 'Medium Risk', 'Low Risk'])
        elif method == 'quartile':
            groups = pd.qcut(dii, 4, labels=['Q1 (Highest Risk)', 'Q2', 'Q3', 'Q4 (Lowest Risk)'])
        else:
            raise ValueError(f"Unknown method: {method}")
        
        groups.name = 'Risk_Group'
        return groups
    
    def analyze(
        self,
        expression_data: pd.DataFrame,
        survival_time: pd.Series,
        survival_event: pd.Series,
        aneuploidy_scores: pd.Series,
        filter_instability: bool = True
    ) -> Dict:
        """
        Run full FRANCIS survival analysis.
        
        Parameters
        ----------
        expression_data : pd.DataFrame
            Gene expression matrix (samples x genes)
        survival_time : pd.Series
            Overall survival time in days
        survival_event : pd.Series
            Event indicator (1 = death, 0 = censored)
        aneuploidy_scores : pd.Series
            Aneuploidy scores for each sample
        filter_instability : bool
            Whether to filter to intermediate instability regime
            
        Returns
        -------
        dict
            Analysis results including p-value, hazard ratio, KM data
        """
        # Align all data
        common_samples = (
            expression_data.index
            .intersection(survival_time.index)
            .intersection(survival_event.index)
            .intersection(aneuploidy_scores.index)
        )
        
        if len(common_samples) == 0:
            raise ValueError("No common samples across input data")
        
        # Filter to intermediate instability regime
        if filter_instability:
            intermediate_samples = self.filter_intermediate_instability(
                aneuploidy_scores.loc[common_samples]
            )
            samples = common_samples.intersection(intermediate_samples)
        else:
            samples = common_samples
        
        if len(samples) < 20:
            raise ValueError(f"Insufficient samples after filtering: {len(samples)}")
        
        # Compute DII
        dii = self.compute_dii(expression_data.loc[samples])
        
        # Stratify
        risk_groups = self.stratify(dii)
        
        # Prepare survival data
        T = survival_time.loc[samples]
        E = survival_event.loc[samples]
        
        # Log-rank test
        high_risk_mask = risk_groups == 'High Risk'
        lr_result = logrank_test(
            T[high_risk_mask], T[~high_risk_mask],
            E[high_risk_mask], E[~high_risk_mask]
        )
        
        # Kaplan-Meier fits
        kmf_high = KaplanMeierFitter()
        kmf_low = KaplanMeierFitter()
        
        kmf_high.fit(T[high_risk_mask], E[high_risk_mask], label='High Risk (Low DII)')
        kmf_low.fit(T[~high_risk_mask], E[~high_risk_mask], label='Low Risk (High DII)')
        
        # Store results
        self.results = {
            'cancer_type': self.cancer_type,
            'triplet': self.triplet,
            'n_samples': len(samples),
            'n_events': int(E.sum()),
            'n_high_risk': int(high_risk_mask.sum()),
            'n_low_risk': int((~high_risk_mask).sum()),
            'p_value': lr_result.p_value,
            'test_statistic': lr_result.test_statistic,
            'dii': dii,
            'risk_groups': risk_groups,
            'kmf_high': kmf_high,
            'kmf_low': kmf_low,
            'survival_time': T,
            'survival_event': E,
        }
        
        return self.results
    
    def summary(self) -> str:
        """Return formatted summary of analysis results."""
        if self.results is None:
            return "No analysis run yet. Call analyze() first."
        
        r = self.results
        return f"""
FRANCIS Survival Analysis Results
=================================
Cancer Type: {r['cancer_type']}
Biomarker Triplet: {' + '.join(r['triplet'])}

Sample Size: {r['n_samples']}
Death Events: {r['n_events']}
High Risk Group: {r['n_high_risk']}
Low Risk Group: {r['n_low_risk']}

Log-Rank P-Value: {r['p_value']:.6f}
{'*** SIGNIFICANT ***' if r['p_value'] < 0.05 else ''}
"""

    def plot_km(self, ax=None, title=None):
        """
        Plot Kaplan-Meier survival curves.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to plot on
        title : str, optional
            Plot title
        """
        import matplotlib.pyplot as plt
        
        if self.results is None:
            raise ValueError("No analysis run yet. Call analyze() first.")
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
        
        r = self.results
        
        r['kmf_high'].plot_survival_function(ax=ax, ci_show=True, color='red', linestyle='--')
        r['kmf_low'].plot_survival_function(ax=ax, ci_show=True, color='blue', linestyle='-')
        
        if title is None:
            title = f"{r['cancer_type']}: Survival by DII ({' + '.join(r['triplet'])})"
        
        ax.set_title(title)
        ax.set_xlabel('Time (Days)')
        ax.set_ylabel('Survival Probability')
        ax.legend(loc='lower left')
        
        # Add p-value annotation
        ax.text(
            0.95, 0.95,
            f"Log-rank p = {r['p_value']:.4f}",
            transform=ax.transAxes,
            ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )
        
        return ax


def validate_all_cancers(
    expression_data: pd.DataFrame,
    survival_time: pd.Series,
    survival_event: pd.Series,
    aneuploidy_scores: pd.Series,
    cancer_labels: pd.Series
) -> pd.DataFrame:
    """
    Run FRANCIS validation across all 7 cancer types.
    
    Parameters
    ----------
    expression_data : pd.DataFrame
        Combined gene expression matrix
    survival_time : pd.Series
        Overall survival time
    survival_event : pd.Series
        Event indicator
    aneuploidy_scores : pd.Series
        Aneuploidy scores
    cancer_labels : pd.Series
        Cancer type labels (TCGA codes)
        
    Returns
    -------
    pd.DataFrame
        Validation results for all cancers
    """
    results = []
    
    for cancer_type in VALIDATED_TRIPLETS.keys():
        # Filter to this cancer
        cancer_mask = cancer_labels == cancer_type
        cancer_samples = cancer_labels[cancer_mask].index
        
        if len(cancer_samples) < 20:
            print(f"Skipping {cancer_type}: insufficient samples ({len(cancer_samples)})")
            continue
        
        try:
            model = FRANCISSurvival(cancer_type)
            r = model.analyze(
                expression_data.loc[cancer_samples],
                survival_time.loc[cancer_samples],
                survival_event.loc[cancer_samples],
                aneuploidy_scores.loc[cancer_samples]
            )
            
            results.append({
                'Cancer': cancer_type,
                'Triplet': ' + '.join(r['triplet']),
                'N': r['n_samples'],
                'Deaths': r['n_events'],
                'P-Value': r['p_value'],
                'Significant': r['p_value'] < 0.05
            })
            
            print(f"{cancer_type}: p={r['p_value']:.6f} (n={r['n_samples']}, deaths={r['n_events']})")
            
        except Exception as e:
            print(f"Error processing {cancer_type}: {e}")
    
    return pd.DataFrame(results)


if __name__ == '__main__':
    # Example usage
    print("FRANCIS Survival Analysis Module")
    print("=" * 40)
    print(f"Validated cancer types: {list(VALIDATED_TRIPLETS.keys())}")
    print(f"Universal marker: CXCL13 (present in 7/7 triplets)")
    print(f"Intermediate instability regime: aneuploidy {ANEUPLOIDY_LOW}-{ANEUPLOIDY_HIGH}")
