# FRANCIS

Catching cancer before it tips.

## What this is

Cancer doesn't progress in a straight line. It hits tipping points, moments where the system becomes unstable and tips into something worse. Metastasis. Treatment failure. Collapse.

By the time standard biomarkers catch it, you're already there.

This code catches the approach to the tipping point. Not the event, the warning signs before it.

## The idea

Biological systems show specific signatures when they're approaching a critical transition:

1. Variance goes up. Things get wobbly before they break.
2. Internal correlations spike. A small group of genes start moving together.
3. External correlations drop. That group decouples from everything else.

This is called critical slowing down. It's physics. It happens in ecosystems, climate, finance, and it happens in tumors.

I built this to compute a Dynamical Instability Index (DII) from longitudinal biomarker data. When DII spikes, you're approaching a tipping point. That's your intervention window.

## What it actually detects

Pre-metastatic tipping points. Immunometabolic collapse, when immune checkpoints and Warburg metabolism decouple. The window where intervention can still work.

## The markers

Immune checkpoint: PDCD1, LAG3, CTLA4, CD274

Warburg metabolism: LDHA, HK2, PDK1, PKM

Inflammatory: IL6, TNF, STAT3, NFKB1

These aren't random. They're the axes that correlate across cancers. When that correlation structure breaks down, you're approaching the tipping point.

## Run it

```
git clone https://github.com/FRANCIS-bio/cancer-dynamics-model.git
cd cancer-dynamics-model
pip install -r requirements.txt
python dnb_cancer.py
```

## Usage

```
from dnb_cancer import detect_critical_transition, compute_dii

my_data = load_your_data()  # shape: (genes, timepoints)

genes, score, is_critical = detect_critical_transition(my_data)
dii, components = compute_dii(my_data)

if is_critical:
    print("Approaching tipping point")
```

Works on any cancer. Just plug in your expression data.

## The math

DII(t) = w1 × Entropy(t) + w2 × Divergence(t) + w3 × TrajectoryInstability(t)

Entropy is complexity breakdown. Divergence is how far current state is from healthy baseline. Trajectory instability is whether nearby states are diverging.

When these spike together, the system is losing stability.

## Validated against real data

TCGA Pancreatic (PAAD): t=30.0, p=2.09e-98

TCGA Ovarian (OV): t=23.7, p=2.37e-84

TCGA Prostate (PRAD): t=18.7, p=1.46e-60

TCGA Breast (BRCA): t=27.6, p=2.63e-130

TCGA Melanoma (SKCM): t=51.2, p=4.13e-285

Correlation: r=0.84 to 0.89 across all 5 cancer types

Total: 2,670 tumor samples, 981 normal samples

This isn't theoretical. It replicates.

## Files

dnb_cancer.py is the engine

multicancer_validation.py validates across 5 TCGA cancers

immunometabolic_markers.json has marker configs

demo_cancer_dnb.py runs the visualization

tcga_analysis.py does detailed TCGA analysis

requirements.txt has dependencies

## References

Chen et al. (2012). Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. Scientific Reports.

Scheffer et al. (2009). Early-warning signals for critical transitions. Nature.

## License

Apache 2.0. Use it. Build on it. Credit appreciated.

## Patent

U.S. Provisional Application No. 63/968,064

Systems and Methods for Predicting and Controlling Physiological Collapse Using Dynamical Biomarker Indices

Filed January 26, 2026

The code is open for research. If you want to build a commercial product around the method, let's talk.

## Contact

Open an issue or reach out through GitHub.
