# FRANCIS

Detecting cancer tipping points before they happen.

## What this is

Code that implements Dynamical Network Biomarker (DNB) theory to detect when a biological system is approaching a critical transition. Built for cancer, but the math works for any complex system approaching a tipping point.

## Status

Code: Working. Math is correct.

Theory: Based on published research (Chen et al. 2012, cited 300+ times).

Testing: Simulated data only. Clinical validation on real patient data not yet done.

## How it works

The code computes a Dynamical Instability Index (DII) from biomarker time series:

DII = Entropy + Divergence + Trajectory Instability

When DII spikes, the system is approaching a tipping point. That's the window for intervention.

DNB theory says that before a system tips, a group of genes will show:
1. Increased variance
2. Increased correlation with each other
3. Decreased correlation with everything else

This is called critical slowing down. The code detects it.

## What the simulation showed

```
Tipping score: 0.351 (below threshold, system stable)
DII: 0.328
Top markers: IDO1, CD274 (PD-L1), HIF1A, LDHA
```

These are real immunotherapy targets. The code identified them correctly from simulated data designed to mimic healthy dynamics.

## What it doesn't show yet

Real patient validation. The simulation proves the code runs. It doesn't prove the method detects real cancer tipping points in real patients.

That requires longitudinal data (multiple measurements over time from the same patients as they progress). That validation is next.

## Run it

```
git clone https://github.com/FRANCIS-bio/cancer-dynamics-model.git
cd cancer-dynamics-model
pip install -r requirements.txt
python dnb_cancer.py
```

## Files

dnb_cancer.py — the engine (DNB scores, DII computation)

immunometabolic_markers.json — marker panel configs

demo_cancer_dnb.py — visualization demo

requirements.txt — dependencies

## The science

DNB theory: Chen L, et al. (2012). Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. Scientific Reports.

Critical transitions: Scheffer M, et al. (2009). Early-warning signals for critical transitions. Nature.

Both are real, peer-reviewed, widely cited.

## License

Apache 2.0. Use it. Build on it.

## Patent

U.S. Provisional Patent Application No. 63/968,064

Systems and Methods for Predicting and Controlling Physiological Collapse Using Dynamical Biomarker Indices

Filed January 26, 2026

The patent protects the method. The code is open for research. Commercial licensing available.

## What's next

1. Validate on TCGA longitudinal proxy data (stage progression as pseudo-time)
2. Test on liquid biopsy time series if available
3. Publish results when real validation is complete

## Contact

Open an issue or reach out through GitHub.
