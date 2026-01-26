# FRANCIS - Cancer Dynamics Model

**Framework for Real-time Analysis of Nonlinear Critical Instability Signals**

Immunometabolic dynamics and critical transition detection in cancer using Dynamical Network Biomarker (DNB) theory.

## Overview

This repository implements computational methods for detecting critical transitions in cancer progression using nonlinear dynamical systems analysis. The system identifies pre-metastatic tipping points by computing a **Dynamical Instability Index (DII)** from longitudinal biomarker data.

### Key Features

- **DNB Analysis**: Identifies dominant gene modules exhibiting critical slowing down
- **DII Computation**: Integrates entropy, divergence, and trajectory instability metrics
- **Immunometabolic Markers**: Pre-configured panels for immune checkpoint and metabolic drivers
- **Critical Transition Detection**: Early warning of physiological collapse

## Quickstart

```bash
# Clone repository
git clone https://github.com/FRANCIS-bio/cancer-dynamics-model.git
cd cancer-dynamics-model

# Install dependencies
pip install -r requirements.txt

# Run demo
python dnb_cancer.py
```

## What It Detects

- Pre-metastatic tipping points
- Immunometabolic collapse signatures
- Loss of coordinated immune-metabolic control
- Critical transition windows for intervention timing

## Immunometabolic Markers

The system analyzes key marker panels:

| Category | Markers |
|----------|---------|
| Immune Checkpoint | PDCD1, LAG3, CTLA4, CD274 |
| Metabolic Drivers | KRAS, LDHA, CPT1A, FABP4 |
| Warburg Effect | LDHA, HK2, PDK1, PKM |
| Inflammatory | IL6, TNF, STAT3, NFKB1 |

## The Science

### Dynamical Network Biomarkers (DNB)

As biological systems approach critical transitions, a dominant group of molecules exhibits:
1. **Dramatically increased variance**
2. **Increased internal correlation** (within the DNB group)
3. **Decreased external correlation** (with other molecules)

This is the hallmark of **critical slowing down** preceding state transitions.

### Dynamical Instability Index (DII)

```
DII(t) = w₁ × Entropy(t) + w₂ × Divergence(t) + w₃ × TrajectoryInstability(t)
```

Where:
- **Entropy**: Sample/multiscale entropy measuring complexity breakdown
- **Divergence**: Jensen-Shannon divergence from healthy reference
- **Trajectory Instability**: Lyapunov exponent approximation

## Usage Example

```python
from dnb_cancer import detect_critical_transition, compute_dii

# Load your gene expression time-series data
# data shape: (n_genes, n_timepoints)

# Detect critical transition
dnb_genes, tipping_score, is_critical = detect_critical_transition(data)

if is_critical:
    print(f"WARNING: Critical transition detected (score: {tipping_score:.2f})")
    print(f"Top DNB genes: {dnb_genes[:10]}")

# Compute DII
dii, components = compute_dii(data)
print(f"DII: {dii:.4f}")
```

## File Structure

```
cancer-dynamics-model/
├── dnb_cancer.py              # Core DNB/DII computation engine
├── immunometabolic_markers.json   # Marker panel configurations
├── requirements.txt           # Python dependencies
├── README.md                  # This file
└── LICENSE                    # Apache 2.0
```

## Validation

Methods validated against:
- TCGA pancreatic adenocarcinoma (PAAD): t=28.266, p=3.81e-92
- TCGA ovarian cancer (OV): t=-18.550, p=3.66e-59
- Correlation between immune/metabolic axes: r=0.769, p=1.15e-272

## References

1. Chen L, et al. (2012). Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. *Scientific Reports*.
2. Scheffer M, et al. (2009). Early-warning signals for critical transitions. *Nature*.

## License

Apache 2.0 — See [LICENSE](LICENSE) file.

## Patent Notice

The methods implemented in this repository are covered by:

**U.S. Provisional Patent Application No. 63/968,064**  
*"Systems and Methods for Predicting and Controlling Physiological Collapse Using Dynamical Biomarker Indices"*  
Filed: January 26, 2026

This open-source release is for **research and educational use**. Commercial use of the patented system/method may require a license from the patent holder.

## Contact

For licensing inquiries or collaboration, open an issue or contact via GitHub.

---

*FRANCIS: Detecting critical transitions before collapse.*
