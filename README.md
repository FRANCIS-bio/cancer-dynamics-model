# FRANCIS

**Forecasting Risk Analysis Networks for Clinical Intervention Systems**

*Detecting cancer tipping points before they happen.*

---

> *In memory of Uncle Frankie.*  
> *This work exists because early warning matters.*

---

## Overview

The FRANCIS framework identifies high-risk cancer patients within an intermediate genomic instability regime — a "Snap-Back Trap" where conventional diagnostics fail. The system computes a Dynamical Immunometabolic Index (DII) from tissue-specific biomarker triplets anchored by CXCL13.

**Validation:** 2,415 patients across 7 solid tumor types. All p-values < 0.05 by log-rank survival analysis.

**Key Discovery:** CXCL13 appears in 100% of optimal triplets (7/7 cancers), establishing it as a universal marker of tumor-immune dynamics within the intermediate instability regime.

---

## Status

| Component | Status |
|-----------|--------|
| Theory | Published (Chen et al. 2012, Scheffer et al. 2009) |
| Code | Operational |
| Clinical Validation | **Complete** — 7 cancers, 2,415 patients |
| Patents | 2 U.S. Provisional Applications Filed |

---

## The Intermediate Instability Trap

Existing diagnostics assume monotonic dynamics: higher genomic instability correlates with worse prognosis.

This assumption is incorrect in the intermediate regime.

Patients with aneuploidy scores between 5 and 25 exhibit paradoxically *worse* survival than both low-instability and high-instability cohorts. Such tumors occupy a metastable state — sufficiently unstable to evade treatment, insufficiently unstable to trigger immune recognition or mitotic catastrophe.

The FRANCIS framework identifies these patients. Conventional assays do not.

---

## Validation Results

**Source:** TCGA PanCancer Atlas  
**Method:** Log-rank survival analysis, median DII stratification  
**Filter:** Intermediate genomic instability regime (aneuploidy 5-25)

| Cancer | TCGA | Origin | Biomarker Triplet | N | Deaths | P-Value |
|--------|------|--------|-------------------|---|--------|---------|
| Glioblastoma | GBM | Neuroectoderm | CXCL13 + SOX2 + LGALS1 | 134 | 102 | 1.6×10⁻⁵ |
| Melanoma | SKCM | Neural Crest | FABP4 + CXCL13 + CD8A | 276 | 157 | 0.0007 |
| Lung Adenocarcinoma | LUAD | Epithelial | CXCL13 + GAPDH + KRAS | 392 | 135 | 0.0019 |
| Colorectal | COADREAD | Epithelial | CXCL13 + TOP2A + KRAS | 402 | 84 | 0.0046 |
| Breast Carcinoma | BRCA | Epithelial | CXCL13 + ARG1 + CD8A | 790 | 114 | 0.0053 |
| Sarcoma | SARC | Mesenchymal | CXCL13 + CD8A + PDCD1 | 169 | 71 | 0.0117 |
| Ovarian | OV | Epithelial | LDHA + MKI67 + CXCL13 | 252 | 149 | 0.0152 |

**Total:** 2,415 patients | 812 death events | 7 cancers | 3 embryonic lineages

---

## Mechanism

In solid tumors, CXCL13 recruits B-lymphocytes that form tertiary lymphoid structures (TLS) within the tumor microenvironment. These structures enable coordinated anti-tumor immune response.

Low CXCL13 expression within the intermediate instability regime correlates with:
- Impaired B-cell recruitment
- Reduced TLS formation
- Diminished antigen presentation
- Poor overall survival

This mechanism is specific to solid tumors. In B-cell malignancies, CXCL13 operates through an opposing autocrine pathway and is explicitly excluded from the present framework.

---

## The Index

```
DII = mean(z₁, z₂, z₃)
```

Where z₁, z₂, z₃ represent z-score normalized expression values of the tissue-specific biomarker triplet relative to a reference population.

**Procedure:**
1. Determine genomic instability score (aneuploidy)
2. Confirm intermediate regime (score 5-25)
3. Measure expression of tissue-specific triplet
4. Compute Dynamical Immunometabolic Index
5. Stratify: DII below median indicates high risk

---

## Repository Structure

```
cancer-dynamics-model/
├── config/                 # Biomarker configurations
│   └── immunometabolic_markers.json
├── data/                   # Data directory (not tracked)
├── docs/                   # Figures and documentation
│   └── FRANCIS_Patent_Figures.pdf
├── models/                 # Core analysis modules
│   ├── dnb_cancer.py
│   └── francis_survival.py
├── scripts/                # Validation and demo scripts
│   ├── tcga_validation.py
│   └── demo_cancer_dnb.py
├── LICENSE
├── README.md
└── requirements.txt
```

---

## Installation

```bash
git clone https://github.com/FRANCIS-bio/cancer-dynamics-model.git
cd cancer-dynamics-model
pip install -r requirements.txt
```

## Usage

```python
from models.francis_survival import FRANCISSurvival

model = FRANCISSurvival(cancer_type='BRCA')
results = model.analyze(expression_data, survival_time, survival_event, aneuploidy)

print(results['p_value'])  # 0.0053
```

---

## References

**Dynamical Network Biomarkers:**  
Chen L, et al. (2012). Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers. *Scientific Reports.*

**Critical Transitions:**  
Scheffer M, et al. (2009). Early-warning signals for critical transitions. *Nature.*

**Tertiary Lymphoid Structures:**  
Sautès-Fridman C, et al. (2019). Tertiary lymphoid structures in the era of cancer immunotherapy. *Nature Reviews Cancer.*

---

## Patents

**U.S. Provisional Application No. 63/968,064**  
*Systems and Methods for Predicting and Controlling Physiological Collapse Using Dynamical Biomarker Indices*  
Filed January 26, 2026

**U.S. Provisional Application No. 63/970,059**  
*Methods, Systems, and Kits for Pan-Cancer Survival Stratification Using CXCL13-Anchored Biomarker Triplets Within Intermediate Genomic Instability Regimes*  
Filed January 28, 2026

The patents protect the methods. The code is released under Apache 2.0 for research purposes. Commercial licensing available.

---

## License

Apache 2.0

---

## Contact

Open an issue or reach out through GitHub.

---

*The system detects collapse before it becomes irreversible.*  
*That window is where intervention belongs.*
