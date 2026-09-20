# Directed brain hierarchy and transdiagnostic anhedonia

Analysis code for the manuscript *"[Manuscript title]"*.

This repository contains the custom code used to derive transdiagnostic anhedonia
subgroups, estimate subject-level generative effective connectivity (GEC) from
resting-state fMRI, decompose that connectivity into trophic levels, and produce
every statistic and figure reported in the paper.

## Data

Analyses use the **Transdiagnostic Connectome Project** (TCP; NDA Study #3397,
https://doi.org/10.15154/w59r-ef41), distributed through the NIMH Data Archive
within NDA Collection 3552. The data are available to qualified investigators who
hold an executed NDA Data Use Certification. No participant data are redistributed here.

This repository is **code only**. In particular, the subject-level GEC matrices
produced by stage 4 are not included, because of their size (~190 MB); they are
available from the corresponding author on reasonable request. Stages 1–4 can be run
from the TCP source data; stage 5 requires the stage-4 outputs.

## Pipeline

Stages run in order. Each consumes the previous stage's output.

| Stage | Language | What it does |
|---|---|---|
| `01_phenotype/` | R | Multiple imputation (MICE) of SHAPS, TEPS, DASS, MADRS and QIDS item responses, excluding participants with >60% missingness within a scale. |
| `02_clustering/` | Python | Polychoric factor analysis of the 37 anhedonia items, Horn's parallel analysis, principal axis factoring with varimax rotation, and Gaussian mixture clustering into low- and high-anhedonia subgroups. |
| `03_parcellation/` | Python | Parcellation of HCP minimally-preprocessed, ICA-FIX–denoised CIFTI resting-state runs into 232 regions (Schaefer-200 7-network cortical + Tian Scale-II subcortical), concatenated run1-AP, run1-PA, run2-AP, run2-PA. |
| `04_gec_matlab/` | MATLAB | Hopf whole-brain model fitting. Estimates each subject's directed GEC matrix against empirical FC and time-shifted covariance, then computes trophic levels and trophic coherence. |
| `05_analysis/` | Python | Covariate-adjusted OLS models, cortico-subcortical and frontostriatal flow asymmetry, dimensional severity correlations, tiered BH-FDR correction, and the undirected functional-connectivity control analysis. |
| `06_figures/` | Python / MATLAB | Cortical and subcortical surface renderings of regional trophic levels. |

### Stage detail

**01_phenotype** — `impute_shaps_mice.Rmd` is the primary scale; the remaining four
scripts follow the same MICE procedure per scale. `post_imputation_processing.ipynb`
assembles the imputed item tables.

**02_clustering** — `anhedonia_factor_analysis_clustering.ipynb` produces the four-factor
solution and the GMM cluster assignments that serve as the grouping variable for all
imaging analyses. `clustering_validation.ipynb` runs bootstrap stability, permutation
distinctness and the published factor-loading figure.

**04_gec_matlab** — run `compute_peak_frequencies.m` before `estimate_gec_trophic.m`.
The latter fits a group-level solution first, then seeds each subject's fit from it.
`global_measure_stats.m` and `regional_trophic_stats.m` run the MATLAB-side group tests.

**05_analysis** — `trophic_levels_main_analysis.ipynb` is the main post-analysis notebook.
`multiple_comparisons_fdr.ipynb` implements the tiered FDR scheme described in the Methods.

## Requirements

**Python** ≥ 3.10 — numpy, scipy, pandas, matplotlib, seaborn, statsmodels,
scikit-learn, nibabel, factor_analyzer, semopy, ptitprince, palettable.
See `environment.yml`.

**R** ≥ 4.2 — mice, dplyr, tidyr, tidyverse, ggplot2, MASS.

**MATLAB** — R2021b or later, with the Parallel Computing Toolbox (the subject loops use
`parfor`) and the Statistics and Machine Learning Toolbox (`ranksum`).

### External dependencies not included here

These are third-party tools under their own licenses. Install separately and add to the
MATLAB path for stage 6.

- [BrainNet Viewer](https://www.nitrc.org/projects/bnv/) — subcortical ROI rendering
- [GIfTI toolbox](https://github.com/gllmflndn/gifti) and `xmltree` — surface file I/O
- `subtightplot` and `othercolor` (MATLAB File Exchange) — figure layout and colormaps

### Atlases and connectome

Not redistributed here. The Schaefer-200 7-network and Tian Scale-II atlases are publicly
available from their original sources. The anatomical scaffold for the Hopf model is the
normative connectome distributed with [Lead-DBS](https://www.lead-dbs.org/).

## Configuration

No absolute paths are hard-coded. All filesystem locations live in two files at the
repository root, `config.py` (Python) and `config.m` (MATLAB), and each can be overridden
by an environment variable.

| Variable | What it points at |
|---|---|
| `TCP_DATA` | `tcp_parcellations/data` — parcellated time series, atlases, structural connectome, cluster `.mat` files |
| `TCP_PHENO` | MICE-imputed item-level questionnaire tables (stage 01 output) |
| `TCP_GEC` | Stage 04 MATLAB outputs and the demographics table |
| `TCP_PROJ` | Project root, for demographics and supplementary inputs |
| `TCP_FIGS` | Where figures are written |
| `TCP_RENDER_UTILS` | Third-party MATLAB rendering toolboxes (stage 06) |

R scripts in stage 01 read `TCP_PHENO_RAW` for the directory holding the raw phenotype CSVs.

Notebooks pick this up with a bootstrap cell at the top:

```python
import os, sys
sys.path.insert(0, os.path.abspath('..'))   # repo root
from config import TCP_DATA, PHENO, GEC, PROJ, FIGS
```

MATLAB scripts call `cfg = config();` and read `cfg.TCP_DATA`, `cfg.CLUSTERS`, `cfg.SC_FILE`,
`cfg.OUTPUT` and `cfg.RENDER_UTILS`.

## Notebook state

Notebooks are stored without executed outputs. Stages 01 through 04 can be run from the TCP
source data; stage 05 additionally requires the stage-04 GEC matrices, which are not
redistributed here.

## License

GPL-3.0. See `LICENSE`.
