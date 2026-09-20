"""Filesystem locations for the analysis pipeline.

Every path below points at data that is NOT distributed with this repository
(see README). Set them for your own machine, either by editing the defaults
here or by exporting the matching environment variable.

    export TCP_DATA=/path/to/tcp_parcellations/data
    export TCP_PHENO=/path/to/phenotype/imputed
    export TCP_GEC=/path/to/trophic_coherence_analysis
    export TCP_PROJ=/path/to/TCP_thesis_project
    export TCP_FIGS=/path/to/figure/output
"""
import os

# Parcellated time series, atlases, structural connectome, cluster .mat files
TCP_DATA = os.environ.get("TCP_DATA", "/path/to/tcp_parcellations/data")

# MICE-imputed item-level questionnaire tables (stage 01 output)
PHENO = os.environ.get("TCP_PHENO", "/path/to/phenotype/imputed")

# Stage 04 MATLAB outputs (per-group results_Ceff_*.mat) and demographics
GEC = os.environ.get("TCP_GEC", "/path/to/trophic_coherence_analysis")

# Project root (demographics, supplementary inputs, figure output)
PROJ = os.environ.get("TCP_PROJ", "/path/to/TCP_thesis_project")

# Where figures are written
FIGS = os.environ.get("TCP_FIGS", "./figures")

ATLAS_LABELS = os.path.join(
    TCP_DATA, "atlases",
    "Schaefer2018_200Parcels_7Networks_order_Tian_Subcortex_S2_label.txt")
