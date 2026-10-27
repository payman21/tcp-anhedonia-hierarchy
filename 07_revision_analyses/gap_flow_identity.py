"""Supplementary Methods 2: the hierarchical gap and the flow-asymmetry index are mathematically tied.

Summing the trophic-level equations (Lambda h = v) over one compartment X gives
    sum_{i in X, j in Y} (W_ij + W_ji)(h_i - h_j) = B - A
with A = X->Y and B = Y->X boundary weight, valid whenever X and Y partition the network.
This script checks the identity and how close the reported (unweighted) gap is to the flow index.
"""
import os, sys
import numpy as np
import scipy.io as sio
from scipy.stats import pearsonr, spearmanr, ttest_ind

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))   # repo root
from config import ATLAS_LABELS, GEC_OUTPUT

dl = sio.loadmat(f'{GEC_OUTPUT}/low_anhedonia/results_Ceff_low_anhedonia.mat')
dh = sio.loadmat(f'{GEC_OUTPUT}/high_anhedonia/results_Ceff_high_anhedonia.mat')
C = np.concatenate([dl['Ceff_LOW'], dh['Ceff_HIGH']]).astype(float)   # Ceff[target, source]
HL = np.concatenate([dl['hierarchicallevels_LOW'], dh['hierarchicallevels_HIGH']]).astype(float)
grp = np.r_[np.zeros(len(dl['Ceff_LOW'])), np.ones(len(dh['Ceff_HIGH']))]  # 0 = LOW, 1 = HIGH

lines = open(ATLAS_LABELS).read().strip().splitlines()
lab = [lines[i] for i in range(0, len(lines), 2)]

# --- Compartments (same definitions as 05_analysis/trophic_levels_main_analysis.ipynb) ---
sub = np.arange(0, 32)
cor = np.arange(32, 232)
mpfc = [i for i, l in enumerate(lab) if 'Default_PFC' in l and 'pCun' not in l]
ofc = [i for i, l in enumerate(lab) if 'Limbic_OFC' in l or 'Cont_OFC' in l]
acc = [i for i, l in enumerate(lab) if 'SalVentAttn_Med' in l or 'Cont_Cing' in l]
dlpfc = [i for i, l in enumerate(lab) if 'Cont_PFCl' in l]
pfc = dlpfc + mpfc + ofc + acc
striatum = [8, 9, 24, 25, 14, 15, 30, 31, 12, 13, 28, 29]


def cohen_d(v):
    a, b = v[grp == 1], v[grp == 0]
    return (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)


def analyse(X, Y, name):
    raw, fi, wsum, wg, ug = [], [], [], [], []
    for s in range(len(C)):
        Ce, h = C[s], HL[s]
        W = Ce.T                                  # W[i, j] = i -> j (as in compute_hierarchy)
        A = Ce[np.ix_(Y, X)].sum()                # X -> Y (top-down)
        B = Ce[np.ix_(X, Y)].sum()                # Y -> X (bottom-up)
        w = W[np.ix_(X, Y)] + W[np.ix_(Y, X)].T   # symmetric cross-edge weight
        d = h[X][:, None] - h[Y][None, :]
        raw.append(A - B)
        fi.append((A - B) / (A + B))
        wsum.append((w * d).sum())
        wg.append((w * d).sum() / w.sum())
        ug.append(h[X].mean() - h[Y].mean())
    raw, fi, wsum, wg, ug = map(np.array, (raw, fi, wsum, wg, ug))

    print(f'\n===== {name} =====')
    print(f'max |weighted level sum + raw (A-B)|  = {np.abs(wsum + raw).max():.1e}')
    print(f'max |weighted gap + normalised index| = {np.abs(wg + fi).max():.1e}')
    print(f'r(weighted gap, flow index)           = {pearsonr(wg, fi)[0]:+.6f}')
    r = pearsonr(ug, fi)[0]
    print(f'r(unweighted gap, flow index)         = {r:+.4f} (Spearman {spearmanr(ug, fi)[0]:+.4f}, R2 = {r * r:.3f})')
    print(f'r(unweighted gap, raw A-B)            = {pearsonr(ug, raw)[0]:+.4f}')
    print(f'sign agreement (gap vs -index)        = {np.mean(np.sign(ug) == np.sign(-fi)) * 100:.1f}%')
    for nm, v in [('flow index', fi), ('weighted gap', wg), ('unweighted gap', ug)]:
        t = ttest_ind(v[grp == 1], v[grp == 0])
        print(f'  {nm:15s} LOW {v[grp == 0].mean():+.4f}  HIGH {v[grp == 1].mean():+.4f}  '
              f'd = {cohen_d(v):+.3f}  t = {t.statistic:+.2f}  p = {t.pvalue:.1e}')


analyse(cor, sub, 'Cortex vs subcortex (full bipartition)')
analyse(pfc, striatum, 'PFC vs striatum (partial: not a bipartition)')

rest = np.setdiff1d(np.arange(232), pfc)
share = np.array([(C[s][np.ix_(striatum, pfc)].sum() + C[s][np.ix_(pfc, striatum)].sum())
                  / (C[s][np.ix_(rest, pfc)].sum() + C[s][np.ix_(pfc, rest)].sum()) for s in range(len(C))])
print(f'\nPFC<->striatum share of all PFC boundary weight: {share.mean() * 100:.1f}%')
