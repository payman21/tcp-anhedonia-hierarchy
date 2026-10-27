"""Supplementary Table 5: split-half reliability and within-half replication of the main measures.

Uses the full-data GECs (stage 04) and the split-half GECs written by run_split_half.py (one GEC on
run-01 AP+PA, one on run-02 AP+PA, per participant; same group seeds and intrinsic frequencies as the
main analysis). Reports Pearson r between halves, the Spearman-Brown-corrected r, ICC(A,1), and the
covariate-adjusted group effect fitted within each half and on the full data.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import scipy.io as sio
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))   # repo root
from config import ATLAS_LABELS, CLUSTERS, GEC_OUTPUT, DEMOS, SPLIT_HALF

lines = open(ATLAS_LABELS).read().strip().splitlines()
lab = [lines[i] for i in range(0, len(lines), 2)]
idx = lambda f: [i for i, l in enumerate(lab) if f(l)]
SUB, COR = np.arange(32), np.arange(32, 232)
mPFC = idx(lambda l: 'Default_PFC' in l and 'pCun' not in l)
OFC = idx(lambda l: 'Limbic_OFC' in l or 'Cont_OFC' in l)
ACC = idx(lambda l: 'SalVentAttn_Med' in l or 'Cont_Cing' in l)
dlPFC = idx(lambda l: 'Cont_PFCl' in l)
NAc, aPUT = [8, 9, 24, 25], [12, 28]
PFC = dlPFC + mPFC + OFC + ACC
STR = NAc + [14, 15, 30, 31] + [12, 13, 28, 29]
EDGES = {'vmPFC_NAc': (mPFC, NAc), 'OFC_NAc': (OFC, NAc), 'ACC_NAc': (ACC, NAc),
         'dlPFC_NAc': (dlPFC, NAc), 'vmPFC_aPUT': (mPFC, aPUT), 'ACC_aPUT': (ACC, aPUT)}


def troph(Ce):
    A = Ce.T
    d, dl = A.sum(0), A.sum(1)
    L = np.diag(d + dl) - A - A.T
    L[0, 0] = 0
    h = np.linalg.solve(L, d - dl)
    h -= h.min()
    F0 = np.sum(A * (h[:, None] - h[None, :] - 1) ** 2) / A.sum()
    return h, 1 - F0


def measures(Ce):
    """Ceff[target, source]."""
    h, tc = troph(Ce)
    m = {'gap_cs': h[COR].mean() - h[SUB].mean(), 'gap_fs': h[PFC].mean() - h[STR].mean(), 'TC': tc}
    for nm, X, Y in [('Q_cs', COR, SUB), ('Q_fs', PFC, STR)]:
        td, bu = Ce[np.ix_(Y, X)].sum(), Ce[np.ix_(X, Y)].sum()
        m[nm] = (td - bu) / (td + bu)
    for e, (P, S) in EDGES.items():
        td, bu = Ce[np.ix_(S, P)].mean(), Ce[np.ix_(P, S)].mean()
        m[e] = (td - bu) / (td + bu)
    return m


def icc_a1(a, b):
    """ICC(A,1): two-way, absolute agreement, single measure."""
    Y = np.c_[a, b]
    n, k = Y.shape
    gm = Y.mean()
    msr = k * np.sum((Y.mean(1) - gm) ** 2) / (n - 1)
    msc = n * np.sum((Y.mean(0) - gm) ** 2) / (k - 1)
    sse = np.sum((Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm) ** 2)
    mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)


# ---------------------------------------------------------------- data
rows, full = [], {}
for g in ['low', 'high']:
    F = sio.loadmat(f'{GEC_OUTPUT}/{g}_anhedonia/results_Ceff_{g}_anhedonia.mat')[f'Ceff_{g.upper()}'].astype(float)
    ids = sio.loadmat(f'{CLUSTERS}/{g}_anhedonia.mat')[f'{g}_anhedonia']
    for k in range(len(ids)):
        full[ids[k, 1].flat[0]] = (g, F[k])
for p in sorted(glob.glob(f'{SPLIT_HALF}/*.npz')):
    z = np.load(p, allow_pickle=True)
    sid = str(z['sid'])
    r = {'id': sid, 'grp': int(str(z['group']) == 'high')}
    for tag, C in [('h1', z['half1']), ('h2', z['half2']), ('full', full[sid][1])]:
        for k, v in measures(C.astype(float)).items():
            r[f'{k}_{tag}'] = v
    rows.append(r)
df = pd.DataFrame(rows)
dm = pd.read_excel(DEMOS, header=1).set_index('subjectkey')
dm['age'] = dm['Age'].replace(999, np.nan)
dm['pat'] = (dm['Primary_Dx_Payman'] != 999).astype(int)
df = df.join(dm[['age', 'sex', 'pat']], on='id')
dc = df.dropna(subset=['age'])
dc = dc[dc.sex.isin(['M', 'F'])]
print(f'N = {len(df)} (low {(df.grp == 0).sum()}, high {(df.grp == 1).sum()}); covariate models N = {len(dc)}')

# ---------------------------------------------------------------- reliability
MEAS = ['gap_cs', 'Q_cs', 'gap_fs', 'Q_fs', 'TC'] + list(EDGES)
print('\n=== Split-half reliability (run-01 vs run-02) ===')
print(f'  {"measure":11s} {"r":>6s} {"SB r":>6s} {"ICC(A,1)":>8s} {"r low":>6s} {"r high":>6s} {"r half~full":>11s}')
for m in MEAS:
    a, b = df[f'{m}_h1'], df[f'{m}_h2']
    r = np.corrcoef(a, b)[0, 1]
    rl = np.corrcoef(a[df.grp == 0], b[df.grp == 0])[0, 1]
    rh = np.corrcoef(a[df.grp == 1], b[df.grp == 1])[0, 1]
    rf = np.mean([np.corrcoef(df[f'{m}_{t}'], df[f'{m}_full'])[0, 1] for t in ['h1', 'h2']])
    print(f'  {m:11s} {r:6.3f} {2 * r / (1 + r):6.3f} {icc_a1(a, b):8.3f} {rl:6.3f} {rh:6.3f} {rf:11.3f}')

# ---------------------------------------------------------------- within-half replication
print('\n=== Group effect replicated in each half: OLS ~ grp + age + sex + patient [β (p)] ===')
for m in MEAS:
    s = []
    for t in ['h1', 'h2', 'full']:
        f = smf.ols(f'{m}_{t} ~ grp + age + sex + pat', dc).fit()
        s.append(f'{t}: {f.params["grp"]:+.4f} ({f.pvalues["grp"]:.2g})')
    print(f'  {m:11s} ' + ' | '.join(s))
