"""Demo: trophic levels and coherence from the group-level GEC matrices.

Runs on the files in group_level_data/ only, so it needs no access to the TCP source data.
Python port of compute_hierarchy() in 04_gec_matlab/estimate_gec_trophic.m.
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DATA = os.path.join(ROOT, 'group_level_data')
SUB, COR = np.arange(0, 32), np.arange(32, 232)   # 0-based: Tian subcortex, Schaefer cortex


def compute_hierarchy(Ceff):
    A = Ceff.T                                    # A[i, j] = i -> j
    d, delta = A.sum(0), A.sum(1)
    u, v = d + delta, d - delta
    Lam = np.diag(u) - A - A.T
    Lam[0, 0] = 0                                 # regularisation, as in the MATLAB code
    h = np.linalg.solve(Lam, v)
    h -= h.min()
    F0 = (A * (h[:, None] - h[None, :] - 1) ** 2).sum() / A.sum()
    return 1 - F0, h


regional = pd.read_csv(os.path.join(DATA, 'regional_trophic_levels_group_means.csv'))
print(f'{"group":6s} {"coherence":>10s} {"cortex":>8s} {"subcortex":>10s} {"gap":>8s} {"flow index":>11s}')
for grp in ['low', 'high']:
    C = np.loadtxt(os.path.join(DATA, f'Ceffgroup_{grp}_anhedonia.csv'), delimiter=',')
    tc, h = compute_hierarchy(C)
    down = C[np.ix_(SUB, COR)].sum()              # cortex -> subcortex (Ceff[target, source])
    up = C[np.ix_(COR, SUB)].sum()                # subcortex -> cortex
    print(f'{grp:6s} {tc:10.4f} {h[COR].mean():8.4f} {h[SUB].mean():10.4f} '
          f'{h[COR].mean() - h[SUB].mean():8.4f} {(down - up) / (down + up):+11.4f}')

d = regional['diff_high_minus_low'].to_numpy()
print(f'\nSubject-averaged regional levels (high - low): cortex {d[COR].mean():+.4f}, '
      f'subcortex {d[SUB].mean():+.4f}')
