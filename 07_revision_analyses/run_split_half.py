"""Fit a GEC separately on run-1 (AP+PA) and run-2 (AP+PA) for every subject.

Seeds and intrinsic frequencies are the same group-level ones used in the main analysis (stage 04
outputs, read-only). Resumable: one npz per subject in SPLIT_HALF (subject-level, never committed).
Usage: python run_split_half.py [n_per_group_for_pilot] [workers]"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'; os.environ['OPENBLAS_NUM_THREADS'] = '1'; os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
import numpy as np, scipy.io as sio
from multiprocessing import Pool
from hopf_lna import *
from config import SPLIT_HALF

HALVES = {'half1': RUNS[:2], 'half2': RUNS[2:]}
SC = load_sc(); MASK = update_mask(SC)
G = {}
for g in ['low', 'high']:
    G[g] = dict(seed=sio.loadmat(f'{GEC_OUTPUT}/{g}_anhedonia/results_Ceff_{g}_anhedonia.mat')[f'Ceffgroup_{g.upper()}'],
                f=sio.loadmat(f'{GEC_OUTPUT}/empirical_{g}_anhedonia.mat')['f_diff'],
                ids=sio.loadmat(f'{CLUSTERS}/{g}_anhedonia.mat')[f'{g}_anhedonia'])

def job(args):
    g, k = args
    sid = G[g]['ids'][k, 1].flat[0]; path = f'{SPLIT_HALF}/{g}_{k:03d}_{sid}.npz'
    if os.path.exists(path): return path
    out = {}
    for h, runs in HALVES.items():
        FC, CV = empirical(load_runs(sid, runs))
        C, rf, rc, it = fit(G[g]['seed'], FC, CV, G[g]['f'], MASK)
        out[h] = C; out[h + '_fit'] = [rf, rc, it]
    np.savez(path, sid=sid, group=g, **out)
    return path

if __name__ == '__main__':
    os.makedirs(SPLIT_HALF, exist_ok=True)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else None
    w = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    jobs = [(g, k) for g in ['low', 'high'] for k in range(len(G[g]['ids']))]
    if n:   # pilot: first n of each group
        jobs = [(g, k) for g, k in jobs if k < n]
    with Pool(w) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs)):
            print(i + 1, '/', len(jobs), r, flush=True)
