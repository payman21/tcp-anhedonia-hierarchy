"""Python port of the linear-Hopf GEC fit in 04_gec_matlab/estimate_gec_trophic.m.

Used for the split-half refits. Defaults reproduce the MATLAB settings; imported by
run_split_half.py and analyze_split_half.py."""
import os, sys
import numpy as np
from scipy.signal import butter, filtfilt, detrend
from scipy.linalg import solve_continuous_lyapunov, expm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))   # repo root
from config import TCP_DATA, CLUSTERS, SC_FILE, GEC_OUTPUT

N, TAU, SIGMA, TR = 232, 3, 0.01, 0.8
EPS_FC, EPS_FCTAU, MAXC, A_BIF = 4e-4, 1e-4, 0.2, -0.02
fnq = 1 / (2 * TR)
BFILT, AFILT = butter(2, [0.01 / fnq, 0.99], btype='band')

RUNS = ['restAP_run-01', 'restPA_run-01', 'restAP_run-02', 'restPA_run-02']  # original concatenation order

def load_sc():
    import h5py
    with h5py.File(SC_FILE) as f:
        return np.array(f['SC']).T   # h5py reads MATLAB v7.3 transposed

def hom_index():
    hom = np.arange(N)
    hom[:16] += 16; hom[16:32] -= 16; hom[32:132] += 100; hom[132:] -= 100
    return hom

def update_mask(SC):
    C = SC / SC.max() * MAXC
    mask = C > 0
    mask[np.arange(N), hom_index()] = True
    return mask

def load_runs(sid, runs):
    """Parcellated time series (N x T) for one subject, runs concatenated in the given order."""
    return np.concatenate([np.loadtxt(f'{TCP_DATA}/parcellated_rsfmri/{sid}/parcellated_task-{r}_bold_Atlas_MSMAll_hp2000_clean.csv',
                                      delimiter='\t').T for r in runs], 1)

def empirical(ts):
    """ts: N x T. Returns FCemp, COVtau_emp (as process_subject_empirical)."""
    ts = ts.astype(float).copy()
    for i in range(N):
        r = detrend(ts[i] - np.nanmean(ts[i]))
        ts[i] = filtfilt(BFILT, AFILT, r, padlen=3 * (max(len(AFILT), len(BFILT)) - 1))
    ts2 = ts[:, 9:-10]            # MATLAB 10:end-10
    FC = np.corrcoef(ts2)
    Xt, Xtau = ts2[:, :-TAU], ts2[:, TAU:]
    Xt = Xt - Xt.mean(1, keepdims=True); Xtau = Xtau - Xtau.mean(1, keepdims=True)
    raw = Xt @ Xtau.T / Xt.shape[1]
    sd = 1 / np.sqrt(np.diag(np.cov(ts2)))
    return FC, raw * np.outer(sd, sd)

def hopf_model(C, f_diff):
    wo = 2 * np.pi * np.asarray(f_diff).ravel()
    Axx = A_BIF * np.eye(N) - np.diag(C.sum(1)) + C
    A = np.block([[Axx, -np.diag(wo)], [np.diag(wo), Axx]])
    K = solve_continuous_lyapunov(A, -(SIGMA ** 2) * np.eye(2 * N))
    d = np.sqrt(np.diag(K))
    FC = (K / np.outer(d, d))[:N, :N]
    CV = K[:N, :N]
    return FC, CV, K, A

def fit(C_init, FCemp, COVemp, f_diff, mask, max_iter=5000, verbose=False):
    Cn = C_init.copy(); old = 1e5
    tri = np.tril_indices(N, -1)
    for it in range(1, max_iter + 1):
        FCs, CV, K, A = hopf_model(Cn, f_diff)
        Ks = (expm(TAU * TR * A) @ K)[:N, :N]
        sd = 1 / np.sqrt(np.diag(CV))
        COVs = Ks * np.outer(sd, sd)
        dFC, dCOV = FCemp - FCs, COVemp - COVs
        if it % 100 == 0:
            err = np.mean(dFC ** 2) + np.mean(dCOV ** 2)
            if verbose: print(it, err)
            if (old - err) / err < 0.001 or old < err:
                break
            old = err
        Cn[mask] += (EPS_FC * dFC + EPS_FCTAU * dCOV)[mask]
        Cn[Cn < 0] = 0
        Cn = Cn / Cn.max() * MAXC
    return Cn, np.corrcoef(FCemp[tri], FCs[tri])[0, 1], np.corrcoef(COVemp[tri], COVs[tri])[0, 1], it

def troph(Ceff):
    A = Ceff.T
    d, dl = A.sum(0), A.sum(1)
    Lam = np.diag(d + dl) - A - A.T; Lam[0, 0] = 0
    g = np.linalg.solve(Lam, d - dl); return g - g.min()

SUB, COR = np.arange(32), np.arange(32, 232)
def measures(Ceff):
    A = Ceff.T
    td = A[np.ix_(COR, SUB)].sum(); bu = A[np.ix_(SUB, COR)].sum()
    h = troph(Ceff)
    return dict(td=td, bu=bu, Q=(td - bu) / (td + bu), gap=h[COR].mean() - h[SUB].mean())
