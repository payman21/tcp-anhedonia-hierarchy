"""Covariate-adjusted group means against zero, and standardized effect sizes d, for the
cortex-subcortex and PFC-striatum gaps and flow indices.

Uses the same data, definitions and OLS covariates as 05_analysis/trophic_levels_main_analysis.ipynb.
Adjusted means come from refitting the model with one coefficient per group and covariates centred
at their sample means. d is the covariate-adjusted group coefficient divided by the pooled
within-group standard deviation of the outcome.
"""
import os, sys
import numpy as np, pandas as pd, scipy.io as sio, statsmodels.formula.api as smf

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))   # repo root
from config import ATLAS_LABELS, CLUSTERS, GEC_OUTPUT, DEMOS

dl = sio.loadmat(f'{GEC_OUTPUT}/low_anhedonia/results_Ceff_low_anhedonia.mat')
dh = sio.loadmat(f'{GEC_OUTPUT}/high_anhedonia/results_Ceff_high_anhedonia.mat')
mat_low = sio.loadmat(f'{CLUSTERS}/low_anhedonia.mat')
mat_high = sio.loadmat(f'{CLUSTERS}/high_anhedonia.mat')
ids_low = [mat_low['low_anhedonia'][i, 1].flat[0] for i in range(mat_low['low_anhedonia'].shape[0])]
ids_high = [mat_high['high_anhedonia'][i, 1].flat[0] for i in range(mat_high['high_anhedonia'].shape[0])]
demos = pd.read_excel(DEMOS, header=1)
demos['Age_clean'] = demos['Age'].replace(999, np.nan)
demos['is_patient'] = demos['Primary_Dx_Payman'].apply(lambda x: 0 if x == 999 else 1)
lines = open(ATLAS_LABELS).read().strip().splitlines()
lab = [lines[i] for i in range(0, len(lines), 2)]

sub = np.arange(0, 32); cor = np.arange(32, 232)
mpfc = [i for i, l in enumerate(lab) if 'Default_PFC' in l and 'pCun' not in l]
ofc = [i for i, l in enumerate(lab) if 'Limbic_OFC' in l or 'Cont_OFC' in l]
acc = [i for i, l in enumerate(lab) if 'SalVentAttn_Med' in l or 'Cont_Cing' in l]
dlp = [i for i, l in enumerate(lab) if 'Cont_PFCl' in l]
pfc = dlp + mpfc + ofc + acc
st = [8, 9, 24, 25] + [14, 15, 30, 31] + [12, 13, 28, 29]

def measures(C, H):
    def fi(X, Y):
        A = np.array([c[np.ix_(Y, X)].sum() for c in C]); B = np.array([c[np.ix_(X, Y)].sum() for c in C])
        return (A - B) / (A + B)
    return {'ctx_sub_gap': H[:, cor].mean(1) - H[:, sub].mean(1),
            'pfc_str_gap': H[:, pfc].mean(1) - H[:, st].mean(1),
            'ctx_sub_flow': fi(cor, sub), 'pfc_str_flow': fi(pfc, st)}

mL = measures(dl['Ceff_LOW'].astype(float), dl['hierarchicallevels_LOW'].astype(float))
mH = measures(dh['Ceff_HIGH'].astype(float), dh['hierarchicallevels_HIGH'].astype(float))
df = pd.concat([pd.DataFrame({'subjectkey': ids_low, 'anhedonia_group': 'LOW', **mL}),
                pd.DataFrame({'subjectkey': ids_high, 'anhedonia_group': 'HIGH', **mH})], ignore_index=True)
df = df.merge(demos[['subjectkey', 'Age_clean', 'sex', 'is_patient']], on='subjectkey', how='left')
df['sex'] = df['sex'].replace('O', np.nan)
df = df.dropna(subset=['Age_clean', 'sex', 'is_patient']).copy()          # listwise deletion as in OLS
df['anhedonia_group'] = pd.Categorical(df['anhedonia_group'], categories=['LOW', 'HIGH'])
df['male'] = (df['sex'] == 'M').astype(float)
print(f'N = {len(df)} (LOW {sum(df.anhedonia_group=="LOW")}, HIGH {sum(df.anhedonia_group=="HIGH")})\n')

def fit_cellmeans(d, y, centre):
    d = d.copy()
    for v in ['Age_clean', 'male', 'is_patient']:
        d[v + '_c'] = d[v] - centre[v]
    return smf.ols(f'{y} ~ 0 + anhedonia_group + Age_clean_c + male_c + is_patient_c', data=d).fit()

def pooled_sd(v, g):
    lo, hi = v[g == 'LOW'], v[g == 'HIGH']
    return np.sqrt(((len(lo) - 1) * lo.var() + (len(hi) - 1) * hi.var()) / (len(lo) + len(hi) - 2))

centre = df[['Age_clean', 'male', 'is_patient']].mean()
for y in ['ctx_sub_gap', 'pfc_str_gap', 'ctx_sub_flow', 'pfc_str_flow']:
    chk = smf.ols(f'{y} ~ anhedonia_group + Age_clean + C(sex) + is_patient', data=df).fit()
    b = chk.params['anhedonia_group[T.HIGH]']
    print(f'=== {y}   group beta = {b:+.5f}, p = {chk.pvalues["anhedonia_group[T.HIGH]"]:.2g}, '
          f'd = {b / pooled_sd(df[y], df.anhedonia_group):+.2f}')
    r = fit_cellmeans(df, y, centre)
    for g in ['LOW', 'HIGH']:
        k = f'anhedonia_group[{g}]'; lo, hi = r.conf_int().loc[k]
        print(f'   {g:4s} adj. mean {r.params[k]:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]  p(≠0) = {r.pvalues[k]:.3g}')
    print()
