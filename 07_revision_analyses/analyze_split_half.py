"""Split-half reliability of the cortico-subcortical flow index, trophic gap and the two flow arms,
and replication of the group effect in each half. Reads the fits written by run_split_half.py."""
import glob, numpy as np, pandas as pd, scipy.io as sio, statsmodels.formula.api as smf
from scipy import stats
from hopf_lna import *
from config import SPLIT_HALF, DEMOS

full = {}
for g in ['low', 'high']:
    F = sio.loadmat(f'{GEC_OUTPUT}/{g}_anhedonia/results_Ceff_{g}_anhedonia.mat')[f'Ceff_{g.upper()}']
    ids = sio.loadmat(f'{CLUSTERS}/{g}_anhedonia.mat')[f'{g}_anhedonia']
    for k in range(len(ids)): full[ids[k, 1].flat[0]] = (g, F[k])
rows = []
for p in sorted(glob.glob(f'{SPLIT_HALF}/*.npz')):
    z = np.load(p, allow_pickle=True); sid = str(z['sid']); g = str(z['group'])
    r = {'id': sid, 'grp': int(g == 'high')}
    for tag, C in [('h1', z['half1']), ('h2', z['half2']), ('full', full[sid][1])]:
        for k, v in measures(C).items(): r[f'{k}_{tag}'] = v
    r['fitFC_h1'], r['fitFC_h2'] = z['half1_fit'][0], z['half2_fit'][0]
    rows.append(r)
df = pd.DataFrame(rows)
dm = pd.read_excel(DEMOS, header=1).set_index('subjectkey')
dm['age'] = dm['Age'].replace(999, np.nan); dm['pat'] = (dm['Primary_Dx_Payman'] != 999).astype(int)
df = df.join(dm[['age', 'sex', 'pat']], on='id')
print(f'n subjects = {len(df)}  (low {int((df.grp==0).sum())}, high {int((df.grp==1).sum())})\n')

def sb(r): return 2 * r / (1 + r)
print('SPLIT-HALF RELIABILITY (half1 vs half2)   pooled r | Spearman rho | Spearman-Brown | within-low r | within-high r')
for m in ['Q', 'gap', 'td', 'bu']:
    a, b = df[f'{m}_h1'], df[f'{m}_h2']
    r = np.corrcoef(a, b)[0, 1]; rho = stats.spearmanr(a, b)[0]
    rl = np.corrcoef(a[df.grp == 0], b[df.grp == 0])[0, 1]; rh = np.corrcoef(a[df.grp == 1], b[df.grp == 1])[0, 1]
    print(f'  {m:4s} {r:6.3f} | {rho:6.3f} | {sb(r):6.3f} | {rl:6.3f} | {rh:6.3f}')
print('\nHALF vs FULL-DATA value (r):', {m: (round(np.corrcoef(df[f"{m}_h1"], df[f"{m}_full"])[0,1],3), round(np.corrcoef(df[f"{m}_h2"], df[f"{m}_full"])[0,1],3)) for m in ['Q', 'gap']})
print('mean model fit (FC corr): half1 %.3f half2 %.3f' % (df.fitFC_h1.mean(), df.fitFC_h2.mean()))

print('\nGROUP EFFECT (high - low), OLS ~ group + age + sex + patient   [beta, p]')
d = df.dropna(subset=['age']); d = d[d.sex.isin(['M', 'F'])]
for m in ['Q', 'gap', 'td', 'bu']:
    s = []
    for tag in ['h1', 'h2', 'full']:
        f = smf.ols(f'{m}_{tag} ~ grp + age + sex + pat', d).fit()
        s.append(f'{tag}: {f.params["grp"]:+.4f} (p={f.pvalues["grp"]:.3g})')
    print(f'  {m:4s}', ' | '.join(s))
df.to_csv(f'{SPLIT_HALF}/split_half_measures.csv', index=False)   # subject-level, stays outside the repo
