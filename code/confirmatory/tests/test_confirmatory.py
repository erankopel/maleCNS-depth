"""Integration tests of the confirmatory scripts on synthetic data (never the connectome).

A. Default synthetic CNS: run prep_neuropil and h1 to h5 and check every scored number, every verdict
   and every descriptive quantity against independent dense computations; check that reruns are
   bit-identical, that an existing result is never overwritten, and that the unscored-only mode must
   reproduce the saved scored result exactly (C6).
B. High-confidence variant: the H4 robust bounds are nontrivial and equal dense values up to r = 32.
C. Decoupled variant: removing the slowmode_check.py set separates brain and cord; H1 applies C1.
D. Corrupted variant: one synapse-partner row removed; H4 stops at the integrity pass.
Run: python test_confirmatory.py [workdir]   (about 10 minutes on 2 cores)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as feather
import pyarrow.ipc as ipc
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, CODE)
import make_synthetic                                   # noqa: E402
import common                                            # noqa: E402

work = sys.argv[1] if len(sys.argv) > 1 else tempfile.mkdtemp(prefix='malecns_test_')
FAIL = []


def check(cond, msg):
    print(('ok    ' if cond else 'FAIL  ') + msg)
    if not cond:
        FAIL.append(msg)


def run(script, data, res, **extra):
    env = dict(os.environ, MALECNS_DATA=data, MALECNS_RESULTS=res, OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1', **extra)
    return subprocess.run([sys.executable, os.path.join(CODE, script + '.py')], env=env,
                          capture_output=True, text=True)


def run_ok(script, data, res, **extra):
    rc = run(script, data, res, **extra)
    print(f'{script:<15s} exit {rc.returncode}')
    if rc.returncode != 0:
        print(rc.stdout[-3000:], rc.stderr[-3000:])
        raise SystemExit(1)
    return rc


def load(res, name):
    return json.load(open(os.path.join(res, name + '.json')))


class Ref:
    """Independent dense reference objects for one synthetic dataset."""

    def __init__(self, data):
        d = np.load(os.path.join(data, 'neuron_edges.npz'))
        self.body_ids = d['body_ids']
        self.N = N = len(self.body_ids)
        self.W = np.zeros((N, N)); np.add.at(self.W, (d['pre'], d['post']), d['w'].astype(float))
        ann = feather.read_table(os.path.join(data, 'body-annotations.feather')).to_pandas().set_index('bodyId')
        self.sc = ann.reindex(self.body_ids)['superclass'].to_numpy()
        ntab = feather.read_table(os.path.join(data, 'body-nt.feather')).to_pandas().set_index('body')
        self.nt = ntab.reindex(self.body_ids)['consensus_nt'].to_numpy()
        self.W5 = np.where(self.W >= 5, self.W, 0.0)
        self.G1 = self.scc(self.W, np.arange(N))
        self.G5 = self.scc(self.W5, np.arange(N))

    @staticmethod
    def labels(Wm, nodes):
        _, lab = connected_components(sp.csr_matrix(Wm[np.ix_(nodes, nodes)] > 0), directed=True,
                                      connection='strong')
        return lab

    def scc(self, Wm, nodes):
        lab = self.labels(Wm, nodes)
        return nodes[lab == np.bincount(lab).argmax()]

    @staticmethod
    def chain(Wm, nodes):
        Ws = Wm[np.ix_(nodes, nodes)]
        return Ws / Ws.sum(axis=1, keepdims=True)

    def lam2(self, Wm, nodes):
        return float(np.sort(np.abs(la.eigvals(self.chain(Wm, nodes))))[::-1][1])

    def confidences(self, data, G):
        inG = -np.ones(self.N, int); inG[G] = np.arange(len(G))
        syn = feather.read_table(os.path.join(data, 'syn-partners.feather'),
                                 columns=['body_pre', 'body_post', 'conf_pre', 'conf_post']).to_pandas()
        pos = pd.Series(np.arange(self.N), index=self.body_ids)
        g = []
        for col in ('body_pre', 'body_post'):
            ix = pos.reindex(syn[col]).to_numpy()
            gg = np.full(len(ix), -1)
            okk = ~np.isnan(ix)
            gg[okk] = inG[ix[okk].astype(int)]
            g.append(gg)
        ok = (g[0] >= 0) & (g[1] >= 0)
        cm = np.minimum(syn.conf_pre.to_numpy(), syn.conf_post.to_numpy())[ok]
        return g[0][ok], cm

    def robust(self, data, witness, levels, rmax):
        gi, cm = self.confidences(data, self.G1)
        Wrow = np.bincount(gi, minlength=len(self.G1)).astype(float)
        E = np.stack([1 - np.bincount(gi[cm >= np.float32(float(c))], minlength=len(self.G1)) / Wrow
                      for c in levels], axis=1)
        P1 = self.chain(self.W, self.G1)
        ia, ib = np.triu_indices(len(witness), 1)
        X = np.zeros((len(ia), len(self.G1)))
        X[np.arange(len(ia)), witness[ia]] += 1; X[np.arange(len(ia)), witness[ib]] -= 1
        S = np.zeros((len(ia), len(levels))); out = {}
        for r in range(1, rmax + 1):
            S += np.abs(X) @ E
            X = X @ P1
            out[r] = (0.5 * np.abs(X).sum(axis=1)[:, None] - S).max(axis=0) - 1e-9 * r
        return E, out


# ============================================================================== A. default dataset
DATA, RES = os.path.join(work, 'A', 'data'), os.path.join(work, 'A', 'results')
print('A. synthetic dataset:', make_synthetic.make(DATA, seed=0))
for s in ('prep_neuropil', 'h1_gap', 'h2_floor', 'h3_depth', 'h4_confidence', 'h5_sign'):
    run_ok(s, DATA, RES)
R = {h: load(RES, h) for h in ('h1', 'h2', 'h3', 'h4', 'h5')}
R1s = load(RES, 'h1_sensitivity')
ref = Ref(DATA)

# ---- H1
l2 = ref.lam2(ref.W5, ref.G5)
check(abs(l2 - R['h1']['rows']['full']['lam2']) < 1e-6, f'H1 |lam2| of the w>=5 SCC: dense {l2:.6f}')
for key, classes, src in (('neck_script', common.NECK_SCRIPT, R['h1']['rows']),
                          ('motor_script', common.MOTOR_SCRIPT, R['h1']['rows']),
                          ('neck_all', common.NECK_ALL, R1s['rows']), ('neck_text', common.NECK_TEXT, R1s['rows'])):
    nodes = ref.scc(ref.W5, ref.G5[~np.isin(ref.sc[ref.G5], classes)])
    l2k = ref.lam2(ref.W5, nodes)
    check(abs(l2k - src[key]['lam2']) < 1e-6 and len(nodes) == src[key]['scc_size'],
          f'H1 knock-out {key}: dense |lam2| {l2k:.6f}, SCC {len(nodes)}')
c = R['h1']['criteria']; rows = R['h1']['rows']
exp = [0.90 <= rows['full']['lam2'] <= 0.97, rows['neck_script']['gap'] < 0.002,
       abs(rows['motor_script']['lam2'] - rows['full']['lam2']) < 0.01]
check([c['i'], c['ii'], c['iii']] == exp and R['h1']['verdict'] == ('PASS' if all(exp) else 'FAIL'),
      f'H1 criteria and verdict follow from the numbers ({R["h1"]["verdict"]})')

# ---- H3
P5 = ref.chain(ref.W5, ref.G5)
src = np.array(R['h3']['witness_rows'])
pi5, _, _ = common.stationary(sp.csr_matrix(P5), tol=1e-13, maxit=10000)
hubs = np.argsort(-pi5)[:500]
X = np.eye(len(ref.G5))[src]; Y = np.eye(len(ref.G5))[:, hubs]; worst = 0.0
full3 = load(RES, 'h3_full')
for r in range(1, 129):
    X = X @ P5
    Y = P5 @ Y
    if r in common.CHECKPOINTS:
        lo = max(0.5 * np.abs(X[a] - X[b]).sum() for a in range(len(src)) for b in range(a + 1, len(src))) - 1e-9 * r
        up = min(1.0, 1 - Y.min(axis=0).sum() + 1e-9 * r)
        row = next(x for x in full3['rows'] if x['r'] == r)
        worst = max(worst, abs(lo - row['lower']), abs(up - row['upper']))
check(worst < 1e-9, f'H3 lower and upper bounds at all checkpoints equal the dense values (max diff {worst:.1e})')
lo32 = next(x['lower'] for x in full3['rows'] if x['r'] == 32)
first = next((x['r'] for x in full3['rows'] if x['lower'] < 0.05), None)
v3 = 'PASS' if (lo32 > 0.20 and first is not None and 48 <= first <= 96) else 'FAIL'
check(R['h3']['verdict'] == v3 and R['h3']['first_checkpoint_below_0_05'] == first,
      f'H3 verdict follows from the bounds ({v3}, first below 0.05 at {first})')

# ---- H4
levels = ('0.70', '0.51', '0.95')
E, rob = ref.robust(DATA, np.array(R['h4']['witness_rows']), levels, 16)
full4 = load(RES, 'h4_full')
for k, cl in enumerate(levels):
    check(abs(E[:, k].mean() - R['h4']['eps_stats'][cl]['mean']) < 1e-12, f'H4 mean eps_i at c = {cl}: {E[:, k].mean():.6f}')
worst = max(abs(rob[r][k] - next(x for x in full4['rows'] if x['r'] == r)['robust_lower'][cl])
            for r in (1, 2, 4, 8, 16) for k, cl in enumerate(levels))
check(worst < 1e-9, f'H4 robust lower bounds, r = 1..16, three levels, equal the dense values (max diff {worst:.1e})')
lb = R['h4']['robust_lower_32_primary']
check(R['h4']['verdict'] == ('PASS' if lb > 0.10 else 'FAIL'), f'H4 verdict follows from the r = 32 bound ({lb:.4f})')

# ---- prep_neuropil
syn = feather.read_table(os.path.join(DATA, 'syn-partners.feather'), columns=['body_pre', 'body_post', 'primary_post']).to_pandas()
npz = np.load(os.path.join(DATA, 'neuron_neuropil.npz'))
vocab = list(npz['vocab'])
code = syn.primary_post.cat.codes.to_numpy()
pos = pd.Series(np.arange(ref.N), index=ref.body_ids)
C = np.zeros((ref.N, len(vocab)), np.int64)
for col in ('body_pre', 'body_post'):
    ix = pos.reindex(syn[col]).to_numpy(); okk = ~np.isnan(ix)
    np.add.at(C, (ix[okk].astype(int), code[okk]), 1)
prim = np.where(C.sum(1) > 0, C.argmax(1), -1)
check(np.array_equal(C, npz['counts']) and np.array_equal(prim, npz['primary']),
      'prep_neuropil counts and primary neuropils equal the brute force')

# ---- H2 and H5 from dense eigenvectors
G = ref.G1; bid = ref.body_ids[G]
lopeb = np.isin(prim[G], [vocab.index('LOP(L)'), vocab.index('LOP(R)'), vocab.index('EB')])


def dense_modes(conv, floor):
    s = np.array([(common.SIGN_BASE if conv == 'base' else common.SIGN_POSPISIL).get(x, 0) for x in ref.nt[G]], float)
    As = ref.W[np.ix_(G, G)] * (s != 0)[:, None]
    inw = As.sum(0)
    den = np.maximum(inw, floor) if floor is not None else np.where(inw > 0, inw, 1.0)
    ev, V = la.eig((s[:, None] * As / den[None, :]).T)
    o = np.lexsort((-ev.imag, -np.abs(ev)))
    p = np.abs(V[:, o]) ** 2
    return ev[o], p / p.sum(0, keepdims=True), inw


ev, p, inw = dense_modes('base', None)
F = float(np.median(inw))
check(F == R['h2']['F'], f'H2 floor F = median signed-eligible input = {F}')
ev, p, _ = dense_modes('base', F)
pr = 1 / (p ** 2).sum(0)
got = R['h2']['primary']['modes']
check(np.allclose(np.abs(ev[:20]), [m['modulus'] for m in got], atol=1e-7), f'H2 top-20 moduli equal the dense spectrum')
check(np.allclose(pr[:20], [m['pr'] for m in got], rtol=1e-4), 'H2 participation ratios equal the dense eigenvectors')
check(np.allclose(p[lopeb, :20].sum(0), [m['lop_eb_share'] for m in got], atol=1e-6), 'H2 LOP+EB power shares equal the dense eigenvectors')
check(R['h2']['primary']['n_modes_pr_gt_20'] == int((pr[:20] > 20).sum()), 'H2 (b) count equals the dense count')
a2, b2 = abs(ev[0]) < 0.85, int((pr[:20] > 20).sum()) >= 10
check(R['h2']['verdict'] == ('PASS' if a2 and b2 else 'FAIL'), f'H2 verdict follows from the numbers ({R["h2"]["verdict"]})')
rho = {}
for conv in ('base', 'pospisil'):
    ev, p, _ = dense_modes(conv, None)
    union = sorted(set(int(bid[i]) for k in range(6) for i in common.carrier_set(p[:, k], 0.9)))
    rho[conv] = abs(ev[0])
    check(union == R['h5'][conv]['carrier_union_body_ids'] and abs(abs(ev[0]) - R['h5'][conv]['rho']) < 1e-7,
          f'H5 {conv}: spectral radius and carrier union equal the dense eigenvectors')
same = R['h5']['base']['carrier_union_body_ids'] == R['h5']['pospisil']['carrier_union_body_ids']
check(R['h5']['verdict'] == ('PASS' if abs(rho['base'] - rho['pospisil']) < 0.02 and same else 'FAIL'),
      f'H5 verdict follows from the numbers ({R["h5"]["verdict"]})')

# ---- C6 mechanics: determinism, no overwrite, unscored-only mode
RES2 = os.path.join(work, 'A', 'results_rerun')
run_ok('h2_floor', DATA, RES2)
same = all(common._strip(load(RES, n)) == common._strip(load(RES2, n)) for n in ('h2', 'h2_sensitivity'))
check(same, 'H2 rerun from scratch is bit-identical (seeded ARPACK, one thread)')
rc = run('h2_floor', DATA, RES)
check(rc.returncode != 0 and 'REFUSING TO OVERWRITE' in rc.stdout + rc.stderr, 'an existing scored result is never overwritten')
os.remove(os.path.join(RES, 'h2_sensitivity.json'))
rc = run('h2_floor', DATA, RES, MALECNS_UNSCORED_ONLY='1')
check(rc.returncode == 0 and 'reproduced exactly' in rc.stdout and os.path.exists(os.path.join(RES, 'h2_sensitivity.json')),
      'unscored-only mode reproduces the scored result exactly and completes the unscored part')

# ============================================================================== B. high confidence
DATA, RES = os.path.join(work, 'B', 'data'), os.path.join(work, 'B', 'results')
print('\nB. high-confidence dataset:', make_synthetic.make(DATA, seed=1, high_conf=True))
run_ok('h4_confidence', DATA, RES)
refB = Ref(DATA)
h4 = load(RES, 'h4'); full4 = load(RES, 'h4_full')
E, rob = refB.robust(DATA, np.array(h4['witness_rows']), ('0.70', '0.95'), 32)
worst = max(abs(rob[r][k] - next(x for x in full4['rows'] if x['r'] == r)['robust_lower'][cl])
            for r in (1, 2, 4, 8, 16, 32) for k, cl in enumerate(('0.70', '0.95')))
check(worst < 1e-9, f'H4 robust lower bounds up to r = 32 equal the dense values (max diff {worst:.1e})')
check(rob[1][0] > 0.5 and rob[2][0] > 0.2, f'H4 bounds are nontrivial in this regime (r=1: {rob[1][0]:.3f}, r=2: {rob[2][0]:.3f})')
check(h4['verdict'] == ('PASS' if h4['robust_lower_32_primary'] > 0.10 else 'FAIL'), f'H4 verdict consistent ({h4["verdict"]})')

# ============================================================================== C. decoupled
DATA, RES = os.path.join(work, 'C', 'data'), os.path.join(work, 'C', 'results')
print('\nC. decoupled dataset:', make_synthetic.make(DATA, seed=2, residual_neck=False))
run_ok('h1_gap', DATA, RES)
refC = Ref(DATA)
h1 = load(RES, 'h1')
rest = refC.G5[~np.isin(refC.sc[refC.G5], common.NECK_SCRIPT)]
lab = refC.labels(refC.W5, rest)
side = common.side_of(refC.sc[rest])
two = any(((side[lab == l] == 'brain').any() and (side[lab == l] == 'cord').any()) for l in np.unique(lab))
check(not two, 'dense check: after the neck removal no SCC holds both brain and cord neurons')
check(h1['rows']['neck_script'].get('decoupled') is True and h1['criteria']['ii'] is True,
      'H1 labels the row decoupled and scores (ii) as passed (C1)')

# ============================================================================== D. corrupted
DATA, RES = os.path.join(work, 'D', 'data'), os.path.join(work, 'D', 'results')
shutil.copytree(os.path.join(work, 'A', 'data'), DATA)
t = ipc.open_file(pa.memory_map(os.path.join(DATA, 'syn-partners.feather'))).read_all()
inG = ref.body_ids[ref.G1]
drop = int(np.flatnonzero(np.isin(t.column('body_pre').to_numpy(), inG) & np.isin(t.column('body_post').to_numpy(), inG))[0])
t = pa.concat_tables([t.slice(0, drop), t.slice(drop + 1)])
os.remove(os.path.join(DATA, 'syn-partners.feather'))
with ipc.new_file(os.path.join(DATA, 'syn-partners.feather'), t.schema) as wr:
    for b in t.to_batches(max_chunksize=20000):
        wr.write_batch(b)
rc = run('h4_confidence', DATA, RES)
check(rc.returncode == 2 and os.path.exists(os.path.join(RES, 'h4_integrity_failure.json'))
      and not os.path.exists(os.path.join(RES, 'h4.json')), 'H4 stops at the integrity pass when one row is missing')

print('\nALL INTEGRATION CHECKS PASSED' if not FAIL else f'\n{len(FAIL)} CHECK(S) FAILED')
raise SystemExit(1 if FAIL else 0)
