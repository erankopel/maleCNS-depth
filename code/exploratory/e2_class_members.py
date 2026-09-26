"""E2, exploratory (after the confirmatory run; not pre-registered): is the H4 failure a property of the
certificate or of the data?

H4's certificate bounds tau(P'^r) for every P' in the confidence class at once, by a first-order
argument, and fails from r = 8 at c = 0.70. Here we compute exactly, for specific chains P', how much
the depth numbers move:
  M(c)   all synapse-partner rows with min(conf_pre, conf_post) < c removed (a vertex of the class at
         level c; a row that would lose all its weight keeps its original weights, which the class
         allows), for c = 0.60, 0.70, 0.80, 0.90;
  B(k)   a random reconstruction: each synapse kept with probability equal to the midpoint of its
         confidence bin (not a member of the class; a natural noise model), five draws.
For each chain: the witness distances d_r(a, b) of the H4 witness set (100 rows, seed 7) at
r = 1 ... 64, their change against P, and |lambda_2| of the giant SCC at weight >= 1 and >= 5.
Needs data/edge_conf_hist.npy from e0_edge_confidence.py. ASCII, no em dashes.
"""
import json
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'confirmatory'))
from common import (DATA, log, load_edges, giant_scc, row_stochastic, stationary, pairwise_tv, witness_sources,
                    second_modulus, NotDetermined, _json_default)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exploratory', 'results', '')
CHECK = [1, 2, 4, 8, 16, 32, 48, 64]
MID = np.array([0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975])
LOWER = np.array([0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95])

pre, post, w, body_ids = load_edges()
N = len(body_ids)
keep, A = giant_scc(pre, post, w, N, thresh=1)
A.sort_indices()
n = len(keep)
hist = np.load(DATA + 'edge_conf_hist.npy')
ekeys = np.load(DATA + 'edge_conf_keys.npy')
rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(A.indptr))
assert (ekeys == rows * n + A.indices).all() and (hist.sum(axis=1) == A.data).all()
P0 = row_stochastic(A)
pi0, _, _ = stationary(P0, tol=1e-13, maxit=10000)
srcs = witness_sources(n, pi0, 90, 10, seed=7)            # the H4 witness set
iu = np.triu_indices(len(srcs), 1)


def chain_from(wp):
    """Row-normalised chain on G from edge weights wp; rows left without weight keep w (C4 class rule)."""
    Wp = np.bincount(rows, weights=wp, minlength=n)
    empty = Wp == 0
    wp = np.where(empty[rows], A.data, wp)
    Ap = sp.csr_matrix((wp.astype(np.float64), A.indices.copy(), A.indptr.copy()), shape=(n, n))
    return Ap, int(empty.sum())


def witness_d(P):
    X = np.zeros((len(srcs), n)); X[np.arange(len(srcs)), srcs] = 1.0
    PT = P.T.tocsr(); out = {}
    for r in range(1, CHECK[-1] + 1):
        X = (PT @ X.T).T
        if r in CHECK:
            out[r] = pairwise_tv(X)[iu]
    return out


def lam2_scc(Ap, thresh):
    m = Ap.data >= thresh
    B = sp.csr_matrix((Ap.data * m, Ap.indices.copy(), Ap.indptr.copy()), shape=(n, n)); B.eliminate_zeros()
    _, lab = connected_components(B, directed=True, connection='strong')
    big = np.flatnonzero(lab == np.bincount(lab).argmax())
    try:
        l2, _, _ = second_modulus(row_stochastic(B[big][:, big]), k=8, ncv=48, tol=1e-8, maxiter=6000)
    except NotDetermined:
        l2 = None
    return int(len(big)), l2


def run(label, wp):
    t0 = time.time()
    Ap, empty = chain_from(wp)
    P = row_stochastic(Ap)
    tvd = 0.5 * np.abs(P - P0).sum(axis=1).A1
    d = witness_d(P)
    s1, l1 = lam2_scc(Ap, 1)
    s5, l5 = lam2_scc(Ap, 5)
    res = {'label': label, 'synapses': int(wp.sum()), 'rows_kept_whole': empty,
           'row_tv_to_P_mean': float(tvd.mean()), 'row_tv_to_P_pi_weighted': float(pi0 @ tvd),
           'scc_w1': s1, 'lam2_w1': l1, 'scc_w5': s5, 'lam2_w5': l5,
           'max_d': {r: float(d[r].max()) for r in CHECK}}
    if 'P' in D0:
        res['max_abs_change_d'] = {r: float(np.abs(d[r] - D0['P'][r]).max()) for r in CHECK}
        res['change_of_max_d'] = {r: float(d[r].max() - D0['P'][r].max()) for r in CHECK}
    log(f'{label:<34s} synapses {wp.sum():>11,.0f}  row TV mean {tvd.mean():.3f}  lam2 w>=1 {l1:.6f} '
        f'(SCC {s1:,})  w>=5 {l5:.6f} (SCC {s5:,})  max d_32 {d[32].max():.4f}'
        + (f'  max |change d_32| {res["max_abs_change_d"][32]:.4f}' if 'P' in D0 else '') + f'  [{time.time() - t0:.0f} s]')
    return res, d


# Each member's result is appended to a JSON-lines file as soon as it exists, so an interrupted run
# resumes where it stopped (the members are deterministic: fixed data and fixed seeds).
os.makedirs(OUT, exist_ok=True)
JL = OUT + 'e2_class_members.jsonl'
DP = OUT + 'e2_P_witness_d.npz'
done = {}
if os.path.exists(JL):
    for line in open(JL):
        rec = json.loads(line); done[rec['label']] = rec
D0 = {}
if os.path.exists(DP):
    z = np.load(DP); D0['P'] = {int(k[1:]): z[k] for k in z.files}


def record(label, wp):
    if label in done:
        log(f'{label}: already computed, skipped')
        return
    res, d = run(label, wp)
    if label.startswith('P '):
        np.savez(DP, **{f'r{r}': d[r] for r in CHECK})
        D0['P'] = d
    with open(JL, 'a') as f:
        f.write(json.dumps(res, default=_json_default) + '\n')


record('P (all synapses, minconf 0.5)', A.data.astype(np.float64))
if 'P' not in D0:
    raise SystemExit('reference distances of P missing')
for c in (0.60, 0.70, 0.80, 0.90):
    record(f'M({c:.2f}): rows below {c:.2f} removed', hist[:, LOWER >= c - 1e-9].sum(axis=1).astype(np.float64))
for k in range(5):
    label = f'B({k + 1}): kept with prob. = confidence'
    if label in done:
        log(f'{label}: already computed, skipped'); continue
    rng = np.random.default_rng(101 + k)
    wp = np.zeros(len(A.data))
    for b in range(10):
        wp += rng.binomial(hist[:, b], MID[b])
    record(label, wp)
recs = [json.loads(l) for l in open(JL)]
json.dump({'witness_rows': srcs.tolist(), 'checkpoints': CHECK, 'members': recs},
          open(OUT + 'e2_class_members.json', 'w'), indent=1, default=_json_default)
log('saved ' + OUT + 'e2_class_members.json')
