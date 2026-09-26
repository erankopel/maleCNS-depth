"""E2c, exploratory (after the confirmatory run; not pre-registered): what happens to the 11-neuron set of
E1 in the chains M(c) of E2 at weight >= 5?

For M(c) (every synapse below confidence c removed; rows left empty keep their weights) and for the
random reconstructions B(k) of E2 (the same draws: default_rng(101 + k), binomial per confidence bin),
threshold at weight >= 5, giant SCC: how many of the 11 neurons it contains, their stationary mass in
the SCC, and the synapses they send and receive on connections of >= 5 synapses. ASCII, no em dashes.
"""
import json
import os
import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'confirmatory'))
from common import DATA, log, load_edges, giant_scc, row_stochastic, stationary, _json_default

OUT = os.path.join(HERE, '..', '..', 'exploratory', 'results', '')
TRAP = [800719, 804848, 806145, 807560, 808551, 808845, 811415, 813264, 819906, 820114, 913939]
LOWER = np.round(np.arange(0.50, 0.96, 0.05), 2)

pre, post, w, body_ids = load_edges()
N = len(body_ids)
keep, A = giant_scc(pre, post, w, N, thresh=1)
del pre, post, w
A.sort_indices()
n = len(keep)
g_ids = body_ids[keep]
tr = np.flatnonzero(np.isin(g_ids, TRAP))
assert len(tr) == 11
hist = np.load(DATA + 'edge_conf_hist.npy', mmap_mode='r')
rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(A.indptr))
MID = np.array([0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975])


def members():
    for c in (0.50, 0.60, 0.70, 0.80, 0.90):
        k = int(np.flatnonzero(np.isclose(LOWER, c))[0])
        hi = np.zeros(A.nnz)
        for e0 in range(0, A.nnz, 2_000_000):               # contiguous blocks of the histogram
            hi[e0:e0 + 2_000_000] = np.asarray(hist[e0:e0 + 2_000_000, k:]).sum(axis=1)
        yield f'M({c:.2f})' if c > 0.5 else 'P', c, hi
    full = np.load(DATA + 'edge_conf_hist.npy')          # B(k) as in e2_class_members.py
    for k in range(5):
        rng = np.random.default_rng(101 + k)
        wp = np.zeros(A.nnz)
        for b in range(10):
            wp += rng.binomial(full[:, b], MID[b])
        yield f'B({k + 1})', None, wp
    del full


out = []
for label, c, hi in members():
    Hr = np.bincount(rows, weights=hi, minlength=n)
    wp = np.where((Hr == 0)[rows], A.data, hi)
    m = wp >= 5
    B = sp.csr_matrix((wp * m, A.indices.copy(), A.indptr.copy()), shape=(n, n)); B.eliminate_zeros()
    _, lab = connected_components(B, directed=True, connection='strong')
    big = lab == np.bincount(lab).argmax()
    inset = np.zeros(n, bool); inset[tr] = True
    out_syn = float(B[tr].sum()); out_syn_ext = float(B[tr][:, ~inset].sum())
    in_syn_ext = float(B[~inset][:, tr].sum())
    rec = {'member': label, 'c': c, 'scc_w5': int(big.sum()), 'trap_in_scc': int(big[tr].sum()),
           'trap_out_synapses_w5_to_rest': out_syn_ext, 'trap_in_synapses_w5_from_rest': in_syn_ext,
           'trap_out_synapses_w5_total': out_syn}
    if big[tr].any():
        idx = np.flatnonzero(big)
        Pb = row_stochastic(B[idx][:, idx])
        pib, _, _ = stationary(Pb, tol=1e-12, maxit=20000)
        pos = np.searchsorted(idx, tr[big[tr]])
        rec['trap_mass_in_scc'] = float(pib[pos].sum())
    out.append(rec)
    log(json.dumps(rec))
json.dump(out, open(OUT + 'e2c_trap_in_members.json', 'w'), indent=1, default=_json_default)
log('saved ' + OUT + 'e2c_trap_in_members.json')
