"""E1b, exploratory (after the confirmatory run; not pre-registered): the weight >= 5 chain without the
near-absorbing set found in E1.

E1 found that the slow mode of the weight >= 5 chain (lambda_2 = 0.9965) sits on 11 mesothoracic
efferent and associated neurons that hold 12.8 percent of the stationary mass, and the H3 witness pair
that sets the certified lower bound from r = 4 on has one member in that set. Here the chain is
recomputed (i) without those 11 neurons and (ii) without every efferent superclass (vnc_efferent,
cb_efferent, efferent_ascending, efferent_descending), each time on the giant SCC of what remains:
|lambda_2|, the next smallest-conductance sweep cut, and the certified witness bounds of H3
(100 witness rows, seed 7; 500 hub columns) at the H3 checkpoints. ASCII, no em dashes.
"""
import json
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'confirmatory'))
from common import (log, load_edges, giant_scc, row_stochastic, stationary, top_eigs, pairwise_tv,
                    witness_sources, annotations, SLACK, CHECKPOINTS, _json_default)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'exploratory', 'results', '')
pre, post, w, body_ids = load_edges()
N = len(body_ids)
ann = annotations(body_ids, cols=('superclass', 'type', 'somaNeuromere'))
trap = np.isin(body_ids, json.load(open(OUT + 'e1_slow_mode.json'))['w5']['cut_body_ids'])
efferent = np.isin(ann['superclass'], ['vnc_efferent', 'cb_efferent', 'efferent_ascending', 'efferent_descending'])
keep5, _ = giant_scc(pre, post, w, N, thresh=5)
in5 = np.zeros(N, bool); in5[keep5] = True


def run(label, removed):
    t0 = time.time()
    keep, A = giant_scc(pre, post, w, N, thresh=5, subset=in5 & ~removed)
    P = row_stochastic(A); n = P.shape[0]
    vals, V, _ = top_eigs(P, k=6, ncv=48, tol=1e-8, maxiter=6000, vectors=True)
    pi, it, _ = stationary(P, tol=1e-13, maxit=20000)
    v2 = V[:, 1].real
    o = np.argsort(v2, kind='stable'); rank = np.empty(n, np.int64); rank[o] = np.arange(n)
    rows = np.repeat(np.arange(n), np.diff(P.indptr)); f = pi[rows] * P.data
    ri, rj = rank[rows], rank[P.indices]; fwd = ri < rj
    F = np.cumsum(np.bincount(ri[fwd] + 1, weights=f[fwd], minlength=n + 1)
                  - np.bincount(rj[fwd] + 1, weights=f[fwd], minlength=n + 1))[1:n]
    piS = np.cumsum(pi[o])[:n - 1]; phi = F / np.minimum(piS, 1 - piS)
    c = int(np.argmin(phi)) + 1
    S = np.zeros(n, bool); S[o[:c]] = True
    small = S if pi[S].sum() <= 0.5 else ~S
    u, cnt = np.unique((ann['superclass'][keep][small]).astype(str), return_counts=True)
    ut, ct = np.unique((ann['type'][keep][small]).astype(str), return_counts=True)
    srcs = witness_sources(n, pi, 90, 10, seed=7)
    hubs = np.argsort(-pi)[:500]
    X = np.zeros((len(srcs), n)); X[np.arange(len(srcs)), srcs] = 1.0
    Y = np.zeros((n, len(hubs))); Y[hubs, np.arange(len(hubs))] = 1.0
    PT = P.T.tocsr(); bounds = {}
    for r in range(1, CHECKPOINTS[-1] + 1):
        X = (PT @ X.T).T
        Y = P @ Y
        if r in CHECKPOINTS:
            bounds[r] = (float(pairwise_tv(X).max() - SLACK * r), float(min(1.0, 1 - Y.min(axis=0).sum() + SLACK * r)))
    res = {'label': label, 'removed_in_scc5': int((in5 & removed).sum()), 'scc': int(n), 'lam2': float(abs(vals[1])),
           'lam3': float(abs(vals[2])), 'top1pct_mass': float(np.sort(pi)[::-1][:n // 100].sum()),
           'cut_size': int(small.sum()), 'cut_mass': float(pi[small].sum()), 'cut_conductance': float(phi[c - 1]),
           'cut_superclasses': {a: int(b) for a, b in sorted(zip(u, cnt), key=lambda x: -x[1])[:8]},
           'cut_types': {a: int(b) for a, b in sorted(zip(ut, ct), key=lambda x: -x[1])[:10]},
           'bounds': bounds}
    log(f'{label}: removed {res["removed_in_scc5"]}, SCC {n:,}, |lam2| {res["lam2"]:.6f}, |lam3| {res["lam3"]:.6f}; '
        f'next cut {res["cut_size"]:,} neurons (mass {res["cut_mass"]:.4f}, conductance {res["cut_conductance"]:.2e}) '
        f'{res["cut_superclasses"]}; bounds {bounds}  [{time.time() - t0:.0f} s]')
    return res


out = [run('without the 11-neuron set of E1', trap), run('without every efferent superclass', efferent)]
json.dump(out, open(OUT + 'e1b_without_trap.json', 'w'), indent=1, default=_json_default)
log('saved ' + OUT + 'e1b_without_trap.json')
