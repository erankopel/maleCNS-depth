"""E3, exploratory (after the confirmatory run; not pre-registered): null models for the slow modes of the
synapse-flow chain.

Each null keeps every edge's source and weight and re-draws its target by permuting targets within a
group of edges (a stub-matching null): out-degree, out-strength and the multiset of in-stubs are kept;
duplicate (source, target) pairs are merged by adding weights.
  N1  one group (configuration null);
  N2  groups = (superclass of source, superclass of target): superclass-to-superclass synapse counts kept;
  N3  groups = (cell type of source, cell type of target): type-to-type synapse counts kept (neurons
      without a type form their own group).
For each null, weight threshold (1 and 5) and seed: |lambda_2| of the giant SCC, its size, the stationary
mass of the top 1 percent of neurons, and the smallest-conductance sweep cut along the slow right
eigenvector. ASCII, no em dashes.
"""
import json
import os
import sys
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'confirmatory'))
from common import (log, load_edges, row_stochastic, stationary, top_eigs, annotations, NotDetermined,
                    _json_default)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exploratory', 'results', '')
pre, post, w, body_ids = load_edges()
N = len(body_ids)
ann = annotations(body_ids, cols=('superclass', 'type'))
sc_code = np.unique(ann['superclass'].astype(str), return_inverse=True)[1]
typ = np.where(ann['type'].astype(str) == '?', np.char.add('body', body_ids.astype(str)), ann['type'].astype(str))
ty_code = np.unique(typ, return_inverse=True)[1]


def slow_stats(pre_, post_, w_):
    A = sp.csr_matrix((w_.astype(np.float64), (pre_, post_)), shape=(N, N))
    _, lab = connected_components(A, directed=True, connection='strong')
    big = np.flatnonzero(lab == np.bincount(lab).argmax())
    P = row_stochastic(A[big][:, big])
    n = P.shape[0]
    try:
        vals, V, _ = top_eigs(P, k=6, ncv=48, tol=1e-8, maxiter=6000, vectors=True)
    except NotDetermined:
        return {'scc': int(n), 'lam2': None}
    pi, _, _ = stationary(P, tol=1e-12, maxit=20000)
    v2 = V[:, 1].real
    o = np.argsort(v2, kind='stable'); rank = np.empty(n, np.int64); rank[o] = np.arange(n)
    rows = np.repeat(np.arange(n), np.diff(P.indptr)); f = pi[rows] * P.data
    ri, rj = rank[rows], rank[P.indices]; fwd = ri < rj
    F = np.cumsum(np.bincount(ri[fwd] + 1, weights=f[fwd], minlength=n + 1)
                  - np.bincount(rj[fwd] + 1, weights=f[fwd], minlength=n + 1))[1:n]
    piS = np.cumsum(pi[o])[:n - 1]
    phi = F / np.minimum(piS, 1 - piS)
    k = int(np.argmin(phi)) + 1
    size = min(k, n - k)
    return {'scc': int(n), 'lam2': float(abs(vals[1])), 'lam3': float(abs(vals[2])),
            'top1pct_mass': float(np.sort(pi)[::-1][:n // 100].sum()),
            'cut_size': int(size), 'cut_conductance': float(phi[k - 1])}


def null_edges(m, groups, seed):
    rng = np.random.default_rng(seed)
    p, q, ww = pre[m], post[m].copy(), w[m]
    if groups is None:
        q = q[rng.permutation(len(q))]
    else:
        g = groups[m]
        o0 = np.argsort(g, kind='stable')                 # edges grouped, original order inside a group
        o1 = np.lexsort((rng.random(len(g)), g))          # same groups, random order inside a group
        q2 = q.copy()
        q2[o0] = q[o1]                                    # each edge takes the target of a random group mate
        q = q2
    key = p.astype(np.int64) * N + q
    u, inv = np.unique(key, return_inverse=True)
    return (u // N).astype(np.int64), (u % N).astype(np.int64), np.bincount(inv, weights=ww)


results = []
for thresh in (1, 5):
    m = w >= thresh
    t0 = time.time()
    obs = slow_stats(pre[m], post[m], w[m])
    obs.update(null='observed', thresh=thresh, seed=None)
    results.append(obs)
    log(f'w>={thresh} observed: {obs}  [{time.time() - t0:.0f} s]')
    groups = {'N1': None,
              'N2': sc_code[pre].astype(np.int64) * (sc_code.max() + 1) + sc_code[post],
              'N3': ty_code[pre].astype(np.int64) * (ty_code.max() + 1) + ty_code[post]}
    for name, g in groups.items():
        for s in range(3):
            t0 = time.time()
            r = slow_stats(*null_edges(m, g, 1000 + 10 * s + thresh))
            r.update(null=name, thresh=thresh, seed=1000 + 10 * s + thresh)
            results.append(r)
            log(f'w>={thresh} {name} seed {r["seed"]}: {r}  [{time.time() - t0:.0f} s]')
os.makedirs(OUT, exist_ok=True)
json.dump(results, open(OUT + 'e3_null_models.json', 'w'), indent=1, default=_json_default)
log('saved ' + OUT + 'e3_null_models.json')
