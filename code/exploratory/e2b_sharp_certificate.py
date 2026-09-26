"""E2b, exploratory (after the confirmatory run; not pre-registered): a sharper certificate for H4.

H4 charged row i the radius of Lemma 1, eps_i(c) = 1 - W_hi,i(c) / W_i, the fraction of its synapses
below confidence c. That radius is reached only when the low-confidence synapses of a row go to
partners that receive no high-confidence synapse from it. Here the robust lower bound of H4 is
recomputed, with the same witness rows (100, seed 7), checkpoints and slack 1e-9 r, for three radii:
  eps_i(c)    the radius of H4 (a check: this column must reproduce results/h4_full.json);
  eps*_i(c)   a certified upper bound on the exact largest change of row i over the class
              (radius.py: exact enumeration for rows with at most 12 partners, the fractional-knapsack
              relaxation otherwise);
  delta_i(c)  the actual change of row i in the vertex M(c) of E2. NOT a certificate: it shows how much
              of the slack of Lemma 2 remains when the radius is exact for one member of the class.
Levels c = 0.55, 0.60, ..., 0.80, the bin edges of data/edge_conf_hist.npy (from e0_edge_confidence.py);
with the argument 'low', levels c = 0.51 ... 0.54 from data/edge_conf_hist_low.npy
(e0b_low_confidence_bins.py), checkpoints to 48, output e2b_sharp_certificate_low.json.
ASCII, no em dashes.
"""
import json
import os
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'confirmatory'))
sys.path.insert(0, HERE)
from common import (DATA, log, load_edges, giant_scc, row_stochastic, stationary, pairwise_tv, witness_sources,
                    SLACK, U, _json_default)
from radius import certified_radius

OUT = os.path.join(os.environ.get('MALECNS_EXPL_OUT', os.path.join(HERE, '..', '..', 'exploratory', 'results')), '')
LOW = len(sys.argv) > 1 and sys.argv[1] == 'low'
LEVELS = [0.51, 0.52, 0.53, 0.54] if LOW else [0.55, 0.60, 0.65, 0.70, 0.75, 0.80]
LOWER = np.round(np.arange(0.50, 0.96, 0.05), 2)            # lower edges of the ten confidence bins
LOWER_LOW = np.array([0.50, 0.51, 0.52, 0.53, 0.54])         # lower edges of the five bins below 0.55
CHECK = [1, 2, 4, 8, 16, 32, 48] if LOW else [1, 2, 4, 8, 16, 32, 48, 64]
BLOCK = 20000                                                # rows per block for the radius computation
NAME = 'e2b_sharp_certificate_low.json' if LOW else 'e2b_sharp_certificate.json'

pre, post, w, body_ids = load_edges()
N = len(body_ids)
keep, A = giant_scc(pre, post, w, N, thresh=1)
del pre, post, w
A.sort_indices()
n = len(keep)
indptr = A.indptr.astype(np.int64)
hist = np.load(DATA + 'edge_conf_hist.npy', mmap_mode='r')
hist_low = np.load(DATA + 'edge_conf_hist_low.npy', mmap_mode='r') if LOW else None
ekeys = np.load(DATA + 'edge_conf_keys.npy', mmap_mode='r')
wA = A.data.astype(np.float64)
for b0 in range(0, n, BLOCK):                                # the histogram is in the CSR order of G
    b1 = min(n, b0 + BLOCK); e0, e1 = indptr[b0], indptr[b1]
    rr = np.repeat(np.arange(b0, b1, dtype=np.int64), np.diff(indptr[b0:b1 + 1]))
    assert (np.asarray(ekeys[e0:e1]) == rr * n + A.indices[e0:e1]).all()
    assert (np.asarray(hist[e0:e1]).sum(axis=1) == A.data[e0:e1]).all()
log(f'G: {n:,} neurons, {A.nnz:,} edges; confidence histogram matches edge by edge')

P = row_stochastic(A)
pi, iters, delta_pi = stationary(P, tol=1e-13, maxit=10000)
log(f'stationary distribution: {iters} iterations, final L1 change {delta_pi:.1e}')

# ---------------------------------------------------------------- radii, block by block
t0 = time.time()
cols, names, stats = [], [], {}
for c in LEVELS:
    k = int(np.flatnonzero(np.isclose(LOWER_LOW if LOW else LOWER, c))[0])
    eps = np.zeros(n); eps_star = np.zeros(n); dlt = np.zeros(n)
    for b0 in range(0, n, BLOCK):
        b1 = min(n, b0 + BLOCK); e0, e1 = indptr[b0], indptr[b1]
        rr = np.repeat(np.arange(b1 - b0, dtype=np.int64), np.diff(indptr[b0:b1 + 1]))
        if LOW:                                                           # synapses with conf >= c
            hb = (np.asarray(hist[e0:e1, 1:]).sum(axis=1) + np.asarray(hist_low[e0:e1, k:]).sum(axis=1)).astype(np.float64)
        else:
            hb = np.asarray(hist[e0:e1, k:]).sum(axis=1).astype(np.float64)
        wb = wA[e0:e1]
        es, ee = certified_radius(rr, wb, hb, b1 - b0)
        eps[b0:b1], eps_star[b0:b1] = ee, es
        # the vertex M(c): every synapse below c removed; a row left empty keeps its weights (C4)
        Wb = np.bincount(rr, weights=wb, minlength=b1 - b0)
        Hb = np.bincount(rr, weights=hb, minlength=b1 - b0)
        empty = Hb == 0
        wp = np.where(empty[rr], wb, hb)
        Wp = np.where(empty, Wb, Hb)
        dlt[b0:b1] = 0.5 * np.bincount(rr, weights=np.abs(wb / Wb[rr] - wp / Wp[rr]), minlength=b1 - b0)
    assert (eps_star <= eps + 1e-15).all() and (dlt <= eps + 1e-12).all()
    assert (dlt <= eps_star + 1e-9).all(), 'a member of the class moved a row by more than eps*'
    stats[f'{c:.2f}'] = {nm: {'mean': float(v.mean()), 'median': float(np.median(v)), 'pi_weighted_mean': float(pi @ v),
                              'max': float(v.max())} for nm, v in (('eps', eps), ('eps_star', eps_star), ('delta_M', dlt))}
    stats[f'{c:.2f}']['rows_eps_star_below_eps'] = int((eps_star < eps - 1e-12).sum())
    s = stats[f'{c:.2f}']
    log(f'c = {c:.2f}: mean eps {s["eps"]["mean"]:.4f}  eps* {s["eps_star"]["mean"]:.4f}  delta(M) '
        f'{s["delta_M"]["mean"]:.4f}   (pi-weighted {s["eps"]["pi_weighted_mean"]:.4f} / '
        f'{s["eps_star"]["pi_weighted_mean"]:.4f} / {s["delta_M"]["pi_weighted_mean"]:.4f})  [{time.time() - t0:.0f} s]')
    cols += [eps, eps_star, dlt]; names += [f'eps_{c:.2f}', f'eps_star_{c:.2f}', f'delta_M_{c:.2f}']
E = np.column_stack(cols)
del cols

# ---------------------------------------------------------------- witness propagation (as in H4)
D = int(max(np.diff(P.indptr).max(), np.bincount(P.indices, minlength=n).max()))
assert D * U * 10 < SLACK
srcs = witness_sources(n, pi, 90, 10, seed=7)
K = len(srcs)
X = np.zeros((K, n)); X[np.arange(K), srcs] = 1.0
S = np.zeros((K, K, E.shape[1]))
PT = P.T.tocsr()
iu = np.triu_indices(K, 1)
rows_out, t0 = [], time.time()
for r in range(1, CHECK[-1] + 1):
    for a in range(K - 1):
        S[a, a + 1:, :] += np.abs(X[a][None, :] - X[a + 1:]) @ E
    X = (PT @ X.T).T
    if r not in CHECK:
        continue
    d = pairwise_tv(X)[iu]
    best = (d[:, None] - S[iu] - SLACK * r).max(axis=0)
    rec = {'r': r, 'nominal_lower': float(d.max() - SLACK * r), 'robust_lower': dict(zip(names, map(float, best)))}
    rows_out.append(rec)
    ca, cb = (f'{LEVELS[-1]:.2f}', f'{LEVELS[0]:.2f}') if LOW else ('0.70', '0.55')
    log(f'r = {r:3d}  nominal {rec["nominal_lower"]:.4f}  c={ca}: eps {best[names.index("eps_" + ca)]:+.4f} '
        f'eps* {best[names.index("eps_star_" + ca)]:+.4f} delta {best[names.index("delta_M_" + ca)]:+.4f}   '
        f'c={cb}: eps {best[names.index("eps_" + cb)]:+.4f} eps* {best[names.index("eps_star_" + cb)]:+.4f}   '
        f'[{time.time() - t0:.0f} s]')
json.dump({'levels': LEVELS, 'checkpoints': CHECK, 'witness_rows': [int(x) for x in srcs], 'radius_stats': stats,
           'rows': rows_out}, open(OUT + NAME, 'w'), indent=1, default=_json_default)
log('saved ' + OUT + NAME)
