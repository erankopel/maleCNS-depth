"""H3 (certified depth), PREREG Sec. 5, with DEVIATIONS.md C3 and C6.

Chain P on the weight >= 5 giant SCC. Witness rows: 90 random (numpy default_rng(7)) + 10 of
largest stationary mass. Lower bound on tau(P^r): max over witness pairs of
d_r(a, b) = (1/2) || (e_a - e_b) P^r ||_1, minus 1e-9 r. Upper bound (descriptive): Markov's
1 - sum_k min_i (P^r)_{ik} over the 500 columns of largest stationary mass, plus 1e-9 r.
Criteria (frozen): lower bound at r = 32 exceeds 0.20, and the first checkpoint
(1, 2, 4, 8, 16, 32, 48, 64, 96, 128) at which the lower bound drops below 0.05 lies in [48, 96].
Pass iff both. The verdict is saved as soon as it is determined; the loop then continues to r = 128
for the descriptive rows.
"""
import sys
import time
import numpy as np
from common import (log, load_edges, giant_scc, row_stochastic, stationary, pairwise_tv,
                    witness_sources, save_result, expect, SLACK, U, CHECKPOINTS)

pre, post, w, body_ids = load_edges()
keep, A = giant_scc(pre, post, w, len(body_ids), thresh=5)
expect(len(body_ids), 'G5', len(keep))
P = row_stochastic(A)
del A
n = P.shape[0]
D = int(max(np.diff(P.indptr).max(), np.bincount(P.indices, minlength=n).max()))
log(f'weight >= 5 giant SCC: n = {n:,}, nnz = {P.nnz:,}, max in- or out-degree D = {D:,}; '
    f'D u = {D * U:.2e} per hop (slack {SLACK:.0e} per hop)')
assert D * U * 10 < SLACK
pi, iters, delta = stationary(P, tol=1e-13, maxit=10000)
log(f'stationary distribution: {iters} iterations, final L1 change {delta:.1e}'
    + ('' if delta < 1e-13 else '  (tolerance not met; used only to choose witness rows, C3)'))

srcs = witness_sources(n, pi, 90, 10, seed=7)
hubs = np.argsort(-pi)[:500]
K, H = len(srcs), len(hubs)
X = np.zeros((K, n)); X[np.arange(K), srcs] = 1.0
Y = np.zeros((n, H)); Y[hubs, np.arange(H)] = 1.0
PT = P.T.tocsr()
rows, t0, saved = [], time.time(), False
common_info = {'hypothesis': 'H3', 'n': n, 'nnz': int(P.nnz), 'max_degree': D,
               'stationary_iterations': iters, 'stationary_delta': delta,
               'witness_rows': [int(x) for x in srcs], 'witness_distinct': int(len(np.unique(srcs)))}
for r in range(1, CHECKPOINTS[-1] + 1):
    X = (PT @ X.T).T
    Y = P @ Y
    if r not in CHECKPOINTS:
        continue
    dm = pairwise_tv(X)
    a, b = np.unravel_index(dm.argmax(), dm.shape)
    lo = float(dm.max() - SLACK * r)
    floor = float(Y.min(axis=0).sum())
    up = float(min(1.0, 1.0 - floor + SLACK * r))
    rows.append({'r': r, 'lower': lo, 'upper': up, 'hub_floor': floor,
                 'pair_body_ids': [int(body_ids[keep[srcs[a]]]), int(body_ids[keep[srcs[b]]])]})
    log(f'r = {r:3d}   lower {lo:.4f}   upper {up:.4f}   (hub floor {floor:.2e})   [{time.time() - t0:.0f} s]')
    first = next((x['r'] for x in rows if x['lower'] < 0.05), None)
    if not saved and r >= 32 and (first is not None or r == CHECKPOINTS[-1]):
        lo32 = next(x['lower'] for x in rows if x['r'] == 32)
        c1 = lo32 > 0.20
        c2 = first is not None and 48 <= first <= 96
        verdict = 'PASS' if (c1 and c2) else 'FAIL'
        save_result('h3', dict(common_info, rows=list(rows), first_checkpoint_below_0_05=first,
                               criteria={'lower32_gt_0.20': c1, 'first_below_0.05_in_48_96': c2},
                               verdict=verdict))
        log(f'\nlower bound at r = 32: {lo32:.4f} > 0.20: {c1}')
        log(f'first checkpoint with lower bound < 0.05: {first} (in [48, 96]): {c2}')
        log(f'H3 {verdict}\n')
        saved = True
save_result('h3_full', dict(common_info, scored=False, rows=rows))
sys.exit(0)
