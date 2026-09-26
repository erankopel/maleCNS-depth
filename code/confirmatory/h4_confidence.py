"""H4 (robustness needs confidences), PREREG Sec. 5, with DEVIATIONS.md D4 and C4.

Chain P on G (weight >= 1 giant SCC), nominal weights w = minconf-0.5 counts.
Perturbation class at confidence level c: any w' with w_hi(c) <= w' <= w edgewise, where w_hi(c)
counts the synapse-partner rows of the edge with min(conf_pre, conf_post) >= c; P' = row
normalisation of w'. For every such P', row i moves by at most eps_i(c) = 1 - W_i(c) / W_i in
total variation (W = row sums). Certified robust lower bound (Sec. 4.4 formula with per-row eps,
t = 0 included; DEVIATIONS.md D4):
    tau(P'^r) >= max over witness pairs (a, b) of
                 d_r(a, b) - sum_{t=0}^{r-1} sum_i |x_t(i)| eps_i(c)  -  1e-9 r,
    x_t = (e_a - e_b) P^t,   d_r = (1/2) ||x_r||_1.
Witness rows: 90 random (default_rng(7)) + 10 of largest stationary mass. Checkpoints to 128.
Prediction (frozen): the robust lower bound on tau(P^32) at c = 0.7 exceeds 0.10. Pass iff so.
If not, the script reports the largest c on the grid at which it does (PREREG H4 fail clause).

Pass 1 checks that the flat-connectome weights equal the syn-partners counts edge by edge; it decodes
only the two body columns, and any mismatch stops the script before pass 2 (C4). Pass 2 stops if a
confidence value is missing, not finite or outside [0, 1]. The verdict is saved as soon as the r = 32
bound exists; the loop then continues to r = 128 for the descriptive rows (C6).
"""
import sys
import time
import numpy as np
from common import (log, load_edges, giant_scc, row_stochastic, stationary, pairwise_tv, witness_sources,
                    save_result, expect, open_synpartners, SLACK, U, CHECKPOINTS)

C_PRIMARY = 0.70
GRID = np.array([0.51, 0.52, 0.53, 0.54, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95])
DELTA_UNIFORM = 0.05                                  # the uniform class of Sec. 4.4, for comparison
EPS_UNIFORM = 2 * DELTA_UNIFORM / (1 - DELTA_UNIFORM)
GRID32 = GRID.astype(np.float32)                    # confidences are stored as float32
ip = int(np.flatnonzero(np.isclose(GRID, C_PRIMARY))[0])

pre, post, w, body_ids = load_edges()
N = len(body_ids)
keep, A = giant_scc(pre, post, w, N, thresh=1)
A.sort_indices()
n = len(keep)
inv = -np.ones(N, dtype=np.int64)
inv[keep] = np.arange(n)
log(f'G: {n:,} neurons, {A.nnz:,} edges, {int(A.data.sum()):,} synapses')
expect(N, 'G', n)
rowsA = np.repeat(np.arange(n, dtype=np.int64), np.diff(A.indptr))
ekeys = rowsA * n + A.indices.astype(np.int64)
assert (np.diff(ekeys) > 0).all()
del rowsA

rd1 = open_synpartners(['body_pre', 'body_post'])                           # pass 1: no confidence column
rd2 = open_synpartners(['body_pre', 'body_post', 'conf_pre', 'conf_post'])   # pass 2
NB = rd1.num_record_batches


def to_g(body):
    i = np.searchsorted(body_ids, body)
    np.clip(i, 0, N - 1, out=i)
    g = inv[i]
    g[body_ids[i] != body] = -1
    return g


# ---------------------------------------------------------------- pass 1: edge-level integrity
t0 = time.time()
ecount = np.zeros(A.nnz, dtype=np.int64)
buf, nbuf, rows_g, missing = [], 0, 0, 0
for bi in range(NB):
    b = rd1.get_batch(bi)
    assert b.schema.names == ['body_pre', 'body_post']
    gi = to_g(b.column('body_pre').to_numpy(zero_copy_only=False))
    gj = to_g(b.column('body_post').to_numpy(zero_copy_only=False))
    ok = (gi >= 0) & (gj >= 0)
    key = gi[ok] * n + gj[ok]
    pos = np.searchsorted(ekeys, key)
    np.clip(pos, 0, len(ekeys) - 1, out=pos)
    found = ekeys[pos] == key
    missing += int((~found).sum())
    rows_g += int(ok.sum())
    buf.append(pos[found]); nbuf += int(found.sum())
    if nbuf > 30_000_000:
        ecount += np.bincount(np.concatenate(buf), minlength=A.nnz); buf, nbuf = [], 0
if buf:
    ecount += np.bincount(np.concatenate(buf), minlength=A.nnz)
mismatch = int((ecount != A.data.astype(np.int64)).sum())
log(f'pass 1 ({time.time() - t0:.0f} s): syn-partners rows inside G {rows_g:,}; rows on no G edge {missing:,}; '
    f'edges whose weight differs from the syn-partners count {mismatch:,}')
if missing or mismatch:
    log('INTEGRITY FAILURE: stopping before any confidence value is read (DEVIATIONS.md C4).')
    save_result('h4_integrity_failure', {'rows_in_G': rows_g, 'rows_on_no_edge': missing,
                                          'edges_mismatched': mismatch})
    sys.exit(2)
del ecount, ekeys

# ---------------------------------------------------------------- pass 2: per-row confidence counts
t0 = time.time()
L = len(GRID)
hist = np.zeros(n * (L + 1), dtype=np.int64)       # bin b = number of grid levels <= min(conf)
below_half, bad_conf, buf, nbuf = 0, 0, [], 0
for bi in range(NB):
    b = rd2.get_batch(bi)
    gi = to_g(b.column('body_pre').to_numpy(zero_copy_only=False))
    gj = to_g(b.column('body_post').to_numpy(zero_copy_only=False))
    ok = (gi >= 0) & (gj >= 0)
    bad_conf += b.column('conf_pre').null_count + b.column('conf_post').null_count
    cm = np.minimum(b.column('conf_pre').to_numpy(zero_copy_only=False),
                    b.column('conf_post').to_numpy(zero_copy_only=False))[ok]   # float32, as stored
    bad_conf += int((~np.isfinite(cm)).sum() + (cm < 0).sum() + (cm > 1).sum())
    below_half += int((cm < np.float32(0.5)).sum())
    bins = np.searchsorted(GRID32, cm, side='right')      # comparisons in float32 (C4)
    buf.append(gi[ok] * (L + 1) + bins); nbuf += int(ok.sum())
    if nbuf > 30_000_000:
        hist += np.bincount(np.concatenate(buf), minlength=n * (L + 1)); buf, nbuf = [], 0
if buf:
    hist += np.bincount(np.concatenate(buf), minlength=n * (L + 1))
if bad_conf:
    log(f'CONFIDENCE FAILURE: {bad_conf:,} missing, non-finite or out-of-range values; stopping (C4).')
    save_result('h4_confidence_failure', {'bad_confidence_values': bad_conf})
    sys.exit(3)
hist = hist.reshape(n, L + 1)
W = np.asarray(A.sum(axis=1)).ravel()
assert (hist.sum(axis=1) == W).all()
# W_c[:, k] = rows of neuron i with min conf >= GRID[k]  (bins k+1 .. L)
Wc = np.cumsum(hist[:, ::-1], axis=1)[:, ::-1][:, 1:]
E = 1.0 - Wc / W[:, None]                           # eps_i(c), one column per grid level
E = np.hstack([E, np.full((n, 1), EPS_UNIFORM)])    # last column: uniform class, for comparison
log(f'pass 2 ({time.time() - t0:.0f} s): rows with min(conf) < 0.5: {below_half:,}')

P = row_stochastic(A)
del A
pi, iters, delta = stationary(P, tol=1e-13, maxit=10000)
log(f'stationary distribution: {iters} iterations, final L1 change {delta:.1e}')
eps_stats = {f'{c:.2f}': {'mean': float(E[:, k].mean()), 'median': float(np.median(E[:, k])),
                          'pi_weighted_mean': float(pi @ E[:, k]), 'max': float(E[:, k].max())}
             for k, c in enumerate(GRID)}
log(f'eps_i at c = 0.70: mean {eps_stats["0.70"]["mean"]:.4f}, median {eps_stats["0.70"]["median"]:.4f}, '
    f'stationary-weighted mean {eps_stats["0.70"]["pi_weighted_mean"]:.4f}')

# ---------------------------------------------------------------- witness propagation
D = int(max(np.diff(P.indptr).max(), np.bincount(P.indices, minlength=n).max()))
assert D * U * 10 < SLACK
srcs = witness_sources(n, pi, 90, 10, seed=7)
K = len(srcs)
X = np.zeros((K, n)); X[np.arange(K), srcs] = 1.0
S = np.zeros((K, K, E.shape[1]))                    # S[a, b, c] = sum_{t<r} sum_i |x_t(i)| eps_i(c), a < b
PT = P.T.tocsr()
out, t0, saved = [], time.time(), False
info = {'hypothesis': 'H4', 'grid': GRID.tolist(), 'c_primary': C_PRIMARY, 'rows_in_G': rows_g,
        'rows_min_conf_below_0.5': below_half, 'eps_stats': eps_stats, 'stationary_iterations': iters,
        'stationary_delta': delta, 'max_degree': D, 'witness_rows': [int(x) for x in srcs],
        'witness_distinct': int(len(np.unique(srcs)))}
for r in range(1, CHECKPOINTS[-1] + 1):
    for a in range(K - 1):                          # add the term of x_{r-1} (t = r - 1, t = 0 first)
        S[a, a + 1:, :] += np.abs(X[a][None, :] - X[a + 1:]) @ E
    X = (PT @ X.T).T
    if r not in CHECKPOINTS:
        continue
    dm = pairwise_tv(X)
    iu = np.triu_indices(K, 1)
    d = dm[iu]
    nominal = float(d.max() - SLACK * r)
    robust = d[:, None] - S[iu] - SLACK * r            # pairs x (grid + uniform)
    best = robust.max(axis=0)
    out.append({'r': r, 'nominal_lower': nominal,
                'robust_lower': {f'{c:.2f}': float(best[k]) for k, c in enumerate(GRID)},
                'robust_lower_uniform_0.05': float(best[-1])})
    log(f'r = {r:3d}  nominal {nominal:.4f}  robust c=0.70 {best[ip]:+.4f}  '
        f'c=0.55 {best[4]:+.4f}  c=0.51 {best[0]:+.4f}  uniform {best[-1]:+.4f}   [{time.time() - t0:.0f} s]')
    if r == 32 and not saved:
        lb = out[-1]['robust_lower'][f'{C_PRIMARY:.2f}']
        verdict = 'PASS' if lb > 0.10 else 'FAIL'
        ok_levels = [c for c in GRID if out[-1]['robust_lower'][f'{c:.2f}'] > 0.10]
        c_star = float(max(ok_levels)) if ok_levels else None
        save_result('h4', dict(info, rows=list(out), robust_lower_32_primary=lb,
                               largest_level_certified_0_10_at_32=c_star, verdict=verdict))
        log(f'\nrobust lower bound on tau(P^32) at c = 0.70: {lb:+.4f} (> 0.10: {lb > 0.10})')
        log(f'largest grid level with robust lower bound on tau(P^32) > 0.10: {c_star}')
        log(f'H4 {verdict}\n')
        saved = True
save_result('h4_full', dict(info, scored=False, rows=out))
sys.exit(0)
