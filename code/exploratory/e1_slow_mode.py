"""E1, exploratory (after the confirmatory run; not pre-registered): anatomy of the slow mode of the
synapse-flow chain on the weight >= 5 giant SCC, where H1 found |lambda_2| = 0.9965.

1. Right and left eigenvectors of lambda_2 (ARPACK, seeded as in C7).
2. Sweep cut along the right eigenvector: the prefix set of smallest conductance
   phi(S) = F(S -> S^c) / min(pi(S), pi(S^c)), F = stationary probability flow per hop.
3. Composition of the smaller side: superclass, soma neuromere, side, primary neuropil, cell types.
4. The same cut in the weight >= 1 chain: escape rates and the synapses that cross it.
The same analysis is repeated for the weight >= 1 chain (lambda_2 = 0.945) for comparison.
Run 2 (26 Sep): corrected synapse counts across the cut in step 4 (see the comment there).
ASCII, no em dashes.
"""
import json
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'confirmatory'))
from common import (DATA, log, load_edges, giant_scc, row_stochastic, stationary, top_eigs, annotations,
                    superclass_of, side_of, _json_default)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exploratory', 'results', '')
pre, post, w, body_ids = load_edges()
N = len(body_ids)
ann_all = annotations(body_ids, cols=('superclass', 'class', 'type', 'somaSide', 'somaNeuromere', 'rootSide',
                                      'statusLabel'))
side_all = side_of(ann_all['superclass'])
npz = np.load(DATA + 'neuron_neuropil.npz')
vocab = np.array(npz['vocab'])
prim_all = np.where(npz['primary'] >= 0, vocab[np.clip(npz['primary'], 0, None)], 'none')


def counts(values, mask, top=12):
    u, c = np.unique(values[mask].astype(str), return_counts=True)
    o = np.argsort(-c)[:top]
    return {str(u[i]): int(c[i]) for i in o}


def analyse(thresh):
    keep, A = giant_scc(pre, post, w, N, thresh=thresh)
    P = row_stochastic(A)
    n = P.shape[0]
    pi, it, dl = stationary(P, tol=1e-13, maxit=20000)
    vals, V, info = top_eigs(P, k=6, ncv=48, tol=1e-8, maxiter=6000, vectors=True)
    valsL, U, infoL = top_eigs(P.T.tocsr(), k=6, ncv=48, tol=1e-8, maxiter=6000, vectors=True)
    lam2 = vals[1]
    v2, u2 = V[:, 1], U[:, 1]
    log(f'\n=== weight >= {thresh}: n = {n:,}; lambda_1 = {vals[0]:.6f}, lambda_2 = {lam2:.6f}, '
        f'lambda_3 = {vals[2]:.6f}; left lambda_2 = {valsL[1]:.6f}')
    log(f'imaginary parts: v2 {np.abs(v2.imag).max():.1e}, u2 {np.abs(u2.imag).max():.1e}')
    v2, u2 = v2.real / np.abs(v2.real).max(), u2.real / np.abs(u2.real).sum()

    # sweep cut along v2: prefix sets in the order of increasing v2
    o = np.argsort(v2, kind='stable')
    rank = np.empty(n, np.int64); rank[o] = np.arange(n)
    rows = np.repeat(np.arange(n), np.diff(P.indptr)); cols = P.indices
    f = pi[rows] * P.data
    ri, rj = rank[rows], rank[cols]
    fwd = ri < rj                                  # edge from an earlier to a later node
    diff = (np.bincount(ri[fwd] + 1, weights=f[fwd], minlength=n + 1)
            - np.bincount(rj[fwd] + 1, weights=f[fwd], minlength=n + 1))
    F = np.cumsum(diff)[1:n]                       # F[k-1] = flow out of the first k nodes, k = 1..n-1
    piS = np.cumsum(pi[o])[:n - 1]
    phi = F / np.minimum(piS, 1 - piS)
    k = int(np.argmin(phi)) + 1
    S = np.zeros(n, bool); S[o[:k]] = True
    small = S if pi[S].sum() <= 0.5 else ~S
    p_out = F[k - 1] / pi[small].sum()
    q_out = F[k - 1] / pi[~small].sum()
    log(f'best sweep cut: |S| = {small.sum():,} neurons, stationary mass {pi[small].sum():.4f}, '
        f'conductance {phi[k - 1]:.2e}; escape per hop {p_out:.2e} (small side), {q_out:.2e} (rest); '
        f'two-block prediction 1 - p - q = {1 - p_out - q_out:.6f} vs lambda_2 {lam2.real:.6f}')
    gidx = keep[small]
    comp = {'superclass': counts(ann_all['superclass'], np.isin(np.arange(N), gidx)),
            'side_of_cns': counts(side_all, np.isin(np.arange(N), gidx)),
            'soma_neuromere': counts(ann_all['somaNeuromere'], np.isin(np.arange(N), gidx)),
            'soma_side': counts(ann_all['somaSide'], np.isin(np.arange(N), gidx)),
            'primary_neuropil': counts(prim_all, np.isin(np.arange(N), gidx)),
            'class': counts(ann_all['class'], np.isin(np.arange(N), gidx)),
            'type': counts(ann_all['type'], np.isin(np.arange(N), gidx), top=20)}
    for key, val in comp.items():
        log(f'  {key:<17s} {val}')
    # left eigenvector: where the slowly decaying deviation of a distribution sits
    pos, neg = u2 > 0, u2 < 0
    lsides = {sd: {'pos_share': float(u2[pos & (side_all[keep] == sd)].sum() / u2[pos].sum()),
                   'neg_share': float(-u2[neg & (side_all[keep] == sd)].sum() / -u2[neg].sum())}
              for sd in ('brain', 'cord', 'neck')}
    log(f'  left eigenvector shares by side: {lsides}')
    top = np.argsort(-np.abs(u2))[:15]
    carriers = [f'{ann_all["type"][keep[i]]}/{ann_all["superclass"][keep[i]]}/{ann_all["somaNeuromere"][keep[i]]}'
                f'/{ann_all["somaSide"][keep[i]]} {u2[i]:+.5f}' for i in top]
    log('  top left-eigenvector carriers: ' + '; '.join(carriers[:8]))
    return {'thresh': thresh, 'n': n, 'lambda': [complex(x) for x in vals], 'lambda_left': [complex(x) for x in valsL],
            'stationary_iterations': it, 'cut_size': int(small.sum()), 'cut_mass': float(pi[small].sum()),
            'conductance': float(phi[k - 1]), 'escape_small': float(p_out), 'escape_rest': float(q_out),
            'two_block_prediction': float(1 - p_out - q_out), 'composition': comp, 'left_sides': lsides,
            'left_top_carriers': carriers, 'cut_body_ids': body_ids[gidx].tolist(), 'arpack': [info, infoL]}


res5 = analyse(5)
res1 = analyse(1)

# the weight >= 5 cut, looked at in the weight >= 1 chain
cut = np.isin(body_ids, res5['cut_body_ids'])
keep1, A1 = giant_scc(pre, post, w, N, thresh=1)
P1 = row_stochastic(A1)
pi1, _, _ = stationary(P1, tol=1e-13, maxit=20000)
S1 = cut[keep1]
rows = np.repeat(np.arange(len(keep1)), np.diff(P1.indptr))
cross_out = S1[rows] & ~S1[P1.indices]
Fw1 = float((pi1[rows] * P1.data)[cross_out].sum())
# synapse counts from A1's own index arrays (run 1 paired A1.data with P1.indices, which store the
# entries of a row in a different order; its counts were wrong, see check_trap_synapses.py)
rowsA = np.repeat(np.arange(len(keep1)), np.diff(A1.indptr))
co, ci = S1[rowsA] & ~S1[A1.indices], ~S1[rowsA] & S1[A1.indices]
wcross = A1.data
syn_out, syn_in = int(wcross[co].sum()), int(wcross[ci].sum())
strong_out, strong_in = int(wcross[co & (wcross >= 5)].sum()), int(wcross[ci & (wcross >= 5)].sum())
log(f'\nthe weight >= 5 cut ({S1.sum():,} neurons in G) in the weight >= 1 chain: stationary mass '
    f'{pi1[S1].sum():.4f}, escape per hop {Fw1 / pi1[S1].sum():.2e}; synapses leaving the set {syn_out:,} '
    f'({strong_out:,} on edges >= 5), entering {syn_in:,} ({strong_in:,} on edges >= 5)')
out = {'w5': res5, 'w1': res1,
       'w5_cut_in_w1_chain': {'size_in_G': int(S1.sum()), 'stationary_mass': float(pi1[S1].sum()),
                              'escape_per_hop': Fw1 / float(pi1[S1].sum()), 'synapses_out': syn_out,
                              'synapses_out_on_edges_ge5': strong_out, 'synapses_in': syn_in,
                              'synapses_in_on_edges_ge5': strong_in}}
os.makedirs(OUT, exist_ok=True)
json.dump(out, open(OUT + 'e1_slow_mode.json', 'w'), indent=1, default=_json_default)
log('saved ' + OUT + 'e1_slow_mode.json')
