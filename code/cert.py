"""Item 4: certified bounds on tau(P^r), the Dobrushin coefficient of the r-hop synapse-flow map.

tau(Q) = (1/2) max_{a,b} || Q_a - Q_b ||_1 is the exact classical analogue of the
entanglement-breaking question: tau(P^r) = 0 iff P^r has rank one (every start is forgotten).

Lower bound (certified):  tau(P^r) >= (1/2) || (e_a - e_b) P^r ||_1 for any witness pair (a,b).
Upper bound (certified):  tau(P^r) <= 1 - sum_k min_i (P^r)_{ik}  for any set of columns k  (Markov 1906).
Perturbation class:       every edge weight w -> w(1+eta), |eta| <= delta.  Then for the perturbed P'
    (P'^r)_{ik} >= ((1-delta)/(1+delta))^r (P^r)_{ik}                       (upper bound degrades by this factor)
    | d_r(P') - d_r(P) | <= eps * sum_{t<r} 2 d_t(P),  eps = 2 delta/(1-delta)   (lower bound degrades by this)
Rounding: each matvec entry is a sum of at most D = max degree nonnegative terms, so relative error
    per hop <= D * 2^-53 < 1.3e-12; we add 1e-9 * r as slack to every printed bound.
"""
import os
DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')), '')
import numpy as np, scipy.sparse as sp, time
d = np.load(DATA + 'neuron_edges.npz'); pre, post, w, body_ids = d['pre'], d['post'], d['w'], d['body_ids']
keep = np.load(DATA + 'scc_keep.npy'); pi = np.load(DATA + 'scc_pi.npy')
inv = -np.ones(len(body_ids), dtype=np.int64); inv[keep] = np.arange(len(keep)); n = len(keep)
m = (inv[pre] >= 0) & (inv[post] >= 0)
A = sp.csr_matrix((w[m].astype(np.float64), (inv[pre[m]], inv[post[m]])), shape=(n, n))
rs = np.asarray(A.sum(axis=1)).ravel(); P = (sp.diags(1.0 / rs) @ A).tocsr(); del A
D = int(np.diff(P.indptr).max()); slack = 1e-9
print(f'giant SCC n={n:,}  nnz={P.nnz:,}  max degree D={D:,}\n')

rng = np.random.default_rng(7)
KW, KH = 40, 60
srcs = np.concatenate([rng.choice(n, KW - 10, replace=False), np.argsort(-pi)[:10]])   # random + 10 hubs as sources
hubs = np.argsort(-pi)[:KH]                                                               # hub columns for the upper bound
checkpoints = [1, 2, 4, 8, 16, 32, 64]
R = checkpoints[-1]

# rows: X_t = E_src P^t  (K x n);  columns: Y_t = P^t E_hub  (n x K)
X = np.zeros((KW, n)); X[np.arange(KW), srcs] = 1.0
Y = np.zeros((n, KH)); Y[hubs, np.arange(KH)] = 1.0
PT = P.T.tocsr()
dsum = np.zeros((KW, KW))           # sum_{t<r} 2 d_t for the perturbation term (pairwise)
rows = []
t0 = time.time()
for r in range(1, R + 1):
    # pairwise distances before the step contribute to the perturbation sum for this step
    if r > 1:
        dsum += 2 * dmat
    X = (PT @ X.T).T                # X P
    Y = P @ Y                       # P Y
    # pairwise TV among witness rows
    dmat = np.zeros((KW, KW))
    for a in range(KW):
        dmat[a] = 0.5 * np.abs(X[a][None, :] - X).sum(axis=1)
    if r in checkpoints:
        lo_pair = np.unravel_index(dmat.argmax(), dmat.shape)
        lo = dmat.max() - slack * r
        up = 1.0 - Y.min(axis=0).sum() + slack * r
        rows.append((r, lo, up, dsum[lo_pair], srcs[lo_pair[0]], srcs[lo_pair[1]], Y.min(axis=0).sum()))
        print(f'r={r:3d}  lower {lo:.4f}  upper {up:.4f}   (hub-column mass floor {Y.min(axis=0).sum():.2e})  '
              f'[{time.time()-t0:.0f}s]', flush=True)

print('\n=== certified bounds on tau(P^r), giant SCC, minconf-0.5 weights ===')
print('  r   lower bound (witness pair)   upper bound (60 hub columns)   robust lower, delta=0.05   robust upper, delta=0.05')
for r, lo, up, ds, a, b, floor in rows:
    delta = 0.05; eps = 2 * delta / (1 - delta)
    lo_rob = max(0.0, lo - eps * ds)
    up_rob = min(1.0, 1.0 - ((1 - delta) / (1 + delta)) ** r * floor + slack * r)
    print(f'{r:4d}   {lo:22.4f}   {up:26.4f}   {lo_rob:22.4f}   {up_rob:22.4f}')
print('\nwitness pair at r=64 (giant-SCC indices -> body IDs):',
      body_ids[keep[rows[-1][4]]], body_ids[keep[rows[-1][5]]])
print('reading: tau(P^r) is the largest TV distance any two starting neurons can retain after r hops;')
print('         P^r is information-breaking to within eps once the upper bound falls below eps.')
