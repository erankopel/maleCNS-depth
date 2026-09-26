"""Signed effective-connectivity map and spectra for male-cns v1.0.

Sign convention (Drosophila CNS, standard in FlyWire-style analyses):
  acetylcholine +1 ; gaba -1 ; glutamate -1 (GluCl-dominated) ; histamine -1 ;
  dopamine / serotonin / octopamine / tyramine / unknown -> 0 (modulatory, excluded).
"""
import os
DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')), '')
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import pyarrow.feather as feather

SIGN = {'acetylcholine': 1, 'gaba': -1, 'glutamate': -1, 'histamine': -1}

d = np.load(DATA + 'neuron_edges.npz')
pre, post, w, body_ids = d['pre'], d['post'], d['w'], d['body_ids']
n = len(body_ids)

nt = feather.read_table(DATA + 'body-nt.feather',
                        columns=['body', 'consensus_nt', 'predicted_nt_confidence']).to_pandas()
nt = nt.set_index('body')
nt = nt.reindex(body_ids)
cons = nt['consensus_nt'].fillna('unknown').to_numpy()
conf = nt['predicted_nt_confidence'].fillna(0).to_numpy()
sign = np.array([SIGN.get(c, 0) for c in cons], dtype=np.int8)

vals, cnts = np.unique(cons, return_counts=True)
order = np.argsort(-cnts)
print('consensus neurotransmitter among the 166,700 neurons:')
for i in order:
    print(f'  {vals[i]:<16s} {cnts[i]:>8,}   sign {SIGN.get(vals[i], 0):+d}')
print(f'excitatory {int((sign>0).sum()):,}  inhibitory {int((sign<0).sum()):,}  '
      f'excluded {int((sign==0).sum()):,}')
print(f'median NT confidence: {np.median(conf[conf>0]):.3f}\n')

mode = sys.argv[1] if len(sys.argv) > 1 else 'unsigned'

if mode == 'unsigned':
    # row-stochastic synapse-flow chain P (presynaptic normalisation), as in analyse.py
    A = sp.csr_matrix((w.astype(np.float64), (pre, post)), shape=(n, n))
    rs = np.asarray(A.sum(axis=1)).ravel(); rs[rs == 0] = 1
    P = (sp.diags(1.0 / rs) @ A).tocsr()
    print('top eigenvalues of P (largest magnitude), ARPACK:')
    ev = sla.eigs(P, k=8, which='LM', tol=1e-7, maxiter=4000, ncv=48,
                  return_eigenvectors=False)
    ev = ev[np.argsort(-np.abs(ev))]
    for e in ev:
        print(f'  |lam| = {abs(e):.6f}   lam = {e.real:+.6f} {e.imag:+.6f}i')
    lam2 = abs(ev[1])
    print(f'\nspectral gap 1-|lam2| = {1-lam2:.6f}')
    for eps in (0.1, 0.01, 0.001):
        print(f'  asymptotic depth for TV factor {eps}: r ~ {np.log(eps)/np.log(lam2):.0f} hops')

    # extend the empirical TV curve
    v = np.full(n, 1.0 / n)
    for _ in range(300):
        vn = P.T @ v; vn /= vn.sum()
        if np.abs(vn - v).sum() < 1e-13: v = vn; break
        v = vn
    pi = v
    rng = np.random.default_rng(0)
    rows = rng.choice(n, 20, replace=False)
    X = np.zeros((20, n)); X[np.arange(20), rows] = 1.0
    print('\nhops  mean TV to stationary (20 random sources)')
    for r in range(1, 201):
        X = X @ P
        X /= X.sum(axis=1, keepdims=True)
        if r in (1, 2, 4, 8, 16, 32, 64, 100, 150, 200):
            tv = 0.5 * np.abs(X - pi[None, :]).sum(axis=1)
            print(f'{r:4d}  {tv.mean():.4f}   (min {tv.min():.4f}, max {tv.max():.4f})')

elif mode == 'signed':
    # restrict to the giant SCC of the unsigned graph, then sign by presynaptic NT,
    # postsynaptic normalisation: each neuron's total (signed-eligible) input weight = 1
    keep = np.load(DATA + 'scc_keep.npy')
    inv = -np.ones(n, dtype=np.int64); inv[keep] = np.arange(len(keep))
    m = (sign[pre] != 0) & (inv[pre] >= 0) & (inv[post] >= 0)
    n = len(keep); sign = sign[keep]
    A = sp.csr_matrix((w[m].astype(np.float64), (inv[pre[m]], inv[post[m]])), shape=(n, n))
    cs = np.asarray(A.sum(axis=0)).ravel(); cs[cs == 0] = 1
    W = (A @ sp.diags(1.0 / cs)).tocsr()              # W[i,j] = w(i->j)/in(j)
    S = sp.diags(sign.astype(np.float64))
    M = (S @ W).T.tocsr()                              # x_{t+1} = M x_t ; M[j,i] = s_i w(i->j)/in(j)
    print(f'signed map: {M.nnz:,} edges;  signed-input fraction check: '
          f'{np.abs(np.asarray(abs(M).sum(axis=1)).ravel()).max():.3f} (should be <= 1)')
    print('top eigenvalues of the signed map M (largest magnitude), ARPACK:')
    ev = sla.eigs(M, k=10, which='LM', tol=1e-7, maxiter=4000, ncv=60,
                  return_eigenvectors=False)
    ev = ev[np.argsort(-np.abs(ev))]
    for e in ev:
        print(f'  |lam| = {abs(e):.6f}   lam = {e.real:+.6f} {e.imag:+.6f}i')
    print(f'\nspectral radius rho(M) = {abs(ev[0]):.6f}   '
          f'ratio |lam2/lam1| = {abs(ev[1])/abs(ev[0]):.6f}')
    # norm decay of signed propagation from random unit impulses
    rng = np.random.default_rng(1)
    X = np.zeros((n, 20)); X[rng.choice(n, 20, replace=False), np.arange(20)] = 1.0
    print('\nhops  mean ||M^r e||_1   mean ||M^r e||_2   mean net (sum) sign balance')
    for r in range(1, 65):
        X = M @ X
        if r in (1, 2, 4, 8, 16, 32, 64):
            l1 = np.abs(X).sum(axis=0).mean(); l2 = np.linalg.norm(X, axis=0).mean()
            bal = (X.sum(axis=0) / np.abs(X).sum(axis=0)).mean()
            print(f'{r:4d}  {l1:14.5f}  {l2:14.5f}  {bal:+.4f}')
