"""Spectral gap of the synapse-flow chain restricted to the giant strongly connected component."""
import os
DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')), '')
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components
d = np.load(DATA + 'neuron_edges.npz'); pre, post, w = d['pre'], d['post'], d['w']; n = len(d['body_ids'])
A = sp.csr_matrix((w.astype(np.float64), (pre, post)), shape=(n, n))
ncomp, lab = connected_components(A, directed=True, connection='strong')
giant = np.bincount(lab).argmax(); keep = np.flatnonzero(lab == giant)
A = A[keep][:, keep]; m = A.shape[0]
print(f'giant SCC: {m:,} neurons, {A.nnz:,} internal edges, {int(A.data.sum()):,} synapses')
rs = np.asarray(A.sum(axis=1)).ravel(); assert (rs > 0).all()
P = (sp.diags(1.0 / rs) @ A).tocsr()
# period check: gcd of cycle lengths via eigenvalues on unit circle; ARPACK for top-8 by magnitude
ev = sla.eigs(P, k=8, which='LM', tol=1e-8, maxiter=6000, ncv=48, return_eigenvectors=False)
ev = ev[np.argsort(-np.abs(ev))]
print('top eigenvalues of P on the giant SCC:')
for e in ev: print(f'  |lam| = {abs(e):.6f}   lam = {e.real:+.6f} {e.imag:+.6f}i')
lam2 = abs(ev[1]); print(f'\nspectral gap 1-|lam2| = {1-lam2:.6f}')
for eps in (0.1, 0.01, 0.001): print(f'  asymptotic depth for TV factor {eps}: r ~ {np.log(eps)/np.log(lam2):.0f} hops')
# stationary distribution and its concentration
v = np.full(m, 1.0/m)
for _ in range(2000):
    vn = P.T @ v; vn /= vn.sum()
    if np.abs(vn - v).sum() < 1e-13: v = vn; break
    v = vn
pi = np.sort(v)[::-1]; cum = np.cumsum(pi)
print(f'\nstationary mass: top 1% of neurons hold {cum[m//100-1]:.3f}, top 10% hold {cum[m//10-1]:.3f}; '
      f'entropy {-(v*np.log2(v)).sum():.2f} bits of {np.log2(m):.2f}')
np.save(DATA + 'scc_keep.npy', keep); np.save(DATA + 'scc_pi.npy', v)
