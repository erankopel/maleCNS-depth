"""Structure of the male-CNS neuron graph, and a first look at mixing depth."""
import os
DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')), '')
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components, breadth_first_order

d = np.load(DATA + 'neuron_edges.npz')
pre, post, w = d['pre'], d['post'], d['w']
n = len(d['body_ids'])
print(f'neurons {n:,}   edges {len(pre):,}   synapses {int(w.sum()):,}\n')


def build(thresh):
    m = w >= thresh
    A = sp.csr_matrix((w[m].astype(np.float32), (pre[m], post[m])), shape=(n, n))
    return A


for thresh in (1, 5):
    A = build(thresh)
    nnz = A.nnz
    outdeg = np.diff(A.indptr)
    indeg = np.bincount(A.indices, minlength=n)
    ncomp, labels = connected_components(A, directed=True, connection='strong')
    sizes = np.bincount(labels)
    big = sizes.max()
    print(f'--- weight >= {thresh} ---')
    print(f'edges {nnz:,}  synapses {int(A.data.sum()):,}')
    print(f'out-degree  median {np.median(outdeg):.0f}  mean {outdeg.mean():.1f}  max {outdeg.max():,}')
    print(f'in-degree   median {np.median(indeg):.0f}  mean {indeg.mean():.1f}  max {indeg.max():,}')
    print(f'isolated (no in and no out): {int(((outdeg == 0) & (indeg == 0)).sum()):,}')
    print(f'strongly connected components: {ncomp:,}   largest {big:,} '
          f'({100*big/n:.1f}% of neurons)\n')
    if thresh == 5:
        A5, lab5, big5 = A, labels, sizes.argmax()

# hop reachability from random sources in the w>=5 graph
rng = np.random.default_rng(0)
core = np.flatnonzero(lab5 == big5)
srcs = rng.choice(core, 25, replace=False)
A1 = build(1)
Ab = (A1 > 0).astype(np.int8).tocsr()
reach = []
for s in srcs:
    order, _ = breadth_first_order(Ab, int(s), directed=True, return_predecessors=True)
    reach.append(len(order))
print(f'--- forward reachability (w>=1, 25 random core sources) ---')
print(f'neurons reachable: median {np.median(reach):,.0f} '
      f'({100*np.median(reach)/n:.1f}%)  min {min(reach):,}  max {max(reach):,}\n')

# mixing depth of the row-stochastic synapse-flow map
A = build(1)
rs = np.asarray(A.sum(axis=1)).ravel()
dangling = rs == 0
rs[dangling] = 1.0
P = sp.diags((1.0 / rs).astype(np.float32)) @ A
P = P.tocsr()

# stationary distribution by power iteration (ignoring dangling mass)
v = np.full(n, 1.0 / n, dtype=np.float64)
for _ in range(200):
    v_new = P.T @ v
    s = v_new.sum()
    if s <= 0:
        break
    v_new /= s
    if np.abs(v_new - v).sum() < 1e-12:
        v = v_new
        break
    v = v_new
pi = v
print('--- mixing of the row-stochastic synapse-flow map (w>=1) ---')
print(f'stationary mass on top 100 neurons: {np.sort(pi)[-100:].sum():.4f}')

rows = rng.choice(core, 5, replace=False)
print('hops  mean TV distance to stationary  mean support size')
dists = np.zeros((len(rows), n))
for i, r in enumerate(rows):
    dists[i, r] = 1.0
for r in range(1, 21):
    dists = dists @ P
    s = dists.sum(axis=1, keepdims=True)
    s[s == 0] = 1
    dists = dists / s
    tv = 0.5 * np.abs(dists - pi[None, :]).sum(axis=1)
    supp = (dists > 1e-12).sum(axis=1)
    if r <= 12 or r % 4 == 0:
        print(f'{r:4d}  {tv.mean():.4f}  {supp.mean():12,.0f}')
