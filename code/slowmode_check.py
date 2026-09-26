import os
DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')), '')
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, pyarrow.feather as feather
d = np.load(DATA + 'neuron_edges.npz'); pre, post, w, body_ids = d['pre'], d['post'], d['w'], d['body_ids']
keep = np.load(DATA + 'scc_keep.npy'); pi = np.load(DATA + 'scc_pi.npy')
inv = -np.ones(len(body_ids), dtype=np.int64); inv[keep] = np.arange(len(keep)); n = len(keep)
ann = feather.read_table(DATA + 'body-annotations.feather', columns=['bodyId','superclass']).to_pandas().astype(object)
sc = ann.set_index('bodyId').reindex(body_ids[keep])['superclass'].fillna('?').to_numpy()
m = (inv[pre] >= 0) & (inv[post] >= 0); p2, q2, w2 = inv[pre[m]], inv[post[m]], w[m].astype(np.float64)
A = sp.csr_matrix((w2, (p2, q2)), shape=(n, n))
outw = np.asarray(A.sum(axis=1)).ravel(); inw = np.asarray(A.sum(axis=0)).ravel()
for g in ['vnc_motor', 'vnc_intrinsic', 'ascending_neuron', 'descending_neuron', 'cb_intrinsic', 'ol_intrinsic']:
    s = sc == g
    print(f'{g:<18s} n={s.sum():>6,}  stationary mass {pi[s].sum():.3f}  median in-syn {np.median(inw[s]):7.0f}  median out-syn {np.median(outw[s]):6.0f}')
def lam2(mask, label):
    idx = np.flatnonzero(mask); B = A[idx][:, idx]
    from scipy.sparse.csgraph import connected_components
    nc, lab = connected_components(B, directed=True, connection='strong'); big = np.bincount(lab).argmax(); idx2 = np.flatnonzero(lab == big)
    B = B[idx2][:, idx2]; rs = np.asarray(B.sum(axis=1)).ravel(); P = (sp.diags(1/rs) @ B).tocsr()
    ev = sla.eigs(P, k=3, which='LM', tol=1e-7, maxiter=6000, ncv=40, return_eigenvectors=False); ev = ev[np.argsort(-np.abs(ev))]
    print(f'{label:<45s} SCC {B.shape[0]:>7,}  lam2 = {abs(ev[1]):.4f}  (gap {1-abs(ev[1]):.4f})')
lam2(np.ones(n, bool), 'full giant SCC (reference)')
lam2(sc != 'vnc_motor', 'motor neurons removed')
lam2(~np.isin(sc, ['ascending_neuron','descending_neuron','sensory_ascending']), 'ascending + descending neurons removed')
lam2(~np.isin(sc, ['vnc_motor','vnc_intrinsic','vnc_sensory','vnc_efferent','vnc_endocrine','vnc_tbc','vnc_sensory_tbc']), 'brain only (VNC removed)')
