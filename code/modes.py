"""Items 1 and 2: localisation of the leading signed modes; anatomy of the slow unsigned mode."""
import os
DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')), '')
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import pyarrow.feather as feather

SIGN = {'acetylcholine': 1, 'gaba': -1, 'glutamate': -1, 'histamine': -1}
BRAIN = {'ol_intrinsic', 'cb_intrinsic', 'visual_projection', 'visual_centrifugal',
         'ol_sensory', 'cb_sensory', 'endocrine'}
VNC = {'vnc_intrinsic', 'vnc_sensory', 'vnc_motor'}
BRIDGE = {'ascending_neuron', 'descending_neuron'}
OL = {'ol_intrinsic', 'ol_sensory'}          # optic lobe proper
OLBRIDGE = {'visual_projection', 'visual_centrifugal'}

d = np.load(DATA + 'neuron_edges.npz')
pre, post, w, body_ids = d['pre'], d['post'], d['w'], d['body_ids']
keep = np.load(DATA + 'scc_keep.npy'); pi = np.load(DATA + 'scc_pi.npy')
n_all = len(body_ids); inv = -np.ones(n_all, dtype=np.int64); inv[keep] = np.arange(len(keep)); n = len(keep)

ann = feather.read_table(DATA + 'body-annotations.feather',
                         columns=['bodyId', 'superclass', 'class', 'type', 'somaSide', 'statusLabel']).to_pandas()
ann = ann.set_index('bodyId').reindex(body_ids[keep]).astype(object)
sc = ann['superclass'].fillna('?').to_numpy(); typ = ann['type'].fillna('?').to_numpy()
side = ann['somaSide'].fillna('?').to_numpy(); status = ann['statusLabel'].fillna('?').to_numpy()
print('superclasses in giant SCC:', dict(zip(*np.unique(sc, return_counts=True))), '\n')

nt = feather.read_table(DATA + 'body-nt.feather', columns=['body', 'consensus_nt']).to_pandas()
nt = nt.set_index('body').reindex(body_ids[keep])['consensus_nt'].fillna('unknown').to_numpy()
sign = np.array([SIGN.get(c, 0) for c in nt], dtype=np.float64)

m = (inv[pre] >= 0) & (inv[post] >= 0)
p2, q2, w2 = inv[pre[m]], inv[post[m]], w[m].astype(np.float64)

# ---------------- item 1: signed map, leading eigenvectors ----------------
ms = sign[p2] != 0
A = sp.csr_matrix((w2[ms], (p2[ms], q2[ms])), shape=(n, n))
cs = np.asarray(A.sum(axis=0)).ravel(); cs[cs == 0] = 1
M = (sp.diags(sign) @ (A @ sp.diags(1.0 / cs))).T.tocsr(); del A
vals, vecs = sla.eigs(M, k=10, which='LM', tol=1e-8, maxiter=6000, ncv=60)
o = np.argsort(-np.abs(vals)); vals, vecs = vals[o], vecs[:, o]
print('=== item 1: leading modes of the signed map (giant SCC) ===')
print('mode  lam        PR (eff. #neurons)  mass in top-10  top carriers (type, side, NT)')
for k in range(len(vals)):
    v = np.abs(vecs[:, k]) ** 2; v /= v.sum()
    pr = 1.0 / (v ** 2).sum(); top = np.argsort(-v)[:10]
    carriers = ', '.join(f'{typ[i]}/{side[i]}/{nt[i][:3]}' for i in top[:4])
    print(f'{k:3d}  {vals[k].real:+.4f}   {pr:9.1f}          {v[top].sum():.3f}          {carriers}')
np.savez(DATA + 'signed_modes.npz', vals=vals, vecs=vecs)
del M, vecs

# ---------------- item 2: unsigned chain, slow mode anatomy ----------------
A = sp.csr_matrix((w2, (p2, q2)), shape=(n, n))
rs = np.asarray(A.sum(axis=1)).ravel()
P = (sp.diags(1.0 / rs) @ A).tocsr(); del A
valsL, vecsL = sla.eigs(P.T.tocsr(), k=3, which='LM', tol=1e-8, maxiter=6000, ncv=40)
o = np.argsort(-np.abs(valsL)); valsL, vecsL = valsL[o], vecsL[:, o]
u = vecsL[:, 1].real; u /= np.abs(u).sum()   # left eigenvector for lam2 (slow deviation of a distribution)
print(f'\n=== item 2: slow mode lam2 = {valsL[1].real:.5f} (left eigenvector) ===')
pos, neg = u > 0, u < 0
groups = {'brain (OL+CB)': np.isin(sc, list(BRAIN)), 'VNC': np.isin(sc, list(VNC)),
          'asc/desc bridge': np.isin(sc, list(BRIDGE)), 'optic lobe proper': np.isin(sc, list(OL)),
          'OL<->CB bridge': np.isin(sc, list(OLBRIDGE)),
          'central brain intrinsic': sc == 'cb_intrinsic'}
print('group                     #neurons   share of + mass   share of - mass   stationary mass')
for name, g in groups.items():
    print(f'{name:<24s} {g.sum():>9,}   {u[g & pos].sum()/u[pos].sum():14.3f}   '
          f'{-u[g & neg].sum()/-u[neg].sum():15.3f}   {pi[g].sum():14.3f}')
# concentration of the slow mode by cell type
v = np.abs(u); top = np.argsort(-v)[:15]
print('\ntop carriers of the slow mode: type/superclass/side  (signed weight)')
for i in top: print(f'  {typ[i]:<14s} {sc[i]:<18s} {side[i]:<3s} {u[i]:+.5f}')

# stationary flow across the candidate bottlenecks (escape probability per hop)
def flow(gA, gB):
    mask = gA[p2] & gB[q2]
    return (pi[p2[mask]] * P[p2[mask], q2[mask]].A1).sum() / pi[gA].sum()
B, V, R = groups['brain (OL+CB)'], groups['VNC'], groups['asc/desc bridge']
brainside, vncside = B | R, V
print('\nper-hop stationary escape probabilities:')
print(f'  brain(+bridge) -> VNC : {flow(brainside, vncside):.4f}    VNC -> brain(+bridge): {flow(vncside, brainside):.4f}')
OLp, CBp = groups['optic lobe proper'], ~groups['optic lobe proper']
print(f'  optic lobe -> rest    : {flow(OLp, CBp):.4f}    rest -> optic lobe   : {flow(CBp, OLp):.4f}')
print(f'  two-block prediction for lam2 (1 - p - q): brain/VNC {1-flow(brainside,vncside)-flow(vncside,brainside):.4f}, '
      f'OL/rest {1-flow(OLp,CBp)-flow(CBp,OLp):.4f}')
