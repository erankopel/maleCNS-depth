"""Shared helpers for the confirmatory scripts h1 to h5 (PREREG Sec. 5; DEVIATIONS.md C-entries).

Every numerical choice here reproduces the exploratory scripts of PREREG Sec. 6 unless an entry in
DEVIATIONS.md says otherwise. Nothing in this file reads a synapse-confidence value.
ASCII, no em dashes.
"""
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.abspath(os.environ.get('MALECNS_DATA', os.path.join(HERE, '..', '..', 'data'))), '')
RESULTS = os.path.join(os.path.abspath(os.environ.get('MALECNS_RESULTS', os.path.join(HERE, '..', '..', 'results'))), '')

# Sign conventions. BASE is PREREG Sec. 2 (spectra.py, modes.py); POSPISIL is the sensitivity
# convention of Sec. 2 and H5 (dopamine +, serotonin and octopamine -).
SIGN_BASE = {'acetylcholine': 1, 'gaba': -1, 'glutamate': -1, 'histamine': -1}
SIGN_POSPISIL = {'acetylcholine': 1, 'dopamine': 1, 'gaba': -1, 'glutamate': -1, 'histamine': -1,
                 'serotonin': -1, 'octopamine': -1}

# Superclass groups (DEVIATIONS.md D5 and C1). NECK_SCRIPT is the set removed by slowmode_check.py.
NECK_SCRIPT = ['ascending_neuron', 'descending_neuron', 'sensory_ascending']
NECK_TEXT = ['ascending_neuron', 'descending_neuron']          # the 3,160 of the PREREG 4.2 text
NECK_ALL = ['ascending_neuron', 'descending_neuron', 'sensory_ascending', 'sensory_ascending_tbc',
            'descending_neuron_tbc', 'efferent_ascending', 'efferent_descending', 'sensory_descending']
MOTOR_SCRIPT = ['vnc_motor']
MOTOR_ALL = ['vnc_motor', 'cb_motor']
VNC_SCRIPT = ['vnc_motor', 'vnc_intrinsic', 'vnc_sensory', 'vnc_efferent', 'vnc_endocrine', 'vnc_tbc',
              'vnc_sensory_tbc']

SLACK = 1e-9          # rounding slack per hop (PREREG Sec. 6)
U = 2.0 ** -53        # unit roundoff of float64
CHECKPOINTS = [1, 2, 4, 8, 16, 32, 48, 64, 96, 128]
ARPACK_SEED = 11      # starting vectors of every scored ARPACK call (C7)
ARPACK_SEED_CHECK = 12  # starting vectors of the unscored cross-check calls (C7)

# Structural facts of the real inputs (PREREG Sec. 2, DEVIATIONS.md C2); checked only on the real data.
REAL_N = 166_700
EXPECTED = {'G': 165_314, 'G5': 157_821, 'F': 342.0}


class NotDetermined(Exception):
    """ARPACK did not converge after the escalation of C7."""


def log(msg):
    print(msg, flush=True)


def is_real(n_neurons):
    return n_neurons == REAL_N


def expect(n_neurons, name, value):
    """Stop if a structural fact of the real data differs from the frozen value (C6)."""
    if is_real(n_neurons) and value != EXPECTED[name]:
        raise SystemExit(f'STRUCTURAL CHECK FAILED: {name} = {value}, frozen value {EXPECTED[name]}')
    if is_real(n_neurons):
        log(f'structural check {name} = {value}: ok')


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def environment():
    import pandas
    import pyarrow
    return {'python': sys.version.split()[0], 'numpy': np.__version__, 'scipy': scipy.__version__,
            'pandas': pandas.__version__, 'pyarrow': pyarrow.__version__,
            'platform': platform.platform(), 'cpu_count': os.cpu_count(),
            'threads': {k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')},
            'date_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}


def _strip(x):
    """Drop run-dependent fields (environment, wall times) before comparing results."""
    if isinstance(x, dict):
        return {k: _strip(v) for k, v in x.items() if k not in ('environment', 'wall_s')}
    if isinstance(x, list):
        return [_strip(v) for v in x]
    return x


def save_result(name, obj):
    """Write results/<name>.json. An existing result is never overwritten (C6). With
    MALECNS_UNSCORED_ONLY=1 (a rerun that only completes the unscored part after an interruption), an
    existing result must be reproduced exactly, number for number, or the run stops."""
    os.makedirs(RESULTS, exist_ok=True)
    path = RESULTS + name + '.json'
    if os.path.exists(path):
        if os.environ.get('MALECNS_UNSCORED_ONLY') != '1':
            raise SystemExit(f'REFUSING TO OVERWRITE {path} (C6: a script that has produced a result is not rerun)')
        new = _strip(json.loads(json.dumps(dict(obj), default=_json_default)))
        old = _strip(json.load(open(path)))
        if new != old:
            raise SystemExit(f'{path}: the rerun did NOT reproduce the saved result exactly; stopping (C6)')
        log(f'{path}: saved result reproduced exactly (C6); not rewritten')
        return path
    obj = dict(obj)
    obj['environment'] = environment()
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=2, default=_json_default)
    os.replace(tmp, path)
    log(f'saved {path}')
    return path


def _json_default(x):
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (complex, np.complexfloating)):
        return [float(np.real(x)), float(np.imag(x))]
    raise TypeError(type(x))


# ----------------------------------------------------------------------------- graph objects

def load_edges():
    d = np.load(DATA + 'neuron_edges.npz')
    return d['pre'], d['post'], d['w'], d['body_ids']


def giant_scc(pre, post, w, n, thresh=1, subset=None, return_labels=False):
    """Giant strongly connected component of the graph with edges w >= thresh.

    subset: optional boolean mask over the n neurons; only neurons in it are kept before the SCC is
    taken (knock-outs). Returns (keep, A): keep are neuron indices (sorted), A is the float64 CSR
    weight matrix restricted to keep, rows = presynaptic. With return_labels, also returns the
    neuron indices of the subset and the SCC label of each.
    """
    m = w >= thresh
    if subset is not None:
        m &= subset[pre] & subset[post]
    A = sp.csr_matrix((w[m].astype(np.float64), (pre[m], post[m])), shape=(n, n))
    if subset is not None:
        idx = np.flatnonzero(subset)
        A = A[idx][:, idx]
    else:
        idx = np.arange(n)
    ncomp, lab = connected_components(A, directed=True, connection='strong')
    big = np.bincount(lab).argmax()
    loc = np.flatnonzero(lab == big)
    if return_labels:
        return idx[loc], A[loc][:, loc].tocsr(), idx, lab
    return idx[loc], A[loc][:, loc].tocsr()


def row_stochastic(A):
    rs = np.asarray(A.sum(axis=1)).ravel()
    if not (rs > 0).all():
        raise ValueError('row with zero out-weight inside an SCC')
    return (sp.diags(1.0 / rs) @ A).tocsr()


def stationary(P, tol=1e-13, maxit=10000):
    """Power iteration as in spectra_scc.py (uniform start, L1 stopping rule). If it has not met the
    tolerance after maxit iterations the last vector is returned and the caller reports it (C3)."""
    m = P.shape[0]
    v = np.full(m, 1.0 / m)
    PT = P.T.tocsr()
    for it in range(1, maxit + 1):
        vn = PT @ v
        vn /= vn.sum()
        delta = np.abs(vn - v).sum()
        v = vn
        if delta < tol:
            return v, it, delta
    return v, maxit, delta


def top_eigs(M, k, ncv, tol=1e-8, maxiter=6000, vectors=False, seed=ARPACK_SEED):
    """ARPACK, which='LM', seeded starting vector, sorted by decreasing modulus and then by decreasing
    imaginary part (so a conjugate pair is always in the same order).
    Escalation (C7): on non-convergence retry once with ncv = min(n - 1, 2 ncv) and maxiter = 20000;
    if that also fails, raise NotDetermined."""
    n = M.shape[0]
    info = {'k': k, 'ncv': ncv, 'tol': tol, 'maxiter': maxiter, 'seed': seed, 'retry': False}
    t0 = time.time()
    try:
        out = sla.eigs(M, k=k, which='LM', tol=tol, maxiter=maxiter, ncv=ncv,
                       return_eigenvectors=vectors, rng=np.random.default_rng(seed))
    except sla.ArpackNoConvergence:
        ncv2 = int(min(n - 1, 2 * ncv))
        info.update(retry=True, ncv_retry=ncv2, maxiter_retry=20000)
        log(f'  ARPACK did not converge (k={k}, ncv={ncv}); retrying with ncv={ncv2}, maxiter=20000')
        try:
            out = sla.eigs(M, k=k, which='LM', tol=tol, maxiter=20000, ncv=ncv2,
                           return_eigenvectors=vectors, rng=np.random.default_rng(seed))
        except sla.ArpackNoConvergence:
            raise NotDetermined(f'ARPACK did not converge after the retry (k={k}, ncv={ncv2})')
    info['wall_s'] = round(time.time() - t0, 1)
    vals = out[0] if vectors else out
    o = np.lexsort((-np.imag(vals), -np.abs(vals)))
    if vectors:
        return vals[o], out[1][:, o], info
    return vals[o], None, info


def second_modulus(P, k=8, ncv=48, tol=1e-8, maxiter=6000, seed=ARPACK_SEED):
    vals, _, info = top_eigs(P, k, ncv, tol, maxiter, seed=seed)
    return float(abs(vals[1])), vals, info


# ----------------------------------------------------------------------------- annotations

def superclass_of(body_ids):
    import pyarrow.feather as feather
    ann = feather.read_table(DATA + 'body-annotations.feather', columns=['bodyId', 'superclass']).to_pandas()
    return ann.set_index('bodyId').reindex(body_ids)['superclass'].astype(object).fillna('?').to_numpy()


def side_of(superclasses):
    """'cord' for the VNC superclasses of slowmode_check.py, 'neck' for NECK_ALL, 'brain' otherwise (C1)."""
    out = np.full(len(superclasses), 'brain', dtype=object)
    out[np.isin(superclasses, VNC_SCRIPT)] = 'cord'
    out[np.isin(superclasses, NECK_ALL)] = 'neck'
    return out


def annotations(body_ids, cols=('type', 'somaSide', 'superclass', 'statusLabel')):
    import pyarrow.feather as feather
    ann = feather.read_table(DATA + 'body-annotations.feather', columns=['bodyId'] + list(cols)).to_pandas()
    ann = ann.set_index('bodyId').reindex(body_ids).astype(object)
    return {c: ann[c].fillna('?').to_numpy() for c in cols}


def consensus_nt(body_ids):
    import pyarrow.feather as feather
    nt = feather.read_table(DATA + 'body-nt.feather', columns=['body', 'consensus_nt']).to_pandas()
    return nt.set_index('body').reindex(body_ids)['consensus_nt'].astype(object).fillna('unknown').to_numpy()


def signs(nt, convention):
    table = SIGN_BASE if convention == 'base' else SIGN_POSPISIL
    return np.array([table.get(c, 0) for c in nt], dtype=np.float64)


# ----------------------------------------------------------------------------- signed map

def signed_map(A, s, floor=None):
    """M[j, i] = s_i w(i->j) / den_j with den_j = in_j (modes.py; zero input -> 1) or
    den_j = max(in_j, floor) (H2). in_j counts only presynaptic partners with s_i != 0.
    A: CSR on the SCC, rows = presynaptic. Returns (M, in_j)."""
    As = (sp.diags((s != 0).astype(np.float64)) @ A).tocsr()
    As.eliminate_zeros()
    inw = np.asarray(As.sum(axis=0)).ravel()
    if floor is None:
        den = inw.copy()
        den[den == 0] = 1.0
    else:
        den = np.maximum(inw, float(floor))
    M = (sp.diags(s) @ (As @ sp.diags(1.0 / den))).T.tocsr()
    return M, inw


def mode_stats(vecs):
    """Power distribution p = |v|^2 / sum |v|^2 of each right eigenvector; participation ratio 1 / sum p^2."""
    p = np.abs(vecs) ** 2
    p /= p.sum(axis=0, keepdims=True)
    pr = 1.0 / (p ** 2).sum(axis=0)
    return p, pr


def carrier_set(p, share=0.90):
    """Smallest set of neurons holding >= share of one mode's power (indices)."""
    o = np.argsort(-p, kind='stable')
    c = np.cumsum(p[o])
    k = int(np.searchsorted(c, share) + 1)
    return o[:k]


# ----------------------------------------------------------------------------- Dobrushin bounds

def pairwise_tv(X):
    """d(a, b) = (1/2) || X_a - X_b ||_1 for all rows of X (K x n)."""
    K = X.shape[0]
    D = np.zeros((K, K))
    for a in range(K - 1):
        D[a, a + 1:] = 0.5 * np.abs(X[a][None, :] - X[a + 1:]).sum(axis=1)
    return D + D.T


def witness_sources(n, pi, n_random, n_top, seed):
    """As cert.py: n_random distinct random rows (default_rng(seed)), then the n_top rows of largest
    stationary mass. A top row may repeat a random one; the number of distinct rows is logged."""
    rng = np.random.default_rng(seed)
    src = np.concatenate([rng.choice(n, n_random, replace=False), np.argsort(-pi)[:n_top]])
    log(f'witness rows: {len(src)} ({len(np.unique(src))} distinct)')
    return src


def open_synpartners(fields):
    """Memory-mapped syn-partners reader that decodes only the named columns (C4)."""
    import pyarrow as pa
    import pyarrow.ipc as ipc
    src = pa.memory_map(DATA + 'syn-partners.feather')
    names = ipc.open_file(src).schema.names
    return ipc.open_file(pa.memory_map(DATA + 'syn-partners.feather'),
                         options=ipc.IpcReadOptions(included_fields=[names.index(f) for f in fields]))
