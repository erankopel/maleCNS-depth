"""E4, exploratory (after the confirmatory run; not pre-registered): the same synapse-flow chain on the
FlyWire whole-brain connectome v783 (female, brain only; Dorkenwald et al. 2024, Schlegel et al. 2024;
data doi:10.5281/zenodo.10676866, CC-BY) and on the male CNS restricted to the brain.

Inputs in data/flywire/: proofread_connections_783.feather (md5 f48f972d262323a102aed49af1396b8a),
proofread_root_ids_783.npy (md5 e0e6c19732fd8c7a4e39a2d170105421) and the annotation table
Supplemental_file1_neuron_annotations.tsv from github.com/flyconnectome/flywire_annotations.
FlyWire connections are summed over neuropils per (pre, post) pair; the neuron set is the proofread
root ids. Male brain only: the VNC superclasses of slowmode_check.py removed.
For weight >= 1 and >= 5: giant SCC, |lambda_2|, |lambda_3|, stationary concentration, the
smallest-conductance sweep cut along the slow right eigenvector (with its composition), and certified
witness lower bounds on tau(P^r) (100 witness rows, seed 7) at r = 8, 16, 32, 64.
ASCII, no em dashes.
"""
import json
import os
import sys
import time
import numpy as np
import pandas as pd
import pyarrow.feather as feather
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'confirmatory'))
from common import (DATA, log, load_edges, row_stochastic, stationary, top_eigs, pairwise_tv, witness_sources,
                    annotations, VNC_SCRIPT, SLACK, _json_default)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'exploratory', 'results', '')
FW = DATA + 'flywire/'


def analyse(label, pre, post, w, n, labels):
    out = {'label': label}
    for thresh in (1, 5):
        t0 = time.time()
        m = w >= thresh
        A = sp.csr_matrix((w[m].astype(np.float64), (pre[m], post[m])), shape=(n, n))
        _, lab = connected_components(A, directed=True, connection='strong')
        big = np.flatnonzero(lab == np.bincount(lab).argmax())
        P = row_stochastic(A[big][:, big])
        k = P.shape[0]
        vals, V, _ = top_eigs(P, k=6, ncv=48, tol=1e-8, maxiter=6000, vectors=True)
        pi, it, _ = stationary(P, tol=1e-13, maxit=20000)
        ps = np.sort(pi)[::-1]
        v2 = V[:, 1].real
        o = np.argsort(v2, kind='stable'); rank = np.empty(k, np.int64); rank[o] = np.arange(k)
        rows = np.repeat(np.arange(k), np.diff(P.indptr)); f = pi[rows] * P.data
        ri, rj = rank[rows], rank[P.indices]; fwd = ri < rj
        F = np.cumsum(np.bincount(ri[fwd] + 1, weights=f[fwd], minlength=k + 1)
                      - np.bincount(rj[fwd] + 1, weights=f[fwd], minlength=k + 1))[1:k]
        piS = np.cumsum(pi[o])[:k - 1]
        phi = F / np.minimum(piS, 1 - piS)
        c = int(np.argmin(phi)) + 1
        S = np.zeros(k, bool); S[o[:c]] = True
        small = S if pi[S].sum() <= 0.5 else ~S
        p, q = F[c - 1] / pi[small].sum(), F[c - 1] / pi[~small].sum()
        comp = pd.Series(labels[big][small]).value_counts().head(10).to_dict()
        srcs = witness_sources(k, pi, 90, 10, seed=7)
        X = np.zeros((len(srcs), k)); X[np.arange(len(srcs)), srcs] = 1.0
        PT = P.T.tocsr(); lower = {}
        for r in range(1, 65):
            X = (PT @ X.T).T
            if r in (8, 16, 32, 64):
                lower[r] = float(pairwise_tv(X).max() - SLACK * r)
        res = {'n_scc': int(k), 'edges': int(P.nnz), 'lam2': float(abs(vals[1])), 'lam3': float(abs(vals[2])),
               'top1pct_mass': float(ps[:k // 100].sum()), 'top10pct_mass': float(ps[:k // 10].sum()),
               'entropy_bits': float(-(pi * np.log2(pi)).sum()), 'cut_size': int(small.sum()),
               'cut_mass': float(pi[small].sum()), 'cut_conductance': float(phi[c - 1]),
               'two_block_prediction': float(1 - p - q), 'cut_composition': {str(a): int(b) for a, b in comp.items()},
               'certified_lower': lower}
        out[f'w{thresh}'] = res
        log(f'{label}, w>={thresh}: SCC {k:,}, |lam2| {res["lam2"]:.6f}, |lam3| {res["lam3"]:.6f}, top 1% mass '
            f'{res["top1pct_mass"]:.3f}; cut {res["cut_size"]:,} neurons, mass {res["cut_mass"]:.4f}, '
            f'conductance {res["cut_conductance"]:.2e}, two-block {res["two_block_prediction"]:.6f}; '
            f'composition {res["cut_composition"]}; lower bounds {lower}  [{time.time() - t0:.0f} s]')
    return out


WHICH = sys.argv[1] if len(sys.argv) > 1 else 'both'
os.makedirs(OUT, exist_ok=True)
# ---------------------------------------------------------------- FlyWire v783
def run_flywire():
    ids = np.sort(np.load(FW + 'proofread_root_ids_783.npy').astype(np.int64))
    t = feather.read_table(FW + 'proofread_connections_783.feather', columns=['pre_pt_root_id', 'post_pt_root_id', 'syn_count'])
    a, b, s = (t.column(c).to_numpy() for c in ('pre_pt_root_id', 'post_pt_root_id', 'syn_count'))
    ia, ib = np.searchsorted(ids, a), np.searchsorted(ids, b)
    ok = (ia < len(ids)) & (ib < len(ids))
    ok[ok] &= (ids[ia[ok]] == a[ok]) & (ids[ib[ok]] == b[ok])
    key = ia[ok].astype(np.int64) * len(ids) + ib[ok]
    u, inv = np.unique(key, return_inverse=True)
    wf = np.bincount(inv, weights=s[ok]).astype(np.int64)
    pf, qf = (u // len(ids)).astype(np.int64), (u % len(ids)).astype(np.int64)
    log(f'FlyWire v783: {len(ids):,} proofread neurons, {len(u):,} connected pairs, {int(wf.sum()):,} synapses '
        f'({int((~ok).sum()):,} rows outside the proofread set dropped)')
    fa = pd.read_csv(FW + 'Supplemental_file1_neuron_annotations.tsv', sep='\t', low_memory=False,
                     usecols=['root_id', 'super_class', 'cell_type'])
    fa = fa.set_index('root_id').reindex(ids)
    flabels = (fa['super_class'].fillna('?').astype(str) + ':' + fa['cell_type'].fillna('?').astype(str)).to_numpy()
    fly = analyse('FlyWire v783 (female brain)', pf, qf, wf, len(ids), flabels)
    json.dump(fly, open(OUT + 'e4_flywire_fw.json', 'w'), indent=1, default=_json_default)
    log('saved ' + OUT + 'e4_flywire_fw.json')


# ---------------------------------------------------------------- male CNS, brain only
def run_male():
    pre, post, w, body_ids = load_edges()
    ann = annotations(body_ids, cols=('superclass', 'type'))
    brain = ~np.isin(ann['superclass'], VNC_SCRIPT)
    idx = np.flatnonzero(brain); pos = -np.ones(len(body_ids), np.int64); pos[idx] = np.arange(len(idx))
    m = brain[pre] & brain[post]
    mlabels = (ann['superclass'].astype(str) + ':' + ann['type'].astype(str))[idx]
    male = analyse('male CNS, brain only', pos[pre[m]], pos[post[m]], w[m], len(idx), mlabels)
    json.dump(male, open(OUT + 'e4_flywire_male.json', 'w'), indent=1, default=_json_default)
    log('saved ' + OUT + 'e4_flywire_male.json')


if WHICH in ('both', 'flywire'):
    run_flywire()
if WHICH in ('both', 'male'):
    run_male()
