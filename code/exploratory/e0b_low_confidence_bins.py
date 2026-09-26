"""Exploratory (post-confirmatory, not pre-registered). Finer confidence bins below 0.55 for E2b.

For every edge of G counts the syn-partners rows whose min(conf_pre, conf_post) falls in each of the
five bins [0.50, 0.51), [0.51, 0.52), [0.52, 0.53), [0.53, 0.54), [0.54, 0.55) (float32 comparisons,
as in H4). Output: data/edge_conf_hist_low.npy (int32, edges x 5, edges in the CSR order of G); its
row sums equal the first column of data/edge_conf_hist.npy. ASCII, no em dashes.
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'confirmatory'))
from common import DATA, log, load_edges, giant_scc, open_synpartners

EDGES = np.array([0.51, 0.52, 0.53, 0.54], dtype=np.float32)
pre, post, w, body_ids = load_edges()
N = len(body_ids)
keep, A = giant_scc(pre, post, w, N, thresh=1)
A.sort_indices()
n = len(keep)
inv = -np.ones(N, dtype=np.int64); inv[keep] = np.arange(n)
ekeys = np.repeat(np.arange(n, dtype=np.int64), np.diff(A.indptr)) * n + A.indices.astype(np.int64)
E = len(ekeys)
hist = np.zeros(E * 5, dtype=np.int32)
rd = open_synpartners(['body_pre', 'body_post', 'conf_pre', 'conf_post'])


def to_g(body):
    i = np.searchsorted(body_ids, body); np.clip(i, 0, N - 1, out=i)
    g = inv[i]; g[body_ids[i] != body] = -1
    return g


buf, nbuf, t0 = [], 0, time.time()


def flush():
    global buf, nbuf
    if buf:
        u, c = np.unique(np.concatenate(buf), return_counts=True)
        hist[u] += c.astype(np.int32)
    buf, nbuf = [], 0


for bi in range(rd.num_record_batches):
    b = rd.get_batch(bi)
    gi = to_g(b.column('body_pre').to_numpy(zero_copy_only=False))
    gj = to_g(b.column('body_post').to_numpy(zero_copy_only=False))
    ok = (gi >= 0) & (gj >= 0)
    key = gi[ok] * n + gj[ok]
    pos = np.searchsorted(ekeys, key)
    cm = np.minimum(b.column('conf_pre').to_numpy(zero_copy_only=False),
                    b.column('conf_post').to_numpy(zero_copy_only=False))[ok]
    low = cm < np.float32(0.55)
    buf.append(pos[low] * 5 + np.searchsorted(EDGES, cm[low], side='right')); nbuf += int(low.sum())
    if nbuf > 25_000_000:
        flush()
    if (bi + 1) % 1000 == 0:
        log(f'  batch {bi + 1}/{rd.num_record_batches}  [{time.time() - t0:.0f} s]')
flush()
hist = hist.reshape(E, 5)
coarse = np.load(DATA + 'edge_conf_hist.npy', mmap_mode='r')
assert (np.load(DATA + 'edge_conf_keys.npy', mmap_mode='r') == ekeys).all()
for e0 in range(0, E, 2_000_000):
    assert (hist[e0:e0 + 2_000_000].sum(axis=1) == np.asarray(coarse[e0:e0 + 2_000_000, 0])).all()
np.save(DATA + 'edge_conf_hist_low.npy', hist)
log(f'saved per-edge counts below 0.55: {E:,} edges x 5 bins; synapses per bin: {hist.sum(axis=0).tolist()}')
