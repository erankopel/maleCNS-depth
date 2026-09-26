"""Synthetic stand-in for the male-CNS inputs, used only to test the confirmatory scripts.

Writes, in the same formats as the real data: neuron_edges.npz, body-annotations.feather,
body-nt.feather and syn-partners.feather (several record batches, dictionary-encoded neuropil,
float32 confidences, rows on non-neuron bodies to be filtered out). A brain block and a cord block
are joined only through ascending/descending-type neurons, as in the real CNS.
Options: residual_neck=False replaces the two small neck classes (efferent_ascending,
sensory_descending) by intrinsic neurons, so that removing the slowmode_check.py set decouples brain
and cord (tests the C1 rule); high_conf=True draws confidences close to 1, so that the H4 robust
bounds are not trivially zero (tests C4 at the scored point).
Usage: python make_synthetic.py OUTDIR [seed]
"""
import os
import sys
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as feather
import pyarrow.ipc as ipc


def make(out, seed=0, residual_neck=True, high_conf=False):
    rng = np.random.default_rng(seed)
    os.makedirs(out, exist_ok=True)
    small = (['efferent_ascending'] * 10 + ['sensory_descending'] * 10 if residual_neck
             else ['vnc_intrinsic'] * 10 + ['cb_intrinsic'] * 10)
    classes = (['cb_intrinsic'] * 1300 + ['ol_intrinsic'] * 400 + ['cb_sensory'] * 100 +
               ['cb_motor'] * 30 + ['vnc_intrinsic'] * 700 + ['vnc_motor'] * 100 +
               ['vnc_sensory'] * 50 + ['ascending_neuron'] * 140 + ['descending_neuron'] * 120 +
               ['sensory_ascending'] * 40 + small)
    N = len(classes)
    sc = np.array(classes, dtype=object)
    rng.shuffle(sc)
    body_ids = np.sort(rng.choice(10 ** 12, N, replace=False)).astype(np.int64)
    brain = np.flatnonzero(np.isin(sc, ['cb_intrinsic', 'ol_intrinsic', 'cb_sensory', 'cb_motor']))
    cord = np.flatnonzero(np.isin(sc, ['vnc_intrinsic', 'vnc_motor', 'vnc_sensory']))
    up = np.flatnonzero(np.isin(sc, ['ascending_neuron', 'sensory_ascending', 'efferent_ascending']))
    down = np.flatnonzero(np.isin(sc, ['descending_neuron', 'sensory_descending']))
    role = np.zeros(N, dtype=np.int8)                # 0 brain, 1 cord, 2 up (to brain), 3 down (to cord)
    role[cord], role[up], role[down] = 1, 2, 3
    pools = [np.concatenate([brain, down]), np.concatenate([cord, up]),
             np.concatenate([brain, brain, down]), np.concatenate([cord, cord, up])]
    targets = {i: pools[role[i]] for i in range(N)}
    pre, post, w = [], [], []
    for i in range(N):
        d = 5 + rng.poisson(20)
        t = rng.choice(targets[i], d, replace=True)
        t = t[t != i]
        pre.append(np.full(len(t), i)); post.append(t)
        w.append(1 + rng.geometric(0.22, len(t)))
    pre, post, w = np.concatenate(pre), np.concatenate(post), np.concatenate(w)
    key = pre.astype(np.int64) * N + post
    uk, inv = np.unique(key, return_inverse=True)
    w = np.bincount(inv, weights=w).astype(np.int64)
    pre, post = (uk // N).astype(np.int32), (uk % N).astype(np.int32)
    np.savez_compressed(os.path.join(out, 'neuron_edges.npz'), pre=pre, post=post, w=w.astype(np.int32),
                        body_ids=body_ids)

    types = np.array([f'T{k}' for k in rng.integers(0, 300, N)], dtype=object)
    side = rng.choice(['L', 'R'], N)
    ann = pd.DataFrame({'bodyId': np.concatenate([body_ids, body_ids.max() + 1 + np.arange(50)]),
                        'superclass': list(sc) + [None] * 50,
                        'type': list(types) + [None] * 50,
                        'somaSide': list(side) + [None] * 50,
                        'statusLabel': ['Traced'] * N + [None] * 50})
    feather.write_feather(ann, os.path.join(out, 'body-annotations.feather'))
    nts = rng.choice(['acetylcholine', 'gaba', 'glutamate', 'histamine', 'dopamine', 'serotonin',
                      'octopamine', 'unclear'], N, p=[0.6, 0.13, 0.17, 0.04, 0.02, 0.015, 0.015, 0.01])
    feather.write_feather(pd.DataFrame({'body': body_ids, 'consensus_nt': nts,
                                        'predicted_nt_confidence': rng.uniform(0.5, 1, N)}),
                          os.path.join(out, 'body-nt.feather'))

    vocab = ['<unspecified>', 'EB', 'LOP(L)', 'LOP(R)', 'SMP(L)', 'SMP(R)', 'LegNp(T1)(L)']
    home = rng.integers(0, len(vocab), N)
    rows_pre = np.repeat(body_ids[pre], w)
    rows_post = np.repeat(body_ids[post], w)
    rows_postidx = np.repeat(post, w)
    R = len(rows_pre)
    extra = R // 20                                 # rows touching bodies that are not neurons
    fake = body_ids.max() + 1000 + rng.integers(0, 10 ** 6, extra)
    ep = rng.choice(body_ids, extra)
    swap = rng.random(extra) < 0.5
    rows_pre = np.concatenate([rows_pre, np.where(swap, fake, ep)])
    rows_post = np.concatenate([rows_post, np.where(swap, ep, fake)])
    roi = np.where(rng.random(R) < 0.8, home[rows_postidx], rng.integers(0, len(vocab), R))
    roi = np.concatenate([roi, rng.integers(0, len(vocab), extra)]).astype(np.int16)
    T = len(rows_pre)
    if high_conf:
        conf_pre = (1.0 - 0.5 * rng.beta(1, 60, T)).astype(np.float32)
        conf_post = (1.0 - 0.5 * rng.beta(1, 60, T)).astype(np.float32)
    else:
        conf_pre = rng.uniform(0.5, 1.0, T).astype(np.float32)
        conf_post = rng.uniform(0.5, 1.0, T).astype(np.float32)
    tie = rng.random(T) < (0.002 if high_conf else 0.05)   # exact float32 0.7 values test the >= comparison
    conf_post[tie] = np.float32(0.7)
    o = rng.permutation(T)
    table = pa.table({
        'x_pre': pa.array(np.zeros(T, np.int32)), 'y_pre': pa.array(np.zeros(T, np.int32)),
        'z_pre': pa.array(np.zeros(T, np.int32)), 'body_pre': pa.array(rows_pre[o]),
        'conf_pre': pa.array(conf_pre[o]), 'x_post': pa.array(np.zeros(T, np.int32)),
        'y_post': pa.array(np.zeros(T, np.int32)), 'z_post': pa.array(np.zeros(T, np.int32)),
        'body_post': pa.array(rows_post[o]), 'conf_post': pa.array(conf_post[o]),
        'primary_post': pa.DictionaryArray.from_arrays(pa.array(roi[o]), pa.array(vocab))})
    with ipc.new_file(os.path.join(out, 'syn-partners.feather'), table.schema) as wr:
        for b in table.to_batches(max_chunksize=20000):
            wr.write_batch(b)
    return {'N': N, 'edges': len(w), 'synapses': int(w.sum()), 'rows': T}


if __name__ == '__main__':
    info = make(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    print(info)
