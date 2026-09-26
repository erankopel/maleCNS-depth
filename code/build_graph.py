"""Build a compact neuron-level edge list for male-cns v1.0.

Streams the 1.1 GB flat connectome (152M segment-to-segment rows) in Arrow
record batches, keeps only edges whose two endpoints are annotated neurons
(non-null superclass), and stores them as int32 index pairs + int32 weight.
"""
import os
DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')), '')
import numpy as np
import pyarrow as pa
import pyarrow.feather as feather

ANN = DATA + 'body-annotations.feather'
EDGES = DATA + 'connectome-weights.feather'
OUT = DATA + 'neuron_edges.npz'

ann = feather.read_table(ANN, columns=['bodyId', 'superclass', 'class']).to_pandas()
neurons = ann[ann['superclass'].notna()].copy()
body_ids = np.sort(neurons['bodyId'].to_numpy().astype(np.int64))
print(f'annotated bodies      : {len(ann):,}')
print(f'neurons (superclass)  : {len(body_ids):,}')

src = pa.memory_map(EDGES)
rd = pa.ipc.open_file(src)
n_batches = rd.num_record_batches

pre_chunks, post_chunks, w_chunks = [], [], []
total_rows = 0
total_syn_all = 0

for i in range(n_batches):
    b = rd.get_batch(i)
    pre = b.column(0).to_numpy(zero_copy_only=False)
    post = b.column(1).to_numpy(zero_copy_only=False)
    w = b.column(2).to_numpy(zero_copy_only=False)
    total_rows += len(pre)
    total_syn_all += int(w.sum())

    ipre = np.searchsorted(body_ids, pre)
    np.clip(ipre, 0, len(body_ids) - 1, out=ipre)
    ok_pre = body_ids[ipre] == pre

    ipost = np.searchsorted(body_ids, post)
    np.clip(ipost, 0, len(body_ids) - 1, out=ipost)
    ok_post = body_ids[ipost] == post

    m = ok_pre & ok_post
    if m.any():
        pre_chunks.append(ipre[m].astype(np.int32))
        post_chunks.append(ipost[m].astype(np.int32))
        w_chunks.append(w[m].astype(np.int32))
    if (i + 1) % 400 == 0:
        print(f'  batch {i+1}/{n_batches}  kept so far '
              f'{sum(len(c) for c in pre_chunks):,}', flush=True)

pre = np.concatenate(pre_chunks); del pre_chunks
post = np.concatenate(post_chunks); del post_chunks
w = np.concatenate(w_chunks); del w_chunks

print(f'\nsegment-to-segment rows : {total_rows:,}')
print(f'synapses in all rows    : {total_syn_all:,}')
print(f'neuron-to-neuron edges  : {len(pre):,}')
print(f'synapses on those edges : {int(w.sum()):,}')
print(f'weight>=5 edges         : {int((w >= 5).sum()):,}')

np.savez_compressed(OUT, pre=pre, post=post, w=w, body_ids=body_ids)
print(f'saved {OUT}')
