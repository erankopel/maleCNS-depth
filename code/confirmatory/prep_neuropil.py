"""Per-neuron neuropil profile from syn-partners, for the descriptive part (c) of H2 (DEVIATIONS.md C2).

Reads only body_pre, body_post and primary_post; never a confidence column. Every synapse-partner
row adds one count to (body_pre, neuropil) and one to (body_post, neuropil), where neuropil is the
row's primary_post. A neuron's primary neuropil is the one holding most of its counts (ties: first
in the file's vocabulary order).
Output: data/neuron_neuropil.npz with counts (166,700 x 146, int32), vocabulary, primary index.
"""
import numpy as np
from common import DATA, log, load_edges, open_synpartners

_, _, _, body_ids = load_edges()
n = len(body_ids)
rd = open_synpartners(['body_pre', 'body_post', 'primary_post'])    # no confidence column is decoded
vocab = rd.get_batch(0).column('primary_post').dictionary.to_pylist()
R = len(vocab)
counts = np.zeros(n * R, dtype=np.int64)
buf, buffered, rows_total, rows_in_n = [], 0, 0, 0


def flush():
    global buf, buffered
    if buf:
        counts[:] += np.bincount(np.concatenate(buf), minlength=n * R)
    buf, buffered = [], 0


for bi in range(rd.num_record_batches):
    b = rd.get_batch(bi)
    assert b.schema.names == ['body_pre', 'body_post', 'primary_post']
    col = b.column('primary_post')
    d = col.dictionary.to_pylist()
    codes = col.indices.to_numpy(zero_copy_only=False).astype(np.int64)
    if d != vocab:                                   # remap a batch-local dictionary
        remap = np.array([vocab.index(x) if x in vocab else -1 for x in d])
        if (remap < 0).any():
            raise ValueError(f'batch {bi}: neuropil not in vocabulary')
        codes = remap[codes]
    rows_total += b.num_rows
    for name in ('body_pre', 'body_post'):
        body = b.column(name).to_numpy(zero_copy_only=False)
        i = np.searchsorted(body_ids, body)
        np.clip(i, 0, n - 1, out=i)
        ok = body_ids[i] == body
        if name == 'body_pre':
            rows_in_n += int(ok.sum())
        buf.append(i[ok] * R + codes[ok])
        buffered += int(ok.sum())
    if buffered > 40_000_000:
        flush()
    if (bi + 1) % 500 == 0:
        log(f'  batch {bi + 1}/{rd.num_record_batches}')
flush()
counts = counts.reshape(n, R).astype(np.int32)
primary = counts.argmax(axis=1)
primary[counts.sum(axis=1) == 0] = -1
log(f'rows {rows_total:,}; rows with an annotated presynaptic neuron {rows_in_n:,}')
log(f'neurons with no synapse in the table: {(primary < 0).sum():,}')
for roi in ('LOP(L)', 'LOP(R)', 'EB'):
    log(f'  primary neuropil {roi:<7s}: {(primary == vocab.index(roi)).sum():,} neurons')
np.savez_compressed(DATA + 'neuron_neuropil.npz', counts=counts, vocab=np.array(vocab), primary=primary,
                    body_ids=body_ids)
log('saved ' + DATA + 'neuron_neuropil.npz')
