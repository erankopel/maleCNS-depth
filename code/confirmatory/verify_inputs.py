"""Gate before the confirmatory runs: every input must match the value frozen in DEVIATIONS.md C6.
Exits non-zero on any mismatch, so run_confirmatory.sh stops before h1 to h5."""
import hashlib
import sys
import numpy as np
from common import DATA, log

MD5 = {  # PREREG Sec. 1, plus syn-partners (md5 reported by Google Cloud Storage, REPRODUCTION_2026-09-26.md)
    'body-annotations.feather': '50a7718770c57220f160ba4f431ab89e',
    'body-nt.feather': '3d842b12fe5c49eefade528d7dd24a1f',
    'connectome-weights.feather': 'f30e9dcca25cfd021bf1e7b3d975599e',
    'syn-partners.feather': '58efcf712f8c4d4de5f2ad51e97def76',
}
GRAPH_CONTENT_SHA256 = 'cd6bbde8abf155fa746aaa64afc34f5c5b2454ac7867426a2a4972fc2566b0e8'


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


bad = 0
for name, want in MD5.items():
    got = md5(DATA + name)
    log(f'{"ok  " if got == want else "BAD "} md5 {name}: {got}')
    bad += got != want
d = np.load(DATA + 'neuron_edges.npz')
h = hashlib.sha256()
for k in ('pre', 'post', 'w', 'body_ids'):
    a = np.ascontiguousarray(d[k]); h.update(k.encode()); h.update(str(a.dtype).encode()); h.update(a.tobytes())
ok = h.hexdigest() == GRAPH_CONTENT_SHA256
log(f'{"ok  " if ok else "BAD "} content sha256 of neuron_edges.npz (rebuilt by build_graph.py): {h.hexdigest()}')
bad += not ok
sys.exit(1 if bad else 0)
