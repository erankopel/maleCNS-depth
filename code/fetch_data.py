"""Download the three male-cns v1.0 flat-connectome inputs and verify them against the md5
checksums frozen in PREREG_maleCNS_depth_2026-09-19.md, Section 1. Public bucket, no token.
Files already present with the right size and md5 are not downloaded again."""
import hashlib
import os
import sys
import urllib.request

DATA = os.path.join(os.environ.get('MALECNS_DATA', os.path.join(os.path.dirname(
    os.path.abspath(__file__)), '..', 'data')), '')
BASE = 'https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/'
FILES = [  # local name, remote name, md5 (PREREG Section 1), bytes
    ('body-annotations.feather', 'body-annotations-male-cns-v1.0-minconf-0.5.feather',
     '50a7718770c57220f160ba4f431ab89e', 14483314),
    ('body-nt.feather', 'body-neurotransmitters-male-cns-v1.0.feather',
     '3d842b12fe5c49eefade528d7dd24a1f', 43282834),
    ('connectome-weights.feather', 'connectome-weights-male-cns-v1.0-minconf-0.5.feather',
     'f30e9dcca25cfd021bf1e7b3d975599e', 1051241946),
]


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


os.makedirs(DATA, exist_ok=True)
bad = 0
for local, remote, want, size in FILES:
    path = DATA + local
    if os.path.exists(path) and os.path.getsize(path) == size and md5(path) == want:
        print(f'ok (already present)  {local}')
        continue
    print(f'downloading {remote} ({size / 1e6:.0f} MB) ...', flush=True)
    tmp = path + '.part'
    urllib.request.urlretrieve(BASE + remote, tmp)
    got = md5(tmp)
    if got != want:
        print(f'CHECKSUM MISMATCH  {local}: got {got}, frozen {want}; kept as {tmp}')
        bad += 1
        continue
    os.replace(tmp, path)
    print(f'ok (downloaded)       {local}')
print('all inputs match the pre-registration' if not bad else f'{bad} file(s) do NOT match')
sys.exit(1 if bad else 0)
