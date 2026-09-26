"""Test of radius.py (exploratory): on random rows with 1 to 15 partners the certified radius is never
below the exact sup of the row change over the confidence class (brute force over all partner sets),
never above the radius of Lemma 1, and exact for rows with at most 12 partners. ASCII, no em dashes."""
import itertools, sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from radius import sharp_radius, exact_small, certified_radius
rng = np.random.default_rng(2)
rows, ws, his = [], [], []
n = 2500
for i in range(n):
    d = int(rng.integers(1, 16))
    w = rng.integers(1, 12, size=d).astype(float)
    mode = rng.integers(3)
    if mode == 0: hi = np.array([rng.integers(0, int(x) + 1) for x in w]).astype(float)
    elif mode == 1: hi = np.where(rng.random(d) < 0.7, w, 0.0)
    else: hi = np.floor(w * rng.uniform(0.6, 1.0))
    rows += [i] * d; ws += list(w); his += list(hi)
rows = np.array(rows); ws = np.array(ws); his = np.array(his)
perm = rng.permutation(len(rows)); rows, ws, his = rows[perm], ws[perm], his[perm]
cr, e = certified_radius(rows, ws, his, n)
rel, _ = sharp_radius(rows, ws, his, n)
bad = 0; exact_rows = 0; gaps = []
for i in range(n):
    m = rows == i; w = ws[m]; hi = his[m]; W = w.sum(); d = len(w)
    if hi.sum() == 0: brute = 1.0
    else:
        brute = 0.0
        for S in itertools.product([0, 1], repeat=d):
            S = np.array(S, bool)
            if not S.any(): continue
            brute = max(brute, w[S].sum() / W - hi[S].sum() / (hi[S].sum() + W - w[S].sum()))
    if cr[i] < brute - 1e-12 or cr[i] > e[i] + 1e-15 or rel[i] < brute - 1e-12: bad += 1
    if d <= 12:
        exact_rows += 1
        assert abs(cr[i] - min(brute + 1e-12, e[i])) < 1e-9 or hi.sum() == 0, (i, cr[i], brute)
    else:
        gaps.append(cr[i] - brute)
print('violations', bad, '| exact rows', exact_rows, '| relaxation rows', len(gaps), 'max gap %.4f mean gap %.4f' % (max(gaps), np.mean(gaps)))
