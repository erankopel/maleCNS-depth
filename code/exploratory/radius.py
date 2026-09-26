"""Sharper certified row radius for the confidence class (exploratory, post-confirmatory).

Row i has targets j with weights w_j and high-confidence parts hi_j (0 <= hi_j <= w_j); the class
allows any w' with hi <= w' <= w. For every member, TV(w/W, w'/W') = max_S (p(S) - p'(S)), and for a
fixed set S the smallest p'(S) over the box is attained at w' = hi on S and w on S^c, so

    sup over the class of TV = max_S g(S),   g(S) = w(S)/W - hi(S) / (W - l(S)),   l = w - hi.

For a fixed removed weight C = l(S) this is a knapsack problem. Its fractional relaxation is solved by
taking targets in increasing order of q_j = hi_j / w_j: the value per unit of removed weight,
(1/W - q_j/y) / (1 - q_j) with y = W - C <= W, decreases in q_j for every y. Hence

    eps*_i = max over the greedy path of  G(x) = sum_j x_j w_j / W - sum_j x_j hi_j / (W - sum_j x_j l_j)

is an upper bound on the sup (it can exceed the sup, never fall below it). On each segment of the path
(prefix k, fraction phi of the next target) G is concave in phi, with stationary point
phi* = (Y - sqrt(W (h Y + l H) / w)) / l, Y = W - L_k, H = hi of the prefix. Always
eps*_i <= C/W <= eps_i = L_i / W, the radius of Lemma 1. A row without high-confidence synapses may be
emptied, and then any distribution is allowed: eps*_i = 1. For rows with at most DMAX targets the
sup is computed exactly by enumerating the sets S (exact_small), which removes the looseness of the
relaxation where a single target carries much of the row. ASCII, no em dashes.
"""
import numpy as np


def sharp_radius(rows, w, hi, n, margin=1e-12):
    """rows: row index of every edge (any order); w > 0 and 0 <= hi <= w per edge.
    Returns (eps_star, eps), both of length n; eps_star <= eps."""
    rows = np.asarray(rows, dtype=np.int64)
    w = np.asarray(w, dtype=np.float64)
    hi = np.asarray(hi, dtype=np.float64)
    assert (w > 0).all() and (hi >= 0).all() and (hi <= w).all()
    W = np.bincount(rows, weights=w, minlength=n)
    Whi = np.bincount(rows, weights=hi, minlength=n)
    safeW = np.where(W > 0, W, 1.0)
    eps = np.where(W > 0, 1.0 - Whi / safeW, 0.0)

    order = np.lexsort((hi / w, rows))                   # by row, then q = hi/w ascending
    r_, w_, h_ = rows[order], w[order], hi[order]
    l_ = w_ - h_
    start = np.flatnonzero(np.r_[True, r_[1:] != r_[:-1]])
    seg = np.repeat(np.arange(len(start)), np.diff(np.r_[start, len(r_)]))

    def excl(x):                                           # exclusive prefix sum within each row
        c = np.cumsum(x) - x
        return c - c[start][seg]

    A, H, Lk = excl(w_), excl(h_), excl(l_)
    Wr = W[r_]
    Y = Wr - Lk                                            # weight left after removing the prefix
    with np.errstate(divide='ignore', invalid='ignore'):
        root = np.sqrt(Wr * (h_ * Y + l_ * H) / w_)
        phi = np.where(l_ > 0, (Y - root) / np.where(l_ > 0, l_, 1.0), 0.0)
    phi = np.clip(np.nan_to_num(phi, nan=0.0), 0.0, 1.0)
    best = np.full(len(r_), -np.inf)
    for f in (phi, np.zeros_like(phi), np.ones_like(phi)):
        den = Y - f * l_
        ok = den > 0
        G = np.where(ok, (A + f * w_) / Wr - (H + f * h_) / np.where(ok, den, 1.0), -np.inf)
        best = np.maximum(best, G)
    eps_star = np.zeros(n)
    eps_star[r_[start]] = np.maximum(np.maximum.reduceat(best, start), 0.0)
    eps_star = np.where(Whi > 0, eps_star, np.where(W > 0, 1.0, 0.0))
    eps_star = np.minimum(eps_star + margin, eps)          # rounding margin; never above Lemma 1
    return eps_star, eps


DMAX = 12


def exact_small(rows, w, hi, n, dmax=DMAX):
    """Exact sup of the row change for rows with at most dmax targets (NaN for the other rows)."""
    rows = np.asarray(rows, dtype=np.int64)
    w = np.asarray(w, dtype=np.float64)
    hi = np.asarray(hi, dtype=np.float64)
    order = np.argsort(rows, kind='stable')
    r_, w_, h_ = rows[order], w[order], hi[order]
    deg = np.bincount(r_, minlength=n)
    out = np.full(n, np.nan)
    start = np.r_[0, np.cumsum(deg)[:-1]]
    for d in range(1, dmax + 1):
        R = np.flatnonzero(deg == d)
        if len(R) == 0:
            continue
        idx = start[R][:, None] + np.arange(d)[None, :]
        wd, hd = w_[idx], h_[idx]                      # rows x d
        Wd = wd.sum(axis=1, keepdims=True)
        masks = ((np.arange(1, 2 ** d)[:, None] >> np.arange(d)[None, :]) & 1).astype(np.float64)
        best = np.full(len(R), -np.inf)
        for m0 in range(0, len(masks), 512):            # chunks of subsets
            M = masks[m0:m0 + 512]
            wS, hS = wd @ M.T, hd @ M.T                 # rows x subsets
            den = hS + (Wd - wS)
            ok = den > 0
            g = np.where(ok, wS / Wd - hS / np.where(ok, den, 1.0), -np.inf)
            best = np.maximum(best, g.max(axis=1))
        out[R] = best
    return out


def certified_radius(rows, w, hi, n, dmax=DMAX, margin=1e-12):
    """eps*: exact sup for rows with <= dmax targets, the relaxation for the others; and eps (Lemma 1)."""
    rel, eps = sharp_radius(rows, w, hi, n, margin=margin)
    ex = exact_small(rows, w, hi, n, dmax)
    Whi = np.bincount(np.asarray(rows, dtype=np.int64), weights=np.asarray(hi, dtype=np.float64), minlength=n)
    use = ~np.isnan(ex) & (Whi > 0)
    out = rel.copy()
    out[use] = np.minimum(np.maximum(ex[use], 0.0) + margin, eps[use])
    return out, eps
