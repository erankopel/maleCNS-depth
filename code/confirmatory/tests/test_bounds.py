"""Unit tests of the bounds used by cert.py, h3_depth.py and h4_confidence.py, on small random chains.

Every claim is checked against brute force: exact tau(P^r) from dense powers, and thousands of
chains P' sampled from the perturbation classes (interior points and vertices).
Also demonstrates that the robust bound without the t = 0 term (code/cert.py, DEVIATIONS.md D4)
is violated by chains inside the class. Synthetic data only. Run: python test_bounds.py
"""
import numpy as np

rng = np.random.default_rng(1)
FAIL = []


def check(cond, msg):
    print(('ok    ' if cond else 'FAIL  ') + msg)
    if not cond:
        FAIL.append(msg)


def random_weights(n, density=0.35, wmax=20):
    w = (rng.random((n, n)) < density) * rng.integers(1, wmax + 1, (n, n))
    np.fill_diagonal(w, 0)
    for i in range(n):                                     # a Hamiltonian cycle keeps it irreducible
        w[i, (i + 1) % n] = max(w[i, (i + 1) % n], 1)
    return w.astype(float)


def stoch(w):
    rs = w.sum(axis=1, keepdims=True)
    P = np.where(rs > 0, w / np.where(rs > 0, rs, 1), 1.0 / w.shape[0])  # empty row: any distribution
    return P


def tau(Q):
    n = Q.shape[0]
    return max(0.5 * np.abs(Q[a] - Q[b]).sum() for a in range(n) for b in range(n))


def trajectories(P, a, b, R):
    x = np.zeros(P.shape[0]); x[a] += 1; x[b] -= 1
    xs = [x]
    for _ in range(R):
        x = x @ P
        xs.append(x)
    return xs


n, R = 25, 8
w = random_weights(n)
P = stoch(w)

# ---- 1. certified two-sided bounds on tau(P^r) against the exact value
witness = rng.choice(n, 8, replace=False)
hubs = rng.choice(n, 6, replace=False)
Q = np.eye(n)
ok = True
for r in range(1, R + 1):
    Q = Q @ P
    exact = tau(Q)
    lower = max(0.5 * np.abs(Q[a] - Q[b]).sum() for a in witness for b in witness)
    upper = 1 - Q[:, hubs].min(axis=0).sum()
    ok &= lower <= exact + 1e-12 and exact <= upper + 1e-12
check(ok, 'witness lower bound <= exact tau(P^r) <= Markov upper bound, r = 1..8')

# ---- 2. per-row eps for the confidence class [w_hi, w] (h4_confidence.py)
w_hi = rng.binomial(w.astype(int), 0.75).astype(float)
eps = 1 - w_hi.sum(axis=1) / w.sum(axis=1)
pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
xs = {p: trajectories(P, *p, R) for p in pairs}
lb_ok = {p: [0.5 * np.abs(xs[p][r]).sum() - sum(np.abs(xs[p][t]) @ eps for t in range(r))
             for r in range(R + 1)] for p in pairs}
lb_bug = {p: [0.5 * np.abs(xs[p][r]).sum() - sum(np.abs(xs[p][t]) @ eps for t in range(1, r))
              for r in range(R + 1)] for p in pairs}


def sample_class(kind):
    if kind == 'interior':
        return w_hi + np.floor(rng.random((n, n)) * (w - w_hi + 1))
    return np.where(rng.random((n, n)) < 0.5, w_hi, w)          # a random vertex of the box


viol_ok = viol_bug = tv_bad = 0
for s in range(3000):
    wp = sample_class('interior' if s % 2 else 'vertex')
    Pp = stoch(wp)
    tv = 0.5 * np.abs(Pp - P).sum(axis=1)
    tv_bad += int((tv > eps + 1e-12).sum())
    for p in pairs[::7]:
        x = trajectories(Pp, *p, R)
        for r in range(1, R + 1):
            d = 0.5 * np.abs(x[r]).sum()
            viol_ok += d < lb_ok[p][r] - 1e-12
            viol_bug += d < lb_bug[p][r] - 1e-12
check(tv_bad == 0, 'row total-variation change <= eps_i = 1 - W_hi,i / W_i for 3000 sampled chains')
check(viol_ok == 0, f'robust lower bound with the t = 0 term holds for every sampled chain (violations: {viol_ok})')
check(viol_bug > 0, f'the bound without the t = 0 term (cert.py) is violated by chains in the class '
                    f'({viol_bug} violations): the D4 correction is needed')

# ---- 3. the uniform class C(delta) of PREREG Sec. 4.4
delta = 0.05
eps_u = 2 * delta / (1 - delta)
viol_lo = viol_up = viol_cert = tv_bad = 0
Q = np.linalg.matrix_power(P, 6)
floor6 = Q.min(axis=0).sum()
for s in range(2000):
    eta = rng.uniform(-delta, delta, (n, n))
    if s % 2 == 0:
        eta = np.sign(eta) * delta
    Pp = stoch(w * (1 + eta))
    tv_bad += int((0.5 * np.abs(Pp - P).sum(axis=1) > eps_u + 1e-12).sum())
    Qp = np.linalg.matrix_power(Pp, 6)
    viol_up += tau(Qp) > 1 - ((1 - delta) / (1 + delta)) ** 6 * floor6 + 1e-12
    for p in pairs[::11]:
        x = trajectories(Pp, *p, R)
        for r in range(1, R + 1):
            d = 0.5 * np.abs(x[r]).sum()
            lb = 0.5 * np.abs(xs[p][r]).sum() - eps_u * sum(2 * 0.5 * np.abs(xs[p][t]).sum() for t in range(r))
            lb_cert = 0.5 * np.abs(xs[p][r]).sum() - eps_u * sum(2 * 0.5 * np.abs(xs[p][t]).sum() for t in range(1, r))
            viol_lo += d < lb - 1e-12
            viol_cert += d < lb_cert - 1e-12
check(tv_bad == 0, 'uniform class: row total-variation change <= 2 delta / (1 - delta)')
check(viol_lo == 0, 'uniform class: robust lower bound with t = 0 holds for every sampled chain')
check(viol_up == 0, 'uniform class: robust upper bound at r = 6 holds for every sampled chain')
check(viol_cert > 0, f'uniform class: cert.py\'s bound (no t = 0 term) is violated by sampled chains ({viol_cert} violations)')

print('\nALL TESTS PASSED' if not FAIL else f'\n{len(FAIL)} TEST(S) FAILED')
raise SystemExit(1 if FAIL else 0)
