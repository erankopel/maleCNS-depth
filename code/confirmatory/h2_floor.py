"""H2 (localisation is a normalisation effect), PREREG Sec. 5, with DEVIATIONS.md C2, C6 and C7.

Signed map on G (weight >= 1 giant SCC) with a per-neuron input floor:
  M_ji = s_i w_ij / max(in_j, F),   in_j = input weight of j from partners with s_i != 0,
  F = median of in_j over all neurons of G (C2; 342 synapses).
ARPACK computes 24 eigenpairs (right eigenvectors); the scored modes are the first 20 in the order of
decreasing modulus, then decreasing imaginary part (C2).
Criteria (frozen): (a) spectral radius < 0.85; (b) at least 10 of the 20 modes have participation
ratio > 20 neurons. Pass iff (a) and (b).
(c), descriptive: whether any of the 20 modes has > 50 percent of its power on neurons whose primary
neuropil (prep_neuropil.py) is LOP(L), LOP(R) or EB (power on the three neuropils counted together).
Sensitivity (not scored, computed after the scored result is saved): F = 700 (the approximate value
the PREREG states) and F = median total input weight.
"""
import sys
import numpy as np
from common import (DATA, log, load_edges, giant_scc, consensus_nt, signs, signed_map, top_eigs,
                    mode_stats, annotations, save_result, expect, NotDetermined)

pre, post, w, body_ids = load_edges()
n = len(body_ids)
keep, A = giant_scc(pre, post, w, n, thresh=1)
log(f'G: {len(keep):,} neurons, {A.nnz:,} edges')
expect(n, 'G', len(keep))
bid = body_ids[keep]
nt = consensus_nt(bid)
s = signs(nt, 'base')
ann = annotations(bid, cols=('type', 'somaSide'))
npz = np.load(DATA + 'neuron_neuropil.npz', allow_pickle=False)
assert (npz['body_ids'] == body_ids).all()
vocab = list(npz['vocab'])
prim = npz['primary'][keep]
lop = np.isin(prim, [vocab.index('LOP(L)'), vocab.index('LOP(R)')])
eb = prim == vocab.index('EB')
log(f'neurons of G with primary neuropil LOP: {lop.sum():,}; EB: {eb.sum():,}')


def run(label, floor):
    M, inw = signed_map(A, s, floor=floor)
    log(f'\n--- {label}: floor F = {floor:.1f} ---')
    try:
        vals, vecs, info = top_eigs(M, k=24, ncv=60, tol=1e-8, maxiter=10000, vectors=True)
    except NotDetermined as e:
        log(f'NOT DETERMINED: {e}')
        return {'label': label, 'floor': float(floor), 'not_determined': str(e)}
    vals, vecs = vals[:20], vecs[:, :20]
    p, pr = mode_stats(vecs)
    modes = []
    log('mode  lambda               |lambda|   PR      LOP+EB power   top carriers (type/side/NT)')
    for k in range(len(vals)):
        top = np.argsort(-p[:, k])[:4]
        share = float(p[lop | eb, k].sum())
        car = ', '.join(f'{ann["type"][i]}/{ann["somaSide"][i]}/{nt[i][:3]}' for i in top)
        log(f'{k:3d}  {vals[k].real:+.5f}{vals[k].imag:+.5f}i  {abs(vals[k]):.5f}  {pr[k]:7.1f}  '
            f'{share:12.3f}   {car}')
        modes.append({'lambda': complex(vals[k]), 'modulus': float(abs(vals[k])), 'pr': float(pr[k]),
                      'lop_share': float(p[lop, k].sum()), 'eb_share': float(p[eb, k].sum()),
                      'lop_eb_share': share, 'top10_body_ids': bid[np.argsort(-p[:, k])[:10]].tolist(),
                      'top4': car})
    return {'label': label, 'floor': float(floor), 'rho': float(abs(vals[0])),
            'n_modes_pr_gt_20': int((pr > 20).sum()),
            'any_mode_lop_eb_gt_half': bool(any(m['lop_eb_share'] > 0.5 for m in modes)),
            'modes': modes, 'arpack': info}


M0, inw = signed_map(A, s, floor=None)
F = float(np.median(inw))
F_total = float(np.median(np.asarray(A.sum(axis=0)).ravel()))
del M0
log(f'median signed-eligible input of G: F = {F:.1f} synapses; median total input: {F_total:.1f}')
expect(n, 'F', F)

primary = run('primary (C2): F = median signed-eligible input', F)
if 'not_determined' in primary:
    a = b = None
    verdict = 'NOT DETERMINED'
else:
    a = primary['rho'] < 0.85
    b = primary['n_modes_pr_gt_20'] >= 10
    verdict = 'PASS' if (a and b) else 'FAIL'
    log(f'\n(a) spectral radius {primary["rho"]:.5f} < 0.85: {a}')
    log(f'(b) modes with PR > 20 among the top 20: {primary["n_modes_pr_gt_20"]} (>= 10): {b}')
    log(f'(c) descriptive: some top-20 mode with > 50% power on LOP/EB neurons: '
        f'{primary["any_mode_lop_eb_gt_half"]}')
save_result('h2', {'hypothesis': 'H2', 'F': F, 'F_total_input': F_total,
                   'criteria': {'a': a, 'b': b, 'c_descriptive': primary.get('any_mode_lop_eb_gt_half')},
                   'verdict': verdict, 'primary': primary,
                   'n_lop_neurons': int(lop.sum()), 'n_eb_neurons': int(eb.sum())})
log(f'H2 {verdict}')

sens = [run('sensitivity: F = 700, the value stated in the PREREG', 700.0),
        run('sensitivity: F = median total input', F_total)]
save_result('h2_sensitivity', {'hypothesis': 'H2', 'scored': False, 'rows': sens})
sys.exit(0)
