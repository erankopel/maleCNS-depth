"""H5 (sign convention), PREREG Sec. 5, with DEVIATIONS.md C5 and C7. A sensitivity result.

Repeats Sec. 4.3 (modes.py item 1: signed map on G, postsynaptic normalisation over partners with
nonzero sign, 10 leading modes) under the base convention and under Pospisil's (dopamine +,
serotonin and octopamine -).
Criteria (frozen): the spectral radius changes by less than 0.02, and the top-six carriers are
unchanged. Carriers (C5): for each of the six leading modes, the smallest set of neurons holding
>= 90 percent of the mode's power; "unchanged" means the union of the six sets, as body IDs, is
identical under both conventions. Pass iff both. Modes are right eigenvectors; a conjugate pair has
one power distribution, so the union does not depend on which member of a pair is listed first.
"""
import sys
import numpy as np
from common import (log, load_edges, giant_scc, consensus_nt, signs, signed_map, top_eigs,
                    mode_stats, carrier_set, annotations, save_result, expect, NotDetermined)

pre, post, w, body_ids = load_edges()
keep, A = giant_scc(pre, post, w, len(body_ids), thresh=1)
expect(len(body_ids), 'G', len(keep))
bid = body_ids[keep]
nt = consensus_nt(bid)
ann = annotations(bid, cols=('type', 'somaSide'))


def run(convention):
    s = signs(nt, convention)
    M, _ = signed_map(A, s, floor=None)
    vals, vecs, info = top_eigs(M, k=10, ncv=60, tol=1e-8, maxiter=6000, vectors=True)
    p, pr = mode_stats(vecs)
    sets = [carrier_set(p[:, k], 0.90) for k in range(6)]
    union = sorted(set(int(bid[i]) for st in sets for i in st))
    log(f'\n--- {convention}: {int((s > 0).sum()):,} +, {int((s < 0).sum()):,} -, {int((s == 0).sum()):,} excluded; '
        f'{M.nnz:,} edges ---')
    for k in range(len(vals)):
        top = np.argsort(-p[:, k])[:4]
        car = ', '.join(f'{ann["type"][i]}/{ann["somaSide"][i]}/{nt[i][:3]}' for i in top)
        extra = f'   90% set: {len(sets[k])} neurons' if k < 6 else ''
        log(f'{k:3d}  {vals[k].real:+.5f}{vals[k].imag:+.5f}i  PR {pr[k]:6.1f}   {car}{extra}')
    return {'convention': convention, 'rho': float(abs(vals[0])), 'vals': [complex(v) for v in vals],
            'pr': pr.tolist(), 'carrier_union_body_ids': union,
            'carrier_sets_body_ids': [[int(bid[i]) for i in st] for st in sets],
            'carrier_types': sorted(set(str(ann['type'][i]) for st in sets for i in st)), 'arpack': info}


try:
    base = run('base')
    posp = run('pospisil')
except NotDetermined as e:
    log(f'NOT DETERMINED: {e}')
    save_result('h5', {'hypothesis': 'H5', 'verdict': 'NOT DETERMINED', 'not_determined': str(e)})
    log('H5 NOT DETERMINED (sensitivity result)')
    sys.exit(0)
d_rho = abs(posp['rho'] - base['rho'])
same = base['carrier_union_body_ids'] == posp['carrier_union_body_ids']
c1, c2 = d_rho < 0.02, same
verdict = 'PASS' if (c1 and c2) else 'FAIL'
log(f'\nspectral radius: base {base["rho"]:.6f}, Pospisil {posp["rho"]:.6f}, change {d_rho:.6f} (< 0.02): {c1}')
log(f'top-six carrier union unchanged ({len(base["carrier_union_body_ids"])} vs '
    f'{len(posp["carrier_union_body_ids"])} neurons): {c2}')
save_result('h5', {'hypothesis': 'H5', 'criteria': {'rho_change_lt_0.02': c1, 'carriers_unchanged': c2},
                   'rho_change': d_rho, 'verdict': verdict, 'base': base, 'pospisil': posp})
log(f'H5 {verdict} (sensitivity result)')
sys.exit(0)
