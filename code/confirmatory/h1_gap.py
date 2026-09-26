"""H1 (neck-connective gap), PREREG Sec. 5, with DEVIATIONS.md C1, C6 and C7.

Criteria (frozen): on the weight >= 5 giant SCC (157,821 neurons)
  (i)   |lambda_2| lies in [0.90, 0.97];
  (ii)  removing ascending and descending neurons gives 1 - |lambda_2| < 0.002;
  (iii) removing all motor neurons moves |lambda_2| by less than 0.01.
Removal sets are those of slowmode_check.py (C1). For (ii): if the giant SCC of what remains still
holds at least half of the remaining brain-side and half of the remaining cord-side neurons, (ii) is
scored on its |lambda_2|, as in slowmode_check.py; otherwise brain and cord no longer form one
communicating class and (ii) is scored as passed, labelled "decoupled" (C1).
Pass iff (i), (ii) and (iii). The scored result is saved before any sensitivity row is computed.
"""
import sys
import numpy as np
from common import (log, load_edges, giant_scc, row_stochastic, second_modulus, superclass_of, side_of,
                    save_result, expect, NotDetermined, ARPACK_SEED_CHECK, NECK_SCRIPT, NECK_TEXT,
                    NECK_ALL, MOTOR_SCRIPT, MOTOR_ALL, VNC_SCRIPT)

pre, post, w, body_ids = load_edges()
n = len(body_ids)
sc_all = superclass_of(body_ids)
side_all = side_of(sc_all)

keep5, A5 = giant_scc(pre, post, w, n, thresh=5)
log(f'weight >= 5 giant SCC: {len(keep5):,} neurons, {A5.nnz:,} edges')
expect(n, 'G5', len(keep5))
in5 = np.zeros(n, bool)
in5[keep5] = True


def composition(idx):
    s = side_all[idx]
    return {k: int((s == k).sum()) for k in ('brain', 'cord', 'neck')}


def lam2_row(label, removed_classes, check_decoupling=False):
    mask = in5.copy()
    if removed_classes:
        present = [c for c in removed_classes if (in5 & (sc_all == c)).any()]
        mask &= ~np.isin(sc_all, removed_classes)
    else:
        present = []
    removed = int(in5.sum() - mask.sum())
    keep, A, ridx, lab = giant_scc(pre, post, w, n, thresh=5, subset=mask, return_labels=True)
    rest, giant = composition(ridx), composition(keep)
    two_sided = (giant['brain'] >= 0.5 * rest['brain']) and (giant['cord'] >= 0.5 * rest['cord'])
    log(f'{label}: removed {removed:,} (classes present: {present}); remaining {rest}; giant SCC '
        f'{len(keep):,} = {giant}; two-sided: {two_sided}')
    row = {'label': label, 'removed_classes': removed_classes, 'classes_present': present,
           'removed_neurons': removed, 'remaining_composition': rest, 'scc_size': int(len(keep)),
           'scc_composition': giant, 'scc_edges': int(A.nnz), 'two_sided': bool(two_sided)}
    try:
        P = row_stochastic(A)
        l2, vals, info = second_modulus(P, k=8, ncv=48, tol=1e-8, maxiter=6000)
        row.update(lam2=l2, gap=1 - l2, lam1=complex(vals[0]), top8=[complex(v) for v in vals], arpack=info)
        log(f'    |lam1| = {abs(vals[0]):.6f}   |lam2| = {l2:.6f}   gap {1 - l2:.6f}')
    except NotDetermined as e:
        row.update(lam2=None, gap=None, not_determined=str(e))
        log(f'    NOT DETERMINED: {e}')
    if check_decoupling and not two_sided:
        row['decoupled'] = True
        row['sides'] = {}
        for sd in ('brain', 'cord'):
            labs, cnt = np.unique(lab[side_all[ridx] == sd], return_counts=True)
            if len(labs) == 0:
                continue
            sub = np.zeros(n, bool)
            sub[ridx[lab == labs[cnt.argmax()]]] = True
            k2, A2 = giant_scc(pre, post, w, n, thresh=5, subset=sub)
            try:
                l2s, _, _ = second_modulus(row_stochastic(A2), k=8, ncv=48, tol=1e-8, maxiter=6000)
            except NotDetermined:
                l2s = None
            row['sides'][sd] = {'scc_size': int(len(k2)), 'lam2': l2s}
            log(f'    decoupled; largest {sd}-side SCC {len(k2):,}, |lam2| = {l2s}')
    return row


rows ={'full': lam2_row('weight >= 5 giant SCC (reference)', []),
        'neck_script': lam2_row('ascending + descending removed (slowmode_check set)', NECK_SCRIPT,
                                check_decoupling=True),
        'motor_script': lam2_row('motor neurons removed (vnc_motor)', MOTOR_SCRIPT)}
for key in ('neck_script', 'motor_script'):
    if not rows[key]['classes_present']:
        raise SystemExit(f'STRUCTURAL CHECK FAILED: no neuron of {rows[key]["removed_classes"]} in the SCC')

l2 = rows['full']['lam2']
nd = [k for k, r in rows.items() if r.get('lam2') is None and not r.get('decoupled')]
if nd:
    verdict, c1, c2, c3 = 'NOT DETERMINED', None, None, None
else:
    c1 = 0.90 <= l2 <= 0.97
    c2 = True if rows['neck_script'].get('decoupled') else rows['neck_script']['gap'] < 0.002
    c3 = abs(rows['motor_script']['lam2'] - l2) < 0.01
    verdict = 'PASS' if (c1 and c2 and c3) else 'FAIL'
    log(f'\n(i)   |lam2| = {l2:.6f} in [0.90, 0.97]: {c1}')
    if rows['neck_script'].get('decoupled'):
        log('(ii)  brain and cord decoupled after neck removal: scored as passed (C1)')
    else:
        log(f'(ii)  1 - |lam2| after neck removal = {rows["neck_script"]["gap"]:.6f} < 0.002: {c2}')
    log(f'(iii) shift after motor removal = {abs(rows["motor_script"]["lam2"] - l2):.6f} < 0.01: {c3}')
save_result('h1', {'hypothesis': 'H1', 'criteria': {'i': c1, 'ii': c2, 'iii': c3}, 'verdict': verdict,
                   'rows': rows})
log(f'H1 {verdict}')

# ---------------------------------------------------------------- unscored rows
log('\n--- sensitivity rows (not scored) ---')
sens = {'neck_text': lam2_row('ascending_neuron + descending_neuron only (PREREG 4.2 text)', NECK_TEXT,
                              check_decoupling=True),
        'neck_all': lam2_row('every neck-crossing superclass removed', NECK_ALL, check_decoupling=True),
        'motor_all': lam2_row('vnc_motor and cb_motor removed', MOTOR_ALL),
        'brain_only': lam2_row('brain only (VNC superclasses removed, as script)', VNC_SCRIPT)}
checks = {}
for key, classes in (('full', []), ('neck_script', NECK_SCRIPT), ('motor_script', MOTOR_SCRIPT)):
    mask = in5 & ~np.isin(sc_all, classes) if classes else in5
    _, A = giant_scc(pre, post, w, n, thresh=5, subset=mask)
    try:
        l2b, _, info = second_modulus(row_stochastic(A), k=8, ncv=64, tol=1e-8, maxiter=6000,
                                      seed=ARPACK_SEED_CHECK)
    except NotDetermined:
        l2b, info = None, None
    checks[key] = {'lam2_ncv64': l2b, 'arpack': info}
    log(f'cross-check {key}: |lam2| with ncv = 64, seed {ARPACK_SEED_CHECK}: {l2b}')
save_result('h1_sensitivity', {'hypothesis': 'H1', 'scored': False, 'rows': sens, 'arpack_checks': checks})
sys.exit(0)
