"""Writes paper/tab_members_rows.tex (rows of Table 4) from exploratory/results/e2_class_members.json and
e2c_trap_in_members.json, so that the table cannot drift from the results. ASCII, no em dashes."""
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
E = json.load(open(os.path.join(ROOT, 'exploratory', 'results', 'e2_class_members.json')))
fc = os.path.join(ROOT, 'exploratory', 'results', 'e2c_trap_in_members.json')
C = {r['member']: r for r in json.load(open(fc))}
lines = []
for m in E['members']:
    lab = m['label'].split(':')[0].strip()
    if lab.startswith('P '):
        lab = '$P$'
    else:
        lab = '$' + lab.replace('(', '(').replace(')', ')') + '$'
    dlt = '' if 'max_abs_change_d' not in m else f"{m['max_abs_change_d']['32']:.3f}"
    w5 = f"{m['lam2_w5']:.4f}"
    key = m['label'].split(':')[0].strip()
    key = 'P' if key.startswith('P ') else key
    trap = str(C[key]['trap_in_scc'])                   # neurons of the 11-neuron set in the w >= 5 SCC
    lines.append(f"{lab} & {m['synapses'] / 1e6:.2f} & {m['row_tv_to_P_mean']:.3f} & {m['lam2_w1']:.4f} & {w5} & "
                 f"{trap} & {m['max_d']['32']:.3f} & {dlt} \\\\")
    if lab == '$P$':
        lines.append(r'\addlinespace')
    if lab == '$M(0.90)$':
        lines.append(r'\addlinespace')
open(os.path.join(ROOT, 'paper', 'tab_members_rows.tex'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
