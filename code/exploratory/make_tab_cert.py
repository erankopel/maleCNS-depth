"""Writes paper/tab_cert_rows.tex (rows of the certificate table of Section 4.2) from
exploratory/results/e2b_sharp_certificate.json and e2_class_members.json. ASCII, no em dashes."""
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
B = json.load(open(os.path.join(ROOT, 'exploratory', 'results', 'e2b_sharp_certificate.json')))
E = json.load(open(os.path.join(ROOT, 'exploratory', 'results', 'e2_class_members.json')))
RS = [2, 4, 8, 16, 32]
rows = {x['r']: x for x in B['rows']}
exact = {m['label'].split(':')[0].strip(): m for m in E['members']}


def f(v):
    return f'{v:.3f}' if v >= 0 else f'$-{abs(v):.3f}$'


out = []
for c in ('0.55', '0.60', '0.70'):
    st = B['radius_stats'][c]
    for key, name in (('eps', r'$\varepsilon$ (Lemma~\ref{lem:row}, H4)'), ('eps_star', r'$\varepsilon^*$ (Lemma~\ref{lem:sharp})'),
                      ('delta_M', r'row change of $M(c)$, not certified')):
        lab = f'$c={c}$' if key == 'eps' else ''
        vals = ' & '.join(f(rows[r]['robust_lower'][f'{key}_{c}']) for r in RS)
        out.append(f"{lab} & {name} & {st[key]['mean']:.3f} & {vals} \\\\")
    m = exact.get(f'M({c})')
    if m is not None:
        vals = ' & '.join(f"{m['max_d'][str(r)]:.3f}" for r in RS)
        out.append(f"& exact witness bound of $M(c)$ & & {vals} \\\\")
    if c != '0.70':
        out.append(r'\addlinespace')
vals = ' & '.join(f"{rows[r]['nominal_lower']:.3f}" for r in RS)
out.append(r'\addlinespace')
out.append(f"nominal $P$ & & & {vals} \\\\")
open(os.path.join(ROOT, 'paper', 'tab_cert_rows.tex'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
