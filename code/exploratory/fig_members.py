"""Figure 2 of the paper (exploratory): the certificate of H4 against exact chains built from the
confidence data. For each chain the witness lower bound max_{a,b} d_r(a, b) over the 4,950 pairs of
the H4 witness rows; for the confidence class at c = 0.51, 0.60, 0.70 the certified lower bound of
Lemma 2 (results/h4_full.json), drawn at zero where no certificate is obtained.
Reads exploratory/results/e2_class_members.json (or the .jsonl while E2 runs); writes
paper/fig_members.pdf and .png. ASCII, no em dashes."""
import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
BLUE, BLUE_FILL, ORANGE = '#2a78d6', '#2a78d6', '#eb6834'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'

h4 = json.load(open(os.path.join(ROOT, 'results', 'h4_full.json')))
fj = os.path.join(ROOT, 'exploratory', 'results', 'e2_class_members.json')
if os.path.exists(fj):
    members = json.load(open(fj))['members']
else:
    members = [json.loads(l) for l in open(fj + 'l')]
P = [m for m in members if m['label'].startswith('P ')][0]
others = [m for m in members if not m['label'].startswith('P ')]
R = sorted(int(k) for k in P['max_d'])
p_line = [P['max_d'][str(r)] for r in R]
lo = [min(m['max_d'][str(r)] for m in others) for r in R]
hi = [max(m['max_d'][str(r)] for m in others) for r in R]

cert = {}
for c in ('0.51', '0.60', '0.70'):
    rr = [x['r'] for x in h4['rows'] if x['r'] <= R[-1]]
    cert[c] = (rr, [max(x['robust_lower'][c], 0.0) for x in h4['rows'] if x['r'] <= R[-1]])

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'font.family': 'serif'})
fig, ax = plt.subplots(figsize=(5.6, 3.9))
ax.fill_between(R, lo, hi, color=BLUE_FILL, alpha=0.18, linewidth=0,
                label=f'{len(others)} chains built from the confidence data (range)')
ax.plot(R, p_line, color=BLUE, lw=1.6, marker='o', ms=4, label=r'$P$ (all synapses)')
styles = {'0.51': ('-', '^'), '0.60': ('--', 's'), '0.70': (':', 'D')}
for c, (rr, v) in cert.items():
    ls, mk = styles[c]
    ax.plot(rr, v, color=ORANGE, lw=1.4, ls=ls, marker=mk, ms=3.5,
            label=f'certified for the whole class, $c={c}$')
ax.set_xscale('log', base=2)
ax.set_xticks(R); ax.set_xticklabels([str(r) for r in R])
ax.set_ylim(-0.02, 1.03); ax.set_xlim(0.9, 72)
ax.set_xlabel('synaptic steps $r$')
ax.set_ylabel(r'lower bound on $\tau$ of the $r$-step chain')
ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
h, l = ax.get_legend_handles_labels()
order = [1, 0, 2, 3, 4]                                   # P first, then the band, then the certificates
fig.legend([h[i] for i in order], [l[i] for i in order], frameon=False, loc='lower center', ncol=2,
           fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=(0, 0.17, 1, 1))
out = os.path.join(ROOT, 'paper', 'fig_members.pdf')
fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
print('saved', out, '| members:', [m['label'].split(':')[0] for m in members])
