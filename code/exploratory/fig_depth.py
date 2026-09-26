"""Figure 1 of the paper: certified two-sided bounds on the Dobrushin coefficient tau(P^r) of the male CNS
synapse-flow chain. Reads results/h3_full.json, results/h4_full.json and
exploratory/results/e1b_without_trap.json; writes paper/fig_depth.pdf. ASCII, no em dashes."""
import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'

h3 = json.load(open(os.path.join(ROOT, 'results', 'h3_full.json')))
h4 = json.load(open(os.path.join(ROOT, 'results', 'h4_full.json')))
e1b = json.load(open(os.path.join(ROOT, 'exploratory', 'results', 'e1b_without_trap.json')))[0]

r_w1 = [x['r'] for x in h4['rows']]; lo_w1 = [max(x['nominal_lower'], 0) for x in h4['rows']]
r_w5 = [x['r'] for x in h3['rows']]; lo_w5 = [x['lower'] for x in h3['rows']]; up_w5 = [x['upper'] for x in h3['rows']]
r_t = sorted(int(k) for k in e1b['bounds']); lo_t = [e1b['bounds'][str(k)][0] for k in r_t]
up_t = [e1b['bounds'][str(k)][1] for k in r_t]

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'font.family': 'serif'})
fig, ax = plt.subplots(figsize=(5.6, 3.5))
ax.fill_between(r_w5, lo_w5, up_w5, color=ORANGE, alpha=0.12, linewidth=0)
ax.fill_between(r_t, lo_t, up_t, color=AQUA, alpha=0.12, linewidth=0)
ax.plot(r_w1, lo_w1, color=BLUE, lw=1.6, marker='o', ms=4, label=r'$w\geq1$ (lower bound)')
ax.plot(r_w5, lo_w5, color=ORANGE, lw=1.6, marker='s', ms=4, label=r'$w\geq5$ (lower and upper bound)')
ax.plot(r_w5, up_w5, color=ORANGE, lw=1.2, ls='--')
ax.plot(r_t, lo_t, color=AQUA, lw=1.6, marker='^', ms=4, label=r'$w\geq5$ without the 11-neuron set')
ax.plot(r_t, up_t, color=AQUA, lw=1.2, ls='--')
ax.set_xscale('log', base=2)
ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 128]); ax.set_xticklabels(['1', '2', '4', '8', '16', '32', '64', '128'])
ax.set_ylim(0, 1.03); ax.set_xlim(0.9, 150)
ax.set_xlabel('synaptic steps $r$')
ax.set_ylabel(r'certified bounds on $\tau(P^r)$')
ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.annotate('0.669', (128, lo_w5[-1]), xytext=(4, -2), textcoords='offset points', color=INK2, fontsize=8)
ax.annotate('0.107', (128, lo_t[-1]), xytext=(4, -2), textcoords='offset points', color=INK2, fontsize=8)
ax.annotate('0.001', (128, lo_w1[-1]), xytext=(4, 2), textcoords='offset points', color=INK2, fontsize=8)
ax.legend(frameon=False, loc='lower left', fontsize=8)
fig.tight_layout()
out = os.path.join(ROOT, 'paper', 'fig_depth.pdf')
fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
print('saved', out)
