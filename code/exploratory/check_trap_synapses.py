"""Check (exploratory): synapses crossing the boundary of the 11-neuron set of E1, from the raw edge list
of the flat connectome restricted to G (weight >= 1 giant SCC), with no sparse-matrix index order
involved. Written because e1_slow_mode.py paired the weights of one CSR matrix with the column indices
of another (A1.data with P1.indices), which is valid only if both store each row in the same order.
ASCII, no em dashes."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'confirmatory'))
from common import DATA, log, load_edges, giant_scc
TRAP = [800719, 804848, 806145, 807560, 808551, 808845, 811415, 813264, 819906, 820114, 913939]
pre, post, w, body_ids = load_edges()
N = len(body_ids)
keep, A = giant_scc(pre, post, w, N, thresh=1)
ing = np.zeros(N, bool); ing[keep] = True
t = np.isin(body_ids, TRAP); assert t.sum() == 11 and ing[t].all()
m = ing[pre] & ing[post]
out_ = m & t[pre] & ~t[post]; in_ = m & ~t[pre] & t[post]; inside = m & t[pre] & t[post]
res = {'synapses_out': int(w[out_].sum()), 'synapses_out_on_edges_ge5': int(w[out_ & (w >= 5)].sum()),
       'synapses_in': int(w[in_].sum()), 'synapses_in_on_edges_ge5': int(w[in_ & (w >= 5)].sum()),
       'synapses_inside': int(w[inside].sum()), 'synapses_inside_on_edges_ge5': int(w[inside & (w >= 5)].sum()),
       'edges_out': int(out_.sum()), 'edges_out_ge5': int((out_ & (w >= 5)).sum())}
# the same through the CSR matrix of G, pairing A's own data and indices
A.sort_indices()
rows = np.repeat(np.arange(len(keep)), np.diff(A.indptr)); S = t[keep]
co = S[rows] & ~S[A.indices]; ci = ~S[rows] & S[A.indices]
res['csr_check'] = [int(A.data[co].sum()), int(A.data[co & (A.data >= 5)].sum()), int(A.data[ci].sum()),
                    int(A.data[ci & (A.data >= 5)].sum())]
log(json.dumps(res))
json.dump(res, open(os.path.join(HERE, '..', '..', 'exploratory', 'results', 'check_trap_synapses.json'), 'w'), indent=1)
