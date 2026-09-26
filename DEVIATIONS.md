# Deviations, corrections and clarifications

`PREREG_maleCNS_depth_2026-09-19.md` is frozen and is never edited. Every departure from it, every
error found in it, and every definition it leaves open is recorded here, with the date the entry
was written and the state of the analysis on that date. ASCII, no em dashes.

## Entries written 26 Sep 2026, before registration and before any H1 to H5 computation

D1. Repository. PREREG Sec. 6 says the scripts go to "the GitHub repo under a new fly/ directory".
They are published instead in this dedicated repository, github.com/erankopel/maleCNS-depth, with
the same file layout as the 19 Sep bundle. No analysis is affected.

D2. Registration date. The PREREG text was drafted and frozen in a private chat session on
19 Sep 2026 and first made public by this repository's first commit and its Zenodo record on
26 Sep 2026. The registration date to cite is 26 Sep 2026. Between the two dates nothing from
H1 to H5 was computed (REGISTRATION.md).

D3. Compute location. PREREG Sec. 5 says each hypothesis is run "on the PC". The runs are made
instead on a cloud Linux container (2 vCPU Intel Xeon 2.1 GHz, 7.8 GB RAM; Python 3.11.15,
numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, pyarrow 25.0.1), because the sandbox through which the
author's PC is reached limits each command to 3 minutes and does not keep background jobs. Inputs
are md5-checked against PREREG Sec. 1 before every run, and every run's log and environment are
committed. The 26 Sep reproduction of Sec. 4 on this container matches the 19 Sep numbers
(REPRODUCTION_2026-09-26.md).

D4. Error in code/cert.py: the robust lower bound omits the t = 0 term. PREREG Sec. 4.4 states
|d_r(P') - d_r(P)| <= eps sum_{t<r} 2 d_t(P), and the sum includes t = 0, where d_0(a, b) = 1 for
any witness pair a != b: the perturbation already acts on the first step, from x_0 = e_a - e_b.
cert.py accumulates only t = 1, ..., r-1, so every robust lower bound it prints is too high by
2 eps = 0.2105 (delta = 0.05, eps = 2 delta / (1 - delta)). Values from the 26 Sep reproduction
(unchanged cert.py) and corrected, max(0, printed - 2 eps):

    r     lower bound   robust lower, printed   robust lower, corrected
    1     1.0000        1.0000                  0.789
    2     1.0000        0.7895                  0.579
    4     0.9896        0.3583                  0.148
    8     0.8794        0.0000                  0.000
    16    0.6231        0.0000                  0.000
    32    0.2663        0.0000                  0.000
    64    0.0428        0.0000                  0.000

Not affected: the non-robust lower bounds, the upper bounds and the robust upper bounds. The
conclusion carried into H4 (a uniform per-edge perturbation class cannot be certified through
beyond a few hops) is unchanged and becomes stronger. The H4 script implements the bound exactly
as stated in Sec. 4.4, including t = 0.

D5. Error in PREREG Sec. 4.2: the knock-out count. Sec. 4.2 reports "remove ascending +
descending neurons (3,160), lambda_2 = 0.9997". 3,160 is the size of the ascending_neuron (1,846)
plus descending_neuron (1,314) group of modes.py. slowmode_check.py, which produced
lambda_2 = 0.9997, removes in addition the superclass sensory_ascending (533 neurons in G): 3,693
neurons were removed, and the giant SCC of the rest has 161,597 neurons. Five smaller
neck-crossing superclasses (descending_neuron_tbc 2, sensory_ascending_tbc 2, efferent_ascending
8, efferent_descending 4, sensory_descending 12; 28 neurons in G) were not removed, so brain and
cord stay weakly coupled in that knock-out. How H1 treats this is fixed in an entry written
before the confirmatory runs.

D6. Prior work (PREREG Sec. 3), checked against full texts on 26 Sep 2026. These corrections
change the positioning of the paper, not any hypothesis.

(a) Brain-plus-cord connectomes. Sec. 3 calls a spectral-gap or mixing-depth analysis of a
brain-plus-cord chain "impossible before this dataset: FlyWire has no VNC and MANC has no brain".
That is wrong. BANC, a female brain-and-nerve-cord connectome, preceded this dataset (preprint
August 2025, doi:10.1101/2025.07.31.667571; Bates et al., "Distributed control circuits across a
brain-and-cord connectome", Nature, published 8 Jun 2026, doi:10.1038/s41586-026-10735-w). BANC
analyses the whole graph with an unsigned, postsynaptically normalised "influence" model whose
matrix is rescaled so that "its largest real eigenvalue is 0.99", with Laplacian spectral
clustering and with betweenness centrality. As far as we found, neither BANC nor any other study
reports the spectral gap, the mixing depth or the stationary distribution of a synapse-weighted
walk on either brain-plus-cord connectome, and no study reports certified bounds for any
connectome-derived operator.

(b) Lin et al., Nature 634, 153-165 (2024), doi:10.1038/s41586-024-07968-y. Their random walk is
built on the unweighted "0-1 adjacency matrix" (transition probability 1/out-degree), not on
synapse counts, and they show the eigenvalue spectra of the normalised Laplacians of the forward
and reverse walks ("The gaps between eigenvalues indicate the conductance properties of the
graph."). A spectral look at a connectome random walk is therefore not new. What remains new is
the synapse-weighted chain on a whole CNS, the slow mode tied to the neck connective, and
certified two-sided bounds on the mixing depth.

(c) Male CNS paper: Berg et al., "Sexual dimorphism in the complete Drosophila male central
nervous system connectome", Cell 189, 5504-5526.e15 (2026), doi:10.1016/j.cell.2026.08.015
(preprint doi:10.1101/2025.10.09.680999). The full preprint contains no eigenvalue, spectral-gap,
random-walk or mixing analysis. Its depth notions are maximum flow with input-normalised
capacities (Dinic's method), layers from the pseudo-inverse of the Laplacian of the flow graph,
and the probabilistic traversal layers of Schlegel et al., eLife 10, e66018 (2021). Its wording
on feedback is that neck-connective neurons "participate in the majority of feedback
connections" (superclass level) and that "feedback connections have a larger inhibitory
proportion". The published Cell text could not be retrieved on 26 Sep 2026; the check rests on
the preprint and the Cell metadata.

## Entries written 26 Sep 2026, after registration (commit cd16354, Zenodo doi:10.5281/zenodo.22980248) and before any H1 to H5 computation

These entries fix every definition that the PREREG leaves open. They were written together with the
confirmatory scripts in `code/confirmatory/`, which were run only on synthetic data before the
freeze: `tests/test_bounds.py` checks every bound against brute force on sampled chains, and
`tests/test_confirmatory.py` checks every scored number, verdict and descriptive quantity of h1 to h5
against independent dense computations on four synthetic data sets (default, high confidence,
brain and cord decoupled, one synapse row missing). An independent review of the scripts was made
before the freeze and its findings are addressed below. The release made from the commit that adds
these entries freezes them.

D7. Scripts and real-data reads after registration. PREREG Sec. 5 says each hypothesis is run "with
the scripts in Section 6". H1 to H4 need variants of those scripts (the 19 Sep README says so), and
they are the scripts in `code/confirmatory/`; they reuse the Sec. 6 numerics except where an entry
here says otherwise. The only computations on the real data between registration and the freeze
were: the input-weight statistics of G used in C2 (median signed-eligible input 342, mean 730;
median total input 349, mean 750); the content hash of the rebuilt graph (C6); and, during the
review, metadata only: the syn-partners schema, the dictionary and null counts of body_pre,
body_post and primary_post (no confidence column decoded), and the graph's self-loop count (101)
and duplicate-edge count (0). No eigenvalue, bound or confidence statistic was computed.

C1. H1.
- Graph: edges with w >= 5 among the neurons N; its giant SCC (157,821 neurons, checked); P
  row-normalised.
- lambda_2 is the eigenvalue of second-largest modulus of P, from ARPACK with which = 'LM', k = 8,
  ncv = 48, tol = 1e-8, maxiter = 6000 (the setting of spectra_scc.py), seeded as in C7. Criterion
  (i) is on |lambda_2|.
- Removal sets are those of slowmode_check.py: for (ii) the superclasses ascending_neuron,
  descending_neuron and sensory_ascending; for (iii) the superclass vnc_motor (the "674 motor
  neurons" of Sec. 4.2). Neurons are removed from the w >= 5 giant SCC; each removal set must
  remove at least one neuron.
- Brain side and cord side: the cord side is the VNC superclass list of slowmode_check.py; the neck
  is the eight neck-crossing superclasses listed below; everything else is the brain side. Before
  any eigenvalue, h1 prints how many brain-side, cord-side and neck neurons remain and how many lie
  in the giant SCC of what remains.
- Scoring of (ii). If the giant SCC of what remains holds at least half of the remaining
  brain-side neurons and at least half of the remaining cord-side neurons, (ii) is scored on its
  |lambda_2|, exactly as slowmode_check.py computes it. Otherwise brain and cord no longer form one
  communicating class: that is the limiting case of what (ii) predicts (a chain that splits into
  blocks has lambda_2 = 1), so (ii) is scored as passed and labelled "decoupled". The |lambda_2| of
  the largest brain-side and cord-side SCCs is reported alongside.
- Why this removal set for (ii): it is the one that produced the exploratory 0.9997 on which the
  0.002 threshold was set (D5). With the 28 remaining neck-crossing neurons in place, (ii) measures
  how weakly brain and cord are coupled through them.
- Reported, not scored, after the scored result is saved: ascending_neuron + descending_neuron only
  (the 3,160 of the Sec. 4.2 text); every neck-crossing superclass removed (the three above plus
  sensory_ascending_tbc, descending_neuron_tbc, efferent_ascending, efferent_descending and
  sensory_descending); vnc_motor and cb_motor removed; brain only; and a cross-check of the three
  scored rows with ncv = 64 (outside the Sec. 6 range 40 to 60, hence unscored).

C2. H2.
- G is the w >= 1 giant SCC (165,314 neurons, checked), with the signs of Sec. 2.
- in_j is the total weight of edges into j from presynaptic neurons with nonzero sign: the
  denominator of the unfloored map of spectra.py and modes.py.
- F is the median of in_j over all neurons of G: 342 synapses (checked). The PREREG adds "(about 700
  synapses)". That figure is wrong for the stated statistic: the median is 342 (signed-eligible
  input) or 349 (total input), while 700 is close to the mean (730 signed-eligible, 750 total) and to
  the median of the cb_intrinsic superclass alone (698, printed in the Sec. 4 logs). The defining
  word, median, governs.
- Denominator max(in_j, F) for every neuron j of G.
- ARPACK (which = 'LM', ncv = 60, tol = 1e-8, maxiter = 10000, seeded as in C7) computes 24 right
  eigenvectors; the modes are ordered by decreasing modulus and then by decreasing imaginary part,
  and the scored modes are the first 20 (computing 24 keeps modes near the cut from being missed; a
  conjugate pair split by the cut contributes its first member).
- (a) The spectral radius is the modulus of the first mode; (a) passes if it is below 0.85.
  (b) The participation ratio of a mode is 1 / sum_i p_i^2 with p_i = |v_i|^2 / sum_k |v_k|^2;
  (b) passes if at least 10 of the 20 modes have PR > 20.
- (c), descriptive: a neuron counts if its primary neuropil is LOP(L), LOP(R) or EB; a mode's power
  on the three neuropils is counted together (the single-neuropil shares are reported too). A
  neuron's primary neuropil is the category of the file's primary_post vocabulary that holds most
  of the syn-partners rows in which the neuron is the presynaptic or the postsynaptic partner;
  ties go to the earlier category, and the "<unspecified>" category competes like any other
  (prep_neuropil.py; no confidence column is decoded).
- Reported, not scored, after the scored result is saved: the same computation with F = 700 and
  with F = median total input.

C3. H3.
- Chain as in C1. Stationary distribution by power iteration from the uniform vector, as in
  spectra_scc.py, stopping at an L1 change below 1e-13 (at most 10,000 iterations; if the tolerance
  is not met the last vector is used and reported: it only chooses witness rows, and every pair
  gives a valid bound).
- Witness rows: numpy default_rng(7), 90 distinct rows from rng.choice(n, 90, replace=False), then
  the 10 rows of largest stationary mass (the recipe of cert.py with 90 + 10). A top row may repeat
  a random row; the number of distinct rows is reported.
- Lower bound at each checkpoint: the maximum over witness pairs of d_r, minus 1e-9 r. Upper bound,
  reported only: the Markov bound over the 500 columns of largest stationary mass (Sec. 7), plus
  1e-9 r.
- Criterion 2: the first checkpoint of Sec. 6 (1, 2, 4, 8, 16, 32, 48, 64, 96, 128) at which the
  lower bound is below 0.05 must be 48, 64 or 96. If the bound is still at least 0.05 at 128,
  criterion 2 fails. The verdict is saved as soon as both criteria are determined; the loop then
  continues to r = 128 for the descriptive rows.
- Rounding: D is the larger of the maximum in-degree and the maximum out-degree of P (cert.py used
  the out-degree only; the row and column propagations need both). After r hops each propagated
  row has an L1 rounding error of at most about r (D + 1) u (below 2e-12 r), and the final sums over
  n terms add at most about 2 n u (below 4e-11), so the total stays far inside the 1e-9 r slack for
  every r >= 1.

C4. H4.
- Graph and chain: G and P of Sec. 2, the chain of cert.py, whose method H4 invokes (H4 itself
  names no threshold). Witness rows: the H3 recipe (90 + 10, seed 7) applied to G; checkpoints to
  128 (the Sec. 7 sizes).
- Nominal weights: the minconf-0.5 flat weights. Pass 1 decodes only the body_pre and body_post
  columns and checks that for every edge of G the weight equals the number of syn-partners rows with
  that (body_pre, body_post). Any mismatch stops the script before a confidence value is read; the
  handling is then decided in a new, dated entry here before the script is run again.
- Pass 2 decodes the confidence columns and stops if any value is missing, not finite, or outside
  [0, 1].
- Class at level c: w_hi(c) counts the rows of an edge with min(conf_pre, conf_post) >= c, compared
  in float32, the stored type. Every P' obtained by row-normalising some w' with
  w_hi(c) <= w' <= w entrywise is in the class; a row with no remaining weight may take any
  distribution.
- Per-row radius: eps_i(c) = 1 - W_i(c) / W_i, with W the row sums. Removing mass R from a row of
  mass W moves the normalised row by at most R / W in total variation, so eps_i(c) bounds
  TV(P_i, P'_i) for every P' in the class (eps_i = 1 for a row whose weight can vanish).
- Robust lower bound: the Sec. 4.4 argument with per-row eps and the t = 0 term (D4):
  LB_r = max over witness pairs of [ d_r(a, b) - sum_{t=0}^{r-1} sum_i |x_t(i)| eps_i(c) ] - 1e-9 r.
  The maximum is taken over each pair's robust value. cert.py used the pair with the largest
  nominal d_r; every pair gives a valid bound, so maximising the robust value is also valid and
  never smaller. The rounding error of the eps sums is below 5e-11 r; with the propagation error the
  total stays below 1e-10 r, inside the 1e-9 r slack.
- Scored: c = 0.70 and r = 32; H4 passes if LB > 0.10. The verdict is saved as soon as the r = 32
  bound exists; the loop then continues to r = 128 for the descriptive rows.
- Fail clause ("reports at which confidence level it is"): the largest c on the grid 0.51, 0.52,
  0.53, 0.54, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 with LB at r = 32 above 0.10, or
  "none".
- Reported for comparison: the corrected uniform-class bound C(0.05) with the same 100 witnesses.

C5. H5.
- Pospisil convention: acetylcholine and dopamine +1; GABA, glutamate, histamine, serotonin and
  octopamine -1; tyramine, unclear and unknown 0 (their output edges dropped), as in Sec. 2.
- Map, normalisation and eigen-solver as in modes.py item 1 (ARPACK, k = 10, ncv = 60, tol = 1e-8,
  seeded as in C7) on G (165,314 neurons, checked), under both conventions; modes are right
  eigenvectors ordered as in C2.
- "The spectral radius changes by less than 0.02" is |rho_Pospisil - rho_base| < 0.02.
- "The top-six carriers are unchanged": for each of the six leading modes, the carrier set is the
  smallest set of neurons holding at least 90 percent of the mode's power. The union of the six sets,
  as body IDs, must be identical under the two conventions. (A conjugate pair has one power
  distribution, so the union does not depend on which member of a pair is listed first.)
- As the PREREG says, H5 is reported as a sensitivity result, not as a hypothesis test.

C6. Run protocol.
- The code that runs must match `code/confirmatory/MANIFEST.sha256` of the frozen release, and the
  package versions must match C8; otherwise the runner stops before any computation.
- Before any hypothesis, verify_inputs.py checks the md5 sums of Sec. 1, the syn-partners md5
  58efcf712f8c4d4de5f2ad51e97def76, and the content hash
  cd6bbde8abf155fa746aaa64afc34f5c5b2454ac7867426a2a4972fc2566b0e8 of the graph rebuilt by
  build_graph.py; each h-script checks the structural facts marked "checked" above. Any mismatch
  stops the run.
- run_confirmatory.sh executes prep_neuropil.py and h1 to h5 once, in that order. Every attempt
  writes its logs to a new directory; results are never overwritten. All logs, results and verdicts
  are committed, whatever they say.
- Each script saves and prints its verdict as soon as its scored numbers exist, before any unscored
  computation. A hypothesis counts as run once its scored result has been saved, and it is not run
  again. If a script stops after saving its scored result, the runner completes only the unscored
  part, in a mode that recomputes everything and stops unless the saved scored result is reproduced
  exactly, number for number. If a script stops before saving its scored result, the error, the fix
  and the rerun are recorded here with dates, and any scored number already printed in the failed
  attempt's log must be reproduced exactly by the rerun.
- Runs are bit-reproducible: ARPACK starting vectors are seeded (C7) and BLAS runs on one thread
  (C8).

C7. Eigen-solver. Every ARPACK call takes its starting vector from numpy default_rng(11) (default_rng(12)
for the unscored cross-checks of C1), passed as the rng argument of scipy.sparse.linalg.eigs. If
ARPACK raises a non-convergence error, the call is repeated once with ncv doubled (at most n - 1)
and maxiter = 20000. If that also fails, the quantity is "not determined": the criterion that needs
it is reported as not determined and the hypothesis as NOT DETERMINED, not as passed or failed.
Wall times of every call are reported.

C8. Environment: Python 3.11.15, numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, pyarrow 25.0.1
(`code/confirmatory/requirements-confirmatory.txt`), with OPENBLAS_NUM_THREADS = OMP_NUM_THREADS =
MKL_NUM_THREADS = 1. This is the environment on which the Sec. 4 reproduction matched; the runner
checks it before starting.
