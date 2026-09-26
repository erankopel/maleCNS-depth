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
