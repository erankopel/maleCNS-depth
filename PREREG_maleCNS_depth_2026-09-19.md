# Pre-registration: certified information-breaking depth of the male Drosophila CNS connectome

Draft 1, prepared 19 Sep 2026 (Claude, from the 19 Sep 2026 chat session). ASCII, no em dashes.
Status: exploratory phase closed at Section 4; confirmatory hypotheses in Section 5 are frozen
before any of the PC-phase computations in Section 7 are run.

## 1. Data and provenance

Male CNS connectome v1.0 (FlyEM / Cambridge / MRC LMB / Google Research), flat-connectome
release, minimum synapse confidence 0.5, CC-BY. Public GCS bucket, no token needed.
gs://flyem-male-cns/v1.0/connectome-data/flat-connectome/

Files used (md5, size):
- body-annotations-male-cns-v1.0-minconf-0.5.feather   50a7718770c57220f160ba4f431ab89e   14,483,314 B
- body-neurotransmitters-male-cns-v1.0.feather          3d842b12fe5c49eefade528d7dd24a1f   42 MB
- connectome-weights-male-cns-v1.0-minconf-0.5.feather  f30e9dcca25cfd021bf1e7b3d975599e   1,051,241,946 B

Deferred to the PC: syn-partners (6.8 GB, per-synapse conf_pre / conf_post), syn-points (12.7 GB),
skeletons, neo4j dump. neuPrint API (male-cns:v1.0) requires an account token; not used.

## 2. Objects (fixed)

- Neuron set N: bodies with non-null superclass in the annotation table. |N| = 166,700
  (211,577 annotated bodies; the rest are glia, orphans, "unimportant" fragments, unlabelled).
- Edge table: 151,856,684 segment-to-segment rows carrying 311,833,243 synapses; restricted to
  N x N: 25,582,938 directed edges, 124,177,617 synapses. At weight >= 5: 6,242,118 edges,
  89,860,280 synapses.
- Giant strongly connected component G (weight >= 1): 165,314 neurons (99.2 percent),
  25,545,360 internal edges. At weight >= 5 the giant SCC has 157,821 neurons (94.7 percent).
- Synapse-flow chain P on G: P_ij = w_ij / sum_k w_ik (presynaptic normalisation; row-stochastic).
- Signed map M on G: M_ji = s_i w_ij / sum_k w_kj (postsynaptic normalisation), s_i = +1 for
  acetylcholine, -1 for GABA, glutamate, histamine, 0 (edges dropped) for dopamine, serotonin,
  octopamine, tyramine, unclear, unknown. Consensus neurotransmitter covers 162,982 of 166,700
  neurons (103,720 +, 59,262 -, 3,718 excluded), median prediction confidence 0.936.
  Sensitivity convention (Pospisil et al. 2024): dopamine +, serotonin and octopamine -.
- Dobrushin coefficient tau(Q) = (1/2) max_{a,b} ||Q_a - Q_b||_1. tau(P^r) = 0 iff P^r has rank
  one; this is the classical analogue of the entanglement-breaking index (the composite map
  forgets its input). Witness distance d_r(a,b) = (1/2) ||(e_a - e_b) P^r||_1 <= tau(P^r).
- Perturbation class C(delta): every edge weight w -> w (1 + eta), |eta| <= delta, independently.

## 3. Prior work this must be positioned against (checked 19 Sep 2026)

- Lin et al., Nature 634 (2024), FlyWire female brain (no VNC): degree and motif statistics,
  giant SCC 93.3 percent at threshold 5, shortest paths, rich club (30 percent of neurons),
  forward and reverse random-walk stationary distributions ("attractor" and "repeller" neurons).
  Our stationary distribution and SCC numbers are the male-CNS counterpart and are not new.
- Pospisil et al., Nature 634, 201 (2024): eigendecomposition of the signed synapse-count
  matrix W (threshold 5, scaled for stability, no input normalisation). Eigenvalues decay slowly
  (1000th about 1/10 of the first); eigenvectors extremely sparse (top 10 modes: about 50 neurons
  for 75 percent of power); first mode localised to the lobula plate (VCH, DCH, LPi15, Am1);
  eigenvalue angles cluster at 0, 90, 180 degrees; numerical sensitivity analysis to biological
  variability, measurement error and saturation (Extended Data Fig. 9). Our localisation result
  (Section 4.2) replicates this in the male CNS under a different normalisation and is not new.
- Male CNS paper, Cell (Sept 2026): superclass connectivity is principally feedforward; neck
  connective (ascending and descending) neurons account for most feedback; end-to-end
  maximum-flow analysis; DN/AN clustering by sensorimotor flow. No spectral or mixing analysis
  reported (abstract and snippets checked; full text to be read on the PC).
- Not found: any certified or interval-arithmetic treatment of connectome-derived operators;
  any spectral-gap or mixing-depth analysis of a brain-plus-cord chain (impossible before this
  dataset: FlyWire has no VNC and MANC has no brain).

What would be new: (i) the spectral gap of the whole-CNS chain identified with neck-connective
exchange, as a number; (ii) certified two-sided bounds on tau(P^r) with a stated perturbation
class; (iii) the finding that a uniform per-edge perturbation model destroys certification within
eight hops, so robustness must come from per-synapse confidences (Section 4.4 and H4).

## 4. Exploratory results already obtained (19 Sep 2026; not confirmatory)

4.1 Unsigned chain on G. Second eigenvalue lambda_2 = 0.94503 (ARPACK, tol 1e-8), spectral gap
0.0550. Empirical TV to stationary from 20 random sources: 0.77 (r=4), 0.57 (8), 0.33 (16),
0.11 (32), 0.018 (64), 0.0045 (100). Stationary mass: top 1 percent of neurons hold 0.422, top
10 percent hold 0.854; entropy 13.96 bits of 17.33. VNC holds 0.653 of stationary mass with
20,109 neurons; 674 motor neurons alone hold 0.168 (median 2,729 input, 15 output synapses).
Whole-graph (not SCC-restricted) run gave a spurious 0.9986 from a small closed component.

4.2 Slow mode anatomy. Left eigenvector for lambda_2: 0.921 of positive mass in VNC, 0.938 of
negative mass in brain; top carriers leg motor neurons (+) and APL (-). Knock-outs: remove all
674 motor neurons, lambda_2 = 0.9430; remove ascending + descending neurons (3,160),
lambda_2 = 0.9997 (gap 0.0003); brain only, lambda_2 = 0.9256. Naive two-block Markov
prediction 1 - p - q = 0.78 fails (hub-weighted escape probabilities 0.144 and 0.071).

4.3 Signed map on G. Spectral radius 0.9151. Leading eigenvalues in near-equal +/- pairs
(-0.9151, +0.9137, -0.9122, +0.9098). Participation ratio of the top six modes 2.0 to 2.4
neurons, > 0.97 of mode power on 10 neurons; carriers R7/R8 photoreceptor pairs (histamine,
reciprocal, same sign) with left/right duplicates. Modes 7 and 9: PR 76 and 97 (R7d/R8d family;
ER3 ring neurons, GABA). Random unit impulse: ||M^r e||_1 = 0.86 (r=1), 0.18 (4), 0.016 (8),
2.4e-4 (16). Interpretation: the spectral radius is a two-neuron motif inflated by postsynaptic
normalisation of low-input sensory cells, not a global property.

4.4 Certified bounds on tau(P^r) (40 witness rows: 30 random + 10 top-stationary; 60 hub
columns; slack 1e-9 r for rounding, max degree 11,201):
  r     lower (witness)   upper (hubs)
  4     0.9896            1.0000
  8     0.8794            0.9997
  16    0.6231            0.9936
  32    0.2663            0.9697
  64    0.0428            0.9434
Robust versions under C(0.05) with the Lipschitz argument |d_r(P') - d_r(P)| <= eps sum_{t<r}
2 d_t(P), eps = 2 delta/(1 - delta): lower bound is 0.79 at r=2, 0.36 at r=4, zero from r=8.
Upper bound under C(0.05): no better than 0.9987 at any r. Conclusion carried into H4: uniform
per-edge uncertainty cannot be certified through; per-synapse confidences are required.

## 5. Confirmatory hypotheses (frozen 19 Sep 2026)

Each is run once, on the PC, with the scripts in Section 6, seeds as fixed there, and reported
whatever the outcome. Pass and fail criteria are stated before running.

H1 (neck-connective gap). On the weight >= 5 giant SCC (157,821 neurons), lambda_2 lies in
[0.90, 0.97]; removing ascending and descending neurons gives 1 - lambda_2 < 0.002; removing
all motor neurons moves lambda_2 by less than 0.01. Fail: any of the three false.

H2 (localisation is a normalisation effect). With a per-neuron input floor F = median input
weight of G (about 700 synapses) in the postsynaptic normalisation, i.e. M_ji = s_i w_ij /
max(in_j, F): (a) the spectral radius falls below 0.85; (b) the top 20 modes have participation
ratio > 20 neurons in at least 10 cases; (c) at least one of the top 20 modes has > 50 percent
of its power in lobula plate or ellipsoid body neurons (Pospisil's first and 45th eigencircuits).
Fail: (a) or (b) false; (c) is descriptive and reported either way.

H3 (certified depth). On the weight >= 5 giant SCC, with 100 witness rows (90 random + 10
top-stationary, seed 7): the certified lower bound on tau(P^32) exceeds 0.20, and the first
checkpoint at which the lower bound drops below 0.05 lies in [48, 96] hops. Fail: either false.

H4 (robustness needs confidences). Build per-edge intervals from syn-partners: an edge's weight
interval is [w_hi-only, w] where w_hi-only counts synapses with min(conf_pre, conf_post) >= 0.7.
Propagate as in cert.py with per-row eps_i computed from the row's low-confidence fraction.
Prediction: the certified lower bound on tau(P^32) remains above 0.10 under this class.
Fail: it does not; then the paper states that the depth result is not certifiably robust to
reconstruction uncertainty at the 0.7 level, and reports at which confidence level it is.

H5 (sign convention). Repeating 4.3 with Pospisil's convention changes the spectral radius by
less than 0.02 and leaves the top-six carriers unchanged. Fail: either false; reported as a
sensitivity result, not as a hypothesis test in the paper.

## 6. Fixed analysis choices

Thresholds: weight >= 1 and weight >= 5, both reported. Neuron set and SCC as in Section 2.
Seeds: numpy default_rng(0) for reachability and mixing sources, (1) for signed impulses,
(7) for witness rows. ARPACK: which='LM', tol 1e-7 or 1e-8, ncv 40 to 60. Checkpoints for
tau bounds: 1, 2, 4, 8, 16, 32, 48, 64, 96, 128. Rounding slack 1e-9 per hop (bound: relative
error per hop <= D 2^-53 < 1.3e-12 for D = 11,201). Scripts: build_graph.py, analyse.py,
spectra.py, spectra_scc.py, modes.py, slowmode_check.py, cert.py (this session's outputs;
to be committed to the GitHub repo under a new fly/ directory before H1 to H5 are run).

## 7. Deferred to the PC (heavy)

- syn-partners table: per-edge confidence intervals; per-neuron proofreading status
  (statusLabel) as a second perturbation class (Reviewed vs Roughly traced vs Prelim).
- Witness and hub set sizes 100 and 500; checkpoints to 128 hops.
- Full leading eigenvector sets (top 100) for both maps; anatomical rendering via navis.
- Read the Cell paper's methods in full for any spectral or flow analysis (Section 3 caveat).
- Comparison values on FlyWire v783 (brain only) with the identical pipeline, for H1's
  "brain only" row and for H2's Pospisil check.
- Later: Fish Fire and Wire pairing with ZAPBench activity (structure vs dynamics), once the
  draft connectome is downloadable.

## 8. Reporting commitments

All five hypotheses reported with their pre-stated criteria, pass or fail. Exploratory
results (Section 4) labelled as such. No hypothesis added after 19 Sep 2026 is labelled
confirmatory. Any deviation from Section 6 is listed in a "Deviations" paragraph.
The normalisation artefact behind the signed spectral radius (4.3) is reported as a caution
for the linear-effectome literature, not as a result about photoreceptors.
