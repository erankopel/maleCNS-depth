# Exploratory analyses after the confirmatory run (26 September 2026)

**Not pre-registered.** These analyses were made after the single confirmatory run of H1-H5
(`RESULTS_2026-09-26.md`) to understand its outcome. They do not change any confirmatory verdict
(H1 fail, H2 pass, H3 fail, H4 fail, H5 pass). Code: `code/exploratory/` (see its README);
results: `exploratory/results/`; logs: `exploratory/logs/`; write-up: Section 4 of `paper/main.pdf`.
Same data, environment and one-BLAS-thread setting as the confirmatory run. ASCII, no em dashes.

## E1: what carries the slow mode (why H1 (i) and H3 failed)

- At weight >= 5 the slow mode (|lambda_2| = 0.996524) sits on 11 neurons of the cord: seven efferent
  neurons with somata on the mesothoracic midline (EN00B001, EN00B008, EN00B011 x2, EN00B015 x2,
  mesVUM-MJ), IN03B088 x2, IN19A061 and the motor neuron MNad21. They hold 12.8 % of the stationary mass;
  a walker leaves them with probability 3.5e-3 per step; 1 - p - q = 0.99598 reproduces lambda_2.
- At weight >= 1 the same set holds 0.4 % of the mass and is left with probability 0.29 per step. It
  sends 299 synapses to the rest of the CNS (43 on connections of >= 5 synapses) and receives 21,846
  (19,378 on such connections): the threshold removes most of the way out and little of the way in.
- The H3 witness pair that attains the bound from r = 4 has one member in the set.
- Without the set, lambda_2 = 0.9806 and the next metastable set is again small and in the cord
  (5 VNC interneurons and 1 motor neuron); certified 0.674 <= tau(P^32) <= 0.863 and
  0.107 <= tau(P^128) <= 0.511. Removing all efferent superclasses (24 neurons) gives the same.
- At weight >= 1 no small set dominates: the smallest-conductance cut is 727 neurons of the central
  complex and anterior visual pathway (conductance 0.116).
- Correction: run 1 of `e1_slow_mode.py` miscounted the synapses across the set (see the README);
  the numbers above are from run 2 and `check_trap_synapses.py`.

## E2: is the H4 failure the certificate or the data?

- Exact chains built from the confidence data (`e2_class_members.json`): M(c) removes every synapse
  below confidence c (a vertex of the class at level c); B(k) keeps each synapse with probability
  equal to the midpoint of its confidence bin (five draws). At c = 0.70 no witness distance at r = 32
  changes by more than 0.017 (max d_32: 0.312 -> 0.296); B(k): at most 0.016; M(0.90), which removes
  half of all synapses: 0.075. The certified H4 bound at c = 0.70 is negative from r = 8.
- |lambda_2| at weight >= 1 moves by less than 0.008 in every chain. At weight >= 5 it follows the
  11-neuron set (`e2c_trap_in_members.json`): in M(0.80) and M(0.90) the set leaves the giant SCC
  and |lambda_2| is 0.940 and 0.952; in the random reconstructions it is 0.994 and 0.996 where the set
  survives with its mass and 0.967 to 0.976 where it is broken up.
- E2b (`e2b_sharp_certificate*.json`): the exact largest row change over the class (Lemma 3 of the
  paper; `radius.py`) is on average 0.089 at c = 0.70 against the H4 radius 0.144. With it the
  certified bound at c = 0.70 is 0.621 at r = 4 (H4: 0.396) and still negative at r = 8; at c = 0.55
  it certifies tau(P'^16) >= 0.106 (H4: none). Evaluated with the actual row changes of M(0.70), the
  first-order bound still gives only 0.026 at r = 8, where M(0.70) has 0.934: most of the slack is in
  the step-by-step accumulation of Lemma 2, not in the radius. At r = 32 the best certificate, at
  c = 0.51, rises from 0.071 to 0.096, still short of the pre-registered 0.10; from c = 0.52 on there is
  none. The H4 column of E2b reproduces `results/h4_full.json` exactly at all ten levels computed.

## E3: null models (`e3_null_models.json`)

- Configuration null N1: |lambda_2| about 0.15 (w >= 1) and 0.23 (w >= 5). Superclass-block null N2:
  0.928 and 0.91 to 0.94. Cell-type-block null N3: 0.9446 at w >= 1 (observed 0.9450, with a cut of
  714 to 716 neurons against 727), and 0.998, 0.996, 0.953 at w >= 5 (a trap of 11 to 13 neurons in two
  of three draws). The slow modes are a property of cell-type wiring.

## E4: the FlyWire female brain (v783) and the male brain alone (`e4_flywire_*.json`)

- The slowest well-populated subsystem is the anterior visual pathway and central complex in both
  brains (cuts of 891 to 1,495 neurons, conductance 0.07 to 0.11).
- Certified lower bounds on tau(P^32): 0.285 and 0.261 (FlyWire, w >= 1 and >= 5), 0.112 and 0.165
  (male brain), against 0.312 for the whole male CNS at w >= 1.
- FlyWire at w >= 1 has |lambda_2| = 0.99999, set by seven optic-lobe neurons at the lamina (four
  R1-6 photoreceptors, one L2, two unlabelled) with 1e-4 of the stationary mass: a boundary set, like
  the male cord set, that dominates lambda_2 but not the depth profile.
