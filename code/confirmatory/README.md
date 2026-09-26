# Confirmatory scripts (PREREG Sec. 5)

Written and tested on synthetic data only, before any of H1 to H5 was computed on the connectome.
Every definition the PREREG leaves open is fixed in DEVIATIONS.md, entries C1 to C8.
ASCII, no em dashes.

| file | role |
|---|---|
| `common.py` | shared graph, eigen-solver and bound helpers |
| `verify_inputs.py` | gate: md5 of the four inputs and content hash of the rebuilt graph (C6) |
| `prep_neuropil.py` | per-neuron primary neuropil from syn-partners, for H2(c) (C2) |
| `h1_gap.py` | H1, neck-connective gap (C1) |
| `h2_floor.py` | H2, input floor in the signed normalisation (C2) |
| `h3_depth.py` | H3, certified depth on the weight >= 5 giant SCC (C3) |
| `h4_confidence.py` | H4, robustness under per-synapse confidence classes (C4, D4) |
| `h5_sign.py` | H5, sign-convention sensitivity (C5) |
| `run_confirmatory.sh` | runs everything once, in order; each attempt logs to a new logs/confirmatory_<date>_attempt<k>/ (C6) |
| `requirements-confirmatory.txt` | pinned environment (C8) |
| `tests/test_bounds.py` | brute-force checks of every bound on sampled chains; shows the cert.py t = 0 error |
| `tests/test_confirmatory.py` | runs all scripts on four synthetic data sets and checks every scored number, verdict and descriptive quantity against dense computations, plus the C6 mechanics (bit-identical reruns, no overwrite, unscored-only mode, integrity stop) |
| `tests/make_synthetic.py` | the synthetic CNS (same file formats as the real inputs) |
| `MANIFEST.sha256` | sha256 of every file the runner executes; checked before any computation (C6) |

Run: `MALECNS_DATA=/abs/path/to/data bash code/confirmatory/run_confirmatory.sh` (one to two hours on
2 cores and 8 GB, one BLAS thread; needs the three PREREG inputs plus
syn-partners-male-cns-v1.0-minconf-0.5.feather saved as data/syn-partners.feather).
Tests: `python code/confirmatory/tests/test_bounds.py` and
`python code/confirmatory/tests/test_confirmatory.py` (about 3 minutes, synthetic data only).
