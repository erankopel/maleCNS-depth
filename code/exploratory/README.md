# Post-confirmatory exploratory analyses (26 September 2026)

**Not pre-registered.** Everything in this directory was written and run *after* the single
confirmatory run of H1-H5 (results in `results/`, verdicts in `RESULTS_2026-09-26.md`), to understand
its outcome. None of it changes a confirmatory verdict. These scripts are distinct from:

- `code/*.py`: the exploratory phase of 19 September 2026 on which the pre-registration was based;
- `code/confirmatory/`: the frozen confirmatory scripts (release v0.2.0-frozen,
  doi:10.5281/zenodo.22981189), which nothing here modifies.

The scripts import helper functions from `code/confirmatory/common.py` but never write to `results/`.
Outputs go to `exploratory/results/`, logs to `exploratory/logs/`, figures and generated table rows
to `paper/`. A summary is in `EXPLORATORY_2026-09-26.md` at the top of the repository.

| Script | Question | Output (`exploratory/results/`) |
|---|---|---|
| `e0_edge_confidence.py` | per-edge synapse counts of G in ten confidence bins (input to E2) | `data/edge_conf_hist.npy`, `data/edge_conf_keys.npy` (not in git, 1.2 GB) |
| `e0b_low_confidence_bins.py` | the same, five bins below 0.55 (input to E2b low) | `data/edge_conf_hist_low.npy` (not in git) |
| `e1_slow_mode.py` | which neurons carry the slow mode at w>=5 and at w>=1 (eigenvectors, sweep cut) | `e1_slow_mode.json` |
| `check_trap_synapses.py` | independent count of the synapses crossing the 11-neuron set | `check_trap_synapses.json` |
| `e1b_without_trap.py` | the w>=5 chain without the 11-neuron set, and without all efferent superclasses | `e1b_without_trap.json` |
| `e2_class_members.py` | is the H4 failure the certificate or the data? exact chains M(c), B(k) | `e2_class_members.json` (+ `.jsonl`, `e2_P_witness_d.npz`) |
| `e2b_sharp_certificate.py` | the H4 bound with the exact row radius (Lemma 3) and with the row changes of M(c) | `e2b_sharp_certificate.json`, `e2b_sharp_certificate_low.json` (argument `low`) |
| `radius.py`, `test_radius.py` | the exact row radius, and its test against brute force | |
| `e2c_trap_in_members.py` | the 11-neuron set at w>=5 in the chains of E2 | `e2c_trap_in_members.json` |
| `e3_null_models.py` | stub-matching null models N1-N3 for lambda_2 | `e3_null_models.json` |
| `e4_flywire.py` | the same chain on FlyWire v783 and on the male brain only | `e4_flywire_fw.json`, `e4_flywire_male.json` |
| `fig_depth.py`, `fig_members.py` | Figures 1 and 2 of the paper | `paper/fig_depth.pdf`, `paper/fig_members.pdf` |
| `make_tab_members.py`, `make_tab_cert.py` | rows of Tables 4 and 5 of the paper, written from the result files | `paper/tab_members_rows.tex`, `paper/tab_cert_rows.tex` |

Run each script from the repository root with the environment of the confirmatory run:

```
env MALECNS_DATA=$PWD/data OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 code/exploratory/e1_slow_mode.py
```

E2 and E2c need the output of `e0_edge_confidence.py`, E2b low that of `e0b_low_confidence_bins.py`;
`e4_flywire.py` needs the FlyWire v783 files in `data/flywire/` (md5 sums in its docstring) and takes
`flywire`, `male` or `both`. Peak memory is up to about 4 GB per script; run them one at a time on a
machine with 8 GB (two at once were killed for lack of memory; see the `*_oom.log` files).

## Corrections made during the exploratory phase

- `e1_slow_mode.py`, run 1 (`e1_slow_mode_run1.json`, `.log`): the synapse counts across the
  11-neuron set (step 4) paired the weights of the CSR matrix `A1` with the column indices of the
  row-normalised matrix `P1`, which stores the entries of a row in a different order; the four
  counts were wrong (617, 378, 13,366, 10,711). Found by `check_trap_synapses.py`, which counts from
  the raw edge list; run 2 (`e1_slow_mode.json`) uses `A1`'s own indices and gives 299, 43, 21,846,
  19,378. Nothing else in the output changed.
- `e2_class_members.py`: attempt 1 stopped on a shared index array (fixed by copying it), attempt 2
  was killed for lack of memory; the script was made resumable and attempt 3 completed.
- `e4_flywire.py`: attempt 1 was killed for lack of memory while E2 ran; it was split per dataset.
