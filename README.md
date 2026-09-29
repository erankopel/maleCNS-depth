# maleCNS-depth

**Pre-registered spectral and certified mixing analysis of the male *Drosophila* central nervous
system connectome**

- **Paper:** arXiv:2609.33054 [q-bio.NC], https://arxiv.org/abs/2609.33054
- **Archive:** Zenodo, all versions doi:10.5281/zenodo.22980247 (resolves to the latest release)
- **Author:** Eran Kopel, Tel Aviv University, ORCID https://orcid.org/0000-0003-4657-8636
- **Licence:** code MIT (`LICENSE`); the connectome data keep their own licences (CC BY)

## What the study is

The study asks after how many synaptic steps a random walk on the connectome has forgotten where it
started. The walk runs on the male CNS connectome (male-cns v1.0, 166,700 neurons, 124 million
synapses between neurons), and a walker moves to a postsynaptic partner in proportion to synapse
number. The measure is certified bounds on the Dobrushin coefficient tau(P^r) of the r-step chain.

Five hypotheses, H1 to H5, were registered publicly with numerical pass and fail criteria. The
confirmatory scripts were tested on synthetic data only and frozen. Each hypothesis was then run
once on the connectome. H2 and H5 held; H1, H3 and H4 failed.

The analyses made afterwards to understand the failures are exploratory. They are kept in separate
folders and labelled as not pre-registered.

## Order of events, and how to check it

| When (2026, UTC) | Event | Record |
|---|---|---|
| 19 Sep | Exploratory phase. The pre-registration is drafted and frozen, but not yet public. | `PREREG_maleCNS_depth_2026-09-19.md`, sha256 in `CHECKSUMS.sha256` |
| 26 Sep, 15:56 | Registration made public: release `v0.1.0-prereg`, commit cd16354 | doi:10.5281/zenodo.22980248 |
| 26 Sep | Every definition the pre-registration left open is fixed, before any confirmatory computation. | `DEVIATIONS.md`, entries D1 to D7 and C1 to C8 |
| 26 Sep, 17:18 | Confirmatory scripts frozen: release `v0.2.0-frozen`, commit 7eafd18 | doi:10.5281/zenodo.22981189 |
| 26 Sep, 17:23 to 18:13 | The single confirmatory run (attempt 1): no error, no rerun | `logs/confirmatory_2026-09-26_attempt1/`, `RESULTS_2026-09-26.md` |
| 26 to 27 Sep | Exploratory analyses (not pre-registered) and the manuscript | `EXPLORATORY_2026-09-26.md`, `code/exploratory/`, `paper/` |
| 28 Sep | Release `v0.4.0-exploratory` | doi:10.5281/zenodo.23010800 |

The registration date is 26 Sep 2026, not 19 Sep, because the text was not public before then;
`REGISTRATION.md` explains this.

The release and Zenodo timestamps do not depend on the git history. The integrity checks:

```
sha256sum -c CHECKSUMS.sha256                       # the registered pre-registration and 19 Sep scripts
sha256sum -c code/confirmatory/MANIFEST.sha256      # the frozen confirmatory code
(cd results && sha256sum -c SHA256SUMS)             # the confirmatory result files
```

## Verdicts against the frozen criteria

| | Criterion (pre-registration Sec. 5, `DEVIATIONS.md` C1 to C5) | Value | Verdict |
|---|---|---|---|
| H1 | abs(lambda_2) of the weight >= 5 giant SCC in [0.90, 0.97]; knock-out parts (ii) and (iii) | 0.996524; (ii) and (iii) met | **fail** |
| H2 | floored signed map: spectral radius < 0.85 and at least 10 of the top 20 modes with participation ratio > 20 | 0.83484; 19 of 20 | **pass** |
| H3 | certified lower bound on tau(P^32) > 0.20, and the bound first below 0.05 at r in [48, 96] | 0.9229; still 0.6695 at r = 128 | **fail** |
| H4 | certified lower bound on tau(P'^32) over the confidence class at c = 0.70 > 0.10 | no certificate (-0.520); best level c = 0.51 gives 0.071 | **fail** |
| H5 | spectral radius changes by < 0.02 under the alternative sign convention, same carriers of the top six modes | 0.000015; same six neurons | **pass** (sensitivity) |

`RESULTS_2026-09-26.md` gives the numbers behind each verdict. `EXPLORATORY_2026-09-26.md`
summarises what the exploratory analyses found:

- the neurons that carry the slow mode;
- chains rebuilt from the synapse confidences;
- null models;
- the FlyWire comparison.

## Corrections a reader should know about

- **Before registration: three errors in the pre-registration text** (`DEVIATIONS.md` D4 to D6).
  They were found before registration and are corrected throughout.
  - The robust bounds of its Sec. 4.4 omitted a term.
  - A knock-out count was misstated.
  - Its prior-work section overlooked the BANC connectome.
- **Exploratory: a miscount in one script.** Run 1 of `code/exploratory/e1_slow_mode.py`
  miscounted the synapses crossing the eleven-neuron set; run 2 fixed it. Both runs are kept, and
  `code/exploratory/README.md` explains the error.
- **Exploratory: failed attempts are kept.** They are in `exploratory/logs/` (out of memory, or
  stopped by an error).

## Where things are

| Path | Contents |
|---|---|
| `PREREG_maleCNS_depth_2026-09-19.md` | the pre-registration (unchanged since 19 Sep) |
| `REGISTRATION.md` | what the registration release contains, and the records of each step |
| `DEVIATIONS.md` | every correction and every definition fixed before the run, plus the post-run entry R1 |
| `REPRODUCTION_2026-09-26.md` | the pre-registration's exploratory numbers, reproduced with the unchanged 19 Sep scripts before the confirmatory run |
| `RESULTS_2026-09-26.md` | the confirmatory verdicts and the numbers behind them |
| `EXPLORATORY_2026-09-26.md` | summary of the post-confirmatory exploratory analyses (not pre-registered) |
| `code/*.py` | the eight scripts of the 19 Sep exploratory phase, including `fetch_data.py` |
| `code/confirmatory/` | the frozen confirmatory scripts, their tests, the manifest, the pinned requirements and the runner |
| `code/exploratory/` | the post-confirmatory scripts (see its README) |
| `results/` | confirmatory result files (JSON) with `SHA256SUMS` |
| `logs/repro_2026-09-26/`, `logs/confirmatory_2026-09-26_attempt1/` | logs of the reproduction and of the one confirmatory run |
| `exploratory/results/`, `exploratory/logs/` | outputs and logs of the exploratory scripts |
| `paper/` | LaTeX sources, the compiled `main.pdf`, the figures, and `make_arxiv.py` (builds the single-file arXiv source) |
| `tools/` | runners and helpers; `run_repro.sh` and `fetch_synpartners.sh` contain the paths of the container they ran in |
| `README.txt` | the README of the 19 Sep bundle. It is kept byte-identical because `CHECKSUMS.sha256`, part of the registration record, covers it. Its run order applies to the 19 Sep scripts. |
| `CITATION.cff`, `.zenodo.json`, `LICENSE` | citation and archive metadata, licence |

## Data

The data are not in git; `data/` is created on download.

**Male CNS connectome v1.0** (FlyEM project team at HHMI Janelia and collaborators, CC BY),
flat-connectome release at minimum synapse confidence 0.5, from https://male-cns.janelia.org/download/.

- `code/fetch_data.py` downloads the three input tables (1.1 GB) from the public bucket and checks
  them against the md5 sums frozen in the pre-registration.
- H4 also needs the synapse-partner table, 6.8 GB, md5 58efcf712f8c4d4de5f2ad51e97def76:

  ```
  https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/syn-partners-male-cns-v1.0-minconf-0.5.feather
  ```

  Save it as `data/syn-partners.feather`.

**FlyWire v783**, used only for the exploratory comparison E4. The connections and root ids come
from doi:10.5281/zenodo.10676866, and the annotations from
https://github.com/flyconnectome/flywire_annotations. The md5 sums are in the docstring of
`code/exploratory/e4_flywire.py`.

## Reproducing

**Environment.** The pinned environment is in `code/confirmatory/requirements-confirmatory.txt`:
Python 3.11.15, numpy 2.4.4, scipy 1.17.1, pandas 3.0.2 and pyarrow 25.0.1, with one BLAS thread.
The confirmatory runner stops if the versions differ. The original runs used a Linux container with
2 vCPUs and 7.8 GB of RAM.

1. **Tests.** These use synthetic data only and take about 3 minutes.

   ```
   python3 code/confirmatory/tests/test_bounds.py
   python3 code/confirmatory/tests/test_confirmatory.py
   ```

2. **Confirmatory analysis.** This takes one to two hours and downloads the inputs if they are
   absent. The runner never overwrites results, and it skips any hypothesis whose result file
   exists, so point it at a new folder:

   ```
   MALECNS_DATA=$PWD/data MALECNS_RESULTS=$PWD/results_rerun bash code/confirmatory/run_confirmatory.sh
   (cd results_rerun && sha256sum -c ../results/SHA256SUMS)
   ```

   - **Same environment and hardware:** in the tests, reruns reproduced the result files bit for bit.
   - **Other hardware:** the numbers should agree to the digits reported, even if the files differ
     in the last bits.

3. **Exploratory analyses.** Run them one at a time. Some need about 4 GB of memory, and
   `code/exploratory/README.md` gives the order and the inputs. For example:

   ```
   env MALECNS_DATA=$PWD/data OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 code/exploratory/e1_slow_mode.py
   ```

4. **Paper.**

   ```
   cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
   ```

## Citing

- **Paper:** Kopel E. 2026. Pre-registered spectral and certified mixing analysis of the male
  Drosophila central nervous system connectome. arXiv:2609.33054, doi:10.48550/arXiv.2609.33054.
- **Code and records:** see `CITATION.cff`. Cite either the Zenodo concept DOI
  10.5281/zenodo.22980247 or the DOI of the version you used.

## Use of AI

As stated in the paper, the study and the manuscript were prepared with the help of Claude
(Anthropic). It was used to write and run the scripts, to prepare the pre-registration record, the
deviation log and the notes, to check references and to draft text. The author designed the study,
approved every pre-registration decision, directed and checked the work, verified the results, and
takes full responsibility for the content.
