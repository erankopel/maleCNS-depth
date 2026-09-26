# Registration record

ASCII, no em dashes.

## What this repository's first commit is

The first commit publishes `PREREG_maleCNS_depth_2026-09-19.md` byte-identical to the file
drafted and frozen on 19 Sep 2026 (sha256
34bee5603441cc33c87693861fdf786f5da071917fe94ffd534f565abe4287c4; see `CHECKSUMS.sha256`),
together with the eight scripts of the 19 Sep bundle (`code/`) and the bundle README.

The pre-registration was not public before this commit. Its registration date is therefore the
date of this commit and of the Zenodo record made from the first release (tag `v0.1.0-prereg`),
26 Sep 2026, and not 19 Sep 2026. PREREG Sec. 5 and 8 are read with that date.

## State of the analysis at registration (26 Sep 2026)

- None of the confirmatory hypotheses H1 to H5 (PREREG Sec. 5) has been computed, in any form,
  on any version of the data.
- The only computations made so far are the exploratory runs of PREREG Sec. 4 (19 Sep 2026) and
  their reproduction on 26 Sep 2026 with the unchanged scripts (`REPRODUCTION_2026-09-26.md`,
  logs in `logs/repro_2026-09-26/`).
- `data/syn-partners.feather` (needed for H4) has been downloaded and checked against the md5
  that Google Cloud Storage reports for the object. Only its schema, row count and neuropil
  vocabulary have been read; no confidence value has been looked at.
- Errors and open definitions found before registration are listed in `DEVIATIONS.md`. Further
  entries are added there, dated, before the confirmatory scripts are run.

## What comes next, in order

1. Confirmatory scripts for H1 to H5 and the definitions the PREREG leaves open, tested on
   synthetic graphs only, committed and released (second tag) before any of them is run on the
   connectome.
2. Each hypothesis run once; logs, outputs and pass/fail against the frozen criteria committed.
3. Exploratory work (PREREG Sec. 7 and anything added later) labelled as exploratory.
