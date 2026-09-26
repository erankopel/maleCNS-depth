maleCNS-depth: certified information-breaking depth of the male Drosophila CNS connectome
Exploratory phase run 19 Sep 2026 in the claude.ai chat container; bundle prepared 26 Sep 2026.
ASCII, no em dashes.

Suggested place on the PC: Documents/Arxiv/ (the zip unpacks to a folder maleCNS-depth/).

CONTENTS
  PREREG_maleCNS_depth_2026-09-19.md  Frozen pre-registration (H1 to H5). Byte-identical to the
                                      19 Sep file. Do not edit; deviations go in a separate note.
  code/                               The seven analysis scripts of 19 Sep plus fetch_data.py.
  data/                               Not shipped. Created by fetch_data.py (1.1 GB of inputs,
                                      plus about 115 MB of intermediate files from the scripts).
  CHECKSUMS.sha256                    sha256 of every shipped file ("sha256sum -c CHECKSUMS.sha256").

ONE CHANGE FROM THE 19 SEP SCRIPTS
  The container path /home/claude/fly/ is replaced by DATA, which defaults to ../data relative to
  code/ and can be overridden with the environment variable MALECNS_DATA. No analysis line changed.
  Tested 26 Sep: all scripts compile and reach their data; fetch_data.py downloaded and verified.

RUN ORDER  (Python 3 with numpy, scipy, pandas, pyarrow; each step is minutes on 1 CPU / 3 GB)
  1. python code/fetch_data.py         3 inputs from the public bucket, md5 checked vs PREREG Sec. 1
  2. python code/build_graph.py        -> data/neuron_edges.npz (166,700 neurons, 25,582,938 edges)
  3. python code/analyse.py            SCCs, degrees, reachability, first TV curve (to 20 hops)
  4. python code/spectra_scc.py        -> data/scc_keep.npy, data/scc_pi.npy; lambda_2 = 0.94503
  5. python code/spectra.py unsigned   whole-graph eigenvalues (spurious 0.9986 from a closed
                                       component) and the 200-hop TV curve of PREREG 4.1
  6. python code/spectra.py signed     signed map on the giant SCC; spectral radius 0.9151
  7. python code/modes.py              PREREG 4.2 and 4.3 -> data/signed_modes.npz
  8. python code/slowmode_check.py     knock-outs: motor neurons, neck connective, VNC
  9. python code/cert.py               certified bounds on tau(P^r) to r = 64 (PREREG 4.4)
  Steps 6 to 9 need step 4. Step 5 needs only step 2.

NEXT (PC PHASE; PREREG Sections 5 and 7)
  - Before running anything for H1 to H5: commit PREREG and code/ to the GitHub repo (PREREG
    Sec. 6 names a new fly/ directory) and push a tag. The frozen date only counts once it has a
    public timestamp; nothing from H1 to H5 has been run as of 26 Sep 2026.
  - H1 to H4 need variants of these scripts: weight >= 5 giant SCC (H1, H3), an input floor in
    the signed normalisation (H2), 100 witness rows and checkpoints to 128 hops (H3),
    per-synapse confidence intervals from syn-partners (H4).
  - Heavy downloads: syn-partners (6.8 GB) and syn-points (12.7 GB) from
    https://male-cns.janelia.org/download/ . neuPrint needs an account token; not required.
  - Read the male CNS Cell paper (Sept 2026) in full for any spectral or flow analysis
    (PREREG Sec. 3 caveat).
