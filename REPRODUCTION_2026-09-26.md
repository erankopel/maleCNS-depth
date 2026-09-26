# Reproduction of the PREREG Sec. 1, 2 and 4 numbers, 26 Sep 2026

The eight scripts of the 19 Sep bundle were run unchanged, in the README order (steps 1 to 9), on a
fresh download of the inputs. Runner: `tools/run_repro.sh`; logs: `logs/repro_2026-09-26/`.
Environment (`00_environment.log`): cloud Linux container, 2 vCPU Intel Xeon 2.1 GHz, 7.8 GB RAM,
Python 3.11.15, numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, pyarrow 25.0.1. Total wall time 10 min;
peak memory 2.3 GB (modes.py). Exploratory only; nothing from H1 to H5 was computed.
ASCII, no em dashes.

Every number below matches the PREREG to the digits printed there. The only differences are the
sign of the slow-mode eigenvector (arbitrary; see 4.2) and two errors in the PREREG text that the
logs expose (DEVIATIONS.md D4 and D5).

## Sec. 1, inputs

| file | md5 frozen in PREREG | 26 Sep download |
|---|---|---|
| body-annotations (14,483,314 B) | 50a7718770c57220f160ba4f431ab89e | match |
| body-neurotransmitters (43,282,834 B) | 3d842b12fe5c49eefade528d7dd24a1f | match |
| connectome-weights (1,051,241,946 B) | f30e9dcca25cfd021bf1e7b3d975599e | match |

Also downloaded, for H4: syn-partners-male-cns-v1.0-minconf-0.5.feather, 6,777,179,098 B, md5
58efcf712f8c4d4de5f2ad51e97def76 (matches the md5 reported by Google Cloud Storage; object
generation 1780494942562468, last modified 3 Jun 2026), 311,833,243 rows (equal to the synapse
total of the connectome-weights table), 4,759 record batches. Columns: x_pre, y_pre, z_pre,
body_pre (int64), conf_pre (float32), x_post, y_post, z_post, body_post (int64), conf_post
(float32), primary_post (dictionary of 146 neuropil names, including LOP(L), LOP(R) and EB).
No confidence value was read.

## Sec. 2, objects

| quantity | PREREG | reproduced | log |
|---|---|---|---|
| annotated bodies | 211,577 | 211,577 | 02 |
| neurons N (non-null superclass) | 166,700 | 166,700 | 02 |
| segment-to-segment rows | 151,856,684 | 151,856,684 | 02 |
| synapses in those rows | 311,833,243 | 311,833,243 | 02 |
| N x N edges / synapses | 25,582,938 / 124,177,617 | 25,582,938 / 124,177,617 | 02 |
| weight >= 5 edges / synapses | 6,242,118 / 89,860,280 | 6,242,118 / 89,860,280 | 02, 03 |
| giant SCC G (w >= 1) | 165,314 (99.2%), 25,545,360 edges | 165,314 (99.2%), 25,545,360 edges | 03, 04 |
| giant SCC at w >= 5 | 157,821 (94.7%) | 157,821 (94.7%) | 03 |
| signed neurons + / - / excluded | 103,720 / 59,262 / 3,718 | 103,720 / 59,262 / 3,718 | 05 |
| median NT confidence | 0.936 | 0.936 | 05 |

## Sec. 4.1, unsigned chain on G

| quantity | PREREG | reproduced | log |
|---|---|---|---|
| lambda_2 | 0.94503 | 0.945027 | 04 |
| spectral gap | 0.0550 | 0.054973 | 04 |
| TV at r = 4, 8, 16, 32, 64, 100 (20 sources, whole graph) | 0.77, 0.57, 0.33, 0.11, 0.018, 0.0045 | 0.7704, 0.5694, 0.3303, 0.1145, 0.0182, 0.0045 | 05 |
| stationary mass, top 1% / top 10% | 0.422 / 0.854 | 0.422 / 0.854 | 04 |
| stationary entropy | 13.96 of 17.33 bits | 13.96 of 17.33 bits | 04 |
| VNC stationary mass (20,109 neurons) | 0.653 | 0.653 | 07 |
| 674 motor neurons: mass, median in / out synapses | 0.168, 2,729 / 15 | 0.168, 2,729 / 15 | 08 |
| whole-graph spurious eigenvalue | 0.9986 | 0.998614 (plus four eigenvalues of modulus 1) | 05 |

## Sec. 4.2, slow-mode anatomy

| quantity | PREREG | reproduced | log |
|---|---|---|---|
| slow mode: share of one sign in VNC / of the other in brain | 0.921 / 0.938 | 0.921 / 0.938 | 07 |
| top carriers | leg motor neurons and APL, opposite signs | same | 07 |
| motor neurons removed | 0.9430 | 0.9430 | 08 |
| ascending + descending removed | 0.9997 | 0.9997 (3,693 removed, SCC 161,597; see D5) | 08 |
| brain only | 0.9256 | 0.9256 | 08 |
| escape probabilities brain to VNC / VNC to brain | 0.144 / 0.071 | 0.1443 / 0.0708 | 07 |
| two-block prediction 1 - p - q | 0.78 | 0.7849 | 07 |

The eigenvector returned by ARPACK has the opposite overall sign to the 19 Sep run (here the brain
side is positive and the motor neurons negative). An eigenvector's sign is arbitrary; the shares
are identical.

## Sec. 4.3, signed map on G

| quantity | PREREG | reproduced | log |
|---|---|---|---|
| spectral radius | 0.9151 | 0.915062 | 06, 07 |
| leading eigenvalues | -0.9151, +0.9137, -0.9122, +0.9098 | -0.915062, +0.913732, -0.912235, +0.909779 | 06 |
| participation ratio, top six modes | 2.0 to 2.4 | 2.0 to 2.4 | 07 |
| mode power on the top 10 neurons, top six modes | > 0.97 | 0.970 to 0.994 | 07 |
| carriers | R7/R8 photoreceptor pairs | R7/R8 pairs (R7_unclear/R8_unclear, R7p/R8p) | 07 |
| modes 7 and 9 (1-based) | PR 76 and 97; R7d/R8d, ER3 | PR 76.1 and 96.6; R7d/R8d, ER3d_b/ER3p_a | 07 |
| mean L1 norm of M^r e at r = 1, 4, 8, 16 | 0.86, 0.18, 0.016, 2.4e-4 | 0.859, 0.185, 0.0159, 2.4e-4 | 06 |

## Sec. 4.4, certified bounds on tau(P^r), G, 40 witnesses, 60 hub columns

| r | lower, PREREG | lower, reproduced | upper, PREREG | upper, reproduced |
|---|---|---|---|---|
| 4 | 0.9896 | 0.9896 | 1.0000 | 1.0000 |
| 8 | 0.8794 | 0.8794 | 0.9997 | 0.9997 |
| 16 | 0.6231 | 0.6231 | 0.9936 | 0.9936 |
| 32 | 0.2663 | 0.2663 | 0.9697 | 0.9697 |
| 64 | 0.0428 | 0.0428 | 0.9434 | 0.9434 |

Robust bounds under C(0.05) as printed by cert.py: lower 0.7895 (r = 2), 0.3583 (r = 4), zero from
r = 8; upper no better than 0.9987. These match the PREREG, but the robust lower bounds are wrong
because cert.py omits the t = 0 term; the corrected values are 0.579 (r = 2), 0.148 (r = 4) and
zero from r = 8 (DEVIATIONS.md D4). Witness pair at r = 64: body IDs 533585 and 801768.
