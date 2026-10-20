# Validating and optimising the pharmacophore screen against the μ-opioid receptor

**Run 2026-10-20.** Design and pre-registered metric split: [`STUDY.md`](../STUDY.md).
Machine-readable results: [`selection_results.json`](selection_results.json) (selection),
[`sweep_results.json`](sweep_results.json) (evaluation).

The first implementation of this screen failed its own reproducibility criterion and
showed no robust added value over a one-line sequence rule
([`screen_validation.md`](screen_validation.md)). This study fixes the defects behind
that result and asks whether the method, run properly, does the thing it claims to do.

---

## 1. Three defects, two of them silent

All three were found by measurement, and the two silent ones had been corrupting results
without raising an error.

**Threads.** RDKit's `numThreads=0` means every core the *machine* has, not every core
the *job* was allocated. Eight 8-core array tasks on a 64-core node each spawned 64
threads and drove the node to a measured CPULoad of 319. This is the mirror image of the
`OMP_NUM_THREADS=1` fault that cost this project an earlier campaign.

**Silent under-delivery of conformers.** Distance-geometry embedding quietly returns far
fewer conformers than requested for flexible peptides. Asking for 50 on the N-terminal
window of dynorphin A:

| window | rotatable bonds | conformers returned | with `useRandomCoords` |
|---|---|---|---|
| 5 | 19 | 48 | 50 |
| 8 | 37 | **21** | 50 |
| 10 | 44 | **2** | 50 |
| 12 | 57 | **1** | 50 |

No error is raised. The first screen ran at window 5 and got 48 of 50, which is why this
never surfaced; any longer window would have scored one- and two-member ensembles as
though they were full ones. `useRandomCoords` is now the default and an ensemble below
half the requested size is refused rather than scored.

**Whole-peptide conformers are not viable.** A 17-residue peptide has 81 rotatable
bonds; default embedding returns *zero* conformers for it, and with random coordinates
20 take 73 s, so 500 would cost about half an hour for one peptide. Scoring a whole
digestion peptide against a 5-residue template was never going to work, for reasons of
conformational space rather than cost alone.

## 2. Reproducibility is solved

The blocking defect was not the conformer count alone — it was the **best-of-ensemble
statistic**. The maximum is the noisiest order statistic available, and it was the
default.

| | first run | this study |
|---|---|---|
| conformers | 50 | 400 |
| statistic | max | mean / top-decile |
| seed-to-seed Spearman | — (73% decile churn) | **0.945 – 0.980** |
| median seed jitter ÷ score spread | ~0.16 | **0.006 – 0.05** |

Switching from `max` to `mean` or the top-decile mean matters more than any other single
change: at matched conformer counts it is the difference between failing and passing the
gate.

## 3. What the selection stage chose

Selection used only reproducibility (S1, gate ρ ≥ 0.95) and actives-vs-scrambles
classification (S2), computed on a 72-peptide set containing **only** actives and
scrambles — the evaluation peptides are not in it, so the evaluation metrics could not
have been computed at this stage even by accident.

| window | conformers | metric | statistic | S1 ρ | S2 |
|---|---|---|---|---|---|
| 5 | 400 | tanimoto | top-decile | 0.973 | 0.762 |
| 8 | 400 | tversky_ref | mean | 0.960 | 0.691 |
| 12 | 150 | tanimoto | top-decile | 0.987 | 0.545 |

Ranked by the original rule the winner is window 5. But window 5 collapses the
11-member `YGGFL` family and the 6-member `YGPFL` family to one point each, so the
study's central test is undefined under it. All three windows were therefore carried
into evaluation — an amendment written and committed **before any evaluation metric was
computed**.

## 4. Evaluation

Each window evaluated only at the metric and statistic its selection stage chose.

| | win 5, 400 | win 8, 400 | **win 8, 800** | win 12, 150 |
|---|---|---|---|---|
| S1 reproducibility | 0.968 | 0.945 *(fails gate)* | **0.983** | 0.980 |
| S2 vs scrambles | 0.772 | 0.691 | 0.681 | 0.550 |
| **E1 potency pair** | 1.000 *(2 comparisons — not powered)* | 1.000 *(40 comp.)* | **0.803 [0.561, 0.970]** | **0.348** |
| E2 ρ vs pKi | 0.316 | 0.551 | **0.496** | 0.051 |
| E3 vs Tyr-matched decoys | 0.706 [0.581, 0.817] | 0.726 [0.626, 0.820] | **0.734 [0.632, 0.827]** | 0.668 [0.538, 0.792] |
| E4 endomorphin percentile | 0.671 | 0.976 | **0.988** | 1.000 |
| E5 positional rule | 1.000 | 0.946 | **0.946** | 0.909 |

**The gate-passing configuration is window 8 at 800 conformers**, and it is the one the
conclusions below rest on. Window 8 at 400 conformers is shown because the difference
between the two columns is the study's most important methodological result (§4.2).

### 4.1 The potency pair: the method does the thing it claims to do

**At window 8 with a converged ensemble the screen separates the `YGGFL` family from the
`YGPFL` family at AUROC 0.803, 95% interval [0.561, 0.970] over 40 comparisons.** The
interval excludes chance.

These two families differ by a single **Gly3→Pro3** substitution and about three orders
of magnitude of measured affinity (pKi 7.39–8.71 against 4.71–5.51). They are
length-matched, nearly composition-matched, and every member of both begins with Tyr, so
no sequence or composition rule distinguishes them. Proline at position 3 restricts the
backbone: the difference is conformational, which is precisely what a pharmacophore
method is supposed to read and what the cheaper methods cannot.

This matters because **this project's co-folding pipeline fails the same contrast** —
iLIS scores 0.409, ranking the inactive series *above* the active one. This is the first
method tested here that gets it right.

### 4.2 The perfect score was a sampling artefact, and the gate caught it

At 400 conformers this same contrast scored **AUROC 1.000 with a degenerate [1.000,
1.000] interval**. That configuration fails the reproducibility gate (ρ = 0.945). Doubling
the ensemble to 800 conformers raises reproducibility to 0.983 and drops E1 to 0.803.

The perfect separation was therefore partly an artefact of under-sampling: with too few
conformers the best-matching conformer found for each peptide is effectively a draw from
a long tail, and the draw happened to order the two families cleanly. More sampling gives
a better-estimated and *less* flattering number.

This is the same failure mode as the first study, where a single seed produced an
apparent effect that a second seed removed — and it is the second time a pre-registered
reproducibility criterion has caught an inflated result in this project. **An effect
measured at a configuration that fails the gate should not be quoted**, however clean it
looks. E2 moves the same way, 0.551 to 0.496.

Note that E3 does *not* behave this way: it is 0.726 at 400 conformers and 0.734 at 800,
i.e. stable under better sampling, which is a reason to trust it more than E1's point
estimate.

### 4.3 Added value over the sequence baseline

E3 compares actives against Tyr-starting, length-matched library decoys, so the
"starts with Tyr" baseline is constant across the comparison and cannot contribute.

**Window 8 at 800 conformers reaches 0.734 with a 95% interval of [0.632, 0.827]** —
above 0.70 and excluding 0.5, and stable across conformer counts. The first run's equivalent measurement was 0.603–0.669 with an interval
that included chance. The improvement is attributable to the fixes in §1–2, not to a
change of test.

### 4.4 There is an optimum window, and it is interior

Window 12 **fails** the potency pair (0.348, worse than chance) and has essentially no
potency correlation (ρ = 0.051), despite having the most comparisons available and the
best reproducibility. At 57 rotatable bonds, 150 conformers do not cover the accessible
space, so a well-behaved but meaningless score results — reproducible noise.

Window 5 cannot be evaluated on the central test at all.

The useful window is therefore **interior**: long enough to retain the residues that
distinguish family members, short enough that the conformational space can actually be
sampled. Window 8 is that compromise here.

### 4.5 Potency ranking: a genuine partial positive

E2 at the gated window-8 configuration is **ρ = 0.496** against measured pKi across 36 peptides spanning 4.43 log
units. This project's prior work found potency ranking failing on seven separate metrics
over the same dynamic range, so a moderate positive correlation is a real and
unexpected result. It is reported as a secondary finding: E2 was pre-registered as an
evaluation metric and not used to choose anything, but 36 peptides in a few sequence
families is not a basis for claiming calibrated affinity prediction.

## 5. Honest caveats

- **Window 8 at 400 conformers fails the reproducibility gate** (0.945 against 0.95),
  which it had cleared at 0.960 on the shorter selection set. The gate was not relaxed;
  the ensemble was doubled to 800, which restores ρ = 0.983 — and lowers E1 (§4.2).
- **E1's 40 comparisons come from 8 and 5 distinct windows**, not 40 independent
  molecules. Family members share sequence by construction, so this is one well-posed
  chemical contrast measured many ways, not 40 independent tests.
- **The actives are concentrated in a few families.** `YGGFL`, `YGPFL` and the
  endomorphins account for most of the set, so the benchmark tests a narrow region of
  peptide space well and the rest hardly at all.
- **E4 is percentile against the hard decoy background**, not against the whole library,
  so it is not comparable to the first run's figure.
- Nothing here is experimentally validated. Every number is computational.

## 6. Production configuration and cost

Measured on the 235-peptide benchmark, 8 cores:

| configuration | core-hours (benchmark) | extrapolated to 14,595 library targets |
|---|---|---|
| window 5, 400 conformers | 13.0 | ~810 |
| window 8, 400 conformers | 35.3 | ~2,200 |
| window 12, 150 conformers | 29.3 | ~1,800 |
| **window 8, 800 conformers** | **75.8** | **~4,700** |

Library peptides have median length 15, so nearly all are truncated to the full window;
the benchmark contains many peptides shorter than 8, so the library figure is a lower
bound.

**Window 8 at 800 conformers is the configuration to scale.** It is the only
configuration that both passes the reproducibility gate and can be evaluated on the
study's central test, which it passes. At roughly 4,700 core-hours for the library it
costs about 90× the first run's 52, which bought a result that did not reproduce.

At 48 concurrent 8-core tasks that is approximately 12 hours of wall time.
