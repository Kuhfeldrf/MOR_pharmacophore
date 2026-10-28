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

**Measured:** 2,800 core-hours, 8.2 h wall over 300 array blocks (298 non-empty; two
blocks received zero peptides because 49 × 298 > 14,595, and the empty-table guard
refused them rather than writing a header-only file). The benchmark extrapolation was
68% high because the benchmark contains many peptides shorter than the window, which
embed faster per peptide than a truncated 8-mer. All 14,571 library peptides scored.

---

## 7. The full library screen, and the confound it exposed

Ranking the library at the validated configuration produced a result that the
pre-registered benchmark had no way to catch. **This section is post-hoc** — everything
in it was written after the library ranking existed — and is labelled as such throughout.
Machine-readable: [`library_bare_report.json`](library_bare_report.json),
[`composition_control.json`](composition_control.json),
[`library_bare_shortlist.json`](library_bare_shortlist.json).

### 7.1 No known active reaches tier 1

Tiers were anchored to the Tyr-matched decoy background, not to library deciles. Edges
were 0.3998 / 0.3787 / 0.3751 (decoy p99/p95/p90) against a seed jitter of 0.0013, so
**every edge clears the jitter** — the tier-safety check passes, which was the question
§6 was built to answer.

It is the contents that are wrong:

| | |
|---|---|
| tier 1 | 22 peptides — **8 distinct windows**, 2 Tyr-starting |
| tier 2 | 133 peptides, 75 windows |
| tier 3 | 201 peptides, 105 windows |
| untiered | 14,370 |

Best endogenous ligand, endomorphin-1, scores 0.3875 — the **99.7th library percentile**,
but *below* the tier-1 cut. Leu-enkephalin is at 97.3%, the dynorphins at 95.4%, the
template's own met-enkephalin at 95.6%. So 22 library peptides outrank every known
opioid peptide, and the eight windows they come from are `EFHSHEFH`, `FHSHEFHS`,
`HSHEFHSH`, `SHEFHSHE`, `FHWKGHEH`, `EKIWHHTF`, `YKYDHNFV`, `YSTYHTHG` — one His-rich
repeat family and two Tyr/His windows.

### 7.2 The score is close to a staircase in aromatic count

Spearman ρ between the score and simple composition of the scored window, over all
14,571 library windows:

| | ρ |
|---|---|
| **aromatic residues (F/W/Y/H)** | **+0.909** |
| His count | +0.434 |
| matched features | +0.350 |
| starts with Tyr | +0.240 |
| window length | −0.069 |
| charged (K/R/D/E/H) | −0.006 |

The strata barely overlap:

| aromatics | distinct windows | median | max |
|---|---|---|---|
| 0 | 2,469 | 0.2840 | 0.3020 |
| 1 | 2,077 | 0.3382 | 0.3519 |
| 2 | 703 | 0.3583 | 0.3775 |
| 3 | 144 | 0.3766 | 0.3906 |
| 4 | 19 | 0.3965 | 0.4066 |
| 5 | 5 | 0.4098 | 0.4182 |

Every endogenous ligand sits at 2 or 3 aromatics. **Every tier-1 window sits at 4 or 5** —
strata of 19 and 5 windows, containing no known active. His supplies an aromatic ring, a
donor and an acceptor at one residue, so a His-rich repeat accumulates template coverage
without presenting the message motif at all.

This does not touch E1. The `YGGFL` and `YGPFL` families are exactly aromatic-matched, so
§4.1's 0.803 is a conformational result and stands. E3 is partly affected: its decoys are
Tyr-matched, which fixes the first residue but not the total aromatic count.

### 7.3 A remedy that failed its own control, and is therefore not used

Percentile stratification — score converted to its rank within the library windows of the
same aromatic count — was implemented ([`14_composition_control.py`](../scripts/14_composition_control.py))
with E1 as its invariance control, and **rejected on two independent grounds**:

| | raw | stratified |
|---|---|---|
| E1 potency pair | 0.803 | **0.712** |
| E2 ρ vs pKi | 0.497 | 0.356 |
| E3 vs Tyr-matched decoys | 0.733 | 0.643 |
| E4 endomorphin percentile | 0.988 | 0.940 |
| E5 positional rule | 0.946 | 0.959 |

E1 should have been *unchanged* and moved by 0.09, because at window 8 the family members
do not all share an aromatic count (`YGGFLRRI` has 2, `YGGFLRKY` has 3), so stratifying
splits a contrast it was meant to leave alone. And it promoted zero-aromatic peptides
(`LTTLKNLRETAE`, raw 0.3019) to the 99.96th percentile for being the best of a stratum in
which nothing matches the template. Aromatic count is not a nuisance variable here: the
template's anchors *are* the Tyr1 and Phe4 aromatics of the message.

Recording the rejected attempt rather than deleting it, because the reason it failed is
the clearest statement of what the confound actually is.

### 7.4 What is reported instead

Rank only within the strata the endogenous ligands occupy (2 and 3 aromatics), against a
bar set by the weakest endogenous ligand in the *same* stratum
([`15_shortlist.py`](../scripts/15_shortlist.py)):

| stratum | bar | library windows clearing it |
|---|---|---|
| 2 aromatics | dynorphin A, 0.3623 | 355 of 703 |
| 3 aromatics | endomorphin-1, 0.3875 | 8 of 144 |

**363 peptides, 212 distinct windows, 72 Tyr-starting** —
[`library_bare_shortlist.tsv`](library_bare_shortlist.tsv). The 52 peptides in the 4- and
5-aromatic strata are listed separately and flagged as uncalibratable
([`library_bare_high_aromatic_flagged.tsv`](library_bare_high_aromatic_flagged.tsv)) rather than
ranked alongside.

The 2-aromatic bar is weak: half that stratum clears the weakest endogenous ligand. That
is the library-scale restatement of E3 = 0.734 — real separation, not a sharp filter. The
3-aromatic stratum is the selective one, 8 of 144 clearing endomorphin-1, and those eight
are the part of this output worth anything downstream.

### 7.5 What this changes about the framework

The study's conclusion in §4 is unaltered: with composition held fixed the method reads
conformation, at AUROC 0.803 on a three-log-unit contrast that this project's co-folding
pipeline gets backwards. What §7 adds is that **the pre-registered benchmark could not
have detected the library-scale confound**, because every one of its decoy sets controls
composition by construction — scrambles are composition-matched, Tyr-decoys are
Tyr-matched and length-matched. A benchmark built entirely from matched controls measures
the thing you care about and is blind to the thing that dominates an unmatched library.

The framework needs a sixth evaluation metric, pre-registered for the next receptor:
**E6, the composition confound** — ρ between score and aromatic count across the screened
library, with a stated ceiling above which raw cross-library tiers are not reported. At
+0.909 this screen is far above any ceiling worth setting.

Window 8 remains the configuration. The raw tiers do not survive; the within-stratum
shortlist does.

---

## 8. What the score actually reads

Three changes, all post-hoc, all driven by §7. Machine-readable:
[`reselection_e6.json`](reselection_e6.json), [`epimer_test.json`](epimer_test.json),
[`effect_sizes.json`](effect_sizes.json).

### 8.1 The confound cannot be fixed by choice of metric

All nine metric × statistic combinations, scored on the confound and the pre-registered
metrics together:

| metric / statistic | E6 (arom ρ) | S1 | S2 | E1 | E2 | E3 |
|---|---|---|---|---|---|---|
| tversky_ref top10 | +0.910 | 0.883 | 0.711 | 0.879 | 0.550 | 0.679 |
| **tversky_ref mean** | **+0.913** | **0.983** | **0.681** | **0.803** | **0.496** | **0.734** |
| tanimoto mean | +0.722 | 0.989 | 0.573 | 0.470 | 0.121 | 0.675 |
| tanimoto top10 | +0.730 | 0.934 | 0.593 | 0.515 | 0.146 | 0.614 |
| tversky_db top10 | +0.334 | 0.988 | 0.555 | 0.470 | 0.087 | 0.609 |
| tversky_db mean | +0.310 | 0.998 | 0.534 | 0.439 | 0.074 | 0.639 |

**E6 and the evaluation metrics are collinear.** Every metric that weakens the confound
weakens E1 and E2 in proportion, and the two combinations that clear a ρ = 0.6 ceiling
sit at chance on the study's central test. The ceiling committed in STUDY.md §7 was
therefore withdrawn the same day: a gate that eliminates every working configuration is
badly posed, and E6 now constrains the *comparison the score may make* rather than which
score is used.

The mechanism is not a defect. TVERSKY_REF measures coverage of the template's
pharmacophore; the template's anchors are the Tyr1 and Phe4 aromatics of the message;
Gly3→Pro changes how well those can be presented. Tracking aromatic coverage and
separating the potency pair are one property.

### 8.2 E7: the screen does not resolve stereochemistry

Since no metric can decouple composition from shape, the benchmark has to. An L/D epimer
pair is the limiting case — identical atoms, identical aromatic count, identical
constitution, differing only in the handedness of one α-carbon — so any separation is
necessarily conformational. `conformers.py` now builds these (lowercase = D).

| pair | D | L | Δ | direction |
|---|---|---|---|---|
| dermorphin | 0.3826 | 0.3817 | +0.0009 | correct |
| [D-Ala²]Met-enkephalinamide | 0.3673 | 0.3691 | −0.0018 | wrong |
| [D-Ala²]Leu-enkephalinamide | 0.3662 | 0.3719 | −0.0057 | wrong |
| deltorphin I *(δ-selective, excluded)* | 0.3758 | 0.3750 | +0.0008 | — |
| deltorphin II *(δ-selective, excluded)* | 0.3691 | 0.3699 | −0.0008 | — |

**1 of 3 μ-relevant pairs in the predicted direction, exact two-sided p = 1.00.** The
deltorphins are δ-selective, carry no directional prediction against a μ template, and
are excluded from the sign test rather than counted in whichever direction helps.

The direction for L-Ala² comes from qualitative SAR, not a matched pKi pair, so this is
scored as a sign test and no interval is invented for a quantity with no measured scale
here.

### 8.3 Effect sizes against the noise floor — the reporting change that matters most

An AUROC says how consistently two groups are ordered and nothing about how far apart.
This project has been misled by that gap twice: 73% decile churn in the first study, and
E1 = 1.000 at 400 conformers in this one. So every contrast is now reported as a margin
in units of the measured seed jitter (0.0013 at this configuration):

| contrast | AUROC | margin | in jitter | status |
|---|---|---|---|---|
| S2 actives vs scrambles | 0.680 | +0.0078 | 5.9× | resolved |
| **E1 potency pair** | **0.803** | **+0.0043** | **3.2×** | **resolved** |
| **E3 vs Tyr-matched decoys** | **0.733** | **+0.0119** | **9.0×** | **resolved** |
| E5 positional rule | 0.946 | +0.0041 | 3.1× | resolved |
| **E7 epimer pairs** | — | **−0.0018** | **−1.4×** | **noise-floor limited** |

E1 survives this test at 3.2× the jitter, and **E3 is the strongest contrast in the
study at 9.0×** — which is consistent with §4.2, where E3 was the metric stable under
better sampling while E1 was not. E7 is the only contrast below the floor, and it points
the wrong way.

### 8.4 The conclusion this forces

The three contrasts form a ladder, and reading them together is more informative than any
one of them:

| what changes | example | effect |
|---|---|---|
| aromatic count | library-wide | ρ = 0.909, dominates |
| feature composition at fixed aromatic count | Gly3→Pro (Pro adds a lipophilic ring, removes a backbone donor) | E1 = 0.803, 3.2× jitter |
| nothing but 3D arrangement | L/D epimer | p = 1.00, below the noise floor |

**The score tracks pharmacophore feature composition, with no resolvable contribution
from 3D arrangement.** That is a narrower claim than §4.1 made and it supersedes it.
E1 is real and is not the aromatic confound — Gly3→Pro is aromatic-matched — but it is
also not evidence of shape reading, because Pro changes the feature inventory. The one
contrast that isolates shape returns nothing.

This is a negative result for the method *as a 3D shape-matching method*, and it stands
regardless of settings, because it was measured at the configuration that passes every
gate. What the screen is, demonstrably, is a well-calibrated pharmacophore-feature
composition filter: E3 at 9× the noise floor against Tyr- and length-matched decoys is a
real and useful capability, and it is not what "pharmacophore matching" is usually taken
to claim.

Practical consequence for the library output: the within-stratum shortlist of §7.4 stands
as the deliverable, and should be read as "peptides whose feature inventory matches the
opioid message template", not as "peptides predicted to adopt the message conformation".
The 8 peptides clearing endomorphin-1 in the 3-aromatic stratum remain the part worth
following up, on that narrower basis.

*§8.4's conclusion is further narrowed by §9.2 below: the feature-composition claim
itself does not survive a matched benchmark.*

---

## 9. Checked against the methodology literature

Everything to §8 was decided on this project's own measurements. Checking those choices
against the virtual-screening methodology literature found two metric errors, and
vindicated one choice. Machine-readable:
[`screening_metrics.json`](screening_metrics.json).

### 9.1 AUROC is the wrong summary for a screen

Truchon & Bayly (*J Chem Inf Model* 2007, 47:488–508) show ROC AUC is equivalent to an
average of the actives' ranks, so early recognition is cancelled by late recognition —
the wrong summary for a method whose entire purpose is the top of a list. They propose
**BEDROC** with α = 20, putting 80% of the weight in the top 8%.

AUROC is not wrong everywhere, and the distinction matters: **E1 is a two-sample
comparison**, where AUROC is the Mann-Whitney *U* statistic and is exactly right. The
error is confined to the screening-shaped questions. Those, computed properly:

| ranking all 14,595 scored peptides | AUROC | BEDROC α=20 | null p99 | EF 1% | EF 5% |
|---|---|---|---|---|---|
| 24 reference ligands | 0.915 | **0.322** | 0.136 | 8.3 | 9.2 |
| 7 endogenous + template only | 0.964 | **0.517** | 0.214 | 14.3 | 20.0 |

Both clear a permutation null (computed rather than assumed, because 36 actives is below
the ~50 Truchon & Bayly recommend for stable BEDROC). Best-ranked reference ligand is
43rd of 14,595. **Read alone this is the study's most favourable result** — and §9.2 is
why it cannot be read alone.

A second lesson from the same source: BEDROC on the *benchmark* contrasts is not
meaningful. With 36 actives against 85 decoys, N = 121 and α = 20 concentrates all weight
in the top ~10 ranks. E3's BEDROC is 0.610 against a null p99 of 0.612 — i.e. not
distinguishable from random. That is a statement about the metric being misapplied at
that sample size, not a result.

### 9.2 The decoys were not matched on the property that drives the score — and E3 does not survive matching

`tyr_hard` matches the actives on Tyr-start and length but **not on aromatic count**, the
property correlating with the score at ρ +0.909. The benchmark-design literature names
this failure mode: **artificial enrichment**, performance attributable to the difference
in chemical space between actives and decoys rather than to binding (*Bioinformatics*
2021, 37:2134; and *PLOS ONE* 2019 e0220113 on DUD-E's hidden bias, which distinguishes
analogue bias from decoy bias). DUD-E's own principle is property similarity with
topological dissimilarity; this benchmark satisfies it for length and the N-terminal
residue and violates it for the one property that matters.

Recomputing E3 within matched aromatic strata — actives against Tyr-starting library
decoys of the *same* aromatic count, distinct windows only:

| aromatic count | actives | matched decoys | AUROC |
|---|---|---|---|
| 1 | 2 | 118 | 0.011 |
| **2** | **26** | **86** | **0.513** |
| 3 | 7 | 24 | 0.833 |
| 4 | 1 | 5 | 0.900 |
| **pooled (comparison-weighted)** | | | **0.490** |

**E3 falls from 0.733 to 0.490 — chance.** The stratum holding 26 of the 36 actives sits
at 0.513. The 0.733 was artificial enrichment from aromatic count, exactly as the decoy
literature predicts, and the 9.0× noise-floor margin of §8.3 measured that confound
rather than a capability.

This is the controlled version of what script 14 reached for and got wrong: the fix is a
**matched benchmark**, not a score transform. It required no new computation — the
library screen had already scored every Tyr-starting library peptide at the same
configuration and seed as the benchmark run.

The one stratum with signal is 3 aromatics, AUROC 0.833 on 7 actives against 24 decoys.
That is thin, and it is also exactly where §7.4's 8 shortlisted peptides sit, so the two
results are at least internally consistent. It is the only surviving quantitative basis
for the shortlist.

### 9.3 TVERSKY_REF was the right choice, and the confound is its documented flip side

A large virtual-screening benchmark on pharmacophore-triplet descriptors found Tversky
weighted toward the **query** features "by far the most successful and the least
failure-prone" of ten dissimilarity scores, with Tanimoto failing to defend its
privileged status, and performance insensitive to tuning α. Tanimoto carries a
documented bias toward smaller compounds and stops decreasing as the second molecule
grows. So selecting TVERSKY_REF empirically landed on the literature's recommendation —
and the §7.2 confound is the known cost of the asymmetry, not a bug: ignoring excess
candidate features is exactly what makes it robust to size mismatch *and* what lets a
His-rich repeat accumulate coverage. Caveat: that benchmark used 2-D pharmacophore
fingerprints, not 3-D overlap.

### 9.4 Mean-over-conformers is non-standard but defensible

No head-to-head benchmark of max / mean / Boltzmann-weighted conformer aggregation turned
up. Max-over-conformers is described as standard practice for shape methods, so this
pipeline's `mean` is the minority choice. Supporting it: one study aggregating across
*pharmacophore models* (not conformers) found mean gave slightly better enrichment than
max, and this project's own S1 measurement is direct — ρ = 0.537 at `max` against 0.983
at `mean`. A separate finding that raising conformer counts increases recall but reduces
precision matches §4.2, where doubling the ensemble lowered E1 from 1.000 to 0.803.

The interpretive cost stands and is now sharper: `mean` asks how well a peptide's
accessible ensemble matches the template on average, not whether it can adopt a matching
pose. The second question is what binding requires; §8.2 shows this pipeline cannot
resolve it at any setting.

### 9.5 Revised conclusion

Three statements, in decreasing order of confidence:

1. **Uncontrolled, the screen enriches known opioid peptides from the library** —
   BEDROC 0.322–0.517 above a permutation null, EF(1%) 8–14×, best reference ligand 43rd
   of 14,595 (§9.1).
2. **That enrichment is not demonstrably more than aromatic-group content.** With
   aromatic count matched, E3 is 0.490 (§9.2). The screen has not been shown to add value
   over counting F/W/Y/H in the N-terminal window.
3. **It does not read 3-D arrangement at all** — epimer pairs, below the noise floor and
   in the wrong direction (§8.2).

The honest summary of this pipeline is that it is a reproducible, well-characterised
implementation whose discriminating power has not been separated from aromatic-residue
content, with one underpowered exception (3-aromatic stratum, AUROC 0.833, n = 7).

What would settle it: an aromatic-matched active set large enough to power the 2- and
3-aromatic strata. That is a data-curation task, not a compute task, and it is the only
remaining question worth spending on.
