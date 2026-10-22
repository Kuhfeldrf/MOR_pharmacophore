# Validation and optimisation study: pharmacophore matching against the μ-opioid receptor

**Written 2026-10-20, before the sweep is run.** This supersedes `PLAN.md` as the active
document. `PLAN.md` described a first implementation whose purpose was to demonstrate the
pipeline; this describes a study whose purpose is to establish whether the method works
and, if it does, with what settings. The first run's result
([`docs/screen_validation.md`](docs/screen_validation.md)) is the starting point, not
something to be defended.

Engineering and statistics decisions, with rationale and status: [`DECISIONS.md`](DECISIONS.md).

Scope: one receptor, OPRM1 (P35372). Establish a validated framework and optimised
settings, then scale to the full library.

---

## 1. What the first run established

| | |
|---|---|
| Works | Template retrieves its own ligand and, more importantly, retrieves endomorphin-1 at the ~99th percentile from **8F7R, a structure held out of template construction** |
| Works | The positional rule separates peptides starting with the opioid message from those merely spanning it, in both seeds |
| **Broken** | At 50 conformers the score is not reproducible enough to tier: 73% decile churn, median \|Δ\| 0.021 against a p10–p90 range of 0.134 |
| **Not shown** | Any seed-robust added value over the free "starts with Tyr" rule |

The reproducibility defect is blocking: an effect smaller than 0.021 cannot be resolved,
and the first run's apparent C3 signal disappeared on the second seed. Nothing else is
interpretable until that is fixed, so convergence is parameter one.

---

## 2. The benchmark

The first run's positives (21, mostly not in the library, pKi spanning 0.4 log units)
were too few and too narrow. This study uses the ChEMBL-derived affinity set already
curated in this project: **36 peptides with measured pKi spanning 4.71–9.14, i.e. 4.43
log units**, 31 Tyr-starting, lengths 4–17, 23 C-terminally amidated.

### 2.1 The potency-pair control — the centrepiece

Two families in that set differ by a **single Gly3→Pro substitution** and about three
orders of magnitude in affinity:

| Family | n | pKi range |
|---|---|---|
| `YGGFL…` | 11 | 7.39 – 8.71 |
| `YGPFL…` | 6 | 4.71 – 5.51 |

They are length-matched, nearly composition-matched, and every member of both families
begins with Tyr. **No sequence rule separates them.** Proline at position 3 restricts the
backbone, so the difference is conformational — exactly what a pharmacophore method
claims to read and what a sequence or composition method cannot.

This project has already measured that co-folding **fails** this contrast: iLIS scores
0.409, ranking the inactive series above the active one, while measured affinity
separates them at AUROC 1.000.

This is the hard test. It is the reason to run the study at all.

### 2.2 Decoy sets

| Set | Construction | Tests |
|---|---|---|
| `scramble` | one composition-matched permutation per active, fixed seed, identity excluded | basic discrimination; controls amino-acid content |
| `tyr_hard` | Tyr-starting library peptides, length-matched to the actives, no documented activity | added value over "starts with Tyr", which is constant across this set |
| `positional` | the 74 library peptides that span YPFP with a blocked N-terminus, vs the 4 that start with it | the positional rule, isolated from motif presence |

### 2.3 Honest handling of ties

Peptides sharing an N-terminal window are one molecule to a windowed screen. In the
affinity set, window 5 collapses 36 peptides into **20** distinct windows, and the
11-member `YGGFL` family into **one**. A windowed screen therefore cannot even attempt
§2.1. This tension between tractability and resolution is the central thing the sweep
must measure, not something to be assumed away.

---

## 3. Pre-registered separation of selection from evaluation

The risk in a parameter sweep is choosing settings on the same numbers that are then
reported as the result. To prevent that, the metrics are split **before any sweep runs**:

**Selection metrics — may be used to choose the configuration:**

- **S1 Reproducibility.** Seed-to-seed agreement: Spearman ρ between two independent
  conformer seeds, and median \|Δ\| relative to the spread of the scored set.
  *Gate: a configuration is ineligible unless ρ ≥ 0.95.*
- **S2 Classification.** AUROC, actives vs `scramble`.

**Evaluation metrics — reported for the chosen configuration, never used to choose it:**

- **E1 The potency-pair control** (§2.1), AUROC `YGGFL` family vs `YGPFL` family.
- **E2 Potency ranking.** Spearman ρ between score and pKi across all 36 actives.
- **E3 Added value over the sequence baseline.** AUROC, actives vs `tyr_hard`.
- **E4 Held-out structure retrieval.** Endomorphin-1 percentile, template from 8F7Q only.
- **E5 Positional rule.** AUROC, starts-with-YPFP vs spans-YPFP.

If the selected configuration fails E1, the method does not do the thing it claims to do,
and that is the finding. E2 is expected to be hard: this project has already shown
potency ranking failing on seven metrics with this much dynamic range, so a negative E2
is a confirmation, not a surprise, and only a positive E2 would be news.

---

## 4. The sweep

Conformer generation is the only expensive axis; scoring metric and ensemble statistic
are read from the same score table and cost nothing extra.

| Axis | Values | Why |
|---|---|---|
| `n_conformers` | 50, 200, 500 | the blocking defect; 50 is known to be too few |
| `window` | 5, 8, full peptide | window 5 cannot attempt E1 at all; full is the honest comparator |
| `metric` | TANIMOTO, TVERSKY_REF, TVERSKY_DB | free; TVERSKY_REF is asymmetric and expected to suit size-mismatched candidates |
| `ensemble statistic` | max, mean, top-10% mean | free; max is the noisiest order statistic and is the current default |

9 conformer configurations × ~190 benchmark peptides, each run twice with independent
seeds for S1. Scoring variants are expanded post hoc from the stored tables.

**Pre-registered tie-break:** among configurations passing the S1 gate, choose the
highest S2; ties broken by lower cost.

### 4.1 Amendment, 2026-10-20, made on S1/S2 only and before any evaluation metric was computed

Two facts from stage A forced a change, and both were visible from the selection
metrics alone:

1. **Full-peptide conformers are not viable.** A 17-mer has 81 rotatable bonds and
   default ETKDGv3 returns *zero* conformers for it; with random starting coordinates it
   returns 20 in 73 s, so 500 would cost ~30 min for one peptide. `window = full` is
   replaced by `window = 12`. More importantly, embedding silently under-delivers at
   window >= 8 without `useRandomCoords` (50 requested -> 21, 2, 1 at windows 8, 10, 12),
   which is now fixed and guarded.
2. **The selection rule picks a window under which the central test is undefined.**
   Ranked by S1 then S2, the best configuration uses window 5 — but window 5 collapses
   the 11-member `YGGFL` family and the 6-member `YGPFL` family to one point each, so E1
   becomes a single comparison with no error bars. Window 8 gives 30 family-pair
   comparisons and window 12 gives 54.

**Amended rule:** the selection metrics choose `n_conformers`, `metric` and `statistic`
*within each window*, and all three windows are carried into evaluation. The window is a
structural choice about which questions are answerable, not a tunable parameter, so it is
reported rather than optimised away. The globally best configuration under the original
rule is still reported as primary.

This amendment was written before any E metric had been computed, and the evaluation
metrics remain untouched by selection. Cost is recorded but is not a selection criterion
on its own, because the benchmark is small and the library screen is only ~52 CPU-hours
at the cheapest setting.

---

## 5. Scaling decision

The full library is run **only after** a configuration passes the S1 gate. The cost of
the chosen configuration is measured on the benchmark and extrapolated before submission,
as it was for the first run. If the chosen window is `full`, the cost model changes
substantially — the library's median length is 15, not 5 — and that estimate is made
before, not after, the array is submitted.

---

## 6. What would make this a negative result

Stated now so it cannot be softened later:

- No configuration reaches ρ ≥ 0.95 → the score is not reproducible at achievable cost,
  and the method is not usable for ranking regardless of any AUROC.
- The chosen configuration fails E1 → the method cannot separate a three-log-unit
  affinity difference that is conformational in origin, which is the specific capability
  it is supposed to have over sequence and composition methods.
- E3 at chance → the screen reduces to "starts with Tyr", as every structure-based tier
  tested in this project so far has.

Any of these is reported as the study's result.

---

## 7. Post-hoc addition, 2026-10-21: E6, the composition confound

**Added after the full library screen, from its result, and therefore post-hoc for this
receptor.** Pre-registered here for the next one.

None of E1–E5 could detect what the library ranking showed. Every decoy set above
controls composition by construction — `scramble` is composition-matched, `tyr_hard` is
Tyr- and length-matched — so all of them hold the confound fixed and are blind to it. On
the unmatched library, ρ between the score and the aromatic-residue count of the scored
window is **+0.909**, the raw tier 1 contains no known active, and 22 library peptides
outrank every endogenous opioid peptide
([`docs/optimisation_study.md`](docs/optimisation_study.md) §7).

- **E6 Composition confound.** Spearman ρ between score and the number of aromatic
  (F/W/Y/H) residues in the scored window, across the screened library.

### 7.1 Correction, same day: E6 is a reporting constraint, not a gate

E6 was first written here with a *ceiling of ρ = 0.6, above which a metric is
ineligible*. **That was wrong and is withdrawn.** Testing all nine metric × statistic
combinations against E6 and E1–E5 together
([`scripts/16_reselect_with_e6.py`](scripts/16_reselect_with_e6.py),
[`docs/reselection_e6.json`](docs/reselection_e6.json)) shows E6 and the evaluation
metrics are **collinear**: TVERSKY_REF/mean has E6 +0.913 and E1 0.803, TANIMOTO/top10
has E6 +0.730 and E1 0.515, TVERSKY_DB/mean has E6 +0.310 and E1 0.439. The only
combinations that clear a 0.6 ceiling sit at chance on the study's central test.

The mechanism explains why, and it is not a defect to be tuned away. TVERSKY_REF
measures how much of the *template's* pharmacophore a candidate covers, the template's
anchors are the Tyr1 and Phe4 aromatics of the opioid message, and a Gly3→Pro
substitution changes how well those anchors can be presented. "Tracks aromatic coverage"
and "separates the potency pair" are one property, not two that trade off. A gate that
eliminates every working configuration is a badly posed gate.

**E6 is therefore a constraint on how output is reported, not on which metric is used:**
above ρ = 0.6, raw cross-library tiers are not reported, and ranking is confined to the
composition strata the known actives occupy with the bar set by the weakest known active
in the same stratum. The score stays as selected; what changes is the comparison it is
allowed to make. A score can be a valid within-stratum ranking and an invalid
cross-stratum one at the same time, and E6 measures exactly that gap.

Full percentile stratification was tried as an alternative remedy and rejected: it failed
its own E1 invariance control and promoted peptides with no aromatic groups at all. The
rejection is kept in the record (`scripts/14_composition_control.py`) because why it
failed is the clearest available statement of what the confound is.

### 7.2 What can actually break the collinearity

No choice of metric can, because composition and shape are confounded in the *benchmark*,
not in the scoring function: of the 23 single- and double-substitution pairs with ≥1 log
pKi gap in the active set, only 9 are aromatic-matched and those collapse to about two
substitution types. The remaining 14 all gain or lose an aromatic residue, so E6 would
answer them for free.

**L/D epimer pairs are the limiting case and the way out.** An epimer pair has identical
atoms, identical aromatic count and identical everything a composition rule can see,
differing only in the handedness of one alpha carbon — so aromatic count is held exactly
fixed and any separation is necessarily conformational. `conformers.py` now builds them
(lowercase = D, so `YaFGYPS` is [D-Ala²]dermorphin), validated to invert the alpha carbon
S→R at identical molecular formula and skeleton, with D-Gly correctly refused.

**E7, pre-registered for the next receptor: epimer discrimination.** AUROC over
naturally-occurring D-residue opioid peptides against their L diastereomers. This is the
one test that can distinguish "reads the shape" from "counts the aromatics", and it is
the test that should have been the centrepiece from the start.

A matched-control benchmark measures the capability you care about and cannot see what
dominates an unmatched screen. Both are needed.

---

## 8. Metric corrections from the methodology literature, 2026-10-21

Everything above was decided on this project's own measurements. Checked against the
virtual-screening methodology literature, two of those choices are wrong and one is
confirmed. Recorded here as amendments; full evidence in
[`docs/optimisation_study.md`](docs/optimisation_study.md) §9.

**8.1 Screening contrasts are scored with BEDROC, not AUROC.** ROC AUC equals an average
of the actives' ranks, so early recognition is cancelled by late recognition (Truchon &
Bayly, *J Chem Inf Model* 2007, 47:488–508). BEDROC at α = 20 replaces it **for
screening-shaped questions only** — E3 and the library ranking. **E1 keeps AUROC**: it is
a two-sample comparison, where AUROC is the Mann-Whitney *U* and is correct. Because this
study has 36 actives against the ~50 recommended for stable BEDROC, a permutation null is
computed rather than assumed, and BEDROC is not reported for benchmark-sized contrasts at
all (N = 121 with α = 20 puts all weight in the top ~10 ranks).

**8.2 Decoys must be matched on the properties that drive the score.** `tyr_hard` matched
Tyr-start and length but not aromatic count, which is textbook **artificial enrichment**
(*Bioinformatics* 2021, 37:2134; *PLOS ONE* 2019 e0220113). Matched E3 is **0.490 against
0.733 unmatched**, so the original figure measured the confound. From now on a decoy set
must be matched on the pharmacophore feature inventory — not just length and the
N-terminal residue — and E6 is the check that it was.

**8.3 TVERSKY_REF is confirmed.** Tversky weighted toward query features outperformed
nine alternatives in a large pharmacophore-descriptor screening benchmark, with Tanimoto
carrying a documented small-molecule bias. The §7 confound is the known cost of the
asymmetry, not a defect. No change.

**8.4 The resulting pre-registration for the next receptor** is E1 (AUROC), E2, E3
(BEDROC + aromatic-matched decoys), E4, E5, E6 (reporting constraint), E7 (epimer sign
test) — with the decoy set matched on feature inventory before anything is run, and the
active set checked for sufficient actives per stratum to power E3 within strata. The
present study cannot power the 2-aromatic stratum, which is where 26 of its 36 actives
sit.

---

## 9. The steric bound: exclusion volumes, and E8 as the confound-resistant test

§7 and §8 end in a bad place. Every positive number in the study was readable as
aromatic-residue count, and matched-decoy E3 (0.490) said the method had no demonstrated
signal beyond composition. §9 is the resolution, and it came from the structure paper
rather than from more tuning.

**9.1 The defect was mis-stated.** Wang et al. (*Cell* 2023, 186:413–427) report that the
YGGF message occupies a "highly aromatic pocket" conserved across all four opioid
receptors, and that synthetic YGGF alone activates all four. Aromatic preference is
therefore **real chemistry**, not an artefact — our ρ = +0.909 was reading something true.
What was wrong is that the pocket is a chamber of *fixed size*, and our score had no upper
bound: it rewarded aromatic coverage without limit, which is how a His-rich repeat
(`EFHSHEFH`, 5 aromatics) outscored every endogenous ligand (2–3). The correct statement
is **the score captured the pocket's aromatic preference and omitted its steric bound**,
which is why no choice of similarity metric helped (§8.1 — every metric that weakened the
correlation weakened E1 proportionally) and why exclusion volumes are the right fix.

**9.2 Implemented and measured.** 86 receptor-derived exclusion spheres over 32 contact
residues. On the 14-peptide diagnostic the defect **inverts**: message peptides minus
His-rich windows goes −0.0444 → **+0.0677**, with per-peptide penalties ordered as the
structure requires (met-enkephalin, whose window *is* the template ligand, penalised least
at −0.074; His-rich repeats most, −0.19 to −0.21).

**9.3 E8, and why it is the test that should have been designed first.** The dietary
peptide set ([`docs/dietary_screen.md`](docs/dietary_screen.md)) contains 14 agonists and
3 antagonists with published IC50 values we did not curate. Its antagonists are **more**
aromatic than its agonists (2.67 vs 2.00), so the confound predicts the *opposite* of the
structural hypothesis. That makes it the first discriminating test in the study.

| | agonist-vs-antagonist AUROC (μ assays) | exact p |
|---|---|---|
| p5, no exclusion | 0.528 | 0.47 |
| p8, no exclusion | 0.444 | 0.63 |
| **p8 + exclusion** | **0.917** | **0.0154** |

Leave-one-antagonist-out gives 0.875 / 0.875 / 1.000, so it is not one peptide. ρ vs
aromatic count flips +0.149 → −0.261 and S1 rises to 0.994.

**9.4 Two negative results worth as much as the positive one.** (a) Extending the template
from 5 to 8 residues to capture the load-bearing E8–R213(ECL2) hydrogen bond makes E8
*worse* on its own (0.528 → 0.444): more reference features without a bound gives a
candidate more ways to accumulate overlap. A bigger template is not automatically a better
one. (b) **Potency ranking still fails** — E9 reaches only ρ = +0.16, and the two weakest
agonists in the set rank first and second while met-enkephalin, 250× more potent, ranks
third. The method discriminates *binding mode*, not affinity, and must be described that
way.

**9.5 Pre-registration amended.** E8 replaces E3 as the headline evaluation metric: E3 is
confound-limited by construction, whereas E8's confound is adversarial. E3 is retained as
a reported secondary. Two independent sets now support the steric bound; neither is the
pre-registered 235-peptide benchmark, which is what Stage C exists to settle.
