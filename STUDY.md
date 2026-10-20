# Validation and optimisation study: pharmacophore matching against the μ-opioid receptor

**Written 2026-10-20, before the sweep is run.** This supersedes `PLAN.md` as the active
document. `PLAN.md` described a first implementation whose purpose was to demonstrate the
pipeline; this describes a study whose purpose is to establish whether the method works
and, if it does, with what settings. The first run's result
([`docs/screen_validation.md`](docs/screen_validation.md)) is the starting point, not
something to be defended.

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
highest S2; ties broken by lower cost. Cost is recorded but is not a selection criterion
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
