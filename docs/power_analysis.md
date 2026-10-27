# What the dietary set can and cannot resolve

**Run 2026-10-27**, `scripts/28_power_analysis.py`, seed 20261027. Machine-readable:
[`power_analysis.json`](power_analysis.json).

Every validation statistic in this project is an agonist-vs-antagonist AUROC on the
dietary set, and the antagonist class has **three** members. This document measures what
that buys. It exists because [`dual_template_result.md`](dual_template_result.md) was
decided by a comparison this set cannot make, and because the same limit applies to any
future template change.

## The evidence floor

The exact permutation test enumerates every assignment of labels, so the set's resolution
is fixed by combinatorics before any score is computed.

| Pool | Label assignments | Minimum one-sided p | AUROC granule |
|---|---|---|---|
| **12 ago × 3 ant (primary)** | **455** | **0.00220** | **0.0278** |
| 14 ago × 3 ant (secondary) | 680 | 0.00147 | 0.0238 |
| 12 ago × 5 ant | 6,188 | 0.00016 | 0.0167 |
| 12 ago × 8 ant | 125,970 | 0.00001 | 0.0104 |
| 12 ago × 12 ant | 2,704,156 | <0.00001 | 0.0069 |

No result on the primary pool can ever report p < 0.0022, however perfectly it separates
the classes. The observed values sit above that floor: `A8` p = 0.0088, `B8` p = 0.0352.
One discordant pair moves the AUROC by a full 0.028, which is the scale at which the
regenerated-runs discrepancy in `dual_template_result.md` lives.

## The interval, not the point estimate, is the result

| | AUROC | Hanley–McNeil SE | Percentile bootstrap (`23_dietary_screen.py`) | **Logit 95%** |
|---|---|---|---|---|
| `A8` | 0.944 | 0.061 | [0.78, 1.00] | **[0.635, 0.994]** |
| `B8` | 0.861 | 0.104 | [0.56, 1.00] | **[0.529, 0.972]** |

The symmetric Wald interval on `A8` is [0.825, **1.064**] — it breaks the bound, which is
the failure Nicholls (2014) works through explicitly and fixes with the logit transform
used here. The project's percentile bootstrap cannot exceed 1.0 but is correspondingly
optimistic, because at three antagonists many resamples separate perfectly; its own
docstring already calls it coarse. Both are reported above rather than one being dropped.

**`A8`'s AUROC is consistent with anything from 0.64 to 0.99.** That is the honest
reading, and it applies equally to this project's headline positive result — the
exclusion-volume effect — and to its negative one.

**Adding agonists cannot fix this.** The Hanley–McNeil variance goes as
1/n_pos + 1/n_neg, so the class of three dominates the uncertainty no matter how many
agonists are added. Only antagonists buy precision.

## Can a working template be detected at all?

Power of the one-sided permutation test to reject AUROC = 0.5 at α = 0.05, 12 agonists
fixed, 600 simulations per cell:

| n antagonists | true AUROC 0.75 | 0.85 | 0.95 |
|---|---|---|---|
| **3 (current)** | **0.39** | **0.62** | 0.93 |
| 5 | 0.48 | 0.74 | 0.98 |
| 8 | 0.57 | 0.90 | 1.00 |
| 12 | 0.70 | 0.95 | 1.00 |

At three antagonists a genuinely good template (true 0.85) is missed nearly two times in
five, and a merely decent one (0.75) is missed more often than it is found. A negative
result on this set is therefore weak evidence of absence.

## Can two templates be told apart? No.

The observed `A8` − `B8` gap is +0.083, three discordant pairs. A paired bootstrap over
peptides (20,000 resamples, agonists and antagonists resampled independently) gives:

> **Δ AUROC = +0.083, 95% CI [0.000, +0.333], P(Δ ≤ 0) = 0.316.**

Nicholls (2016) prescribes exactly this — for dependent systems the error bar belongs on
the difference per system — so the non-result is not an artefact of the wrong test.

Power to resolve a **true** gap between two templates whose score vectors correlate at
ρ = 0.9, which is the regime T3 measured:

| n antagonists | true gap 0.05 | 0.10 | 0.15 |
|---|---|---|---|
| **3 (current)** | **0.01** | **0.06** | **0.14** |
| 5 | 0.03 | 0.19 | 0.32 |
| 8 | 0.08 | 0.29 | 0.56 |
| 12 | 0.14 | 0.48 | 0.74 |
| 20 | 0.22 | 0.65 | 0.90 |

**This is the result that matters for design.** At three antagonists, the chance of
resolving even a large 0.15 AUROC gap between two correlated templates is 0.14. The
dual-template design did not fail because endomorphin-1 was the wrong second template; it
failed because **no template comparison is adjudicable on this set**. T1 and T2 were
decided by a class of three.

It also explains which pre-registered tests returned usable answers. T3 (ρ between
templates) and T6 (ρ vs aromatic count) are computed over all 20 peptides and were
decisive. T1 and T2 are gated on 3 antagonists and were not. The tests that worked are
the ones not bottlenecked by the minority class.

## What the literature says, and what it does not

- **Hanley & McNeil 1982**, *Radiology* 143:29, [10.1148/radiology.143.1.7063747](https://doi.org/10.1148/radiology.143.1.7063747) — the SE formula; variance scales as 1/n_pos + 1/n_neg.
- **Nicholls 2014**, *JCAMD* 28:887, [10.1007/s10822-014-9753-z](https://doi.org/10.1007/s10822-014-9753-z) — use DeLong with primary data; effective degrees of freedom ≈ n_small − 1; "about twenty samples before the inaccuracy in the prefactor of 1.96 is less than 10%"; works an AUC = 0.9 / ten-actives example where the symmetric upper limit exceeds one, and prescribes the logit transform.
- **Nicholls 2016**, *JCAMD* 30:103, [10.1007/s10822-016-9904-5](https://doi.org/10.1007/s10822-016-9904-5) — "the correct way to construct an error bar for a dependent system is on the difference between A and B per system"; extends DeLong to correlated AUCs. Does **not** discuss permutation tests, so the exact p used here is defensible but not prescribed by this literature.
- **Feng et al. 2017**, *Stat Methods Med Res* 26:2603, [10.1177/0962280215602040](https://doi.org/10.1177/0962280215602040) — 29 CI methods; the discrepancy between them grows as the true AUC rises and n falls, i.e. exactly this regime. Reason to state which method produced an interval.
- **Kottas et al. 2014**, *BMC Med Res Methodol* 14:26, [10.1186/1471-2288-14-26](https://doi.org/10.1186/1471-2288-14-26) — for small samples recommends the continuity-corrected Wald, with coverage comparable to logit.
- **Do not swap AUROC for an early-recognition metric to compensate.** Zhao et al. 2009, *BMC Bioinformatics* 10:225, [10.1186/1471-2105-10-225](https://doi.org/10.1186/1471-2105-10-225), report the "seesaw effect": a metric that overemphasises early recognition loses statistical power as the active set shrinks. Empereur-Mot et al. 2015, *J Cheminform* 7:52, [10.1186/s13321-015-0100-8](https://doi.org/10.1186/s13321-015-0100-8), find ROC AUC, BEDROC, AUAC and average-rank closely related despite appearing distinct. Enrichment factor is not even computable here: with 15 compounds EF@1% is undefined and EF@10% is one compound.
- **There is no minimum-n convention.** No paper in cheminformatics or biostatistics was found stating a hard minimum number of actives or inactives for a usable AUC CI. The defensible statement is that the CI must be reported, not that some n is sufficient.
- **There is no conventional model-redundancy correlation threshold.** T3's ρ ≥ 0.9 has no precedent in the VS literature: inter-model score correlation is usually not measured at all. The nearest conceptual standard is Yang et al. 2005's "distinctiveness" precondition ([10.1021/ci050034w](https://doi.org/10.1021/ci050034w)), which attaches no threshold; the nearest operational control is geometric deduplication (Polishchuk et al. 2019, [10.3390/ijms20235834](https://doi.org/10.3390/ijms20235834)). The |r| > 0.9 convention exists for **descriptor** collinearity in QSAR, which is an analogy, not a precedent. T3 should be presented as a pre-registered choice made in the absence of convention — which is what it was.
- **The consensus mechanism requires independence.** Rank/score fusion averages errors across estimates assumed independent (Wang & Wang 2001, [10.1021/ci010025x](https://doi.org/10.1021/ci010025x)); at ρ = 0.917 that mechanism does not apply, so T2's failure is mechanistically expected rather than surprising. Even in ensemble ML, where diversity measures are formalised, their prospective usefulness is contested (Kuncheva & Whitaker 2003, [10.1023/A:1022859003006](https://doi.org/10.1023/A:1022859003006)).

**Provenance note.** Truchon & Bayly 2007, Zhao et al. 2009 and Yang et al. 2005 were
read at abstract level only — publisher full text was inaccessible — and their internals
above come from secondary sources that quote them. The DOIs were resolved; the quoted
internals should be checked against the full text before being relied on in a manuscript.

## Consequence for the study

The evaluation set, not the template, is the binding constraint. Three things follow:

1. **Report the interval, not the point estimate**, on every AUROC in this project, and
   report which method produced it.
2. **Any further template work is unvalidatable until the antagonist class grows.** Not
   "weakly validated" — unvalidatable: power 0.01–0.14 against the gaps at issue.
3. **Only antagonists help.** More agonists, more seeds, more conformers and better
   templates all leave 1/n_neg untouched.
