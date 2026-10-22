# Dietary opioid peptide screen (E8, E9)

**Run.** SLURM 241252, 6 configurations (3 templates × 2 seeds), 20 peptides,
800 conformers, 8-residue window. Analysis: `scripts/23_dietary_screen.py`.
Machine-readable output: [`dietary_screen.json`](dietary_screen.json).

This is the only evaluation in the study whose labels we did not curate. The potencies
are published IC50 values from the primary pharmacology literature, the peptides are
dietary — the actual domain of the digestion library — and the set contains both
agonists and antagonists, which the ChEMBL benchmark does not.

---

## 1. Why E8 is the test that matters

Every earlier positive result in this project turned out to be readable as aromatic-residue
count ([`optimisation_study.md`](optimisation_study.md) §7.2; E3 collapsed 0.733 → 0.490
under aromatic-matched decoys). E8 is immune to that, because **the confound points the
wrong way**:

| | mean aromatic residues in the scored window |
|---|---|
| agonists (n=14) | **2.00** |
| antagonists (n=3) | **2.67** |

A score that merely counted aromatic groups would rank the antagonists **higher**. So an
agonist-over-antagonist result here cannot be produced by the confound — it has to come
from geometry. That was stated in the script before the scores existed, and is recorded in
the JSON under `predictions_registered_before_scoring`.

## 2. Result

Pre-registered primary analysis: published peptides, literature-supported termini,
μ-dominant assays only. Agonists vs antagonists, AUROC.

| template | AUROC | exact p (1-sided) | leave-one-antagonist-out | margin / measured jitter | E6 ρ vs aromatic count | S1 |
|---|---|---|---|---|---|---|
| p5 (message only) | 0.528 | 0.47 | 0.458 / 0.500 / 0.625 | +2.5× | +0.149 | 0.962 |
| p8 (adds E8–R213 ECL2) | 0.444 | 0.63 | 0.333 / 0.417 / 0.583 | +0.2× | +0.149 | 0.942 |
| **p8 + exclusion α=1.0** | **0.917** | **0.0154** | **0.875 / 0.875 / 1.000** | **+30.8×** | **−0.261** | **0.994** |

Secondary, all published peptides regardless of assay bias (n=14/3): p5 0.452, p8 0.381,
**p8excl 0.833 (exact p = 0.046)**.

Three things in that table are worth stating separately.

**Exclusion volumes take E8 from chance to 0.917.** The two configurations without a
steric bound are indistinguishable from a coin flip (p = 0.47, p = 0.63), and p8 base is
*below* chance. This is the first result in the study that the aromatic confound cannot
explain, and it is the direct confirmation of the mechanism argued from Wang et al. in
[`structural_basis.md`](structural_basis.md) §4: the pocket is a bounded aromatic chamber,
and our score was missing the bound, not the preference.

**It is not carried by one peptide.** With 3 antagonists that was the obvious failure
mode. Dropping each antagonist in turn leaves AUROC at 0.875, 0.875 and 1.000 — every
subset still well above chance. The bootstrap CI [0.67, 1.00] is reported but is coarse at
this n; the exact permutation p is the honest statistic.

**The confound sign flips.** ρ(score, aromatic count) goes +0.149 → **−0.261**. The score
stops rewarding aromatic content without limit, which is precisely the intended effect and
not something that was tuned for — α=1.0 was chosen on the 14-peptide diagnostic, before
this set was scored.

### Extending the template is not a substitute for bounding it

p8 adds the E8–R213(ECL2) hydrogen bond that Wang et al. show is load-bearing, and on its
own it makes E8 **worse** (0.528 → 0.444). More reference features without a steric bound
gives a candidate more ways to accumulate overlap. Only the two together work. This is a
point against the intuition that a larger template is automatically a better one.

## 3. What this does not show: E9 still fails

Spearman ρ between score and pIC50, μ-dominant assays:

| template | all published | μ-only | μ agonists only |
|---|---|---|---|
| p5 | +0.081 | −0.025 | +0.119 |
| p8 | −0.002 | −0.111 | −0.056 |
| p8excl | +0.221 | +0.164 | +0.140 |

Exclusion volumes improve this too, but **+0.16 over 15 peptides is not a potency
model**. The prediction registered beforehand was "weak positive at best", and that is
what happened.

The ranking makes the failure concrete. Under p8excl the top two are α-lactorphin
(IC50 50 µM) and β-lactorphin (160 µM) — the two *weakest* agonists in the set — while
met-enkephalin, more potent by a factor of 250, ranks third. The score separates
**agonist from antagonist** but carries almost no information about *how strong* an
agonist is.

That is a coherent and defensible position for a shape-complementarity score to be in:
binding mode is geometric, and potency depends on affinity, efficacy, assay tissue and
proteolytic stability, none of which this method represents. It should be stated that way
rather than presented as a potency predictor.

Caveats applied in code, not just noted: the two αs1-casein exorphins were measured on MVD
(δ-dominant) and are excluded from the μ-restricted analysis; casoxin-6 and lactoferroxin
A are methyl esters in their literature form and are scored as free acids and flagged;
the bovine lactoferrin f319–324 fragment has no published IC50. IC50 values span three
log units across different tissues and decades and are an ordering, not a calibration.

## 4. Mechanism: which peptides moved

The effect is concentrated where the structural argument predicts it. Casoxin A
(`YPSYGLNY`, 8 residues, 3 aromatics) is the clearest case — rank **3 → 17** of 20 when
the steric bound is added. It is the longest and most aromatic antagonist, exactly the
profile that accumulated unbounded overlap before. Casoxin C (`YIPIQYVLSR`, truncated to
the window) goes 10 → 18. Met-enkephalin, whose own window is the template ligand, moves
*up*, 8 → 3.

Conversely the peptides that rise are the short `Y..F` agonists — α- and β-lactorphin
(`YGLF`, `YLLF`) and the β-casomorphins — which fit a bounded chamber.

## 5. Consequences

1. **The dietary set is the better-matched evaluation**, as predicted in
   [`structural_basis.md`](structural_basis.md) §6: most of these peptides are `YP…`
   rather than `YGGF…`. It should be reported as a primary result, not an appendix.
2. **E8 becomes the headline evaluation metric**, displacing E3. E3 is confound-limited by
   construction; E8's confound points the other way.
3. **D2.4 (exclusion volumes) should move to HELD** once Stage C confirms it on the
   pre-registered 235-peptide benchmark. Two independent sets now agree — the 14-peptide
   diagnostic and these 20 — but neither is the pre-registered one.
4. **Any claim about potency ranking must be dropped.** The score discriminates binding
   mode, not affinity (§3).
5. **The library shortlist must be rebuilt under exclusion volumes.** It was generated
   under the unbounded score, which §2 and §4 show systematically favours over-aromatic,
   over-long windows — the same profile as casoxin A.

### Open: the endomorphin binding mode

Most of this set is `YP…`, and endomorphin-1 is `YPWF`. Wang et al. show endomorphin binds
β-endorphin's pocket *differently* — its W3 overlaps β-endorphin's F4, and its F4 occupies
an extended TM1/2/7 subpocket β-endorphin does not use. So these peptides may be better
matched to **8F7R** than to our 8F7Q template, and the E8 = 0.917 result may be an
*underestimate* of what the right template would give. Testing that costs the E4
held-out-structure validation, which is a deliberate trade-off and still
[open for decision](../DECISIONS.md).
