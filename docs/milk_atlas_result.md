# The milk atlas screen: the positive controls do not clear the bar

**Run 2026-10-27** (`242501_[0-7]`, `242502`). The 34 motif-bearing peptides released from
milk proteins at a free N-terminal Tyr ([`milk_atlas.md`](milk_atlas.md)) scored against
both μ-opioid templates, each **with and without** exclusion volumes — four arms, two
conformer seeds, 800 conformers, window 8, `tversky_ref_mean`.

Predictions M1–M4 were written in `scripts/27_milk_atlas_analysis.py` and
[`milk_atlas.md`](milk_atlas.md) §5 before any atlas score existed. Machine-readable
results: [`milk_atlas_result.json`](milk_atlas_result.json).

## Result

| Arm | M1 known-vs-novel AUROC (≥ 0.70) | M3 ρ vs aromatic count | M4 seed ρ (≥ 0.95) |
|---|---|---|---|
| `a8bare` | 0.421  p=0.78  LOO 0.407–0.426 | **+0.805** | +0.961 pass |
| `a8excl` | 0.663  p=0.055  LOO 0.644–0.667 | **−0.012** | +0.970 pass |
| `b8bare` | 0.491  p=0.54  LOO 0.481–0.500 | **+0.805** | +0.964 pass |
| `b8excl` | 0.625  p=0.11  LOO 0.604–0.659 | **+0.310** | +0.980 pass |

| # | Prediction | Outcome |
|---|---|---|
| **M1** | AUROC(known, novel) ≥ 0.70 | **fail** in all four arms; best 0.663 |
| **M2** | the two aromatic-rich novels lose more under exclusion than the knowns | **pass**, both templates |
| **M3** | ρ vs aromatic count weaker with exclusion than without | **pass**, decisively |
| **M4** | seed ρ ≥ 0.95 | **pass**, 0.961–0.980 |

## M1 failed, and M1 was the one that licensed everything else

The pre-registration is unambiguous about what M1 is for:

> **M1** PRIMARY. […] This is the positive-control test, and it is the one that licenses
> reading anything into the novel scores at all.

The best arm reaches 0.663 against a required 0.70. **So the novel candidates from this
atlas are not licensed by this screen**, and no entry in the ranking below should be
carried forward as a candidate on the strength of its score.

That is the pre-registered consequence and it stands. Three qualifications belong with it,
none of which rescues the prediction:

**The `n` is inflated by pseudo-replication.** M1 was pre-registered as 5 known sites
vs 5 novel sites. It ran as n = 19/15, because each *site* enters as several length
variants of the same released peptide — β-casomorphin-7 contributes `YPFPGPI`,
`YPFPGPIPNS`, `YPFPGPIPNSLPQ` and `YPFPGPIPNSLPQNIPPLTQT`, which are one observation of
one site, not four. The permutation p of 0.055 and the narrow LOO range both treat them as
independent, so they are optimistic. The effective n is 10, and at n = 5 vs 5 this test
could not have resolved anything short of perfect separation.

**"Novel" does not mean inactive.** The negative class is *motif-bearing but never
assayed*, not *measured non-binder*. A novel peptide ranking high is not demonstrably a
false positive. This makes M1 a conservative test of the screen — but it also means M1
could not have been passed by a screen that was merely correct, so it is a weak instrument
in both directions. It was still the instrument chosen in advance.

**The top of the ranking is nonetheless all positive controls once exclusion volumes are
on.** `b8excl`'s entire top 8 are known opioid sites — β-casomorphin-7 and its extensions,
then the human β-casomorphin-5 region — and `a8excl`'s top 3 are the conserved
α-lactorphin site and β-casomorphin-7. The global AUROC is dragged down by the mid-table,
not by the knowns failing to surface. Retrieval at the top works; separation across the
whole set does not.

## M2 and M3: the exclusion volumes are confirmed on an independent set

This is the substantive positive result, and it replicates
[`dual_template_result.md`](dual_template_result.md) T6 on a different peptide set:

**Without exclusion volumes, the screen measures aromatic content.** ρ = **+0.805** for
both templates, and the bare arms' top five are the two aromatic-rich novel candidates —
`YPTYHTHGRYVPP` (4 aromatics) and `YPYAVWYYPQ` (5) — exactly the documented failure mode
from [`optimisation_study.md`](optimisation_study.md) §7.2, reproduced here unprompted.

**With exclusion volumes, it does not.** ρ falls to **−0.012** (template A) and **+0.310**
(template B), and the knowns displace the aromatic-rich novels at the top.

**M2 identifies those two candidates as confound artefacts, which is what it was for.**
`milk_atlas.md` §5 pre-committed to judging them by the with-vs-without comparison rather
than by rank. Under exclusion they lose **−0.106** and **−0.115** (template A) and
**−0.170** and **−0.186** (template B), against set medians of −0.084 and −0.136. They lose
more than the known controls under both templates. Per the rule written in advance, they
are **not believed**: they ranked first by over-filling the chamber with aromatic groups,
and the steric bound removes them.

That the comparison could be made at all is the reason the bare arms were added to the job
— the documentation had claimed this comparison while the job ran only the "with" arm.

## Verdict

The atlas itself remains sound and is unaffected by this: it recovered six known opioid
families unprompted, at their correct mature-protein positions, from motif definitions
alone. What failed is the **screen's** ability to rank the atlas, and with it the case for
carrying any novel site forward.

Nothing here measures binding, bioavailability, or survival of further digestion. Per
[`structural_basis.md`](structural_basis.md) §5 both templates are pan-opioid message
matchers, so none of this supports a μ-selectivity claim.

**What would change the answer** is not more compute on this set. M1 is limited by 10
sites, 5 of them positive. The fix is the same as D6.6 item 3: more assayed actives. A
screen cannot be validated against a positive-control class of five.
