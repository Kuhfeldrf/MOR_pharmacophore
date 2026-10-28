# Two templates, one measurement: the dual-template design does not hold

**Run 2026-10-27.** The 20-peptide dietary set scored against two experimental μ-opioid
templates — **A** = 8F7Q/β-endorphin(1–31), extent 8 (`p8excl`), and **B** =
8F7R/endomorphin-1, complete ligand (`b8excl`) — both with α-1 exclusion volumes, 800
conformers, window 8, `tversky_ref_mean`, two independent conformer seeds. Tests T1–T6
were written in [`dual_template_design.md`](dual_template_design.md) §6 before any of
these scores existed.

Machine-readable results: [`dual_template_result.json`](dual_template_result.json).
Single-template numbers for both templates: [`dietary_screen.json`](dietary_screen.json).

## Result

Seed-averaged on the primary pool (`E8_mu_assays_only`, 12 agonists / 3 antagonists), as
§6 requires, with per-seed values reported because a test holds only if it also fails to
reverse in either seed.

| # | Test | Prediction, stated in advance | Seed-averaged | Seed 1 / Seed 2 | |
|---|---|---|---|---|---|
| **T1** | `B8` AUROC ≥ 0.75 **and** ≥ `A8` | B beats A | **B 0.861** vs **A 0.944** | 0.861 / 0.917 vs 0.944 / 0.944 | **fail** |
| **T2** | Consensus ≥ max(`A8`, `B8`) | consensus wins | **0.861** max-norm, **0.903** rank-fusion vs max 0.944 | 0.861 / 0.917; 0.903 / 0.944 | **fail** |
| **T3** | ρ(`A8`, `B8`) < 0.9 | templates see different things | **+0.917** | +0.926 / +0.912 | **fail** |
| **T4** | `B8` ρ vs pIC50 in (0, 0.6) | weak positive | **+0.061** | +0.036 / +0.155 | pass |
| **T5** | DAMGO ≥ 90th pct of library | top decile for ≥ 1 template | **not computed** | — | — |
| **T6** | `B8` ρ vs aromatic count < 0.909 | B less composition-driven | **−0.116** | −0.087 / −0.173 | pass |

Every AUROC above has a leave-one-antagonist-out range entirely above 0.5 (T1 B8:
0.792, 0.792, 1.000; T1 A8: 0.917, 0.917, 1.000), so none is withheld under the §6 rule.
`B8` on the secondary pool (`published, supported termini`) is the exception worth
naming: AUROC 0.738, exact one-sided p = 0.122 — not significant, and below `A8`'s 0.857
on the same pool.

## What this means, stated plainly

**The dual-template design is unjustified.** §9, first bullet, written in advance:

> If **T3** shows ρ ≥ 0.9 between the two templates' score vectors, the second screen is
> measuring the same thing as the first and the dual-template design is unjustified
> regardless of how T1 comes out.

ρ = +0.917. Template B re-measures what template A already measured. The consensus of §7
is therefore not a consensus of two views, and `E4′` would have been a percentile against
a distribution nobody should quote.

**The ~52 CPU-hour library screen was not run.** `scripts/26_t3_gate.py` refused it
(`GATE FAIL: |rho| = 0.917 >= 0.9`) and the chain reported
`LIBRARY SCREEN NOT SUBMITTED - T3 gate failed`. This is the gate working as designed:
the cheap screen came first precisely so this could be decided for the price of a few
minutes, by a number rather than by whoever was at the terminal.

**The conclusion does not rest on where the T3 threshold was drawn.** ρ = 0.917 is close
to 0.9, and the margin would be easy to argue. It does not need arguing: T1 and T2 fail
independently of T3. `B8` did not reach `A8` on either seed or either pool, and neither
consensus method beat the better single template — rank fusion only tied it, on one seed.
At a threshold of 0.95 the design still fails. The complete-ligand advantage predicted in
§3 did not appear.

**But "A beats B" is not a supported claim either, and T1 was never equipped to make
one.** T1 was pre-registered as a comparison of point estimates with no test attached.
A paired bootstrap over peptides puts the observed Δ AUROC of +0.083 at a 95% CI of
**[0.000, +0.333]**, with 32% of resamples placing A at or below B. The honest reading is
that `B8` failed to beat `A8` as predicted, **not** that `A8` is the better template. See
[`power_analysis.md`](power_analysis.md): on 12 agonists × 3 antagonists, the probability
of resolving even a 0.15 AUROC gap between two templates correlated at ρ = 0.9 is 0.11.
This set cannot adjudicate templates, and T1 and T2 were decided by a class of three.

## What passed, and it is worth stating

**T6 is a real positive.** The aromatic-composition confound that once let a His-rich
repeat outscore every endogenous ligand at ρ = +0.909
([`optimisation_study.md`](optimisation_study.md) §7.2) is now *negative* for both
templates: −0.231 for A, −0.116 for B. The exclusion volumes (D2.4) did what they were
added to do, and this is the clearest evidence yet that they work. T4 passed as the
confirmation it was pre-registered as, not as a discovery.

## A reproducibility note on `A8`

§6 pre-registered `A8` as **0.917**, its published value. It is reported above as
**0.944**, and the reason is recorded here rather than left as a discrepancy.

Template A's dietary runs were lost when the scratch allocation was reaped and were
regenerated for this analysis. Most of the set reproduces exactly: 14 of 20 peptides for
seed 20261020 and 13 of 20 for seed 61453 are identical to four decimal places, and mean
|Δ| is 0.0005 and 0.0028 against a documented seed jitter of 0.00133–0.00225.

Two peptides in seed 61453 did not. β-lactorphin moved **+0.0277** and casoxin A
**+0.0213**, roughly ten times jitter, with conformer counts unchanged (793→793,
800→800) — so this is not an ensemble-size effect. Either way it is a ~0.02 effect on 2 of
20 peptides, and it is recorded here rather than averaged away.

**Resolved 2026-10-28 as far as the evidence allows; see
[`determinism_check.md`](determinism_check.md).** Non-determinism is ruled out: conformer
generation and align-it's alignment search are both bit-for-bit identical on repeated runs
over these same peptides, and this run reproduces the regenerated production scores to
full precision two days later (β-lactorphin 0.190446 against 0.1904, casoxin A 0.134792
against 0.1348). The surviving explanation is that the original and regenerated runs used
different library builds — the environment was rebuilt from scratch between them — which
**cannot be confirmed**, because the original run recorded only align-it's version and not
RDKit's. The numbers are sound within an environment and must not be compared across a
rebuild. `04_screen.py` now records the full library stack so this is answerable next time.

That is enough to move the AUROC, because on 12 agonists × 3 antagonists the statistic is
granular: **33/36 → 34/36, one discordant pair.** A single shift of this size flips one
pair and moves the headline by a full 0.028. That granule is the right scale against which
to read every AUROC in this document — including A's lead over B, 3/36 = 0.083, three
times it. No verdict in the table changes under either value of `A8`: `B8` at 0.861 sits
below A at 0.917 and at 0.944 alike.

**Regenerating rather than restoring these runs removed a confound.** The surviving
copies were produced before `a89b9a2`; template B's runs were produced after it. Scoring
A-from-the-old-files against B-from-the-new-code would have made T1 and T3 — the two
tests this document turns on — a comparison across two code versions as well as two
templates. Both arms are now post-fix.

## T5 was not computed, and cannot be salvaged from what survives

`E4′` needs the library scored against the same template. The raw library scores
(`results/library_w8`) went with the reaped allocation; what survives is the committed
[`library_bare_ranked.tsv`](library_bare_ranked.tsv), 14,571 peptides at the same window
8, 800 conformers and `TVERSKY_REF`.

It is tempting and it is wrong. `jobs/library_screen.sbatch` scored against the **bare**
`mor_template.phar`, while every holdout run here used `mor_template_excl_a1.phar`.
Those are different distributions, so a DAMGO percentile taken across them would be
exactly the error D4.2 warns about. T5 stays open, and it requires new compute that the
T3 gate has just refused to authorise.

## What this does not claim

The ceiling in [`structural_basis.md`](structural_basis.md) §5 stands. Template B is
still a message-region template with no ECL2 contact — endomorphin-1 is four residues and
does not reach ECL2 — and nothing here tests δ/κ/NOPR, so no result in this document
licenses a μ-selectivity claim. The failure reported here is a failure of the
*dual-template* design, not evidence against template A, which remains the better of the
two on every test that distinguished them.
