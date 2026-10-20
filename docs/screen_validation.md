# Pharmacophore screen against the μ-opioid receptor: first run

**Run 2026-10-20.** 14,595 targets (14,571 digestion peptides + 24 reference ligands)
scored against a pharmacophore template read from the experimental μ-opioid receptor
complex 8F7Q. Criteria C1–C5 were written in [`PLAN.md`](../PLAN.md) §3 before the run.

Machine-readable results: [`screen_validation.json`](screen_validation.json).

## Result

The screen was run twice, with independent conformer seeds. **Both runs are reported.**
The second seed is not a replicate for averaging — it is the C5 stability check, and it
changes how two of the other criteria read.

| # | Criterion | Threshold | Seed 1 | Seed 2 | |
|---|---|---|---|---|---|
| **C1** | Template recovers its own ligand (met-enkephalin) | top decile | 0.458, **93.1st** pct | 0.469, **94.8th** pct | **pass** |
| **C2** | Template recovers an independent ligand (endomorphin-1, from held-out 8F7R) | top decile | 0.503, **98.8th** pct | 0.549, **99.9th** pct | **pass** |
| **C3** | **Added value over the free sequence baseline** | AUROC ≥ 0.70, CI excluding 0.5 | 0.669, CI [0.528, 0.799] | **0.603, CI [0.484, 0.714]** | **fail** |
| **C4** | Positional rule is load-bearing | median difference > 0 | AUROC 0.973, Δmedian 0.061 | **AUROC 0.716**, Δmedian 0.055 | **pass, weakly** |
| **C5** | Conformer stability across seeds | ≤10% decile churn | **73.0% churn** | | **fail** |

Scores are Align-it TVERSKY_REF over a 50-conformer ensemble per peptide, best
conformer taken, N-terminal 5-residue window (PLAN.md D3).

## What passed, and it is worth stating plainly

**The template works as a template.** It retrieves the ligand it was built from (C1) and,
more importantly, retrieves **endomorphin-1 at the 98.8th percentile** (C2) — a peptide
whose structure, 8F7R, was deliberately held out of template construction. A template
built from one experimental pose generalising to an independently determined ligand is a
real result, not plumbing.

**The positional rule is load-bearing, though more weakly than one seed suggested
(C4).** Of the 78 library peptides containing the β-casomorphin message YPFP, only 4
begin with it; the other 74 span it with the Tyr α-amine locked in a peptide bond. The
screen separates the two groups in both seeds — identical motif, identical parent
protein, only the position differing, and no sequence-similarity method separates that
pair. But the effect size is seed-sensitive: **AUROC 0.973 on seed 1 and 0.716 on
seed 2**, on only 4 positives. The direction is reproducible; the magnitude is not, and
the seed-1 figure should not be quoted on its own.

**The global separation is in the expected direction.** Tyr-starting library peptides
have median 0.475 against 0.358 for the rest.

## What failed: C3, the criterion that matters

**AUROC 0.669 against a pre-registered threshold of 0.70, CI [0.528, 0.799].**

C3 asks the only question that distinguishes this screen from a one-line rule: *within
the Tyr-starting subset, where "starts with Tyr" is constant by construction, does 3D
pharmacophore matching separate documented actives from the rest?* It does not, at the
bar set in advance.

On seed 1 alone this looked like weak evidence of added value, because the interval
excluded 0.5 (lower bound 0.528). **The second seed removes that reading**: AUROC falls
to 0.603 with CI [0.484, 0.714], which includes 0.5. There is no seed-robust evidence
that the 3D score adds discrimination beyond the sequence rule.

The test is also badly powered: 21 positives against 639 negatives that are only **219
distinct N-terminal windows**, so ties are pervasive and the effective n is far below the
nominal one.

**The criterion stands as written.** C3 fails, in both seeds, and the apparent signal in
the first was not stable. The bar was set before the data and is not being moved now.
Had only one seed been run, this run would have been reported as a borderline near-miss;
it is not one. **That is what C5 was pre-registered for.**

## C5: the scores are not stable enough to tier

**73.0% of peptides change decile between the two seeds**, against a pre-registered
ceiling of 10%. C5 fails, and it fails for a reason that matters more than the number.

The per-peptide score barely moves: median |Δ| between seeds is **0.021**, max 0.174.
The problem is that the library's score distribution is extremely compressed — the whole
interdecile range (p10 0.310 to p90 0.444) is **0.134**. A 0.021 jitter is about one
sixth of that range, so it reshuffles decile membership wholesale even though no
peptide's score meaningfully changed.

Two consequences:

1. **Decile tiering is not a usable construct here.** Tiers built on this distribution
   describe the conformer seed as much as the peptide. Any tiering must be defined on a
   scale coarser than the seed jitter, or the ensemble must be converged much further.
2. **It bounds what C3 could have shown.** C4's effect (Δmedian ≈ 0.055–0.061) is roughly
   3× the seed jitter and survives; C3's separation is weaker than that and does not.
   A score whose reproducibility is 0.021 cannot resolve differences near 0.02.

The direct fix is more conformers. 50 per 5-residue window was chosen before measurement
and is evidently too few; the ensemble spread (median 0.144) is itself larger than the
interdecile range of the library, which says the best-of-ensemble statistic is still
sampling-limited. Re-running at a higher conformer count, with convergence checked rather
than assumed, is the first item on the list below.

## A false-positive cluster worth recording

The second through twelfth ranked library peptides are a single family sharing the
N-terminal window `EGIHA` (scoring 0.5597) — **no N-terminal Tyr, no opioid message**.
They outrank met-enkephalin.

This is informative in both directions. It shows the score is not merely a Tyr detector,
which is what C3 was built to probe. It also shows the template's positive-charge and
donor/acceptor features can be satisfied by a chemically unrelated N-terminus. The
free N-terminal amine supplies the POSC feature for *every* peptide, so that feature
carries no discriminating information on its own — only its geometric relationship to
the aromatic and hydrophobic features does.

## Known limits of this run

- **The window collapses the library.** 14,595 targets are 5,183 distinct N-terminal
  windows; peptides differing only beyond residue 5 are one molecule to this screen and
  score identically. This is defensible (D3) but it caps resolution and power.
- **Conformers are generated free, not bound.** Met-enkephalin scores 0.46 against a
  template derived from its own bound pose, because ETKDGv3 ensembles do not reproduce
  the bound conformation. Every score is therefore a lower bound on achievable overlap.
- **PLIP finds no salt bridge in the deposited complex.** The peptide N-terminus sits
  **4.55 Å** from Asp149 OD2 in 8F7Q (3.22 Å cryo-EM), too long and too ambiguous for a
  salt-bridge call, where this project's own MD tightened the same contact to 2.6–3.0 Å.
  The template's POSC feature comes from protonation at pH 7.4 (D4), not from a PLIP
  salt-bridge assignment.
- **21 positives is too few**, and most are not in the digestion library at all — the
  enkephalins and dynorphins come from proenkephalin and prodynorphin, not milk.
  β-casomorphin-7 is the one documented food-derived active that is genuinely a library
  member.
- **Three documented actives do not start with Tyr** (α_s1-casein exorphins 6 and 7,
  casoxin D), so they are excluded from C3's Tyr-restricted comparison by construction
  while remaining real counterexamples to the sequence baseline.

## What this does and does not support

It supports: the pipeline runs end to end on ORCA, CPU-only, at measured cost; the
template generalises to a held-out experimental ligand (C2, ~99th percentile in both
seeds, the most robust result here); the positional rule is real and measurable in
direction.

It does not support: using this score as a triage tier, or quoting any single-seed
effect size from this run. On the evidence here the screen
has not been shown to add enough to the free sequence rule to rank candidates, which is
the same conclusion the co-folding and interface-rescoring tiers reached in
`MOR_peptide_binding`, reached again by a different route and with the bar set in
advance.

## Cost

| | |
|---|---|
| Benchmark | 100 peptides, 8 cores, 2 min 39 s (123 s conformers, 30 s scoring) |
| Full screen | 146 array tasks × 8 cores, **≈52 CPU-hours**, ~15 min wall at 48 concurrent |
| Run twice (two seeds) | ≈104 CPU-hours total |
| GPU use | none |

The narrative's budget for motif matching, fingerprinting and model fitting is under
800 CPU-hours. One receptor costs ≈52, which is the first measured number to check that
figure against.

## Next, in order

1. **Converge the conformer ensemble.** C5 is the blocking defect: at 50 conformers the
   score's reproducibility (0.021) is comparable to the effects being measured. Raise the
   count until seed-to-seed churn meets the 10% bar, and check convergence rather than
   assuming it. Nothing else in this list is interpretable until this is fixed.
2. **Power C3 properly.** Expand the positive set beyond 21 — the binding constraint is
   documented actives, not compute. 52 CPU-hours per full screen means compute is not
   what limits this.
3. **Replace decile tiering.** Define tiers on a scale coarser than the seed jitter, from
   the negative distribution, as PLAN.md §4.5 intended.
4. **Weight the template by PLIP evidence.** Features engaging no receptor residue
   currently count as much as those that do; the PLIP interaction set is already
   recorded alongside the template and is not yet used to weight it. The POSC feature in
   particular is present on every peptide and so carries no discriminating information by
   itself.
5. **Test the whole-peptide alternative to D3** on the Tyr-starting subset, where it is
   affordable, to confirm the window is not itself the limitation.
