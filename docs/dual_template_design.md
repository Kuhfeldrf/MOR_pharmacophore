# Dual-template screen: design and pre-registration

**Written 2026-10-26, before any score against the new template was computed.** Every
prediction in §6 is recorded here first so that it cannot be adjusted after the fact. The
one quantity in this document that *was* computed before writing it is the aromatic
composition of the dietary set per window (§5.2) — it is an input to the design, not a
result, and it is what forced the window choice.

Companion documents: [`structural_basis.md`](structural_basis.md) (what the source paper
says), [`dietary_screen.md`](dietary_screen.md) (the single-template result this extends),
[`../DECISIONS.md`](../DECISIONS.md) (engineering record), [`../STUDY.md`](../STUDY.md)
(the parent study's pre-registration).

---

## 1. The question this answers

The pipeline's template has always come from **8F7Q**, the mOR–β-endorphin–Gi complex. The
choice was never argued; it was inherited. Table S2 of the source paper, now extracted to
[`../data/reference/wang2023_tableS2.csv`](../data/reference/wang2023_tableS2.csv), says it
was the wrong default.

### 1.1 How many endogenous ligands does mOR have?

More than one, which settles the question of whether to run multiple screens.

| Endogenous μOR ligand | Precursor | μOR pKi | mOR-bound structure |
|---|---|---|---|
| Endomorphin-1 | unidentified | **8.73** | **8F7R** |
| Leu-enkephalin | PENK / PDYN | 7.97 | none |
| Met-enkephalin | PENK | 7.40 | none |
| β-Endorphin (1-31) | POMC | 7.11 | **8F7Q** |
| Dynorphin A (1-13) | PDYN | 6.02 | 8F7W (κOR, not μ) |

Five endogenous ligands bind μOR; **exactly two have an experimentally determined
mOR-bound pose.** A structure-based pharmacophore needs a pose, so the number of
templatable endogenous ligands is two, not five and not one. That is the whole argument for
this design: two screens, because there are two poses — not because two is a nice number.

The enkephalins are the notable absence. They are higher-affinity than β-endorphin at μOR
and have no mOR complex in this deposition set, so they can only ever be *candidates* here,
never templates.

### 1.2 Why 8F7Q alone is a bad default

Three reasons, all from Table S2 and all independent of each other:

1. **Affinity.** Endomorphin-1 binds μOR at pKi 8.73 (Ki 1.86 nM) against β-endorphin's
   7.11 (77.6 nM) — **42-fold tighter**. β-endorphin is the 4th-best of the five.
2. **Selectivity.** β-endorphin is effectively non-selective μ vs δ (7.11 vs 6.95, ΔpKi
   0.16 = 1.4-fold). Endomorphin-1 is the only μ-exclusive peptide in the set: δ, κ and
   NOPR are all `N.A.`, i.e. ≥2.4 log units of selectivity. Since
   [`structural_basis.md`](structural_basis.md) §5 established that a YGGF message template
   **cannot be μ-selective by construction** (YGGF alone activates all four receptors),
   the only route to a selective template is a ligand whose selectivity is encoded in the
   message region itself. Endomorphin-1 is that ligand: its W3/F4 occupy an extended
   TM1/2/7 pocket β-endorphin does not use (source paper Fig. 2F–H).
3. **Completeness of the pose.** This is the strongest practical argument and it is new in
   this document. See §3.

---

## 2. The two templates

| | **T_A** | **T_B** |
|---|---|---|
| Source | 8F7Q | **8F7R** |
| Ligand | β-endorphin (1-31) | **endomorphin-1** |
| μOR pKi | 7.11 | **8.73** |
| Selectivity | μ≈δ | **μ-exclusive** |
| Residues resolved in deposition | 21 of 31 | **4 of 4** |
| Residues in template | 8 (a truncation) | **5 — the entire ligand incl. amide cap** |
| Binding mode | YGGF message | **YPWF, extended TM1/2/7 pocket** |
| Complex file | `data/8f7q_complex_p8.pdb` | `data/8f7r_complex_p4.pdb` |
| Exclusion volumes | yes (α=1.0) | yes (α=1.0), to be built |

Both templates get exclusion volumes. This is not optional: `structural_basis.md` §4
establishes that the recognition chamber's *steric bound* is the part the score was
missing, and [`dietary_screen.md`](dietary_screen.md) measured the consequence — AUROC
0.917 with exclusion volumes against 0.444 without, on the same template and same set.
A template without them is known-broken and is not run here.

---

## 3. T_B needs no truncation decision, and that matters

8F7Q resolves 21 of β-endorphin's 31 residues, and the pipeline uses 8 of those 21. Two
truncations sit between the molecule and the template, and
[`../DECISIONS.md`](../DECISIONS.md) D2.6 records that getting the first one wrong (5
instead of 8) silently removed every ECL2 contact.

8F7R resolves **all four** residues of endomorphin-1 plus its C-terminal amide cap
(`YPWF` + `NH2`, chain P residues 1–5, verified in
`data/8f7r_complex_p4.provenance.json`). The template is the complete ligand. There is no
extent parameter to choose, no address region left unmodelled, and no possibility of the
D2.6 class of error.

**Coordinate precision.** 8F7R contains two independent copies of the complex (chains R/P
and M/Q). After superposing the receptors (287 Cα, RMSD 0.001 Å — the copies are
coordinate-identical), the two endomorphin-1 copies differ by **0.635 Å** heavy-atom RMSD
(per-residue 0.43/0.50/0.72/0.74 Å). That is an empirical floor on how precisely this pose
is known, and it is comfortably inside Align-it's feature-matching tolerance, so the
template does not depend on which copy is used. Copy **R/P** is primary;
`data/8f7r_complex_p4_copyMQ.pdb` is retained so the control can be re-run.

---

## 4. Re-pointing the held-out test (E4 → E4′)

Promoting 8F7R to a template spends the study's only held-out structure: C2/E4 was
"endomorphin-1 retrieved at the 98.8th percentile from held-out 8F7R", the strongest
independent evidence in [`screen_validation.md`](screen_validation.md). It cannot remain a
held-out test once it is a template.

**Replacement: 8EFQ, the mOR–DAMGO–Gi complex.** It qualifies on every count that matters:
never used anywhere in this project, deposited by a *different* paper (so independent of
both templates), and DAMGO is the highest-affinity ligand in Table S2 (pKi 8.88, Ki 1.32
nM). `data/8efq_complex_damgo.pdb` is built with provenance.

**One obstacle, and how it is handled.** DAMGO is Tyr-D-Ala-Gly-N-MePhe-Gly-ol. The
candidate builder takes one-letter sequences; it already handles D-residues via the
lowercase convention (`conformers.parse_sequence`) and C-terminal amides via `c_term`, but
it cannot express **N-methyl-Phe** or the **C-terminal glycinol**. So the target table
gains an optional `smiles` column, and a target carrying one is built from SMILES instead
of from sequence. The conformer protocol is otherwise identical — same ensemble size, same
seed, same embedding parameters — so the resulting percentile is comparable with the
library distribution rather than being a special case.

DAMGO is a synthetic analogue, not an endogenous ligand. It is used here *only* as a
held-out retrieval target, never as a template and never in the endogenous-ligand analysis.

---

## 5. Run matrix, and why the window is 8 for both

### 5.1 The runs

The dietary set (20 canonical sequences, 14 agonists / 6 antagonists, 19 with published
IC50) is the evaluation set, for the reasons given in `dietary_screen.md` §1: externally
measured potencies, the actual domain of the library, and both agonists and antagonists.

| Run | Template | Window | Conformers | Seeds | Role |
|---|---|---|---|---|---|
| `A8` | T_A + excl | 8 | 800 | 20261020, 61453 | **already run** — `results/dietary/p8excl_s*` |
| `B8` | T_B + excl | 8 | 800 | 20261020, 61453 | **primary for T_B** |
| `B4` | T_B + excl | 4 | 800 | 20261020, 61453 | sensitivity only — see §5.2 |

`A8` and `B8` share window, conformer count, seeds, metric and exclusion-volume treatment.
**The only thing that differs is which template**, which is what makes the comparison
interpretable.

### 5.2 Why not window 4, when T_B's ligand is 4 residues?

Matching the candidate window to the template extent is the project's own stated principle
(`22_prepare_complex.py` docstring: "Matching the reference extent to the candidate window
makes the comparison like-for-like"), and it is why T_A is screened at 8. Applied to T_B it
would give window 4.

These figures are computed on the **actual E8 pools**, which is the correction that
matters here. E8 does not score all 20 dietary peptides: `E8_all_published` keeps the 17
with `published=yes` and `terminal_supported=yes` (14 agonists / **3** antagonists), and
`E8_mu_assays_only` further restricts to μ-dominant assays (12 / **3**). The headline
0.917 is the μ pool. An earlier draft of this section computed the aromatic composition
over all 20 peptides and 6 antagonists, concluded that the confound *reverses* at window
4, and was wrong; the pool-correct numbers are below.

| | window 4 | window 8 |
|---|---|---|
| Aromatics, agonists — `clean` pool (14) | 1.79 | 2.00 |
| Aromatics, antagonists — `clean` pool (3) | 2.00 | 2.67 |
| Antagonist − agonist, `clean` | **+0.21** | **+0.67** |
| Antagonist − agonist, `mu` (12 / 3) | **+0.08** | **+0.67** |
| Confound direction | opposes, but barely | **opposes** |
| Distinct scored windows, `mu` pool | **9 of 15** | 15 of 15 |

Two reasons for window 8, one weaker than the first draft claimed and one decisive:

1. **The confound protection thins to nothing.** The reason the 0.917 result was worth
   reporting is that antagonists in this set are *more* aromatic than agonists, so a score
   that merely counted aromatic groups would rank them **higher** — the confound opposes
   the hypothesis, and a win cannot be explained by composition. That protection does not
   reverse at window 4, but on the μ pool it shrinks from **+0.67 to +0.08 aromatic
   residues**, which is indistinguishable from no protection at all. Window 4 would
   therefore not *invert* the safeguard, it would *remove* it.
2. **Six of fifteen μ-pool peptides become degenerate.** Window 4 collapses
   β-casomorphin-4/5/7 onto `YPFP`, three soymorphins and two human β-casomorphins onto
   `YPFV`, the two αs1-casein exorphins onto `RYLG`, and two lactoferrin peptides onto
   `YLGS`. Identical window means identical score, which destroys the within-series
   potency ladder the ranking test depends on. This is the decisive objection.

So window 8 is primary for both templates. `B4` is still run, because the extent-matching
principle deserves a measurement rather than an argument, but it is **pre-registered as
weak evidence at best**: with the confound margin at +0.08 a positive B4 result is close to
what a composition artefact alone would produce, and 6 of 15 peptides carry no independent
information. It can be read as *agreement* with B8 (reassuring) or *disagreement* (a flag
that window sensitivity exceeds template choice), not as support for T_B on its own.

**Power, stated plainly.** Every E8 number in this project rests on **3 antagonists**. The
leave-one-antagonist-out range is reported for exactly this reason, and no AUROC here —
existing or new — should be read as more than suggestive.

This leaves T_B with a 5-feature-source template scored against an 8-residue window — an
extent mismatch, which is the mechanism behind the ρ = +0.909 aromatic correlation. The
exclusion volumes are what is expected to contain it, and §6 T3 is the test of whether
they do.

---

## 6. Pre-registered predictions

Primary metric throughout: `tversky_ref_mean` over the 800-conformer ensemble, the
statistic selected in D4.2/D6.8.

**Pool and seed handling, fixed here to match the existing E8 exactly**, so that T1 is
comparable with the published 0.917 rather than merely resembling it: the primary pool is
`E8_mu_assays_only` (12 agonists / 3 antagonists), the secondary is `E8_all_published`
(14 / 3), and scores are **averaged across the two seeds** before the AUROC is taken —
which is what `23_dietary_screen.py` does. Per-seed AUROCs are additionally reported, and
a test counts as holding only if it holds in the seed-averaged statistic **and** does not
reverse in either individual seed.

| # | Test | Prediction, stated in advance | Falsified if |
|---|---|---|---|
| **T1** | `B8` agonist-vs-antagonist AUROC, μ pool, seed-averaged | ≥ 0.75, and point estimate ≥ `A8`'s **0.917** | < 0.75, or clearly below `A8` |
| **T2** | Consensus (§7) AUROC, μ pool | ≥ max(`A8`, `B8`) | below both single templates |
| **T3** | Paired `B8` vs `A8`, same window, same seeds | `B8` ≥ `A8`; Spearman ρ between the two score vectors < 0.9 (the templates see *different* things) | ρ ≥ 0.9, which would mean T_B adds nothing and the second screen is redundant |
| **T4** | `B8` Spearman ρ vs pIC50, μ-restricted subset | weak positive, 0 < ρ < 0.6. A confirmation, not a discovery — E2 gave 0.496 and E9 is harder | ρ < 0, which would contradict the agonist result |
| **T5** | `E4′`: DAMGO percentile vs the library distribution, both templates | top decile (≥ 90th) for at least one template | below the 90th for both |
| **T6** | Aromatic-count correlation, `B8` | ρ(score, aromatic count) < 0.909, i.e. the complete-ligand template with exclusion volumes is *less* composition-driven than the truncated one | ρ ≥ 0.909 |

Because every AUROC rests on 3 antagonists, each is reported with its exact one-sided
permutation p and its full leave-one-antagonist-out range. **A result whose LOO range
straddles 0.5 is not claimed**, whatever the point estimate.

**What the whole design cannot show.** T_B is still a message-region template with no ECL2
contact (endomorphin-1 is four residues; it does not reach ECL2). The μ-selectivity
argument of §1.2 is a *mechanistic hypothesis* drawn from the source paper's Fig. 2F–H,
not something this screen measures — nothing here tests δ/κ/NOPR. The ceiling in
`structural_basis.md` §5 stands, and no result from this design licenses a selectivity
claim.

---

## 7. The consensus statistic

Two templates give two score vectors over the same peptides. Pre-registered combination,
fixed before scoring:

1. **Rank fusion (primary).** Per template, rank peptides by `tversky_ref_mean`; the
   consensus score is the mean of the two ranks. Scale-free, so it needs no calibration
   between templates whose score distributions differ.
2. **Max-normalised score (secondary).** Per template, min-max scale scores to [0,1]
   across the scored set; take the per-peptide maximum. This encodes "a peptide is
   interesting if it matches *either* endogenous binding mode", which is the biologically
   motivated reading — a dietary peptide need only mimic one natural ligand.

Both are reported. Rank fusion is primary because max-of-two is upward-biased by
construction (two chances to score well), and that bias is exactly what the antagonists
also get, so it is not self-cancelling.

---

## 8. Execution

Scoring requires Align-it, which builds only against the conda toolchain on ORCA
(`env/build_alignit.sbatch`); it is not available on the workstation. Split accordingly:

**Done locally** (no Align-it needed): structures downloaded; all three complexes built
with provenance; two-copy precision control (§3); aromatic/degeneracy analysis (§5.2);
SMILES target support; this document.

**On ORCA**: `sbatch jobs/dual_template.sbatch` — builds T_B and its exclusion volumes,
then runs `B8` and `B4` at both seeds. Then `python scripts/24_dual_template_consensus.py`
for §6 and §7.

---

## 9. What would make this design wrong

Stated plainly, so the result is readable either way:

- If **T3** shows ρ ≥ 0.9 between the two templates' score vectors, the second screen is
  measuring the same thing as the first and the dual-template design is unjustified
  regardless of how T1 comes out.
- If **T6** fails and T_B is *more* aromatic-correlated than T_A, then the complete-ligand
  advantage of §3 is outweighed by the extent mismatch of §5.2, and the right response is
  to revisit the window rather than the template.
- If **T1** passes but **T5** fails for both templates, the templates work on the dietary
  set but generalise to no independent structure, which would point at the dietary set
  rather than at the templates.
