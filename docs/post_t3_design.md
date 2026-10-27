# After T3: what replaces the dual-template design

**Written 2026-10-27, before the runs it specifies.** Pre-registration in the sense this
project uses it: every threshold below is fixed here, and the gate scripts read these
numbers rather than a judgement made after seeing them.

[`dual_template_result.md`](dual_template_result.md) reports ρ(A8, B8) = +0.917 and
[`dual_template_design.md`](dual_template_design.md) §9 says plainly what that means —
the dual-template design is unjustified. This document decides what follows, and the
central claim is that **the template was never the binding constraint.**

---

## 1. Two independent reasons the design failed, neither of them "endomorphin-1"

**Statistical.** [`power_analysis.md`](power_analysis.md): on 12 agonists × 3 antagonists,
the power to resolve a true 0.15 AUROC gap between two templates correlated at ρ = 0.9 is
**0.14**, and 0.01–0.06 for the smaller gaps actually at issue. T1 and T2 were decided by
a class of three. The observed A-vs-B gap of +0.083 has a paired bootstrap CI of
[0.000, +0.333]. No template comparison is adjudicable on this set — so the dual-template
experiment could not have returned a usable answer whichever structure supplied T_B.

**Mechanistic.** Truncation to the message epitope is expected to compress score
separation, because the shared anchor is what the window keeps and the discriminating
address is what it discards. Anchor-residue truncation is established practice
(Rajamani et al. 2004, *PNAS* 101:11287, [10.1073/pnas.0401942101](https://doi.org/10.1073/pnas.0401942101)
— anchors bury ≥100 Å² and already adopt near-bound conformations unbound), and the
opioid message–address model is its oldest instance (Chavkin & Goldstein 1981, *PNAS*
78:6543, [10.1073/pnas.78.10.6543](https://doi.org/10.1073/pnas.78.10.6543): the
N-terminal four residues carry agonist activity, the C-terminal domain is the address).
The cost of keeping only the message is that actives and decoys come to look alike. ρ =
+0.917 between two **message-region** templates is what that looks like in a number.

Both reasons predict the same thing: **another short template will be redundant again.**

---

## 2. The structural ceiling in `structural_basis.md` §5 is wrong, and the correction is one residue

§5 states that reaching the selectivity-determining region requires spanning the "address
region proper" at β-endorphin residues ~13–21, that a 21-residue reference cannot be
conformationally enumerated, and that this is "a structural limit on the method, not a
tuning problem."

**Measured from the 8F7Q coordinates** (chain P vs receptor chain R, 4.5 Å heavy-atom
cutoff; chain Q agrees). Full per-residue table in the commit that adds this document:

| β-endorphin residue | Regions contacted |
|---|---|
| 1–5 | pocket core, TM7, TM5/TM6 tips, ECL1, ECL2 floor (Phe4↔Cys219 4.1, Thr220 3.7) |
| **6** | **ECL2** — Thr6↔Arg213 3.7 |
| **7** | **TM6 tip** — Ser7↔Lys305 (6.58) 4.2 |
| **8** | **ECL2** — Glu8↔Gln214 3.1, Arg213 3.3 |
| **9** | **ECL2** — Lys9↔Asp218 3.1; ECL1 Gly133 3.5 |
| 10–21 | **ECL1 and TM1 tip only.** Residues 10, 16, 17, 18, 20, 21 contact nothing at 4.5 Å |

Every ECL2 and TM6-tip contact β-endorphin makes lies at residues 4–9. Residues 13–21 —
the "address proper" §5 declares unreachable — barely touch the receptor at all, and half
of them do not touch it. **Not modelling them costs almost nothing.** Lys305(6.58) is the
residue whose mutation changes dynorphin potency 100-fold, and Ser7 reaches it at extent
7.

So the selectivity-relevant geometry is **not** out of reach. Extent 8, already built,
spans the Glu8↔Arg213/Gln214 salt bridge and the TM6 tip. **Extent 9 adds the last one,
Lys9↔Asp218, and then there is nothing further to add.** §5's pessimism should be read as
superseded by this table; its underlying observation (a 21-mer is not enumerable) is true
and now irrelevant.

---

## 3. Options considered and rejected

| Option | Verdict | Basis |
|---|---|---|
| Second short-peptide template — 9WSW/9WSX (endomorphin-1, 2.8 Å), 9ODN (DAMGO, 2.4 Å), 8Y72 (2.65 Å) | **reject** | All 18 non-8F7Q MOR–peptide structures are 4–5-mers stopping at the ECL2 floor; same orthosteric volume, so §1's two reasons apply unchanged. Better resolution does not buy reach. 8F7Q is the only MOR structure resolving past ligand residue 5 |
| Window 11–13, following the dynorphin message–address literature | **reject** | Already measured: [`../DECISIONS.md`](../DECISIONS.md) D3.2, window 12 gives E1 0.348 and ρ = 0.051 against pKi — "reproducible noise". Not an embedding failure; conformer yield is 100/100 at windows 8–13 (measured 2026-10-27). It is a **sampling** failure at 44–48 rotatable bonds |
| Cross-receptor selectivity axis from Wang Table S2 | **reject** | `data/reference/wang2023_tableS2.csv` holds 8 peptides, only 6 with ≥2 receptors measured — a weaker axis than the set it would replace |
| Sequence/descriptor ML for agonist vs antagonist | **reject** | The 2025 GPCR-peptide benchmark's best model used ~124 actives at 1:10; descriptor-based agonist/antagonist classifiers used ~4,600 ligands for AUC 0.795 (Oh et al. 2022). GPCRVS (Latek et al. 2025, [10.3390/ijms26052160](https://doi.org/10.3390/ijms26052160)) refuses peptides ≥7 residues outright. No sequence model was found that classifies 8–30-residue GPCR peptide ligands by functional class |
| Whole-peptide scoring | **reject, unchanged** | D3.3 holds: a 17-mer returns zero conformers under default ETKDGv3 |

**Independent support for the window, worth recording:** GPCRVS's authors truncate ChEMBL
peptide entries to **6 residues** to merge them with small molecules, and state that
peptide agonists of opioid, chemokine and class B receptors tend to bind through their
N-terminal short fragments. A separate group reached the same windowing decision and
published it.

---

## 4. The design

**Extent 9 / window 9, single template, from 8F7Q.** One template, because there is only
one MOR structure that reaches past residue 5 and §1 says a second short one adds nothing.
Extent 9, because §2 shows it captures the last selectivity-relevant contact and nothing
beyond it is worth having. Exclusion volumes at α = 1.0, mandatory — `dietary_screen.md`
measured AUROC 0.917 with against 0.444 without, and a template without them is
known-broken (D2.4).

And **the validation moves off the minority class.** The tests that decided the
dual-template question — T3's ρ between templates, T6's ρ against aromatic count — are
computed over all 20 peptides. The tests that failed to decide anything, T1 and T2, rest
on 3 antagonists. That is not a coincidence, and the new gates are chosen accordingly.

### Pre-registered tests

| # | Test | Threshold, fixed now | Falsified if | Gated on 3 antagonists? |
|---|---|---|---|---|
| **W1** | **Sampling adequacy.** Seed-to-seed Spearman ρ of dietary scores at window 9, 800 conformers, two seeds | **ρ ≥ 0.98**, matching window 8's measured 0.983 | ρ < 0.98 — then 800 conformers do not cover 31–37 rotatable bonds, window 8 is the ceiling **for sampling reasons**, and this design is over | no |
| **W2** | **The extent-9 template contains the new feature.** An ECL2-directed feature corresponding to Lys9↔Asp218 is present in the built `.phar` | present, and absent from the extent-8 template | absent — then extent 9 is extent 8 with extra atoms and buys nothing | no |
| **W3** | **Non-redundancy with extent 8.** Spearman ρ between extent-9 and extent-8 dietary score vectors | **ρ < 0.9**, the same threshold T3 used | ρ ≥ 0.9 — extent 9 re-measures extent 8, and the honest conclusion is that window 8 is the whole method | no |
| **W4** | **Composition confound does not return.** ρ(score, aromatic count) at window 9 | < 0.909, as T6 | ≥ 0.909 — the exclusion volumes stop working at this extent | no |
| **W5** | **Agonist/antagonist discrimination**, μ pool, seed-averaged | **reported with a logit CI and an exact permutation p, and claimed only if the leave-one-antagonist-out range stays above 0.5** | — | **yes — and therefore not a gate** |

W5 is deliberately not a gate. At n_antagonist = 3 its power against a true AUROC of 0.85
is 0.62 and its minimum attainable p is 0.0022, so it can neither confirm nor refute this
design. It is reported because withholding it would be worse, and it is explicitly not
load-bearing. **W1 and W3 are the decision**: W1 says whether window 9 is computable at
all, W3 says whether it is a different measurement from window 8.

W3 is the test this design is most likely to fail, and failing it is a publishable
result: it would say the method saturates at the message epitope and that the ECL2
contacts, though geometrically present in the template, do not change what the score
sees.

---

## 5. What is required before any AUROC claim, and it is not more compute

Hanley–McNeil variance goes as 1/n_pos + 1/n_neg, so the class of three dominates
regardless of everything else. From [`power_analysis.md`](power_analysis.md):

| n antagonists | power, true AUROC 0.85 | power to resolve a 0.15 template gap |
|---|---|---|
| 3 (current) | 0.62 | 0.14 |
| 8 | 0.90 | 0.56 |
| 12 | 0.95 | 0.74 |
| 20 | — | 0.90 |

**Eight antagonists** is the point where a single-template claim becomes properly powered;
**twelve to twenty** is what adjudicating two templates needs. More agonists, more seeds,
more conformers and better templates leave 1/n_neg untouched. Whether that many MOR
antagonist peptides with quantitative data exist in the public record is a separate
question and is being surveyed; **until it is answered, no AUROC-based claim in this
project should be strengthened**, including the exclusion-volume result.

## 6. What this design cannot show

Unchanged from `structural_basis.md` and restated so it is not quietly dropped. Nothing
here tests δ/κ/NOPR, so no result licenses a μ-selectivity claim — extent 9 makes the
ECL2 contacts *representable*, which is not the same as demonstrating selectivity.
Separately, the literature offers no precedent for what this project attempts: no
published study demonstrates agonist/antagonist discrimination by pharmacophore or shape
for peptide ligands at any opioid receptor, and the GPCR work that succeeds at functional
discrimination pairs shape with interaction-pattern scoring rather than using shape alone
(de Graaf & Rognan 2008, [10.1021/jm800710x](https://doi.org/10.1021/jm800710x); Kooistra
et al. 2016, [10.1038/srep28288](https://doi.org/10.1038/srep28288)). The closest
published analogue of this whole workflow is Poli et al. 2019, *Molecules* 24:3872
([10.3390/molecules24213872](https://doi.org/10.3390/molecules24213872)) — a
receptor-based pharmacophore screen over a tetrapeptide library at MOR that yielded one
confirmed inverse agonist from three synthesised. That is the realistic yardstick, and
`src/mor_pharmacophore/plip_step.py` already profiles receptor–peptide interactions, so a
shape-plus-interaction arm is reachable without new infrastructure. It is deliberately out
of scope here: W1 and W3 decide whether this geometry is worth scoring at all first.
