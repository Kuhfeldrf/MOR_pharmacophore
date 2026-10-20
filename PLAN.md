# Pharmacophore motif screen against the μ-opioid receptor — execution plan

**Written 2026-10-20, before any screening run.** Acceptance criteria in §3 are
pre-registered: they are written here so that a result can fail them. This follows the
`oprm1-opes-activation` §10 pattern, where criteria written before the data rejected the
first landscape.

Scope: implement and run, end to end on ORCA, the pharmacophore motif-matching screen
described in the Huang Complex seed-fund narrative (§ *Technical approach*, Aim 1,
second screen). One receptor — OPRM1 / P35372 — against the digestion peptide library.
MOR is the calibration anchor, not the target: the deliverable is a method that can be
pointed at the 83 receptors with a characterised endogenous ligand.

CPU only. This stage places no demand on the GPU partition, which is what the narrative
claims for it (under 800 CPU-hours for motif matching, fingerprinting and model fitting).

---

## 1. What exists and what is missing

`fdde167` is a tested skeleton, not a screen. Per the README's own status table:

| Component | State | This plan |
|---|---|---|
| PLIP interaction reading (8F7Q chain P) | working, tested | §4.1 — extend to emit a template file |
| RDKit feature encoding | working, tested | §4.2 — add conformer ensembles |
| Free N-terminus rule | working, tested | §4.5 — wire into tiering |
| Align-it scoring | **wrapper stub**; `alignit_path()` returns None | §4.0, §4.4 — build it, then implement |
| Tiered scoring | not started | §4.5 |
| Library ingestion | not started | §4.3 |
| Validation against controls | not started | §5 |

The binding gap is Align-it. It is **not on conda-forge**, is GPL-3.0, and its upstream
(`OliverBScott/align-it`) was last pushed 2022-01-26. It must be built from source
against OpenBabel 3. Everything downstream of §4.4 is blocked on that build, so it is
stage zero and it is the single largest schedule risk (§7).

---

## 2. The scientific question, stated so it can fail

The sibling repositories have established, by measurement, that **no structure-based
score tested so far beats the free sequence rule "does the peptide start with Tyr?"**
(AUROC 0.875; `MOR_peptide_binding`, branch `pose-generation`,
`docs/pose_generation_comparison.md`). Co-folding confidence, interface rescoring,
MM/GBSA and free docking have each failed that comparison.

The pharmacophore screen is vulnerable to the same critique in a sharper form, because
**its own positional rule is nearly that baseline.** The free-N-terminal-Tyr requirement
and "starts with Tyr" select almost the same peptides. If the screen is scored as a
whole, the sequence rule will do the work and the 3D matching will look successful
without contributing anything.

> **The question this screen must answer is therefore not "does it find opioid
> peptides?" It is: does 3D pharmacophore matching add discrimination *beyond* the free
> sequence rule?**

That question is only answerable **within** the Tyr-starting subset, where the baseline
is constant by construction. Everything in §3 and §5 is built to ask it that way.

---

## 3. Pre-registered acceptance criteria

Evaluated after the run described in §6. A screen that fails C1–C3 is reported as a
negative technical-readiness result with a quantified baseline, exactly as the
pose-generation experiment was.

| # | Criterion | Threshold | Why this number |
|---|---|---|---|
| **C1** | **Template recovery.** The pharmacophore read from 8F7Q must score the ligand it was read from above the library median. | positive control in top decile | A template that cannot retrieve its own ligand is not a template. Weakest possible check; failing it voids the run. |
| **C2** | **Independent-ligand recovery.** Endomorphin (from 8F7R, a structure not used to build the template) must score in the top decile. | top decile | Tests that the template generalises past the one pose it was derived from. 8F7R is an independent experimental determination. |
| **C3** | **Added value over the baseline.** Within the 640 Tyr-starting library peptides, pharmacophore score must separate known opioid-active peptides from the rest. | AUROC ≥ 0.70, bootstrap 95% CI excluding 0.5 | This is §2's question. The sequence baseline is constant across this subset, so any separation here is attributable to the 3D matching. |
| **C4** | **Positional rule is load-bearing.** The 74 peptides that *span* YPFP with a blocked N-terminus must score below the 4 that *start* with it. | median difference, reported with CI | A principled control (§5.2) that isolates the N-terminus rule from motif presence. |
| **C5** | **Conformer stability.** Re-running conformer generation with a different random seed must not change the tier assignment of more than 10% of peptides. | ≤10% tier churn | Cross-seed agreement was the only signal that separated a good pose from a bad one in the Chai-1 work. Same lesson applied here. |

**C3 is the criterion that matters.** C1 and C2 are plumbing checks. A run that passes
C1, C2, C4 and C5 but fails C3 means the screen reduces to the sequence rule — a real
and reportable finding, and one the narrative would need to reflect.

---

## 4. Pipeline design

Seven stages, one script each, Snakemake-driven, matching the `MOR_peptide_binding`
layout so the two are operable by the same person. Every stage writes a provenance
record; stubs exit non-zero with a `STUBBED` marker rather than writing empty output.

### 4.0 Stage 0 — Align-it build *(critical path)*

Build from source on ORCA against OpenBabel 3.2.1 (conda-forge, has `linux-aarch64`,
which matters for the later Vera/ARM port). Deliverables:

- `env/build_alignit.sh` — pinned commit, CMake flags, recorded in `docs/provenance.md`.
- A smoke test: align a molecule to itself and assert the score is maximal, plus align
  two known-dissimilar molecules and assert it is not. `tests/test_alignit.py` currently
  skips when the binary is absent; it stops skipping once this lands.
- **Fallback, decided now rather than under pressure:** if the build cannot be made to
  work, substitute RDKit's own `rdMolAlign` / feature-map scoring
  (`rdkit.Chem.Features.FeatMaps`). It is weaker — no pharmacophore-typed Gaussian
  overlap — but it is pure Python, already installed, and keeps every other stage
  honest. The substitution would be recorded as a deviation, not hidden.

### 4.1 Stage 1 — Template extraction

Input `data/8f7q_complex.pdb` (already in repo: receptor chain R + chain P residues 1–5,
the YGGFM message segment of β-endorphin). Existing `plip_step.py` reads the
interactions. New work: convert that interaction set into a positioned pharmacophore
template — each feature typed (hydrogen-bond donor/acceptor, hydrophobic, aromatic,
positive ionisable) and placed at the coordinate it occupies in the experimental pose.

Output: `results/01_template/mor_template.phar` (Align-it format) plus a human-readable
TSV of features with the receptor residue each one engages.

**Decision D1. The template is read from the experimental pose only.** Not from a
predicted or minimised one. The pose-generation work found that all three predicted
generators rotate the Tyr phenol into an H-bond the experimental structure does not
make; building the template from a predicted pose would bake that error into the screen.

**Decision D2. Template structure is 8F7Q; 8F7R is held out.** 8F7Q is already in the
repo and is the structure the met-enkephalin pose overlay was validated against. 8F7R
(endomorphin) is reserved as the independent recovery check in C2 and is not used to
build or tune the template.

### 4.2 Stage 2 — Candidate conformers

A pharmacophore alignment compares 3D point sets, so each candidate needs conformers.
`features.py` currently embeds a single conformer per sequence (`EmbedMolecule`,
`randomSeed=7`). That is adequate for a 5-mer and inadequate for the library.

**Decision D3 — the N-terminal window.** Score the **N-terminal 5 residues**, not the
whole peptide.

The library's median length is 15 and its maximum 46. The template is 5 residues. A
15-mer has far more conformational freedom than any ensemble can cover, and the extra
residues contribute pharmacophore points that have nothing to match in the template,
which dilutes the score. This is the same failure mode already measured in this project
as ipTM dilution: *"ipTM is known to be diluted when one chain is large and mostly
off-interface, and 281 residues against 5 is close to the worst case."*

The window is justified by the binding model, not only by cost: the opioid message is
an N-terminal motif anchored by the free α-amine at the Tyr, and residues beyond the
message are not what the template describes. Carrying the whole peptide is recorded as
the alternative to test once the windowed screen has a result — if the windowed screen
works, whole-peptide scoring is the obvious next experiment, and if it does not, the
window was not the reason.

Conformers: ETKDGv3, N=50 per window, MMFF-minimised, RMSD-pruned at 0.5 Å. Two seeds
for the C5 stability check.

### 4.3 Stage 3 — Library ingestion

Source: `combine data sets simplifed.xlsx` (Dallas lab, 5 populated sheets). Profiled
2026-10-20:

| Sheet | Peptides (≥2 aa, standard residues) | Tyr-starting |
|---|---:|---:|
| `abbott` | 2,803 | 172 |
| `Preterm_gastric` | 10,758 | 509 |
| `pre-term_intestinal` | 5,804 | 180 |
| `term_intestinial` | 3,349 | 103 |
| **union (deduplicated)** | **14,571** | **640** |

This matches the narrative's "14,500 peptides" and is the provenance for that figure.
One non-standard entry (`B`) and one empty sheet are dropped, with the drop recorded.

The sheet of origin is retained per peptide: `abbott` carries intensity columns
distinguishing human/bovine and gastric/intestinal, and the narrative prioritises
peptides present in the same digestive phase as the receptor. Phase is carried as
metadata from the start rather than rejoined later.

Output: `data/library/peptides.tsv` — sequence, length, source sheets, digestive phase,
`starts_with_Y`, and a stable `peptide_id`.

### 4.4 Stage 4 — Align-it scoring

Each candidate conformer ensemble against the template. Align-it reports Tanimoto-like
pharmacophore overlap; the per-peptide score is the best over its ensemble, with the
ensemble spread retained as the convergence readout (C5).

Array job over peptide blocks. Embarrassingly parallel: no block depends on another.

### 4.5 Stage 5 — Tiering and the positional rule

The narrative specifies graded tiers rather than one cutoff, *"because partial matches
recover known ligands that strict matching misses"*, and a positional requirement.

- Tiers from the score distribution of the **negative** set (§5), not from an absolute
  cutoff, so tiers mean the same thing when the method is pointed at another receptor.
- The free N-terminus rule (`nterm_rule.py`) is applied as a **separate reported flag,
  not folded into the score.** This is what keeps C3 answerable: if the rule is mixed
  into the score, its contribution can no longer be separated from the 3D matching's.

### 4.6 Stage 6 — Validation

Computes C1–C5 and writes `docs/screen_validation.md`, including the bootstrap CIs and
the baseline comparison. Written to be readable as a negative result.

### 4.7 Stage 7 — Candidate output

Ranked candidate list per tier, with digestive phase and the interactions each candidate
reproduces. This is the artefact that feeds co-folding and, later, synthesis and the
PDSP assays.

---

## 5. Control sets

Three sets, all fixed before the run. Two of them fall out of the library itself and are
the strongest part of this design.

### 5.1 Positive controls

The IUPHAR/BPS-curated endogenous ligands of OPRM1 (6 peptides, pKi 8.1–8.5) plus
food-derived opioid peptides from `MOR_peptide_binding/data/reference/`.

**β-casomorphin-7 (YPFPGPI) is present in the digestion library itself** — a documented
food-derived MOR ligand appearing in real digestion data, so at least one positive is a
genuine library member rather than a spiked-in reference.

**The enkephalins and dynorphins are absent from the library**, as expected: they come
from proenkephalin and prodynorphin, not milk. This is a useful property — the template
ligand's exact sequence is not in the screening pool, so the screen cannot succeed by
retrieving a copy of its own template.

**pKi spans only 0.4 log units (8.1–8.5), so these are binary positives and cannot
support a potency-ranking claim.** The potency question was already answered negatively
with far more range (4.43 log units, 88% power, seven metrics) in
`docs/validated_screening_framework.md`. This screen does not re-ask it.

### 5.2 The positional control — motif present, N-terminus blocked

Of the library's 78 peptides containing **YPFP** (the β-casomorphin message), only **4**
begin with it. The other **74 span the motif with the Tyr amine locked in a peptide
bond** — e.g. `AQTQSLVYPFPGPI`, which contains β-casomorphin-7 in full but is not it.

These 74 are a near-perfect negative set: identical motif, identical source protein,
wrong position. They isolate the positional rule from motif presence (C4), and no
sequence-similarity method can separate them from the positives.

### 5.3 The baseline-matched negative set

The 640 Tyr-starting library peptides minus the known actives. Holding the sequence
baseline constant, this is the set on which C3 is evaluated. Composition-matched
scrambles of the positives are retained as a second negative set for continuity with the
existing reference data.

---

## 6. Compute plan on ORCA

CPU only. ORCA partitions available: `short` (4 h), `normal` (24 h), `long` (7 d).

| Stage | Shape | Estimate |
|---|---|---|
| 0 Align-it build | 1 job, 1 node | < 1 h, one time |
| 1 Template | 1 job, serial | minutes |
| 2 Conformers | array, 146 blocks × 100 peptides, 8 CPU each | the cost driver; benchmark first |
| 3 Ingestion | login node | < 1 min |
| 4 Align-it scoring | array over the same blocks | minutes per block |
| 5–7 Tiering, validation, output | serial | minutes |

**No stage-2 estimate is given here on purpose.** Conformer generation for a 5-residue
window is fast, but 14,571 of them has not been timed, and this project has already been
burned once by an unmeasured throughput assumption (the 1-thread/20-hour fault that cost
the met-enkephalin campaign). **Stage 2 is benchmarked on 100 peptides and the real
number recorded before the full array is submitted.** The narrative's "under 800
CPU-hours" is the figure to check against.

Scratch is at **94% (3.8 TB free)**. This screen's output is small — scores and
templates, not trajectories — but conformer SDFs for 14,571 peptides are not nothing.
Conformers are written to a scratch working directory, consumed by stage 4, and not
retained; only scores, the template and the ranked output are kept.

Run root: `/scratch/kuhfeldr-Kuhfeld_temp/MOR_pharmacophore`.

---

## 7. Risks

| Risk | Severity | Mitigation |
|---|---|---|
| **Align-it will not build** (unmaintained since 2022, OpenBabel 3 API drift) | blocks §4.4 onward | Stage 0 is first. RDKit FeatMaps fallback is decided in advance (§4.0). |
| **The screen reduces to "starts with Tyr"** | the scientific risk | C3 is designed to detect exactly this, within the Tyr-starting subset. A failure is reportable, not fatal. |
| Conformer coverage inadequate even for a 5-mer window | moderate | C5 cross-seed tier churn detects it. |
| Template over-fit to one experimental pose | moderate | C2 holds out 8F7R. |
| Scratch fills | moderate | Conformers are transient; only scores retained. |
| ARM/Vera port | deferred | Not this plan. OpenBabel and RDKit both have `linux-aarch64` builds, which is the encouraging part; Align-it is the unknown and is another reason the fallback matters. |

---

## 8. Milestones

| # | Deliverable | Gate |
|---|---|---|
| M1 | Align-it built and smoke-tested, or fallback adopted and recorded | `tests/test_alignit.py` no longer skips |
| M2 | Template extracted from 8F7Q, features tabulated | C1 passes on the template's own ligand |
| M3 | Library ingested, control sets frozen and committed | counts reproduce §4.3 exactly |
| M4 | Stage 2 benchmarked on 100 peptides | real CPU-hour estimate recorded |
| M5 | Full screen run | C1–C5 computed |
| M6 | `docs/screen_validation.md` written | states plainly whether C3 passed |

M1–M3 are this week's work and are what "initiate the process" means here. M5 waits on
M4's benchmark.

---

## 9. What this plan deliberately does not claim

- It does not claim the screen will work. §2 names the way it most plausibly fails.
- It does not re-ask the potency-ranking question (§5.1).
- It does not extend to the other 82 receptors. One receptor, one template, measured.
- It does not touch the GPU partition.
