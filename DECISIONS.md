# Computational decision log

Every non-obvious choice in this pipeline, why it was made, what evidence supports it, and
whether it survived. Written so that any number in
[`docs/optimisation_study.md`](docs/optimisation_study.md) can be traced to the decisions
that produced it, including the decisions that were later wrong.

**Status vocabulary**

| | |
|---|---|
| **HELD** | made, tested, still in force |
| **REVISED** | changed after evidence; both versions recorded |
| **WITHDRAWN** | made, found wrong, removed; kept here because the reasoning is instructive |
| **OPEN** | identified, not yet acted on |

Evidence is either **measured** (a number from this project, with the script that produced
it) or **literature** (a citation). Decisions marked *measured* were tested on this
system; decisions marked *literature* were taken on precedent and should be treated as
weaker until tested here.

Design intent and pre-registered criteria live in [`STUDY.md`](STUDY.md); this file is the
engineering and statistics record.

---

## 1. Source data and library construction

### D1.1 Peptide library from the digestion workbook — HELD
**Decision.** `combine data sets simplifed.xlsx` → `data/library/peptides.tsv`, 14,571
unique sequences. `peptide_id = "pep_" + sha256(sequence)[:10]`.
**Rationale.** A content-derived ID is stable across re-ingests and makes a changed
sequence a changed ID rather than a silent mutation of an existing record.
**Evidence.** `scripts/01_library.py` writes `provenance.json` with the SHA-256 of the
source workbook, so the library can be proved to match a specific file.

### D1.2 Deduplicate by sequence, not by source row — HELD
**Decision.** Identical sequences from different phases/matrices collapse to one peptide,
with the phases/matrices retained as metadata.
**Rationale.** The screen scores chemistry, not provenance. Scoring the same sequence
twice would double-count it in every downstream distribution.

### D1.3 CTOP excluded — HELD
**Decision.** Skipped as `NON_CANONICAL`.
**Rationale.** Contains residues outside the 20-letter alphabet, which
`Chem.MolFromSequence` cannot build. Recorded as an explicit skip rather than a silent
drop.

---

## 2. Receptor and template

### D2.1 Template from 8F7Q, held-out structure 8F7R — REVISED by D2.7
**Decision.** Build the pharmacophore from 8F7Q chain P only; never use 8F7R in template
construction, and use it solely for held-out retrieval (E4).
**Rationale.** A template validated on the structure it was built from proves nothing.
**Evidence.** *Measured* — endomorphin-1 retrieved at the ~99th percentile from 8F7R
(`docs/screen_validation.md`).
**Superseded.** The premise — that 8F7Q is the right template and 8F7R is only a test set
— was never argued; 8F7Q was inherited. Table S2 of the source paper says it is the worse
of the two choices (D2.7). The held-out role moves to 8EFQ (D2.9).

### D2.7 A second template from 8F7R, and two screens rather than one — OPEN, pre-registered
**Decision.** Add template **B** from 8F7R (mOR–endomorphin-1–Gi) alongside template A
(8F7Q, β-endorphin, extent 8), screen the dietary set against both, and report each
template separately plus a consensus. Design and predictions:
[`docs/dual_template_design.md`](docs/dual_template_design.md).

**Rationale — why two, and why these two.** Five endogenous peptides bind μOR
(endomorphin-1, Leu- and Met-enkephalin, β-endorphin, dynorphin A) but only **two** have an
experimentally determined mOR-bound pose, and a structure-based pharmacophore needs a pose.
So the number of screens follows from the number of poses, not from preference. The
enkephalins are higher-affinity than β-endorphin at μOR and have no mOR complex in this
deposition set, so they can only ever be candidates here.

**Evidence that A alone was the wrong default.** *Literature*, from
[`data/reference/wang2023_tableS2.csv`](data/reference/wang2023_tableS2.csv):
endomorphin-1 binds μOR at pKi 8.73 against β-endorphin's 7.11 — **42-fold tighter** — and
is the only μ-exclusive peptide in the set (δ/κ/NOPR all `N.A.`), where β-endorphin is
non-selective μ vs δ (ΔpKi 0.16). *Measured*, this repo: 8F7R resolves **all 4** residues
of endomorphin-1 plus its amide cap, so template B involves no truncation decision at all,
where template A is a truncation (8 of 21 resolved, of 31 total) and D2.6 records what
getting that truncation wrong already cost. The two copies of the complex in 8F7R agree to
**0.635 Å** ligand RMSD after receptor superposition (287 Cα, 0.001 Å), so the pose does not
depend on which copy is used.

**Cost, accepted deliberately.** Promoting 8F7R spends the only held-out structure; see
D2.9 for the replacement.

**What this cannot do.** Template B is still a message-region template — endomorphin-1 is
four residues and does not reach ECL2 — so the μ-selectivity argument is a mechanistic
hypothesis from the source paper's Fig. 2F–H, not something the screen measures. The
pan-opioid ceiling of `docs/structural_basis.md` §5 stands.

### D2.8 Template B screened at window 8, not at its own extent of 4 — OPEN, pre-registered
**Decision.** Screen template B at candidate window **8**, the same window as template A,
despite its ligand being 4 residues. Also run window 4 as a sensitivity check, flagged as
weak evidence at best.
**Rationale.** Extent-matching is this project's own principle (D3.1,
`22_prepare_complex.py`), and it points at window 4. Two measured facts override it.
**Evidence.** *Measured*, computed before any score existed, on the actual E8 pools:
(i) the aromatic-confound margin that makes E8 interpretable — antagonists being *more*
aromatic than agonists, so composition alone would rank them higher — collapses from
**+0.67 to +0.08** aromatic residues on the μ pool at window 4, removing the safeguard
rather than reversing it; (ii) window 4 collapses **6 of 15** μ-pool peptides onto shared
windows (β-casomorphin-4/5/7 → `YPFP`; three soymorphins and two human β-casomorphins →
`YPFV`), which destroys the within-series potency ladder. An earlier draft of the design
doc claimed the confound *reverses* at window 4; that was computed over all 20 peptides
instead of the 15-peptide μ pool and was wrong. Corrected in place.
**Consequence.** Template B carries a 4-residue reference against an 8-residue window — an
extent mismatch, the mechanism behind ρ = +0.909 (`docs/optimisation_study.md` §7.2).
Exclusion volumes are expected to contain it, and T6 is the test of whether they do.

### D2.9 Held-out structure re-pointed from 8F7R to 8EFQ — OPEN, pre-registered
**Decision.** E4 becomes **E4′**: held-out retrieval of **DAMGO** from 8EFQ (mOR–DAMGO–Gi),
replacing endomorphin-1 from 8F7R.
**Rationale.** 8F7R cannot be both a template and a held-out test. Dropping the held-out
test instead would discard the study's strongest independent evidence; 8EFQ restores it.
**Evidence.** *Literature* — 8EFQ is used nowhere else in this project, was deposited by a
different paper than 8F7Q/8F7R (so it is independent of both templates), and DAMGO is the
highest-affinity ligand in Table S2 (pKi 8.88, Ki 1.32 nM).
**Implementation cost, and why it was paid rather than worked around.** DAMGO is
Tyr-D-Ala-Gly-N-MePhe-Gly-ol. The candidate builder handles D-residues (lowercase
convention) and amide caps (`c_term`) but has no way to express N-methyl-Phe or the
C-terminal glycinol, so the target table gained an optional `smiles` column
(`conformers.build_from_smiles`). Scoring DAMGO from its deposited coordinates would have
been less code and was rejected: it would have made DAMGO the one target whose conformers
were not generated like everyone else's, and a percentile against the library distribution
cannot tolerate that. Verified: `C26H35N5O6`, MW 513.60, stereocentres S/R/S = L-Tyr /
D-Ala / L-N-MePhe, 50/50 conformers embedded (`tests/test_conformers_smiles.py`).
**Caveat.** DAMGO is synthetic. It is a held-out retrieval target only — never a template,
never part of an endogenous-ligand analysis.

### D2.2 PLIP for interaction detection — HELD
**Decision.** Derive features from PLIP-detected receptor–ligand interactions rather than
from the ligand alone.
**Rationale.** Makes the template a model of *the interaction*, not of the ligand's own
feature set.

### D2.3 Protonate at pH 7.4 before typing — REVISED (originally omitted)
**Decision.** `obabel -p 7.4` on the complex before feature assignment.
**Why it changed.** *Measured* — the deposited PDB has no hydrogens, so the free
N-terminal amine typed as a neutral donor and the template was built with **zero charge
features**. Protonation recovered the `POSC` feature that corresponds to the Tyr1
ammonium–D147 salt bridge, the single most conserved contact in the opioid message.
**Status.** Held since. A warning fires if `POSC` is absent, so this cannot regress
silently.

### D2.4 No exclusion volumes — REVISED; exclusion volumes added, two independent sets agree
**Decision as made.** The template contains 16 purely attractive features (2 AROM, 6
HDON, 5 HACC, 2 LIPO, 1 POSC) and no exclusion volumes.
**Why this is wrong.** Align-it's score table has a dedicated exclusion-volume overlap
column, which in every run of this pipeline is identically zero: **nothing penalises a
candidate for bulk or for features the pocket cannot accommodate.** This is the
mechanistic origin of the aromatic confound of §7 — a His-rich 8-mer accumulates template
coverage with no counter-pressure.
**Literature.** Both relevant precedents include them. The 4DKL MOR tetrapeptide screen
used "LigandScout default volumes derived from the receptor structure"
([PMC6865014](https://pmc.ncbi.nlm.nih.gov/articles/PMC6865014/)); apo2ph4 keeps all
exclusion volumes, removing only duplicates
([PMC9832483](https://pmc.ncbi.nlm.nih.gov/articles/PMC9832483/)).
**Predicted effect, stated before testing.** Receptor-derived exclusion spheres should
lower scores for over-large, feature-dense windows more than for message-like ones, i.e.
reduce E6 *without* the proportional loss of E1 seen when the metric was changed
(§8.1 of the study). If E6 falls and E1 holds, the confound was an artefact of an
incomplete template rather than a property of the method. If both fall together, the
collinearity is intrinsic. **Either outcome is informative, which is why this is the next
thing to run.**

**Measured outcome — the prediction held, on two independent sets.** Implemented in
`scripts/20_template_exclusion.py` (86 spheres over 32 contact residues including ECL2,
placed on receptor heavy atoms within 5.0 Å of ligand, deduped at 1.5 Å, cleared 2.0 Å
from attractive features) and plumbed through `alignit.score(with_exclusion=True)`.

Note first a trap that cost a run: **`align_it` parses EXCL features and silently ignores
them unless `--withExclusion` is passed.** A template carrying exclusion volumes scored
without the flag returns numerically identical results, so the change can appear to do
nothing. Recorded in the `alignit.py` docstring.

| evidence | without exclusion | with exclusion (α=1.0) |
|---|---|---|
| 14-peptide diagnostic (job 241229): message peptides − His-rich windows | **−0.0444** (His-rich win) | **+0.0677** |
| dietary E8 AUROC, agonist vs antagonist, μ-only (job 241252) | 0.444 (chance, exact p=0.63) | **0.917, exact p=0.0154** |
| dietary E6 ρ vs aromatic count | +0.149 | **−0.261** |
| dietary S1 seed reproducibility | 0.942 | **0.994** |

E6 fell and performance *rose* — the outcome the prediction labelled "the confound was an
artefact of an incomplete template rather than a property of the method". The dietary
result is the stronger of the two because its confound points the wrong way (antagonists
are *more* aromatic, 2.67 vs 2.00), so it cannot be produced by composition; see
[`docs/dietary_screen.md`](docs/dietary_screen.md).

α=1.0 over α=0.5: α=0.5 gave a larger gap (+0.0943) but crushed scores to 0.03–0.17, and a
compressed range against an unchanged *absolute* jitter costs resolution.

**Cost.** `--withExclusion` scoring is ~2.5× slower (559 s → 1374–1494 s per block at the
p5 template; scoring is ~20% of a run). At the p8 template it is far worse — 6921 s for a
single 14-peptide block, which overran the 4 h `short` partition on job 241247.

**Still open:** neither set is the pre-registered 235-peptide benchmark. Stage C (job
241259) is what moves this to HELD. The seed jitter was re-measured rather than assumed:
scores compress under exclusion (median 0.228 → 0.154) but the absolute jitter does **not**
(0.00105 either way), so margins in jitter units remain comparable. `23_dietary_screen.py`
now measures jitter per configuration instead of taking it from a flag.

### D2.5 No mandatory features — OPEN
**Decision as made.** Pure similarity score; no feature is required to match.
**Literature.** The 4DKL tetrapeptide screen made the positive-ionizable feature
*mandatory* plus ≥2 others, retaining 28,070 of ~198,000 peptides as a filter before any
ranking ([PMC6865014](https://pmc.ncbi.nlm.nih.gov/articles/PMC6865014/)).
**Why it matters here.** Requiring `POSC` (the Tyr1 ammonium anchor) would enforce the
message motif instead of letting accumulated feature count substitute for it — a
filter-then-rank design rather than rank-only.
**Cost correction.** This is *not* free re-analysis, as first recorded here. Align-it's
score table reports only a *count* of matched features, not their identities, so enforcing
a geometric match on a specific feature needs the aligned-pharmacophore output
(`align_it -o`) and therefore a rescreen. What *is* free is the weaker version — requiring
the candidate to *possess* a positive-charge feature at all, i.e. an unblocked
N-terminus — but that is close to the free-N-terminus flag already carried separately
(D6.11) and is not the same test.

---

## 3. Peptide preparation

### D3.1 N-terminal window truncation — HELD, with a documented cost
**Decision.** Score only the N-terminal *w* residues; `window = 8` in production.
**Rationale.** Opioid activity requires a free N-terminal Tyr and the message motif is
N-terminal; whole peptides are also conformationally intractable (D3.3).
**Cost, measured and material.** Peptides sharing a window are one molecule to the screen.
14,595 targets collapse to 5,183 distinct windows; the 640 Tyr-starting library peptides
collapse to 234 windows. Everything past residue 8 is invisible to the score. This is
recorded in every output file via the `window` column so no ranking can be read without
it.

### D3.2 Window 8 specifically — HELD
**Evidence.** *Measured* — window 5 collapses both E1 families to one point each, so the
central test is undefined; window 12 fails E1 at 0.348 with ρ = 0.051 against pKi, i.e.
reproducible noise at 57 rotatable bonds. The usable window is interior.

### D3.3 Whole-peptide scoring abandoned — HELD
**Evidence.** *Measured* — a 17-mer has 81 rotatable bonds and default ETKDGv3 returns
**zero** conformers; with random coordinates 20 take 73 s, so 500 would cost ~30 min for
one peptide.

### D3.4 Protonation at pH 7.4 — HELD
Applied to candidates as well as the template, for consistency of charge features.

### D3.5 C-terminal amidation of truncated windows — HELD
**Decision.** Cap the artificial C-terminus of a truncated window as an amide; leave
genuinely free acids free; honour known amidation.
**Rationale.** Truncation creates a carboxylate that does not exist in the parent peptide.
Leaving it would add a spurious negative charge and acceptor at a position that is, in the
real molecule, a peptide bond.

---

## 4. Conformer generation

### D4.1 ETKDGv3 with `useRandomCoords=True` — REVISED (originally default coords)
**Why it changed.** *Measured, and this was a silent corruption.* Asking for 50 conformers
on the N-terminal window of dynorphin A returned 48 / 21 / 2 / 1 at windows 5 / 8 / 10 /
12, **with no error raised.** The first screen ran at window 5 and got 48 of 50, which is
why it never surfaced. With random starting coordinates all return 50.
**Guard.** An ensemble below `MIN_CONFORMER_FRACTION` (0.5) of the requested size now
raises rather than being scored.

### D4.2 800 conformers in production — HELD
**Evidence.** *Measured* — window 8 at 400 conformers fails the S1 gate at ρ = 0.945;
800 gives ρ = 0.983. The doubling also lowered E1 from 1.000 to 0.803, which is the
correct direction (§6.2).
**Literature tension — sharper than first recorded, and still unresolved** (McNutt et
al., *J Chem Inf Model* 2023, 63:6598, read in full 2026-10-22). Their recommendation is
not a cost argument and it is specifically about pharmacophore screening:

> "in screening tasks, generating larger ensembles can decrease performance … due to
> increasing the number of false positives … a reduced conformational ensemble of less
> than 25 conformers … is likely sufficient … We emphasize that this recommendation is
> **not primarily motivated by the reduced computational demands** … but by the potential
> **decrease in accuracy** that arises when generating more conformers, which **increases
> the false positive rate faster than the true positive rate**."

Taken at face value this says 800 conformers is actively harmful, and that the pipeline's
single most expensive decision is wrong by a factor of >30.

**Why it is probably not a contradiction — the mechanism depends on D5.5.** Conventional
pharmacophore screening is a *hit-or-miss* match, equivalent to the `max` statistic: one
conformer matching is enough. Under `max`, each additional conformer is another chance at
a spurious match, so the false-positive rate grows with ensemble size — exactly the
failure they describe. This pipeline scores the ensemble **mean**, where additional
conformers improve an average rather than adding draws to a maximum, and the measured
consequence is the opposite: reproducibility *rises* with ensemble size (ρ 0.945 → 0.983
from 400 → 800). Their result and ours are consistent if the mechanism is `max`-specific.

**Tested, and the reconciliation above is wrong.** `scripts/21_ensemble_size_statistic.py`
over 10 / 25 / 50 / 100 (and separately 400 / 800) conformers at window 8, both seeds:

| statistic | 10 → 100 conformers | 400 → 800 |
|---|---|---|
| max | dS1 **+0.095**, dS2 **+0.050**, dE3 **+0.092** | dS1 −0.031, dS2 −0.026, dE3 +0.025 |
| mean | dS1 **+0.302**, dS2 +0.031, dE3 +0.047 | dS1 +0.038, dS2 −0.012, dE3 +0.002 |
| top10 | dS1 +0.301, dS2 +0.021, dE3 **+0.108** | dS1 +0.004, dS2 −0.008, dE3 +0.001 |

In the 10–100 regime **everything improves, `max` included**, so the penalty is not
`max`-specific and the proposed mechanism is not what is going on. Across the whole range
under `mean`: S1 0.566 → 0.868 → 0.945 → 0.983 and E3 0.693 → 0.740 → 0.740 → 0.742, with
S2 peaking near 400 and declining by 0.012 at 800.

**The actual resolution is molecular flexibility, which was the caveat first recorded here
and then over-corrected away from on reading their emphatic phrasing.** Conformer count
needed to cover accessible space scales with rotatable bonds. McNutt et al. benchmark
drug-like ligands; an 8-residue peptide window has ~37 rotatable bonds, for which 25
conformers is severe undersampling — measured S1 = 0.576 at 25, against a 0.95 gate. Their
recommendation is regime-specific and does not transfer to flexible peptides.

**What survives of their warning:** S2 does decline slightly beyond ~400 conformers under
every statistic, so the false-positive effect may exist weakly at the top of our range. It
is an order of magnitude smaller than the reproducibility gain, so 800 stands — but 400 is
defensible and cheaper, and the only metric that prefers 800 is S1.

**Lesson recorded:** the first instinct ("that was for drug-like molecules") was right. The
error was inventing a mechanism to reconcile a conflict that was really a domain mismatch,
and then running a test (400 vs 800) that could not address the question because it did not
span the regime in dispute.

**Two of their findings support choices already made:** selecting the lowest-energy
conformers does *not* improve retrieval of bioactive conformations, while rmsd-based
diversity filtering does (supports D4.3's `pruneRmsThresh`, and they note stricter
thresholds matter more for pharmacophore search); and ETKDGv3's built-in geometric and
knowledge-based constraints "need less refinement" than ML-generated conformers
(supports D4.1).

### D4.3 `pruneRmsThresh = 0.5`, MMFF minimisation — HELD
Standard settings; not independently tested here. Marked as *untested* rather than
validated.

### D4.4 Thread count from the SLURM allocation — HELD
**Why.** *Measured* — RDKit's `numThreads=0` means every core on the *machine*, not every
core in the *allocation*. Eight 8-core tasks on a 64-core node each spawned 64 threads and
drove measured CPULoad to 319. `allocated_threads()` reads `SLURM_CPUS_PER_TASK`, then
`OMP_NUM_THREADS`, then `sched_getaffinity`. This is the mirror image of the
`OMP_NUM_THREADS=1` fault that cost this project an earlier campaign.

### D4.5 Fixed, recorded seeds; two seeds per configuration — HELD
Reproducibility (S1) is a *selection gate*, so it has to be measurable; one seed cannot
measure it. Every run records its seed in the output path.

### D4.6 D-amino-acid support via α-carbon inversion — HELD
**Decision.** Lowercase = D (`YaFGYPS` is [D-Ala²]dermorphin). Implemented by inverting
the chiral tag of the alpha carbon identified through RDKit's PDB residue info, since
`Chem.MolFromSequence` builds all-L.
**Rationale.** Required for E7 (§6.5); no other contrast available here holds composition
*exactly* fixed.
**Evidence.** *Measured* — inverts residue 2 from S to R at identical molecular formula
and identical stereo-free skeleton; D-glycine correctly refused.

---

## 5. Scoring

### D5.1 Align-it v1.0.4 — HELD
Built from source with the conda toolchain (conda-forge OpenBabel needs a newer CXXABI
than RHEL's libstdc++ provides). Version is recorded in every run log.

### D5.2 Score-table column layout verified empirically — HELD
**Why.** The columns were initially *guessed* and the `overlap` index was wrong. Corrected
against actual v1.0.4 output. An 11-column layout is now asserted in code.
**Lesson recorded.** Parsing an external binary's output format from assumption is a
silent-corruption risk of the same class as D4.1.

### D5.3 TVERSKY_REF as the similarity metric — HELD, and literature-confirmed
**Decision.** Reference-weighted Tversky, i.e. what fraction of the *template's* features
the candidate covers.
**Evidence, measured.** Selected on S1 + S2 before the evaluation metrics were computed.
**Evidence, literature.** Independently the recommended choice: Tversky weighted toward
query features was "by far the most successful and the least failure-prone" of ten
dissimilarity scores in a large pharmacophore-descriptor screening benchmark, with
Tanimoto failing to defend its privileged status and performance insensitive to tuning α
([summary](https://www.researchgate.net/publication/237015604_Do_Not_Hesitate_to_Use_Tversky-and_Other_Hints_for_Successful_Active_Analogue_Searches_with_Feature_Count_Descriptors)).
Tanimoto additionally carries a documented bias toward smaller molecules and stops
decreasing as the second molecule grows
([shape-similarity review](https://pmc.ncbi.nlm.nih.gov/articles/PMC6068280/)).
**Known cost, not a defect.** The asymmetry that makes it robust to size mismatch is
exactly what lets a feature-dense candidate accumulate coverage. The right remedy is
exclusion volumes (D2.4), not a different metric — see D5.4.

**Support qualified 2026-10-22, on reading the primary paper** (Horvath, Marcou & Varnek,
*J Chem Inf Model* 2013). Two corrections to how this entry read the evidence:

1. **The descriptors are topological, not 3-D.** ISIDA *colored fragment / pharmacophore
   triplet count* vectors over ~0.75 M similarity-search runs on 165 ChEMBL targets. The
   headline result — query-biased Tversky the best and least failure-prone of ten scores,
   Tanimoto failing to defend its status, performance insensitive to α (0.9 and 0.7 both
   beat Tanimoto) — is a 2-D feature-count result. It is **weaker support for a 3-D
   overlap score than this entry originally implied**, and is a conceptual import rather
   than a transferred benchmark.
2. **Their conditional may point the other way for this target.** They separate
   *growable* queries, which admit active analogues carrying additional features and are
   best served by α > 0.5, from *conservative* queries that allow no growth, whose
   "specific application domain" is **candidate-biased Tversky at α < 0.5**. The opioid
   message pocket is sterically tight — 49 pocket-lining exclusion spheres from 18
   residues (D2.4) — which makes it a conservative query being scored with a strongly
   query-biased metric.

This does not reverse D5.3 on its own, because the measured alternative here was poor:
TVERSKY_DB (the candidate-biased direction) gave E1 = 0.439, i.e. chance (D5.4). But that
measurement was made **without exclusion volumes**, so it tested the asymmetry while the
"no growth" constraint was absent from the template entirely. The two are confounded, and
D2.4 is the experiment that separates them. Retest the metric asymmetry once exclusion
volumes are in (now item 5 in §9).

### D5.4 Metric NOT changed to reduce the confound — HELD (an explicit decision not to act)
**Evidence.** *Measured* — across all nine metric × statistic combinations, the confound
and the evaluation metrics are collinear: tversky_ref/mean E6 +0.913 with E1 0.803;
tanimoto/top10 E6 +0.730 with E1 0.515; tversky_db/mean E6 +0.310 with E1 0.439. Every
metric that weakens the confound weakens the signal proportionally, and the only
combinations clearing a ρ = 0.6 ceiling sit at chance on E1.
**Conclusion.** Composition and shape are confounded in the *benchmark*, not in the
scoring function. `scripts/16_reselect_with_e6.py`.

### D5.5 Mean over the conformer ensemble, not max — HELD, with the interpretive cost stated
**Evidence, measured.** `max` is the noisiest order statistic available: seed-to-seed
ρ = 0.537 at `max` against 0.983 at `mean`. This single change mattered more than any
other in the pipeline — at matched conformer counts it is the difference between failing
and passing the S1 gate.
**Literature position.** `mean` is the *minority* choice; max-over-conformers is commonly
described as standard for shape methods. The one comparison found favoured mean, but
aggregated across pharmacophore *models* rather than conformers
([Molecules 25:385](https://doi.org/10.3390/molecules25020385)), so it does not transfer
cleanly.
**Interpretive cost, explicit.** `mean` asks how well a peptide's accessible ensemble
matches the template *on average*, not whether it can adopt a matching pose. The second
question is closer to what binding requires. It is also the one that is not reproducible
at achievable cost here — and E7 (§6.5) indicates this pipeline cannot resolve pose-level
differences at any setting.

### D5.6 Per-conformer score tables retained — HELD
All ensemble statistics and all three metrics are recomputable post hoc from
`scores.tab` without rescreening. This is why the matched-decoy E3 (§6.6) and the
nine-combination comparison (D5.4) cost no new compute. Storage cost: 1.2 GB for the
library screen, kept on scratch and never committed.

---

## 6. Benchmark design and statistics

### D6.1 Selection metrics separated from evaluation metrics, pre-registered — HELD
**Decision.** S1 (reproducibility, gate ρ ≥ 0.95) and S2 (actives vs scrambles) may
choose the configuration; E1–E5 are reported and never used to choose.
**Enforced in code, not only in prose.** `scripts/11_select.py` is given a selection set
containing *only* actives and scrambles, so the evaluation metrics cannot be computed at
that stage even by accident.
**Why it earned its place.** *Measured* — the gate caught an inflated result twice: the
first study's tier structure (73% decile churn) and E1 = 1.000 at 400 conformers, which
fell to 0.803 once the ensemble converged.

### D6.2 Amendment to the selection rule, made before any E metric was computed — HELD
**Decision.** Choose `n_conformers`/metric/statistic *within* each window; carry all
windows into evaluation.
**Why.** The original rule selected window 5, under which E1 is undefined. Recorded in
`STUDY.md` §4.1 with an explicit timestamp relative to the evaluation metrics, because an
amendment made after seeing evaluation numbers would not be worth anything.

### D6.3 The potency-pair control as the centrepiece — HELD
`YGGFL…` (n = 11, pKi 7.39–8.71) vs `YGPFL…` (n = 6, pKi 4.71–5.51): one Gly3→Pro
substitution, ~3 log units, length- and aromatic-matched, all Tyr-starting. No sequence
rule separates them.

### D6.4 E1 kept as AUROC — HELD (deliberately not changed to BEDROC)
**Rationale.** E1 is a two-sample comparison, where AUROC *is* the Mann-Whitney *U*
statistic and is the correct choice. The early-recognition critique (D6.8) applies to
screening-shaped questions, not to this one. Changing it would have been cargo-culting a
correct criticism into the wrong place.

### D6.5 E7, epimer discrimination, as a sign test — HELD
**Decision.** L/D epimer pairs; score by exact two-sided binomial sign test, not AUROC.
**Rationale.** The literature direction for the L-Ala² substitution is qualitative SAR,
not a matched pKi pair, so there is no measured scale on which to build an interval.
Inventing one would overstate the evidence.
**Scope control.** δ-selective deltorphins carry no directional prediction against a μ
template and are reported but **excluded** from the test, rather than counted in whichever
direction helps.
**Result.** 1 of 3 μ-relevant pairs correct, p = 1.00.

### D6.6 Decoys must be matched on the properties that drive the score — REVISED
**Decision as made.** `tyr_hard` = Tyr-starting library peptides, length-matched.
**Why it was wrong.** Not matched on aromatic count, the property correlating with the
score at ρ +0.909. The benchmark literature names this **artificial enrichment**:
performance attributable to the difference in chemical space between actives and decoys
rather than to binding
([Bioinformatics 37:2134](https://academic.oup.com/bioinformatics/article/37/15/2134/6126797)),
and DUD-E's hidden bias is the canonical case
([PLOS ONE 14:e0220113](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0220113)).
DUD-E's own principle is property similarity with topological dissimilarity; this
benchmark met it for length and the N-terminal residue and violated it for the one
property that mattered.
**Effect, measured.** E3 falls from **0.733 to 0.490** when aromatic count is matched; the
stratum holding 26 of the 36 actives sits at 0.513. `scripts/19_screening_metrics.py`.
**Now required.** A decoy set must be matched on the pharmacophore feature inventory
before anything is run, and E6 is the check that it was.

### D6.7 Percentile stratification of the score — WITHDRAWN
**Decision as made.** Convert each score to its rank within the library windows of the
same aromatic count.
**Why withdrawn.** Failed its own invariance control: E1 should have been unchanged and
moved 0.803 → 0.712, because family members do not share an aromatic count at window 8,
so stratifying split a contrast it was meant to leave alone. It also promoted
zero-aromatic peptides (raw 0.3019) to the 99.96th percentile for topping a stratum in
which nothing matches the template.
**Kept in the repo** (`scripts/14_composition_control.py`) because the reason it failed is
the clearest statement of what the confound is: aromatic count is not a nuisance variable
when the template's anchors *are* aromatics.
**Superseded by** D6.6 — the fix is a matched benchmark, not a score transform.

### D6.8 BEDROC for screening-shaped contrasts — REVISED (originally AUROC everywhere)
**Why it changed.** *Literature* — ROC AUC is equivalent to an average of the actives'
ranks, so early recognition is cancelled by late recognition, which is the wrong summary
for a method whose purpose is the top of a list (Truchon & Bayly, *J Chem Inf Model* 2007,
47:488–508;
[overview](https://www.researchgate.net/publication/6517236_Evaluating_Virtual_Screening_Methods_Good_and_Bad_Metrics_for_the_Early_Recognition_Problem)).
BEDROC at α = 20 puts 80% of the weight in the top 8%.
**Scope.** Applied to the library ranking and E3 only; **not** to E1 (D6.4).
**Sample-size caveat, acted on.** A **permutation null is computed rather than assumed**,
and BEDROC is *not reported* for benchmark-sized contrasts.

**Attribution corrected 2026-10-22, on reading the primary paper.** This entry previously
said BEDROC's variance "is only controlled near ~50 actives". Truchon & Bayly make no such
threshold claim — n = 50 is an illustrative figure in their worked example. What they
actually give is stronger and usable directly:

- An analytical upper bound on the metric's standard deviation,
  **STD_max = 1 / sqrt(8n)**, offered explicitly "for planning evaluations". For this
  study: n = 24 reference ligands → 0.072; n = 7 endogenous → 0.134. The observed BEDROC
  values sit ~3.7–3.8 of those units above their permutation-null medians.
- A **comparability condition, α·R_a ≪ 1** (R_a = n/N), required for BEDROC to be free of
  the saturation effect and for cross-study comparison to be meaningful. This is the
  precise version of what this entry previously hand-waved as "all weight falls in the
  top ~10 ranks":

| contrast | n / N | α·R_a | verdict |
|---|---|---|---|
| library, 24 ref ligands | 24 / 14,595 | 0.033 | ✓ reportable |
| library, 7 endogenous | 7 / 14,595 | 0.010 | ✓ reportable |
| E3, unmatched decoys | 36 / 121 | 5.95 | ✗ saturated |
| E3, aromatic-matched | 36 / 269 | 2.68 | ✗ saturated |

So the decision not to report BEDROC for E3 was right, but for a reason I can now state
exactly rather than approximately. Truchon & Bayly additionally criticise EF for not
discriminating ranking quality *before* the threshold and for having high variance with a
maximum that depends on R_a — so the EF values reported alongside BEDROC in §9.1 of the
study are secondary, not primary.

### D6.9 Effect sizes reported in units of the measured seed jitter — HELD
**Decision.** Every contrast reports `margin / S1_median_abs_delta` alongside its AUROC,
and is flagged **noise-floor limited** below 2×.
**Why.** An AUROC says how consistently two groups are ordered and nothing about how far
apart. This project was misled by that gap twice (D6.1). Measured at the gated
configuration: S2 5.9×, E1 3.2×, E3 9.0×, E5 3.1×, **E7 −1.4× (below the floor, wrong
direction)**.
**Status.** A reporting change only, applied to numbers already computed; it chooses
nothing.

### D6.10 Tiers anchored to a decoy background, not library deciles — HELD
**Rationale.** Library deciles describe the conformer seed as much as the peptide, because
the library score distribution is compressed relative to the seed jitter. Every tier edge
is additionally checked against the measured jitter and reported as *unsafe* if closer to
its neighbour than the jitter.
**Outcome.** All edges cleared the jitter — and the tiers were still wrong, for the
unrelated reason in D6.6. A passing safety check on a badly posed comparison is still a
badly posed comparison.

### D6.11 The free-N-terminal-Tyr rule kept separate from the score — HELD
Never folded into the ranking; carried as its own column so the contribution of the free
one-line sequence rule stays separable from the 3D matching's.

### D6.12 Consensus over two templates is rank fusion, not score averaging — OPEN, pre-registered
**Decision.** The two-template consensus is the **mean rank** across templates (primary),
with per-peptide **max of min-max-scaled scores** reported secondary.
**Rationale.** The two templates' score distributions are not on a common scale: Tversky
reference similarity is normalised by the reference, and template A carries 27 features
from 8 residues against template B's complete 4-residue ligand. Averaging raw scores would
silently weight whichever template has the wider spread. Rank fusion needs no calibration.
Max-of-two is reported because it encodes the biologically motivated reading — a dietary
peptide need only mimic *one* natural binding mode — but it is secondary because it is
upward-biased by construction (two chances to score well) and the antagonists receive the
same bias, so it does not cancel.
**Evidence.** *Literature* for rank fusion in consensus virtual screening; untested here.

### D6.13 Template-B statistics must reuse the existing E8 pools and seed-averaging — HELD
**Decision.** `24_dual_template_consensus.py` imports its pool predicates and statistics
from `23_dietary_screen.py` instead of reimplementing them, and reports the seed-averaged
AUROC as primary.
**Rationale.** A near-duplicate of a statistic is how two numbers that should agree stop
agreeing. E8 does not score all 20 dietary peptides: the published **0.917** is the
μ-dominant-assay pool, **12 agonists against 3 antagonists**, seed-averaged. A first draft
of the new script computed the AUROC over all 20 peptides and 6 antagonists and returned
0.774 for template A — a number that looks comparable to 0.917 and is not.
**Evidence.** *Measured* — after reusing the pool definitions, the new script reproduces
template A's published 0.917 (μ pool, n=12/3, exact p=0.0154) exactly, which is what makes
T1 a comparison rather than a coincidence.
**Standing caveat.** Every E8 number in this project rests on **3 antagonists**. The exact
one-sided permutation p and the full leave-one-antagonist-out range are reported for each,
and a result whose LOO range straddles 0.5 is not claimed whatever its point estimate.

---

## 7. Infrastructure

### D7.1 Everything runs as a SLURM batch job — HELD
A benchmark was once run on the ORCA login node for 12+ minutes. Killed and resubmitted.
No analysis step, including pure re-analysis of score tables, now runs interactively.

### D7.2 Block-parallel arrays with per-block output — HELD
Each block writes its own `.tsv` and `scores.tab`, so a failed block can be rerun alone.
**Measured benefit.** Two of 300 library blocks "failed"; both had received zero peptides
(49 × 298 > 14,595) and the empty-table guard refused to write a header-only file. The
distinction between an empty block and a lost block was visible immediately.

### D7.3 Refuse to write an empty score table — HELD
A zero-row output that looks like a successful run is worse than a failure.

### D7.4 Trajectories and score tables stay on scratch — HELD
Only derived, human-readable results are committed. Grant-internal documents are excluded
via `.git/info/exclude`.

---

## 8. Where this sits against the published work

Calibration, so the numbers above can be read against the field rather than in isolation.

| | AUROC |
|---|---|
| This screen, library retrieval, **uncontrolled** | 0.915 |
| Best ML model for MOR ligands, held-out split ([µORScreen, 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13488679/)) | 0.714 |
| Shape-overlap screening, mean over DUD-E | 0.58 |
| Shape-overlap screening, mean over LIT-PCBA (debiased) | 0.54 |
| **This screen, E3, aromatic-matched** | **0.490** |

Two things follow. First, an uncontrolled 0.915 exceeding the best published MOR model by
0.20 is itself evidence of residual bias rather than of a better method. Second, the
matched 0.490 is *close to* what shape-based methods generally achieve on debiased
benchmarks — so this is a field-wide difficulty showing up here, not a defect unique to
this implementation. Neither observation makes the matched result good.

On template provenance: no published MOR pharmacophore appears to be derived from a
cryo-EM peptide complex. The existing models are built from the 4DKL X-ray structure with
a morphinan antagonist
([PMC6865014](https://pmc.ncbi.nlm.nih.gov/articles/PMC6865014/),
[PMC4757823](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4757823/)) or are ligand-based
on morphinans
([model](https://link.springer.com/article/10.1186/1471-2210-10-S1-A4)). A peptide-complex
template is therefore not redundant with the literature — which is a point in favour of
the approach, independent of how this validation turned out.

On validation rigour: the closest methodological precedent — a receptor-based
pharmacophore from 4DKL screening ~198,000 tetrapeptides, 28,070 passing the filter, 3
assayed, 1 active — validated its pharmacophore on **two compounds** (Tyr-Gly-Gly-Phe
retained, rubiscolin-5 discarded) with no decoy set, no enrichment factor and no
retrospective screen
([PMC6865014](https://pmc.ncbi.nlm.nih.gov/articles/PMC6865014/)). The validation in this
repository is considerably more thorough than the published comparator, and it is the
thoroughness that produced the negative result. That asymmetry is worth stating plainly in
any write-up: a less careful version of this same pipeline would have reported a clean
success.

---

## 9. Open decisions, in priority order

| # | Change | Basis | Cost |
|---|---|---|---|
| 1 | **Add receptor-derived exclusion volumes** (D2.4) | literature; mechanistically targets the confound | one template rebuild + one benchmark rescreen |
| 2 | **Make `POSC` a mandatory feature** (D2.5) | literature; filter-then-rank | rescreen with `-o`; see cost correction in D2.5 |
| 3 | **Expand the active set to power the 2- and 3-aromatic strata** (D6.6) | measured; 26 of 36 actives in a stratum at AUROC 0.513 | curation, not compute |
| 4 | Adopt AVE-style debiasing for decoy selection | literature; [LIT-PCBA](https://www.researchgate.net/publication/340614114_LIT-PCBA_An_unbiased_dataset_for_machine_learning_and_virtual_screening) uses asymmetric validation embedding for exactly this | moderate |
| 5 | **Retest the Tversky asymmetry once exclusion volumes are in** (D5.3) | literature; TVERSKY_DB was tested with the "no growth" constraint absent from the template, so metric and constraint were confounded | re-analysis of the D2.4 rescreen |
| 6 | **Test whether the ensemble-size penalty is `max`-specific** (D4.2) | literature; McNutt's false-positive mechanism should apply to `max` and not to `mean`, which would reconcile 800 conformers with a <25 recommendation | re-analysis only, data already on disk |
| 7 | **Run the dual-template screen** (D2.7–D2.9, D6.12) — built and pre-registered, not yet scored | measured + literature; β-endorphin is the 4th-best of five endogenous μOR ligands and 42× weaker than endomorphin-1, and 8F7Q is a truncation where 8F7R is not | one template build + 4 dietary runs (`jobs/dual_template.sbatch`, minutes); the E4′ percentile additionally needs a full library screen (`jobs/dual_template_library.sbatch`, ~52 CPU-h), which should **not** be run until T3 shows the templates are non-redundant |

Items 1 and 2 are the ones that could change the study's conclusion. Item 3 is what would
make the conclusion trustworthy either way.
