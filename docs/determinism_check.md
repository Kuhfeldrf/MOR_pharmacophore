# The screen is deterministic; the discrepancy was the environment

**Run 2026-10-28** (job `242919`, superseding `242912` — see §4).
Machine-readable: [`determinism_check.json`](determinism_check.json).

## 1. The question

When template A's dietary runs were regenerated after the scratch allocation was reaped,
27 of 40 peptide-scores reproduced to four decimal places and two did not:
**β-lactorphin +0.0277** and **casoxin A +0.0213**, roughly ten times the documented seed
jitter of 0.00133–0.00225, with conformer counts unchanged
([`dual_template_result.md`](dual_template_result.md), "A reproducibility note on `A8`").

Three candidate causes, needing three different fixes:

| | Cause | Consequence if true |
|---|---|---|
| **C1** | conformer generation is not deterministic given a fixed seed | the seed does not pin the ensemble; every reported score is one draw |
| **C2** | align-it's alignment search is not deterministic given fixed input | same, and the S1 seed-stability gates measure run-to-run noise too |
| **C3** | neither — the two runs used different library builds | the numbers are sound; the provenance was not |

C1 and C2 are testable directly. C3 is not testable retrospectively, and that is itself
the finding: `04_screen.py` recorded **only align-it's version**, so the original run's
RDKit — the library that embeds the conformers — is unrecoverable. The regenerated run's
`block0000.json` has no `environment` key at all, because it predates the fix below.

## 2. C1 and C2: clean

`scripts/31_determinism_check.py` runs each step twice on identical input in one
environment, over the two peptides that moved plus Met-enkephalin, which reproduced
exactly, as a control — a test that calls everything non-deterministic is not measuring
anything. For C2 **one** protonated SDF is scored twice, so a difference cannot be
conformer variation leaking in.

Environment: rdkit 2024.09.6, align-it v1.0.4, Open Babel 3.2.1, Python 3.11.18.

| peptide | `c_term` | C1 conformers | C1 coordinate digest | C2 max \|Δ\| |
|---|---|---|---|---|
| casoxin A | free | 800 = 800 | identical | 0.000000 |
| β-lactorphin | amide | 793 = 793 | identical | 0.000000 |
| Met-enkephalin | free | 797 = 797 | identical | 0.000000 |

C1 hashes every atom coordinate of every conformer at full precision — stricter than the
four-decimal comparison that surfaced the discrepancy. Both are bit-for-bit identical.

## 3. C3 is the surviving explanation, and now has positive evidence

Elimination alone would be weak. The stronger result is that **this run reproduces the
regenerated production run to full precision**, two days later, on all three peptides
including both anomalies:

| peptide | regenerated production run (`p8excl_s61453`) | this run, independent | conformers |
|---|---|---|---|
| β-lactorphin | 0.1904 | 0.190446 | 793 = 793 |
| casoxin A | 0.1348 | 0.134792 | 800 = 800 |
| Met-enkephalin | 0.1847 | 0.184666 | 797 = 797 |

So the pipeline is reproducible *within* this environment, across weeks and separate jobs.
The only difference remaining between the original and regenerated runs is the stack that
executed them, which was rebuilt from scratch in between.

**This is as far as the evidence goes.** C3 is the only surviving candidate and it is
consistent with everything measured, but it is not confirmed: confirming it would need the
original RDKit version, which was never recorded. Stated as elimination plus consistency,
not as a demonstrated cause.

**What it means for the reported numbers.** The regenerated `A8` of 0.944 is a sound
measurement of the current stack, and the pre-registered 0.917 is a sound measurement of
the previous one. Neither is a mistake; they are not the same experiment, and the project
could not tell because it did not write down which stack it was using. Scores are
comparable within an environment and must not be compared across a rebuild.

## 4. Two bugs found in the course of this, both in my own work

**The check built the wrong molecule for β-lactorphin.** Job `242912` hardcoded
`c_term="auto"` while the production targets table marks `diet_beta_lactorphin` as
`amide`, so it embedded a free acid where the screen embeds an amide: **788** conformers
against production's **793**. Casoxin A and Met-enkephalin are `free`, which is what
`auto` resolves to, so those two were right by coincidence. The C1/C2 verdict was
unaffected — identical input twice is identical whatever the molecule — but β-lactorphin
is one of the two peptides the check exists to investigate, and it was investigated as a
molecule the screen never scored. Sequence and `c_term` are now read from the targets
table by `id` with the same normalisation guard `04_screen.py` applies. The 788-conformer
output is kept as `determinism_check.job242912_wrong_cterm.json` rather than deleted,
being the evidence for the bug it was replaced over.

**`environment_record()` raised when align-it was absent**, the opposite of its stated
contract, caught by `tests/test_environment_record.py` on the way in. Each lookup is now
independently guarded: by the time that block is written the screen has already spent its
CPU, so a missing tool costs one field, not the run.

## 5. The fix that matters

Every block JSON now records **rdkit, align-it, obabel, python and platform**. An unknown
version is written as `null` rather than omitted, because an absent key cannot be
distinguished from an older run that never had one.

A seed is not provenance if the stack that consumed it is unrecorded.

## 6. Reproducing it

```bash
sbatch jobs/determinism.sbatch     # idempotent on docs/determinism_check.json
```
