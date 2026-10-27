# Extent 9 re-measures extent 8: the method saturates at the message epitope

**Run 2026-10-27**, jobs 242786 (screen) and 242787 (gate). Pre-registered in
[`post_t3_design.md`](post_t3_design.md) §4 before these runs existed. Machine-readable:
[`extent9_gates.json`](extent9_gates.json).

## Result

| # | Test | Threshold, fixed in advance | Measured | |
|---|---|---|---|---|
| **W1** | seed-to-seed Spearman at window 9 | ≥ 0.98 | **+0.9774** | **fail** |
| **W2** | extent-9 template carries a feature extent 8 lacks | POSC count must rise | **2 vs 1** (30 vs 27 features) | **pass** |
| **W3** | Spearman vs extent 8 | < 0.9 | **+0.9880** | **fail** |
| **W4** | ρ(score, aromatic count) | \|ρ\| < 0.909 | **−0.2312** | **pass** |

Gate exit 1, as designed. Job 242787 is recorded by SLURM as `FAILED` because the gate
returns non-zero when the design does not survive; that is the mechanism working, not a
crash.

## What this means

**The feature is in the template, and it changes nothing.** W2 passed exactly as
predicted: the extent-9 template gains one POSC and two HDON features over extent 8 —
Lys9's side-chain amine and its backbone donors, the amine that salt-bridges Asp218 in
ECL2 at 3.1 Å ([`8f7q_contacts.tsv`](8f7q_contacts.tsv)). The geometry that
[`post_t3_design.md`](post_t3_design.md) §2 showed was reachable is now demonstrably
represented in the pharmacophore.

And the score does not care. ρ(extent 9, extent 8) = **+0.988** over all 20 peptides.
That is not a near miss on the 0.9 threshold; it is the two templates producing the same
ranking. Adding the ECL2 contact to the reference changes the geometry of the template
without changing what the scoring function sees.

**Both available directions are now closed, empirically.**

| Direction tried | Result |
|---|---|
| Add a second, different short template (8F7R/endomorphin-1) | ρ = **0.917** — [`dual_template_result.md`](dual_template_result.md) |
| Extend the one long-peptide template to its last ECL2 contact | ρ = **0.988** — this document |

The first failure could be blamed on the choice of second structure. The second cannot:
extent 9 is the same ligand, the same deposition, the same receptor copy, one residue
longer, with the new feature verified present. **The method saturates at the message
epitope.** What the Align-it TVERSKY_REF similarity reads is the N-terminal aromatic
anchor; address-region features can be added to the template and they do not move the
ranking.

This is the mechanism `post_t3_design.md` §1 predicted from the literature — truncation
to the message epitope compresses score separation because the shared anchor is what is
kept — now measured directly rather than inferred from a correlation between two
different templates.

## W1, honestly

W1 failed at 0.9774 against a threshold of 0.98. The margin is small, window 8's measured
value is 0.983, and it would be easy to argue 0.9774 is "near enough". **It was
pre-registered at 0.98 and it missed**; relitigating the threshold after seeing the number
is precisely what pre-registration exists to prevent, so it is recorded as a failure.

W1 is also moot. Had W3 passed, a marginal W1 would have justified one more experiment —
re-run at 1600 conformers and see whether seed agreement recovers, since 31–37 rotatable
bonds at 800 conformers is the suspected cause. With W3 at 0.988 there is nothing to
recover: better sampling of a measurement that duplicates extent 8 is still a duplicate.

## One caveat on the inputs

These runs used casoxin A as `YPSYGLNY`, which
[`evaluation_set_provenance.md`](evaluation_set_provenance.md) establishes is one residue
too long. W1, W3 and W4 are rank statistics over all 20 peptides, so a single peptide's
sequence cannot account for ρ = 0.988 — the conclusion is not sensitive to it. W5 was not
computed here and is not a gate; every AUROC remains subject to that correction and to
the three-antagonist power ceiling in [`power_analysis.md`](power_analysis.md).

## Where this leaves the study

Not with another template. Three independent lines now say the same thing — the
statistical one (no template comparison is adjudicable at three antagonists), the
mechanistic one (message truncation compresses separation), and this empirical one (the
ECL2 feature is inert to the score). The structural inventory adds that there is nothing
left to try: 8F7Q is the only MOR structure resolving past ligand residue 5, and extent 9
exhausts its selectivity-relevant contacts.

Two avenues remain, and both are stated in `post_t3_design.md` §6 rather than being new:

- **Scoring, not geometry.** The GPCR literature that succeeds at functional
  discrimination pairs shape with interaction-pattern scoring rather than using shape
  alone (de Graaf & Rognan 2008, [10.1021/jm800710x](https://doi.org/10.1021/jm800710x);
  Kooistra et al. 2016, [10.1038/srep28288](https://doi.org/10.1038/srep28288)). This
  result — a template feature that is present but inert — is a direct argument for that
  pairing, and `src/mor_pharmacophore/plip_step.py` already profiles receptor–peptide
  interactions.
- **Data, not method.** Eight to twelve antagonists would make an AUROC claim properly
  powered, and [`evaluation_set_provenance.md`](evaluation_set_provenance.md) §5
  establishes that this requires new measurement or primary-literature curation, because
  the public databases do not contain these peptides.

Nothing in either avenue is licensed by this document. What is settled is narrower and
worth stating plainly: **extending the template's reach does not work, and the reason is
not that the geometry was out of reach.**
