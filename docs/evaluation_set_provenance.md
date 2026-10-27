# Provenance audit of the dietary evaluation set

**Written 2026-10-27.** Three findings about the three antagonists that every AUROC in
this project rests on, each checked against the primary source rather than against a
review or a compilation. One of them is a sequence error and is **deliberately not fixed
here** — see §4.

Prompted by [`power_analysis.md`](power_analysis.md), which showed the evaluation set is
the binding constraint on the whole study. If three antagonists decide everything, their
provenance deserves more scrutiny than the rest of the data put together.

---

## 1. Casoxin A's sequence is wrong: `YPSYGLNY` should be `YPSYGLN`

**Primary source.** Chiba H, Tani F, Yoshikawa M, "Opioid antagonist peptides derived
from κ-casein", *J Dairy Res* 1989;56(3):363–366,
[10.1017/s0022029900028818](https://doi.org/10.1017/s0022029900028818). The abstract
states: *"A synthetic bovine κ-casein peptide (35-41), Tyr-Pro-Ser-Tyr-Gly-Leu-Asn
(casoxin A)"* — **seven residues**, `YPSYGLN`.

**Independent confirmation.** UniProt P02668 (bovine CSN3) has a 21-residue signal
peptide, so mature numbering 35–41 is precursor 56–62, and that heptamer is exactly
`YPSYGLN`. Precursor position 63 is another Tyr, so `YPSYGLNY` is mature **35–42** — the
published peptide plus one residue.

**This repository already has it right elsewhere.** [`milk_atlas.md`](milk_atlas.md) §
lists "casoxin A `YPSYGLN`". The evaluation set disagrees with the project's own
documentation, which is what identifies it as an error rather than a recorded choice.

**Why it matters more than one residue.** The extra residue is a Tyr, so casoxin A's
aromatic count is 3 where it should be 2, and the antagonist aromatic mean is **2.67 where
it should be 2.33** (agonists 2.00). That is the project's documented confound direction:
`23_dietary_screen.py` argues E8 is unusually robust *because* the antagonists are more
aromatic than the agonists, so the score cannot be passing by counting aromatics. The
correction weakens that argument without overturning it — 2.33 is still above 2.00.

Chemically the correction is clean: at window 8 a 7-mer is untruncated and keeps a free
acid, exactly as the 8-mer did, so no amidation decision changes (D3.5).

## 2. The antagonist potencies are threshold concentrations, not IC50 values

The same abstract: casoxin A *"showed opioid antagonist activity **at** 200 µM in the
guinea pig ileum assay"*; casoxin B *"also showed opioid antagonist activity **at** 100
µM"*; casoxin C *"was **active at** 5 µM"*. These are the lowest concentrations at which
activity was observed. **The paper reports no IC50, no Ki and no dose-response fit, runs
no binding assay, and resolves no receptor subtype.**

The data file records them as `ic50_uM` = 200 / 100 / 5 and derives `pIC50` = 3.699 /
4.000 / 5.301. A threshold and a half-maximal inhibitory concentration are different
quantities, and the pIC50 column implies a measurement that was never made.

**What this does and does not touch.** The AUROC tests (E8, T1, T2) use only the
Agonist/Antagonist *label*, so they are unaffected. The potency-correlation tests (E9,
T4) use the values: `E9_mu_assays_only` at n = 15 includes these three pseudo-pIC50s,
while `E9_mu_agonists_only` at n = 12 does not. The n = 12 variant is the one that rests
entirely on measured IC50s.

**Decided:** annotated in place; the values and the tests are unchanged, so existing
numbers stay comparable. The `assay` field now reads `GPI (threshold active concn, not
IC50)` for the three casoxins. Note separately that assigning `mu_dominant` is an
inference from the tissue — guinea-pig ileum is µ-dominant as a preparation — and not a
subtype determination reported by Chiba et al.

## 3. The citation pointed at the wrong paper — fixed

Both [`../data/dietary/dietary_peptides.tsv`](../data/dietary/dietary_peptides.tsv) and
[`../data/reference/reference_peptides.tsv`](../data/reference/reference_peptides.tsv)
carried `PMID:2760234` for Chiba et al. 1989. That PMID is *"Comparison of intralobar
non-parenchyma, subcapsular non-parenchyma, and liver capsule thickness"*, Ryoo & Buschmann,
*J Clin Pathol* 1989 — an unrelated paper. The correct PMID is **2760302**. Corrected in
both files; it changes no result.

## 4. What was deliberately not changed, and why

The sequence in §1 is **not** corrected in this commit. `reference_peptides.tsv` records
its provenance as *"Table1_Updated_GPI_IC50.xlsx (NIH R21 … Jun 2026); values transcribed
verbatim 2026-09-24"*, and both the extra residue and the IC50 framing of §2 are
consistent with having been transcribed faithfully from that spreadsheet. The error is
therefore **upstream of this repository**, and it may affect other work drawing on the
same file.

Correcting it here would also move every AUROC, T3's ρ, T4/E9 and T6 at once, and
invalidate the runs in flight. So the evidence is recorded first and the spreadsheet is to
be checked with its author before the evaluation set changes. The files still carry
`YPSYGLNY`; this document is the record that they are known to be wrong.

**Rows that would change, and only these five:** `data/targets.tsv` (`ref__casoxin_A`),
`data/reference/reference_peptides.tsv` (`casoxin_A`), `data/controls/control_sets.json`,
`data/dietary/dietary_targets.tsv` and `data/dietary/dietary_peptides.tsv`
(`diet_casoxin_a`). Everything else matching `YPSYGLNY` is correct and must be left alone:
the library hits are genuine digestion products (`LSRYPSYGLNY`,
`SRYPSYGLNYYQQKPVAL`), and `data/milk_atlas/spanning_blocked.tsv` holds an 8-residue
*window* at that site whose own annotation column already reads "casoxin A YPSYGLN".

---

## 5. The power ceiling cannot be lifted by querying databases

[`power_analysis.md`](power_analysis.md) concluded that only more antagonists can buy
precision. A survey of the public bioactivity sources (queried 2026-10-27) says that route
is largely closed.

**This evaluation set is not in the bioactivity databases.** Matching all 24 non-library
rows of `data/targets.tsv` by InChIKey against human and rat ChEMBL (CHEMBL233/CHEMBL270),
BindingDB (P35372/P33535) and GtoPdb: 9 of 24 are present. Of the 17 food-derived
positives, only **two** — β-casomorphin-5 and -7, bovine — and only at **rat** MOR. **All
four casoxins and both lactoferroxins, i.e. every antagonist in this set, are absent from
ChEMBL and BindingDB entirely**, as are casomorphin/casoxin/lactoferroxin/lactorphin/
rubiscolin/soymorphin/exorphin as name searches. The granularity constraint is structural,
not a curation shortfall that querying harder would fix.

**Antagonist peptides do exist, but not on a poolable scale.** Across all MOR orthologs,
ChEMBL holds ~101 distinct peptide molecules with a true antagonist potency
(Ke/Kb/pA2/pKb) — but ~92% are guinea-pig-ileum Ke values from unnamed analogue series,
and only three carry a preferred name. GtoPdb's expert-curated MOR antagonist peptide set
is **exactly three**: CTOP, CTAP, zyklophin. A Ke is an agonist-independent equilibrium
constant derived from a shift in an agonist curve; an agonist IC50 is a system-dependent
midpoint (Neubig et al. 2003, *Pharmacol Rev* 55:597,
[10.1124/pr.55.4.4](https://doi.org/10.1124/pr.55.4.4), which also notes pA2 = pKB only
when the Schild slope is indistinguishable from 1). They are not on a common scale, so
adding them to this set would trade a power problem for a commensurability problem.

**Pooling noise is itself contested, and both sides should be cited.** Landrum & Riniker
2024, *JCIM* 64:1560, [10.1021/acs.jcim.4c00049](https://doi.org/10.1021/acs.jcim.4c00049),
find that for compounds measured against the same target in different ChEMBL assays ~65%
of paired points differ by >0.3 log and 27% by >1 log, with Kendall's τ ≈ 0.51.
Kalliokoski et al. 2013, *PLoS ONE* 8:e61007, measure σ(pIC50) ≈ 0.68 in public ChEMBL
against ≈0.22 intra-lab, yet conclude mixing adds "only a moderate amount of noise". The
contradiction is real and worth reporting as such. Either way a per-point uncertainty near
0.5 log units is the scale against which this project's 0.028 AUROC granule should be
read.

**There is no peptide benchmark to fall back on.** The only reusable µOR virtual-screening
benchmark is GPCR-Bench's OPRM set (82 actives / 4,800 decoys, built on PDB 4DKL, the
*mouse* receptor; Weiss et al. 2016,
[10.1021/acs.jcim.5b00660](https://doi.org/10.1021/acs.jcim.5b00660)), and its
construction rule — MW 250–500, 0–7 rotatable bonds — **excludes every opioid peptide,
DAMGO included, by design**. DUD-E has no opioid target and caps ligands at MW < 600;
LIT-PCBA covers OPRK1 only. LEADS-PEP, PepBDB and Propedia are pose-prediction sets with
no screening library and no GPCR.

**One operational constraint worth recording:** BIOPEP-UWM, the only food-peptide
compilation with usable numbers (157 opioid records, of which 7 are coded
`opioid antagonist`), prohibits automated extraction in its
[terms of use](https://biochemia.uwm.edu.pl/en/terms-of-use-2/) and requires a written
agreement for a machine-readable copy. Manual queries only — it must not be scraped.

### Consequence

Expanding the antagonist class to the 8–12 needed for a properly powered AUROC
([`power_analysis.md`](power_analysis.md)) would require **new measurements, or primary
literature curation, not a database query**. Until then the conclusion of
[`post_t3_design.md`](post_t3_design.md) §5 stands and hardens: no AUROC-based claim in
this project should be strengthened, including the exclusion-volume result.
