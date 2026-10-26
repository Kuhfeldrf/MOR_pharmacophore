# Opioid-motif atlas over the milk proteome

Which peptides carrying an opioid motif are actually released from milk proteins, with the
motif at their own N-terminus? Built by
[`../scripts/25_milk_motif_atlas.py`](../scripts/25_milk_motif_atlas.py); outputs in
[`../data/milk_atlas/`](../data/milk_atlas/). Run 2026-10-26, ~4 s, no cluster needed.

This replaces the sequence-level pre-screen that README referred to in the sibling
`milk-opioid-motif-atlas` repository, which is not present on this machine. The two agree
closely on the number that matters — that work reported 11 motif sites released as
N-termini, this one finds **12**.

---

## 1. Why released-at-the-N-terminus is the only population that counts

The positional rule (`mor_pharmacophore.nterm_rule`, D6.11): opioid binding needs a **free
α-amine** on the anchoring Tyr, because that amine makes the D149(3.32) salt bridge the
source paper calls "highly conserved" and whose mutation "nearly abolished downstream
signaling" ([`structural_basis.md`](structural_basis.md) §3). A peptide that merely
*contains* `YPFP` has that amine locked in a peptide bond and cannot anchor.

So the scan produces three populations, and only the first is a candidate:

| Population | Meaning | n |
|---|---|---|
| **`released_nterminal`** | observed peptide **starts exactly at** a motif site | **35 peptides, 12 sites** |
| `spanning_blocked` | observed peptide contains a site but starts upstream | **1,767 peptides** |
| not observed | motif site never seen in any observed peptide | remainder of 235 |

**The positional rule discards 98% of motif-bearing observed peptides** (1,767 spanning
against 35 released). That ratio is the quantitative case for the rule: motif presence is
nearly uninformative on its own, which is exactly what C4 measured on `YPFP` and what no
sequence-similarity method captures.

## 2. Inputs

| | |
|---|---|
| Proteomes | 429 human milk proteins + 494 bovine (UniProt FASTA) |
| Observed peptides | **14,571**, MS-identified in infant gastric and intestinal digesta — `data/library/peptides.tsv` (13,043 human milk, 1,071 bovine+human, 457 bovine) |
| Mapping | exact substring match, peptide → parent protein and 1-based start |
| Unmapped | 18 of 14,571 (0.1%) — reported, not dropped silently |

The observed peptides are **experimentally derived**, not an in-silico digest: they are
real MS identifications from digesta, which is what makes "could be released" a
measurement here rather than a prediction.

## 3. Motif definitions, and where they come from

Derived from this repository's own reference peptides, not invented for this script.

**Tier 1 — canonical** (235 sites; the only tier that becomes a screening target):

| Class | Regex | Reference peptides that motivate it | Sites |
|---|---|---|---|
| `casomorphin_YP` | `YP[FWYVLIM]` | β-casomorphins `YPFP`/`YPFPG`/`YPFPGPI`, human `YPFV`, soymorphin `YPFVVN`, neocasomorphin `YPVEPF`, endomorphin-1 `YPWF`, casoxin B `YPYY` | 159 |
| `lactorphin_YxxF` | `Y[GL][LGSF][FY]` | α-lactorphin `YGLF`, β-lactorphin `YLLF`, lactoferroxin A `YLGSGY` | 51 |
| `casoxin_YPxY` | `YP[STAG][YFW]` | casoxin A `YPSYGLN`, casoxin-6 `SRYPSY` | 14 |
| `message_GGF` | `[YF]GGF` | enkephalins `YGGFM`/`YGGFL`, β-endorphin, dynorphin A, N/OFQ `FGGF` (F anchor) | 11 |

**Tier 2 — aromatic pair** (`[YF].{1,2}[FWY]`, 4,437 sites): the mechanistic
generalisation. The source paper finds that YGGF, YPW, Y[D-Ala]F and FGGF all occupy one
"highly aromatic pocket" (§4), i.e. what the receptor reads is an aromatic anchor plus a
second aromatic two or three residues later. **Reported for breadth, never screened.** At
4,437 sites it is permissive by construction and a tier-2 hit is a hypothesis, not a
prediction.

Only **16 of the 235** tier-1 sites fall in a major milk protein. A motif in a trace
signalling protein is not a dietary opioid source in any practical sense, so the
`major_milk_protein` flag is carried on every row.

## 4. Result: the atlas recovers the known dietary opioids unprompted

The 35 released peptides cluster at 12 sites, in **β-casein (24), αs1-casein (11),
α-lactalbumin (2), MUC1 (1)** and two incidental actin hits. Six of the twelve sites are
*known* opioid peptides with published activity, recovered without the script being told
they exist:

| Site | Released peptide (shortest observed) | Known as | Published activity |
|---|---|---|---|
| β-casein P02666 @75 | **`YPFPGPI`** | **β-casomorphin-7** | GPI 57 / 5.14 / 3.2 µM; μ K_D 2 µM |
| β-casein P05814 @66 | `YPFVEPIP` | human β-casomorphin (`YPFVE`…) | μ K_D 0.4–1.6 µM |
| β-casein P02666 @129 | `YPVEPFTES` | neocasomorphin-6 (`YPVEPF`) | GPI 59 µM; μ IC50 92 µM |
| α-lactalbumin P00709 @69 | `YGLFQISNK` | α-lactorphin (`YGLF`) | in the dietary set |
| αs1-casein P02662 @106 | `YLGYLEQLL` | αs1-casein f91–96 (`YLGYLE`) | IC50 5.2 µM (Loukas 1983) |
| αs1-casein P47710-3 @137 | `YPYAVWYYPQ` | — | novel |

`YPFPGPI` is β-casomorphin-7 exactly, observed as a released N-terminus. That the scan
finds the single best-characterised dietary opioid peptide, at its correct position in
bovine β-casein, with a free Tyr, is the positive control for the whole approach.

**The αs1-casein hit is worth a second look.** The canonical exorphin in the dietary set is
`RYLGYLE` (f90–96), which begins with **Arg** — it carries the motif but *fails the
positional rule*. What is actually released and observed here is `YLGYLEQLL` starting at
f91, i.e. the **free-Tyr form**, which satisfies the rule. Loukas 1983 measured that form
separately (`YLGYLE`, IC50 5.2 µM). So the atlas did not merely re-find a known peptide; it
identified the form of it that the binding model says should work.

### Genuinely novel candidates

Five sites, released and motif-bearing, with no published opioid activity:

| Site | Peptide | Class | Note |
|---|---|---|---|
| αs1-casein P47710-3 @137 | `YPYAVWYYPQ` | `casomorphin_YP` | **Most interesting.** Human αs1-casein, 4 aromatics in 10 residues (Y,W,Y,Y). High aromatic content is the known confound (ρ = +0.909), so this one needs the exclusion-volume template specifically — see §5. |
| β-casein P05814 @56 | `YPSFQPQPLI` | `casoxin_YPxY` | human β-casein; `YPSF` resembles casoxin A's `YPSY` |
| β-casein P05814 @211 | `YPVTQPLAPVH` | `casomorphin_YP` | human β-casein C-terminal region |
| MUC1 P15941 @1209 | `YPTYHTHGRYVPP` | `casoxin_YPxY` | mucin; 4 aromatics, same caveat as above |
| ACTA1 P68133 @71 | `YPIEHGIITNWDDMEK` | `casomorphin_YP` | actin — not a milk protein; almost certainly an incidental match |

## 5. How to read the screen when it comes back

The screening target table is `data/milk_atlas/milk_atlas_targets.tsv` (35 peptides),
screened against both templates by `jobs/milk_atlas_screen.sbatch`.

**Check the positive controls first.** β-casomorphin-7, the human β-casomorphin series,
neocasomorphin-6, α-lactorphin and the αs1 exorphin are all in this set with published
activity. If the screen does not rank them above the novel sites, the result says something
about the screen, not about the novel sites.

**Two novel candidates are aromatic-rich and must not be believed on score alone.**
`YPYAVWYYPQ` (4 aromatics) and `YPTYHTHGRYVPP` (4 aromatics) sit exactly where this
project's documented confound lives: score correlates with aromatic count at ρ = +0.909,
and a His-rich repeat once outscored every endogenous ligand
([`optimisation_study.md`](optimisation_study.md) §7.2). Exclusion volumes are the fix and
both templates carry them (D2.4), but for these two peptides the exclusion-volume
comparison — scored with against scored without — is the result, not the raw rank.

**What this cannot tell you.** Nothing here measures binding, bioavailability or survival
of further digestion. The atlas selects what is worth scoring; the screen ranks geometric
match; neither is an affinity. And per `structural_basis.md` §5 the templates are
pan-opioid message matchers, so none of this supports a μ-selectivity claim.

## 6. Reproducing it

```bash
python scripts/25_milk_motif_atlas.py \
    --fasta human=/path/to/HumanMilkProteinDatabase_v2.fasta \
    --fasta bovine=/path/to/CowMilkProteinDatabase_v2.fasta \
    --library data/library/peptides.tsv --out data/milk_atlas
sbatch jobs/milk_atlas_screen.sbatch     # needs both templates built
```

The FASTA paths are inputs rather than committed data because they are shared reference
databases maintained outside this repository; `summary.json` records the protein counts so
a rerun against a different release is detectable.
