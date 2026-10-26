# Structural basis: what Wang et al. 2023 says, and what it means for this pipeline

**Source.** Wang Y, Zhuang Y, DiBerto JF, Zhou XE, Schmitz GP, Yuan Q, Jain MK, Liu W,
Melcher K, Jiang Y, Roth BL, Xu HE. *Structures of the entire human opioid receptor
family.* **Cell** 186:413–427 (19 Jan 2023). doi:10.1016/j.cell.2022.12.026

This is the paper that deposited the structure this pipeline's template is built from, so
it is the authority on what that structure does and does not contain. Every claim below
is theirs; the commentary under each is the consequence for our design, and is ours.

---

## 1. The depositions, and which one we use

| PDB | Complex | Role here |
|---|---|---|
| **8F7Q** | mOR – **β-endorphin(1–31)** – Gi | **our template source** |
| **8F7R** | mOR – endomorphin-1 – Gi | held out, used only for E4 retrieval |
| 8F7S | dOR – [D-Ala²]deltorphin – Gi | not used |
| 8F7W | kOR – dynorphin A(1–13) – Gi | not used |
| 8F7X | NOPR – N/OFQ – Gi | not used |

Prior mOR structures they compare against, for orientation: 8EF5 (fentanyl), 8EF6
(morphine), 5C1M (BU72), 4DKL (β-FNA, inactive), 8EFQ (DAMGO).

**Consequence.** All existing published mOR pharmacophore models derive from **4DKL**, an
*inactive-state, morphinan-antagonist* structure, or are ligand-based on morphinans. A
template from an *active-state, peptide-bound* structure is therefore not redundant with
the literature. That is a genuine argument for this approach, independent of how our
validation turned out.

## 2. 8F7Q resolves 21 peptide residues, not 5 — and we used 5

Their description: residues **M5–A21** are visible in the density, divided into the
hydrophilic **M5–T12** stretch, which "interact[s] mainly with ECL2" and serves as a
linker, and **P13–A21**, which "formed an α-helix turn at the extracellular end of the
TM1/2 cleft". The deposited chain P is 21 residues: `YGGFMTSEKSQTPLVTLFKNA`.

The inherited `data/8f7q_complex.pdb` carried **5**, plus 287 of 292 receptor residues.
The template therefore encoded the message motif alone and contained no ECL2 contact at
all. Corrected in `scripts/22_prepare_complex.py`, which rebuilds from the deposited mmCIF
and records provenance; see [`../DECISIONS.md`](../DECISIONS.md) D2.6.

Rebuilding at three extents reproduces their interaction map:

| extent | ligand | features | PLIP salt bridges |
|---|---|---|---|
| 5 | YGGFM | 16 | **none** |
| **8** | YGGFMTSE | **27** | **ARG213** |
| 21 | YGGFMTSEKSQTPLVTLFKNA | 62 | ARG213, ASP218 |

## 3. The specific interactions, and whether our template has them

Their Figure 2 interaction set for β-endorphin–mOR:

| Interaction | In our template? |
|---|---|
| **Y1 primary amine – D149(3.32) salt bridge, ≤4.5 Å**; D149 stabilised by H-bonds to Q126(2.60) and Y328(7.43). "Highly conserved"; D149³·³²A, Q126²·⁶⁰A or Y328⁷·⁴³A "nearly abolished downstream signaling" | `POSC` feature **yes** — but only because Align-it scores the *protonated* ligand (D2.3). **PLIP never reports this salt bridge at any extent**, because PLIP runs on the unprotonated complex and types the free N-terminal amine as neutral. Our recorded interaction evidence understates the single most important contact in the family. |
| G2 main-chain carbonyl – W320(7.35) H-bond | implicit in HACC features |
| **E8 side chain – R213(ECL2)/Q214(ECL2) H-bond** | **only at extent ≥ 8.** Absent from the template we actually screened with. |
| **K9 side chain – D218(ECL2) H-bond** | **only at extent 21.** Absent. |
| Y1 side chain hydrophobic network with Y150 and neighbours | `LIPO`/`AROM` yes |
| **Y1 participates in a *water-mediated* polar network with Y3.33 and H6.52** | **no.** The template contains no water positions, so a bridging-water interaction cannot be represented. |

Mutating the ECL2-contacting residues "showed compromised mOR activation by
β-endorphin", so these are load-bearing, not incidental.

## 4. The recognition chamber is "highly aromatic" — and it is bounded

Their conserved-chamber finding, which is the most important single paragraph for this
project:

> the YGGF motif from β-endorphin and dynorphin "adopt a similar conformation and occupy a
> **highly aromatic pocket** with nearly the same residue composition from TM2/3/5/6/7"…
> "we observed that **YPW of endomorphin, Y[D-Ala]F of deltorphin, and FGGF of NFQ** all
> share such highly aromatic pocket"… "the synthesized **YGGF peptide alone can activate
> all the four types of ORs**".

**This reframes our central problem.** The ρ = +0.909 correlation between our score and
aromatic-residue count
([`optimisation_study.md`](optimisation_study.md) §7.2) is **not purely an artefact** —
aromatic complementarity is the physical basis of opioid-motif recognition, so a method
that rewards aromatic groups is reading something real.

What makes our version wrong is that the pocket is a **chamber of fixed size**. Aromatic
groups beyond what it accommodates cannot contribute to binding, and may prevent it.
Our score had no upper bound: it rewarded aromatic coverage without limit, which is why
a His-rich repeat (`EFHSHEFH`, 5 aromatics) outscored every endogenous ligand (2–3).

The correct statement of the defect is therefore sharper than "the score is confounded by
composition": **the score captures the pocket's aromatic preference and omits its steric
bound.** Exclusion volumes are exactly the representation of that bound, which is why
D2.4 is the right fix and why changing the similarity metric was not
(§8.1: every metric that weakened the correlation weakened E1 proportionally).

## 5. Selectivity lives in ECL2 and TM6/7 — not in the message

> "we propose that the **ECL2 and extracellular ends of TM6/7 serve as important
> selectivity filters** of ORs to distinct opioid peptides."

Supporting evidence they report: mOR K305(6.58)E **increased dynorphin potency 100-fold**;
replacing kOR's ECL2 with mOR's **decreased** dynorphin potency and efficacy; kOR's bulky
Y313(7.36) is predicted to clash with endomorphin and deltorphin; ECL3 side chains in dOR
clash with dynorphin's C-terminus. ECL2 composition differs sharply by receptor — kOR and
NOPR carry ~30% acidic residues in ECL2.

**Consequence, and it is a hard ceiling.** Taken with §4 — YGGF alone activates *all four*
receptors — a message-motif template **cannot be selective by construction**. Our screen
is a pan-opioid message matcher, not an mOR-specific one. Any claim of mOR selectivity
from this pipeline is unsupported, and the collapse of E3 under aromatic matching
(§9.2, 0.733 → 0.490) is what that limitation looks like in a number.

Reaching selectivity would require the template to span the ECL2-contacting region, i.e.
extent ≥ 8 for the E8–R213 contact and ≥ 9 for K9–D218. Extent 8 is now built. But the
*address* region proper (their §"C-terminal sequence variations address receptor
selectivity") sits at residues ~13–21, and a 21-residue reference cannot be covered by a
tractable candidate conformer ensemble — a 17-mer has 81 rotatable bonds and returns zero
conformers (§1). **This is a structural limit on the method, not a tuning problem.**

## 6. Endomorphin-1 binds differently from β-endorphin

Their Figure 2F–H: endomorphin's **W3 overlaps F4 of β-endorphin**, and endomorphin's
**F4 occupies an extended hydrophobic pocket formed by TM1, TM2 and TM7** that β-endorphin
does not use. Mutations in that extended pocket affect β-endorphin and endomorphin
*differently*.

**Consequence.** Our template is the β-endorphin binding mode. A peptide that binds the
endomorphin way is being scored against the wrong geometry — the same subpocket is
occupied by a residue one position earlier. This is a concrete argument for a
**two-template (consensus) screen**, taking the second template from 8F7R.

The cost of doing that: 8F7R is currently the held-out structure for E4, and using it for
a template forfeits the only held-out-structure validation in the study. That is a
trade-off to decide deliberately, not silently. Note that most dietary opioid peptides are
`YP…` (casomorphins, soymorphins, lactorphins) rather than `YGGF…`, and endomorphin is
`YPWF` — so the dietary library may be better matched to the **8F7R** mode than to our
8F7Q template.

## 7. N/OFQ shows Y1 is not universally required

N/OFQ begins with **phenylalanine**, not tyrosine, and is fully potent at NOPR. Their
reading: "Y1 is not necessary for the high potency of NFQ to NOPR", although Y1 is a
universal feature of the *classical* endogenous peptides and participates in the
water-mediated network of §3.

**Consequence.** Our free-N-terminal-Tyr flag (D6.11) is a reasonable prior for mOR but is
not a physical requirement of the receptor family. It is correctly kept out of the score
and carried as a separate column.

## 8. Actionable items this paper generates

| # | Item | Status |
|---|---|---|
| 1 | Rebuild the complex from the deposition; peptide extent an explicit choice | **done**, `22_prepare_complex.py` |
| 2 | Template at extent 8, matching the candidate window and capturing E8–R213(ECL2) | **done**, `results/template_p8` |
| 3 | Exclusion volumes to encode the chamber's steric bound (§4) | **done and validated**, 86 spheres / 32 residues incl. ECL2; see [`../DECISIONS.md`](../DECISIONS.md) D2.4 |
| 4 | Dietary-peptide screen against measured IC50 (§6 makes this the better-matched set) | **done** — [`dietary_screen.md`](dietary_screen.md). E8 agonist-vs-antagonist AUROC **0.917** (exact p=0.015) with exclusion volumes vs 0.444 without; confirms §4 |
| 5 | Record that PLIP evidence omits the D149 salt bridge; protonate before PLIP | open |
| 6 | Second template from 8F7R for the endomorphin mode; costs the E4 held-out test | **decided — proceed.** Table S2 settled it: endomorphin-1 is 42× tighter at μOR and the only μ-exclusive peptide, and 8F7R resolves its ligand completely (4/4 + amide cap) where 8F7Q is a truncation. The E4 cost is paid by re-pointing the held-out test to 8EFQ/DAMGO, not by dropping it. Design and pre-registration: [`dual_template_design.md`](dual_template_design.md); decisions D2.7–D2.9 |
| 7 | Water-mediated Y1 network cannot be represented without ordered waters | open; may be out of scope for Align-it |
| 8 | Stop any claim of mOR-selectivity; the message template is pan-opioid (§5) | **applied** to the write-up |
| 9 | Extract their Figure 1A / Table S2 pEC50 values across all four receptors as additional affinity data | **done** — `data/reference/wang2023_tableS2.csv`, 8 ligands × 4 receptors, tidy long form with pKi/pEC50/Emax and all SEMs, transcribed from the supplemental PDF and cross-checked against `Known_Opioid_Peptides_v11.xlsx` (0 value mismatches) |
