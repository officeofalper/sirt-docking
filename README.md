# Was the Docking Wrong, or Was the Question Impossible? Induced Fit in SIRT2

> **Disclaimer.** This repository runs **standard, published docking software on public crystal structures**. Nothing here is a new scientific finding, and none of it was done in a laboratory. The SIRT2 induced-fit effect it measures was described by the crystallographers who solved the structures (Rumpf et al., 2015). The purpose is to turn one limitation named in my review — that rigid-receptor docking cannot represent conformational change — into something I measured myself, with a number attached.

**Author:** Alper Öztürk

**Related paper:**

Öztürk, A. (2026). AI-Assisted Drug Discovery: Predicting Protein–Ligand Interactions and Identifying Candidate Therapeutic Molecules. *Next Frontier For Life Sciences and AI*, 10(1). https://doi.org/10.62802/nflsai.110

**Companion projects:** [binding-affinity-leakage](https://github.com/officeofalper/binding-affinity-leakage) (does the model learn binding or memorise the benchmark?), [cancer-flux-balance](https://github.com/officeofalper/cancer-flux-balance) and [tme-glucose-competition](https://github.com/officeofalper/tme-glucose-competition). That project asked whether a reported accuracy is real. This one asks something narrower and more physical: whether the input a docking program is given can support a correct answer at all.

| | Version 1: `redock_crossdock.py` | Version 2: `effort_sweep.py`, `ensemble.py`, `pocket_fit.py` | Version 3: `v3_newstructures.py` |
|---|---|---|---|
| Question | Does docking find the right pose when the receptor conformation is wrong? | Can anything fix it — and why does it fail? | Is the problem the missing ligand, or the closed pocket? |
| Design | One ligand, five SIRT2 conformations, five random seeds | Search effort 8→64; six ensembles; the crystal pose scored in every receptor | Two further structures that break the tie between those two variables |
| Output | Figures 1–2 | Figures 3–6 | Figure 7 |

---

## 1. Connection to my research

| What my review says | How this repository tests it |
|---|---|
| Docking treats the receptor as rigid, so it cannot represent induced fit; this is listed as a central limitation of structure-based virtual screening (Section 3). | SirReal2 is docked into the structure it came from and into four other conformations of the same protein. Re-docking reaches 0.35 Å; the apo structures reach 8–11 Å and never succeed. |
| Scoring functions are the weak link: they rank poses imperfectly and correlate only loosely with measured affinity (Section 3.2). | The score falls by only 3.5 kcal/mol between the structure where docking works and the ones where it fails — not enough to flag the failure. The score-versus-RMSD funnel that exists for the correct receptor disappears for the wrong one. |
| Ensemble docking and flexible-receptor methods are proposed as remedies (Section 3.4). | Six ensembles are compared. Ensemble docking can only interpolate among the conformations it is given, and version 3 shows that adding a structure can make it *worse*: a receptor that scores better can be the one that is more wrong. |
| Hybrid physics-plus-machine-learning approaches are presented as the direction of travel (Section 6). | Not tested here. This repository is about what the physics-based half is and is not doing. |

**What is mine and what is not.** The structures, the software, and the induced-fit phenomenon are other people's work, cited in Section 7. The experimental design, every number in Section 3, the pocket-fit measurement in Section 3.4 and the limitations in Section 4 are my own.

---

## 2. Method

### The ligand

SirReal2, PDB chemical component **3TE**, 29 heavy atoms. It is a selective SIRT2 inhibitor that works by opening a pocket that does not exist in the unbound enzyme — the reason it is the right test case here.

```
Cc1cc(C)nc(SCC(=O)Nc2ncc(Cc3cccc4ccccc34)s2)n1
```

The crystal coordinates carry no bond orders, so the ligand is rebuilt with `AssignBondOrdersFromTemplate` against that SMILES string. For docking, its 3D conformer is **generated from scratch** (`embed=True`), so the program is not handed the answer it is supposed to find.

### The receptors

Every structure is superposed onto 4RMH by Cα Kabsch alignment, then reduced to **the protein chain alone**: zinc, NAD⁺, ADP-ribose, the substrate peptide and all solvent are removed. That is deliberate — it leaves the conformation of the protein backbone and side chains as the only thing that differs between receptors.

| PDB | Resolution | What is bound in the crystal | Role here |
|---|---|---|---|
| **4RMH** | 1.42 Å | SirReal2 + acetyl-lysine peptide | the reference; re-docking control |
| **4RMG** | 1.88 Å | SirReal2 + NAD⁺ | same ligand, different crystal form |
| **1J8F** | 1.70 Å | nothing | apo |
| **3ZGO** | 1.63 Å | nothing | apo, independent crystal |
| **3ZGV** | 2.27 Å | ADP-ribose (reaction product) | a third conformational state |
| **5MAT** | 2.07 Å | 7KJ, a thienopyrimidinone | ligand bound, pocket open — but opened by a different chemotype |
| **9FDR** | 1.25 Å | nothing | no ligand, pocket open — deposited as an apo form with the selectivity pocket opened |

The last two were added in version 3; versions 1 and 2 used the first five. One docking box serves all five: centred on the crystal ligand, padded by 8 Å. This tells the program **where** the site is, which is the normal assumption in re-docking studies; it does not tell it the pose.

### Docking

AutoDock Vina 1.2.5. Receptors protonated at pH 7.4 and converted with Open Babel (`-xr`, rigid); ligands prepared with Meeko. Version 1: exhaustiveness 16, 9 modes, seeds 1–5. Version 2 sweep: exhaustiveness 8, 16, 32, 64, seeds 1–3.

RMSD to the crystal pose is **symmetry-aware** (`rdMolAlign.CalcRMS`), computed with no realignment because both molecules are already in the same frame. A pose under **2 Å** counts as correct, the usual convention.

### Does the pipeline work at all?

Two independent checks, both run before any result below was trusted:

1. Vina's own bundled example (imatinib into 1IEP) re-docks to **0.36 Å**.
2. This repository's own module, re-docking SirReal2 into 4RMH from a regenerated conformer, reaches **0.72 Å** on the first attempt.

---

## 3. Results

### 3.1 Version 1 — re-docking and cross-docking

Five seeds each, exhaustiveness 16. "Correct" means the **top-scored** pose is under 2 Å.

| Receptor | State | Score (kcal/mol) | Top-pose RMSD (median) | Correct |
|---|---|---|---|---|
| 4RMH | holo, self | −12.19 ± 0.02 | **0.35 Å** | **100%** |
| 4RMG | holo, other crystal | −10.91 ± 0.47 | 0.90 Å | 60% |
| 1J8F | apo | −8.71 ± 0.05 | **11.21 Å** | **0%** |
| 3ZGO | apo | −8.52 ± 0.28 | 8.38 Å | 0% |
| 3ZGV | ADP-ribose | −10.12 ± 0.38 | 9.68 Å | 0% |

Two things matter more than the failure itself.

**The score does not warn you.** From the structure where docking is essentially perfect to the one where it is 11 Å wrong, the score falls by 3.5 kcal/mol. −8.7 kcal/mol is a result most people would keep. In Figure 2B the holo structures show a funnel — better score, better pose — and for the apo structures that funnel is simply absent: every pose sits at 8–12 Å with a respectable score.

**Even the same conformational state is not a guarantee.** 4RMG is the same protein with the same ligand bound, from a different crystal form, and it still fails one run in five.

### 3.2 Where the protein moves

Pocket residues are defined geometrically: every residue with an atom within 5 Å of the crystal ligand. That gives 22 residues, among them Ile93, Ala135, Leu138, Tyr139, Pro140, Phe143, Phe190 and Leu206 — which is, without being told, the selectivity-pocket lining that Rumpf et al. (2015) describe.

| Receptor | Whole protein | Pocket | Largest movers |
|---|---|---|---|
| 4RMG (holo) | 0.72 Å | 0.71 Å | 96, 104 |
| 1J8F (apo) | 1.53 Å | **2.12 Å** | **140, 96, 139, 138** |
| 3ZGO (apo) | 1.51 Å | **2.12 Å** | **140, 96, 139, 138** |
| 3ZGV (ADPR) | 2.37 Å | 2.20 Å | 104, 206, 139 |

In the apo structures the pocket moves **more than the protein as a whole** — the change is concentrated at the binding site, which is what induced fit means. Between the two holo structures it is not. Two independently solved apo crystals give the same mean shift to two decimal places and the same four largest movers, so this is a property of the unbound state, not of one crystal.

### 3.3 Version 2a — searching harder does not help

Exhaustiveness is the number of independent Monte-Carlo runs Vina performs. Raising it is the standard first response to a bad docking result.

| Receptor | exh 8 | exh 16 | exh 32 | exh 64 | Success |
|---|---|---|---|---|---|
| 4RMH (holo) | 0.35 Å | 0.35 Å | 0.35 Å | 0.37 Å | 100% throughout |
| 1J8F (apo) | 11.13 Å | 11.13 Å | 11.22 Å | 11.20 Å | 0% throughout |
| 3ZGO (apo) | 8.45 Å | 8.45 Å | 8.38 Å | 8.38 Å | 0% throughout |

Eight times the search, **7.5× the compute**, and the apo results move by 0.07 Å — within run-to-run noise. Meanwhile the holo structure is already perfect at exhaustiveness 8, half the default.

Put the two levers side by side (Figure 6):

| On the apo structures | RMSD improvement | Cost |
|---|---|---|
| exhaustiveness 8 → 64 | −0.07 Å (1J8F), +0.07 Å (3ZGO) | 7.5× compute |
| choosing the right conformation | **+10.78 Å**, **+8.10 Å** | none |

### 3.4 Version 2b — why: the correct answer does not fit

Before calling this a docking failure, there is a prior question: in those structures, is the correct pose physically available? The test places SirReal2 in its **crystal** pose into each superposed receptor and asks two things — how badly it overlaps the protein, and what Vina's energy function says about that exact pose with no search at all (`--score_only`).

| Receptor | Heavy atoms overlapping protein (<2.5 Å) | Closest atom–atom distance | Energy of the crystal pose |
|---|---|---|---|
| 4RMH (holo) | 0 / 29 | 3.27 Å | **−11.62 kcal/mol** |
| 4RMG (holo) | 0 / 29 | 2.92 Å | −9.30 |
| 3ZGV (ADPR) | 12 / 29 | 0.61 Å | **+122.84** |
| 1J8F (apo) | 26 / 29 | 0.65 Å | **+188.37** |
| 3ZGO (apo) | 27 / 29 | 0.63 Å | **+186.71** |

In the apo structures, 26 to 27 of the ligand's 29 heavy atoms sit inside protein atoms, with closest approaches around 0.63 Å — atoms essentially superimposed. The correct pose is not merely disfavoured there; it is impossible.

This changes what the whole repository is reporting:

> **The search did not make a mistake.** It was asked to place a ligand in a pocket that is not open, it correctly found that it could not, and it returned the best available alternative — 11 Å away, at −8.7 kcal/mol — with no indication that anything was wrong.

3ZGV is the informative middle case: 12 clashing atoms rather than 26, consistent with the product-bound state ordering the site partway.

### 3.5 Version 2c — ensemble docking

Ensemble docking is the usual remedy: dock into several conformations and keep the best-scoring result. No new runs are needed — version 1 already docked the same ligand into all five with the same seeds.

| Ensemble | Structures | RMSD of the pose the score chose | Correct | Which it chose |
|---|---|---|---|---|
| Single: holo (self) | 1 | 0.35 Å | 100% | 4RMH |
| Single: apo | 1 | 11.21 Å | 0% | 1J8F |
| 2 apo | 2 | 11.21 Å | 0% | 1J8F |
| 2 apo + ADPR | 3 | 9.68 Å | 0% | 3ZGV |
| 4 structures, other holo included | 4 | 0.90 Å | 60% | 4RMG |
| All 5 | 5 | 0.35 Å | 100% | 4RMH |

Adding structures does nothing by itself. What the score does well is **choose**: in every ensemble containing an open-pocket structure it picked that structure, and the pose it selected matched the pose an oracle would have selected. But the margin is thin — the worst success (4RMG, −10.91) and the best failure (3ZGV, −10.12) are 0.8 kcal/mol apart, and among two equally wrong apo structures the score picked the worse one.

**The honest reading, and the flaw in this experiment.** 4RMG contains SirReal2 too; it differs from 4RMH in having NAD⁺ rather than the acetyl-lysine peptide. So the only thing that rescued ensemble docking here was adding a crystal of this same complex — which is circular, because in a real prospective project a crystal of your ligand bound to your target is the thing you do not have. What this experiment demonstrates is narrower than "ensemble docking works": it is that ensemble docking **can only interpolate among the conformations it is given**, and none of the ligand-free conformations available here has the pocket open.

Version 3 removes the circularity by adding two structures that do not contain SirReal2.

### 3.6 Version 3 — pocket shape, not ligand presence

In the first five structures, *is a ligand bound?* and *is the pocket open?* always moved together. **5MAT** (a different chemical series bound, pocket induced) and **9FDR** (no ligand, pocket open) break that tie. Same ligand, same box, same settings, same seeds as version 1.

| Receptor | Pocket as deposited | Ligand in crystal | Crystal-pose clashes | Crystal-pose energy | Score | RMSD | Correct |
|---|---|---|---|---|---|---|---|
| 4RMH | open, by SirReal2 | yes | 0 / 29 | −11.6 | −12.19 | **0.35 Å** | 100% |
| 4RMG | open, by SirReal2 | yes | 0 / 29 | −9.3 | −10.91 | 0.90 Å | 60% |
| **9FDR** | **open, ligand-free** | **no** | 7 / 29 | +45.0 | −10.50 | **3.48 Å** | 0% |
| **5MAT** | **open, by 7KJ** | **yes** | 9 / 29 | +56.7 | **−12.45** | **9.67 Å** | 0% |
| 3ZGV | partly open | yes | 12 / 29 | +122.8 | −10.12 | 9.68 Å | 0% |
| 1J8F | closed | no | 26 / 29 | +188.4 | −8.71 | 11.21 Å | 0% |
| 3ZGO | closed | no | 27 / 29 | +186.7 | −8.52 | 8.38 Å | 0% |

**Three findings, in order of importance.**

**1. The prediction was wrong, and that is the useful part.** I expected both new structures to succeed. Neither reaches 2 Å. An open pocket is **not sufficient** — it has to be open in roughly the right shape. 5MAT's pocket was shaped by 7KJ, and 9FDR's by yet another molecule during soaking, so neither accommodates the exact SirReal2 pose (7–9 clashing atoms, against 0 for the SirReal2-bound crystals).

**2. But the 2×2 question is answered.** Ordering the seven structures by how open the pocket is reproduces the docking outcome; ordering them by whether a ligand is bound does not. Across the seven, Spearman ρ = **+0.77** (p = 0.04) for clash count against RMSD, and **−0.29** (p = 0.53) for ligand presence. With n = 7 these are a direction, not a proof — but the single clearest comparison needs no statistics: **ligand-free 9FDR reaches 3.48 Å while ligand-bound 5MAT reaches 9.67 Å**, and 9FDR beats the closed apo structures by nearly 8 Å. So the claim moves from the wrong *apo structures fail* to the supported *closed-pocket structures fail, whether or not a ligand is bound.*

**3. The most alarming result in the whole repository.** 5MAT produces **the best score of any receptor tested** — −12.45 kcal/mol, better than the correct receptor's −12.19 — while being 9.67 Å wrong. Section 3.1 reported that the score falls by only 3.5 kcal/mol when the receptor is wrong, which is not enough to warn you. This is worse: with a plausibly-open but wrongly-shaped pocket, **the score goes up**. The most dangerous receptor is not the obviously closed one; it is the one that looks right.

**The failure mode also changes.** In a closed pocket the correct pose is not in the search space at all. In 5MAT it is reachable — one of its nine modes lands at 3.33 Å — but the scoring function ranks a wrong pose above it. Closed pocket is a search-space problem; open-but-wrong-shape is a scoring problem.

**An alternative explanation, ruled out.** 5MAT superposes onto 4RMH less well than the others (1.61 Å Cα RMSD, against 0.53 Å for 9FDR), so its pocket might simply have fallen outside the shared docking box. It has not: after superposition 5MAT's own ligand 7KJ sits **1.29 Å** from the centre of the SirReal2 crystal pose, with all 36 of its heavy atoms inside the box and a median nearest-atom distance of 1.82 Å to the SirReal2 pose. The box is over the right site; the failure is real.

---

## 4. What this does not show

- **One ligand, one target.** Every number comes from SirReal2 and SIRT2. SIRT2 was chosen *because* its induced fit is dramatic, so the size of the effect here is not an estimate of the average case — it is close to a worst case. Nothing here says how often rigid-receptor docking fails in general.
- **Seven structures is a small n.** The correlations in Section 3.6 are reported with their p-values and should be read as a direction, not as evidence. The comparison that carries the argument is the single contrast between 9FDR and 5MAT, not the rank statistic.
- **"Pocket open" is measured by the answer.** Clash count is computed against the SirReal2 crystal pose, so it measures *openness to this ligand*, not openness in the abstract. That is the right quantity for this question but it is circular as a general descriptor: you cannot compute it prospectively, because it needs the answer. A prospective version would need a pocket-volume or shape measure that does not reference the ligand.
- **5MAT and 9FDR may have been opened by their own ligands.** 9FDR is ligand-free as deposited, but the crystals were obtained by adapting the protein to a different compound. Neither structure is a clean sample of "what the unbound enzyme does on its own".
- **The site is given.** The box is centred on the known ligand. Blind docking, where the site must also be found, is a harder and separate problem.
- **Cofactors are removed from every receptor.** This isolates protein conformation, but it is not how one would dock for a real prediction: SIRT2's pocket forms at an interface where NAD⁺ sits, and a realistic setup would keep it.
- **One program, one scoring function.** A different program might behave differently. The claim that a bad pose scores well enough to pass is a claim about Vina's function, not about scoring functions in general — though the clash energies in Section 3.4 are a physical fact that any function would have to reproduce.
- **`--score_only` numbers are not binding free energies.** +188 kcal/mol is the clash term dominating; read it as "severe overlap", not as a thermodynamic quantity.
- **Crystal structures are models fitted to density,** not photographs. Side-chain positions in a 2.27 Å structure carry real uncertainty. The agreement between two independent apo crystals is the main reason to trust the apo result.

## 5. What I would do next

Version 3 answered the question it was built for and raised a sharper one. Three directions, in the order I would take them:

1. **Is the 5MAT effect general?** It produced the best score in the study while being 9.7 Å wrong. One structure is an anecdote. The test: assemble every SIRT2 co-crystal with a non-SirReal ligand in the selectivity pocket — 5MAR, 5Y5N, 5Y0Z, 5D7P and 5D7Q are candidates reported in the literature — and ask how often an open-but-differently-shaped pocket produces a confidently-scored wrong answer. If it is common, it is a practical warning about receptor choice in virtual screening rather than a curiosity.
2. **A prospective measure of openness.** Section 4 notes that clash count needs the answer in order to be computed. Replacing it with a ligand-free descriptor — pocket volume, or the separation between the residues lining the two walls — would let the same ranking be produced *before* docking, which is when it would be useful.
3. **Does the right receptor rescue enrichment, not just pose?** This repository measures pose accuracy for one known binder. What a screening project cares about is whether actives rank above decoys. That needs a curated active/decoy set, and therefore ChEMBL.

Friedrich et al. (2025) state directly that the selectivity pocket is absent in the previously available SIRT2 apo structures — independent confirmation of Section 3.4, reached from crystallography rather than from docking.

---

## 6. Running it

```bash
pip install -r requirements.txt     # rdkit, meeko, numpy, scipy, matplotlib
# AutoDock Vina 1.2.5 binary, and Open Babel on the PATH
# Put 4RMH, 4RMG, 1J8F, 3ZGO, 3ZGV (and for version 3, 5MAT and 9FDR) in structures/

python prep_receptors.py     # superpose, strip to protein, convert
python redock_crossdock.py   # version 1  (~15 min, 25 docking runs)
python effort_sweep.py       # version 2a (~45 min, 36 runs; resumable)
python pocket_fit.py         # version 2b (seconds)
python ensemble.py           # version 2c (seconds, reuses version 1)
python v2_figures.py         # figures 3-6
python v3_newstructures.py   # version 3  (~6 min; needs 5MAT and 9FDR too)
python v3_figures.py         # figure 7, and the seven-structure table
```

`VINA_BIN` overrides the path to the Vina executable. `effort_sweep.py` saves after every run and skips what it has already done, so it can be interrupted.

Timings are from 2 CPU cores.

---

## 7. References

Structures and the phenomenon:

- Rumpf, T., Schiedel, M., Karaman, B., et al. (2015). Selective Sirt2 inhibition by ligand-induced rearrangement of the active site. *Nature Communications*, 6, 6263. https://doi.org/10.1038/ncomms7263 — 4RMG, 4RMH; the selectivity pocket.
- Finnin, M.S., Donigian, J.R., & Pavletich, N.P. (2001). Structure of the histone deacetylase SIRT2. *Nature Structural Biology*, 8, 621–625. https://doi.org/10.1038/89668 — 1J8F.
- Moniot, S., Schutkowski, M., & Steegborn, C. (2013). Crystal structure analysis of human Sirt2 and its ADP-ribose complex. *Journal of Structural Biology*, 182, 136. https://doi.org/10.1016/j.jsb.2013.02.012 — 3ZGO, 3ZGV.
- Sundriyal, S., Moniot, S., Mahmud, Z., et al. (2017). Thienopyrimidinone based sirtuin-2 (SIRT2)-selective inhibitors bind in the ligand induced selectivity pocket. *Journal of Medicinal Chemistry*, 60, 1928–1945. https://doi.org/10.1021/acs.jmedchem.6b01690 — 5MAT.
- Friedrich, F., Schiedel, M., Swyter, S., et al. (2025). Efficient crystallization of apo Sirt2 for small-molecule soaking and structural analysis of ligand interactions. *Journal of Medicinal Chemistry*, 68, 10771–10780. https://doi.org/10.1021/acs.jmedchem.4c02896 — 9FDR.

Software:

- Eberhardt, J., Santos-Martins, D., Tillack, A.F., & Forli, S. (2021). AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. *Journal of Chemical Information and Modeling*, 61, 3891–3898. https://doi.org/10.1021/acs.jcim.1c00203
- Trott, O., & Olson, A.J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry*, 31, 455–461. https://doi.org/10.1002/jcc.21334
- O'Boyle, N.M., Banck, M., James, C.A., et al. (2011). Open Babel: an open chemical toolbox. *Journal of Cheminformatics*, 3, 33. https://doi.org/10.1186/1758-2946-3-33
- RDKit: open-source cheminformatics. https://www.rdkit.org
- Meeko, Forli Lab. https://github.com/forlilab/Meeko
