"""
Version 1: can docking find SirReal2's binding mode without being given the
rearranged protein?

IMPORTANT: This is a standard docking program (AutoDock Vina) run on public
crystal structures. It is a methods exercise, not a drug-discovery campaign,
and no number here is a new scientific finding.

Background
----------
SirReal2 is a selective SIRT2 inhibitor. Rumpf et al. (2015) showed that it
binds by rearranging the enzyme's active site and opening an adjacent pocket,
so the protein in the SirReal2 complex is NOT in its resting shape.

My review of AI-assisted drug discovery notes that docking usually treats the
receptor as rigid, and that success rates fall for proteins with substantial
induced-fit change. That claim is testable here, because both the rearranged
(holo) and the resting (apo) structures of SIRT2 are in the PDB.

The experiment
--------------
One ligand (SirReal2, taken from 4RMH and given a fresh conformer, so the
docking never starts from the answer) is docked into five receptors that differ
only in protein conformation: two SirReal2-bound, two unbound, one with the
reaction product. Every receptor is superposed onto 4RMH first, and cofactors
and the substrate peptide are stripped from all of them, so conformation is the
only variable. Each run is repeated with several random seeds.

A pose counts as correct at the usual threshold: within 2 A RMSD of the crystal
pose of SirReal2.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt

import docking as D
from prep_receptors import RECEPTORS, build

SIRREAL2_SMILES = "Cc1cc(C)nc(SCC(=O)Nc2ncc(Cc3cccc4ccccc34)s2)n1"
SEEDS = [1, 2, 3, 4, 5]
EXHAUSTIVENESS = 16
N_MODES = 9


def reference_ligand():
    _, lig_pdb = D.split("4RMH", "3TE", chain="A")
    return D.ligand_from_pdb(lig_pdb, template_smiles=SIRREAL2_SMILES)


def run():
    cry = reference_ligand()
    centre, size = D.box_from_ligand(cry, padding=8.0)
    lig = D.prepare_ligand(cry, os.path.join(D.WORK, "sirreal2.pdbqt"), embed=True)
    results = {}
    for pid, chain, kind, note in RECEPTORS:
        rec_pdb, fit_rmsd, _ = build(pid, chain)
        rec = D.prepare_receptor(rec_pdb)
        runs = []
        for seed in SEEDS:
            out = os.path.join(D.WORK, f"dock_{pid}_s{seed}.pdbqt")
            scores = D.dock(rec, lig, centre, size, out,
                            exhaustiveness=EXHAUSTIVENESS, n_modes=N_MODES, seed=seed)
            poses = D.poses_from_pdbqt(out)
            rmsds = D.rmsd_to_crystal(poses, cry)
            runs.append({"seed": seed, "scores": scores[:len(rmsds)], "rmsds": rmsds})
            print(f"  {pid} seed {seed}: top {scores[0]:6.2f} kcal/mol, "
                  f"RMSD {rmsds[0]:5.2f} A (best of {len(rmsds)}: {min(rmsds):5.2f} A)")
        results[pid] = {"kind": kind, "note": note, "fit_rmsd": fit_rmsd, "runs": runs}
    return results, cry


def summarise(results):
    rows = []
    for pid, r in results.items():
        top_rmsd = np.array([x["rmsds"][0] for x in r["runs"]])
        best_rmsd = np.array([min(x["rmsds"]) for x in r["runs"]])
        top_score = np.array([x["scores"][0] for x in r["runs"]])
        rows.append({
            "pdb": pid, "kind": r["kind"],
            "score": top_score.mean(), "score_sd": top_score.std(),
            "top_rmsd_med": float(np.median(top_rmsd)),
            "top_rmsd_min": float(top_rmsd.min()), "top_rmsd_max": float(top_rmsd.max()),
            "best_rmsd_med": float(np.median(best_rmsd)),
            "success": float((top_rmsd < D.SUCCESS_RMSD).mean()),
            "success_any": float((best_rmsd < D.SUCCESS_RMSD).mean()),
        })
    return rows


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
KIND_COLOUR = {"holo": D.COLOURS["holo"], "apo": D.COLOURS["apo"], "other": D.COLOURS["other"]}


def fig_pocket(path):
    """How far the pocket moves between each structure and the reference."""
    cry = reference_ligand()
    lig_xyz = cry.GetConformer().GetPositions()
    ref_lines = D.read_pdb("4RMH")
    pocket = D.pocket_residues(ref_lines, lig_xyz, 5.0, "A")
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    labels, allv, pocv = [], [], []
    for pid, chain, kind, _ in RECEPTORS[1:]:
        fit = open(os.path.join(D.WORK, f"{pid}_fit.pdb")).readlines()
        allsh = D.residue_shift(ref_lines, "A", fit, chain, sorted(D._ca_atoms(ref_lines, "A")))
        posh = D.residue_shift(ref_lines, "A", fit, chain, pocket)
        labels.append(f"{pid}\n({kind})")
        allv.append(np.mean(list(allsh.values())))
        pocv.append(np.mean(list(posh.values())))
        ax[1].plot(pocket, [posh.get(r, np.nan) for r in pocket], "o-", ms=4, lw=1.6,
                   color=KIND_COLOUR[kind], label=f"{pid} ({kind})")
    x = np.arange(len(labels))
    w = 0.36
    ax[0].bar(x - w / 2, allv, w, color=D.COLOURS["grey"], label="whole protein")
    ax[0].bar(x + w / 2, pocv, w, color=[KIND_COLOUR[r[2]] for r in RECEPTORS[1:]],
              label="SirReal2 pocket")
    for i, (a, b) in enumerate(zip(allv, pocv)):
        ax[0].text(i - w / 2, a + 0.05, f"{a:.2f}", ha="center", fontsize=9)
        ax[0].text(i + w / 2, b + 0.05, f"{b:.2f}", ha="center", fontsize=9)
    ax[0].set_xticks(x, labels)
    ax[0].set_ylabel("Mean C$\\alpha$ shift vs 4RMH ($\\AA$)")
    ax[0].set_title("A  How much does the protein move?", loc="left", fontweight="bold")
    ax[0].legend(frameon=False, fontsize=9)
    ax[1].set_xlabel("Residue number (pocket lining)")
    ax[1].set_ylabel("C$\\alpha$ shift vs 4RMH ($\\AA$)")
    ax[1].set_title("B  Residue by residue", loc="left", fontweight="bold")
    ax[1].legend(frameon=False, fontsize=9)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_results(results, path):
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.5))
    names = list(results)
    for i, pid in enumerate(names):
        r = results[pid]
        col = KIND_COLOUR[r["kind"]]
        tops = [x["rmsds"][0] for x in r["runs"]]
        ax[0].scatter([i] * len(tops), tops, s=55, color=col, zorder=3,
                      edgecolors="white", linewidths=1)
        ax[0].plot([i - 0.22, i + 0.22], [np.median(tops)] * 2, color=col, lw=2.5)
    ax[0].axhline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    ax[0].text(len(names) - 0.4, D.SUCCESS_RMSD + 0.25, "2 $\\AA$ success threshold",
               ha="right", fontsize=9)
    ax[0].set_xticks(range(len(names)),
                     [f"{p}\n({results[p]['kind']})" for p in names])
    ax[0].set_ylabel("RMSD of the top-scored pose ($\\AA$)")
    ax[0].set_title("A  Did docking find the right pose?", loc="left", fontweight="bold")
    ax[0].set_ylim(bottom=0)

    for pid in names:
        r = results[pid]
        col = KIND_COLOUR[r["kind"]]
        xs = [v for x in r["runs"] for v in x["rmsds"]]
        ys = [v for x in r["runs"] for v in x["scores"]]
        ax[1].scatter(xs, ys, s=22, alpha=0.6, color=col, edgecolors="none",
                      label=f"{pid} ({r['kind']})")
    ax[1].axvline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    ax[1].set_xlabel("RMSD to the crystal pose ($\\AA$)")
    ax[1].set_ylabel("Vina score (kcal/mol, lower = better)")
    ax[1].set_title("B  Every pose: score against correctness", loc="left", fontweight="bold")
    ax[1].legend(frameon=False, fontsize=8.5)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    print(f"Version 1: docking SirReal2 into {len(RECEPTORS)} SIRT2 conformations, "
          f"{len(SEEDS)} seeds each\n")
    results, cry = run()
    rows = summarise(results)
    print(f"\n{'receptor':<10} {'kind':<6} {'score':>14} {'top RMSD':>10} "
          f"{'range':>12} {'correct':>8}")
    for r in rows:
        print(f"{r['pdb']:<10} {r['kind']:<6} {r['score']:7.2f} +-{r['score_sd']:4.2f} "
              f"{r['top_rmsd_med']:10.2f} {r['top_rmsd_min']:5.1f}-{r['top_rmsd_max']:<5.1f} "
              f"{100*r['success']:7.0f}%")
    json.dump({"rows": rows,
               "results": {k: {kk: vv for kk, vv in v.items()} for k, v in results.items()}},
              open(os.path.join("figures", "v1_results.json"), "w"), indent=1)
    fig_pocket("figures/fig1_pocket.png")
    fig_results(results, "figures/fig2_docking.png")
    print("\nFigures saved in figures/")
