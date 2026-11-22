"""Version 3 figure: what actually predicts whether docking works."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import docking as D

# ligand_bound, pocket state as deposited, short label
INFO = {
    "4RMH": (True,  "SirReal2 ile acik", "holo"),
    "4RMG": (True,  "SirReal2 ile acik", "holo"),
    "5MAT": (True,  "7KJ ile acik",      "induced"),
    "9FDR": (False, "ligandsiz acik",    "apo-open"),
    "3ZGV": (True,  "kismen",            "other"),
    "1J8F": (False, "kapali",            "apo"),
    "3ZGO": (False, "kapali",            "apo"),
}
COL = {"holo": D.COLOURS["holo"], "apo": D.COLOURS["apo"],
       "other": D.COLOURS["other"], "induced": D.COLOURS["warn"],
       "apo-open": "#6a2c91"}


def collect():
    v1 = json.load(open(os.path.join(D.FIG, "v1_results.json")))["results"]
    fit1 = json.load(open(os.path.join(D.FIG, "v2_pocketfit.json")))["rows"]
    v3 = json.load(open(os.path.join(D.FIG, "v3_results.json")))
    fit = {r["pdb"]: r for r in fit1} | {r["pdb"]: r for r in v3["fit"]}
    runs = dict(v1) | v3["results"]
    rows = []
    for pid, f in fit.items():
        top = np.array([x["rmsds"][0] for x in runs[pid]["runs"]])
        best = np.array([min(x["rmsds"]) for x in runs[pid]["runs"]])
        sc = np.array([x["scores"][0] for x in runs[pid]["runs"]])
        bound, pocket, kind = INFO[pid]
        rows.append({"pdb": pid, "kind": kind, "bound": bound, "pocket": pocket,
                     "clash": f["n_hard"], "energy": f["score_only"],
                     "rmsd": float(np.median(top)),
                     "best_rmsd": float(np.median(best)),
                     "score": float(sc.mean()), "score_sd": float(sc.std()),
                     "success": float((top < D.SUCCESS_RMSD).mean())})
    return sorted(rows, key=lambda r: r["clash"])


def figure(rows, path):
    fig, ax = plt.subplots(1, 2, figsize=(13.0, 4.8))

    for r in rows:
        c = COL[r["kind"]]
        m, fc = ("o", c) if r["bound"] else ("s", "none")
        ax[0].scatter(r["clash"], r["rmsd"], s=150, marker=m, facecolors=fc,
                      edgecolors=c, linewidths=2.2, zorder=3)
        off = (-11, -4) if r["pdb"] == "5MAT" else (9, -4)
        ax[0].annotate(r["pdb"], (r["clash"], r["rmsd"]), fontsize=9.5,
                       textcoords="offset points", xytext=off,
                       ha="right" if off[0] < 0 else "left")
    ax[0].axhline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    ax[0].text(27, D.SUCCESS_RMSD + 0.4, "2 $\\AA$ success threshold", ha="right",
               fontsize=9)
    ax[0].set_xlabel("Crystal pose: ligand heavy atoms overlapping protein (of 29)")
    ax[0].set_ylabel("Docking RMSD, top pose ($\\AA$)")
    ax[0].set_title("A  Pocket shape, not ligand presence, tracks the outcome",
                    loc="left", fontweight="bold")
    ax[0].set_ylim(-0.5, 13)
    h = [plt.Line2D([], [], marker="o", ls="", mfc=D.COLOURS["grey"],
                    mec=D.COLOURS["grey"], ms=11, label="kristalde ligand var"),
         plt.Line2D([], [], marker="s", ls="", mfc="none",
                    mec=D.COLOURS["grey"], mew=2.2, ms=11, label="ligandsiz")]
    ax[0].legend(handles=h, frameon=False, fontsize=9, loc="lower right")

    for r in rows:
        c = COL[r["kind"]]
        m, fc = ("o", c) if r["bound"] else ("s", "none")
        ax[1].scatter(r["rmsd"], r["score"], s=150, marker=m, facecolors=fc,
                      edgecolors=c, linewidths=2.2, zorder=3)
        ax[1].annotate(r["pdb"], (r["rmsd"], r["score"]), fontsize=9.5,
                       textcoords="offset points", xytext=(9, -4))
    ax[1].axvline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    # the dangerous point: best score in the whole study, and badly wrong
    worst = min((r for r in rows if r["rmsd"] > D.SUCCESS_RMSD),
                key=lambda r: r["score"])
    ax[1].annotate(f"calismanin en iyi skoru,\n{worst['rmsd']:.1f} $\\AA$ yanlis",
                   (worst["rmsd"], worst["score"]), fontsize=9,
                   textcoords="offset points", xytext=(-22, 34), ha="right",
                   color=D.COLOURS["warn"],
                   arrowprops=dict(arrowstyle="->", color=D.COLOURS["warn"], lw=1.4))
    ax[1].set_xlabel("Docking RMSD, top pose ($\\AA$)")
    ax[1].set_ylabel("Vina score (kcal/mol)")
    ax[1].set_title("B  The score cannot tell you which receptor to trust",
                    loc="left", fontweight="bold")

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    rows = collect()
    print(f"{'PDB':<6} {'cep':<18} {'ligand':<8} {'cakisan':>8} {'kristal poz':>12} "
          f"{'skor':>8} {'RMSD':>8} {'en iyi poz':>11} {'basari':>7}")
    print("-" * 94)
    for r in rows:
        print(f"{r['pdb']:<6} {r['pocket']:<18} {'var' if r['bound'] else 'yok':<8} "
              f"{r['clash']:>5}/29 {r['energy']:>+11.1f} {r['score']:>8.2f} "
              f"{r['rmsd']:>7.2f}A {r['best_rmsd']:>10.2f}A {100*r['success']:>6.0f}%")

    cl = [r["clash"] for r in rows]; rm = [r["rmsd"] for r in rows]
    bd = [int(r["bound"]) for r in rows]
    rho_c = stats.spearmanr(cl, rm)
    rho_b = stats.spearmanr(bd, rm)
    print(f"\nSpearman, cakisma sayisi - RMSD : rho = {rho_c.statistic:+.2f} "
          f"(p = {rho_c.pvalue:.3f})")
    print(f"Spearman, ligand varligi - RMSD : rho = {rho_b.statistic:+.2f} "
          f"(p = {rho_b.pvalue:.3f})")
    print("n = 7; bu korelasyonlar yon gostergesi, kanit degil.")

    figure(rows, os.path.join(D.FIG, "fig7_pocket_vs_outcome.png"))
    json.dump({"rows": rows,
               "spearman_clash": [rho_c.statistic, rho_c.pvalue],
               "spearman_bound": [rho_b.statistic, rho_b.pvalue]},
              open(os.path.join(D.FIG, "v3_summary.json"), "w"), indent=1)
    print("\nfigures/fig7_pocket_vs_outcome.png")
