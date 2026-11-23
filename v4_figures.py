"""Version 4 analysis: ligand similarity against docking outcome, over eight
SIRT2 structures."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import docking as D

# structures whose pocket was shaped by a ligand, plus their Tanimoto to SirReal2
LIGANDED = {
    "4RMH": (1.000, "SirReal2 (self)"),
    "4RMG": (1.000, "SirReal2, other crystal"),
    "5DY4": (0.814, "brominated SirReal"),
    "5DY5": (0.558, "SirReal probe"),
    "5YQL": (0.463, "shared head group"),
    "5MAT": (0.235, "thienopyrimidinone"),
    "9FDT": (0.194, "pyrazole fragment"),
    "5Y0Z": (0.120, "NPD11033"),
}


def collect():
    v1 = json.load(open(os.path.join(D.FIG, "v1_results.json")))["results"]
    v3 = json.load(open(os.path.join(D.FIG, "v3_results.json")))
    v4 = json.load(open(os.path.join(D.FIG, "v4_results.json")))
    runs = dict(v1) | v3["results"] | v4["results"]
    fit = ({r["pdb"]: r for r in json.load(open(os.path.join(D.FIG, "v2_pocketfit.json")))["rows"]}
           | {r["pdb"]: r for r in v3["fit"]} | {r["pdb"]: r for r in v4["fit"]})
    rows = []
    for pid, (tan, label) in LIGANDED.items():
        top = np.array([x["rmsds"][0] for x in runs[pid]["runs"]])
        sc = np.array([x["scores"][0] for x in runs[pid]["runs"]])
        rows.append({"pdb": pid, "tanimoto": tan, "label": label,
                     "rmsd": float(np.median(top)), "score": float(sc.mean()),
                     "score_sd": float(sc.std()), "clash": fit[pid]["n_hard"],
                     "energy": fit[pid]["score_only"],
                     "success": float((top < D.SUCCESS_RMSD).mean())})
    return sorted(rows, key=lambda r: -r["tanimoto"])


def figure(rows, path):
    ok = [r for r in rows if r["success"] > 0.5]
    bad = [r for r in rows if r["success"] <= 0.5]
    cg, cb = D.COLOURS["holo"], D.COLOURS["warn"]
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # A: the headline
    for grp, c in ((ok, cg), (bad, cb)):
        ax[0].scatter([r["tanimoto"] for r in grp], [r["rmsd"] for r in grp],
                      s=150, color=c, zorder=3, edgecolors="white", linewidths=1.5)
    # explicit offsets: 4RMH/4RMG share x, and the three failures sit close
    offA = {"4RMH": (10, -9), "4RMG": (10, 2), "5DY4": (10, -4), "5DY5": (10, -4),
            "5YQL": (-10, -4), "5MAT": (10, -4), "9FDT": (6, -14), "5Y0Z": (2, 9)}
    for r in rows:
        dx, dy = offA[r["pdb"]]
        ax[0].annotate(r["pdb"], (r["tanimoto"], r["rmsd"]), fontsize=9,
                       textcoords="offset points", xytext=(dx, dy),
                       ha="right" if dx < 0 else "left")
    lo, hi = max(r["tanimoto"] for r in bad), min(r["tanimoto"] for r in ok)
    ax[0].axvspan(lo, hi, color=D.COLOURS["grey"], alpha=.14, zorder=1)
    ax[0].text((lo + hi) / 2, 6.2, "no structure\nin this range", ha="center",
               fontsize=8.5, color=D.COLOURS["grey"])
    ax[0].axhline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    ax[0].text(1.0, D.SUCCESS_RMSD + 0.4, "2 $\\AA$", ha="right", fontsize=9)
    ax[0].set_xlabel("Tanimoto similarity of the pocket-shaping ligand to SirReal2")
    ax[0].set_ylabel("Docking RMSD, top pose ($\\AA$)")
    ax[0].set_title("A  A cliff, not a slope", loc="left", fontweight="bold")
    ax[0].set_ylim(-0.5, 11)
    ax[0].set_xlim(0.04, 1.10)

    # B: the mechanism
    for grp, c in ((ok, cg), (bad, cb)):
        ax[1].scatter([r["tanimoto"] for r in grp], [r["clash"] for r in grp],
                      s=150, color=c, zorder=3, edgecolors="white", linewidths=1.5)
    offB = {"4RMH": (10, -9), "4RMG": (10, 2), "5DY4": (10, -4), "5DY5": (10, -4),
            "5YQL": (-10, -4), "5MAT": (10, -4), "9FDT": (10, -4), "5Y0Z": (-10, -4)}
    for r in rows:
        dx, dy = offB[r["pdb"]]
        ax[1].annotate(r["pdb"], (r["tanimoto"], r["clash"]), fontsize=9,
                       textcoords="offset points", xytext=(dx, dy),
                       ha="right" if dx < 0 else "left")
    ax[1].set_xlabel("Tanimoto similarity to SirReal2")
    ax[1].set_ylabel("Crystal-pose heavy atoms overlapping protein (of 29)")
    ax[1].set_title("B  Why: a similar ligand leaves a usable pocket",
                    loc="left", fontweight="bold")
    ax[1].set_ylim(-0.6, 11)

    # C: the warning
    for i, (grp, c, name) in enumerate(((ok, cg, "docking correct"),
                                        (bad, cb, "docking wrong"))):
        v = [r["score"] for r in grp]
        ax[2].scatter([i] * len(v), v, s=150, color=c, zorder=3,
                      edgecolors="white", linewidths=1.5)
        offC = {"4RMH": (12, -4), "4RMG": (12, -4), "5DY4": (12, 2),
                "5DY5": (12, -4), "5YQL": (12, -9), "5MAT": (12, -4),
                "9FDT": (12, -4), "5Y0Z": (12, -4)}
        for r in grp:
            dx, dy = offC[r["pdb"]]
            ax[2].annotate(r["pdb"], (i, r["score"]), fontsize=9,
                           textcoords="offset points", xytext=(dx, dy))
        ax[2].hlines(np.mean(v), i - .22, i + .22, color=c, lw=2.5)
    worst = min(bad, key=lambda r: r["score"])
    ax[2].annotate("best score in the whole\nstudy, 9.7 $\\AA$ wrong",
                   (1, worst["score"]), fontsize=8.5, ha="right", va="center",
                   textcoords="offset points", xytext=(-26, 16),
                   color=cb, arrowprops=dict(arrowstyle="->", color=cb, lw=1.4))
    ax[2].set_xticks([0, 1], ["docking correct", "docking wrong"])
    ax[2].set_xlim(-0.5, 1.6)
    ax[2].set_ylabel("Vina score (kcal/mol)")
    ax[2].set_title("C  The score cannot tell them apart", loc="left",
                    fontweight="bold")

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.995, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    rows = collect()
    print(f"{'PDB':<6} {'ligand':<26} {'Tanimoto':>9} {'cakisan':>8} "
          f"{'kristal poz':>12} {'skor':>8} {'RMSD':>8} {'basari':>7}")
    print("-" * 92)
    for r in rows:
        print(f"{r['pdb']:<6} {r['label']:<26} {r['tanimoto']:>9.3f} {r['clash']:>5}/29 "
              f"{r['energy']:>+11.1f} {r['score']:>8.2f} {r['rmsd']:>7.2f}A "
              f"{100*r['success']:>6.0f}%")

    ok = [r for r in rows if r["success"] > 0.5]
    bad = [r for r in rows if r["success"] <= 0.5]
    print(f"\nBasarili  : T >= {min(r['tanimoto'] for r in ok):.3f}  "
          f"(n={len(ok)}), cakisan atom {max(r['clash'] for r in ok)} veya daha az")
    print(f"Basarisiz : T <= {max(r['tanimoto'] for r in bad):.3f}  "
          f"(n={len(bad)}), cakisan atom {min(r['clash'] for r in bad)} veya daha cok")
    so, sb = [r["score"] for r in ok], [r["score"] for r in bad]
    print(f"\nSkor araligi, dogru pozlar   : {min(so):.2f} .. {max(so):.2f} kcal/mol")
    print(f"Skor araligi, yanlis pozlar  : {min(sb):.2f} .. {max(sb):.2f} kcal/mol")
    print(f"  -> ortusme tam; en iyi skor YANLIS bir pozdan geliyor "
          f"({min(bad, key=lambda r: r['score'])['pdb']}, {min(sb):.2f})")
    t = [r["tanimoto"] for r in rows]
    print(f"\nSpearman T-RMSD : rho = {stats.spearmanr(t, [r['rmsd'] for r in rows]).statistic:+.2f}"
          f"   T-cakisma : rho = {stats.spearmanr(t, [r['clash'] for r in rows]).statistic:+.2f}"
          f"   (n = {len(rows)})")

    figure(rows, os.path.join(D.FIG, "fig8_similarity.png"))
    json.dump({"rows": rows}, open(os.path.join(D.FIG, "v4_summary.json"), "w"), indent=1)
    print("\nfigures/fig8_similarity.png")
