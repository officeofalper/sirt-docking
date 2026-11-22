"""Figures for version 2: search effort, ensemble docking, and which of the
two actually matters."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import docking as D
from ensemble import ENSEMBLES, LABELS, evaluate, load_v1

KIND_COLOUR = {"holo": D.COLOURS["holo"], "apo": D.COLOURS["apo"],
               "other": D.COLOURS["other"]}


def sweep_table():
    rows = json.load(open(os.path.join(D.FIG, "v2_sweep.json")))["rows"]
    by = {}
    for r in rows:
        by.setdefault((r["pdb"], r["exh"]), []).append(r)
    out = {}
    for (pid, exh), rs in by.items():
        tops = np.array([r["rmsds"][0] for r in rs])
        out[(pid, exh)] = {
            "kind": rs[0]["kind"],
            "rmsd_med": float(np.median(tops)),
            "rmsds": tops.tolist(),
            "score": float(np.mean([r["scores"][0] for r in rs])),
            "secs": float(np.mean([r["secs"] for r in rs])),
            "success": float((tops < D.SUCCESS_RMSD).mean()),
        }
    return out


def fig_sweep(tab, path):
    pids = sorted({p for p, _ in tab}, key=lambda p: tab[(p, 8)]["kind"])
    exhs = sorted({e for _, e in tab})
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.5))

    styles = {}
    for pid in pids:                      # the two apo share a colour: vary marker
        k = tab[(pid, exhs[0])]["kind"]
        styles[pid] = ["o-", "s--", "^:"][sum(1 for p in styles
                                              if tab[(p, exhs[0])]["kind"] == k)]
    for pid in pids:
        kind = tab[(pid, exhs[0])]["kind"]
        col = KIND_COLOUR[kind]
        med = [tab[(pid, e)]["rmsd_med"] for e in exhs]
        ax[0].plot(exhs, med, styles[pid], color=col, lw=2, ms=7, zorder=3,
                   label=f"{pid} ({kind})")
        for e in exhs:                                  # every individual run
            pts = tab[(pid, e)]["rmsds"]
            ax[0].scatter([e] * len(pts), pts, s=18, color=col, alpha=0.45,
                          edgecolors="none", zorder=2)
        ax[1].plot(exhs, [tab[(pid, e)]["secs"] for e in exhs], styles[pid],
                   color=col, lw=2, ms=7, label=f"{pid} ({kind})")

    ax[0].axhline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    ax[0].text(exhs[-1], D.SUCCESS_RMSD + 0.4, "2 $\\AA$ success threshold",
               ha="right", fontsize=9)
    ax[0].set_xscale("log", base=2)
    ax[0].set_xticks(exhs, [str(e) for e in exhs])
    ax[0].set_xlabel("Exhaustiveness (independent search runs)")
    ax[0].set_ylabel("RMSD of the top-scored pose ($\\AA$)")
    ax[0].set_title("A  Does searching harder help?", loc="left", fontweight="bold")
    ax[0].set_ylim(bottom=0)
    ax[0].legend(frameon=False, fontsize=9)

    ax[1].set_xscale("log", base=2)
    ax[1].set_xticks(exhs, [str(e) for e in exhs])
    ax[1].set_xlabel("Exhaustiveness")
    ax[1].set_ylabel("Wall-clock time per run (s)")
    ax[1].set_title("B  What it costs", loc="left", fontweight="bold")
    ax[1].legend(frameon=False, fontsize=9)

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# structures that contain SirReal2 itself -- the circularity to flag
SIRREAL_BOUND = {"4RMH", "4RMG"}


def fig_ensemble(path):
    top, seeds = load_v1()
    names = list(ENSEMBLES)
    res = {n: evaluate(top, seeds, ENSEMBLES[n]) for n in names}

    fig, ax = plt.subplots(1, 2, figsize=(13.0, 4.8))
    y = np.arange(len(names))
    vals = [res[n]["rmsd_med"] for n in names]
    cols = [D.COLOURS["holo"] if set(ENSEMBLES[n]) & SIRREAL_BOUND
            else D.COLOURS["apo"] for n in names]
    ax[0].barh(y, vals, color=cols, height=0.62)
    for i, n in enumerate(names):
        r = res[n]
        ax[0].text(max(r["rmsd_med"], D.SUCCESS_RMSD) + 0.3, i,
                   f"{r['rmsd_med']:.2f} $\\AA$  ({r['success']*100:.0f}% basarili)",
                   va="center", fontsize=9)
    ax[0].axvline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    ax[0].set_yticks(y, [LABELS[n] for n in names], fontsize=9)
    ax[0].invert_yaxis()
    ax[0].set_xlabel("RMSD of the pose the score selected ($\\AA$)")
    ax[0].set_xlim(0, 15)
    ax[0].set_title("A  Ensemble docking: only as good as its box",
                    loc="left", fontweight="bold")
    hand = [plt.Rectangle((0, 0), 1, 1, color=D.COLOURS["holo"]),
            plt.Rectangle((0, 0), 1, 1, color=D.COLOURS["apo"])]
    ax[0].legend(hand, ["SirReal2-bagli kristal iceriyor", "icermiyor"],
                 frameon=False, fontsize=8.5, loc="lower right")

    # does the score rank the conformations?
    for pid in {m for mem in ENSEMBLES.values() for m in mem}:
        sc = np.mean([top[(pid, s)][0] for s in seeds])
        rm = np.median([top[(pid, s)][1] for s in seeds])
        kind = "holo" if pid in SIRREAL_BOUND else ("other" if pid == "3ZGV" else "apo")
        ax[1].scatter(rm, sc, s=120, color=KIND_COLOUR[kind], zorder=3,
                      edgecolors="white", linewidths=1.5)
        ax[1].annotate(pid, (rm, sc), textcoords="offset points",
                       xytext=(8, -4), fontsize=9)
    ax[1].axvline(D.SUCCESS_RMSD, color=D.COLOURS["ink"], ls="--", lw=1.2)
    ax[1].set_xlabel("RMSD of the top pose ($\\AA$)")
    ax[1].set_ylabel("Mean Vina score (kcal/mol)")
    ax[1].set_title("B  The score ranks them, but by a thin margin",
                    loc="left", fontweight="bold")

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return res


def fig_pocket_fit(path):
    """Why it fails: the correct pose does not fit in the apo structures."""
    rows = json.load(open(os.path.join(D.FIG, "v2_pocketfit.json")))["rows"]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.5))
    x = np.arange(len(rows))
    cols = [KIND_COLOUR[r["kind"]] for r in rows]
    n_tot = rows[0]["n_lig_atoms"]

    ax[0].bar(x, [r["n_hard"] for r in rows], 0.6, color=cols)
    for i, r in enumerate(rows):
        ax[0].text(i, r["n_hard"] + 0.6, f"{r['n_hard']}/{n_tot}", ha="center",
                   fontsize=9)
        inside = r["n_hard"] > 4
        ax[0].text(i, r["n_hard"] / 2 if inside else 2.4,
                   f"{r['min_dist']:.2f} $\\AA$", ha="center", fontsize=8.5,
                   color="white" if inside else D.COLOURS["grey"])
    ax[0].set_xticks(x, [f"{r['pdb']}\n({r['kind']})" for r in rows])
    ax[0].set_ylabel(f"Ligand heavy atoms overlapping protein (of {n_tot})")
    ax[0].set_title("A  Does the correct pose even fit?", loc="left",
                    fontweight="bold")
    ax[0].set_ylim(0, n_tot + 3)
    ax[0].text(0.02, 0.95, "number inside each bar = closest atom-atom distance",
               transform=ax[0].transAxes, fontsize=8, color=D.COLOURS["grey"],
               va="top")

    sc = [r["score_only"] for r in rows]
    ax[1].bar(x, sc, 0.6, color=cols)
    for i, v in enumerate(sc):
        ax[1].text(i, v + (6 if v > 0 else -10), f"{v:+.1f}", ha="center", fontsize=9)
    ax[1].axhline(0, color=D.COLOURS["ink"], lw=1)
    ax[1].set_xticks(x, [f"{r['pdb']}\n({r['kind']})" for r in rows])
    ax[1].set_ylabel("Vina energy of the crystal pose (kcal/mol)")
    ax[1].set_title("B  Scored without any search", loc="left", fontweight="bold")
    ax[1].text(0.98, 0.5, "positive = severe clash\n(the pose is impossible here)",
               transform=ax[1].transAxes, fontsize=8.5, ha="right",
               color=D.COLOURS["grey"])

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return rows


def fig_effort_vs_receptor(tab, path):
    """Side by side: what 8x the compute buys, and what the right crystal buys."""
    apo = sorted({p for p, _ in tab if tab[(p, 8)]["kind"] == "apo"})
    exhs = sorted({e for _, e in tab})
    lo, hi = exhs[0], exhs[-1]

    gain_effort = [tab[(p, lo)]["rmsd_med"] - tab[(p, hi)]["rmsd_med"] for p in apo]
    cost_effort = [tab[(p, hi)]["secs"] / tab[(p, lo)]["secs"] for p in apo]
    gain_struct = [tab[(p, lo)]["rmsd_med"] - tab[("4RMH", lo)]["rmsd_med"] for p in apo]

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    x = np.arange(len(apo))
    w = 0.36
    ax.bar(x - w / 2, gain_effort, w, color=D.COLOURS["grey"],
           label=f"exhaustiveness {lo} $\\rightarrow$ {hi} "
                 f"({np.mean(cost_effort):.1f}$\\times$ hesaplama)")
    ax.bar(x + w / 2, gain_struct, w, color=D.COLOURS["holo"],
           label="dogru konformasyonu secmek (ayni hesaplama)")
    for i, (a, b) in enumerate(zip(gain_effort, gain_struct)):
        ax.text(i - w / 2, a + 0.2, f"{a:+.2f}", ha="center", fontsize=9)
        ax.text(i + w / 2, b + 0.2, f"{b:+.2f}", ha="center", fontsize=9)
    ax.axhline(0, color=D.COLOURS["ink"], lw=1)
    ax.set_xticks(x, [f"{p}\n(apo)" for p in apo])
    ax.set_ylabel("RMSD improvement ($\\AA$, higher = better)")
    ax.set_title("Search effort against receptor choice", loc="left",
                 fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(D.NOTE, fontsize=9, color=D.COLOURS["grey"], x=0.99, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return {"apo": apo, "gain_effort": gain_effort, "cost_effort": cost_effort,
            "gain_struct": gain_struct}


if __name__ == "__main__":
    os.makedirs(D.FIG, exist_ok=True)
    tab = sweep_table()
    exhs = sorted({e for _, e in tab})
    print(f"{'receptor':<8} {'exh':>4} {'RMSD med':>9} {'basari':>7} "
          f"{'score':>8} {'sure':>7}")
    print("-" * 52)
    for pid in sorted({p for p, _ in tab}):
        for e in exhs:
            r = tab[(pid, e)]
            print(f"{pid:<8} {e:>4} {r['rmsd_med']:>8.2f}A {r['success']*100:>6.0f}% "
                  f"{r['score']:>8.2f} {r['secs']:>6.1f}s")

    fig_sweep(tab, os.path.join(D.FIG, "fig3_exhaustiveness.png"))
    ens = fig_ensemble(os.path.join(D.FIG, "fig4_ensemble.png"))
    fit = fig_pocket_fit(os.path.join(D.FIG, "fig5_pocketfit.png"))
    eff = fig_effort_vs_receptor(tab, os.path.join(D.FIG, "fig6_effort.png"))
    json.dump({"sweep": {f"{p}_{e}": v for (p, e), v in tab.items()},
               "effort_vs_receptor": eff},
              open(os.path.join(D.FIG, "v2_summary.json"), "w"), indent=1)
    print("\nFigures: fig3_exhaustiveness.png fig4_ensemble.png "
          "fig5_pocketfit.png fig6_effort.png")
