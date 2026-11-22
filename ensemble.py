"""
Version 2, part B: does ensemble docking rescue the apo failure?

Ensemble docking is the standard remedy for induced fit: instead of one
receptor conformation, dock into several and keep the best-scoring result.
It needs no new docking runs here -- version 1 already docked the same
ligand into all five conformations with the same five seeds, so the
ensembles can be assembled from that output.

The question each ensemble answers:
  all5        the textbook best case  -- the bound conformation is in the box
  no_self     a holo crystal, but not the one this ligand came from
  apo_only    what you actually have for a new target
  apo_strict  only genuinely ligand-free structures

For each ensemble and seed: dock into every member, take the top pose of
each, and keep the one with the best score -- exactly what a user would do.
"Oracle" is the same choice made by someone who already knows the answer;
the gap between the two is how much the scoring function costs you.
"""
import itertools
import json
import os

import numpy as np

import docking as D

ENSEMBLES = {
    "4RMH only":   ["4RMH"],
    "1J8F only":   ["1J8F"],
    "apo_strict":  ["1J8F", "3ZGO"],
    "apo_only":    ["1J8F", "3ZGO", "3ZGV"],
    "no_self":     ["4RMG", "1J8F", "3ZGO", "3ZGV"],
    "all5":        ["4RMH", "4RMG", "1J8F", "3ZGO", "3ZGV"],
}
LABELS = {
    "4RMH only":  "Tek yapi: holo (kendi kristali)",
    "1J8F only":  "Tek yapi: apo",
    "apo_strict": "Ensemble: 2 apo",
    "apo_only":   "Ensemble: 2 apo + ADPR",
    "no_self":    "Ensemble: 4 yapi, baska holo dahil",
    "all5":       "Ensemble: 5 yapi (dogru konformasyon iceride)",
}


def load_v1():
    d = json.load(open(os.path.join(D.FIG, "v1_results.json")))["results"]
    # top pose per (receptor, seed): its score and its RMSD to the crystal
    top = {}
    for pid, r in d.items():
        for run in r["runs"]:
            top[(pid, run["seed"])] = (run["scores"][0], run["rmsds"][0])
    seeds = sorted({s for _, s in top})
    return top, seeds


def evaluate(top, seeds, members):
    """One ensemble: per seed, pick the member with the best score."""
    picked, oracle, chosen = [], [], []
    for s in seeds:
        cand = [(top[(m, s)][0], top[(m, s)][1], m) for m in members]
        sc, rm, who = min(cand, key=lambda x: x[0])        # best (most negative) score
        picked.append(rm)
        chosen.append(who)
        oracle.append(min(c[1] for c in cand))             # knows the answer already
    picked, oracle = np.array(picked), np.array(oracle)
    return {
        "members": members,
        "rmsd_med": float(np.median(picked)),
        "rmsd_min": float(picked.min()), "rmsd_max": float(picked.max()),
        "success": float((picked < D.SUCCESS_RMSD).mean()),
        "oracle_med": float(np.median(oracle)),
        "oracle_success": float((oracle < D.SUCCESS_RMSD).mean()),
        "chosen": chosen,
        "n_struct": len(members),
    }


def main():
    top, seeds = load_v1()
    out = {name: evaluate(top, seeds, mem) for name, mem in ENSEMBLES.items()}

    print(f"{'ensemble':<34} {'n':>2} {'RMSD med':>9} {'basari':>7} "
          f"{'oracle':>8} {'secilen yapi'}")
    print("-" * 92)
    for name, r in out.items():
        pick = ", ".join(sorted(set(r["chosen"])))
        print(f"{LABELS[name]:<34} {r['n_struct']:>2} {r['rmsd_med']:>8.2f}A "
              f"{r['success']*100:>6.0f}% {r['oracle_med']:>7.2f}A  {pick}")

    # Does the score rank the conformations correctly?
    rank = sorted({m for mem in ENSEMBLES.values() for m in mem},
                  key=lambda m: np.mean([top[(m, s)][0] for s in seeds]))
    print("\nSkora gore yapi sirasi (en iyiden en kotuye):")
    for i, m in enumerate(rank, 1):
        sc = np.mean([top[(m, s)][0] for s in seeds])
        rm = np.median([top[(m, s)][1] for s in seeds])
        print(f"  {i}. {m}  {sc:7.2f} kcal/mol   poz RMSD {rm:6.2f} A")

    json.dump({"ensembles": out, "score_rank": rank},
              open(os.path.join(D.FIG, "v2_ensemble.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
