"""
Version 3: separating two variables that versions 1 and 2 confounded.

In the first five structures, "is a ligand bound?" and "is the selectivity
pocket open?" always moved together. Two further structures break that tie:

  5MAT  a different chemical series (7KJ, a thienopyrimidinone) bound, and
        it induces the same pocket          -> ligand YES, pocket OPEN
  9FDR  no ligand at all, but deposited as an apo form with the selectivity
        pocket open                          -> ligand NO,  pocket OPEN

With these the design is a near-complete 2x2, and the claim can move from
"apo structures fail" to the testable "closed-pocket structures fail,
whether or not a ligand is bound".

Same ligand, same box, same settings and same seeds as version 1.
"""
import json
import os

import numpy as np

import docking as D
from pocket_fit import CLOSE, HARD, receptor_heavy_atoms, score_only
from prep_receptors import build
from redock_crossdock import EXHAUSTIVENESS, N_MODES, SEEDS, reference_ligand

EXTRA = [
    ("5MAT", "A", "induced", "7KJ bound (thienopyrimidinone), pocket induced"),
    ("9FDR", "A", "apo-open", "no ligand bound, pocket reported open"),
]
OUT = os.path.join(D.FIG, "v3_results.json")


def fit_test(cry, centre, size):
    """Cheap and decisive: does the crystal pose physically fit?"""
    lig = np.array([cry.GetConformer().GetAtomPosition(i)
                    for i, a in enumerate(cry.GetAtoms()) if a.GetAtomicNum() > 1])
    cry_pdbqt = D.prepare_ligand(cry, os.path.join(D.WORK, "sirreal2_crystal.pdbqt"),
                                 embed=False)
    rows = []
    for pid, chain, kind, note in EXTRA:
        rec_pdb, rms, _ = build(pid, chain)
        rec = D.prepare_receptor(rec_pdb)
        prot = receptor_heavy_atoms(rec_pdb)
        near = np.linalg.norm(lig[:, None, :] - prot[None, :, :], axis=2).min(axis=1)
        rows.append({"pdb": pid, "kind": kind, "note": note, "fit_rmsd": rms,
                     "n_hard": int((near < HARD).sum()),
                     "n_close": int((near < CLOSE).sum()),
                     "min_dist": float(near.min()),
                     "score_only": score_only(rec, cry_pdbqt, centre, size),
                     "n_lig_atoms": int(len(lig)),
                     "nearest": [round(float(v), 3) for v in near]})
        r = rows[-1]
        print(f"  {pid} ({kind:9s}) hizalama {rms:.2f} A | cakisan {r['n_hard']:2d}/"
              f"{r['n_lig_atoms']} | en kisa {r['min_dist']:4.2f} A | "
              f"kristal poz {r['score_only']:+8.2f} kcal/mol", flush=True)
    return rows


def dock_new(cry, centre, size):
    lig = D.prepare_ligand(cry, os.path.join(D.WORK, "sirreal2.pdbqt"), embed=True)
    out = {}
    for pid, chain, kind, note in EXTRA:
        rec = D.prepare_receptor(build(pid, chain)[0])
        runs = []
        for seed in SEEDS:
            path = os.path.join(D.WORK, f"dock_{pid}_s{seed}.pdbqt")
            scores = D.dock(rec, lig, centre, size, path,
                            exhaustiveness=EXHAUSTIVENESS, n_modes=N_MODES, seed=seed)
            poses = D.poses_from_pdbqt(path)
            rmsds = D.rmsd_to_crystal(poses, cry)
            runs.append({"seed": seed, "scores": scores[:len(rmsds)], "rmsds": rmsds})
            print(f"  {pid} seed {seed}: top {scores[0]:6.2f} kcal/mol, "
                  f"RMSD {rmsds[0]:6.2f} A (best {min(rmsds):6.2f} A)", flush=True)
        out[pid] = {"kind": kind, "note": note, "runs": runs}
    return out


if __name__ == "__main__":
    cry = reference_ligand()
    centre, size = D.box_from_ligand(cry, padding=8.0)
    print("1) Kristal pozu yeni yapilara yerlestir -- cep acik mi?\n")
    fit = fit_test(cry, centre, size)
    print("\n2) SirReal2'yi bu yapilara dockla (5 tohum, exhaustiveness 16)\n")
    res = dock_new(cry, centre, size)
    for pid, r in res.items():
        top = np.array([x["rmsds"][0] for x in r["runs"]])
        sc = np.array([x["scores"][0] for x in r["runs"]])
        print(f"\n{pid} ({r['kind']}): skor {sc.mean():.2f} +-{sc.std():.2f}, "
              f"RMSD medyan {np.median(top):.2f} A, basari "
              f"{100*(top < D.SUCCESS_RMSD).mean():.0f}%")
    json.dump({"fit": fit, "results": res}, open(OUT, "w"), indent=1)
    print(f"\n{OUT}")
