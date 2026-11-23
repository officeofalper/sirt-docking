"""
Version 4: does docking accuracy track the chemical similarity of the ligand
that shaped the pocket?

Version 3 found that an open pocket is not enough -- it has to be open in
roughly the right shape -- and that 5MAT, whose pocket was opened by a
different chemotype, produced the best score in the study while being 9.7 A
wrong. That was one structure, so it was an anecdote.

Five further SIRT2 co-crystals span the similarity range from "SirReal2 plus a
bromine" to "nothing in common". If docking accuracy tracks that similarity,
the result becomes a prospective rule: pick the receptor whose bound ligand
most resembles your compound. Unlike the clash count used in version 2, ligand
similarity can be computed before docking, without knowing the answer.

Same ligand, same box, same settings and same seeds as versions 1 and 3.
"""
import json
import os

import numpy as np

import docking as D
from pocket_fit import CLOSE, HARD, receptor_heavy_atoms, score_only
from prep_receptors import build
from redock_crossdock import EXHAUSTIVENESS, N_MODES, SEEDS, reference_ligand

# pdb, chain, its own ligand, Tanimoto to SirReal2 (Morgan r=2, 1024 bits)
EXTRA = [
    ("5DY4", "A", "5GN", 0.814, "brominated SirReal analogue"),
    ("5DY5", "A", "5GR", 0.558, "SirReal affinity probe, triazole arm"),
    ("5YQL", "A", "A2I", 0.463, "shared head group, different tail"),
    ("9FDT", "A", "A1I", 0.194, "small pyrazole fragment"),
    ("5Y0Z", "A", "8K9", 0.120, "NPD11033, unrelated chemotype"),
]
OUT = os.path.join(D.FIG, "v4_results.json")


def main():
    cry = reference_ligand()
    centre, size = D.box_from_ligand(cry, padding=8.0)
    lig_xyz = np.array([cry.GetConformer().GetAtomPosition(i)
                        for i, a in enumerate(cry.GetAtoms()) if a.GetAtomicNum() > 1])
    cry_pdbqt = D.prepare_ligand(cry, os.path.join(D.WORK, "sirreal2_crystal.pdbqt"),
                                 embed=False)
    lig = D.prepare_ligand(cry, os.path.join(D.WORK, "sirreal2.pdbqt"), embed=True)

    fit, res = [], {}
    if os.path.exists(OUT):                       # resume if interrupted
        prev = json.load(open(OUT)); fit, res = prev["fit"], prev["results"]

    for pid, chain, own, tanimoto, note in EXTRA:
        if pid in res:
            continue
        rec_pdb, rms, _ = build(pid, chain)
        rec = D.prepare_receptor(rec_pdb)

        prot = receptor_heavy_atoms(rec_pdb)
        near = np.linalg.norm(lig_xyz[:, None, :] - prot[None, :, :], axis=2).min(axis=1)
        fit.append({"pdb": pid, "own_ligand": own, "tanimoto": tanimoto, "note": note,
                    "fit_rmsd": rms, "n_hard": int((near < HARD).sum()),
                    "n_close": int((near < CLOSE).sum()),
                    "min_dist": float(near.min()),
                    "score_only": score_only(rec, cry_pdbqt, centre, size),
                    "n_lig_atoms": int(len(lig_xyz)),
                    "nearest": [round(float(v), 3) for v in near]})
        f = fit[-1]
        print(f"  {pid} (T={tanimoto:.2f}) cakisan {f['n_hard']:2d}/{f['n_lig_atoms']}, "
              f"kristal poz {f['score_only']:+8.2f} kcal/mol", flush=True)

        runs = []
        for seed in SEEDS:
            path = os.path.join(D.WORK, f"dock_{pid}_s{seed}.pdbqt")
            scores = D.dock(rec, lig, centre, size, path,
                            exhaustiveness=EXHAUSTIVENESS, n_modes=N_MODES, seed=seed)
            poses = D.poses_from_pdbqt(path)
            rmsds = D.rmsd_to_crystal(poses, cry)
            runs.append({"seed": seed, "scores": scores[:len(rmsds)], "rmsds": rmsds})
            print(f"    seed {seed}: top {scores[0]:6.2f} kcal/mol, "
                  f"RMSD {rmsds[0]:6.2f} A (best {min(rmsds):6.2f} A)", flush=True)
        res[pid] = {"own_ligand": own, "tanimoto": tanimoto, "note": note, "runs": runs}

        top = np.array([r["rmsds"][0] for r in runs])
        sc = np.array([r["scores"][0] for r in runs])
        print(f"  => {pid}: skor {sc.mean():6.2f}, RMSD medyan {np.median(top):6.2f} A, "
              f"basari {100*(top < D.SUCCESS_RMSD).mean():3.0f}%\n", flush=True)
        json.dump({"fit": fit, "results": res}, open(OUT, "w"), indent=1)
    print("done")


if __name__ == "__main__":
    main()
