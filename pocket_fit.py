"""
Version 2, part C: was the search wrong, or was the question impossible?

Everything so far says docking "failed" on the apo structures. But there is
a prior question: in those structures, does the correct answer even exist?

The test puts SirReal2 in its *crystal* pose -- the known right answer --
into each superposed receptor and asks two things:

  1. Geometry: how many ligand heavy atoms sit closer to a protein heavy atom
     than two atoms can physically be (< 3.0 A for a non-bonded contact; < 2.5 A
     is a hard clash).
  2. Vina's own energy for that exact pose (--score_only, no search at all).

If the correct pose clashes badly and scores terribly in the apo structures,
then the search did not make a mistake. It was asked to place a ligand in a
pocket that is not open, and it correctly reported that it could not.
"""
import json
import os
import subprocess

import numpy as np
from rdkit import Chem

import docking as D
from prep_receptors import RECEPTORS, build
from redock_crossdock import reference_ligand

HARD = 2.5      # closer than this: atoms overlapping
CLOSE = 3.0     # closer than this: tighter than a normal non-bonded contact


def receptor_heavy_atoms(pdb_path):
    xyz = []
    for l in open(pdb_path):
        if l.startswith("ATOM") and l[76:78].strip() != "H":
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(xyz)


def score_only(receptor_pdbqt, ligand_pdbqt, centre, size):
    """Vina's energy for one fixed pose, with no conformational search.
    Vina 1.2 still wants a grid box even when it is not searching."""
    r = subprocess.run([D.VINA, "--receptor", receptor_pdbqt,
                        "--ligand", ligand_pdbqt, "--score_only",
                        "--center_x", f"{centre[0]:.3f}",
                        "--center_y", f"{centre[1]:.3f}",
                        "--center_z", f"{centre[2]:.3f}",
                        "--size_x", f"{size[0]:.1f}",
                        "--size_y", f"{size[1]:.1f}",
                        "--size_z", f"{size[2]:.1f}"],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith("Estimated Free Energy of Binding"):
            return float(line.split(":")[1].split()[0])
    raise RuntimeError(r.stdout[-400:] + r.stderr[-400:])


def main():
    cry = reference_ligand()
    centre, size = D.box_from_ligand(cry, padding=8.0)
    lig_xyz = np.array([cry.GetConformer().GetAtomPosition(i)
                        for i, a in enumerate(cry.GetAtoms())
                        if a.GetAtomicNum() > 1])
    # crystal pose, NOT a regenerated conformer
    cry_pdbqt = D.prepare_ligand(cry, os.path.join(D.WORK, "sirreal2_crystal.pdbqt"),
                                 embed=False)

    rows = []
    for pid, chain, kind, note in RECEPTORS:
        rec_pdb, _, _ = build(pid, chain)
        rec = D.prepare_receptor(rec_pdb)
        prot = receptor_heavy_atoms(rec_pdb)
        d = np.linalg.norm(lig_xyz[:, None, :] - prot[None, :, :], axis=2)
        nearest = d.min(axis=1)
        rows.append({
            "pdb": pid, "kind": kind,
            "n_hard": int((nearest < HARD).sum()),
            "n_close": int((nearest < CLOSE).sum()),
            "min_dist": float(nearest.min()),
            "worst_overlap": float(max(0.0, CLOSE - nearest.min())),
            "score_only": score_only(rec, cry_pdbqt, centre, size),
            "n_lig_atoms": int(len(lig_xyz)),
            "nearest": [round(float(v), 3) for v in nearest],
        })
        r = rows[-1]
        print(f"  {pid:5s} ({kind:5s})  cakisan atom: {r['n_hard']:2d}/{r['n_lig_atoms']} "
              f"(<{HARD} A), {r['n_close']:2d} (<{CLOSE} A)  "
              f"en kisa mesafe {r['min_dist']:4.2f} A  "
              f"kristal pozun skoru {r['score_only']:+7.2f} kcal/mol", flush=True)

    json.dump({"rows": rows, "hard": HARD, "close": CLOSE},
              open(os.path.join(D.FIG, "v2_pocketfit.json"), "w"), indent=1)
    return rows


if __name__ == "__main__":
    print("Kristal pozu (dogru cevap) her yapiya yerlestir ve sigip sigmadigina bak\n")
    main()
