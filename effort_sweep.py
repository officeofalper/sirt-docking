"""
Version 2, part A: is the apo failure a search problem?

The obvious first reaction to "docking put the ligand 11 A from where it
belongs" is: the search did not try hard enough. Vina's --exhaustiveness
controls how many independent Monte-Carlo runs it does, so raising it is
the standard answer. This script raises it 8-fold and measures what that
buys, on a receptor where docking works (4RMH, holo) and two where it
fails (1J8F and 3ZGO, both apo).

Writes results incrementally so a long run can be inspected while it goes.
"""
import json
import os
import time

import numpy as np

import docking as D
from prep_receptors import build
from redock_crossdock import SIRREAL2_SMILES, reference_ligand

EXH_LEVELS = [8, 16, 32, 64]
SEEDS = [1, 2, 3]
N_MODES = 9
RECEPTORS = [("4RMH", "A", "holo"), ("1J8F", "A", "apo"), ("3ZGO", "A", "apo")]
OUT = os.path.join(D.FIG, "v2_sweep.json")


def main():
    cry = reference_ligand()
    centre, size = D.box_from_ligand(cry, padding=8.0)
    lig = D.prepare_ligand(cry, os.path.join(D.WORK, "sirreal2.pdbqt"), embed=True)

    rows = []
    if os.path.exists(OUT):                     # resume if interrupted
        rows = json.load(open(OUT))["rows"]
    done = {(r["pdb"], r["exh"], r["seed"]) for r in rows}

    for pid, chain, kind in RECEPTORS:
        rec_pdb, _, _ = build(pid, chain)
        rec = D.prepare_receptor(rec_pdb)
        for exh in EXH_LEVELS:
            for seed in SEEDS:
                if (pid, exh, seed) in done:
                    continue
                out = os.path.join(D.WORK, f"sweep_{pid}_e{exh}_s{seed}.pdbqt")
                t0 = time.time()
                scores = D.dock(rec, lig, centre, size, out,
                                exhaustiveness=exh, n_modes=N_MODES, seed=seed)
                secs = time.time() - t0
                poses = D.poses_from_pdbqt(out)
                rmsds = D.rmsd_to_crystal(poses, cry)
                rows.append({"pdb": pid, "kind": kind, "exh": exh, "seed": seed,
                             "secs": secs,
                             "scores": scores[:len(rmsds)], "rmsds": rmsds})
                json.dump({"rows": rows}, open(OUT, "w"), indent=1)
                print(f"  {pid:5s} exh {exh:3d} seed {seed}: "
                      f"top {scores[0]:6.2f} kcal/mol, RMSD {rmsds[0]:6.2f} A "
                      f"(best {min(rmsds):6.2f} A)  {secs:5.1f}s", flush=True)
    print("done")


if __name__ == "__main__":
    main()
