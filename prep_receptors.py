"""Build one comparable receptor per structure.

Design choice: the receptor is the PROTEIN CHAIN ONLY. Cofactors (NAD+, ADP-ribose)
and the substrate peptide are removed everywhere, so the single thing that differs
between receptors is the conformation of the protein itself. Every structure is
first superposed onto 4RMH so that one docking box serves them all.
"""
import os
import docking as D

REF_ID, REF_CHAIN, LIGAND = "4RMH", "A", "3TE"
RECEPTORS = [
    ("4RMH", "A", "holo", "SirReal2-bound (reference, 1.42 A)"),
    ("4RMG", "A", "holo", "SirReal2-bound, second crystal form (1.88 A)"),
    ("1J8F", "A", "apo", "no ligand bound"),
    ("3ZGO", "A", "apo", "no ligand bound, independent crystal"),
    ("3ZGV", "A", "other", "ADP-ribose product bound"),
]


def protein_only(src_lines, chain):
    return [l for l in src_lines if l.startswith("ATOM") and l[21] == chain
            and l[16] in (" ", "A")]


def build(pdb_id, chain):
    """Superpose onto the reference, then keep protein chain atoms only."""
    os.makedirs(D.WORK, exist_ok=True)
    if pdb_id == REF_ID:
        lines, rms, n = D.read_pdb(pdb_id), 0.0, None
    else:
        fit, rms, n = D.superpose(pdb_id, chain, REF_ID, REF_CHAIN,
                                  os.path.join(D.WORK, f"{pdb_id}_fit.pdb"))
        lines = open(fit).readlines()
    out = os.path.join(D.WORK, f"rec_{pdb_id}.pdb")
    open(out, "w").write("".join(protein_only(lines, chain)) + "END\n")
    return out, rms, n


if __name__ == "__main__":
    for pid, ch, kind, _ in RECEPTORS:
        p, rms, n = build(pid, ch)
        q = D.prepare_receptor(p)
        print(f"{pid} ({kind}): fitted Ca RMSD {rms:.2f} A -> {os.path.basename(q)}")
