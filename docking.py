"""
Shared docking machinery: read a PDB entry, split it into receptor and ligand,
prepare both for AutoDock Vina, run the docking, and measure how far the
predicted pose is from the crystal pose.

IMPORTANT: This repository runs a standard docking program on public crystal
structures. It is a methods exercise, not a drug-discovery campaign, and no
result here is a new scientific finding.

Pipeline
--------
    PDB entry ──┬─> receptor (protein + kept cofactors) ──> PDBQT
                └─> ligand (one HETATM residue)         ──> RDKit ──> PDBQT
                                                               │
                        AutoDock Vina (Trott & Olson, 2010) <──┘
                                   │
                        poses ──> RMSD against the crystal pose

The success criterion used throughout is the usual one in the docking
literature: a pose is "correct" when it lies within 2 A RMSD of the
experimental pose.
"""

import os
import subprocess
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
STRUCT = os.path.join(HERE, "structures")
WORK = os.path.join(HERE, "work")
FIG = os.path.join(HERE, "figures")
VINA = os.environ.get("VINA_BIN", "/home/claude/tools/vina/vina")

SUCCESS_RMSD = 2.0          # the usual threshold for a correct pose
SEED = 42

# Residues that are solvent or crystallisation additives, never the ligand.
JUNK = {"HOH", "DOD", "GOL", "EDO", "SO4", "PO4", "CL", "NA", "K", "MG", "CA",
        "ACT", "DMS", "MPD", "PEG", "TRS", "IMD", "FMT", "NO3", "ZN",
        "P6G", "PGE", "1PE", "2PE", "PG4", "EPE", "BME", "IOD", "BR"}
# Cofactors we keep as part of the receptor rather than treat as the ligand.
COFACTOR = {"NAD", "NAI", "APR", "ADP", "ATP", "CRC", "ZN"}


# ---------------------------------------------------------------------------
# 1. Reading the PDB entry
# ---------------------------------------------------------------------------
def read_pdb(pdb_id):
    path = os.path.join(STRUCT, f"{pdb_id.upper()}.pdb")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Download it from "
            f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb")
    with open(path) as f:
        return f.readlines()


def list_hetero(pdb_id):
    """Every non-solvent HETATM residue in the entry, with its atom count."""
    seen = {}
    for line in read_pdb(pdb_id):
        if not line.startswith("HETATM"):
            continue
        name = line[17:20].strip()
        if name in JUNK:
            continue
        key = (name, line[21], line[22:26].strip())
        seen[key] = seen.get(key, 0) + 1
    return sorted(seen.items(), key=lambda kv: -kv[1])


def split(pdb_id, ligand_resname, chain=None, keep_cofactors=True):
    """Write receptor.pdb and ligand.pdb for one entry.

    The ligand is the residue named `ligand_resname`; everything else that is
    protein (and, if asked, the cofactors) becomes the receptor. Waters and
    crystallisation additives are dropped. Only the first altloc is kept.
    """
    os.makedirs(WORK, exist_ok=True)
    rec, lig = [], []
    for line in read_pdb(pdb_id):
        tag = line[:6].strip()
        if tag not in ("ATOM", "HETATM"):
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        name, ch = line[17:20].strip(), line[21]
        if chain and ch != chain and tag == "ATOM":
            continue
        if name == ligand_resname and (chain is None or ch == chain):
            lig.append(line)
        elif tag == "ATOM":
            rec.append(line)
        elif name in COFACTOR and keep_cofactors:
            rec.append(line)
    if not lig:
        raise ValueError(f"no atoms for ligand {ligand_resname} in {pdb_id}")
    rp = os.path.join(WORK, f"{pdb_id}_receptor.pdb")
    lp = os.path.join(WORK, f"{pdb_id}_{ligand_resname}.pdb")
    open(rp, "w").write("".join(rec) + "END\n")
    open(lp, "w").write("".join(lig) + "END\n")
    return rp, lp


# ---------------------------------------------------------------------------
# 2. Preparation
# ---------------------------------------------------------------------------
def _run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def prepare_receptor(pdb_path):
    """Protonate and convert to PDBQT. -xr marks the receptor as rigid."""
    out = pdb_path.replace(".pdb", ".pdbqt")
    r = _run(["obabel", "-ipdb", pdb_path, "-opdbqt", "-O", out, "-xr", "-p", "7.4"])
    if not os.path.exists(out):
        raise RuntimeError(f"receptor prep failed: {r.stderr[-400:]}")
    return out


def ligand_from_pdb(lig_pdb, template_smiles=None):
    """Rebuild the crystal ligand as a proper RDKit molecule.

    A PDB file carries no bond orders, so without a reference SMILES the
    molecule would come back with every bond single. When the SMILES is given
    (from the PDB chemical component dictionary) it is used as a template.
    """
    mol = Chem.MolFromPDBFile(lig_pdb, removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"could not read {lig_pdb}")
    if template_smiles:
        ref = Chem.MolFromSmiles(template_smiles)
        mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)
    Chem.SanitizeMol(mol)
    return mol


def prepare_ligand(mol, out_path, embed=False):
    """Write a Vina-ready PDBQT. With embed=True the input geometry is thrown
    away and a fresh conformer generated, so the docking cannot start from the
    answer."""
    m = Chem.Mol(mol)
    if embed:
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=SEED)
        AllChem.MMFFOptimizeMolecule(m)
    else:
        m = Chem.AddHs(m, addCoords=True)
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    setups = MoleculePreparation().prepare(m)
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"ligand prep failed: {err}")
    open(out_path, "w").write(s)
    return out_path


# ---------------------------------------------------------------------------
# 3. Docking
# ---------------------------------------------------------------------------
def box_from_ligand(mol, padding=8.0):
    """A docking box centred on the crystal ligand. This is the standard
    'the binding site is known' setup; it does not tell the program the pose."""
    c = mol.GetConformer().GetPositions()
    centre = c.mean(axis=0)
    size = (c.max(axis=0) - c.min(axis=0)) + 2 * padding
    return centre, np.maximum(size, 16.0)


def dock(receptor_pdbqt, ligand_pdbqt, centre, size, out_pdbqt,
         exhaustiveness=16, n_modes=9, seed=SEED):
    cmd = [VINA, "--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt,
           "--center_x", f"{centre[0]:.3f}", "--center_y", f"{centre[1]:.3f}",
           "--center_z", f"{centre[2]:.3f}",
           "--size_x", f"{size[0]:.1f}", "--size_y", f"{size[1]:.1f}",
           "--size_z", f"{size[2]:.1f}",
           "--exhaustiveness", str(exhaustiveness), "--num_modes", str(n_modes),
           "--seed", str(seed), "--out", out_pdbqt]
    r = _run(cmd)
    if not os.path.exists(out_pdbqt):
        raise RuntimeError(f"vina failed: {r.stderr[-500:] or r.stdout[-500:]}")
    scores = [float(l.split()[1]) for l in r.stdout.splitlines()
              if l.strip() and l.split()[0].isdigit() and len(l.split()) >= 2]
    return scores


def poses_from_pdbqt(path):
    """Rebuild docked poses with meeko, which keeps the bond orders that a
    PDBQT round-trip would otherwise destroy."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(path, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    return Chem.RemoveHs(mols[0])


def rmsd_to_crystal(poses, crystal):
    """Symmetry-aware RMSD of every pose against the crystal pose, with no
    realignment: both are already in the same coordinate frame."""
    ref = Chem.RemoveHs(Chem.Mol(crystal))
    out = []
    for i in range(poses.GetNumConformers()):
        try:
            out.append(rdMolAlign.CalcRMS(poses, ref, prbId=i, refId=0))
        except Exception:
            out.append(float("nan"))
    return out


COLOURS = {"holo": "#0072B2", "apo": "#E69F00", "other": "#009E73",
           "ink": "#1d1422", "grey": "#8a8193", "warn": "#CC79A7"}
NOTE = "Standard docking on public structures: a methods exercise, not a new finding"


# ---------------------------------------------------------------------------
# 4. Superposition
# ---------------------------------------------------------------------------
def _ca_atoms(lines, chain):
    """{residue number: coordinates} for the C-alpha atoms of one chain."""
    out = {}
    for l in lines:
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == chain:
            if l[16] not in (" ", "A"):
                continue
            out[int(l[22:26])] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out


def kabsch(P, Q):
    """Rotation+translation that best maps P onto Q (both N x 3)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def superpose(pdb_id, chain, ref_id, ref_chain, out_path):
    """Move a whole structure onto a reference by its C-alpha atoms, matching
    residues by number. Needed so that a box defined on the reference ligand
    lands in the right place on every receptor, and so that RMSD against the
    reference crystal pose is meaningful.

    Returns the path written and the C-alpha RMSD after fitting.
    """
    mob_lines, ref_lines = read_pdb(pdb_id), read_pdb(ref_id)
    mob_ca, ref_ca = _ca_atoms(mob_lines, chain), _ca_atoms(ref_lines, ref_chain)
    common = sorted(set(mob_ca) & set(ref_ca))
    if len(common) < 30:
        raise ValueError(f"only {len(common)} residues in common between "
                         f"{pdb_id}:{chain} and {ref_id}:{ref_chain}")
    P = np.array([mob_ca[i] for i in common])
    Q = np.array([ref_ca[i] for i in common])
    R, t = kabsch(P, Q)
    rmsd = float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum(axis=1).mean()))

    out = []
    for l in mob_lines:
        if l[:6].strip() not in ("ATOM", "HETATM"):
            continue
        if l[16] not in (" ", "A"):
            continue
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        n = R @ xyz + t
        out.append(l[:30] + f"{n[0]:8.3f}{n[1]:8.3f}{n[2]:8.3f}" + l[54:])
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    open(out_path, "w").write("".join(out) + "END\n")
    return out_path, rmsd, len(common)


def pocket_residues(pdb_lines, ligand_xyz, cutoff=5.0, chain=None):
    """Residue numbers with any atom within `cutoff` of the ligand."""
    near = set()
    for l in pdb_lines:
        if not l.startswith("ATOM"):
            continue
        if chain and l[21] != chain:
            continue
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if np.min(np.linalg.norm(ligand_xyz - xyz, axis=1)) <= cutoff:
            near.add(int(l[22:26]))
    return sorted(near)


def residue_shift(lines_a, chain_a, lines_b, chain_b, resnums):
    """Per-residue C-alpha displacement between two superposed structures."""
    ca_a, ca_b = _ca_atoms(lines_a, chain_a), _ca_atoms(lines_b, chain_b)
    out = {}
    for r in resnums:
        if r in ca_a and r in ca_b:
            out[r] = float(np.linalg.norm(ca_a[r] - ca_b[r]))
    return out
