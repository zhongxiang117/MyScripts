
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFreeSASA

# help link: https://rdkit.org/docs/source/rdkit.Chem.rdFreeSASA.html


def compute_vdwsa(mol):
    radii = rdFreeSASA.classifyAtoms(mol)       # this will all be zeros
    asa = rdFreeSASA.CalcSASA(mol, radii)
    return asa


def compute_sasa(mol):
    # van der Waals radii (angstrom)
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
    sasa = rdFreeSASA.CalcSASA(mol, radii)
    return sasa


def compute_ligand_sasa_pocket(prot, lig, vdwsa=False):
    # compute complex SASA
    both = Chem.CombineMols(prot, lig)
    if vdwsa:
        sasa = compute_vdwsa(both)
    else:
        sasa = compute_sasa(both)
    # compute ligand SASA in pocket, ligand is the last component
    comp_lig = Chem.GetMolFrags(both, asMols=True, sanitizeFrags=False)[-1]
    lig_sasa_bound = sum([float(a.GetProp("SASA")) for a in comp_lig.GetAtoms()])
    return lig_sasa_bound


def compute_ligand_sasa_pocket_cutoff(protein, ligand, cutoff=8, vdwsa=False):
    xyz_protein = protein.GetConformer().GetPositions()
    xyz_ligand = ligand.GetConformer().GetPositions()
    # minimum distance between protein atoms and ligand atoms
    diff = xyz_protein[:, np.newaxis, :] - xyz_ligand[np.newaxis, :, :]
    r = np.min(np.linalg.norm(diff,axis=2), axis=1)
    indices = np.argwhere(r > cutoff).flatten()
    mol = Chem.RWMol(protein)
    for idx in sorted(indices, reverse=True):
        mol.RemoveAtom(int(idx))
    return compute_ligand_sasa_pocket(mol.GetMol(), ligand, vdwsa=vdwsa)



