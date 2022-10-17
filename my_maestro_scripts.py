"""My Scripts for Schrodinger Maestro Usage

Prefer syntax `from SCRIPTS import *`, then all customized functions
will be starts with `myfunc_*`
"""

from schrodinger.maestro import maestro

FEATURES = [
    'version 0.1.0  : My Schrodinger Scripts, Oct 17th, 2022',
    'version 0.2.0  : add `title`s functions',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()
__version__ = VERSION

__all__ = [
    'myfunc_get_selected_rows',
    'myfunc_get_selected_titles',
    'myfunc_get_selected_structures',
    'myfunc_get_selected_molecules',
    'myfunc_get_selected_atoms_separated_by_molecule',
    'myfunc_get_selected_atoms_separated_by_entry',
    'myfunc_get_selected_ids',
    'myfunc_get_selected_resnames_lists',
    'myfunc_get_selected_resnames_sets',
    'myfunc_get_selected_resnums',
    'myfunc_get_selected_atom_charges_separated_by_molecule',
    'myfunc_get_selected_atom_charges_separated_by_entry',
]


def myfunc_get_selected_rows():
    """1D: List[row, row, ...]"""
    pt = maestro.project_table_get()
    rows = [i for i in pt.selected_rows]
    print(f'Note: number of selected entries: {len(rows)}')
    return rows


def myfunc_get_selected_titles():
    """1D: List[str, str, ...]"""
    rows = myfunc_get_selected_rows()
    return [i.title for i in rows]


def myfunc_get_selected_structures():
    """1D: List[structure, structure, ...]"""
    rows = myfunc_get_selected_rows()
    return [i.structure for i in rows]


def myfunc_get_selected_molecules():
    """1D: List[mol, mol, ...]"""
    sts = myfunc_get_selected_structures()
    mols = [[m for m in s.molecule] for s in sts]
    print(f'Note: number of total selected molecules: {sum([len(i) for i in mols])}')
    return mols


def myfunc_get_selected_atoms_separated_by_molecule():
    """3D: List[[[[atom, atom, ...], [atom, atom, ...], ...], ...], ...]"""
    mols = myfunc_get_selected_molecules()
    atoms = [[[i for i in m.atom] for m in s] for s in mols]
    total = sum([sum([len(m) for m in s]) for s in atoms])
    print(f'Note: number of total selected atoms: {total}')
    return atoms


def myfunc_get_selected_atoms_separated_by_entry():
    """2D: List[[[atom, atom, ...], ...], ...]"""
    mols = myfunc_get_selected_molecules()
    atoms = [[i for m in s for i in m.atom] for s in mols]
    total = sum([len(s) for s in atoms])
    print(f'Note: number of total selected atoms: {total}')
    return atoms


def myfunc_get_selected_ids():
    """1D: List[str, str, ...]"""
    rows = myfunc_get_selected_rows()
    ids = [i.entry_id for i in rows]
    return ids


def myfunc_get_selected_resnames_lists():
    """1D: List[str, str, ...]"""
    sts = myfunc_get_selected_structures()
    resnames = [[k.pdbres for j in i.molecule for k in j.residue] for i in sts]
    print(f'Note: number of total selected molecules: {sum([len(i) for i in resnames])}')
    return resnames


def myfunc_get_selected_resnames_sets():
    """1D: List[str, str, ...]"""
    resnames = myfunc_get_selected_resnames_lists()
    return [list(set(i)) for i in resnames]


def myfunc_get_selected_resnums():
    """1D: List[int, int, ...]"""
    sts = myfunc_get_selected_structures()
    resnums = [[k.resnum for j in i.molecule for k in j.residue] for i in sts]
    print(f'Note: number of total selected molecules: {sum([len(i) for i in resnums])}')
    return resnums


def myfunc_get_selected_atom_charges_separated_by_molecule():
    """3D: List[[[[float, float, ...], [float, float, ...], ...], ...], ...]"""
    full = myfunc_get_selected_atoms_separated_by_molecule()
    charges = [[[a.partial_charge for a in m] for m in s] for s in full]
    return charges


def myfunc_get_selected_atom_charges_separated_by_entry():
    """2D: List[[[float, float, ...], ...], ...]"""
    full = myfunc_get_selected_atoms_separated_by_entry()
    charges = [[a.partial_charge for a in m] for m in full]
    return charges



