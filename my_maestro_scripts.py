"""My Scripts for Schrodinger Maestro Usage

Prefer to using syntax `from SCRIPTS import *`, then all customized functions
will be starting with `myfunc_*`

Important: All returned entries are sorted by row_number values
"""

from schrodinger.maestro import maestro

FEATURES = [
    'version 0.1.0  : My Schrodinger Scripts, Oct 17th, 2022',
    'version 0.2.0  : add `title`s functions',
    'version 0.3.0  : add selections for workspace',
    'version 0.4.0  : add included functions for workspace',
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
    'myfunc_ws_get_chosen_atoms_ids',
    'myfunc_ws_get_selected_atoms_ids_detail',
    'myfunc_ws_get_selected_atoms',
    'myfunc_ws_get_included_atoms_ids_detail',
    'myfunc_ws_get_included_atoms',
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
    """3D: List[[[atom, atom, ...], [atom, atom, ...], ...], ...]"""
    mols = myfunc_get_selected_molecules()
    atoms = [[[i for i in m.atom] for m in s] for s in mols]
    total = sum([sum([len(m) for m in s]) for s in atoms])
    print(f'Note: number of total selected atoms: {total}')
    return atoms


def myfunc_get_selected_atoms_separated_by_entry():
    """2D: List[[atom, atom, ...], ...]"""
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
    """3D: List[[[float, float, ...], [float, float, ...], ...], ...]"""
    full = myfunc_get_selected_atoms_separated_by_molecule()
    charges = [[[a.partial_charge for a in m] for m in s] for s in full]
    return charges


def myfunc_get_selected_atom_charges_separated_by_entry():
    """2D: List[[float, float, ...], ...]"""
    full = myfunc_get_selected_atoms_separated_by_entry()
    charges = [[a.partial_charge for a in m] for m in full]
    return charges


def myfunc_ws_get_chosen_atoms_ids():
    """1D: List[int, int, ...]"""
    aids = maestro.selected_atoms_get()
    print(f'Note: workspace: number of chosen atoms: {len(aids)}')
    return [i-1 for i in aids]


def _calc_ids(rows,mids):
    # separate to sublist
    totlist = [0]
    for i,e in enumerate(rows):
        totlist.append(totlist[i]+sum([len(m) for m in e]))
    smids = sorted(mids)
    slist = []
    b = 0
    for e in range(len(totlist)-1):
        ls = []
        while b < len(smids):
            v = smids[b]
            if v >= totlist[e+1]:
                break
            b += 1
            ls.append(v-totlist[e])
        slist.append(ls)

    ids = []
    for eid,mols,subs in zip(range(len(rows)),rows,slist):
        srst = sorted(subs)
        srst.append(-1)     # to avoid overflow
        i = 0
        offset = 0
        for mid,atoms in enumerate(mols):
            n = len(atoms)
            for j in range(n):
                if j == srst[i]-offset:
                    ids.append((eid,mid,j))
                    i += 1
            offset += n
    return ids


def myfunc_ws_get_selected_atoms_ids_detail():
    """1D: List[(eid,mid,aid), (eid,mid,aid), ...]"""
    rows = myfunc_get_selected_atoms_separated_by_molecule()
    mids = myfunc_ws_get_chosen_atoms_ids()
    return _calc_ids(rows,mids)


def myfunc_ws_get_selected_atoms():
    """2D: List[[atom, atom, ...], ...]"""
    rows = myfunc_get_selected_atoms_separated_by_molecule()
    mids = myfunc_ws_get_chosen_atoms_ids()
    ids = _calc_ids(rows,mids)
    return [rows[i[0]][i[1]][i[2]] for i in ids]


def myfunc_ws_get_included_atoms_ids_detail():
    """1D: List[(eid,mid,aid), (eid,mid,aid), ...]
    
    Caution:
        By testing, the sequence of the included structures matters,
        to eliminate this influence, the structures are sorted by their
        row_number values, to keep the consistent with the selected functions
    """
    pt = maestro.project_table_get()
    frows = [r for r in pt.included_rows]
    fsts = [i.getStructure() for i in frows]
    num = [i.row_number for i in frows]
    fnum = sorted(range(len(num)),key=lambda t: num[t])
    sts = [fsts[i] for i in fnum]
    rows = [[[a for a in m.atom] for m in r.molecule] for r in sts]
    print(f'Note: number of included entries: {len(sts)}')
    print(f'Note: number of total included molecules: {sum([len(m) for m in rows])}')
    total = sum([sum([len(m) for m in s]) for s in rows])
    print(f'Note: number of total included atoms: {total}')
    mids = myfunc_ws_get_chosen_atoms_ids()
    return _calc_ids(rows,mids)


def myfunc_ws_get_included_atoms():
    """2D: List[[atom, atom, ...], ...]"""
    pt = maestro.project_table_get()
    frows = [r for r in pt.included_rows]
    fsts = [i.getStructure() for i in frows]
    num = [i.row_number for i in frows]
    fnum = sorted(range(len(num)),key=lambda t: num[t])
    sts = [fsts[i] for i in fnum]
    rows = [[[a for a in m.atom] for m in r.molecule] for r in sts]
    print(f'Note: number of included entries: {len(sts)}')
    print(f'Note: number of total included molecules: {sum([len(m) for m in rows])}')
    total = sum([sum([len(m) for m in s]) for s in rows])
    print(f'Note: number of total included atoms: {total}')
    mids = myfunc_ws_get_chosen_atoms_ids()
    ids = _calc_ids(rows,mids)
    return [rows[i[0]][i[1]][i[2]] for i in ids]



