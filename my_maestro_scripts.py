"""My Scripts for Schrodinger Maestro Usage

Prefer to using syntax `from SCRIPTS import *`, then all customized functions
will be starting with `myfunc_*`

Care:
    1) atom's `index' starts from 1
        -> `st.molecule' starts from 1, `st.atom' starts from 1
    2) `[a for a in st.atom] != [a for m in st.molecule for a in m]' (important!!)
    3) by test, `maestro.project_table_get().included_rows' returns reversed
       order of entries, and the operation order matters. It is different
       with `selected_rows', which is sorted by row number.
"""

from schrodinger.maestro import maestro

FEATURES = [
    'version 0.1.0  : My Schrodinger Scripts, Oct 17th, 2022',
    'version 0.2.0  : add `title`s functions',
    'version 0.3.0  : add selections for workspace',
    'version 0.4.0  : add included functions for workspace',
    'version 0.5.0  : make sure properties are copied when getting structures',
    'version 0.6.0  : add more funcs and make their names more clear',
    'version 0.7.0  : add correspondent `included` functions',
    'version 0.8.0  : add `molecule` methods; more workspace specific',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()
__version__ = VERSION

__all__ = [
    'myfunc_get_selected_rows',
    'myfunc_get_included_rows',
    'myfunc_get_selected_titles',
    'myfunc_get_included_titles',
    'myfunc_get_selected_structures',
    'myfunc_get_included_structures',
    'myfunc_get_selected_molecules',
    'myfunc_get_included_molecules',
    'myfunc_get_selected_atoms_separated_by_molecule',
    'myfunc_get_included_atoms_separated_by_molecule',
    'myfunc_get_selected_atoms_separated_by_entry',
    'myfunc_get_included_atoms_separated_by_entry',
    'myfunc_get_selected_entry_ids',
    'myfunc_get_included_entry_ids',
    'myfunc_get_selected_resnames_lists',
    'myfunc_get_included_resnames_lists',
    'myfunc_get_selected_resnames_sets',
    'myfunc_get_included_resnames_sets',
    'myfunc_get_selected_resnums',
    'myfunc_get_included_resnums',
    'myfunc_get_selected_atoms_charges_separated_by_molecule',
    'myfunc_get_included_atoms_charges_separated_by_molecule',
    'myfunc_get_selected_atoms_charges_separated_by_entry',
    'myfunc_get_included_atoms_charges_separated_by_entry',
    'myfunc_ws_get_chosen_atoms_indexes',
    'myfunc_ws_get_chosen_atoms_indexes_detail',
    'myfunc_ws_get_chosen_atoms',
    'myfunc_ws_get_chosen_atoms_center_and_size',
]


def _get_rows(included=False):
    """1D: List[row, row, ...]"""
    pt = maestro.project_table_get()
    if included:
        rows = [i for i in pt.included_rows][-1::-1]        # important!
        print(f'Note: number of included entries: {len(rows)}')
    else:
        rows = [i for i in pt.selected_rows]
        print(f'Note: number of selected entries: {len(rows)}')
    info = [','.join([r.title,r.entry_id]) for r in rows]
    info = [s.replace('"','').replace("'",'') for s in info]
    info = '; '.join(['({:})'.format(s) for s in info])
    print(f'Note: -> key-pairs: (title,entry_id): {info}')
    return rows


def myfunc_get_selected_rows():
    return _get_rows()


def myfunc_get_included_rows():
    return _get_rows(included=True)


def myfunc_get_selected_titles():
    """1D: List[str, str, ...]"""
    rows = myfunc_get_selected_rows()
    return [i.title for i in rows]


def myfunc_get_included_titles():
    """1D: List[str, str, ...]"""
    rows = myfunc_get_included_rows()
    return [i.title for i in rows]


def myfunc_get_selected_structures():
    """1D: List[structure, structure, ...]"""
    rows = myfunc_get_selected_rows()
    #return [i.structure for i in rows]     # deprecated, properties are not copied!
    return [i.getStructure() for i in rows]


def myfunc_get_included_structures():
    """1D: List[structure, structure, ...]"""
    rows = myfunc_get_included_rows()
    return [i.getStructure() for i in rows]


def myfunc_get_selected_molecules():
    sts = myfunc_get_selected_structures()
    mols = [[m for m in s.molecule] for s in sts]
    print(f'Note: number of total molecules: {sum([s.mol_total for s in sts])}')
    return mols


def myfunc_get_included_molecules():
    sts = myfunc_get_included_structures()
    mols = [[m for m in s.molecule] for s in sts]
    print(f'Note: number of total molecules: {sum([s.mol_total for s in sts])}')
    return mols


def _get_atoms_separated_by_molecule(included=False):
    """3D: List[[[atom, atom, ...], [atom, atom, ...], ...], ...]"""
    if included:
        sts = myfunc_get_included_structures()
    else:
        sts = myfunc_get_selected_structures()
    atoms = [[[a for a in m.atom] for m in s.molecule] for s in sts]
    print(f'Note: number of total molecules: {sum([s.mol_total for s in sts])}')
    print(f'Note: number of total atoms: {sum([s.atom_total for s in sts])}')
    return atoms


def myfunc_get_selected_atoms_separated_by_molecule():
    return _get_atoms_separated_by_molecule()


def myfunc_get_included_atoms_separated_by_molecule():
    return _get_atoms_separated_by_molecule(included=True)


def _get_atoms_separated_by_entry(included=False):
    """2D: List[[atom, atom, ...], ...]"""
    if included:
        sts = myfunc_get_included_structures()
    else:
        sts = myfunc_get_selected_structures()
    atoms = [[a for a in s.atom] for s in sts]
    total = sum([len(s) for s in atoms])
    print(f'Note: number of total atoms: {total}')
    return atoms


def myfunc_get_selected_atoms_separated_by_entry():
    return _get_atoms_separated_by_entry()


def myfunc_get_included_atoms_separated_by_entry():
    return _get_atoms_separated_by_entry(included=True)


def myfunc_get_selected_entry_ids():
    """1D: List[str, str, ...]"""
    rows = myfunc_get_selected_rows()
    ids = [i.entry_id for i in rows]
    return ids


def myfunc_get_included_entry_ids():
    """1D: List[str, str, ...]"""
    rows = myfunc_get_included_rows()
    ids = [i.entry_id for i in rows]
    return ids


def myfunc_get_selected_resnames_lists():
    """1D: List[str, str, ...]"""
    sts = myfunc_get_selected_structures()
    resnames = [[k.pdbres for j in i.molecule for k in j.residue] for i in sts]
    print(f'Note: number of total selected residues: {sum([len(i) for i in resnames])}')
    return resnames


def myfunc_get_included_resnames_lists():
    """1D: List[str, str, ...]"""
    sts = myfunc_get_included_structures()
    resnames = [[k.pdbres for j in i.molecule for k in j.residue] for i in sts]
    print(f'Note: number of total included residues: {sum([len(i) for i in resnames])}')
    return resnames


def myfunc_get_selected_resnames_sets():
    """1D: List[str, str, ...]"""
    resnames = myfunc_get_selected_resnames_lists()
    return [list(set(i)) for i in resnames]


def myfunc_get_included_resnames_sets():
    """1D: List[str, str, ...]"""
    resnames = myfunc_get_included_resnames_lists()
    return [list(set(i)) for i in resnames]


def myfunc_get_selected_resnums():
    """1D: List[int, int, ...]"""
    sts = myfunc_get_selected_structures()
    resnums = [[k.resnum for j in i.molecule for k in j.residue] for i in sts]
    print(f'Note: number of total selected residues: {sum([len(i) for i in resnums])}')
    return resnums


def myfunc_get_included_resnums():
    """1D: List[int, int, ...]"""
    sts = myfunc_get_included_structures()
    resnums = [[k.resnum for j in i.molecule for k in j.residue] for i in sts]
    print(f'Note: number of total included residues: {sum([len(i) for i in resnums])}')
    return resnums


def myfunc_get_selected_atoms_charges_separated_by_molecule():
    """3D: List[[[float, float, ...], [float, float, ...], ...], ...]"""
    full = myfunc_get_selected_atoms_separated_by_molecule()
    charges = [[[a.partial_charge for a in m] for m in s] for s in full]
    return charges


def myfunc_get_included_atoms_charges_separated_by_molecule():
    """3D: List[[[float, float, ...], [float, float, ...], ...], ...]"""
    full = myfunc_get_included_atoms_separated_by_molecule()
    charges = [[[a.partial_charge for a in m] for m in s] for s in full]
    return charges


def myfunc_get_selected_atoms_charges_separated_by_entry():
    """2D: List[[float, float, ...], ...]"""
    full = myfunc_get_selected_atoms_separated_by_entry()
    charges = [[a.partial_charge for a in m] for m in full]
    return charges


def myfunc_get_included_atoms_charges_separated_by_entry():
    """2D: List[[float, float, ...], ...]"""
    full = myfunc_get_included_atoms_separated_by_entry()
    charges = [[a.partial_charge for a in m] for m in full]
    return charges


def myfunc_ws_get_chosen_atoms_indexes():
    """1D: List[int, int, ...]"""
    aids = maestro.selected_atoms_get()
    print(f'Note: workspace: number of chosen atoms: {len(aids)}')
    return [i-1 for i in aids]


def myfunc_ws_get_chosen_atoms():
    """2D: List[[atom, atom, ...], ...]"""
    atoms = myfunc_get_included_atoms_separated_by_entry()
    aids = myfunc_ws_get_chosen_atoms_indexes()
    full = [a for st in atoms for a in st]
    return [full[i] for i in aids]


def myfunc_ws_get_chosen_atoms_indexes_detail():
    """1D: List[(eid,mid,aid), (eid,mid,aid), ...]"""
    atoms = myfunc_ws_get_chosen_atoms()
    return [(a.entry_id,a.molecule_number_by_entry-1,a.index-1) for a in atoms]


def _calc_center_and_size(crds):
    # codes may be redundant, but efficient
    if not crds: return [(0.0,0.0,0.0), (0.0,0.0,0.0)]
    xsum = ysum = zsum = 0.0
    xmin = ymin = zmin = 1000000000000000000.
    xmax = ymax = zmax = -1000000000000000000.
    for r in crds:
        xsum += r[0]
        ysum += r[1]
        zsum += r[2]
        xmin = min(xmin,r[0])
        ymin = min(ymin,r[1])
        zmin = min(zmin,r[2])
        xmax = max(xmax,r[0])
        ymax = max(ymax,r[1])
        zmax = max(zmax,r[2])
    n = len(crds)
    return [(xsum/n,ysum/n,zsum/n), (xmax-xmin,ymax-ymin,zmax-zmin)]


def myfunc_ws_get_chosen_atoms_center_and_size():
    """2D: List[ CenterCoord(x,y,z), BoxSize(xlen,ylen,zlen) ]"""
    atoms = myfunc_ws_get_chosen_atoms()
    crds = [a.xyz for a in atoms]
    return _calc_center_and_size(crds)




