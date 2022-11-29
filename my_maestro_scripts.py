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
    'version 0.5.0  : make sure properties are copied when getting structures',
    'version 0.6.0  : add more funcs and make their names more clear',
    'version 0.7.0  : add correspondent `included` functions',
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
#    'myfunc_get_selected_molecules',
    'myfunc_get_selected_atoms_separated_by_molecule',
    'myfunc_get_included_atoms_separated_by_molecule',
    'myfunc_get_selected_atoms_separated_by_entry',
    'myfunc_get_included_atoms_separated_by_entry',
    'myfunc_get_selected_ids',
    'myfunc_get_included_ids',
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
    'myfunc_ws_get_chosen_atoms_ids',
    'myfunc_ws_get_chosen_atoms_ids_detail_from_selected',
    'myfunc_ws_get_chosen_atoms_ids_detail_from_included',
    'myfunc_ws_get_chosen_atoms_from_selected',
    'myfunc_ws_get_chosen_atoms_from_included',
    'myfunc_ws_get_chosen_atoms_center_and_size_from_selected',
    'myfunc_ws_get_chosen_atoms_center_and_size_from_included',
]


def _get_rows(included=False):
    """1D: List[row, row, ...]"""
    pt = maestro.project_table_get()
    if included:
        rows = [i for i in pt.included_rows]
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


# BUG -> Reported
#def myfunc_get_selected_molecules():
#    """1D: List[mol, mol, ...]"""
#    sts = myfunc_get_selected_structures()
#    mols = [[m for m in s.molecule] for s in sts]
#    print(f'Note: number of total selected molecules: {sum([len(i) for i in mols])}')
#    return mols


# Due to BUG
#def myfunc_get_selected_atoms_separated_by_molecule():
#    """3D: List[[[atom, atom, ...], [atom, atom, ...], ...], ...]"""
#    mols = myfunc_get_selected_molecules()
#    atoms = [[[i for i in m.atom] for m in s] for s in mols]
#    total = sum([sum([len(m) for m in s]) for s in atoms])
#    print(f'Note: number of total selected atoms: {total}')
#    return atoms


def _get_atoms_separated_by_molecule(included=False):
    """3D: List[[[atom, atom, ...], [atom, atom, ...], ...], ...]"""
    if included:
        sts = myfunc_get_included_structures()
    else:
        sts = myfunc_get_selected_structures()
    fullatoms = [[a for a in s.atom] for s in sts]
    # Care: using `resnum` rather than `molecule_number` for differentiating
    eatnums = [[a.resnum for a in s] for s in fullatoms]
    entryatoms = []
    molnum = 0
    for mnums,atoms in zip(eatnums,fullatoms):
        mids = {}
        keys = []       # to store sequence
        for i,j in enumerate(mnums):
            if j in mids:
                mids[j].append(i)
            else:
                mids[j] = [i,]
                keys.append(j)
                molnum += 1
        entryatoms.append([[atoms[v] for v in mids[k]] for k in keys])
    print(f'Note: number of total molecules: {molnum}')
    print(f'Note: number of total atoms: {sum([len(i) for i in fullatoms])}')
    return entryatoms


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


def myfunc_get_selected_ids():
    """1D: List[str, str, ...]"""
    rows = myfunc_get_selected_rows()
    ids = [i.entry_id for i in rows]
    return ids


def myfunc_get_included_ids():
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


def myfunc_ws_get_chosen_atoms_ids_detail_from_selected():
    """1D: List[(eid,mid,aid), (eid,mid,aid), ...]"""
    rows = myfunc_get_selected_atoms_separated_by_molecule()
    mids = myfunc_ws_get_chosen_atoms_ids()
    return _calc_ids(rows,mids)


def myfunc_ws_get_chosen_atoms_ids_detail_from_included():
    """1D: List[(eid,mid,aid), (eid,mid,aid), ...]"""
    rows = myfunc_get_included_atoms_separated_by_molecule()
    mids = myfunc_ws_get_chosen_atoms_ids()
    return _calc_ids(rows,mids)


def myfunc_ws_get_chosen_atoms_from_selected():
    """2D: List[[atom, atom, ...], ...]"""
    rows = myfunc_get_selected_atoms_separated_by_molecule()
    mids = myfunc_ws_get_chosen_atoms_ids()
    ids = _calc_ids(rows,mids)
    return [rows[i[0]][i[1]][i[2]] for i in ids]


def myfunc_ws_get_chosen_atoms_from_included():
    """2D: List[[atom, atom, ...], ...]"""
    rows = myfunc_get_included_atoms_separated_by_molecule()
    mids = myfunc_ws_get_chosen_atoms_ids()
    ids = _calc_ids(rows,mids)
    return [rows[i[0]][i[1]][i[2]] for i in ids]


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


def myfunc_ws_get_chosen_atoms_center_and_size_from_selected():
    """2D: List[ CenterCoord(x,y,z), BoxSize(xlen,ylen,zlen) ]"""
    atoms = myfunc_ws_get_chosen_atoms_from_selected()
    crds = [a.xyz for a in atoms]
    return _calc_center_and_size(crds)


def myfunc_ws_get_chosen_atoms_center_and_size_from_included():
    """2D: List[ CenterCoord(x,y,z), BoxSize(xlen,ylen,zlen) ]"""
    atoms = myfunc_ws_get_chosen_atoms_from_included()
    crds = [a.xyz for a in atoms]
    return _calc_center_and_size(crds)




