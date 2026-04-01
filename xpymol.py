"""
My PyMol scripts

'PYMOL_PATH'    : '/my-build-dir/pymol/pymol_path'
'PYMOL_DATA'    : '/my-build-dir/pymol/pymol_path/data'
'PYMOL_SCRIPTS' : '/my-build-dir/pymol/pymol_path/scripts'
"""

USAGE = """
select  :   /object/segi/chain/resi/name
example :   /1m17/A/A/TYR`803/CA        "`"(backtick) means "at"
defines :
    name, resn, resi, resv, chain, segi, elem, alt, q, b, vdw, type,
    partial_charge, formal_charge, elec_radius, text_type, label,
    numeric_type, model*, state*, index*, id (starts from 1), rank, color, ss,
    cartoon, flags


useful:  https://pymolwiki.org/index.php/Selection_Algebra
    select hydrophobes, (resn ala+gly+val+ile+leu+phe+met and !name C+N+O)
    select hydrophilics, (resn arg+lys+his+glu+asp+asn+gln+thr+ser+cys and !name C+N+O)
    select aromatics, (resn phe+tyr+trp+his and !name C+N+O)
    select positive, (resn arg+lys+his and !name C+N+O)
    select negative, (resn glu+asp and !name C+N+O)
    select acid, (resn asp+glu+cgu)
    select basic, (resn arg+lys+his)
    select polar, (resn ser+thr+asn+gln+tyr)
    select nopolar, ('resn met+phe+pro+trp+val+leu+ile+ala')
    select backbone, (name c+n+o+ca)


>>> model = cmd.get_model(sel)
# check `PYMOL_PATH/../../chempy/models.py`
# e.g.: model.get_min_max()

>>> atoms = model.atom      # List[chempy.Atom]
# check `PYMOL_PATH/../../chempy/__init__.py`
# e.g.:  a.resn, a.resi, a.coord, a.ss

    model = cmd.get_model(sel)
    atoms = model.atom
    for a in atoms:
        print(a.resn, a.resi, a.coord, a.ss)


# `stored` is import by default
stored.p = 0
iterate sele, stored.p += 1


set_name old_name, new_name

set cartoon_transparency, 0.5, <sele>

set label_size, 10
label sele, "x"
# "Editing" mode, `ctrl+left_click`, to move label


"""

from openbabel import pybel

import numpy as np
from pymol import cmd
from pymol import stored
from pymol import selector

import math


YES = [1, True, '1', 'T', 'True', 'TRUE', 'Y', 'Yes', 'YES', 'y', 'yes']
NO = [0, False, None, 'NONE', 'None', '0', 'F', 'False', 'FALSE', 'N', 'No', 'NO', 'n', 'no']


def xx_help_grammar(*args, **kws):
    print(USAGE)

cmd.extend('xx_help_grammar', xx_help_grammar)


def xx_get_property(sel='sele'):
    """
    get pdb property, default sel=sele

    PYMOL_PATH/../../pymol/editing.py:

    name, resn, resi, resv, chain, segi, elem, alt, q, b, vdw, type,
    partial_charge, formal_charge, elec_radius, text_type, label,
    numeric_type, model*, state*, index*, id (starts from 1), rank, color, ss,
    cartoon, flags
    """
    if sel not in cmd.get_names('all'): sel = 'all'
    stored.p = []
    cmd.iterate(sel, 'stored.p.append((name, resn, resi, alt, chain, numeric_type, ss))')
    for i in stored.p:
        print(i)

cmd.extend('xx_get_property', xx_get_property)
cmd.auto_arg[0]['xx_get_property'] = cmd.auto_arg[0]['align']   # for auto-completion


def xx_pdb2ss(sel='sele'):
    """
    convert sel to second_structure, default sel=sele
    """
    if sel not in cmd.get_names('all'): sel = 'all'
    stored.ss = ''
    cmd.iterate('%s and (n. CA)' % sel, 'stored.ss += ("%1s" % ss)')
    ss = stored.ss.replace(' ','.')
    print(ss)

cmd.extend('xx_pdb2ss', xx_pdb2ss)
cmd.auto_arg[0]['xx_pdb2ss'] = cmd.auto_arg[0]['align']   # for auto-completion


def xx_get_helix(sel='sele', save2file=True):
    """
    get distance between O[i] and N[i+3|i+4|i+5], save to file or print to stdout
    """
    if sel not in cmd.get_names('all'): sel = 'all'
    model = cmd.get_model(sel)
    reslist = model.get_residues()
    atomlist = model.atom
    resdata = []
    for start,end in reslist:
        nc = oc = i = None
        for atom in atomlist[start:end]:
            if atom.name == 'N':
                i = atom.resi
                nc = atom.coord
            elif atom.name == 'O':
                oc = atom.coord
        if nc and oc and i:
            resdata.append((oc,nc,i))
    if len(resdata) < 6:
        print('Warning: at least 6 residues needed')
        return
    rstlist = []
    for i in range(len(resdata)-5):
        t = resdata[i][2]
        oc = resdata[i][0]
        nc3 = resdata[i+3][1]
        nc4 = resdata[i+4][1]
        nc5 = resdata[i+5][1]
        on3 = pow(sum((j-k)**2 for j,k in zip(oc,nc3)), 0.5)
        on4 = pow(sum((j-k)**2 for j,k in zip(oc,nc4)), 0.5)
        on5 = pow(sum((j-k)**2 for j,k in zip(oc,nc5)), 0.5)
        rstlist.append((t,on3,on4,on5))

    if save2file in YES:
        with open('helix.dat','wt') as f:
            for g in rstlist:
                f.write('{:6}  {:7.3f}  {:7.3f}  {:7.3f}\n'.format(*g))
        print('Note: helix data saved to file: helix.dat')
        print(
            (
            "GNUPLOT> "
            "plot 'helix.dat' using 1:2 title 'd(O(i)-N(i+3))' with lines, "
            " 'helix.dat' using 1:3 title 'd(O(i)-N(i+4))' with lines, "
            " 'helix.dat' using 1:4 title 'd(O(i)-N(i+5))' with lines"
            )
        )
    else:
        for g in rstlist:
            print('{:6}  {:7.3f}  {:7.3f}  {:7.3f}'.format(*g))

cmd.extend('xx_get_helix', xx_get_helix)
cmd.auto_arg[0]['xx_get_helix'] = cmd.auto_arg[0]['align']   # for auto-completion


def xx_get_phi_psi_omega(sel='sele', save2file=True):
    """
    get phi/psi/omega, save to file or print to stdout
    """
    if sel not in cmd.get_names('all'): sel = 'all'
    stored.dict = {}
    cmd.iterate(sel , "stored.dict[int(resi)] = [model, segi, chain, resn, int(resi)]")
    keys = list(stored.dict.keys())
    rstlist = []
    for resi in keys:
        phi = 0
        if resi-1 in keys:
            g = stored.dict[resi-1]
            m = stored.dict[resi]
            s1 = '/{:}/{:}/{:}/{:}`{:}/C'.format(*g)
            s2 = '/{:}/{:}/{:}/{:}`{:}/N'.format(*m)
            s3 = '/{:}/{:}/{:}/{:}`{:}/CA'.format(*m)
            s4 = '/{:}/{:}/{:}/{:}`{:}/C'.format(*m)
            try:
                phi = cmd.get_dihedral(s1, s2, s3, s4, state=0)
            except:
                pass
        psi = 0
        if resi+1 in keys:
            g = stored.dict[resi]
            m = stored.dict[resi+1]
            s1 = '/{:}/{:}/{:}/{:}`{:}/N'.format(*g)
            s2 = '/{:}/{:}/{:}/{:}`{:}/CA'.format(*g)
            s3 = '/{:}/{:}/{:}/{:}`{:}/C'.format(*g)
            s4 = '/{:}/{:}/{:}/{:}`{:}/N'.format(*m)
            try:
                psi = cmd.get_dihedral(s1, s2, s3, s4, state=0)
            except:
                pass
        omege = 0
        if resi-1 in keys:
            g = stored.dict[resi-1]
            m = stored.dict[resi]
            s1 = '/{:}/{:}/{:}/{:}`{:}/CA'.format(*g)
            s2 = '/{:}/{:}/{:}/{:}`{:}/C'.format(*g)
            s3 = '/{:}/{:}/{:}/{:}`{:}/N'.format(*m)
            s4 = '/{:}/{:}/{:}/{:}`{:}/CA'.format(*m)
            try:
                omege = cmd.get_dihedral(s1, s2, s3, s4, state=0)
            except:
                pass
        if phi or psi or omege:
            rstlist.append((resi, phi, psi, omege))

    if save2file in YES:
        with open('phi_psi_omega.dat','wt') as f:
            f.write('resnum    phi       psi      omega\n')
            for g in rstlist:
                f.write('{:6}  {:8.3f}  {:8.3f}  {:8.3f}\n'.format(*g))
        print('Note: result saved to file: phi_psi_omega.dat')
    else:
        for g in rstlist:
            print('{:6}  {:8.3f}  {:8.3f}  {:8.3f}'.format(*g))

cmd.extend('xx_get_phi_psi_omega', xx_get_phi_psi_omega)
cmd.auto_arg[0]['xx_get_phi_psi_omega'] = cmd.auto_arg[0]['align']   # for auto-completion


def xx_get_info_for_pair_fit(sel='sele'):
    """
    selection should only be inside two objects, sequence does not matter,
    outputs can be directly used for `pair_fit`
    """
    stored.p = []
    #cmd.iterate(sel, 'stored.p.append("/%s/%s/%s/%s`%s/%s"%(model,segi,chain,resn,resi,name))')
    cmd.iterate(sel, 'stored.p.append("%s and id %s"%(model,id))')
    if len(stored.p) % 2 != 0:
        print('Warning: number of pair should be even')
        return
    t = len(stored.p) // 2
    f = []
    for k,v in zip(stored.p[:t],stored.p[t:]): f.extend([k,v])
    print(', '.join(f))

cmd.extend('xx_get_info_for_pair_fit', xx_get_info_for_pair_fit)
cmd.auto_arg[0]['xx_get_info_for_pair_fit'] = cmd.auto_arg[0]['align']   # for auto-completion


def xx_get_molinfo(sel=None, solvent_radius=None):
    """
    Formula, Number Atoms, M.W., Exact Mass, Solvent Accessible Surface, Solvent Excluded Surface

    `solvent_radius` can be directly input, or `set solvent_radius, 1.4`
    similar command to calculate surface `get_area`
    """
    def _getinfo(sel):
        stored.p = []
        cmd.iterate(sel, 'stored.p.append(elem)')
        d = {k:[0,0,0] for k in set(stored.p)}
        for c in stored.p:
            d[c][0] += 1
            i = pybel.ob.GetAtomicNum(c)
            d[c][1] += pybel.ob.GetMass(i)
            d[c][2] += pybel.ob.GetExactMass(i)
        info = ''
        mw = mass = 0.0
        n = 0
        for k in ['C','H','O','N','S','P','Cl']:
            if k in d:
                e = d.pop(k)
                n += e[0]
                info += '{:}{:}'.format(k,e[0])
                mw += e[1]
                mass += e[2]
        for k,e in d.items():
            n += e[0]
            info += '{:}{:}'.format(k,e[0])
            mw += e[1]
            mass += e[2]
        old = cmd.get('dot_solvent')
        cmd.set('dot_solvent','on')
        sas = cmd.get_area(sel)
        cmd.set('dot_solvent','off')
        ses = cmd.get_area(sel)
        cmd.set('dot_solvent',old)
        v = '{:}: {:},  NumAtoms: {:},  MW: {:.2f},  ExactMass: {:.3f},  SAS: {:.3f},  SES: {:.3f}'.format(
            sel, info, n, mw, mass, sas, ses
        )
        return v

    old = cmd.get('solvent_radius')
    if solvent_radius is None:
        solvent_radius = cmd.get('solvent_radius')
    else:
        solvent_radius = solvent_radius
    cmd.set('solvent_radius',solvent_radius)
    if sel is None:
        if 'sele' in cmd.get_names('all'):
            sel = 'sele'
        else:
            sel = cmd.get_names()       # not "all"
    info = []
    if isinstance(sel,list):
        for s in sel:
            info.append(_getinfo(s))
    else:
        info.append(_getinfo(sel))
    for i in info: print(i)
    cmd.set('solvent_radius',old)

cmd.extend('xx_get_molinfo', xx_get_molinfo)
cmd.auto_arg[0]['xx_get_molinfo'] = cmd.auto_arg[0]['align']   # for auto-completion


def xx_find_steric_clashes(sel1, sel2=None, distance_cutoff=1.0):
    """
    find steric clashes, if sel2 is input, find clashes between sel1 and sel2
    otherwise, find self clashes after excluding bonds,
    note, sel2 can be skipped by using 0, n, No
    default: sel2=None, distance_cutoff=1.0
    """
    second = sel2 if sel2 else sel1
    pairs = cmd.find_pairs(sel1, second, mode=0, cutoff=float(distance_cutoff))

    if sel2 in NO:
        members = ''
        for i, (a,b) in enumerate(pairs,1):
            sa = '{:} and id {:}'.format(sel1, a[1])
            sb = '{:} and id {:}'.format(second, b[1])
            cmd.distance(f'clash_{i}', sa, sb)
            members += f'clash_{i}  '
        name = cmd.get_unused_name('steric_clashes',False)
        cmd.group(name, members)
    else:
        pass

cmd.extend('xx_find_steric_clashes', xx_find_steric_clashes)
cmd.auto_arg[0]['xx_find_steric_clashes'] = cmd.auto_arg[0]['align']    # for auto-completion
cmd.auto_arg[1]['xx_find_steric_clashes'] = cmd.auto_arg[0]['align']    # for auto-completion


def xx_find_hydrogen_bonds(sel1, sel2, distance_cutoff=3.5, angle_cutoff=100):
    """
    find hydrogen bonds, only geometry is considered but not atom types or force fields
    default: distance_cutoff=3.5, angle_cutoff=100
    """
    pairs = cmd.find_pairs(sel1, sel2, mode=1, cutoff=float(distance_cutoff), angle=float(angle_cutoff))
    members = ''
    use = set()
    for i, (a,b) in enumerate(pairs,1):
        sa = '{:} and id {:}'.format(sel1, a[1])
        sb = '{:} and id {:}'.format(sel2, b[1])
        cmd.distance(f'hbond_{i}', sa, sb)
        members += f'hbond_{i}  '

        at = cmd.get_model(sa).atom[0]
        bt = cmd.get_model(sb).atom[0]
        if at.resi.strip():
            use.add('{:} and resi {:}'.format(sel1, at.resi))
        else:
            use.add(sa)
        if bt.resi.strip():
            use.add('{:} and resi {:}'.format(sel1, bt.resi))
        else:
            use.add(sb)

    name = cmd.get_unused_name('hbonds',False)
    cmd.group(name, members)

    news = ''
    for i,s in enumerate(use,1):
        n = 'r' + str(i)
        cmd.create(n, s, zoom=0)
        news += n + ' '
    name = cmd.get_unused_name('hbonds_frags',False)
    cmd.group(name, news)

cmd.extend('xx_find_hydrogen_bonds', xx_find_hydrogen_bonds)
cmd.auto_arg[0]['xx_find_hydrogen_bonds'] = cmd.auto_arg[0]['align']    # for auto-completion
cmd.auto_arg[1]['xx_find_hydrogen_bonds'] = cmd.auto_arg[0]['align']    # for auto-completion


def xx_find_aromatic_rings(sel, show=True, show_center=False):
    """
    find aromatic rings by using openbabel.pybel
    default: show=True, show_center=False
    """
    rings = []
    pdb = cmd.get_pdbstr(sel)
    mol = pybel.readstring('pdb', pdb)
    obmol = mol.OBMol
    rings = []
    for ring in obmol.GetSSSR():
        aids = list(ring._path)  # atom indices, starts from 1
        atoms = [obmol.GetAtom(i) for i in aids]
        if all([a.IsAromatic() for a in atoms]):

            xyzs = [(t.GetX(), t.GetY(), t.GetZ()) for t in [a.GetVector() for a in atoms]]

            sx = sum([i[0] for i in xyzs])
            sy = sum([i[1] for i in xyzs])
            sz = sum([i[2] for i in xyzs])
            centroid = [sx/len(xyzs), sy/len(xyzs), sz/len(xyzs)]

            u, v, w = xyzs[0], xyzs[1], xyzs[2]
            r = [u[0]-w[0], u[1]-w[1], u[2]-w[2]]
            t = [u[0]-v[0], u[1]-v[1], u[2]-v[2]]
            norm = [r[1]*t[2] - r[2]*t[1],  r[2]*t[0] - r[0]*t[2],  r[0]*t[1] - r[1]*t[0]]
            l = math.sqrt(sum([i*i for i in norm]))
            if l > 0.000001:
                norm = [i/l for i in norm]

            sele = sel + ' and id ' + '+'.join([str(i) for i in aids])

            rings.append({
                'ring_ids': aids,
                'selector': sele,
                'centroid': centroid,
                'norm': norm,
            })

    if show in YES:
        members = ''
        for i,g in enumerate(rings,1):
            cmd.select(f'ring_{i}', g['selector'])
            members += f'ring_{i}  '
        name = cmd.get_unused_name('aromatic_rings',False)
        cmd.group(name, members)

    if show_center in YES:
        members = ''
        for i,g in enumerate(rings,1):
            cmd.pseudoatom(f'cent_{i}', pos=g['centroid'])
            members += f'cent_{i}  '
        name = cmd.get_unused_name('aromatic_centers',False)
        cmd.group(name, members)

    print('XX: number of aromatic rings: {:}:  {:}'.format(sel, len(rings)))

    return rings

cmd.extend('xx_find_aromatic_rings', xx_find_aromatic_rings)
cmd.auto_arg[0]['xx_find_aromatic_rings'] = cmd.auto_arg[0]['align']    # for auto-completion


def xx_find_pi_pi_stacks(sel1, sel2, min_dist=3.5, max_dist=5.5, angle=30.0):
    """
    find pi-pi stacking
    default: min_dist=3.5, max_dist=5.5, angle=30.0
    """
    min_dist = float(min_dist)
    max_dist = float(max_dist)

    rings1 = xx_find_aromatic_rings(sel1, False, False)
    rings2 = xx_find_aromatic_rings(sel2, False, False)

    pis = []
    for r1 in rings1:
        for r2 in rings2:
            c1, n1 = r1['centroid'], r1['norm']
            c2, n2 = r2['centroid'], r2['norm']

            r = [c1[0]-c2[0], c1[1]-c2[1], c1[2]-c2[2]]
            dist = math.sqrt(sum([i*i for i in r]))
            if not (min_dist <= dist and dist <= max_dist): continue

            u = sum([n1[0]*n2[0], n1[1]*n2[1], n1[2]*n2[2]])
            ang = math.degrees(math.acos(min(abs(u), 1.0)))     # ensure [0, 180]
            if ang <= angle or (150 <= ang and ang <= 180):     # stacked or flipped
                pis.append((r1,r2))

    if pis:
        members = ''
        for i, (r1,r2) in enumerate(pis,1):
            cmd.pseudoatom(f'cent1_{i}', pos=r1['centroid'])
            cmd.pseudoatom(f'cent2_{i}', pos=r2['centroid'])
            cmd.distance(f'pp_{i}', f'cent1_{i}', f'cent2_{i}')
            cmd.group(f'pi_pi_{i}', f'cent1_{i}  cent2_{i}  pp_{i}')
            members += f'pi_pi_{i}  '
        name = cmd.get_unused_name('pi_pi_stacks',False)
        cmd.group(name, members)

    print('XX: number of pi-pi stacks: {:}'.format(len(pis)))

cmd.extend('xx_find_pi_pi_stacks', xx_find_pi_pi_stacks)
cmd.auto_arg[0]['xx_find_pi_pi_stacks'] = cmd.auto_arg[0]['align']  # for auto-completion
cmd.auto_arg[1]['xx_find_pi_pi_stacks'] = cmd.auto_arg[0]['align']  # for auto-completion


def xx_split_to_group(group=None, sel=None):
    """
    split selection to group by using predefined keywords
    => chains  residues  ligands  solvents  metals  others  all

    link: https://pymolwiki.org/index.php/Selection_Algebra
    """
    if not group: group = 'all'
    if not sel: sel = 'sele'

    defined = {
        'Chains':'polymer', 'Ligands':'organic', 'Solvents':'solvent', 'Metals':'metals',
    }
    special = {'Hetatoms':'hetatm', 'Others':'others', 'All':'all', 'Residues':'polymer', }
    pdefs = [k.lower() for k in [*defined.keys(), *special.keys()]]
    if group not in pdefs:
        print('XX: not valid: ', group)
        print(' -> ', ' '.join(pdefs))
        return

    sg = group.capitalize()
    kg = defined[sg]
    sg = f'{sel}-{sg}'

    if group == 'others':
        kg = ' and '.join(['not '+v for v in defined.values()])
        sg = f'{sel}-Others'
    elif group == 'all':
        for g in [k.lower() for k in defined.keys()]:
            xx_split_to_group(group=g, sel=sel)
    else:
        u = cmd.get_unused_name()
        cmd.select(u, f'{kg} and {sel}')
        model = cmd.get_model(u)
        atoms = model.atom
        if group == 'chains':
            chains = {}
            for g in model.get_residues():
                ri = atoms[g[0]].resi
                fm = '{:}_{:}-{:}n'.format(atoms[g[0]].resn, ri, g[1]-g[0])
                cmd.create(fm, f'resi {ri} and {sel}', zoom=0)
                rc = atoms[g[0]].chain
                chains.setdefault(rc, []).append(fm)
            if chains:
                gc = []
                for c,rs in chains.items():
                    if not c: c = 'X'
                    print('XX: creating chain: ', f'{sg}-{c}')
                    s = '{:}-{:}n'.format(c, len(rs))
                    cmd.group(s, ' '.join(rs))
                    gc.append(s)
                cmd.group(sg, ' '.join(gc))

                for c,rs in chains.items():
                    if not c: continue
                    ch = 'Chain-{:}-{:}n'.format(c,len(rs))
                    cmd.create(ch, f'chain {c} and {sel}', zoom=0)
                    cmd.group(sg, ch)
            else:
                print('XX: ignoring group: ', sg)
        else:
            names = []
            for g in model.get_residues():
                ri = atoms[g[0]].resi
                fm = '{:}_{:}-{:}n'.format(atoms[g[0]].resn, ri, g[1]-g[0])
                cmd.create(fm, f'resi {ri} and {sel}', zoom=0)
                names.append(fm)
            if names:
                print('XX: creating group: ', sg)
                s = '{:}-{:}n'.format(sg, len(names))
                cmd.group(s, ' '.join(names))
            else:
                print('XX: ignoring group: ', sg)
        cmd.delete(u)

cmd.extend('xx_split_to_group', xx_split_to_group)
cmd.auto_arg[0]['xx_split_to_group'] = [
    lambda: cmd.Shortcut(['chains', 'residues', 'ligands', 'solvents', 'metals', 'others', 'all']),
    '1st',
    ', '
]
cmd.auto_arg[1]['xx_split_to_group'] = cmd.auto_arg[0]['align']  # for auto-completion


def xx_expand_selection_by(sel=None, distance=3.0, expand_type=None, create=None, all_visiable_objects=None):
    """
    expand `sel` by x Angstrom for self or for all GUI visible objects
    """
    if not sel: sel = 'sele'
    distance = float(distance)
    distsq = distance * distance
    if not expand_type:
        expand_type = 'atom'
    elif expand_type.lower() in ['a', 'atom']:
        expand_type = 'atom'
    elif expand_type.lower() in ['r', 'residue']:
        expand_type = 'residue'
    else:
        print('Warning: not valid `expand_type` -- use default `atom` instead')
        expand_type = 'atom'

    if all_visiable_objects in YES:
        names = cmd.get_names(enabled_only=True)
    else:
        stored.names = set()
        cmd.iterate(sel, 'stored.names.add(model)')
        names = [*list(stored.names)]
    print('XX: using: ', '  '.join(names))

    model = cmd.get_model(sel)
    if len(model.atom) == 0:
        print('Fatal: no atoms -- empty selection')
        return

    xyzref = np.array([a.coord for a in model.atom])
    xyzref_exp = xyzref[:, np.newaxis, :]
    gets = []

    info = ''
    if expand_type == 'atom':
        for n in names:
            model = cmd.get_model(n)
            atoms = model.atom

            if len(atoms) == 0:
                print('XX: ignoring: not valid molecule:', n)
                continue

            xyzcmp = np.array([a.coord for a in atoms])
            xyzcmp = xyzcmp[np.newaxis, :, :]
            sdist = np.sum((xyzref_exp - xyzcmp) ** 2, axis=2)
            d = np.where(sdist <= distsq)[1]
            index = [a.index for a in [atoms[i] for i in d]]
            gets.extend([f'(index {i} and {n})' for i in set(index)])
        info = 'XX: numbers atoms: {:}'.format(len(gets))

    elif expand_type == 'residue':
        for n in names:
            model = cmd.get_model(n)
            atoms = model.atom
            for g in model.get_residues():
                xyzcmp = np.array([atoms[i].coord for i in range(g[0],g[1])])
                xyzcmp = xyzcmp[np.newaxis, :, :]
                sdist = np.sum((xyzref_exp - xyzcmp) ** 2, axis=2)
                if np.min(sdist) <= distsq:
                    ri = atoms[g[0]].resi
                    gets.append(f'(resi {ri} and {n})')
        info = 'XX: numbers residues: {:}'.format(len(gets))

    print(info)
    if create in YES:
        n = cmd.get_unused_name(prefix='expd')
        print('XX: creating: ', n)
        cmd.create(n, ' or '.join(gets), zoom=0)
    else:
        n = cmd.get_unused_name(prefix='selc')
        print('XX: selecting: ', n)
        cmd.select(n, ' or '.join(gets))

cmd.extend('xx_expand_selection_by', xx_expand_selection_by)
cmd.auto_arg[0]['xx_expand_selection_by'] = cmd.auto_arg[0]['align']  # for auto-completion
cmd.auto_arg[2]['xx_expand_selection_by'] = [
    lambda: cmd.Shortcut(['atom', 'residue']),
    '3rd',
    ', '
]

import zipimport
plip_pymol_plugin_xz = None
try:
    ip = zipimport.zipimporter('/home/xiang/.pymol/startup/pymol_plip.zip')
    ip.load_module('plip')
    from plip.xpymol import plip_pymol_plugin_xz
except ImportError:
    plip_pymol_plugin_xz = None

if plip_pymol_plugin_xz:
    def xx_plip(obj='', args=''):
        """
Usage:  xx_plip  obj, ddd=1 aap=3 vvv       # use space
check: -h / --help
"""
        ss = obj + args
        if '-h' in ss or '--help' in ss:
            args = '-h'
            pdbstr = None
        else:
            pdbstr = cmd.get_pdbstr(obj)
        plip_pymol_plugin_xz(None,pdbstr,args)

    cmd.extend('xx_plip', xx_plip)
    cmd.auto_arg[0]['xx_plip'] = cmd.auto_arg[0]['align']  # for auto-completion







