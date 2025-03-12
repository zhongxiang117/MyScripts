"""
My PyMol scripts

'PYMOL_PATH'    : '/my-build-dir/pymol/pymol_path'
'PYMOL_DATA'    : '/my-build-dir/pymol/pymol_path/data'
'PYMOL_SCRIPTS' : '/my-build-dir/pymol/pymol_path/scripts'

select  :   /object/segi/chain/resi/name
example :   /1m17/A/A/TYR`803/CA        "`"(backtick) means "at"


useful:
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
"""

from openbabel import pybel

from pymol import cmd
from pymol import stored
from pymol import selector


def xx_get_property_pythonic(sel='sele'):
    """
    get pdb property, default sel=sele

    >>> model = cmd.get_model(sel)
    # check `PYMOL_PATH/../../chempy/models.py`
    # e.g.: model.get_min_max()

    >>> atoms = model.atom      # List[chempy.Atom]
    # check `PYMOL_PATH/../../chempy/__init__.py`
    # e.g.:  a.resn, a.resi, a.coord, a.ss
    """
    if sel not in cmd.get_names('all'): sel = 'all'
    model = cmd.get_model(sel)
    atoms = model.atom
    for a in atoms:
        print(a.resn, a.resi, a.coord, a.ss)


def xx_get_property(sel='sele'):
    """
    get pdb property, default sel=sele

    PYMOL_PATH/../../pymol/editing.py:

    name, resn, resi, resv, chain, segi, elem, alt, q, b, vdw, type,
    partial_charge, formal_charge, elec_radius, text_type, label,
    numeric_type, model*, state*, index*, ID, rank, color, ss,
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


def xx_get_helix(sel='sele',save2file=1):
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

    if save2file:
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


def xx_get_phi_psi_omega(sel='sele',save2file=1):
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

    if save2file:
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
    cmd.iterate(sel, 'stored.p.append("%s and ID %s"%(model,ID))')
    if len(stored.p) % 2 != 0:
        print('Warning: number of pair should be even')
        return
    t = len(stored.p) // 2
    f = []
    for k,v in zip(stored.p[:t],stored.p[t:]): f.extend([k,v])
    print(', '.join(f))

cmd.extend('xx_get_info_for_pair_fit', xx_get_info_for_pair_fit)
cmd.auto_arg[0]['xx_get_info_for_pair_fit'] = cmd.auto_arg[0]['align']   # for auto-completion


def xx_get_molinfo(sel=None,solvent_radius=None):
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





