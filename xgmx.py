#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import Draw

import os
import sys
import collections
import argparse


FEATURES = [
    'version 0.1.0  : GROMACS Topology Parser',
]

VERSION = FEATURES[-1].split()[1]
__version__ = VERSION




class GMXTopParBase:
    """
    Args:
        sections(dict):
            {'directive':List[pars], 'comment':str, 'gro':List[pars], 'xyz':List[float] }

        >> gro: List[[resnum, resname, atomname, x, y, z, vx, vy, vz]] : (x:8.3f, vx:8.4f)
        >> xyz: List[[atomname, x, y, z]]
    """
    def __init__(self,sections=None,*args,**kws):
        self.sections = sections
        self.kws = kws

    def __str__(self):
        keys = list(self.sections.keys())
        for k in [s for s in keys if s.endswith('type')]:
            keys.remove(k)
        for k in ['defaults','systems','molecule']:
            if k in keys: keys.remove(k)

        out = '; GMXTopParser -- xz\n\n'
        if 'comment' in keys:
            keys.remove('comment')
            out += self.sections['comment'] + '\n'

        out += '[moleculetype]\n; Name        nrexcl\n'
        if 'moleculetype' in keys:
            keys.remove('moleculetype')
            for p in self.sections['moleculetype']:
                out += '  {:12s} {:}\n'.format(p[0],p[1])
        else:
            out += '  UNK            3\n'
        
        out += '\n[atoms]\n'
        out += ';   nr      type        resnr   residue    atom     cgnr      charge        mass'
        if 'atoms' in keys:
            keys.remove('atoms')
            for p in self.sections['atoms']:
                l = '{:6} {:14s} {:^6}  {:^8} {:^6}'.format(*p[:5])
                l += '  '.join(['{:10}'.format(i) for i in p[5:]])
                out += l + '\n'
        
        if 'bonds' in keys:
            keys.remove('bonds')
            for g in self._seperate(self.sections['bonds'],2):
                out += '\n[bonds]\n;    i     j     f\n'
                for p in g:
                    l = '{:6}{:6}{:6}'.format(*p[:3])
                    l += '  '.join(['{:10}'.format(i) for i in p[3:]])
                    out += l + '\n'

        if 'angles' in keys:
            keys.remove('angles')
            for g in self._seperate(self.sections['angles'],3):
                out += '\n[angles]\n;    i     j     k     f\n'
                for p in g:
                    l = '{:6}{:6}{:6}{:6}'.format(*p[:4])
                    l += '  '.join(['{:10}'.format(i) for i in p[4:]])
                    out += l + '\n'

        if 'dihedrals' in keys:
            keys.remove('dihedrals')
            for g in self._seperate(self.sections['dihedrals'],4):
                out += '\n[dihedrals]\n;    i     j     k     l     f\n'
                for p in g:
                    l = '{:6}{:6}{:6}{:6}{:6}'.format(*p[:5])
                    l += '  '.join(['{:10}'.format(i) for i in p[5:]])
                    out += l + '\n'

        if 'pairs' in keys:
            keys.remove('pairs')
            for g in self._seperate(self.sections['pairs'],2):
                out += '\n[pairs]\n;    i     j     f\n'
                for p in g:
                    l = '{:6}{:6}{:6}'.format(*p[:3])
                    l += '  '.join(['{:10}'.format(i) for i in p[3:]])
                    out += l + '\n'

        for k in keys:
            out += f'[{k}]\n'
            for p in self.sections[k]:
                out += '    '.join([str(i) for i in p]) + '\n'
            out += '\n'
        return out

    def _seperate(self,pars,idx):
        if not pars: return []
        blocks = list(set([p[idx] for p in pars]))
        out = [[] for i in blocks]
        k = {i:j for i,j in zip(blocks,range(len(blocks)))}
        for p in pars:
            n = k[p[idx]]
            out[n].append(p)
        return out

    def get_mol(self):
        if 'gro' in self.sections:
            mol = self.get_mol_from_gro()
        elif 'xyz' in self.sections:
            mol = self.get_mol_from_xyz()
        else:
            mol = self.get_mol_from_top()
        return mol

    def get_mol_from_top(self):
        if not ('atoms' in self.sections and 'bonds' in self.sections): return
        from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
        mol = Chem.RWMol()
        for g in self.sections['atoms']:
            a = self.guess_elemtype(g[4])
            mol.AddAtom(Chem.Atom(a))
        bonds = []
        for g in self.sections['bonds']:
            mol.AddBond(g[0]-1,g[1]-1)
            bonds.append((g[0]-1,g[1]-1,g[3]*10))   # unit nanometer to angstrom
        mol = mol.GetMol()
        if mol.GetNumAtoms() <= 10:
            n = 20
        elif mol.GetNumAtoms() <= 30:
            n = 50
        else:
            n = 100
        EmbedMultipleConfs(mol,numConfs=n,maxAttempts=500)
        rmsd = []
        for conf in mol.GetConformers():
            tot = 0.0
            for p in bonds:
                d = conf.GetAtomPosition(p[0]).Distance(conf.GetAtomPosition(p[1]))
                tot += (d-p[2]) ** 2
            rmsd.append(tot)
        t = rmsd.index(min(rmsd))
        sel = mol.GetConformer(t)
        new = Chem.RWMol(mol)
        new.RemoveAllConformers()
        new.AddConformer(sel)
        mol = new.GetMol()
        if 'moleculetype' in self.sections:
            name = self.sections['moleculetype'][0][0]
            mol.SetProp('_Name',name)
        mol.SetProp('Title',self.kws['title']) if 'title' in self.kws else mol.SetProp('Title','GMXTop')
        for i,a in enumerate(mol.GetAtoms()):
            a.SetProp('name',self.sections['atoms'][i][4])
            a.SetProp('residue',self.sections['atoms'][i][3])
        return mol

    def get_mol_from_gro(self,section=None):
        if not section:
            if not 'gro' in self.sections: return
            section = self.sections['gro']
        pdb = ''
        for i,g in enumerate(section):
            r = self.guess_residuetype(g[1])
            a = self.guess_elemtype(g[2])
            l = 'ATOM  {:5} {:3}  {:3}  {:4}   {:>8.3f} {:>8.3f} {:>8.3f}'.fromat(
                i+1,a,r[:3],g[1],g[4],g[5],g[6]
            )
            pdb += l + '\n'
        pdb += 'END\n'
        mol = Chem.MolFromPDBBlock(pdb,removeHs=False,sanitize=False)
        if 'moleculetype' in self.sections:
            name = self.sections['moleculetype'][0][0]
            mol.SetProp('_Name',name)
        mol.SetProp('Title',self.kws['title']) if 'title' in self.kws else mol.SetProp('Title','GMXTop')
        for i,a in enumerate(mol.GetAtoms()):
            a.SetProp('name',section[i][2])
        return mol

    def get_mol_from_xyz(self,section=None):
        if not section:
            if not 'xyz' in self.sections: return
            section = self.sections['xyz']
        xyz = '{:}\nGMXTop\n'.format(len(section))
        for g in section:
            a = self.guess_elemtype(g[0])
            l = '{:}  {:}  {:}  {:}'.format(a,g[1],g[2],g[3])
            xyz += l + '\n'
        xyz += '\n\n'
        mol = Chem.MolFromXYZBlock(xyz,removeHs=False,sanitize=False)
        if 'moleculetype' in self.sections:
            name = self.sections['moleculetype'][0][0]
            mol.SetProp('_Name',name)
        mol.SetProp('Title',self.kws['title']) if 'title' in self.kws else mol.SetProp('Title','GMXTop')
        for i,a in enumerate(mol.GetAtoms()):
            a.SetProp('name',section[i][0])
        return mol

    def guess_elemtype(self,s):
        s = s.strip()
        if not s:
            return Chem.Atom(0)
        try:
            atom = Chem.Atom(s[:2].strip())
        except:
            try:
                atom = Chem.Atom(s[0])
            except:
                atom = Chem.Atom(0)
        return atom.GetSymbol()

    def guess_residuetype(self,s):
        s = s.strip()
        if not s:
            return 'UNK'
        if len(s) <= 2:
            return s
        aa = [
            'ASP', 'ALA', 'ARG', 'ASN', 'CYS', 'GLN', 'GLU', 'GLY', 'ILE', 'LEU',
            'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HIS',
        ]
        if s.upper() in aa:
            return s.upper()
        if s[:3].upper() in aa:
            return s[:3].upper()
        return s

    def file_to_pdb(self,file=None):
        mol = self.get_mol()
        if not mol:
            print('Warning: cannot generate PDB file')
            return
        new = self.genfnew(file,ext='.pdb')
        Chem.MolToPDBFile(mol,new)

    def file_to_gro(self,file=None):
        if 'gro' in self.sections:
            grosection = self.sections['gro']
        elif 'atoms' in self.sections:
            mol = self.get_mol_from_gro()
            if not mol:
                print('Warning: cannot generate GRO file')
                return
            grosection = []
            conf = mol.GetConformer()
            g = self.sections['atoms']
            for i,p in enumerate(conf.GetPositions()):
                grosection.append((g[i][2],g[i][3],g[i][4],p[0],p[1],p[2]))
        elif 'xyz' in self.sections:
            grosection = []
            g = self.sections['atoms']
            for i,p in enumerate(self.sections['xyz']):
                grosection.append((g[i][2],g[i][3],g[i][4],p[0],p[1],p[2]))
        else:
            print('Warning: cannot generate GRO file')
            return

        if 'moleculetype' in self.sections:
            out = self.sections['moleculetype'][0][0]
        else:
            out = self.sections['title'] if 'title' in self.kws else 'GMXTopParser'
        out += '\n' + str(len(grosection)) + '\n'
        for g in grosection:
            out += '{:5}{:<5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}'.format(*g[:7])
            if len(g) > 7:
                out += '{:8.4f}{:8.4f}{:8.4f}'.format(g[7],g[8],g[9])
            out += '\n'
        if 'grobox' in self.kws:
            out += self.kws['grobox']
        else:
            xl = min([i[4] for i in grosection])
            xh = max([i[4] for i in grosection])
            yl = min([i[5] for i in grosection])
            yh = max([i[5] for i in grosection])
            zl = min([i[6] for i in grosection])
            zh = max([i[6] for i in grosection])
            l = '{:8.3f}{:8.3}{:8.3f}'.format(xh-xl,yh-yl,zh-zl)
            out += l + '\n'
            self.kws['grobox'] = l

        new = self.genfnew(file,ext='.gro')
        with open(new,'wt') as f:
            f.write(out)

    def file_to_itp(self,file=None):
        new = self.genfnew(file,ext='.itp')
        with open(new,'wt') as f:
            f.write(str(self))

    def genfnew(self,file=None,ext=None):
        if file:
            if os.path.isfile(file): return file
            file = '_'.join(file.split())
        else:
            file = 'gmx_base'
        if ext and not file.endswith(ext):
            file += ext
        i = 1
        while True:
            new = 'gmx-'+str(i)+'-'+file
            if not os.path.isfile(new): break
            i += 1
        print(f'Note: generating file: {new}')
        return new


class GMXTopParser:
    """GROMCAS Topology Parser

    For directive parser:
        they will be maximumly in a format: i j k l f a b c d g h m
        where, [ijkl] are the atom indexes, [f] is the fuction type, [abcdgh] are values,
        [m] is the comments

    > bonds:        ijfab, ijfabc
    > pairs:        ijfab, ijfabcd
    > angles:       ijkfab,ijfabcd
    > dihedrals:    ijklfab, ijklfabcdgh, ijklfabc
    > exclusions:   i

    Args:
        sections(dict): example:
            {"atoms":[line1, line2, line3, ...], "bonds":[line1, line2, line3, ...]}
    """
    supported_miscellaneous = ('default','molecule','system')
    supported_nbtypes = (
        'moleculetype','atomtypes','bondtypes','angletypes','dihedraltypes','constrainttypes',
    )
    supported_directives = (
        'atoms','bonds','pairs','angles','dihedrals','impropers','exclusions',
    )
    def __init__(self,sections=None,*args,**kws):
        self.sections = sections


    def run(self):
        pass
    def parse_atoms(self):
        pass
    def parse_bonds(self):
        pass
    def parse_pairs(self):
        pass
    def parse_pairs_nb(self):
        pass
    def parse_angles(self):
        pass
    def parse_dihedrals(self):
        pass
    def parse_impropers(self):
        pass
    def parse_exclusions(self):
        pass
    def parse_constraints(self):
        pass
    def parse_settles(self):
        pass
    def parse_virtual_sites2(self):
        pass
    def parse_virtual_sites3(self):
        pass
    def parse_virtual_sites4(self):
        pass
    def parse_virtual_sitesn(self):
        pass
    def parse_position_restraints(self):
        pass
    def parse_distance_restraints(self):
        pass
    def parse_dihedral_restraints(self):
        pass
    def parse_orientation_restraints(self):
        pass
    def parse_angle_restraints(self):
        pass






    def parse_block(self,oplists,lines,errmsg=None):
        """
        oplists(str): seperated by ";"
        lines(list) : 1D, List[str]
        """
        if isinstance(oplists,str):
            g = [i for i in oplists.split(';') if i.strip()]
            oplists = [self._get_oplist(i) for i in g]
        elif not isinstance(oplists,list):
            raise NotImplementedError(f'Not Supported: OP {oplists}')
        if not isinstance(lines,list):
            raise NotImplementedError(f'Not Supported: input: {lines}')
        if errmsg is None: errmsg = ''
        results = []
        ns = [len(c) for c in oplists]
        for l in lines:
            p = l.split(';')
            t = p[0].split()
            if not len(t): continue
            if len(t) not in (ns):
                o = errmsg+': '+ l if errmsg else l
                print(f'Warning: length in directives are wrong: {o}')
                return []
            b = oplists[ns.index(len(t))]    # get oplist
            g = []
            for i,v in enumerate(t):
                try:
                    k = b[i](v)
                except:
                    print(f'Warning: cannot parse <{v}> in line: {l}')
                    return []
                else:
                    g.append((k,p[1].strip()))
            results.append(g)
        return results

    def _get_oplist(self,opcode):
        """opcode will be the combinations of char [isr]"""
        if not isinstance(opcode,str):
            raise NotImplementedError(f'Not Supported: OP {opcode}')
        oplist = []
        for v in opcode:
            if v == 'i':
                oplist.append(int)
            elif v == 'r':
                oplist.append(float)
            elif v == 's':
                oplist.append(str)
            else:
                raise NotImplementedError(f'Not Supported: OP {v}')
        return oplist



def parsecmd():
    parser = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         usage=f"""
PROGRAM {VERSION}. General Usage:

 1. explanation
 >>> script
""",
        allow_abbrev=False,
    ),
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION
    )
    parser.add_argument(
        '-f','--file',
        dest='file',
        help='input file'
    )
    g = parser.add_argument_group('options for development or miscellaneous')
    g.add_argument(
        '--show_supported_ftypes',
        action='store_true',
        help='show supported file output types',
    )
    parser.add_argument(
        '--features',
        action='store_true',
        help='show development features'
    )

    if len(sys.argv) == 1:
        parser.print_help()
        return
    w = parser.parse_args(sys.argv[1:])
    if w.features:
        for i in FEATURES: print(i)
        return



