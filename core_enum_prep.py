#!/usr/bin/env python3

from rdkit import Chem

import os
import itertools
import collections

FEATURES = [
    'version 0.1.0  : CoreEnumeration, Sep 19th, 2022',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()


class CoreEnumeration:
    def __init__(self,core_file=None,sub_files=None,filename=None,*args,**kws):
        self.filename = filename if filename else 'mols-core-enum.sdf'
        if os.path.isfile(core_file):
            self.nice = True
            core = Chem.MolFromMolFile(core_file)
            core_bondedinfo = self.calc_dummy_atoms_idx_and_bonded_info(core)
            self.core_dummynum = len(core_bondedinfo)
            self.core_smarts = Chem.MolToSmarts(core)
            if isinstance(sub_files,str): sub_files = [sub_files, ]
            pro_sub_files = []
            for i in sub_files:
                if os.path.isfile(i):
                    pro_sub_files.append(i)
                else:
                    print(f'Warning: not a file: ignoring: {i}')
            subs = [m for f in pro_sub_files for m in Chem.SDMolSupplier(f)]
            subs_bondedinfos = [self.calc_dummy_atoms_idx_and_bonded_info(m) for m in subs]
            # only keep substructures that have one dummy atom
            reflist = [i for i in range(len(subs)) if len(subs_bondedinfos[i])==1]
            subs_bondedinfos = [subs_bondedinfos[i] for i in reflist]
            self.subs_smartss = [Chem.MolToSmarts(subs[i]) for i in reflist]
            if not self.core_dummynum:
                print('Fatal: core_file: no dummy atoms found: exiting')
                self.nice = False
            if not self.subs_smartss:
                print('Fatal: sub_files: no dummy atoms found: exiting')
                self.nice = False
        else:
            print(f'Fatal: not a file: exiting: {core_file}')
            self.nice = False

    def __next__(self):
        if not self.nice: raise StopIteration
        return next(self.run())

    def __iter__(self):
        return self

    def calc_dummy_atoms_idx_and_bonded_info(self,mol):
        """calculate dummy atoms' indexes and their bonded atoms info

        Args:
            mol (file|Chem.Mol):

        Return:
            dict(list(tuple(int,int))):
                {dummy_atom_idx: [(bond_atom_idx, bond_atom_atomicnum), ...], }
        """
        if isinstance(mol,str): mol = Chem.MolFromMolFile(mol)
        if not isinstance(mol,Chem.Mol): return {}
        dummyidx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        bondedinfo = collections.OrderedDict()
        for d in dummyidx:
            bondedinfo[d] = []
            for a in mol.GetAtomWithIdx(d).GetBonds():
                e = a.GetBeginAtom()
                t = e.GetIdx()
                if t == d:
                    e = a.GetEndAtom()
                    t = e.GetIdx()
                bondedinfo[d].append((t,e.GetAtomicNum()))
        return bondedinfo


    def write(self,num=None):
        w = Chem.SDWriter(self.filename)
        for i in mols:
            w.write(i)
        w.close()

        print('DONE')


    def run(self):
        if not self.nice: raise StopIteration
        for n in range(1,self.core_dummynum+1):
            for p in itertools.combinations(range(self.core_dummynum),r=n):
                for g in itertools.permutations(range(n),r=n):
                    # caution: concatenating on sequence
                    fullsmarts = self.core_smarts + '.' + '.'.join([self.subs_smartss[t] for t in g])
                    mol = Chem.MolFromSmarts(fullsmarts)
                    fullbondedinfo = self.calc_dummy_atoms_idx_and_bonded_info(mol)
                    fullitems = list(fullbondedinfo.items())
                    bo = True
                    for i in range(n):
                        beg = fullitems[p[i]]
                        end = fullitems[g[i]+self.core_dummynum]
                        for z in beg[1]:
                            for k in end[1]:
                                if z[1] == 8 and k[1] == 8:     # not allowed for O-O bond
                                    bo = False
                                    break
                            if not bo:
                                break
                        if not bo:
                            break
                    if bo:
                        rwmol = Chem.RWMol(mol)
                        for i in range(n):
                            beg = fullitems[p[i]]
                            end = fullitems[g[i]+self.core_dummynum]
                            rwmol.AddBond(beg[0],end[0],Chem.BondType.SINGLE)
                        # now replace left dummy atom to hydrogen
                        if n < self.core_dummynum:
                            left = [i for i in range(self.core_dummynum) if i not in p]
                            for t in left:
                                hydrogen = Chem.Atom('H')
                                rwmol.ReplaceAtom(fullitems[t][0],hydrogen)

                        # dummy atoms (will always be pair) are connected together
                        bondedinfo = self.calc_dummy_atoms_idx_and_bonded_info(rwmol.GetMol())
                        totalitems = list(bondedinfo.items())
                        visited = set()
                        reflist = []
                        for i in range(len(totalitems)):
                            if i in visited: continue
                            for j in range(i+1,len(totalitems)):
                                if totalitems[i][0] in [t[0] for t in totalitems[j][1]]:
                                    visited.add(j)
                                    reflist.append((i,j))
                                    break
                        for t in reflist:
                            newitems = [totalitems[t[0]], totalitems[t[1]]]
                            # caution: mutual bond,
                            # for example: a is the neighbor b, b is also the neighbor of a
                            beg = [i[0] for i in newitems[0][1] if i[0] != newitems[1][0]]
                            endatomidx = [i[0] for i in newitems[1][1] if i[0] != newitems[0][0]][0]
                            for i in beg:
                                rwmol.AddBond(i,endatomidx,Chem.BondType.SINGLE)

                        # pay attention on sequence, from the end to beginning
                        for i in range(len(totalitems)-1,-1,-1):
                            rwmol.RemoveAtom(totalitems[i][0])

                        yield rwmol.GetMol()


cn = CoreEnumeration(core_file='imp-core-sub-O.sdf',sub_files='imp-sub.sdf')

mols = [next(cn) for i in range(100)]

