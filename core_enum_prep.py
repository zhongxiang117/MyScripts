#!/usr/bin/env python3

from rdkit import Chem

import os
import itertools
import collections

FEATURES = [
    'version 0.1.0  : CoreEnumeration, Sep 19th, 2022',
    'version 0.2.0  : add CoreEnumeration.write method',
    'version 0.3.0  : make write more powerful',
]

VERSION = FEATURES[-1].split(':')[0].replace('version',' ').strip()


class CoreEnumeration:
    def __init__(self,core_file=None,sub_files=None,filename=None,*args,**kws):
        if filename:
            self.filename = filename
            self.filebase = os.path.splitext(filename)[0]
        else:
            self.filename = 'mols-core-enum.sdf'
            self.filebase = 'mols-core-enum'
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
            self.subs_smartss = [Chem.MolToSmarts(subs[i]) for i in reflist]
            self.subs_dummynum = len(self.subs_smartss)
            if not self.core_dummynum:
                print('Fatal: core_file: no dummy atoms found: exiting')
                self.nice = False
            if not self.subs_dummynum:
                print('Fatal: sub_files: no dummy atoms found: exiting')
                self.nice = False
            if self.nice:
                self._it = self._run()
        else:
            print(f'Fatal: not a file: exiting: {core_file}')
            self.nice = False

    def gen_new_file(self,filename=None):
        if not filename: filename = self.filebase + '.sdf'
        if os.path.isfile(filename):
            base,ext = os.path.splitext(filename)
            i = 1
            while True:
                new = f'{base}-{i}{ext}'
                if not os.path.isfile(new): break
                i += 1
        else:
            new = filename
        return new

    def __next__(self):
        if not self.nice: raise StopIteration
        return next(self._it)

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

    def _precondition(self,mols=None,num=None):
        if mols:
            if isinstance(mols,list):
                mols = [m for m in mols if isinstance(m, Chem.Mol)]
            elif isinstance(mols, Chem.Mol):
                mols = [mols, ]
            else:
                print('Fatal: wrong inputs: mols')
                mols = []
        else:
            mols = []
            if self.nice:
                if not num: num = 10
                for i in range(num):
                    try:
                        mols.append(next(self._it))
                    except StopIteration:
                        break
        return mols
    
    def write_sdf_lazyeval(self,num=None,filename=None):
        """write number of "self.mols" to SDF filename, lazy evaluation"""
        if not filename: filename = self.gen_new_file(self.filebase+'.sdf')
        w = Chem.SDWriter(filename)
        if num:
            print(f'Note: writing {num} mols to SDF {filename}')
            for i in range(num):
                try:
                    m = next(self._it)
                except StopIteration:
                    break
                else:
                    w.write(m)
        else:
            print(f'Note: writing all mols to SDF {filename}')
            while True:
                try:
                    m = next(self._it)
                except StopIteration:
                    break
                else:
                    w.write(m)
        w.close()

    def write_smiles_lazyeval(self,num=None,filename=None):
        """write number of "self.mols" to SMILES filename, lazy evaluation"""
        if not filename: filename = self.gen_new_file(self.filebase+'.smiles')
        w = open(filename,'wt')
        if num:
            print(f'Note: writing {num} SMILES to {filename}')
            for i in range(num):
                try:
                    m = next(self._it)
                except StopIteration:
                    break
                else:
                    w.write(Chem.MolToSmiles(m))
                    w.write('\n')
        else:
            print(f'Note: writing all SMILES to {filename}')
            while True:
                try:
                    m = next(self._it)
                except StopIteration:
                    break
                else:
                    w.write(Chem.MolToSmiles(m))
                    w.write('\n')
        w.close()

    def write_smiles(self,mols=None,num=None,filename=None):
        """write mols or number of "self.mols" to SMILES filename"""
        if not filename: filename = self.gen_new_file(self.filebase+'.smiles')
        mols = self._precondition(mols,num)
        if mols:
            print(f'Note: writing {len(mols)} SMILES to {filename}')
            smiles = [Chem.MolToSmiles(m) for m in mols]
            with open(filename,'wt') as f:
                for s in smiles:
                    f.write(s)
                    f.write('\n')
        else:
            print('Fatal: write: no valid inputs: exit')

    def write_sdf(self,mols=None,num=None,filename=None):
        """write mols or number of "self.mols" to SDF filename"""
        if not filename: filename = self.gen_new_file(self.filebase+'.sdf')
        mols = self._precondition(mols,num)
        if mols:
            print(f'Note: writing {len(mols)} mols to SDF {filename}')
            w = Chem.SDWriter(filename)
            for i in mols:
                w.write(i)
            w.close()
        else:
            print('Fatal: write: no valid inputs: exit')
    
    def _run(self):
        """chain all iterators, python>=3.5"""
        for n in range(1,self.core_dummynum+1):
            yield from self.run_subn(n)

    def run_subn(self,n,core_smarts=None,subs_smartss=None):
        """substitute n dummy atoms on core molecule, lazy evaluation"""
        if core_smarts:
            core_dummynum = len(self.calc_dummy_atoms_idx_and_bonded_info(core_smarts))
        else:
            core_smarts = self.core_smarts
            core_dummynum = self.core_dummynum
        if subs_smartss:
            subs_dummynum = len(subs_smartss)
        else:
            subs_smartss = self.subs_smartss
            subs_dummynum = self.subs_dummynum
        for p in itertools.combinations(range(core_dummynum),r=n):
            for g in itertools.permutations(range(subs_dummynum),r=n):
                # caution: concatenating on sequence
                fullsmarts = core_smarts + '.' + '.'.join([subs_smartss[t] for t in g])
                mol = Chem.MolFromSmarts(fullsmarts)
                fullbondedinfo = self.calc_dummy_atoms_idx_and_bonded_info(mol)
                fullitems = list(fullbondedinfo.items())
                bo = True
                for i in range(n):
                    beg = fullitems[p[i]]
                    end = fullitems[i+core_dummynum]   # important: g's are in sequence
                    for z in beg[1]:
                        for k in end[1]:
                            if z[1] == 8 and k[1] == 8:     # not allowed for Oxygen-Oxygen bond
                                bo = False
                                break
                        if not bo:
                            break
                    if not bo:
                        break
                if not bo: continue

                rwmol = Chem.RWMol(mol)
                for i in range(n):
                    beg = fullitems[p[i]]
                    end = fullitems[i+core_dummynum]   # important: g's are in sequence
                    rwmol.AddBond(beg[0],end[0],Chem.BondType.SINGLE)
                # now replace left dummy atom to hydrogen
                if n < core_dummynum:
                    left = [i for i in range(core_dummynum) if i not in p]
                    for t in left:
                        hydrogen = Chem.Atom('H')
                        rwmol.ReplaceAtom(fullitems[t][0],hydrogen)
                    lefthidx = [fullitems[t][0] for t in left]
                else:
                    lefthidx = []

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
                leftdummyidx = [i[0] for i in totalitems]
                for i in sorted(lefthidx+leftdummyidx,reverse=True):
                    rwmol.RemoveAtom(i)

                yield rwmol.GetMol()

    def filter(self,begin_bondedinfo,end_bondedinfo):
        pass

ce = CoreEnumeration(core_file='imp-core-sub-N.sdf',sub_files='imp-sub.sdf')
#ce._run()
ce.write_sdf_lazyeval(num=10000000)



